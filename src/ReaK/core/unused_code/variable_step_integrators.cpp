
/*
 *    Copyright 2011 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */


#include "variable_step_integrators.hpp"

#include <cmath>

namespace ReaK {




//-----------CAdamsBMVar--------------------------------------------------------

/*
    Bashforth Coefficients
    p | k | j->    |  1     2     3     4     5     6      |  CP
    -------------------------------------------------------------------
    1 | 1 | Bj     |  1                                    | 1/2
    2 | 2 | 2Bj    |  3    -1                              | 5/12
    3 | 3 | 12Bj   |  23   -16    5                        | 3/8
    4 | 4 | 24Bj   |  55   -59    37   -9                  | 251/720
    5 | 5 | 720Bj  |  1901 -2774  2616 -1274  251          | 95/288
    6 | 6 | 1440Bj |  4277 -7923  9982 -7298  2877 -475    | 19087/60480

    Yn = Yn-1 + h*sum(j=1->k, Bj*Fn-j);

    Moulton Coefficients
    p | k | j->    |  0     1     2     3     4     5      |  CC
    -------------------------------------------------------------------
    1 | 1 | Bj     |  1                                    | -1/2
    2 | 1 | 2Bj    |  1     1                              | -1/12
    3 | 2 | 12Bj   |  5     8    -1                        | -1/24
    4 | 3 | 24Bj   |  9     19   -5     1                  | -19/720
    5 | 4 | 720Bj  |  251   646  -264   106  -19           | -3/160
    6 | 5 | 1440Bj |  475   1427 -798   482  -173   27     | -863/60480

    Yn = Yn-1 + h*sum(j=0->k, Bj*Fn-j); //note: implicit, hence Yn or Fn is not known but predicted using Bashforth

    Predictor-Corrector Steps:
    Predictor : compute (Yn)0 using Bashforth
    Evaluation : evaluate (Fn)0 using (Yn)0
    Corrector : compute (Yn)1 using Moulton  -> iterate again Eval/Corr by a fixed number of times "Corrections"
    Evaluation : evaluate (Fn)1 using (Yn)1 for the value of Fn-1 of the next time step

    error = (den(CC) * den(CP) * num(CC)) * abs((Yn)C - (Yn)0) / ((num(CP) * den(CC) + num(CC) * den(CP)) * den(CC) * h)
      where num: absolute value of numerator of ()  and den: absolute value of denominator of ()

    Order 1: error = abs((Yn)C - (Yn)0) / (2.0 * h);
    Order 2: error = abs((Yn)C - (Yn)0) / (6.0 * h);
    Order 3: error = abs((Yn)C - (Yn)0) / (10.0 * h);
    Order 4: error = abs((Yn)C - (Yn)0) * 19.0 / (270.0 * h);
    Order 5: error = abs((Yn)C - (Yn)0) * 27.0 / (502.0 * h);
    Order 6: error = abs((Yn)C - (Yn)0) * 863.0 / (19950.0 * h);

    hnew = h * q; for 0.1 < q < 4;
    q = (Tolerance / (2 * error)) ^ (1/Order)
  */

  //TODO: This is only the copy/paste of the fixed-time step one.

/* DEPRECATED
template <class T>
int CAdamsBMVar<T>::Integrate() {
  if ((GetStateRate == NULL) || (_StateCount == 0) || (StepSize == T(0.0)) || ((StepSize > T(0.0)) && (EndTime < StartTime)) || ((StepSize < T(0.0)) && (StartTime < EndTime)) || (Tolerance <= T(0.0)) || (MinStepSize > MaxStepSize)) return -1;

  if(Corrections == 0) Corrections = 1;
  if(Order == 0) Order = 1;
  else if(Order > 6) Order = 6;

  T t = StartTime;
  T E, Emax;
  T* prevY = new T[_StateCount];
  T* Y0 = new T[_StateCount];
  if(Order == 1) {

    GetStateRate(t,GetStateRate_UserData);
    if(OutputState != NULL) OutputState(t,OutputState_UserData);

    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += StepSize * _StateRate[i][0];
        Y0[i] = _State[i][0];
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * _StateRate[i][0];
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs((_State[i][0] - Y0[i])/(T(2.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = Tolerance / (T(2.0) * Emax);
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = Tolerance / (T(2.0) * Emax);
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
      };
    };
  } else if(Order == 2) {
    T* prevF = new T[2*_StateCount];

    GetStateRate(t,GetStateRate_UserData);
    if (OutputState != NULL) OutputState(t,OutputState_UserData);

    // Runge-Kutta-Fehlberg Iterations for the first two points
    T* K = new T[6*_StateCount];

    while(((StepSize > T(0.0)) && (t <= StartTime + T(2.0) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(2.0) * StepSize))) {
      if(fabs(StepSize) < MinStepSize) StepSize *= fabs(MinStepSize / StepSize);
      if(fabs(StepSize) > MaxStepSize) StepSize *= fabs(MaxStepSize / StepSize);

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += (K[6*i] = StepSize * _StateRate[i][0]) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(3.0) * K[6*i] + T(9.0) * (K[6*i+1] = StepSize * _StateRate[i][0])) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(1932.0) * K[6*i] - T(7200.0) * K[6*i+1] + T(7296.0) * (K[6*i+2] = StepSize * _StateRate[i][0])) / T(2197.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * K[6*i] / T(216.0) - T(8.0) * K[6*i+1] + T(3680.0) * K[6*i+2] / T(513.0) - T(845.0) * (K[6*i+3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * K[6*i] / T(27.0) + T(2.0) * K[6*i+1] - T(3544.0) * K[6*i+2] / T(2565.0) + T(1859.0) * K[6*i+3] / T(4104.0) - T(11.0) * (K[6*i+4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        K[6*i+5] = StepSize * _StateRate[i][0];

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(K[6*i] / T(360.0) - T(128.0) * K[6*i+2] / T(4275.0) - T(2197.0) * K[6*i+3] / T(75240.0) + K[6*i+4] / T(50.0) + T(2.0) * K[6*i+5] / T(55.0)) / StepSize;
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {
        t += StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + T(25.0) * K[6*i] / T(216.0) + T(1408.0) * K[6*i+2] / T(2565.0) + T(2197.0) * K[6*i+3] / T(4104.0) - K[6*i+4] / T(5.0);
        GetStateRate(t,GetStateRate_UserData);
        if(OutputState != NULL) OutputState(t,OutputState_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E >= T(4.0)) StepSize *= T(4.0); else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
      } else {
        t -= StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        GetStateRate(t,GetStateRate_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E < T(0.1)) StepSize *= T(0.1); else StepSize *= E;
      };
    };

    //Midpoint backward iterations
    StepSize = -StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      Y0[i] = _State[i][0];
    while(((StepSize > T(0.0)) && (t < StartTime + T(1.1) * StepSize)) || ((StepSize < T(0.0)) && (t > StartTime + T(1.1) * StepSize))) {

      for(unsigned int i=0;i<_StateCount;i++)
      {
        prevY[i] = _State[i][0];
        prevF[2*i+1] = prevF[2*i];
        _State[i][0] += (prevF[2*i] = _StateRate[i][0]) * StepSize / T(2.0);
      };

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + StepSize * _StateRate[i][0];

      t += StepSize / T(2.0);

      GetStateRate(t,GetStateRate_UserData);
    };

    //Back to Forward Iterations
    StepSize = -StepSize;
    t += T(2.0) * StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      _State[i][0] = Y0[i];
    GetStateRate(t,GetStateRate_UserData);
    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        prevF[2*i] = prevF[2*i+1];
        _State[i][0] += StepSize * (T(3.0) * (prevF[2*i+1] = _StateRate[i][0]) - prevF[2*i]) / T(2.0);
        Y0[i] = _State[i][0];
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * (_StateRate[i][0] + prevF[2*i+1]) / T(2.0);
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs((_State[i][0] - Y0[i])/(T(6.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = pow(Tolerance / (T(2.0) * Emax),T(0.5));
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t - T(2.0) * StepSize < StartTime)) || ((StepSize < T(0.0)) && (t - T(2.0) * StepSize > StartTime))) StepSize = (t - StartTime) / T(2.0);
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;

          //Backward Midpoint iterations for 2 points
          StepSize = -StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            Y0[i] = _State[i][0];
          for(unsigned int j=0;j<2;j++) {

            for(unsigned int i=0;i<_StateCount;i++)
            {
              prevY[i] = _State[i][0];
              prevF[2*i+1] = prevF[2*i];
              _State[i][0] += (prevF[2*i] = _StateRate[i][0]) * StepSize / T(2.0);
            };

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + StepSize * _StateRate[i][0];

            t += StepSize / T(2.0);

            GetStateRate(t,GetStateRate_UserData);
          };
          StepSize = -StepSize;
          t += T(2.0) * StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = Y0[i];
          GetStateRate(t,GetStateRate_UserData);
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = pow(Tolerance / (T(2.0) * Emax),T(0.5));
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;

        //Backward Midpoint iterations for 2 points
        StepSize = -StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          Y0[i] = _State[i][0];
        for(unsigned int j=0;j<2;j++) {

          for(unsigned int i=0;i<_StateCount;i++)
          {
            prevY[i] = _State[i][0];
            prevF[2*i+1] = prevF[2*i];
            _State[i][0] += (prevF[2*i] = _StateRate[i][0]) * StepSize / T(2.0);
          };

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + StepSize * _StateRate[i][0];

          t += StepSize / T(2.0);

          GetStateRate(t,GetStateRate_UserData);
        };
        StepSize = -StepSize;
        t += T(2.0) * StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = Y0[i];
        GetStateRate(t,GetStateRate_UserData);
      };
    };
  } else if(Order == 3) {
    T* prevF = new T[3*_StateCount];

    GetStateRate(t,GetStateRate_UserData);
    if (OutputState != NULL) OutputState(t,OutputState_UserData);

    // Runge-Kutta-Fehlberg Iterations for the first 3 points
    {
    T* K = new T[6*_StateCount];

    while(((StepSize > T(0.0)) && (t <= StartTime + T(3.0) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(3.0) * StepSize))) {
      if(fabs(StepSize) < MinStepSize) StepSize *= fabs(MinStepSize / StepSize);
      if(fabs(StepSize) > MaxStepSize) StepSize *= fabs(MaxStepSize / StepSize);

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += (K[6*i] = StepSize * _StateRate[i][0]) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(3.0) * K[6*i] + T(9.0) * (K[6*i+1] = StepSize * _StateRate[i][0])) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(1932.0) * K[6*i] - T(7200.0) * K[6*i+1] + T(7296.0) * (K[6*i+2] = StepSize * _StateRate[i][0])) / T(2197.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * K[6*i] / T(216.0) - T(8.0) * K[6*i+1] + T(3680.0) * K[6*i+2] / T(513.0) - T(845.0) * (K[6*i+3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * K[6*i] / T(27.0) + T(2.0) * K[6*i+1] - T(3544.0) * K[6*i+2] / T(2565.0) + T(1859.0) * K[6*i+3] / T(4104.0) - T(11.0) * (K[6*i+4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        K[6*i+5] = StepSize * _StateRate[i][0];

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(K[6*i] / T(360.0) - T(128.0) * K[6*i+2] / T(4275.0) - T(2197.0) * K[6*i+3] / T(75240.0) + K[6*i+4] / T(50.0) + T(2.0) * K[6*i+5] / T(55.0)) / StepSize;
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {
        t += StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + T(25.0) * K[6*i] / T(216.0) + T(1408.0) * K[6*i+2] / T(2565.0) + T(2197.0) * K[6*i+3] / T(4104.0) - K[6*i+4] / T(5.0);
        GetStateRate(t,GetStateRate_UserData);
        if(OutputState != NULL) OutputState(t,OutputState_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E >= T(4.0)) StepSize *= T(4.0); else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
      } else {
        t -= StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        GetStateRate(t,GetStateRate_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E < T(0.1)) StepSize *= T(0.1); else StepSize *= E;
      };
    };
    };

    //Backward Runge-Kutta order 4 for the first 3 points
    StepSize = -StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      Y0[i] = _State[i][0];
    {
    T* k = new T[3*_StateCount];
    while(((StepSize > T(0.0)) && (t <= StartTime + T(2.1) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(2.1) * StepSize))) {

      for(unsigned int i=0;i<_StateCount;i++)
      {
        prevY[i] = _State[i][0];
        prevF[3*i+2] = prevF[3*i+1];
        prevF[3*i+1] = prevF[3*i];
        _State[i][0] += (k[3*i] = StepSize * (prevF[3*i] = _StateRate[i][0])) / T(2.0);
      };

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

      GetStateRate(t,GetStateRate_UserData);
    };
    };

    //Back To Forward iterations
    StepSize = -StepSize;
    t += T(3.0) * StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      _State[i][0] = Y0[i];
    GetStateRate(t,GetStateRate_UserData);
    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        prevF[3*i] = prevF[3*i+1];
        prevF[3*i+1] = prevF[3*i+2];
        _State[i][0] += StepSize * (T(23.0) * (prevF[3*i+2] = _StateRate[i][0]) - T(16.0) * prevF[3*i+1] + T(5.0) * prevF[3*i]) / T(12.0);
        Y0[i] = _State[i][0];
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * (T(5.0) * _StateRate[i][0] + T(8.0) * prevF[3*i+2] - prevF[3*i+1]) / T(12.0);
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs((_State[i][0] - Y0[i])/(T(10.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = pow(Tolerance / (T(2.0) * Emax),T(1.0/3.0));
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t - T(3.0) * StepSize < StartTime)) || ((StepSize < T(0.0)) && (t - T(3.0) * StepSize > StartTime))) StepSize = (t - StartTime) / T(3.0);
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;

          //Backward Runge-Kutta Order 4 iterations for 3 points
          StepSize = -StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            Y0[i] = _State[i][0];
          {
          T* k = new T[3*_StateCount];
          for(unsigned int j=0;j<3;j++) {

            for(unsigned int i=0;i<_StateCount;i++)
            {
              prevY[i] = _State[i][0];
              prevF[3*i+2] = prevF[3*i+1];
              prevF[3*i+1] = prevF[3*i];
              _State[i][0] += (k[3*i] = StepSize * (prevF[3*i] = _StateRate[i][0])) / T(2.0);
            };

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

            GetStateRate(t,GetStateRate_UserData);
          };
          };
          StepSize = -StepSize;
          t += T(3.0) * StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = Y0[i];
          GetStateRate(t,GetStateRate_UserData);
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = pow(Tolerance / (T(2.0) * Emax),T(1.0/3.0));
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;

        //Backward Runge-Kutta Order 4 iterations for 3 points
        StepSize = -StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          Y0[i] = _State[i][0];
        {
        T* k = new T[3*_StateCount];
        for(unsigned int j=0;j<3;j++) {

          for(unsigned int i=0;i<_StateCount;i++)
          {
            prevY[i] = _State[i][0];
            prevF[3*i+2] = prevF[3*i+1];
            prevF[3*i+1] = prevF[3*i];
            _State[i][0] += (k[3*i] = StepSize * (prevF[3*i] = _StateRate[i][0])) / T(2.0);
          };

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

          GetStateRate(t,GetStateRate_UserData);
        };
        };
        StepSize = -StepSize;
        t += T(3.0) * StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = Y0[i];
        GetStateRate(t,GetStateRate_UserData);
      };
    };
  } else if(Order == 4) {
    T* prevF = new T[4*_StateCount];

    GetStateRate(t,GetStateRate_UserData);
    if (OutputState != NULL) OutputState(t,OutputState_UserData);

    // Runge-Kutta-Fehlberg Iterations for the first 4 points
    {
    T* K = new T[6*_StateCount];

    while(((StepSize > T(0.0)) && (t <= StartTime + T(4.0) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(4.0) * StepSize))) {
      if(fabs(StepSize) < MinStepSize) StepSize *= fabs(MinStepSize / StepSize);
      if(fabs(StepSize) > MaxStepSize) StepSize *= fabs(MaxStepSize / StepSize);

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += (K[6*i] = StepSize * _StateRate[i][0]) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(3.0) * K[6*i] + T(9.0) * (K[6*i+1] = StepSize * _StateRate[i][0])) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(1932.0) * K[6*i] - T(7200.0) * K[6*i+1] + T(7296.0) * (K[6*i+2] = StepSize * _StateRate[i][0])) / T(2197.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * K[6*i] / T(216.0) - T(8.0) * K[6*i+1] + T(3680.0) * K[6*i+2] / T(513.0) - T(845.0) * (K[6*i+3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * K[6*i] / T(27.0) + T(2.0) * K[6*i+1] - T(3544.0) * K[6*i+2] / T(2565.0) + T(1859.0) * K[6*i+3] / T(4104.0) - T(11.0) * (K[6*i+4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        K[6*i+5] = StepSize * _StateRate[i][0];

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(K[6*i] / T(360.0) - T(128.0) * K[6*i+2] / T(4275.0) - T(2197.0) * K[6*i+3] / T(75240.0) + K[6*i+4] / T(50.0) + T(2.0) * K[6*i+5] / T(55.0)) / StepSize;
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {
        t += StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + T(25.0) * K[6*i] / T(216.0) + T(1408.0) * K[6*i+2] / T(2565.0) + T(2197.0) * K[6*i+3] / T(4104.0) - K[6*i+4] / T(5.0);
        GetStateRate(t,GetStateRate_UserData);
        if(OutputState != NULL) OutputState(t,OutputState_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E >= T(4.0)) StepSize *= T(4.0); else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
      } else {
        t -= StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        GetStateRate(t,GetStateRate_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E < T(0.1)) StepSize *= T(0.1); else StepSize *= E;
      };
    };
    };

    //Backward Runge-Kutta order 4 for the first 4 points
    StepSize = -StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      Y0[i] = _State[i][0];
    {
    T* k = new T[3*_StateCount];
    while(((StepSize > T(0.0)) && (t <= StartTime + T(3.1) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(3.1) * StepSize))) {

      for(unsigned int i=0;i<_StateCount;i++)
      {
        prevY[i] = _State[i][0];
        prevF[4*i+3] = prevF[4*i+2];
        prevF[4*i+2] = prevF[4*i+1];
        prevF[4*i+1] = prevF[4*i];
        _State[i][0] += (k[3*i] = StepSize * (prevF[4*i] = _StateRate[i][0])) / T(2.0);
      };

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

      GetStateRate(t,GetStateRate_UserData);
    };
    };

    //Back To Forward iterations
    StepSize = -StepSize;
    t += T(4.0) * StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      _State[i][0] = Y0[i];
    GetStateRate(t,GetStateRate_UserData);
    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        prevF[4*i] = prevF[4*i+1];
        prevF[4*i+1] = prevF[4*i+2];
        prevF[4*i+2] = prevF[4*i+3];
        _State[i][0] += StepSize * (T(55.0) * (prevF[4*i+3] = _StateRate[i][0]) - T(59.0) * prevF[4*i+2] + T(37.0) * prevF[4*i+1] - T(9.0) * prevF[4*i]) / T(24.0);
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * (T(9.0) * _StateRate[i][0] + T(19.0) * prevF[4*i+3] - T(5.0) * prevF[4*i+2] + prevF[4*i+1]) / T(24.0);
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(T(19.0) * (_State[i][0] - Y0[i])/(T(270.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = pow(Tolerance.q / (T(2.0) * Emax),T(0.25));
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t - T(4.0) * StepSize < StartTime)) || ((StepSize < T(0.0)) && (t - T(4.0) * StepSize > StartTime))) StepSize = (t - StartTime) / T(4.0);
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;

          //Backward Runge-Kutta Order 4 iterations for 4 points
          StepSize = -StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            Y0[i] = _State[i][0];
          {
          T* k = new T[3*_StateCount];
          for(unsigned int j=0;j<4;j++) {

            for(unsigned int i=0;i<_StateCount;i++)
            {
              prevY[i] = _State[i][0];
              prevF[4*i+3] = prevF[4*i+2];
              prevF[4*i+2] = prevF[4*i+1];
              prevF[4*i+1] = prevF[4*i];
              _State[i][0] += (k[3*i] = StepSize * (prevF[4*i] = _StateRate[i][0])) / T(2.0);
            };

            t += StepSize / 2.0;
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

            GetStateRate(t,GetStateRate_UserData);
          };
          };
          StepSize = -StepSize;
          t += T(4.0) * StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = Y0[i];
          GetStateRate(t,GetStateRate_UserData);
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = pow(Tolerance / (T(2.0) * Emax),T(0.25));
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;

        //Backward Runge-Kutta Order 4 iterations for 4 points
        StepSize = -StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          Y0[i] = _State[i][0];
        {
        T* k = new T[3*_StateCount];
        for(unsigned int j=0;j<4;j++) {

          for(unsigned int i=0;i<_StateCount;i++)
          {
            prevY[i] = _State[i][0];
            prevF[4*i+3] = prevF[4*i+2];
            prevF[4*i+2] = prevF[4*i+1];
            prevF[4*i+1] = prevF[4*i];
            _State[i][0] += (k[3*i] = StepSize * (prevF[4*i] = _StateRate[i][0])) / T(2.0);
          };

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + (k[3*i + 1] = StepSize * _StateRate[i][0]) / T(2.0);

          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + (k[3*i + 2] = StepSize * _StateRate[i][0]);

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (k[3*i] + T(2.0) * k[3*i + 1] + StepSize * _StateRate[i][0]) / T(6.0) - T(2.0) * k[3*i + 2] / T(3.0);

          GetStateRate(t,GetStateRate_UserData);
        };
        };
        StepSize = -StepSize;
        t += T(4.0) * StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = Y0[i];
        GetStateRate(t,GetStateRate_UserData);
      };
    };
  } else if(Order == 5) {
    T* prevF = new T[5*_StateCount];

    GetStateRate(t,GetStateRate_UserData);
    if(OutputState != NULL) OutputState(t,OutputState_UserData);

    // Runge-Kutta-Fehlberg Iterations for the first 5 points
    {
    T* K = new T[6*_StateCount];

    while(((StepSize > T(0.0)) && (t <= StartTime + T(5.0) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(5.0) * StepSize))) {
      if(fabs(StepSize) < MinStepSize) StepSize *= fabs(MinStepSize / StepSize);
      if(fabs(StepSize) > MaxStepSize) StepSize *= fabs(MaxStepSize / StepSize);

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += (K[6*i] = StepSize * _StateRate[i][0]) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(3.0) * K[6*i] + T(9.0) * (K[6*i+1] = StepSize * _StateRate[i][0])) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(1932.0) * K[6*i] - T(7200.0) * K[6*i+1] + T(7296.0) * (K[6*i+2] = StepSize * _StateRate[i][0])) / T(2197.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * K[6*i] / T(216.0) - T(8.0) * K[6*i+1] + T(3680.0) * K[6*i+2] / T(513.0) - T(845.0) * (K[6*i+3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * K[6*i] / T(27.0) + T(2.0) * K[6*i+1] - T(3544.0) * K[6*i+2] / T(2565.0) + T(1859.0) * K[6*i+3] / T(4104.0) - T(11.0) * (K[6*i+4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        K[6*i+5] = StepSize * _StateRate[i][0];

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(K[6*i] / T(360.0) - T(128.0) * K[6*i+2] / T(4275.0) - T(2197.0) * K[6*i+3] / T(75240.0) + K[6*i+4] / T(50.0) + T(2.0) * K[6*i+5] / T(55.0)) / StepSize;
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {
        t += StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + T(25.0) * K[6*i] / T(216.0) + T(1408.0) * K[6*i+2] / T(2565.0) + T(2197.0) * K[6*i+3] / T(4104.0) - K[6*i+4] / T(5.0);
        GetStateRate(t,GetStateRate_UserData);
        if(OutputState != NULL) OutputState(t,OutputState_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E >= T(4.0)) StepSize *= T(4.0); else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
      } else {
        t -= StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        GetStateRate(t,GetStateRate_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E < T(0.1)) StepSize *= T(0.1); else StepSize *= E;
      };
    };
    };

    //Runge-Kutta order 5 for the first 5 points
    StepSize = -StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      Y0[i] = _State[i][0];
    {
    T* k = new T[5*_StateCount];
    while(((StepSize > T(0.0)) && (t <= StartTime + T(4.1) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(4.1) * StepSize))) {

      for(unsigned int i=0;i<_StateCount;i++)
      {
        prevY[i] = _State[i][0];
        prevF[5*i+4] = prevF[5*i+3];
        prevF[5*i+3] = prevF[5*i+2];
        prevF[5*i+2] = prevF[5*i+1];
        prevF[5*i+1] = prevF[5*i];
        _State[i][0] += (k[5*i] = StepSize * (prevF[5*i] = _StateRate[i][0])) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
    };
    };

    //Back to Forward Iterations
    StepSize = -StepSize;
    t += T(5.0) * StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      _State[i][0] = Y0[i];
    GetStateRate(t,GetStateRate_UserData);
    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        prevF[5*i] = prevF[5*i+1];
        prevF[5*i+1] = prevF[5*i+2];
        prevF[5*i+2] = prevF[5*i+3];
        prevF[5*i+3] = prevF[5*i+4];
        _State[i][0] += StepSize * (T(1901.0) * (prevF[5*i+4] = _StateRate[i][0]) - T(2774.0) * prevF[5*i+3] + T(2616.0) * prevF[5*i+2] - T(1274.0) * prevF[5*i+1] + T(251.0) * prevF[5*i]) / T(720.0);
        Y0[i] = _State[i][0];
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * (T(251.0) * _StateRate[i][0] + T(646.0) * prevF[5*i+4] - T(264.0) * prevF[5*i+3] + T(106.0) * prevF[5*i+2] - T(19.0) * prevF[5*i+1]) / T(720.0);
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(T(27.0) * (_State[i][0] - Y0[i])/(T(502.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = pow(Tolerance / (T(2.0) * Emax),T(0.2));
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t - T(5.0) * StepSize < StartTime)) || ((StepSize < T(0.0)) && (t - T(5.0) * StepSize > StartTime))) StepSize = (t - StartTime) / T(5.0);
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;

          //Backward Runge-Kutta Order 5 iterations for 5 points
          StepSize = -StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            Y0[i] = _State[i][0];
          {
          T* k = new T[5*_StateCount];
          for(unsigned int j=0;j<5;j++) {

            for(unsigned int i=0;i<_StateCount;i++)
            {
              prevY[i] = _State[i][0];
              prevF[5*i+4] = prevF[5*i+3];
              prevF[5*i+3] = prevF[5*i+2];
              prevF[5*i+2] = prevF[5*i+1];
              prevF[5*i+1] = prevF[5*i];
              _State[i][0] += (k[5*i] = StepSize * (prevF[5*i] = _StateRate[i][0])) / T(4.0);
            };

            t += StepSize / T(4.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

            t += StepSize / T(8.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

            t += T(57.0) * StepSize / T(104.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

            t += StepSize / T(13.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

            t -= StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
          };
          };
          StepSize = -StepSize;
          t += T(5.0) * StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = Y0[i];
          GetStateRate(t,GetStateRate_UserData);
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = pow(Tolerance / (T(2.0) * Emax),T(0.2));
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;

        //Backward Runge-Kutta Order 5 iterations for 5 points
        StepSize = -StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          Y0[i] = _State[i][0];
        {
        T* k = new T[5*_StateCount];
        for(unsigned int j=0;j<5;j++) {

          for(unsigned int i=0;i<_StateCount;i++)
          {
            prevY[i] = _State[i][0];
            prevF[5*i+4] = prevF[5*i+3];
            prevF[5*i+3] = prevF[5*i+2];
            prevF[5*i+2] = prevF[5*i+1];
            prevF[5*i+1] = prevF[5*i];
            _State[i][0] += (k[5*i] = StepSize * (prevF[5*i] = _StateRate[i][0])) / T(4.0);
          };

          t += StepSize / T(4.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

          t += StepSize / T(8.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

          t += T(57.0) * StepSize / T(104.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

          t += StepSize / T(13.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

          t -= StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
        };
        };
        StepSize = -StepSize;
        t += T(5.0) * StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = Y0[i];
        GetStateRate(t,GetStateRate_UserData);
      };
    };
  } else if(Order == 6) {
    T* prevF = new T[6*_StateCount];

    GetStateRate(t,GetStateRate_UserData);
    if(OutputState != NULL) OutputState(t,OutputState_UserData);

    // Runge-Kutta-Fehlberg Iterations for the first 6 points
    {
    T* K = new T[6*_StateCount];

    while(((StepSize > T(0.0)) && (t <= StartTime + T(6.0) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(6.0) * StepSize))) {
      if(fabs(StepSize) < MinStepSize) StepSize *= fabs(MinStepSize / StepSize);
      if(fabs(StepSize) > MaxStepSize) StepSize *= fabs(MaxStepSize / StepSize);

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        _State[i][0] += (K[6*i] = StepSize * _StateRate[i][0]) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(3.0) * K[6*i] + T(9.0) * (K[6*i+1] = StepSize * _StateRate[i][0])) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + (T(1932.0) * K[6*i] - T(7200.0) * K[6*i+1] + T(7296.0) * (K[6*i+2] = StepSize * _StateRate[i][0])) / T(2197.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * K[6*i] / T(216.0) - T(8.0) * K[6*i+1] + T(3680.0) * K[6*i+2] / T(513.0) - T(845.0) * (K[6*i+3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * K[6*i] / T(27.0) + T(2.0) * K[6*i+1] - T(3544.0) * K[6*i+2] / T(2565.0) + T(1859.0) * K[6*i+3] / T(4104.0) - T(11.0) * (K[6*i+4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        K[6*i+5] = StepSize * _StateRate[i][0];

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(K[6*i] / T(360.0) - T(128.0) * K[6*i+2] / T(4275.0) - T(2197.0) * K[6*i+3] / T(75240.0) + K[6*i+4] / T(50.0) + T(2.0) * K[6*i+5] / T(55.0)) / StepSize;
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {
        t += StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + T(25.0) * K[6*i] / T(216.0) + T(1408.0) * K[6*i+2] / T(2565.0) + T(2197.0) * K[6*i+3] / T(4104.0) - K[6*i+4] / T(5.0);
        GetStateRate(t,GetStateRate_UserData);
        if(OutputState != NULL) OutputState(t,OutputState_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E >= T(4.0)) StepSize *= T(4.0); else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
      } else {
        t -= StepSize / T(2.0);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        GetStateRate(t,GetStateRate_UserData);
        E = T(0.84) * pow(Tolerance / Emax, T(0.25));
        if(E < T(0.1)) StepSize *= T(0.1); else StepSize *= E;
      };
    };
    };

    //Backward Runge-Kutta order 5 for the first 6 points
    StepSize = -StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      Y0[i] = _State[i][0];
    {
    T* k = new T[5*_StateCount];
    while(((StepSize > T(0.0)) && (t <= StartTime + T(5.1) * StepSize)) || ((StepSize < T(0.0)) && (t >= StartTime + T(5.1) * StepSize))) {

      for(unsigned int i=0;i<_StateCount;i++)
      {
        prevY[i] = _State[i][0];
        prevF[6*i+5] = prevF[6*i+4];
        prevF[6*i+4] = prevF[6*i+3];
        prevF[6*i+3] = prevF[6*i+2];
        prevF[6*i+2] = prevF[6*i+1];
        prevF[6*i+1] = prevF[6*i];
        _State[i][0] += (k[5*i] = StepSize * (prevF[6*i] = _StateRate[i][0])) / T(4.0);
      };

      t += StepSize / T(4.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

      t += StepSize / T(8.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

      t += T(57.0) * StepSize / T(104.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

      t += StepSize / T(13.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

      t -= StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
      for(unsigned int i=0;i<_StateCount;i++)
        _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

      t += StepSize / T(2.0);
      GetStateRate(t,GetStateRate_UserData);
    };
    };

    //Back to Forward iterations
    StepSize = -StepSize;
    t += T(6.0) * StepSize;
    for(unsigned int i=0;i<_StateCount;i++)
      _State[i][0] = Y0[i];
    GetStateRate(t,GetStateRate_UserData);
    while(((StepSize > T(0.0)) && (t < EndTime)) || ((StepSize < T(0.0)) && (t > EndTime))) {

      for(unsigned int i=0;i<_StateCount;i++) {
        prevY[i] = _State[i][0];
        prevF[6*i] = prevF[6*i+1];
        prevF[6*i+1] = prevF[6*i+2];
        prevF[6*i+2] = prevF[6*i+3];
        prevF[6*i+3] = prevF[6*i+4];
        prevF[6*i+4] = prevF[6*i+5];
        _State[i][0] += StepSize * (T(4277.0) * (prevF[6*i+5] = _StateRate[i][0]) - T(7923.0) * prevF[6*i+4] + T(9982.0) * prevF[6*i+3] - T(7298.0) * prevF[6*i+2] + T(2877.0) * prevF[6*i+1] - T(475.0) * prevF[6*i]) / T(1440.0);
        Y0[i] = _State[i][0];
      };

      t += StepSize;
      for(unsigned int j=0;j<Corrections;j++) {
        GetStateRate(t,GetStateRate_UserData);
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i] + StepSize * (T(475.0) * _StateRate[i][0] + T(1427.0) * prevF[6*i+5] - T(798.0) * prevF[6*i+4] + T(482.0) * prevF[6*i+3] - T(173.0) * prevF[6*i+2] + T(27.0) * prevF[6*i+1]) / T(1440.0);
      };

      Emax = 0.0;
      for(unsigned int i=0;i<_StateCount;i++) {
        E = fabs(T(863.0) * (_State[i][0] - Y0[i])/(T(19950.0) * StepSize));
        if(E > Emax) Emax = E;
      };

      if(Emax <= Tolerance) {   //Result Accepted
        GetStateRate(t,GetStateRate_UserData);
        if (OutputState != NULL) OutputState(t,OutputState_UserData);
        if(Emax <= Tolerance / T(4.0)) {
          E = pow(Tolerance / (T(2.0) * Emax),T(1.0/6.0));
          if(E >= T(4.0)) StepSize *= T(4.0);
          else StepSize *= E;
          if(fabs(StepSize) > MaxStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;
          if(((StepSize > T(0.0)) && (t - T(6.0) * StepSize < StartTime)) || ((StepSize < T(0.0)) && (t - T(6.0) * StepSize > StartTime))) StepSize = (t - StartTime) / T(6.0);
          if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;

          //Backward Runge-Kutta Order 5 iterations for 6 points
          StepSize = -StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            Y0[i] = _State[i][0];
          {
          T* k = new T[5*_StateCount];
          for(unsigned int j=0;j<6;j++) {

            for(unsigned int i=0;i<_StateCount;i++)
            {
              prevY[i] = _State[i][0];
              prevF[6*i+5] = prevF[6*i+4];
              prevF[6*i+4] = prevF[6*i+3];
              prevF[6*i+3] = prevF[6*i+2];
              prevF[6*i+2] = prevF[6*i+1];
              prevF[6*i+1] = prevF[6*i];
              _State[i][0] += (k[5*i] = StepSize * (prevF[6*i] = _StateRate[i][0])) / T(4.0);
            };

            t += StepSize / T(4.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

            t += StepSize / T(8.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

            t += T(57.0) * StepSize / T(104.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

            t += StepSize / T(13.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

            t -= StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
            for(unsigned int i=0;i<_StateCount;i++)
              _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

            t += StepSize / T(2.0);
            GetStateRate(t,GetStateRate_UserData);
          };
          };
          StepSize = -StepSize;
          t += T(6.0) * StepSize;
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = Y0[i];
          GetStateRate(t,GetStateRate_UserData);
        };
      } else {                  //Result Rejected
        t -= StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = prevY[i];
        E = pow(Tolerance / (T(2.0) * Emax),T(1.0/6.0));
        if(E <= T(0.1)) StepSize *= T(0.1);
        else StepSize *= E;
        if(((StepSize > T(0.0)) && (t + StepSize > EndTime)) || ((StepSize < T(0.0)) && (t + StepSize < EndTime))) StepSize = EndTime - t;
        if(fabs(StepSize) < MinStepSize) StepSize = fabs(MaxStepSize / StepSize) * StepSize;

        //Backward Runge-Kutta Order 5 iterations for 6 points
        StepSize = -StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          Y0[i] = _State[i][0];
        {
        T* k = new T[5*_StateCount];
        for(unsigned int j=0;j<6;j++) {

          for(unsigned int i=0;i<_StateCount;i++)
          {
            prevY[i] = _State[i][0];
            prevF[6*i+5] = prevF[6*i+4];
            prevF[6*i+4] = prevF[6*i+3];
            prevF[6*i+3] = prevF[6*i+2];
            prevF[6*i+2] = prevF[6*i+1];
            prevF[6*i+1] = prevF[6*i];
            _State[i][0] += (k[5*i] = StepSize * (prevF[6*i] = _StateRate[i][0])) / T(4.0);
          };

          t += StepSize / T(4.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (T(9.0) * (k[5*i + 1] = StepSize * _StateRate[i][0]) - T(5.0) * k[5*i]) / T(32.0);

          t += StepSize / T(8.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] += (T(276165.0) * k[5*i] - T(1250865.0) * k[5*i + 1] + T(1167360.0) * (k[5*i + 2] = StepSize * _StateRate[i][0])) / T(351520.0);

          t += T(57.0) * StepSize / T(104.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + T(439.0) * k[5*i] / T(216.0) - T(8.0) * k[5*i + 1] + T(3680.0) * k[5*i + 2] / T(513.0) - T(845.0) * (k[5*i + 3] = StepSize * _StateRate[i][0]) / T(4104.0);

          t += StepSize / T(13.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] - T(8.0) * k[5*i] / T(27.0) + T(2.0) * k[5*i + 1] - T(3544.0) * k[5*i + 2] / T(2565.0) + T(1859.0) * k[5*i + 3] / T(4104.0) - T(11.0) * (k[5*i + 4] = StepSize * _StateRate[i][0]) / T(40.0);

          t -= StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
          for(unsigned int i=0;i<_StateCount;i++)
            _State[i][0] = prevY[i] + T(16.0) * k[5*i] / T(135.0) + T(6656.0) * k[5*i + 2] / T(12825.0) + T(28561.0) * k[5*i + 3] / T(56430.0) - T(9.0) * k[5*i + 4] / T(50.0) + T(2.0) * StepSize * _StateRate[i][0] / T(55.0);

          t += StepSize / T(2.0);
          GetStateRate(t,GetStateRate_UserData);
        };
        };
        StepSize = -StepSize;
        t += T(6.0) * StepSize;
        for(unsigned int i=0;i<_StateCount;i++)
          _State[i][0] = Y0[i];
        GetStateRate(t,GetStateRate_UserData);
      };
    };
  };
  return 1;
};

//-----------CBDFVar------------------------------------------------------------

template <class T>
int CBDFVar<T>::Integrate() {
  if ((GetStateRate == NULL) || (_StateCount == 0) || (StepSize <= T(0.0)) || (Tolerance <= T(0.0)) || (MinStepSize > MaxStepSize) || (NewtonTolerance <= T(0.0))) return -1;

  if(Order == 0) Order = 1;
  else if(Order > 6) Order = 6;



  return 1;
};*/




};


