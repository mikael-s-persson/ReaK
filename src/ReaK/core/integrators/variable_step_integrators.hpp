/**
 * \file variable_step_integrators.hpp
 *
 * The following library implements numerical methods for integration of systems
 * of ordinary differential equations using variable time steps. The implementations
 * are done as described in the following books:\n\n
 *
 * Burden R.L. and Faires J.D., "Numerical Analysis", 8th Edition, Thomson, 2005.\n\n
 *
 * Ascher U.M. and Petzold L.R., "Computer Methods for Ordinary Differential Equations
 * Differential-Algebraic Equations, Society for Industrial and Applied Mathematics, 1998.\n\n
 *
 * The methods implemented are:\n\n
 *
 *   - Runge-Kutta-Fehlberg order 4-5 (CFEHLBERG45)\n
 *   - Dormand-Prince order 4-5 (CDORMANDPRINCE45)\n
 *   - Adams-Bashforth-Moulton Variable Step (up to order 6) (CADAMSBMVAR)\n
 *   - BDF Variable Step (up to order 6) (Backward Difference Formula with Variable Coefficient Strategy) (CBDFVAR)\n
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date july 2010
 */

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

#ifndef REAK_VARIABLE_STEP_INTEGRATORS_HPP
#define REAK_VARIABLE_STEP_INTEGRATORS_HPP


#include "integrator.hpp"

#include <cmath>

namespace ReaK {

  
/**
 * This class template implements at Runge-Kutta-Fehlberg integrator of order 4-5. This is a variable-step, 
 * explicit integrator of order 4 (with order 5 error estimation). Each integration step entails six evaluations 
 * of the state derivative. Error control is performed and can throw the ReaK::untolerable_integration exception 
 * if the integrator cannot acheive the required tolerance without lowering the time-step below the acceptable minimum.
 * Also basic verification of the integration parameters is done and might throw the ReaK::impossible_integration 
 * exception.
 */
template <class T>
class fehlberg45_integrator : public variable_step_integrator<T> {
  protected:

  public:

    virtual void RK_CALL integrate(double aEndTime);

    /**
     * Default constructor.
     */
    fehlberg45_integrator(const std::string& aName = "") : variable_step_integrator<T>(aName) { };
    
    /**
     * Parametrized constructor.
     * \param aName The name of this integrator object.
     * \param aState The initial state vector that the integrator will work with.
     * \param aStartTime The initial time to which the integrator is set.
     * \param aInitialStepSize The time-step used in the integration to start with (will be variable according to error control).
     * \param aGetStateRate A weak pointer to the object that will compute the state derivatives (see ReaK::state_rate_function).
     * \param aMaxStepSize The maximum time-step to be used during the integration, if error control allows it.
     * \param aMinStepSize The minimum time-step to be reached before declaring the integration untolerable due to error control.
     * \param aTolerance The desired relative error of the integrated state values.
     */
    fehlberg45_integrator(const std::string& aName,
                          const ReaK::vect_n<T>& aState,
                          double aStartTime,
                          double aInitialStepSize,
                          const weak_ptr< state_rate_function<T> >& aGetStateRate,
                          double aMaxStepSize,
                          double aMinStepSize,
                          double aTolerance) :
                          variable_step_integrator<T>(aName,aState,aStartTime,aInitialStepSize,aGetStateRate,aMaxStepSize,aMinStepSize,aTolerance) { };
    /**
     * Default destructor.
     */
    virtual ~fehlberg45_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      variable_step_integrator<T>::save(A,variable_step_integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      variable_step_integrator<T>::load(A,variable_step_integrator<T>::getStaticObjectType()->TypeVersion());
    };
        
    typedef fehlberg45_integrator<T> self;
    typedef variable_step_integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2220001,1,"fehlberg45_integrator",base)
};



//----------CFEHLBERG45-----------------------------------------------------

  /*
       0     |
       1/4   |  1/4
       3/8   |  3/32         9/32
       12/13 |  1932/2197   -7200/2197    7296/2197
       1     |  439/216     -8            3680/513    -845/4104
       1/2   | -8/27         2           -3544/2565    1859/4104   -11/40
       _____________________________________________________________________________________
       O(h^5)|  25/216       0            1408/2565    2197/4104   -1/5          0
       O(h^6)|  16/135       0            6656/12825   28561/56430 -9/50         2/55
  */

template <class T>
void RK_CALL fehlberg45_integrator<T>::integrate(double aEndTime) {
  using std::fabs;
  using std::pow;
  
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (aEndTime > integrator<T>::mTime)) ||
      (variable_step_integrator<T>::mTolerance <= 0.0) ||
      (variable_step_integrator<T>::mMinStepSize > variable_step_integrator<T>::mMaxStepSize))
    throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
  if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  vect_n<T> prevY(integrator<T>::mState.q.size());
  vect_n<T> k1(integrator<T>::mState.q.size());
  vect_n<T> k2(integrator<T>::mState.q.size());
  vect_n<T> k3(integrator<T>::mState.q.size());
  vect_n<T> k4(integrator<T>::mState.q.size());
  vect_n<T> k5(integrator<T>::mState.q.size());
  vect_n<T> k6(integrator<T>::mState.q.size());
  double R, Rmax;
  int worst_DOF;

  func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);

  while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    prevY = integrator<T>::mState;
    integrator<T>::mState += (k1 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(0.25);

    integrator<T>::mTime += integrator<T>::mStepSize * 0.25;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + (k1 * T(3.0) + (k2 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(9.0)) / T(32.0);


    integrator<T>::mTime += integrator<T>::mStepSize * 0.125;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + (k1 * T(1932.0) - k2 * T(7200.0) + (k3 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(7296.0)) / T(2197.0);

    integrator<T>::mTime += 57.0 * integrator<T>::mStepSize / 104.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + k1 * T(439.0 / 216.0) - k2 * T(8.0) + k3 * T(3680.0 / 513.0) - (k4 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(845.0 / 4104.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 13.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY - k1 * T(8.0 / 27.0) + k2 * T(2.0) - k3 * T(3544.0 / 2565.0) + k4 * T(1859.0 / 4104.0) - (k5 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(11.0 / 40.0);

    integrator<T>::mTime -= integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    k6 = integrator<T>::mStateRate * T(integrator<T>::mStepSize);

    Rmax = 0.0;
    worst_DOF = 0;
    for(unsigned int i=0;i<integrator<T>::mState.q.size();++i) {
      R = fabs((k1.q[i] / T(360.0) - T(128.0) * k3.q[i] / T(4275.0) - T(2197.0) * k4.q[i] / T(75240.0) + k5.q[i] / T(50.0) + T(2.0) * k6.q[i] / T(55.0)) / integrator<T>::mStepSize);
      if(R > Rmax) {
        Rmax = R;
        worst_DOF = i;
      };
    };

    if((Rmax > variable_step_integrator<T>::mTolerance) && (fabs(integrator<T>::mStepSize) > variable_step_integrator<T>::mMinStepSize)) {
      if(fabs(integrator<T>::mStepSize) <= variable_step_integrator<T>::mMinStepSize)
        throw untolerable_integration(variable_step_integrator<T>::mTolerance, Rmax, worst_DOF, integrator<T>::mStepSize, integrator<T>::mTime);

      integrator<T>::mTime -= integrator<T>::mStepSize * 0.5;
      integrator<T>::mState = prevY;
      func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
      R = 0.84 * pow(variable_step_integrator<T>::mTolerance / Rmax, 0.25);
      if(R < 0.1)
        integrator<T>::mStepSize *= 0.1;
      else
        integrator<T>::mStepSize *= R;
    } else {
      integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
      integrator<T>::mState = prevY + k1 * T(25.0 / 216.0) + k3 * T(1408.0 / 2565.0) + k4 * T(2197.0 / 4104.0) - k5 * T(0.2);
      func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);

      R = 0.84 * pow(variable_step_integrator<T>::mTolerance / Rmax, 0.25);
      if(R >= 4.0)
        integrator<T>::mStepSize *= 4.0;
      else if(R > 1.0)
        integrator<T>::mStepSize *= R;
      //if(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime + integrator<T>::mStepSize > aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime + integrator<T>::mStepSize < aEndTime)))
        //integrator<T>::mStepSize = aEndTime - integrator<T>::mTime;
    };

    if(fabs(integrator<T>::mStepSize) < variable_step_integrator<T>::mMinStepSize)
      integrator<T>::mStepSize *= fabs(variable_step_integrator<T>::mMinStepSize / integrator<T>::mStepSize);
    if(fabs(integrator<T>::mStepSize) > variable_step_integrator<T>::mMaxStepSize)
      integrator<T>::mStepSize *= fabs(variable_step_integrator<T>::mMaxStepSize / integrator<T>::mStepSize);
  };

};











/**
 * This class template implements at Dormand-Prince integrator of order 4-5. This is a variable-step, 
 * explicit integrator of order 4 (with order 5 error estimation). This integrator is very similar to the 
 * Runge-Kutta-Fehlberg integrator, but has a different set of intermediate points, with, in theory, better
 * numerical stability. Each integration step entails seven evaluations of the state derivative. Error control 
 * is performed and can throw the ReaK::untolerable_integration exception if the integrator cannot acheive the 
 * required tolerance without lowering the time-step below the acceptable minimum. Also basic verification of 
 * the integration parameters is done and might throw the ReaK::impossible_integration exception.
 */
template <class T>
class dormand_prince45_integrator : public variable_step_integrator<T> {
  protected:

  public:

    virtual void RK_CALL integrate(double aEndTime);

    /**
     * Default constructor.
     */
    dormand_prince45_integrator(const std::string& aName = "") : variable_step_integrator<T>(aName) { };
    
    /**
     * Parametrized constructor.
     * \param aName The name of this integrator object.
     * \param aState The initial state vector that the integrator will work with.
     * \param aStartTime The initial time to which the integrator is set.
     * \param aInitialStepSize The time-step used in the integration to start with (will be variable according to error control).
     * \param aGetStateRate A weak pointer to the object that will compute the state derivatives (see ReaK::state_rate_function).
     * \param aMaxStepSize The maximum time-step to be used during the integration, if error control allows it.
     * \param aMinStepSize The minimum time-step to be reached before declaring the integration untolerable due to error control.
     * \param aTolerance The desired relative error of the integrated state values.
     */
    dormand_prince45_integrator(const std::string& aName,
                                const ReaK::vect_n<T>& aState,
                                double aStartTime,
                                double aInitialStepSize,
                                const weak_ptr< state_rate_function<T> > aGetStateRate,
                                double aMaxStepSize,
                                double aMinStepSize,
                                double aTolerance) :
                                variable_step_integrator<T>(aName,aState,aStartTime,aInitialStepSize,aGetStateRate,aMaxStepSize,aMinStepSize,aTolerance) { };
    /**
     * Default destructor.
     */
    virtual ~dormand_prince45_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      variable_step_integrator<T>::save(A,variable_step_integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      variable_step_integrator<T>::load(A,variable_step_integrator<T>::getStaticObjectType()->TypeVersion());
    };
    
    typedef dormand_prince45_integrator<T> self;
    typedef variable_step_integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2220002,1,"dormand_prince45_integrator",base)
};



//----------CDORMANDPRINCE45-----------------------------------------------------

  /*
       0     |
       1/5   |  1/5
       3/10  |  3/40           9/40
       4/5   |  44/45         -56/15          32/9
       8/9   |  19372/6561    -25360/2187     64448/6561    -212/729
       1     |  9017/3168     -355/33        -46732/5247     49/176        -5103/18656
       1     |  35/384         0              500/1113       125/192       -2187/6784      11/84
       ______________________________________________________________________________________________
       O(h^6)|  5170/57600     0              7571/16695     393/640       -92097/339200   187/2100       1/40
       O(h^5)|  35/384         0              500/1113       125/192       -2187/6784      11/84          0
  */

template <class T>
void RK_CALL dormand_prince45_integrator<T>::integrate(double aEndTime) {
  using std::fabs;
  using std::pow;
  
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (aEndTime > integrator<T>::mTime)) ||
      (variable_step_integrator<T>::mTolerance <= 0.0) ||
      (variable_step_integrator<T>::mMinStepSize > variable_step_integrator<T>::mMaxStepSize))
    throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
  if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);
  vect_n<T> prevY(integrator<T>::mState.q.size());
  vect_n<T> k1(integrator<T>::mState.q.size());
  vect_n<T> k2(integrator<T>::mState.q.size());
  vect_n<T> k3(integrator<T>::mState.q.size());
  vect_n<T> k4(integrator<T>::mState.q.size());
  vect_n<T> k5(integrator<T>::mState.q.size());
  vect_n<T> k6(integrator<T>::mState.q.size());
  vect_n<T> k7(integrator<T>::mState.q.size());
  double R, Rmax;
  int worst_DOF;

  func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);

  while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    prevY = integrator<T>::mState;
    integrator<T>::mState += (k1 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) / T(5.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 5.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + (k1 * T(3.0) + (k2 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(9.0)) / T(40.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 10.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + k1 * T(44.0 / 45.0) - k2 * T(56.0 / 15.0) + (k3 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(32.0 / 9.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 2.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + k1 * T(19372.0 / 6561.0) - k2 * T(25360.0 / 2187.0) + k3 * T(64448.0 / 6561.0) - (k4 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(212.0 / 729.0);

    integrator<T>::mTime += 4.0 * integrator<T>::mStepSize / 45.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + k1 * T(9017.0 / 3168.0) - k2 * T(355.0 / 33.0) - k3 * T(46732.0 / 5247.0) + k4 * T(49.0 / 176.0) - (k5 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(5103.0 / 18656.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 9.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = prevY + k1 * T(35.0 / 384.0) + k3 * T(500.0 / 1113.0) + k4 * T(125.0 / 192.0) - k5 * T(2187.0 / 6784.0) + (k6 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(11.0 / 84.0);

    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    k7 = integrator<T>::mStateRate * T(integrator<T>::mStepSize);

    Rmax = 0.0;
    worst_DOF = 0;
    for(unsigned int i=0;i<integrator<T>::mState.q.size();++i) {
      R = fabs((k1.q[i] * T(5170.0 / 57600.0) + k3.q[i] * T(7571.0 / 16695.0) + k4.q[i] * T(393.0 / 640.0) - k5.q[i] * T(92097.0 / 339200.0) + k6.q[i] * T(187.0 / 2100.0) - k7.q[i] * T(39.0 / 40.0)) / T(integrator<T>::mStepSize));
      if(R > Rmax) {
        Rmax = R;
        worst_DOF = i;
      };
    };

    if((Rmax > variable_step_integrator<T>::mTolerance) &&(fabs(integrator<T>::mStepSize) > variable_step_integrator<T>::mMinStepSize)) {
      if(fabs(integrator<T>::mStepSize) <= variable_step_integrator<T>::mMinStepSize)
        throw untolerable_integration(variable_step_integrator<T>::mTolerance, Rmax, worst_DOF, integrator<T>::mStepSize, integrator<T>::mTime);

      integrator<T>::mTime -= integrator<T>::mStepSize;
      integrator<T>::mState = prevY;
      func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
      R = 0.84 * pow(variable_step_integrator<T>::mTolerance / Rmax, 0.25);
      if(R < 0.1)
        integrator<T>::mStepSize *= 0.1;
      else
        integrator<T>::mStepSize *= R;
    } else {
      integrator<T>::mState = prevY + k1 * T(5170.0 / 57600.0) + k3 * T(7571.0 / 16695.0) + k4 * T(393.0 / 640.0) - k5 * T(92097.0 / 339200.0) + k6 * T(187.0 / 2100.0) + k7 / T(40.0);
      func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);

      R = 0.84 * pow(variable_step_integrator<T>::mTolerance / Rmax, 0.25);
      if(R >= 4.0)
        integrator<T>::mStepSize *= 4.0;
      else if(R > 1.0)
        integrator<T>::mStepSize *= R;
      //if(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime + integrator<T>::mStepSize > aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime + integrator<T>::mStepSize < aEndTime)))
        //integrator<T>::mStepSize = aEndTime - integrator<T>::mTime;
    };

    if(fabs(integrator<T>::mStepSize) < variable_step_integrator<T>::mMinStepSize)
      integrator<T>::mStepSize *= fabs(variable_step_integrator<T>::mMinStepSize / integrator<T>::mStepSize);
    if(fabs(integrator<T>::mStepSize) > variable_step_integrator<T>::mMaxStepSize)
      integrator<T>::mStepSize *= fabs(variable_step_integrator<T>::mMaxStepSize / integrator<T>::mStepSize);
  };
};









/*
template <class T>
class CAdamsBMVar : public IVarStepIntegrator<T> {
  protected:

  public:
    virtual unsigned int GetChunkID() { return T::GetChunkID() | CHUNKID_CADAMSBMVAR; };

    virtual int Integrate();

    unsigned int Corrections;
    unsigned int Order;

    virtual unsigned int size() { return IVarStepIntegrator<T>::size() + 2*sizeof(unsigned int); };
    virtual bool WriteBuffer(void** Buf, unsigned int* Count) {
      RETURNNOTIFYWARNING(!IVarStepIntegrator<T>::WriteBuffer(Buf,Count),"Cannot Write CADAMSBMVAR:IVARSTEPINTEGRATOR Buffer!",false);

      RETURNNOTIFYWARNING(!((GlobalWriteBuffer(Corrections,Buf,Count)) && (GlobalWriteBuffer(Order,Buf,Count))),"Cannot Write CADAMSBMVAR Buffer!",false);

      return true;
    };
    virtual bool ReadBuffer(void** Buf, unsigned int* Count) {
      RETURNNOTIFYWARNING(!IVarStepIntegrator<T>::ReadBuffer(Buf,Count),"Cannot Read CADAMSBMVAR:IVARSTEPINTEGRATOR Buffer!",false);

      RETURNNOTIFYWARNING(!((GlobalReadBuffer(&Corrections,Buf,Count)) && (GlobalReadBuffer(&Order,Buf,Count))),"Cannot Read CADAMSBMVAR Buffer!",false);

      return true;
    };

    __fastcall CAdamsBMVar() : IVarStepIntegrator<T>() { Corrections = 1; Order = 3; return; };
    virtual __fastcall ~CAdamsBMVar() {
      if (RKRoot != NULL) RKRoot->MsgQueue->Push(new IMessage(MSGID_WARNING | MSGID_W_NOTIFY,"$DEL: Adams-BM Variable-Step Integrator Deleted",1,this,RKRoot->MCB_ObjectDeleted));
      return;
    };
};

template <class T>
class CBDFVar : public IVarStepIntegrator<T> {
  protected:

  public:
    virtual unsigned int GetChunkID() { return T::GetChunkID() | CHUNKID_CBDFVAR; };

    virtual int Integrate();

    unsigned int Order;
    T NewtonTolerance;

    virtual unsigned int size() { return IVarStepIntegrator<T>::size() + sizeof(unsigned int) + T::size(); };
    virtual bool WriteBuffer(void** Buf, unsigned int* Count) {
      RETURNNOTIFYWARNING(!IVarStepIntegrator<T>::WriteBuffer(Buf,Count),"Cannot Write CBDFVAR:IVARSTEPINTEGRATOR Buffer!",false);

      RETURNNOTIFYWARNING(!((GlobalWriteBuffer(&Order,Buf,Count)) && (NewtonTolerance.WriteBuffer(Buf,Count))),"Cannot Write CBDFVAR Buffer!",false);

      return true;
    };
    virtual bool ReadBuffer(void** Buf, unsigned int* Count) {
      RETURNNOTIFYWARNING(!IVarStepIntegrator<T>::ReadBuffer(Buf,Count),"Cannot Read CBDFVAR:IVARSTEPINTEGRATOR Buffer!",false);

      RETURNNOTIFYWARNING(!((GlobalReadBuffer(&Order,Buf,Count)) && (NewtonTolerance.ReadBuffer(Buf,Count))),"Cannot Read CBDFVAR Buffer!",false);

      return true;
    };

    __fastcall CBDFVar() : IVarStepIntegrator<T>() { NewtonTolerance = T(0.0001); Order = 3; return; };
    virtual __fastcall ~CBDFVar() {
      if (RKRoot != NULL) RKRoot->MsgQueue->Push(new IMessage(MSGID_WARNING | MSGID_W_NOTIFY,"$DEL: BDF Variable-Step Integrator Deleted",1,this,RKRoot->MCB_ObjectDeleted));
      return;
    };
};*/

};

#endif

//------------------------------------------------------------------------------
