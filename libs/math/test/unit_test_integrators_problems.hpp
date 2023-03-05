
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

#ifndef RK_UNIT_TEST_INTEGRATORS_PROBLEMS_HPP
#define RK_UNIT_TEST_INTEGRATORS_PROBLEMS_HPP

#include "ReaK/core/base/defs.hpp"

#include "integrator.hpp"

namespace ReaK {

template <typename T>
class iv_problem : public state_rate_function<T> {
 public:
  virtual ~iv_problem(){};

  virtual double getInitialTime() const = 0;
  virtual vect_n<T> getInitialValue() const = 0;
  virtual double getFinalTime() const = 0;
  virtual vect_n<T> getFinalValue() const = 0;

  typedef iv_problem<T> self;
  typedef state_rate_function<T> base_type;
  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC22FFFF0, 1, "iv_problem", base_type)
};

template <class T>
class HIRES_iv_problem : public iv_problem<T> {
 public:
  virtual ~HIRES_iv_problem(){};

  virtual void computeStateRate(double, const ReaK::vect_n<T>& aState,
                                ReaK::vect_n<T>& aStateRate) {
    aStateRate[0] =
        -1.71 * aState[0] + 0.43 * aState[1] + 8.32 * aState[2] + 0.0007;
    aStateRate[1] = 1.71 * aState[0] - 8.75 * aState[1];
    aStateRate[2] = -10.03 * aState[2] + 0.43 * aState[3] + 0.035 * aState[4];
    aStateRate[3] = 8.32 * aState[1] + 1.71 * aState[2] - 1.12 * aState[3];
    aStateRate[4] = -1.745 * aState[4] + 0.43 * aState[5] + 0.43 * aState[6];
    aStateRate[5] = -280.0 * aState[5] * aState[7] + 0.69 * aState[3] +
                    1.71 * aState[4] - 0.43 * aState[5] + 0.69 * aState[6];
    aStateRate[6] = 280.0 * aState[5] * aState[7] - 1.81 * aState[6];
    aStateRate[7] = -280.0 * aState[5] * aState[7] + 1.81 * aState[6];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(8);
    result[0] = 1.0;
    result[1] = 0.0;
    result[2] = 0.0;
    result[3] = 0.0;
    result[4] = 0.0;
    result[5] = 0.0;
    result[6] = 0.0;
    result[7] = 0.0057;
    return result;
  };

  virtual double getFinalTime() const { return 321.8122; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(8);
    result[0] = 0.7371312573325668e-3;
    result[1] = 0.1442485726316185e-3;
    result[2] = 0.5888729740967575e-4;
    result[3] = 0.1175651343283149e-2;
    result[4] = 0.2386356198831331e-2;
    result[5] = 0.6238968252742796e-2;
    result[6] = 0.2849998395185769e-2;
    result[7] = 0.2850001604814231e-2;
    return result;
  };

  typedef HIRES_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF1, 1, "HIRES_iv_problem",
                              base_type)
};

template <class T>
class Pollution_iv_problem : public iv_problem<T> {
 public:
  virtual ~Pollution_iv_problem(){};

  virtual void computeStateRate(double, const ReaK::vect_n<T>& aState,
                                ReaK::vect_n<T>& aStateRate) {
    const T k1 = 0.35;
    const T k2 = 0.266e2;
    const T k3 = 0.123e5;
    const T k4 = 0.86e-3;
    const T k5 = 0.82e-3;
    const T k6 = 0.15e5;
    const T k7 = 0.13e-3;
    const T k8 = 0.24e5;
    const T k9 = 0.165e5;
    const T k10 = 0.9e4;
    const T k11 = 0.22e-1;
    const T k12 = 0.12e5;
    const T k13 = 0.188e1;
    const T k14 = 0.163e5;
    const T k15 = 0.48e7;
    const T k16 = 0.35e-3;
    const T k17 = 0.175e-1;
    const T k18 = 0.1e9;
    const T k19 = 0.444e12;
    const T k20 = 0.124e4;
    const T k21 = 0.21e1;
    const T k22 = 0.578e1;
    const T k23 = 0.474e-1;
    const T k24 = 0.178e4;
    const T k25 = 0.312e1;

    vect_n<T> r(25);
    r[0] = k1 * aState[0];
    r[1] = k2 * aState[1] * aState[3];
    r[2] = k3 * aState[4] * aState[1];
    r[3] = k4 * aState[6];
    r[4] = k5 * aState[6];
    r[5] = k6 * aState[6] * aState[5];
    r[6] = k7 * aState[8];
    r[7] = k8 * aState[8] * aState[5];
    r[8] = k9 * aState[10] * aState[1];
    r[9] = k10 * aState[10] * aState[0];
    r[10] = k11 * aState[12];
    r[11] = k12 * aState[9] * aState[1];
    r[12] = k13 * aState[13];
    r[13] = k14 * aState[0] * aState[5];
    r[14] = k15 * aState[2];
    r[15] = k16 * aState[3];
    r[16] = k17 * aState[3];
    r[17] = k18 * aState[15];
    r[18] = k19 * aState[15];
    r[19] = k20 * aState[16] * aState[5];
    r[20] = k21 * aState[18];
    r[21] = k22 * aState[18];
    r[22] = k23 * aState[0] * aState[3];
    r[23] = k24 * aState[18] * aState[0];
    r[24] = k25 * aState[19];

    aStateRate[0] = -r[0] - r[9] - r[13] - r[22] - r[23] + r[1] + r[2] + r[8] +
                    r[10] + r[11] + r[21] + r[24];
    aStateRate[1] = -r[1] - r[2] - r[8] - r[11] + r[0] + r[20];
    aStateRate[2] = -r[14] + r[0] + r[16] + r[18] + r[21];
    aStateRate[3] = -r[1] - r[15] - r[16] - r[22] + r[14];
    aStateRate[4] = -r[2] + r[3] + r[3] + r[5] + r[6] + r[12] + r[19];
    aStateRate[5] = -r[5] - r[7] - r[13] - r[19] + r[2] + r[17] + r[17];
    aStateRate[6] = -r[3] - r[4] - r[5] + r[12];
    aStateRate[7] = r[3] + r[4] + r[5] + r[6];
    aStateRate[8] = -r[6] - r[7];
    aStateRate[9] = -r[11] + r[6] + r[8];
    aStateRate[10] = -r[8] - r[9] + r[7] + r[10];
    aStateRate[11] = r[8];
    aStateRate[12] = -r[10] + r[9];
    aStateRate[13] = -r[12] + r[11];
    aStateRate[14] = r[13];
    aStateRate[15] = -r[17] - r[18] + r[15];
    aStateRate[16] = -r[19];
    aStateRate[17] = r[19];
    aStateRate[18] = -r[20] - r[21] - r[23] + r[22] + r[24];
    aStateRate[19] = -r[24] + r[23];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(20);
    result[0] = 0.0;
    result[1] = 0.2;
    result[2] = 0.0;
    result[3] = 0.04;
    result[4] = 0.0;
    result[5] = 0.0;
    result[6] = 0.1;
    result[7] = 0.3;
    result[8] = 0.01;
    result[9] = 0.0;
    result[10] = 0.0;
    result[11] = 0.0;
    result[12] = 0.0;
    result[13] = 0.0;
    result[14] = 0.0;
    result[15] = 0.0;
    result[16] = 0.007;
    result[17] = 0.0;
    result[18] = 0.0;
    result[19] = 0.0;
    return result;
  };

  virtual double getFinalTime() const { return 60.0; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(20);
    result[0] = 0.5646255480022769e-01;
    result[1] = 0.1342484130422339;
    result[2] = 0.4139734331099427e-08;
    result[3] = 0.5523140207484359e-02;
    result[4] = 0.2018977262302196e-06;
    result[5] = 0.1464541863493966e-06;
    result[6] = 0.7784249118997964e-01;
    result[7] = 0.3245075353396018;
    result[8] = 0.7494013383880406e-02;
    result[9] = 0.1622293157301561e-07;
    result[10] = 0.1135863833257075e-07;
    result[11] = 0.2230505975721359e-02;
    result[12] = 0.2087162882798630e-03;
    result[13] = 0.1396921016840158e-04;
    result[14] = 0.8964884856898295e-02;
    result[15] = 0.4352846369330103e-17;
    result[16] = 0.6899219696263405e-02;
    result[17] = 0.1007803037365946e-03;
    result[18] = 0.1772146513969984e-05;
    result[19] = 0.5682943292316392e-04;
    return result;
  };

  typedef Pollution_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF2, 1, "Pollution_iv_problem",
                              base_type)
};

template <class T>
class RingModulator_iv_problem : public iv_problem<T> {
 public:
  virtual ~RingModulator_iv_problem(){};

  virtual void computeStateRate(double aTime, const ReaK::vect_n<T>& aState,
                                ReaK::vect_n<T>& aStateRate) {
    using std::exp;
    using std::sin;

    const T c = 1.6e-8;
    const T cs = 2.0e-12;
    const T cp = 1.0e-8;
    const T r = 25.0e3;
    const T rp = 50.0;
    const T lh = 4.45;
    const T ls1 = 2.0e-3;
    const T ls2 = 5.0e-4;
    const T ls3 = 5.0e-4;
    const T rg1 = 36.3;
    const T rg2 = 17.3;
    const T rg3 = 17.3;
    const T ri = 5.0e1;
    const T rc = 6.0e2;
    const T gamma = 40.67286402e-9;
    const T delta = 17.7493332;
    const T pi = 3.141592653589793238462643383;

    T uin1 = 0.5 * sin(2.0e3 * pi * aTime);
    T uin2 = 2.0 * sin(2.0e4 * pi * aTime);
    T ud1 = aState[2] - aState[4] - aState[6] - uin2;
    T ud2 = -aState[3] + aState[5] - aState[6] - uin2;
    T ud3 = aState[3] + aState[4] + aState[6] + uin2;
    T ud4 = -aState[2] - aState[5] + aState[6] + uin2;

    T qud1 = gamma * (exp(delta * ud1) - 1.0);
    T qud2 = gamma * (exp(delta * ud2) - 1.0);
    T qud3 = gamma * (exp(delta * ud3) - 1.0);
    T qud4 = gamma * (exp(delta * ud4) - 1.0);

    aStateRate[0] = (aState[7] - 0.5 * aState[9] + 0.5 * aState[10] +
                     aState[13] - aState[0] / r) /
                    c;
    aStateRate[1] = (aState[8] - 0.5 * aState[11] + 0.5 * aState[12] +
                     aState[14] - aState[1] / r) /
                    c;
    aStateRate[2] = (aState[9] - qud1 + qud4) / cs;
    aStateRate[3] = (-aState[10] + qud2 - qud3) / cs;
    aStateRate[4] = (aState[11] + qud1 - qud3) / cs;
    aStateRate[5] = (-aState[12] - qud2 + qud4) / cs;
    aStateRate[6] = (-aState[6] / rp + qud1 + qud2 - qud3 - qud4) / cp;
    aStateRate[7] = -aState[0] / lh;
    aStateRate[8] = -aState[1] / lh;
    aStateRate[9] = (0.5 * aState[0] - aState[2] - rg2 * aState[9]) / ls2;
    aStateRate[10] = (-0.5 * aState[0] + aState[3] - rg3 * aState[10]) / ls3;
    aStateRate[11] = (0.5 * aState[1] - aState[4] - rg2 * aState[11]) / ls2;
    aStateRate[12] = (-0.5 * aState[1] + aState[5] - rg3 * aState[12]) / ls3;
    aStateRate[13] = (-aState[0] + uin1 - (ri + rg1) * aState[13]) / ls1;
    aStateRate[14] = (-aState[1] - (rc + rg1) * aState[14]) / ls1;
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const { return vect_n<T>(15, T(0.0)); };

  virtual double getFinalTime() const { return 1.0e-3; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(15);
    result[0] = -0.2339057358486745e-1;
    result[1] = -0.7367485485540825e-2;
    result[2] = 0.2582956709291169;
    result[3] = -0.4064465721283450;
    result[4] = -0.4039455665149794;
    result[5] = 0.2607966765422943;
    result[6] = 0.1106761861269975;
    result[7] = 0.2939904342435596e-6;
    result[8] = -0.2840029933642329e-7;
    result[9] = 0.7267198267264553e-3;
    result[10] = 0.7929487196960840e-3;
    result[11] = -0.7255283495698965e-3;
    result[12] = -0.7941401968526521e-3;
    result[13] = 0.7088495416976114e-4;
    result[14] = 0.2390059075236570e-4;
    return result;
  };

  typedef RingModulator_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF3, 1, "RingModulator_iv_problem",
                              base_type)
};

template <class T>
class AkzoNobel_iv_problem : public iv_problem<T> {
 public:
  virtual ~AkzoNobel_iv_problem(){};

  virtual void computeStateRate(double aTime, const ReaK::vect_n<T>& aState,
                                ReaK::vect_n<T>& aStateRate) {
    const T k = 100.0;
    const T c = 4.0;

    const std::size_t N = 200;
    const T dzeta = 1.0 / T(N);
    const T dzeta2 = dzeta * dzeta;
    T dum = (dzeta - 1.0) * (dzeta - 1.0) / c;
    T alpha = 2.0 * (dzeta - 1.0) * dum / c;
    T beta = dum * dum;

    T phi = 0.0;
    if (aTime < 5.0)
      phi = 2.0;

    aStateRate[0] = (phi - 2.0 * aState[0] + aState[2]) * beta / dzeta2 +
                    alpha * (aState[2] - phi) / (2.0 * dzeta) -
                    k * aState[0] * aState[1];
    aStateRate[1] = -k * aState[0] * aState[1];

    for (std::size_t j = 1; j < N - 1; ++j) {
      T zeta = j * dzeta;
      dum = (zeta - 1.0) * (zeta - 1.0) / c;
      alpha = 2.0 * (zeta - 1.0) * dum / c;
      beta = dum * dum;
      std::size_t i = 2 * j;
      aStateRate[i] =
          (aState[i - 2] - 2.0 * aState[i] + aState[i + 2]) * beta / dzeta2 +
          alpha * (aState[i + 2] - aState[i - 2]) / (2.0 * dzeta) -
          k * aState[i] * aState[i + 1];
      ++i;
      aStateRate[i] = -k * aState[i] * aState[i - 1];
    };

    aStateRate[2 * N - 2] = -k * aState[2 * N - 2] * aState[2 * N - 1];
    aStateRate[2 * N - 1] = -k * aState[2 * N - 2] * aState[2 * N - 1];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(400);

    for (std::size_t i = 0; i < 200; ++i) {
      result[2 * i] = 0.0;
      result[2 * i + 1] = 1.0;
    };

    return result;
  };

  virtual double getFinalTime() const { return 20.0; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(400);
    result[0] = 0.5113983840919909e-005;
    result[1] = 0.1925112884312553e-143;
    result[2] = 0.1027858770570419e-004;
    result[3] = 0.1890518289312031e-142;
    result[4] = 0.1549349862635799e-004;
    result[5] = 0.1774199325357386e-142;
    result[6] = 0.2075835344757462e-004;
    result[7] = 0.5897341137981092e-143;
    result[8] = 0.2607273610116854e-004;
    result[9] = 0.1093527900908030e-143;
    result[10] = 0.3143617475695002e-004;
    result[11] = 0.1188834841626416e-144;
    result[12] = 0.3684813884509626e-004;
    result[13] = 0.9968323236025642e-147;
    result[14] = 0.4230803594492533e-004;
    result[15] = -0.2801994001528093e-146;
    result[16] = 0.4781520853483223e-004;
    result[17] = -0.7337417669341249e-147;
    result[18] = 0.5336893059800053e-004;
    result[19] = -0.1209033101530330e-147;
    result[20] = 0.5896840407836044e-004;
    result[21] = -0.1430357497530360e-148;
    result[22] = 0.6461275518112516e-004;
    result[23] = -0.1063952641824646e-149;
    result[24] = 0.7030103051210320e-004;
    result[25] = 0.7939969136126717e-152;
    result[26] = 0.7603219304985662e-004;
    result[27] = 0.1568246940545520e-150;
    result[28] = 0.8180511794465543e-004;
    result[29] = 0.4074950357924872e-150;
    result[30] = 0.8761858813806752e-004;
    result[31] = 0.5592746648679992e-150;
    result[32] = 0.9347128979692480e-004;
    result[33] = -0.5510388943414421e-151;
    result[34] = 0.9936180755532036e-004;
    result[35] = -0.2724738349250769e-149;
    result[36] = 0.1052886195582220e-003;
    result[37] = -0.9327772452398718e-149;
    result[38] = 0.1112500923002360e-003;
    result[39] = -0.2182885200987554e-148;
    result[40] = 0.1172444752530255e-003;
    result[41] = -0.4041450806475518e-148;
    result[42] = 0.1232698952748828e-003;
    result[43] = -0.5608157478395261e-148;
    result[44] = 0.1293243507959787e-003;
    result[45] = -0.2639662630908699e-148;
    result[46] = 0.1354057057728661e-003;
    result[47] = 0.1801866277537073e-147;
    result[48] = 0.1415116834059119e-003;
    result[49] = 0.8464449882759417e-147;
    result[50] = 0.1476398596134615e-003;
    result[51] = 0.2245234937355967e-146;
    result[52] = 0.1537876562567258e-003;
    result[53] = 0.3359213489153582e-146;
    result[54] = 0.1599523341096154e-003;
    result[55] = -0.3085721171916412e-146;
    result[56] = 0.1661309855680449e-003;
    result[57] = -0.4465322607423735e-145;
    result[58] = 0.1723205270935920e-003;
    result[59] = -0.1970925996866384e-144;
    result[60] = 0.1785176913868402e-003;
    result[61] = -0.6070953121563027e-144;
    result[62] = 0.1847190192862588e-003;
    result[63] = -0.1412011918930335e-143;
    result[64] = 0.1909208513890961e-003;
    result[65] = -0.2378861987352203e-143;
    result[66] = 0.1971193193914910e-003;
    result[67] = -0.2380432473186974e-143;
    result[68] = 0.2033103371458565e-003;
    result[69] = -0.6522557638254663e-145;
    result[70] = 0.2094895914345677e-003;
    result[71] = 0.1784305601809064e-143;
    result[72] = 0.2156525324601176e-003;
    result[73] = -0.1007474781780816e-142;
    result[74] = 0.2217943640531935e-003;
    result[75] = -0.5281511349479423e-142;
    result[76] = 0.2279100336016016e-003;
    result[77] = -0.1117525482975987e-141;
    result[78] = 0.2339942217046434e-003;
    result[79] = -0.1127916494884468e-141;
    result[80] = 0.2400413315594459e-003;
    result[81] = -0.1633306916231411e-142;
    result[82] = 0.2460454780878912e-003;
    result[83] = 0.2708874035585891e-143;
    result[84] = 0.2520004768152150e-003;
    result[85] = -0.2501941069702609e-142;
    result[86] = 0.2578998325140575e-003;
    result[87] = -0.2642308070750020e-141;
    result[88] = 0.2637367276308081e-003;
    result[89] = -0.3684887530751217e-139;
    result[90] = 0.2695040105145025e-003;
    result[91] = -0.3647274179805887e-138;
    result[92] = 0.2751941834723564e-003;
    result[93] = -0.1255641406397419e-137;
    result[94] = 0.2807993906802854e-003;
    result[95] = -0.1694257216823904e-138;
    result[96] = 0.2863114059815211e-003;
    result[97] = -0.1785516142939602e-136;
    result[98] = 0.2917216206117258e-003;
    result[99] = -0.3935939757647002e-135;
    result[100] = 0.2970210308948898e-003;
    result[101] = -0.2514765666933440e-134;
    result[102] = 0.3022002259608294e-003;
    result[103] = -0.7200873856605984e-134;
    result[104] = 0.3072493755423352e-003;
    result[105] = -0.7539683247227422e-134;
    result[106] = 0.3121582179180383e-003;
    result[107] = 0.3738577086039426e-135;
    result[108] = 0.3169160480759169e-003;
    result[109] = -0.2493582962172335e-131;
    result[110] = 0.3215117061821543e-003;
    result[111] = 0.3039632438293726e-130;
    result[112] = 0.3259335664508512e-003;
    result[113] = 0.5321044068586611e-128;
    result[114] = 0.3301695265219917e-003;
    result[115] = -0.1918129324351378e-126;
    result[116] = 0.3342069974681551e-003;
    result[117] = -0.1336929159252586e-124;
    result[118] = 0.3380328945648600e-003;
    result[119] = 0.9521748754010357e-123;
    result[120] = 0.3416336289752354e-003;
    result[121] = 0.1001197393324181e-120;
    result[122] = 0.3449951005170561e-003;
    result[123] = 0.2703860993866771e-119;
    result[124] = 0.3481026916991771e-003;
    result[125] = 0.4365133580297076e-119;
    result[126] = 0.3509412632351946e-003;
    result[127] = 0.4898111237855383e-115;
    result[128] = 0.3534951512648823e-003;
    result[129] = 0.1621439381962246e-112;
    result[130] = 0.3557481665387581e-003;
    result[131] = 0.3003220203772183e-110;
    result[132] = 0.3576835958481664e-003;
    result[133] = 0.5931668289615909e-108;
    result[134] = 0.3592842060126915e-003;
    result[135] = 0.2235590472383775e-105;
    result[136] = 0.3605322507686931e-003;
    result[137] = 0.1025457293602057e-102;
    result[138] = 0.3614094809374544e-003;
    result[139] = 0.3496613568296336e-100;
    result[140] = 0.3618971582890092e-003;
    result[141] = 0.4767073568395508e-098;
    result[142] = 0.3619760735583436e-003;
    result[143] = -0.2410784286794997e-095;
    result[144] = 0.3616265691144918e-003;
    result[145] = -0.9188398110576038e-093;
    result[146] = 0.3608285668302233e-003;
    result[147] = 0.1146623087995081e-089;
    result[148] = 0.3595616017506735e-003;
    result[149] = 0.1649638439865233e-086;
    result[150] = 0.3578048622135169e-003;
    result[151] = 0.1215140240350217e-083;
    result[152] = 0.3555372371311931e-003;
    result[153] = 0.7134490346394154e-081;
    result[154] = 0.3527373712073181e-003;
    result[155] = 0.4502515392738464e-078;
    result[156] = 0.3493837289247301e-003;
    result[157] = 0.7138395988310312e-075;
    result[158] = 0.3454546682115489e-003;
    result[159] = 0.9941693919247076e-071;
    result[160] = 0.3409285247640208e-003;
    result[161] = 0.2012859826753015e-066;
    result[162] = 0.3357837080804970e-003;
    result[163] = 0.3598261520662423e-062;
    result[164] = 0.3299988103392750e-003;
    result[165] = 0.5466580008990664e-058;
    result[166] = 0.3235527293336597e-003;
    result[167] = 0.6945384844951550e-054;
    result[168] = 0.3164248067597393e-003;
    result[169] = 0.7275415527806026e-050;
    result[170] = 0.3085949832350532e-003;
    result[171] = 0.6193143746524996e-046;
    result[172] = 0.3000439715082906e-003;
    result[173] = 0.4219255556214135e-042;
    result[174] = 0.2907534493998412e-003;
    result[175] = 0.2263678154715720e-038;
    result[176] = 0.2807062740884081e-003;
    result[177] = 0.9401607967545219e-035;
    result[178] = 0.2698867194275612e-003;
    result[179] = 0.2968231730793053e-031;
    result[180] = 0.2582807380350103e-003;
    result[181] = 0.6987463944434805e-028;
    result[182] = 0.2458762499428408e-003;
    result[183] = 0.1201641789884051e-024;
    result[184] = 0.2326634596245027e-003;
    result[185] = 0.1477169946829840e-021;
    result[186] = 0.2186352032185982e-003;
    result[187] = 0.1268462422099779e-018;
    result[188] = 0.2037873277440060e-003;
    result[189] = 0.7425015664001834e-016;
    result[190] = 0.1881191040379240e-003;
    result[191] = 0.2886826929895103e-013;
    result[192] = 0.1716336750388461e-003;
    result[193] = 0.7252477041900172e-011;
    result[194] = 0.1543385408702044e-003;
    result[195] = 0.1143390654212691e-008;
    result[196] = 0.1362460820444338e-003;
    result[197] = 0.1096625145716966e-006;
    result[198] = 0.1173741304462833e-003;
    result[199] = 0.6190822732534586e-005;
    result[200] = 0.9774701310627047e-004;
    result[201] = 0.1986273404756002e-003;
    result[202] = 0.7740788649977313e-004;
    result[203] = 0.3489773624098464e-002;
    result[204] = 0.5657119003189305e-004;
    result[205] = 0.3234526094359604e-001;
    result[206] = 0.3643334879766658e-004;
    result[207] = 0.1548747348410801e+000;
    result[208] = 0.2003152841880950e-004;
    result[209] = 0.4026980529594953e+000;
    result[210] = 0.9608297851720770e-005;
    result[211] = 0.6649744834198490e+000;
    result[212] = 0.4215537698495267e-005;
    result[213] = 0.8409284546320647e+000;
    result[214] = 0.1753504402754791e-005;
    result[215] = 0.9314946676956936e+000;
    result[216] = 0.7048158429518009e-006;
    result[217] = 0.9720896201631835e+000;
    result[218] = 0.2760943506466737e-006;
    result[219] = 0.9890204872799944e+000;
    result[220] = 0.1057554501281432e-006;
    result[221] = 0.9957930123519514e+000;
    result[222] = 0.3965142250779033e-007;
    result[223] = 0.9984246531478463e+000;
    result[224] = 0.1455273204279008e-007;
    result[225] = 0.9994229325942358e+000;
    result[226] = 0.5226348147846279e-008;
    result[227] = 0.9997932125999319e+000;
    result[228] = 0.1835610545325733e-008;
    result[229] = 0.9999275409325039e+000;
    result[230] = 0.6301078589385454e-009;
    result[231] = 0.9999751869380269e+000;
    result[232] = 0.2112538351365564e-009;
    result[233] = 0.9999917015131560e+000;
    result[234] = 0.6912550453447044e-010;
    result[235] = 0.9999972914302640e+000;
    result[236] = 0.2205932132514696e-010;
    result[237] = 0.9999991378543379e+000;
    result[238] = 0.6860095639285670e-011;
    result[239] = 0.9999997325855174e+000;
    result[240] = 0.2077324462852526e-011;
    result[241] = 0.9999999192384585e+000;
    result[242] = 0.6120038908594393e-012;
    result[243] = 0.9999999762710279e+000;
    result[244] = 0.1752695518797070e-012;
    result[245] = 0.9999999932230490e+000;
    result[246] = 0.4875001992978682e-013;
    result[247] = 0.9999999981203191e+000;
    result[248] = 0.1315706848908981e-013;
    result[249] = 0.9999999994941428e+000;
    result[250] = 0.3442274192104633e-014;
    result[251] = 0.9999999998680372e+000;
    result[252] = 0.8721783456154470e-015;
    result[253] = 0.9999999999666630e+000;
    result[254] = 0.2137938962858872e-015;
    result[255] = 0.9999999999918528e+000;
    result[256] = 0.5064735930780995e-016;
    result[257] = 0.9999999999980759e+000;
    result[258] = 0.1158284928109727e-016;
    result[259] = 0.9999999999995613e+000;
    result[260] = 0.2554350586347124e-017;
    result[261] = 0.9999999999999036e+000;
    result[262] = 0.5425563935887811e-018;
    result[263] = 0.9999999999999796e+000;
    result[264] = 0.1108623976460997e-018;
    result[265] = 0.9999999999999958e+000;
    result[266] = 0.2176490922739810e-019;
    result[267] = 0.9999999999999992e+000;
    result[268] = 0.4100180074816888e-020;
    result[269] = 0.9999999999999998e+000;
    result[270] = 0.7401919443964595e-021;
    result[271] = 0.1000000000000000e+001;
    result[272] = 0.1278745657114596e-021;
    result[273] = 0.1000000000000000e+001;
    result[274] = 0.2111087049605767e-022;
    result[275] = 0.1000000000000000e+001;
    result[276] = 0.3325632734364699e-023;
    result[277] = 0.1000000000000000e+001;
    result[278] = 0.4991515592566292e-024;
    result[279] = 0.1000000000000000e+001;
    result[280] = 0.7126950428617158e-025;
    result[281] = 0.1000000000000000e+001;
    result[282] = 0.9664740804131475e-026;
    result[283] = 0.1000000000000000e+001;
    result[284] = 0.1242716896959521e-026;
    result[285] = 0.1000000000000000e+001;
    result[286] = 0.1512543532243458e-027;
    result[287] = 0.1000000000000000e+001;
    result[288] = 0.1739533019752215e-028;
    result[289] = 0.1000000000000000e+001;
    result[290] = 0.1886942537979667e-029;
    result[291] = 0.1000000000000000e+001;
    result[292] = 0.1926965705022792e-030;
    result[293] = 0.1000000000000000e+001;
    result[294] = 0.1849021812823421e-031;
    result[295] = 0.1000000000000000e+001;
    result[296] = 0.1663798767415642e-032;
    result[297] = 0.1000000000000000e+001;
    result[298] = 0.1401076830818626e-033;
    result[299] = 0.1000000000000000e+001;
    result[300] = 0.1101818149402153e-034;
    result[301] = 0.1000000000000000e+001;
    result[302] = 0.8074224739509168e-036;
    result[303] = 0.1000000000000000e+001;
    result[304] = 0.5501249196662931e-037;
    result[305] = 0.1000000000000000e+001;
    result[306] = 0.3476859813132770e-038;
    result[307] = 0.1000000000000000e+001;
    result[308] = 0.2033489290876775e-039;
    result[309] = 0.1000000000000000e+001;
    result[310] = 0.1097880013869247e-040;
    result[311] = 0.1000000000000000e+001;
    result[312] = 0.5457825200381417e-042;
    result[313] = 0.1000000000000000e+001;
    result[314] = 0.2491675366427318e-043;
    result[315] = 0.1000000000000000e+001;
    result[316] = 0.1041801880291617e-044;
    result[317] = 0.1000000000000000e+001;
    result[318] = 0.3978066491064419e-046;
    result[319] = 0.1000000000000000e+001;
    result[320] = 0.1383174699098532e-047;
    result[321] = 0.1000000000000000e+001;
    result[322] = 0.4365911791079500e-049;
    result[323] = 0.1000000000000000e+001;
    result[324] = 0.1247057764661705e-050;
    result[325] = 0.1000000000000000e+001;
    result[326] = 0.3212728839963712e-052;
    result[327] = 0.1000000000000000e+001;
    result[328] = 0.7439366703571565e-054;
    result[329] = 0.1000000000000000e+001;
    result[330] = 0.1542770387822259e-055;
    result[331] = 0.1000000000000000e+001;
    result[332] = 0.2854454245592573e-057;
    result[333] = 0.1000000000000000e+001;
    result[334] = 0.4693220411250150e-059;
    result[335] = 0.1000000000000000e+001;
    result[336] = 0.6828458274546624e-061;
    result[337] = 0.1000000000000000e+001;
    result[338] = 0.8752952529541412e-063;
    result[339] = 0.1000000000000000e+001;
    result[340] = 0.9838541433761416e-065;
    result[341] = 0.1000000000000000e+001;
    result[342] = 0.9649177728609193e-067;
    result[343] = 0.1000000000000000e+001;
    result[344] = 0.8213596936190817e-069;
    result[345] = 0.1000000000000000e+001;
    result[346] = 0.6033986647865674e-071;
    result[347] = 0.1000000000000000e+001;
    result[348] = 0.3802531117966294e-073;
    result[349] = 0.1000000000000000e+001;
    result[350] = 0.2042261117698575e-075;
    result[351] = 0.1000000000000000e+001;
    result[352] = 0.9282595096128614e-078;
    result[353] = 0.1000000000000000e+001;
    result[354] = 0.3543587864454877e-080;
    result[355] = 0.1000000000000000e+001;
    result[356] = 0.1126779423370979e-082;
    result[357] = 0.1000000000000000e+001;
    result[358] = 0.2957534367766753e-085;
    result[359] = 0.1000000000000000e+001;
    result[360] = 0.6344600529877694e-088;
    result[361] = 0.1000000000000000e+001;
    result[362] = 0.1100279075462365e-090;
    result[363] = 0.1000000000000000e+001;
    result[364] = 0.1523845293461783e-093;
    result[365] = 0.1000000000000000e+001;
    result[366] = 0.1662696161555950e-096;
    result[367] = 0.1000000000000000e+001;
    result[368] = 0.1407578290673998e-099;
    result[369] = 0.1000000000000000e+001;
    result[370] = 0.9086150803567186e-103;
    result[371] = 0.1000000000000000e+001;
    result[372] = 0.4384339596163745e-106;
    result[373] = 0.1000000000000000e+001;
    result[374] = 0.1545482064392824e-109;
    result[375] = 0.1000000000000000e+001;
    result[376] = 0.3874172613928345e-113;
    result[377] = 0.1000000000000000e+001;
    result[378] = 0.6689452219441953e-117;
    result[379] = 0.1000000000000000e+001;
    result[380] = 0.7655680935317283e-121;
    result[381] = 0.1000000000000000e+001;
    result[382] = 0.5538543899545850e-125;
    result[383] = 0.1000000000000000e+001;
    result[384] = 0.2386173886563501e-129;
    result[385] = 0.1000000000000000e+001;
    result[386] = 0.5664887497790931e-134;
    result[387] = 0.1000000000000000e+001;
    result[388] = 0.6671124967149171e-139;
    result[389] = 0.1000000000000000e+001;
    result[390] = 0.3351973480286951e-144;
    result[391] = 0.1000000000000000e+001;
    result[392] = 0.5684315818559200e-150;
    result[393] = 0.1000000000000000e+001;
    result[394] = 0.2142121793294590e-156;
    result[395] = 0.1000000000000000e+001;
    result[396] = 0.6727117900187205e-164;
    result[397] = 0.1000000000000000e+001;
    result[398] = 0.0000000000000000e+000;
    result[399] = 0.1000000000000000e+001;
    return result;
  };

  typedef AkzoNobel_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF4, 1, "AkzoNobel_iv_problem",
                              base_type)
};

template <class T>
class Pleiades_iv_problem : public iv_problem<T> {
 public:
  virtual ~Pleiades_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {
    using std::pow;

    for (std::size_t i = 0; i < 7; ++i) {
      T sumx = 0.0;
      T sumy = 0.0;
      for (std::size_t j = 0; j < 7; ++j) {
        T rij =
            (aState[i] - aState[j]) * (aState[i] - aState[j]) +
            (aState[i + 7] - aState[j + 7]) * (aState[i + 7] - aState[j + 7]);
        T rij32 = pow(rij, 1.5);
        if (j != i) {
          sumx += j * (aState[j] - aState[i]) / rij32;
          sumy += j * (aState[j + 7] - aState[i + 7]) / rij32;
        };
      };
      aStateRate[i + 14] = sumx;
      aStateRate[i + 21] = sumy;
    };

    for (std::size_t i = 0; i < 14; ++i)
      aStateRate[i] = aState[i + 14];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(28);

    result[0] = 3.0;
    result[1] = 3.0;
    result[2] = -1.0;
    result[3] = -3.0;
    result[4] = 2.0;
    result[5] = -2.0;
    result[6] = 2.0;
    result[7] = 3.0;
    result[8] = -3.0;
    result[9] = 2.0;
    result[10] = 0.0;
    result[11] = 0.0;
    result[12] = -4.0;
    result[13] = 4.0;
    result[14] = 0.0;
    result[15] = 0.0;
    result[16] = 0.0;
    result[17] = 0.0;
    result[18] = 0.0;
    result[19] = 1.75;
    result[20] = -1.5;
    result[21] = 0.0;
    result[22] = 0.0;
    result[23] = 0.0;
    result[24] = -1.25;
    result[25] = 1.0;
    result[26] = 0.0;
    result[27] = 0.0;

    return result;
  };

  virtual double getFinalTime() const { return 3.0; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(28);
    result[0] = 0.3706139143970502;
    result[1] = 0.3237284092057233e1;
    result[2] = -0.3222559032418324e1;
    result[3] = 0.6597091455775310;
    result[4] = 0.3425581707156584;
    result[5] = 0.1562172101400631e1;
    result[6] = -0.7003092922212495;
    result[7] = -0.3943437585517392e1;
    result[8] = -0.3271380973972550e1;
    result[9] = 0.5225081843456543e1;
    result[10] = -0.2590612434977470e1;
    result[11] = 0.1198213693392275e1;
    result[12] = -0.2429682344935824;
    result[13] = 0.1091449240428980e1;
    result[14] = 0.3417003806314313e1;
    result[15] = 0.1354584501625501e1;
    result[16] = -0.2590065597810775e1;
    result[17] = 0.2025053734714242e1;
    result[18] = -0.1155815100160448e1;
    result[19] = -0.8072988170223021;
    result[20] = 0.5952396354208710;
    result[21] = -0.3741244961234010e1;
    result[22] = 0.3773459685750630;
    result[23] = 0.9386858869551073;
    result[24] = 0.3667922227200571;
    result[25] = -0.3474046353808490;
    result[26] = 0.2344915448180937e1;
    result[27] = -0.1947020434263292e1;
    return result;
  };

  typedef Pleiades_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF5, 1, "Pleiades_iv_problem",
                              base_type)
};

template <class T>
class VanDerPol_iv_problem : public iv_problem<T> {
 public:
  virtual ~VanDerPol_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {

    aStateRate[0] = aState[1];
    aStateRate[1] =
        ((1.0 - aState[0] * aState[0]) * aState[1] - aState[0]) / 1.0e-6;
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(2);

    result[0] = 2.0;
    result[1] = 0.0;

    return result;
  };

  virtual double getFinalTime() const { return 2.0; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(2);

    result[0] = 0.1706167732170483e1;
    result[1] = -0.8928097010247975;

    return result;
  };

  typedef VanDerPol_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF6, 1, "VanDerPol_iv_problem",
                              base_type)
};

template <class T>
class VanDerPolMod_iv_problem : public iv_problem<T> {
 public:
  virtual ~VanDerPolMod_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {

    aStateRate[0] = aState[1];
    aStateRate[1] =
        1.0e3 * (1.0 - aState[0] * aState[0]) * aState[1] - aState[0];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(2);

    result[0] = 2.0;
    result[1] = 0.0;

    return result;
  };

  virtual double getFinalTime() const { return 2.0e3; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(2);

    result[0] = 0.1706167732170469e1;
    result[1] = -0.8928097010248125e-3;

    return result;
  };

  typedef VanDerPolMod_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF7, 1, "VanDerPolMod_iv_problem",
                              base_type)
};

template <class T>
class Orego_iv_problem : public iv_problem<T> {
 public:
  virtual ~Orego_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {

    aStateRate[0] =
        77.27 *
        (aState[1] + aState[0] * (1.0 - 8.375e-6 * aState[0] - aState[1]));
    aStateRate[1] = (aState[2] - (1.0 + aState[0]) * aState[1]) / 77.27;
    aStateRate[2] = 0.161 * (aState[0] - aState[2]);
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(3);

    result[0] = 1.0;
    result[1] = 2.0;
    result[2] = 3.0;

    return result;
  };

  virtual double getFinalTime() const { return 360.0; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(3);

    result[0] = 0.1000814870318523e1;
    result[1] = 0.1228178521549917e4;
    result[2] = 0.1320554942846706e3;

    return result;
  };

  typedef Orego_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF8, 1, "Orego_iv_problem",
                              base_type)
};

template <class T>
class Rober_iv_problem : public iv_problem<T> {
 public:
  virtual ~Rober_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {

    aStateRate[0] = -0.04 * aState[0] + 1.0e4 * aState[1] * aState[2];
    aStateRate[1] = 0.04 * aState[0] - 1.0e4 * aState[1] * aState[2] -
                    3.0e7 * aState[1] * aState[1];
    aStateRate[2] = 3.0e7 * aState[1] * aState[1];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(3);

    result[0] = 1.0;
    result[1] = 0.0;
    result[2] = 0.0;

    return result;
  };

  virtual double getFinalTime() const { return 1.0e11; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(3);

    result[0] = 0.2083340149701255e-7;
    result[1] = 0.8333360770334713e-13;
    result[2] = 0.9999999791665050;

    return result;
  };

  typedef Rober_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFF9, 1, "Rober_iv_problem",
                              base_type)
};

template <class T>
class E5_iv_problem : public iv_problem<T> {
 public:
  virtual ~E5_iv_problem(){};

  virtual void computeStateRate(double aTime, const vect_n<T>& aState,
                                vect_n<T>& aStateRate) {

    T prod1 = 7.89e-10 * aState[0];
    T prod2 = 1.1e7 * aState[0] * aState[2];
    T prod3 = 1.13e9 * aState[1] * aState[2];
    T prod4 = 1.13e3 * aState[3];
    aStateRate[0] = -prod1 - prod2;
    aStateRate[1] = prod1 - prod3;
    aStateRate[3] = prod2 - prod4;
    aStateRate[2] = aStateRate[1] - aStateRate[3];
  };

  virtual double getInitialTime() const { return 0.0; };

  virtual vect_n<T> getInitialValue() const {
    vect_n<T> result(4);

    result[0] = 1.76e-3;
    result[1] = 0.0;
    result[2] = 0.0;
    result[3] = 0.0;

    return result;
  };

  virtual double getFinalTime() const { return 1.0e13; };

  virtual vect_n<T> getFinalValue() const {
    vect_n<T> result(4);

    result[0] = 0.1152903278711829e-290;
    result[1] = 0.8867655517642120e-22;
    result[2] = 0.8854814626268838e-22;
    result[3] = 0.0;

    return result;
  };

  typedef E5_iv_problem<T> self;
  typedef iv_problem<T> base_type;
  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC22FFFFA, 1, "E5_iv_problem", base_type)
};
};  // namespace ReaK

#endif
