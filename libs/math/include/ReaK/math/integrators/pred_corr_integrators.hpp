/**
 * \file pred_corr_integrators.hpp
 *
 * The following library implements numerical methods for integration of systems
 * of ordinary differential equations using predictor-corrector strategies. The implementations
 * are done as described in the following books:\n\n
 *
 * Burden R.L. and Faires J.D., "Numerical Analysis", 8th Edition, Thomson, 2005.\n\n
 *
 * Ascher U.M. and Petzold L.R., "Computer Methods for Ordinary Differential Equations
 * Differential-Algebraic Equations", Society for Industrial and Applied Mathematics, 1998.\n\n
 *
 * The methods implemented here:\n\n
 *
 *   - Adams-Bashforth-Moulton order 3 (adamsBM3_integrator)\n
 *   - Adams-Bashforth-Moulton order 5 (adamsBM5_integrator)\n
 *   - Modified Hamming Method (hammind_mod_integrator)\n
 *   - Iterated Modified Hamming Method (hamming_iter_mod_integrator)\n
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

#ifndef REAK_PRED_CORR_INTEGRATORS_HPP
#define REAK_PRED_CORR_INTEGRATORS_HPP

#include "ReaK/math/integrators/integrator.hpp"

#include <cmath>

namespace ReaK {

/**
 * This class template implements at Adams-Bashford-Moulton integrator of order 3. This is a fixed-step, multi-step,
 * predictor-corrector integrator of order 3. Each integration step entails (1 + CorrectionCount) evaluations
 * of the state derivative. Neither error control nor divergence detection is performed.
 * Only basic verification of the integration parameters is done and might throw the ReaK::impossible_integration
 * exception. The multiple steps are initialized using a Runge-Kutta method of the same order.
 */
template <class T>
class adamsBM3_integrator : public integrator<T> {
 protected:
  unsigned int mCorrectionCount;
  vect_n<T> prevY;
  vect_n<T> prevF1;
  vect_n<T> prevF2;
  vect_n<T> prevF3;

  void initializePrevVectors();

 public:
  void setStepSize(double aNewStepSize) override {
    integrator<T>::setStepSize(aNewStepSize);
    this->initializePrevVectors();
  };
  void setTime(double aNewTime) override {
    integrator<T>::setTime(aNewTime);
    this->initializePrevVectors();
  }

  void clearStateVector() override {
    integrator<T>::clearStateVector();
    this->initializePrevVectors();
  }
  void addStateElement(T Element) override {
    integrator<T>::addStateElement(Element);
    this->initializePrevVectors();
  }
  void addStateElements(const ReaK::vect_n<T>& Elements) override {
    integrator<T>::addStateElements(Elements);
    this->initializePrevVectors();
  }

  void setStateRateFunc(
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate) override {
    integrator<T>::setStateRateFunc(aGetStateRate);
    this->initializePrevVectors();
  }

  /**
   * This function allows you to set the number of corrections to perform (in the predictor-corrector scheme).
   * \param aCorrectionCount The new correction count to use.
   */
  void setCorrectionCount(unsigned int aCorrectionCount) {
    mCorrectionCount = aCorrectionCount;
  }

  void integrate(double aEndTime) override;

  /**
   * Default constructor.
   */
  adamsBM3_integrator(const std::string& aName = "")
      : integrator<T>(aName),
        mCorrectionCount(3),
        prevY(),
        prevF1(),
        prevF2(),
        prevF3() {}
  /**
   * Parametrized constructor.
   * \param aName The name of this integrator object.
   * \param aState The initial state vector that the integrator will work with.
   * \param aStartTime The initial time to which the integrator is set.
   * \param aStepSize The time-step used in the integration.
   * \param aCorrectionCount The number of corrections to perform.
   */
  adamsBM3_integrator(
      const std::string& aName, const ReaK::vect_n<T>& aState,
      double aStartTime, double aStepSize,
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate,
      unsigned int aCorrectionCount = 3)
      : integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate),
        mCorrectionCount(aCorrectionCount) {
    this->initializePrevVectors();
  }
  /**
   * Default destructor.
   */
  ~adamsBM3_integrator() override = default;

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    integrator<T>::save(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mCorrectionCount) &
        RK_SERIAL_SAVE_WITH_NAME(prevY) & RK_SERIAL_SAVE_WITH_NAME(prevF1) &
        RK_SERIAL_SAVE_WITH_NAME(prevF2) & RK_SERIAL_SAVE_WITH_NAME(prevF3);
  }
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    integrator<T>::load(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mCorrectionCount) &
        RK_SERIAL_LOAD_WITH_NAME(prevY) & RK_SERIAL_LOAD_WITH_NAME(prevF1) &
        RK_SERIAL_LOAD_WITH_NAME(prevF2) & RK_SERIAL_LOAD_WITH_NAME(prevF3);
  }

  using self = adamsBM3_integrator<T>;
  using base = integrator<T>;

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2210005, 1, "adamsBM3_integrator", base)
};

//----------adamsBM_integrator------------------------------------------------------

/*
  Bashforth Coefficients
  p | k | j->    |  1     2     3     4     5     6      |  C
  -------------------------------------------------------------------
  1 | 1 | Bj     |  1                                    | 1/2
  2 | 2 | 2Bj    |  3    -1                              | 5/12
  3 | 3 | 12Bj   |  23   -16    5                        | 3/8
  4 | 4 | 24Bj   |  55   -59    37   -9                  | 251/720
  5 | 5 | 720Bj  |  1901 -2774  2616 -1274  251          | 95/288
  6 | 6 | 1440Bj |  4277 -7923  9982 -7298  2877 -475    | 19087/60480

  Yn = Yn-1 + h*sum(j=1->k, Bj*Fn-j);

  Moulton Coefficients
  p | k | j->    |  0     1     2     3     4     5      |  C
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

*/

template <class T>
void adamsBM3_integrator<T>::initializePrevVectors() {
  this->prevY.q.resize(integrator<T>::mState.q.size());
  this->prevF1.q.resize(integrator<T>::mState.q.size());
  this->prevF2.q.resize(integrator<T>::mState.q.size());
  this->prevF3.q.resize(integrator<T>::mState.q.size());

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    return;
  }

  double back_time = integrator<T>::mTime;
  vect_n<T> prevY2(integrator<T>::mState.q.size());
  vect_n<T> prevY3(integrator<T>::mState.q.size());

  func_ptr->computeStateRate(back_time, integrator<T>::mState,
                             integrator<T>::mStateRate);

  // Order 1
  this->prevY = integrator<T>::mState -
                integrator<T>::mStateRate * T(integrator<T>::mStepSize);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, this->prevY, this->prevF1);
    this->prevY =
        integrator<T>::mState - this->prevF1 * T(integrator<T>::mStepSize);
  }

  func_ptr->computeStateRate(back_time, this->prevY, this->prevF1);

  // Order 2
  prevY2 = this->prevY - (this->prevF1 * T(3.0) - integrator<T>::mStateRate) *
                             T(integrator<T>::mStepSize * 0.5);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY2, this->prevF2);
    prevY2 = this->prevY -
             (this->prevF1 + this->prevF2) * T(integrator<T>::mStepSize * 0.5);
  }

  func_ptr->computeStateRate(back_time, prevY2, this->prevF2);

  // Order 3
  prevY3 = prevY2 - (this->prevF2 * T(23.0) - this->prevF1 * T(16.0) +
                     integrator<T>::mStateRate * T(5.0)) *
                        T(integrator<T>::mStepSize / 12.0);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY3, this->prevF3);
    prevY3 = prevY2 -
             (this->prevF3 * T(5.0) + this->prevF2 * T(8.0) - this->prevF1) *
                 T(integrator<T>::mStepSize / 12.0);
  }

  func_ptr->computeStateRate(back_time, prevY3, this->prevF3);
}

template <class T>
void adamsBM3_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr)
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);

  if (this->mCorrectionCount == 0) {
    this->mCorrectionCount = 1;
  }

  func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                             integrator<T>::mStateRate);

  while (
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    this->prevY = integrator<T>::mState;
    this->prevF3 = this->prevF2;
    this->prevF2 = this->prevF1;
    integrator<T>::mState +=
        ((this->prevF1 = integrator<T>::mStateRate) * T(23) -
         this->prevF2 * T(16) + this->prevF3 * T(5.0)) *
        T(integrator<T>::mStepSize / 12.0);

    integrator<T>::mTime += integrator<T>::mStepSize;
    for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
      func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                                 integrator<T>::mStateRate);
      integrator<T>::mState =
          this->prevY + (integrator<T>::mStateRate * T(5.0) +
                         this->prevF1 * T(8.0) - this->prevF2) *
                            T(integrator<T>::mStepSize / 12.0);
    }

    func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                               integrator<T>::mStateRate);
  }
}

/**
 * This class template implements at Adams-Bashford-Moulton integrator of order 5. This is a fixed-step, multi-step,
 * predictor-corrector integrator of order 5. Each integration step entails (1 + CorrectionCount) evaluations
 * of the state derivative. Neither error control nor divergence detection is performed.
 * Only basic verification of the integration parameters is done and might throw the ReaK::impossible_integration
 * exception. The multiple steps are initialized using a Runge-Kutta method of the same order.
 */
template <class T>
class adamsBM5_integrator : public integrator<T> {
 protected:
  unsigned int mCorrectionCount;
  vect_n<T> prevY;
  vect_n<T> prevF1;
  vect_n<T> prevF2;
  vect_n<T> prevF3;
  vect_n<T> prevF4;
  vect_n<T> prevF5;

  void initializePrevVectors();

 public:
  void setStepSize(double aNewStepSize) override {
    integrator<T>::setStepSize(aNewStepSize);
    this->initializePrevVectors();
  }
  void setTime(double aNewTime) override {
    integrator<T>::setTime(aNewTime);
    this->initializePrevVectors();
  }

  void clearStateVector() override {
    integrator<T>::clearStateVector();
    this->initializePrevVectors();
  }
  void addStateElement(T Element) override {
    integrator<T>::addStateElement(Element);
    this->initializePrevVectors();
  }
  void addStateElements(const ReaK::vect_n<T>& Elements) override {
    integrator<T>::addStateElements(Elements);
    this->initializePrevVectors();
  }

  void setStateRateFunc(
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate) override {
    integrator<T>::setStateRateFunc(aGetStateRate);
    this->initializePrevVectors();
  }

  /**
   * This function allows you to set the number of corrections to perform (in the predictor-corrector scheme).
   * \param aCorrectionCount The new correction count to use.
   */
  void setCorrectionCount(unsigned int aCorrectionCount) {
    mCorrectionCount = aCorrectionCount;
  }

  void integrate(double aEndTime) override;

  /**
   * Default constructor.
   */
  adamsBM5_integrator(const std::string& aName = "")
      : integrator<T>(aName),
        mCorrectionCount(3),
        prevY(),
        prevF1(),
        prevF2(),
        prevF3(),
        prevF4(),
        prevF5() {}
  /**
   * Parametrized constructor.
   * \param aName The name of this integrator object.
   * \param aState The initial state vector that the integrator will work with.
   * \param aStartTime The initial time to which the integrator is set.
   * \param aStepSize The time-step used in the integration.
   * \param aCorrectionCount The number of corrections to perform.
   */
  adamsBM5_integrator(
      const std::string& aName, const ReaK::vect_n<T>& aState,
      double aStartTime, double aStepSize,
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate,
      unsigned int aCorrectionCount = 3)
      : integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate),
        mCorrectionCount(aCorrectionCount) {
    this->initializePrevVectors();
  }
  /**
   * Default destructor.
   */
  ~adamsBM5_integrator() override = default;

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    integrator<T>::save(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mCorrectionCount) &
        RK_SERIAL_SAVE_WITH_NAME(prevY) & RK_SERIAL_SAVE_WITH_NAME(prevF1) &
        RK_SERIAL_SAVE_WITH_NAME(prevF2) & RK_SERIAL_SAVE_WITH_NAME(prevF3) &
        RK_SERIAL_SAVE_WITH_NAME(prevF4) & RK_SERIAL_SAVE_WITH_NAME(prevF5);
  }
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    integrator<T>::load(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mCorrectionCount) &
        RK_SERIAL_LOAD_WITH_NAME(prevY) & RK_SERIAL_LOAD_WITH_NAME(prevF1) &
        RK_SERIAL_LOAD_WITH_NAME(prevF2) & RK_SERIAL_LOAD_WITH_NAME(prevF3) &
        RK_SERIAL_LOAD_WITH_NAME(prevF4) & RK_SERIAL_LOAD_WITH_NAME(prevF5);
  }

  using self = adamsBM5_integrator<T>;
  using base = integrator<T>;

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2210006, 1, "adamsBM5_integrator", base)
};

template <class T>
void adamsBM5_integrator<T>::initializePrevVectors() {
  this->prevY.q.resize(integrator<T>::mState.q.size());
  this->prevF1.q.resize(integrator<T>::mState.q.size());
  this->prevF2.q.resize(integrator<T>::mState.q.size());
  this->prevF3.q.resize(integrator<T>::mState.q.size());
  this->prevF4.q.resize(integrator<T>::mState.q.size());
  this->prevF5.q.resize(integrator<T>::mState.q.size());

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    return;
  }

  double back_time = integrator<T>::mTime;
  vect_n<T> prevY2(integrator<T>::mState.q.size());
  vect_n<T> prevY3(integrator<T>::mState.q.size());
  vect_n<T> prevY4(integrator<T>::mState.q.size());
  vect_n<T> prevY5(integrator<T>::mState.q.size());

  func_ptr->computeStateRate(back_time, integrator<T>::mState,
                             integrator<T>::mStateRate);

  // Order 1
  this->prevY = integrator<T>::mState -
                integrator<T>::mStateRate * T(integrator<T>::mStepSize);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, this->prevY, this->prevF1);
    this->prevY =
        integrator<T>::mState - this->prevF1 * T(integrator<T>::mStepSize);
  }

  func_ptr->computeStateRate(back_time, this->prevY, this->prevF1);

  // Order 2
  prevY2 = this->prevY - (this->prevF1 * T(3.0) - integrator<T>::mStateRate) *
                             T(integrator<T>::mStepSize * 0.5);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY2, this->prevF2);
    prevY2 = this->prevY -
             (this->prevF1 + this->prevF2) * T(integrator<T>::mStepSize * 0.5);
  }

  func_ptr->computeStateRate(back_time, prevY2, this->prevF2);

  // Order 3
  prevY3 = prevY2 - (this->prevF2 * T(23.0) - this->prevF1 * T(16.0) +
                     integrator<T>::mStateRate * T(5.0)) *
                        T(integrator<T>::mStepSize / 12.0);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY3, this->prevF3);
    prevY3 = prevY2 -
             (this->prevF3 * T(5.0) + this->prevF2 * T(8.0) - this->prevF1) *
                 T(integrator<T>::mStepSize / 12.0);
  }

  func_ptr->computeStateRate(back_time, prevY3, this->prevF3);

  // Order 4
  prevY4 =
      prevY3 - (this->prevF3 * T(55.0) - this->prevF2 * T(59.0) +
                this->prevF1 * T(37.0) - integrator<T>::mStateRate * T(9.0)) *
                   T(integrator<T>::mStepSize / 24.0);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY4, this->prevF4);
    prevY4 = prevY3 - (this->prevF4 * T(9.0) + this->prevF3 * T(19.0) -
                       this->prevF2 * T(5.0) + this->prevF1) *
                          T(integrator<T>::mStepSize / 24.0);
  }

  func_ptr->computeStateRate(back_time, prevY4, this->prevF4);

  // Order 5
  prevY5 = prevY4 - (this->prevF4 * T(1901.0) - this->prevF3 * T(2774.0) +
                     this->prevF2 * T(2616.0) - this->prevF1 * T(1274.0) +
                     integrator<T>::mStateRate * T(251.0)) *
                        T(integrator<T>::mStepSize / 720.0);

  back_time -= integrator<T>::mStepSize;
  for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
    func_ptr->computeStateRate(back_time, prevY5, this->prevF5);
    prevY5 = prevY4 - (this->prevF5 * T(251.0) + this->prevF4 * T(646.0) -
                       this->prevF3 * T(264.0) + this->prevF2 * T(106.0) -
                       this->prevF1 * T(19.0)) *
                          T(integrator<T>::mStepSize / 720.0);
  }
}

template <class T>
void adamsBM5_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  if (this->mCorrectionCount == 0) {
    this->mCorrectionCount = 1;
  }

  func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                             integrator<T>::mStateRate);

  while (
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    this->prevY = integrator<T>::mState;
    this->prevF5 = this->prevF4;
    this->prevF4 = this->prevF3;
    this->prevF3 = this->prevF2;
    this->prevF2 = this->prevF1;
    integrator<T>::mState +=
        ((this->prevF1 = integrator<T>::mStateRate) * T(1901.0) -
         this->prevF2 * T(2774.0) + this->prevF3 * T(2616.0) -
         this->prevF4 * T(1274.0) + this->prevF5 * T(251.0)) *
        T(integrator<T>::mStepSize / 720.0);

    integrator<T>::mTime += integrator<T>::mStepSize;
    for (unsigned int j = 0; j < this->mCorrectionCount; ++j) {
      func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                                 integrator<T>::mStateRate);
      integrator<T>::mState =
          this->prevY + (integrator<T>::mStateRate * T(251.0) +
                         this->prevF1 * T(646.0) - this->prevF2 * T(264.0) +
                         this->prevF3 * T(106.0) - this->prevF4 * T(19.0)) *
                            T(integrator<T>::mStepSize / 720.0);
    }

    func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                               integrator<T>::mStateRate);
  }
}

/**
 * This class template implements at Modified Hamming integrator (order 3). This is a fixed-step, multi-step,
 * predictor-corrector integrator of order 3. Each integration step entails two evaluations
 * of the state derivative. Neither error control nor divergence detection is performed.
 * Only basic verification of the integration parameters is done and might throw the ReaK::impossible_integration
 * exception.
 */
template <class T>
class hamming_mod_integrator : public integrator<T> {
 protected:
  vect_n<T> mY_n_3;
  vect_n<T> mY_n_2;
  vect_n<T> mY_n_1;
  vect_n<T> mYp_n_2;
  vect_n<T> mYp_n_1;
  vect_n<T> mP_n;
  vect_n<T> mP;
  vect_n<T> mC;
  vect_n<T> mM;
  vect_n<T> mMp;

  void initializePrevVectors();

 public:
  void setStepSize(double aNewStepSize) override {
    integrator<T>::setStepSize(aNewStepSize);
    this->initializePrevVectors();
  }
  void setTime(double aNewTime) override {
    integrator<T>::setTime(aNewTime);
    this->initializePrevVectors();
  }

  void clearStateVector() override {
    integrator<T>::clearStateVector();
    this->initializePrevVectors();
  }
  void addStateElement(T Element) override {
    integrator<T>::addStateElement(Element);
    this->initializePrevVectors();
  }
  void addStateElements(const ReaK::vect_n<T>& Elements) override {
    integrator<T>::addStateElements(Elements);
    this->initializePrevVectors();
  }

  void setStateRateFunc(
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate) override {
    integrator<T>::setStateRateFunc(aGetStateRate);
    this->initializePrevVectors();
  }

  void integrate(double aEndTime) override;

  /**
   * Default constructor.
   */
  hamming_mod_integrator(const std::string& aName = "")
      : integrator<T>(aName),
        mY_n_3(),
        mY_n_2(),
        mY_n_1(),
        mYp_n_2(),
        mYp_n_1(),
        mP_n(),
        mP(),
        mC(),
        mM(),
        mMp() {}
  /**
   * Parametrized constructor.
   * \param aName The name of this integrator object.
   * \param aState The initial state vector that the integrator will work with.
   * \param aStartTime The initial time to which the integrator is set.
   * \param aStepSize The time-step used in the integration.
   * \param aGetStateRate A weak pointer to the object that will compute the state derivatives (see
   * ReaK::state_rate_function).
   */
  hamming_mod_integrator(
      const std::string& aName, const ReaK::vect_n<T>& aState,
      double aStartTime, double aStepSize,
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate)
      : integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate) {
    this->initializePrevVectors();
  }
  /**
   * Default destructor.
   */
  ~hamming_mod_integrator() override = default;

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    integrator<T>::save(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mY_n_3) & RK_SERIAL_SAVE_WITH_NAME(mY_n_2) &
        RK_SERIAL_SAVE_WITH_NAME(mY_n_1) & RK_SERIAL_SAVE_WITH_NAME(mYp_n_2) &
        RK_SERIAL_SAVE_WITH_NAME(mYp_n_1) & RK_SERIAL_SAVE_WITH_NAME(mP_n) &
        RK_SERIAL_SAVE_WITH_NAME(mC);
  }
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    integrator<T>::load(A, integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mY_n_3) & RK_SERIAL_LOAD_WITH_NAME(mY_n_2) &
        RK_SERIAL_LOAD_WITH_NAME(mY_n_1) & RK_SERIAL_LOAD_WITH_NAME(mYp_n_2) &
        RK_SERIAL_LOAD_WITH_NAME(mYp_n_1) & RK_SERIAL_LOAD_WITH_NAME(mP_n) &
        RK_SERIAL_LOAD_WITH_NAME(mC);
  }

  using self = hamming_mod_integrator<T>;
  using base = integrator<T>;

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2210008, 1, "hammind_mod_integrator",
                              base)
};

/*
 P_n+1 = Y_n-3 + dt * 4/3 * (2*Yp_n - Yp_n-1 + 2*Yp_n-2)
 M = P_n+1 - 112/121 * (P_n - C)
 C = 1/8 * (9*Y_n - Y_n-2 + dt * 3 * (Mp + 2*Yp_n - Yp_n-1))
 Y_n+1 = C + 9/121 * (P_n+1 - C)
*/

template <class T>
void hamming_mod_integrator<T>::initializePrevVectors() {
  this->mY_n_1.q.resize(integrator<T>::mState.q.size());
  this->mY_n_2.q.resize(integrator<T>::mState.q.size());
  this->mY_n_3.q.resize(integrator<T>::mState.q.size());
  this->mC.q.resize(integrator<T>::mState.q.size());
  this->mM.q.resize(integrator<T>::mState.q.size());
  this->mMp.q.resize(integrator<T>::mState.q.size());
  this->mYp_n_1.q.resize(integrator<T>::mState.q.size());
  this->mYp_n_2.q.resize(integrator<T>::mState.q.size());
  this->mP_n.q.resize(integrator<T>::mState.q.size());
  this->mP.q.resize(integrator<T>::mState.q.size());

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    return;
  }

  double back_time = integrator<T>::mTime;
  vect_n<T> w(integrator<T>::mState.q.size());
  vect_n<T> k1(integrator<T>::mState.q.size());
  vect_n<T> k2(integrator<T>::mState.q.size());
  vect_n<T> k3(integrator<T>::mState.q.size());
  vect_n<T> k4(integrator<T>::mState.q.size());
  vect_n<T> k5(integrator<T>::mState.q.size());

  func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                             integrator<T>::mStateRate);

  w = integrator<T>::mState;
  this->mY_n_1 =
      w -
      (k1 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(0.25);

  back_time -= integrator<T>::mStepSize * 0.25;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);
  this->mY_n_1 -= ((k2 = this->mYp_n_1 * T(integrator<T>::mStepSize)) * T(9.0) -
                   k1 * T(5.0)) /
                  T(32.0);

  back_time -= integrator<T>::mStepSize * 0.125;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);
  this->mY_n_1 -=
      (k1 * T(276165.0) - k2 * T(1250865.0) +
       (k3 = this->mYp_n_1 * T(integrator<T>::mStepSize)) * T(1167360.0)) /
      T(351520.0);

  back_time -= 57.0 * integrator<T>::mStepSize / 104.0;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);
  this->mY_n_1 =
      w - k1 * T(439.0 / 216.0) + k2 * T(8.0) - k3 * T(3680.0 / 513.0) +
      (k4 = this->mYp_n_1 * T(integrator<T>::mStepSize)) * T(845.0 / 4104.0);

  back_time -= integrator<T>::mStepSize / 13.0;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);
  this->mY_n_1 =
      w + k1 * T(8.0 / 27.0) - k2 * T(2.0) + k3 * T(3544.0 / 2565.0) -
      k4 * T(1859.0 / 4104.0) +
      (k5 = this->mYp_n_1 * T(integrator<T>::mStepSize)) * T(11.0 / 40.0);

  back_time += integrator<T>::mStepSize * 0.5;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);
  this->mY_n_1 = w - k1 * T(16.0 / 135.0) - k3 * T(6656.0 / 12825.0) -
                 k4 * T(28561.0 / 56430.0) + k5 * T(9.0 / 50.0) -
                 this->mYp_n_1 * T(2.0 * integrator<T>::mStepSize / 55.0);

  back_time -= integrator<T>::mStepSize * 0.5;
  func_ptr->computeStateRate(back_time, this->mY_n_1, this->mYp_n_1);

  w = this->mY_n_1;
  this->mY_n_2 =
      w - (k1 = this->mYp_n_1 * T(integrator<T>::mStepSize)) * T(0.25);

  back_time -= integrator<T>::mStepSize * 0.25;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);
  this->mY_n_2 -= ((k2 = this->mYp_n_2 * T(integrator<T>::mStepSize)) * T(9.0) -
                   k1 * T(5.0)) /
                  T(32.0);

  back_time -= integrator<T>::mStepSize * 0.125;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);
  this->mY_n_2 -=
      (k1 * T(276165.0) - k2 * T(1250865.0) +
       (k3 = this->mYp_n_2 * T(integrator<T>::mStepSize)) * T(1167360.0)) /
      T(351520.0);

  back_time -= 57.0 * integrator<T>::mStepSize / 104.0;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);
  this->mY_n_2 =
      w - k1 * T(439.0 / 216.0) + k2 * T(8.0) - k3 * T(3680.0 / 513.0) +
      (k4 = this->mYp_n_2 * T(integrator<T>::mStepSize)) * T(845.0 / 4104.0);

  back_time -= integrator<T>::mStepSize / 13.0;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);
  this->mY_n_2 =
      w + k1 * T(8.0 / 27.0) - k2 * T(2.0) + k3 * T(3544.0 / 2565.0) -
      k4 * T(1859.0 / 4104.0) +
      (k5 = this->mYp_n_2 * T(integrator<T>::mStepSize)) * T(11.0 / 40.0);

  back_time += integrator<T>::mStepSize * 0.5;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);
  this->mY_n_2 = w - k1 * T(16.0 / 135.0) - k3 * T(6656.0 / 12825.0) -
                 k4 * T(28561.0 / 56430.0) + k5 * T(9.0 / 50.0) -
                 this->mYp_n_2 * T(2.0 * integrator<T>::mStepSize / 55.0);

  back_time -= integrator<T>::mStepSize * 0.5;
  func_ptr->computeStateRate(back_time, this->mY_n_2, this->mYp_n_2);

  w = this->mY_n_2;
  this->mY_n_3 =
      w - (k1 = this->mYp_n_2 * T(integrator<T>::mStepSize)) * T(0.25);

  back_time -= integrator<T>::mStepSize * 0.25;
  func_ptr->computeStateRate(back_time, this->mY_n_3, this->mC);
  this->mY_n_3 -=
      ((k2 = this->mC * T(integrator<T>::mStepSize)) * T(9.0) - k1 * T(5.0)) /
      T(32.0);

  back_time -= integrator<T>::mStepSize * 0.125;
  func_ptr->computeStateRate(back_time, this->mY_n_3, this->mC);
  this->mY_n_3 -=
      (k1 * T(276165.0) - k2 * T(1250865.0) +
       (k3 = this->mC * T(integrator<T>::mStepSize)) * T(1167360.0)) /
      T(351520.0);

  back_time -= 57.0 * integrator<T>::mStepSize / 104.0;
  func_ptr->computeStateRate(back_time, this->mY_n_3, this->mC);
  this->mY_n_3 =
      w - k1 * T(439.0 / 216.0) + k2 * T(8.0) - k3 * T(3680.0 / 513.0) +
      (k4 = this->mC * T(integrator<T>::mStepSize)) * T(845.0 / 4104.0);

  back_time -= integrator<T>::mStepSize / 13.0;
  func_ptr->computeStateRate(back_time, this->mY_n_3, this->mC);
  this->mY_n_3 = w + k1 * T(8.0 / 27.0) - k2 * T(2.0) +
                 k3 * T(3544.0 / 2565.0) - k4 * T(1859.0 / 4104.0) +
                 (k5 = this->mC * T(integrator<T>::mStepSize)) * T(11.0 / 40.0);

  back_time += integrator<T>::mStepSize * 0.5;
  func_ptr->computeStateRate(back_time, this->mY_n_3, this->mC);
  this->mY_n_3 = w - k1 * T(16.0 / 135.0) - k3 * T(6656.0 / 12825.0) -
                 k4 * T(28561.0 / 56430.0) + k5 * T(9.0 / 50.0) -
                 this->mC * T(2.0 * integrator<T>::mStepSize / 55.0);

  this->mP_n = this->mC;  // this makes sure the initial error sweep is zero.
}

template <class T>
void hamming_mod_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                             integrator<T>::mStateRate);

  while (
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    integrator<T>::mTime += integrator<T>::mStepSize;

    this->mP =
        this->mY_n_3 +
        ((integrator<T>::mStateRate + this->mYp_n_2) * T(2.0) - this->mYp_n_1) *
            T(integrator<T>::mStepSize * 4.0 / 3.0);
    this->mM = this->mP - (this->mP_n - this->mC) * T(112.0 / 121.0);

    func_ptr->computeStateRate(integrator<T>::mTime, this->mM, this->mMp);
    this->mC =
        (integrator<T>::mState * T(9.0) - this->mY_n_2 +
         (this->mMp + integrator<T>::mStateRate * T(2.0) - this->mYp_n_1) *
             T(integrator<T>::mStepSize * 3.0)) *
        T(0.125);

    this->mY_n_3.q.swap(this->mY_n_2.q);
    this->mY_n_2.q.swap(this->mY_n_1.q);
    this->mY_n_1.q.swap(integrator<T>::mState.q);
    integrator<T>::mState = this->mC + (this->mP - this->mC) * T(9.0 / 121.0);

    this->mYp_n_2.q.swap(this->mYp_n_1.q);
    this->mYp_n_1.q.swap(integrator<T>::mStateRate.q);
    func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                               integrator<T>::mStateRate);
  }
}

/**
 * This class template implements at Modified Hamming integrator (order 3). This is a fixed-step, multi-step,
 * implicit integrator of order 3. Each integration step entails at most (1 + MaxIterations) evaluations
 * of the state derivative. Neither error control nor divergence detection is performed.
 * Only basic verification of the integration parameters is done and might throw the ReaK::impossible_integration
 * exception.
 * \todo Add the divergence detection.
 */
template <class T>
class hamming_iter_mod_integrator : public hamming_mod_integrator<T> {
 protected:
  double mTolerance;
  unsigned int mMaxIter;

 public:
  void integrate(double aEndTime) override;

  /**
   * Default constructor.
   */
  hamming_iter_mod_integrator(const std::string& aName = "")
      : hamming_mod_integrator<T>(aName), mTolerance(1E-3), mMaxIter(10) {}
  /**
   * Parametrized constructor.
   * \param aName The name of this integrator object.
   * \param aState The initial state vector that the integrator will work with.
   * \param aStartTime The initial time to which the integrator is set.
   * \param aStepSize The time-step used in the integration.
   * \param aGetStateRate A weak pointer to the object that will compute the state derivatives (see
   * ReaK::state_rate_function).
   * \param aTolerance The relative tolerance at which the iterations (within one integration step) is declared as
   * converged.
   * \param aMaxIter The maximum number of iterations to be performed if the iterations cannot converge.
   */
  hamming_iter_mod_integrator(
      const std::string& aName, const ReaK::vect_n<T>& aState,
      double aStartTime, double aStepSize,
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate,
      double aTolerance = 1E-3, unsigned int aMaxIter = 10)
      : hamming_mod_integrator<T>(aName, aState, aStartTime, aStepSize,
                                  aGetStateRate),
        mTolerance(aTolerance),
        mMaxIter(aMaxIter) {}
  /**
   * Default destructor.
   */
  ~hamming_iter_mod_integrator() override = default;

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    hamming_mod_integrator<T>::save(
        A, hamming_mod_integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mTolerance) &
        RK_SERIAL_SAVE_WITH_NAME(mMaxIter);
  }
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    hamming_mod_integrator<T>::load(
        A, hamming_mod_integrator<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mTolerance) &
        RK_SERIAL_LOAD_WITH_NAME(mMaxIter);
  }

  using self = hamming_iter_mod_integrator<T>;
  using base = hamming_mod_integrator<T>;

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2210009, 1,
                              "hamming_iter_mod_integrator", base)
};

template <class T>
void hamming_iter_mod_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) ||
      (integrator<T>::mState.q.size() == 0) ||
      (integrator<T>::mStepSize == 0.0) ||
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  std::shared_ptr<state_rate_function<T>> func_ptr =
      integrator<T>::mGetStateRate.lock();
  if (!func_ptr) {
    throw impossible_integration(integrator<T>::mTime, aEndTime,
                                 integrator<T>::mStepSize);
  }

  ReaK::vect_n<T> prevM(integrator<T>::mState.q.size());
  // Make sure the state-rate vector is initialized.
  func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                             integrator<T>::mStateRate);

  while (
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) ||
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    integrator<T>::mTime += integrator<T>::mStepSize;

    hamming_mod_integrator<T>::mP =
        hamming_mod_integrator<T>::mY_n_3 +
        ((integrator<T>::mStateRate + hamming_mod_integrator<T>::mYp_n_2) *
             T(2.0) -
         hamming_mod_integrator<T>::mYp_n_1) *
            T(integrator<T>::mStepSize * 4.0 / 3.0);
    prevM = hamming_mod_integrator<T>::mP -
            (hamming_mod_integrator<T>::mP_n - hamming_mod_integrator<T>::mC) *
                T(112.0 / 121.0);

    double ErrorEst = 2.0 * this->mTolerance;
    unsigned int i = 0;
    while ((ErrorEst > this->mTolerance) && (i < this->mMaxIter)) {
      func_ptr->computeStateRate(integrator<T>::mTime, prevM,
                                 hamming_mod_integrator<T>::mMp);
      hamming_mod_integrator<T>::mC =
          (integrator<T>::mState * T(9.0) - hamming_mod_integrator<T>::mY_n_2 +
           (hamming_mod_integrator<T>::mMp +
            integrator<T>::mStateRate * T(2.0) -
            hamming_mod_integrator<T>::mYp_n_1) *
               T(integrator<T>::mStepSize * 3.0)) *
          T(0.125);

      hamming_mod_integrator<T>::mM =
          hamming_mod_integrator<T>::mC +
          (hamming_mod_integrator<T>::mP - hamming_mod_integrator<T>::mC) *
              T(9.0 / 121.0);
      ErrorEst = 0.0;
      for (unsigned int j = 0; j < prevM.q.size(); ++j) {
        ErrorEst += sqr_mag(hamming_mod_integrator<T>::mM.q[j] - prevM.q[j]);
      }
      ErrorEst = std::sqrt(ErrorEst);
      prevM.q.swap(hamming_mod_integrator<T>::mM.q);
      ++i;
    };
    hamming_mod_integrator<T>::mY_n_3.q.swap(
        hamming_mod_integrator<T>::mY_n_2.q);
    hamming_mod_integrator<T>::mY_n_2.q.swap(
        hamming_mod_integrator<T>::mY_n_1.q);
    hamming_mod_integrator<T>::mY_n_1.q.swap(integrator<T>::mState.q);
    integrator<T>::mState.q.swap(prevM.q);

    hamming_mod_integrator<T>::mYp_n_2.q.swap(
        hamming_mod_integrator<T>::mYp_n_1.q);
    hamming_mod_integrator<T>::mYp_n_1.q.swap(integrator<T>::mStateRate.q);
    func_ptr->computeStateRate(integrator<T>::mTime, integrator<T>::mState,
                               integrator<T>::mStateRate);
  }
}

}  // namespace ReaK

#endif
