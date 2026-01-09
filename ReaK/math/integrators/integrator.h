/**
 * \file integrator.h
 *
 * This library declares the basic abstract classes for fixed step and variable
 * step integrator for systems of ODEs. It provides basic capabilities for setting
 * up the state vector and the state's rate of change vector.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_MATH_INTEGRATORS_INTEGRATOR_H_
#define REAK_MATH_INTEGRATORS_INTEGRATOR_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include "ReaK/math/integrators/integration_exceptions.h"

/** Main namespace for ReaK */
namespace ReaK {

/**
 * This class is the function-object interface for a state equation (or state time-derivative computation).
 */
template <class T>
class state_rate_function : public virtual shared_object {
 public:
  /**
   * This function computes the time-derivative of the state vector.
   *
   * \param aTime current integration time
   * \param aState current state vector
   * \param aStateRate holds, as output, the time-derivative of the state vector
   */
  virtual void computeStateRate(double aTime, const ReaK::vect_n<T>& aState,
                                ReaK::vect_n<T>& aStateRate) = 0;

  ~state_rate_function() override = default;

  using self = state_rate_function<T>;
  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2200001, 1, "state_rate_function",
                              shared_object)
};

/**
 * This class is the function-object interface for a state equation (or state time-derivative computation).
 */
template <class T>
class state_rate_function_with_io : public state_rate_function<T> {
 public:
  /**
   * This function computes the output-vector corresponding to a state vector.
   *
   * \param aTime current integration time
   * \param aState current state vector
   * \param aOutput holds, as output, the output-vector
   */
  virtual void computeOutput(double aTime, const ReaK::vect_n<T>& aState,
                             ReaK::vect_n<T>& aOutput) = 0;

  /**
   * This function sets the input-vector.
   *
   * \param aInput current input-vector
   */
  virtual void setInput(const ReaK::vect_n<T>& aInput) = 0;

  ~state_rate_function_with_io() override = default;

  using self = state_rate_function_with_io<T>;
  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2200002, 1,
                              "state_rate_function_with_io",
                              state_rate_function<T>)
};

/**
 * This class is the basis for all fixed time step integrators. It features a
 * state vector with a rate of change vector (specified by vectors of pointers
 * to the variables), and the time span (StartTime, EndTime, and StepSize).
 */
template <class T>
class integrator : public named_object {
 protected:
  /// Holds the current state vector.
  vect_n<T> mState;
  /// Holds the last computed time-derivative of the state vector.
  vect_n<T> mStateRate;

  /// Current integration time.
  double mTime;
  /// Current integration time step.
  double mStepSize;

  /// Pointer to a function-object that computes the time-derivative of the state vector.
  std::weak_ptr<state_rate_function<T>> mGetStateRate;

 public:
  /// Gets the state vector iterator to the current state vector.
  auto getStateBegin() const { return mState.q.begin(); }
  /// Gets the state vector iterator to the current state vector's end element.
  auto getStateEnd() const { return mState.q.end(); }

  /// Set the integration time step to aNewStepSize.
  virtual void setStepSize(double aNewStepSize) { mStepSize = aNewStepSize; }
  /// Sets the integration time to aNewTime.
  virtual void setTime(double aNewTime) { mTime = aNewTime; }
  /// Gets the current integration time.
  virtual double getTime() { return mTime; }

  /// Clears the state vector, makes it zero in size.
  virtual void clearStateVector() {
    mState.q.clear();
    mStateRate.q.clear();
  }
  /// Adds an element to the state vector and initializes its value to Element.
  virtual void addStateElement(T Element) {
    mState.q.push_back(Element);
    mStateRate.q.push_back(T(0.0));
  }
  /// Adds elements to the state vector and initializes their value to Elements.
  virtual void addStateElements(const ReaK::vect_n<T>& Elements) {
    mState.q.insert(mState.q.end(), Elements.q.begin(), Elements.q.end());
    mStateRate.q.insert(mStateRate.q.end(), Elements.q.size(), T(0.0));
  }

  /// Set the function-object pointer to aGetStateRate.
  virtual void setStateRateFunc(
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate) {
    mGetStateRate = aGetStateRate;
  }

  /**
   * Performs the integration.
   * \param aEndTime the time up to which the integration shall go on. Note that it is only guaranteed to exit the
   * function to within one integration time step passed aEndTime.
   * \throw untolerable_integration if the tolerance of the error estimate could not be kept.
   * \throw invalid_state_derivative if the results of the state time-derivative led to a NaN or infinite value.
   * \throw impossible_integration if the integrator is in an invalid state such as expired function-object pointer for
   * example.
   */
  virtual void integrate(double aEndTime) = 0;

  /// Default constructor.
  explicit integrator(const std::string& aName = "")
      : named_object(),
        mState(),
        mStateRate(),
        mTime(),
        mStepSize(),
        mGetStateRate() {
    set_name(aName);
  }
  /// Parametrized constructor, see data members for meaning of parameters.
  integrator(const std::string& aName, const ReaK::vect_n<T>& aState,
             double aStartTime, double aStepSize,
             const std::weak_ptr<state_rate_function<T>>& aGetStateRate)
      : named_object(),
        mState(aState),
        mStateRate(aState.q.size()),
        mTime(aStartTime),
        mStepSize(aStepSize),
        mGetStateRate(aGetStateRate) {
    set_name(aName);
  }
  /// Default virtual destructor.
  ~integrator() override = default;

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, ReaK::named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mState) & RK_SERIAL_SAVE_WITH_NAME(mTime) &
        RK_SERIAL_SAVE_WITH_NAME(mStepSize) &
        RK_SERIAL_SAVE_WITH_NAME(mGetStateRate);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, ReaK::named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mState) & RK_SERIAL_LOAD_WITH_NAME(mTime) &
        RK_SERIAL_LOAD_WITH_NAME(mStepSize) &
        RK_SERIAL_LOAD_WITH_NAME(mGetStateRate);
    mStateRate.q.resize(mState.q.size());
  }

  using self = integrator<T>;

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2210000, 1, "integrator", named_object)
};

/**
 * This class is the basis for all variable time step integrators. It features
 * a maximum and minimum time step size as well as a tolerance for error control.
 */
template <class T>
class variable_step_integrator : public integrator<T> {
 protected:
  /// Maximum allowable integration time step.
  double mMaxStepSize;
  /// Minimum allowable integration time step.
  double mMinStepSize;
  /// Tolerance by which the estimated error is used to assess changes to the time step.
  double mTolerance;

 public:
  /// Sets maximum allowable integration time step to aNewMaxStepSize.
  virtual void setMaxStepSize(double aNewMaxStepSize) {
    mMaxStepSize = aNewMaxStepSize;
  }
  /// Sets minimum allowable integration time step to aNewMinStepSize.
  virtual void setMinStepSize(double aNewMinStepSize) {
    mMinStepSize = aNewMinStepSize;
  };
  /// Sets tolerance to aNewTolerance.
  virtual void setTolerance(double aNewTolerance) {
    mTolerance = aNewTolerance;
  }

  /// Default constructor.
  explicit variable_step_integrator(const std::string& aName = "")
      : integrator<T>(aName),
        mMaxStepSize(1.0),
        mMinStepSize(1E-6),
        mTolerance(1E-4) {}
  /// Parametrized constructor, see data members for meaning of parameters.
  variable_step_integrator(
      const std::string& aName, const ReaK::vect_n<T>& aState,
      double aStartTime, double aInitialStepSize,
      const std::weak_ptr<state_rate_function<T>>& aGetStateRate,
      double aMaxStepSize, double aMinStepSize, double aTolerance)
      : integrator<T>(aName, aState, aStartTime, aInitialStepSize,
                      aGetStateRate),
        mMaxStepSize(aMaxStepSize),
        mMinStepSize(aMinStepSize),
        mTolerance(aTolerance) {}
  /// Default virtual destructor.
  ~variable_step_integrator() override = default;

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    integrator<T>::save(A, integrator<T>::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mMaxStepSize) &
        RK_SERIAL_SAVE_WITH_NAME(mMinStepSize) &
        RK_SERIAL_SAVE_WITH_NAME(mTolerance);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    integrator<T>::load(A, integrator<T>::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mMaxStepSize) &
        RK_SERIAL_LOAD_WITH_NAME(mMinStepSize) &
        RK_SERIAL_LOAD_WITH_NAME(mTolerance);
  }

  using self = variable_step_integrator<T>;
  using base = integrator<T>;

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2220000, 1, "variable_step_integrator",
                              base)
};
}  // namespace ReaK

#endif  // REAK_MATH_INTEGRATORS_INTEGRATOR_H_
