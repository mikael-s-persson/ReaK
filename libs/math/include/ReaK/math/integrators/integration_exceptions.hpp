
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

#ifndef REAK_INTEGRATION_EXCEPTIONS_HPP
#define REAK_INTEGRATION_EXCEPTIONS_HPP

#include <exception>
#include <sstream>

namespace ReaK {

/**
 * Exception thrown whenever the parameters of the integration were invalid.
 */
class impossible_integration : public std::exception {
 public:
  double mTime;
  double mEndTime;
  double mTimeStep;
  std::string message;

  /**
   * Constructor.
   * \param aTime integration time at which the exception occurred.
   * \param aEndTime integration time at which the integration was supposed to end.
   * \param aTimeStep integration time step at the time of the exception.
   */
  impossible_integration(double aTime, double aEndTime, double aTimeStep)
      : mTime(aTime), mEndTime(aEndTime), mTimeStep(aTimeStep) {
    std::stringstream sstr;
    sstr << "Integration is impossible! Error occurred at time " << mTime
         << " towards end time of " << mEndTime << " with time step "
         << mTimeStep << ".";
    message = sstr.str();
  }
  /**
   * Destructor.
   */
  ~impossible_integration() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

/**
 * Exception thrown whenever the integration reached the tolerance without being able to improve the results (reaching a
 * minimum time step).
 */
class untolerable_integration : public std::exception {
 public:
  double mTolerance;
  double mErrorEstimate;
  int mDOF;
  double mTimeStep;
  double mTime;
  std::string message;

  /**
   * Constructor.
   * \param aTolerance tolerance that was not reachable by the integrator.
   * \param eErrorEstimate estimated error at the time of the exception.
   * \param aDOF degree-of-freedom at which the maximum and untolerable estimated error occurred.
   * \param aTimeStep integration time step at the time of the exception.
   * \param aTime integration time at which the exception occurred.
   */
  untolerable_integration(double aTolerance, double aErrorEstimate, int aDOF,
                          double aTimeStep, double aTime)
      : mTolerance(aTolerance),
        mErrorEstimate(aErrorEstimate),
        mDOF(aDOF),
        mTimeStep(aTimeStep),
        mTime(aTime) {
    std::stringstream sstr;
    sstr << "Integration was deemed untolerable! Error occurred at time "
         << mTime << " with current time step of " << mTimeStep
         << " in violation of a tolerance of " << mTolerance
         << " by an estimated error of " << mErrorEstimate
         << " at state element " << mDOF << ".";
    message = sstr.str();
  }
  /**
   * Destructor.
   */
  ~untolerable_integration() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

/**
 * Exception thrown whenever the integration reached an invalid state such as not-a-number (NaN).
 */
class invalid_state_derivative : public std::exception {
 public:
  int mDOF;
  double mTime;
  std::string message;

  /**
   * Constructor.
   * \param aDOF degree-of-freedom that went invalid.
   * \param aTime integration time at which the exception occurred.
   */
  invalid_state_derivative(int aDOF, double aTime) : mDOF(aDOF), mTime(aTime) {
    std::stringstream sstr;
    sstr << "Integration has reached an invalid state derivative! Error "
            "occurred at time "
         << mTime << " with state element " << mDOF << ".";
    message = sstr.str();
  }
  /**
   * Destructor.
   */
  ~invalid_state_derivative() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

}  // namespace ReaK

#endif
