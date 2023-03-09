/**
 * \file mat_num_exceptions.h
 *
 * This library defines a number of exception classes to represent the typical
 * exceptions that could be thrown during matrix numerical methods.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_NUM_EXCEPTIONS_H_
#define REAK_MATH_LIN_ALG_MAT_NUM_EXCEPTIONS_H_

#include <exception>
#include <sstream>

namespace ReaK {

/**
 * This exception signifies an error with a mathematical singularity in a matrix involved
 * in the operation that threw this exception.
 */
class singularity_error : public std::exception {
 public:
  std::string message;  ///< Message string that identifies the singular matrix.

  /**
   * Default constructor.
   * \param aMatrixName the name of the matrix which is singular.
   */
  explicit singularity_error(const std::string& aMatrixName)
      : message(std::string("Singularity detected! For matrix ") +
                aMatrixName) {}
  /**
   * Destructor.
   */
  ~singularity_error() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

/**
 * This exception signifies that the throwing iterative numerical method has reached maximum
 * iteration before having reached any "successful" terminal condition.
 */
class maximum_iteration : public std::exception {
 public:
  std::string message;  ///< Error message string.

  /**
   * Default constructor.
   * \param aMaxIter the iteration limit reached.
   */
  explicit maximum_iteration(unsigned int aMaxIter) {
    std::stringstream ss;
    ss << "Maximum iteration reached at " << aMaxIter;
    message = ss.str();
  }
  /**
   * Destructor.
   */
  ~maximum_iteration() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_NUM_EXCEPTIONS_H_
