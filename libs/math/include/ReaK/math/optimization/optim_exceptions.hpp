/**
 * \file optim_exceptions.hpp
 *
 * This library contains the exception classes which are relevant to the optimization algorithms.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_OPTIM_EXCEPTIONS_HPP
#define REAK_OPTIM_EXCEPTIONS_HPP

#include <ReaK/core/base/defs.hpp>

#include <exception>
#include <utility>

namespace ReaK::optim {

/**
 * This exception signifies that the constraints imposed on the search domain render
 * the optimization infeasible (i.e. the search domain is empty).
 */
class infeasible_problem : public std::exception {
 public:
  std::string message;  ///< Message string.

  /**
   * Default constructor.
   * \param aMessage The error message.
   */
  explicit infeasible_problem(std::string aMessage)
      : message(std::move(aMessage)) {}
  /**
   * Destructor.
   */
  ~infeasible_problem() noexcept override = default;
  ;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

/**
 * This exception signifies that the constraints imposed on the search domain are
 * not sufficient to make it a bounded domain (may be thrown by methods that expect a closed domain).
 */
class unbounded_problem : public std::exception {
 public:
  std::string message;  ///< Message string.

  /**
   * Default constructor.
   * \param aMessage The error message.
   */
  explicit unbounded_problem(std::string aMessage)
      : message(std::move(aMessage)) {}
  /**
   * Destructor.
   */
  ~unbounded_problem() noexcept override = default;
  ;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

/**
 * This exception signifies that the optimization problem is not proper (i.e. some
 * predicate of the problem definition is violated).
 */
class improper_problem : public std::exception {
 public:
  std::string message;  ///< Message string.

  /**
   * Default constructor.
   * \param aMessage The error message.
   */
  explicit improper_problem(std::string aMessage)
      : message(std::move(aMessage)) {}
  /**
   * Destructor.
   */
  ~improper_problem() noexcept override = default;
  ;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

}  // namespace ReaK::optim

#endif
