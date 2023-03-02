/**
 * \file sss_exceptions.hpp
 *
 * This library defines a number of exceptions related to the definition of state-space systems.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_SSS_EXCEPTIONS_HPP
#define REAK_SSS_EXCEPTIONS_HPP

#include <exception>
#include <string>

namespace ReaK::ctrl {

/**
 * This exception is thrown whenever there is an incoherency in the definition of a
 * state-space system, such as mismatching dimensions of system matrices or singular
 * matrix definitions (singular mass-matrix for example).
 */
class system_incoherency : public std::exception {
 public:
  std::string message;  ///< Message string that identifies the singular matrix.

  /**
   * Default constructor.
   * \param aMessage the message corresponding to the incoherency.
   */
  explicit system_incoherency(const std::string& aMessage)
      : message(
            std::string("State space system is incoherent, with message '") +
            aMessage + "'") {}
  /**
   * Destructor.
   */
  ~system_incoherency() noexcept override = default;

  /**
   * Gets the error message.
   * \return c_string of the error message.
   */
  const char* what() const noexcept override { return message.c_str(); }
};

}  // namespace ReaK::ctrl

#endif
