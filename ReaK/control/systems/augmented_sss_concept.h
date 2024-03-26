/**
 * \file augmented_sss_concept.h
 *
 * This library defines the traits class and concept that represent an augmented
 * state-space system. This type of systems is characterized by a set of actual
 * states (dynamic system states) and a set of quasi-constant parameters to be
 * identified.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_CONTROL_SYSTEMS_AUGMENTED_SSS_CONCEPT_H_
#define REAK_CONTROL_SYSTEMS_AUGMENTED_SSS_CONCEPT_H_

#include "ReaK/core/base/defs.h"

#include <concepts>
#include <type_traits>

namespace ReaK::ctrl {

/**
 * This traits class defines the characteristics of a discrete-time state-space system.
 * \tparam DiscreteSystem The discrete-time state-space system type for which the traits are sought.
 */
template <typename T>
struct augmented_sss_traits {
  /** This constant describes the dimensions of the output vector (0 if not known at compile-time). */
  static constexpr std::size_t actual_state_dimensions =
      T::actual_state_dimensions;
};

/**
 * This class template defines the concept for an augmented state-space system as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in augmented_sss_traits, an augmented state-space
 * system should provide a number of valid expressions.
 *
 * Valid expressions:
 *
 * x = sys.get_actual_state_dimensions();  The state-space system (sys) can deliver the dimensions count (x) for the
 *actual (dynamic) states of the system.
 */
template <typename T>
concept AugmentedSystem = requires(const T& sys) {
  { sys.get_actual_state_dimensions() } -> std::integral;
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_AUGMENTED_SSS_CONCEPT_H_
