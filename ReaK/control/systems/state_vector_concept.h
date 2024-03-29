/**
 * \file state_vector_concept.h
 *
 * This library provides the traits class and the concept definition for a
 * state-vector as used within the ReaK::ctrl namespace. A state-vector is
 * an abstraction of the quantity that describes the state of a state-space
 * system (see SSSystemConcept).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_CONTROL_SYSTEMS_STATE_VECTOR_CONCEPT_H_
#define REAK_CONTROL_SYSTEMS_STATE_VECTOR_CONCEPT_H_

#include "ReaK/math/lin_alg/mat_slices.h"
#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include <concepts>
#include <utility>

namespace ReaK::ctrl {

/**
 * This class template is the traits class that defines the traits that a
 * state-vector should have.
 */
template <typename V>
struct state_vector_traits {
  /** This is the type of the state-vector descriptor, usually the same as V. */
  using state_type = int;
  /** This is the type that describes the difference between two state-vectors. */
  using state_difference_type = int;
  /** This is the value-type of the elements of the state-vector. */
  using value_type = int;
  /** This is the type that describes the size of the state-vector. */
  using size_type = int;

  /** This constant describes the dimension of the state-vector (0 if only known at run-time). */
  static constexpr std::size_t dimensions = 0;
};
template <typename V>
concept HasAllStateVectorTraits = requires {
  typename V::state_type;
  typename V::state_difference_type;
  typename V::value_type;
  typename V::size_type;
  V::dimensions;
};
template <HasAllStateVectorTraits V>
struct state_vector_traits<V> {
  /** This is the type of the state-vector descriptor, usually the same as V. */
  using state_type = typename V::state_type;
  /** This is the type that describes the difference between two state-vectors. */
  using state_difference_type = typename V::state_difference_type;
  /** This is the value-type of the elements of the state-vector. */
  using value_type = typename V::value_type;
  /** This is the type that describes the size of the state-vector. */
  using size_type = typename V::size_type;

  /** This constant describes the dimension of the state-vector (0 if only known at run-time). */
  static constexpr std::size_t dimensions = V::dimensions;
};

/**
 * This class template defines the concept that state-vectors should model.
 *
 * Required concepts:
 *
 * the state-vector's state_difference_type should model ReadableVectorConcept.
 *
 * Valid expressions (state_type s, state_difference_type ds, value_type v, size_type sz):
 *
 * ds = diff(s,s);  The difference between state-vectors is obtained by the diff() function.
 *
 * s = add(s,ds);  A state-difference can be added to a state-vector with the add() function.
 *
 * ds = v * ds;  A state-difference is scalable.
 *
 * ds = ds + ds;  State-differences can be added.
 *
 * ds += ds;  State-differences can be added and stored.
 *
 * ds = ds - ds;  State-differences can be subtracted.
 *
 * ds -= ds;  State-differences can be subtracted and stored.
 *
 * ds = -ds;  A state-difference can be negated.
 *
 * ds = unit(ds);  A state-difference can be made into a unit-vector.
 *
 * v = norm_2(ds);  A state-difference can be taken the norm of.
 */
template <typename T>
concept StateVector = requires(std::decay_t<T> s,
                               std::decay_t<decltype(diff(s, s))> ds) {
  { to_vect<double>(ds) } -> ReadableVector;
  from_vect<decltype(ds)>(to_vect<double>(ds));
  ds = diff(s, s);
  s = add(s, ds);
  ds = double{} * ds;
  ds *= double{};
  ds = ds + ds;
  ds += ds;
  ds = ds - ds;
  ds -= ds;
  ds = -ds;
  ds = unit(ds);
  { norm_2(ds) } -> std::convertible_to<double>;
};

template <typename T, unsigned int Size>
struct state_vector_traits<vect<T, Size>> {
  using state_type = vect<T, Size>;
  using state_difference_type = vect<T, Size>;
  using value_type = typename vect_traits<vect<T, Size>>::value_type;
  using size_type = typename vect_traits<vect<T, Size>>::size_type;

  static constexpr std::size_t dimensions = vect_traits<state_type>::dimensions;
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_STATE_VECTOR_CONCEPT_H_
