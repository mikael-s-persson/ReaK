/**
 * \file vect_traits.hpp
 *
 * This library declares the vector traits that are used throughout the ReaK library.
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

#ifndef REAK_VECT_TRAITS_HPP
#define REAK_VECT_TRAITS_HPP

#include <ReaK/core/base/defs.hpp>

namespace ReaK {

/**
 * This trait class defines all the traits of a vector in the ReaK library.
 * \tparam Vector A vector type.
 */
template <typename Vector>
struct vect_traits {
  /** The type of the elements of the vector. */
  using value_type = typename Vector::value_type;
  /** The type of a reference to an element of the vector. */
  using reference = typename Vector::reference;
  /** The type of a const-reference to an element of the vector. */
  using const_reference = typename Vector::const_reference;
  /** The type of a pointer to an element of the vector. */
  using pointer = typename Vector::pointer;
  /** The type of a const-pointer to an element of the vector. */
  using const_pointer = typename Vector::const_pointer;
  /** The type of the allocator used by the vector (void if none). */
  using allocator_type = typename Vector::allocator_type;

  /** The type of a iterator through the vector. */
  using iterator = typename Vector::iterator;
  /** The type of a const-iterator through the vector. */
  using const_iterator = typename Vector::const_iterator;

  /** The type to describe the size of the vector or an index to it. */
  using size_type = typename Vector::size_type;
  /** The type to describe the difference between two indices. */
  using difference_type = typename Vector::difference_type;

  /** The dimension of the vector (0 if not known at compile-time). */
  static constexpr std::size_t dimensions = Vector::dimensions;
};

template <typename Vector, typename = void>
struct vect_value_type {
  using type = double;
};

// Only get value type without relying on vect_traits.
template <typename Vector>
struct vect_value_type<Vector,
                       std::void_t<decltype(std::declval<Vector>()[0])>> {
  using type =
      std::decay_t<decltype(std::declval<std::add_const_t<Vector>>()[0])>;
};

template <typename Vector>
using vect_value_type_t = typename vect_value_type<Vector>::type;

/**
 * This meta-function provides the vector type needed as a destination type for a
 * copy (deep-copy) of a given vector type.
 * \tparam Vector A vector type.
 */
template <typename Vector>
struct vect_copy {
  /** The type of a copy of the given vector type. */
  using type = Vector;
};

template <typename Vector>
using vect_copy_t = typename vect_copy<Vector>::type;

}  // namespace ReaK

#endif
