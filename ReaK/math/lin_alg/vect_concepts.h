/**
 * \file vect_concepts.h
 *
 * This library defines the concepts related to generic vector types used in the ReaK library.
 * These concepts are based on the principle of minimum requirement, they were designed to
 * require only the minimum set of valid expressions that will be used by the algorithms
 * that pertain to generic vector types.
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

#ifndef REAK_MATH_LIN_ALG_VECT_CONCEPTS_H_
#define REAK_MATH_LIN_ALG_VECT_CONCEPTS_H_

#include "ReaK/core/base/defs.h"

#include <concepts>
#include <iterator>
#include <type_traits>
#include <vector>

#include "ReaK/math/lin_alg/vect_traits.h"

namespace ReaK {

/**
 * This concept class defines what makes a vector a readable vector, that is,
 * a vector whose elements can be read.
 *
 * Valid Expressions:
 *
 * e = v[i];   can be indexed to be read.
 *
 * s = v.size();   the size of the vector can be obtained.
 *
 * cit = v.begin();   a const-iterator to the first vector element can be obtained.
 *
 * ++cit;   the const-iterator can be incremented.
 *
 * cit = v.end();   a const-iterator to the one-past-last vector element can be obtained.
 *
 * \tparam Vector The vector type.
 */
template <typename Vector>
concept ReadableVector = requires(const Vector& v) {
  { v[0] } -> std::convertible_to<vect_value_type_t<Vector>>;
  { v.size() } -> std::integral<>;
  { v.begin() } -> std::input_iterator<>;
  { v.end() } -> std::input_iterator<>;
  { v.begin() != v.end() } -> std::convertible_to<bool>;
};

// Legacy
template <typename Vector>
static constexpr bool is_readable_vector_v = ReadableVector<Vector>;

/**
 * This concept class defines what makes a vector a writable vector, that is,
 * a vector whose elements can be written.
 *
 * Valid Expressions (in addition to those of ReadableVector):
 *
 * v[i] = e;   can be indexed to be written.
 *
 * \tparam Vector The vector type.
 */
template <typename Vector>
concept WritableVector = ReadableVector<Vector>&& requires(Vector& v) {
  { v[0] } -> std::assignable_from<vect_value_type_t<Vector>>;
};

// Legacy
template <typename Vector>
static constexpr bool is_writable_vector_v = WritableVector<Vector>;

/**
 * This concept class defines what makes a vector a resizable vector, that is,
 * a vector whose size can be changed at run-time.
 *
 * Valid Expressions:
 *
 * v.resize(s);   can be resized.
 *
 * \tparam Vector The vector type.
 */
template <typename Vector>
concept ResizableVector = ReadableVector<Vector>&& requires(Vector& v, int sz) {
  {v.resize(sz)};
};

// Legacy
template <typename Vector>
static constexpr bool is_resizable_vector_v = ResizableVector<Vector>;

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_CONCEPTS_H_
