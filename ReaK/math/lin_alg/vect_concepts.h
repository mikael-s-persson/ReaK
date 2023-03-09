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

#include <vector>
#include "ReaK/math/lin_alg/vect_traits.h"

#include "boost/concept_check.hpp"

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
struct ReadableVectorConcept {
  Vector v;

  typename vect_traits<Vector>::size_type s;
  vect_value_type_t<Vector> cr;

  typename vect_traits<Vector>::const_iterator it;

  void constraints() {
    cr = v[0];  // can be indexed and given an rvalue
    s = v.size();
    it = v.begin();
    ++it;
    bool b = (it != v.end());
    RK_UNUSED(b);
  }
};

/**
 * This meta-function evaluates whether a Vector class fulfills the ReadableVectorConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
template <typename Vector>
struct is_readable_vector {
  static constexpr bool value = false;
  using type = is_readable_vector<Vector>;
};

template <typename Vector>
static constexpr bool is_readable_vector_v = is_readable_vector<Vector>::value;

template <typename T>
struct is_readable_vector<std::vector<T>> {
  static constexpr bool value = true;
  using type = is_readable_vector<std::vector<T>>;
};

/**
 * This concept class defines what makes a vector a writable vector, that is,
 * a vector whose elements can be written.
 *
 * Valid Expressions (in addition to those of ReadableVectorConcept):
 *
 * v[i] = e;   can be indexed to be written.
 *
 * \tparam Vector The vector type.
 */
template <typename Vector>
struct WritableVectorConcept
    : ReadableVectorConcept<Vector> {  // must also be readable.

  vect_value_type_t<Vector> r;

  void constraints() {
    this->v[0] = r;  // can be indexed and given an lvalue
  }
};

/**
 * This meta-function evaluates whether a Vector class fulfills the WritableVectorConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
template <typename Vector>
struct is_writable_vector {
  static constexpr bool value = false;
  using type = is_writable_vector<Vector>;
};

template <typename Vector>
static constexpr bool is_writable_vector_v = is_writable_vector<Vector>::value;

template <typename T>
struct is_writable_vector<std::vector<T>> {
  static constexpr bool value = true;
  using type = is_writable_vector<std::vector<T>>;
};

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
struct ResizableVectorConcept {
  Vector v;

  typename vect_traits<Vector>::size_type sz;

  BOOST_CONCEPT_USAGE(ResizableVectorConcept) { v.resize(sz); };
};

/**
 * This meta-function evaluates whether a Vector class fulfills the ResizableVectorConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
template <typename Vector>
struct is_resizable_vector {
  static constexpr bool value = false;
  using type = is_resizable_vector<Vector>;
};

template <typename Vector>
static constexpr bool is_resizable_vector_v =
    is_resizable_vector<Vector>::value;

template <typename T>
struct is_resizable_vector<std::vector<T>> {
  static constexpr bool value = true;
  using type = is_resizable_vector<std::vector<T>>;
};

/**
 * This concept class defines what makes a vector a resizable vector, that is,
 * a vector whose size can be changed at run-time.
 *
 * Valid Expressions (in addition to those of ResizableVectorConcept):
 *
 * al = v.get_allocator();   the allocator object can be obtained.
 *
 * \tparam Vector The vector type.
 */
template <typename Vector>
struct DynAllocVectorConcept : ResizableVectorConcept<Vector> {
  Vector v;

  typename vect_traits<Vector>::allocator_type al;

  BOOST_CONCEPT_USAGE(DynAllocVectorConcept) { al = v.get_allocator(); };
};

/**
 * This meta-function evaluates whether a Vector class fulfills the DynAllocVectorConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
template <typename Vector>
struct has_allocator_vector {
  static constexpr bool value = false;
  using type = has_allocator_vector<Vector>;
};

template <typename Vector>
static constexpr bool has_allocator_vector_v =
    has_allocator_vector<Vector>::value;

template <typename T>
struct has_allocator_vector<std::vector<T>> {
  static constexpr bool value = true;
  using type = has_allocator_vector<std::vector<T>>;
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_CONCEPTS_H_
