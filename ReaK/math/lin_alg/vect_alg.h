/**
 * \file vect_alg.h
 *
 * This library provides classes to represent fixed and variable sized vectors. These
 * vectors can be used in linear algebra expressions, with standard vector algebra semantics.
 * Generally, vectors are considered to be column-vectors (but pre-multiplication to a matrix
 * turns assumes them to a row-vector for the duration of the operation). The interface
 * of the vectors are also generally compatible with STL containers (and the underlying
 * container for the variable-size vector is a std::vector class template).
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

#ifndef REAK_MATH_LIN_ALG_VECT_ALG_H_
#define REAK_MATH_LIN_ALG_VECT_ALG_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "ReaK/core/base/defs.h"
#include "ReaK/core/rtti/so_register_type.h"
#include "ReaK/core/rtti/typed_primitives.h"
#include "ReaK/core/serialization/archiver.h"

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_index_iterator.h"
#include "ReaK/math/lin_alg/vect_views.h"

namespace ReaK {

/**
 * This class implements a fixed-index vector component of primitive type.
 * This class is mainly meant to be used with the "vect" class template.
 */
template <typename T, unsigned int Index>
class vect_component {
 public:
  using self = vect_component<T, Index>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;

  value_type q;

  explicit vect_component(const_reference Q = 0) noexcept : q(Q) {}

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  value_type operator[](size_type i) const noexcept {
    return (i == Index ? q : value_type());
  }

  /**
   * Call operator, accessor for read only.
   * \test PASSED
   */
  value_type operator()(size_type i) const noexcept {
    return (i == Index ? q : value_type());
  }

  friend self operator+(const self& lhs, const self& rhs) noexcept {
    return self(lhs.q + rhs.q);
  }

  friend self operator-(const self& lhs, const self& rhs) noexcept {
    return self(lhs.q - rhs.q);
  }

  friend self operator-(const self& lhs) noexcept { return self(-lhs.q); }

  self& operator+=(const self& rhs) noexcept {
    q += rhs.q;
    return *this;
  }

  self& operator-=(const self& rhs) noexcept {
    q -= rhs.q;
    return *this;
  }

  friend self operator*(const self& lhs, const_reference rhs) noexcept {
    return self(lhs.q * rhs);
  }

  friend self operator*(const_reference lhs, const self& rhs) noexcept {
    return self(lhs * rhs.q);
  }

  self& operator*=(const_reference rhs) noexcept {
    q *= rhs;
    return *this;
  }
};

static const vect_component<double, 0> vect_i = vect_component<double, 0>(1.0);
static const vect_component<double, 1> vect_j = vect_component<double, 1>(1.0);
static const vect_component<double, 2> vect_k = vect_component<double, 2>(1.0);

template <typename T, typename Allocator>
class vect_n;  // forward-declaration.

/**
 * This class implements a fixed-size templated vector class which holds components of primitive type.
 */
template <typename T, unsigned int Size>
class vect {
 public:
  using self = vect<T, Size>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using allocator_type = std::allocator<T>;

  using iterator = typename std::array<T, Size>::iterator;
  using const_iterator = typename std::array<T, Size>::const_iterator;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t dimensions = Size;

  /// Components of the vector.
  std::array<T, Size> q = {};

  /**
   * Returns the size of the vector.
   */
  size_type size() const noexcept { return Size; }
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const noexcept { return Size; }
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const noexcept { return Size; }
  /**
   * Resizes the vector.
   */
  void resize(size_type sz, T c = T()) const noexcept {}
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return false; }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(size_type sz) const noexcept {}

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() noexcept { return q.begin(); }
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return q.begin(); }
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() noexcept { return q.end(); }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return q.end(); }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor: sets all to zero.
   * \test PASSED
   */
  vect() noexcept {
    for (auto& v : q) {
      v = value_type{};
    }
  }

  template <typename U, typename Allocator>
  vect(const vect_n<U, Allocator>& V) noexcept;  // NOLINT

  /**
   * Constructor from an array of values of type T.
   * \test PASSED
   */
  explicit vect(const_pointer Q) noexcept {
    for (auto& v : q) {
      v = *Q++;
    }
  }

  vect(const self& V) noexcept = default;

  template <typename U, unsigned int OtherSize>
  explicit vect(const vect<U, OtherSize>& rhs) noexcept {
    static_assert(Size >= OtherSize);
    auto it = q.begin();
    for (auto rhs_it = rhs.begin(); rhs_it < rhs.end(); ++it, ++rhs_it) {
      *it = *rhs_it;
    }
    for (; it < q.end(); ++it) {
      *it = value_type{};
    }
  }

  template <typename U, unsigned int Index>
  explicit vect(const vect_component<U, Index>& rhs) noexcept : vect() {
    static_assert(Size > Index);
    q[Index] = rhs.q;
  }

 private:
  template <typename Arg1>
  static void set_value_impl(iterator& it, const Arg1& a1) noexcept {
    *it++ = value_type(a1);
  }

  template <typename Arg1, typename... Args>
  static void set_value_impl(iterator& it, const Arg1& a1,
                             const Args&... tail) noexcept {
    *it++ = value_type(a1);
    set_value_impl(it, tail...);
  }

 public:
  /**
   * Constructor for Size values.
   * \test PASSED
   */
  template <typename... Args>
  vect(const value_type& a1, const value_type& a2,
       const Args&... args) noexcept {
    static_assert(Size > sizeof...(Args) + 1);
    auto it = q.begin();
    set_value_impl(it, a1, a2, args...);
    for (; it < q.end(); ++it) {
      *it = value_type{};
    }
  }

  explicit vect(const value_type& a1) noexcept {
    static_assert(Size == 1);
    q[0] = a1;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read/write.
   * \test PASSED
   */
  reference operator[](int i) noexcept { return q[i]; }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](int i) const noexcept { return q[i]; }

  /**
   * Sub-vector operator, accessor for read/write.
   * \test PASSED
   */
  vect_ref_view<self> operator[](const std::pair<int, int>& r) {
    return sub(*this)[r];
  }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](const std::pair<int, int>& r) const {
    return sub(*this)[r];
  }

  /**
   * Array indexing operator, accessor for read/write.
   * \test PASSED
   */
  reference operator()(int i) noexcept { return q[i]; }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()(int i) const noexcept { return q[i]; }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self& V) noexcept = default;

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  template <typename Vector>
  self& operator=(const Vector& V) {
    static_assert(is_readable_vector_v<Vector>);
    if (Size != V.size()) {
      throw std::range_error("Vector size mismatch.");
    }
    for (size_type i = 0; i < Size; ++i) {
      q[i] = V[i];
    }
    return *this;
  }

  template <typename U, unsigned int Index>
  self& operator=(const vect_component<U, Index>& V) noexcept {
    static_assert(Size > Index);
    for (pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i) {
      *p_i = value_type();
    }
    q[Index] = V.q;
  }

  template <typename U, unsigned int Index>
  friend self operator+(self lhs,
                        const vect_component<U, Index>& rhs) noexcept {
    static_assert(Size > Index);
    lhs.q[Index] += rhs.q;
    return lhs;
  }

  template <typename U, unsigned int Index>
  self& operator+=(const vect_component<U, Index>& rhs) noexcept {
    static_assert(Size > Index);
    q[Index] += rhs.q;
    return *this;
  }

  template <typename U, unsigned int Index>
  friend self operator+(const vect_component<U, Index>& lhs,
                        self rhs) noexcept {
    static_assert(Size > Index);
    rhs.q[Index] += lhs.q;
    return rhs;
  }

  template <typename U, unsigned int Index>
  friend self operator-(self lhs,
                        const vect_component<U, Index>& rhs) noexcept {
    static_assert(Size > Index);
    lhs.q[Index] -= rhs.q;
    return lhs;
  }

  template <typename U, unsigned int Index>
  self& operator-=(const vect_component<U, Index>& rhs) noexcept {
    static_assert(Size > Index);
    q[Index] -= rhs.q;
    return *this;
  }

  template <typename U, unsigned int Index>
  friend self operator-(const vect_component<U, Index>& lhs,
                        const self& rhs) noexcept {
    static_assert(Size > Index);
    self result = -rhs;
    result.q[Index] += lhs.q;
    return result;
  }
};

namespace rtti {

template <typename T, unsigned int Size>
struct get_type_id<vect<T, Size>> {
  static constexpr unsigned int ID = 0x00000011;
  static constexpr auto type_name = std::string_view{"vect"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const vect<T, Size>&;
  using load_type = vect<T, Size>&;
};

template <typename T, unsigned int Size, typename Tail>
struct get_type_info<vect<T, Size>, Tail> {
  using type =
      type_id<vect<T, Size>,
              typename get_type_info<
                  T, get_type_info<std::integral_constant<unsigned int, Size>,
                                   Tail>>::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<vect<T, Size>>::type_name, lsl_left_bracket,
      get_type_id<T>::type_name, lsl_comma,
      get_type_id<std::integral_constant<unsigned int, Size>>::type_name,
      lsl_right_bracket, get_type_name_tail<Tail>::value>;
};
}  // namespace rtti

namespace serialization {
template <typename T, unsigned int Size>
iarchive& operator>>(iarchive& in, vect<T, Size>& v) {
  for (unsigned int i = 0; i != Size; ++i) {
    in >> v[i];
  }
  return in;
}

template <typename T, unsigned int Size>
iarchive& operator&(iarchive& in,
                    const std::pair<std::string, vect<T, Size>&>& v) {
  for (unsigned int i = 0; i != Size; ++i) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    in& RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), v.second[i]);
  }
  return in;
}

template <typename T, unsigned int Size>
oarchive& operator<<(oarchive& out, const vect<T, Size>& v) {
  for (unsigned int i = 0; i != Size; ++i) {
    out << v[i];
  }
  return out;
}

template <typename T, unsigned int Size>
oarchive& operator&(oarchive& out,
                    const std::pair<std::string, const vect<T, Size>&>& v) {
  for (unsigned int i = 0; i != Size; ++i) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    out& RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), v.second[i]);
  }
  return out;
}
}  // namespace serialization

template <typename T, unsigned int Size>
struct is_readable_vector<vect<T, Size>> {
  static constexpr bool value = true;
  using type = is_readable_vector<vect<T, Size>>;
};

template <typename T, unsigned int Size>
struct is_writable_vector<vect<T, Size>> {
  static constexpr bool value = true;
  using type = is_writable_vector<vect<T, Size>>;
};

template <typename T, unsigned int Size>
struct is_resizable_vector<vect<T, Size>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<vect<T, Size>>;
};

template <typename T, unsigned int Size>
struct has_allocator_vector<vect<T, Size>> {
  static constexpr bool value = false;
  using type = has_allocator_vector<vect<T, Size>>;
};

/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template <typename T, typename... Args>
vect<T, sizeof...(Args) + 1> make_vect(const T& Q1,
                                       const Args&... args) noexcept {
  return {Q1, args...};
}

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T, Size> diff(const vect<T, Size>& v1, const vect<T, Size>& v2) noexcept {
  return v1 - v2;
}

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T, Size> add(const vect<T, Size>& v1, const vect<T, Size>& v2) noexcept {
  return v1 + v2;
}

/*******************************************************************************
                         Special Vector Products / Operators
*******************************************************************************/

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
T operator%(const vect<T, 2>& V1, const vect<T, 2>& V2) noexcept {
  return V1[0] * V2[1] - V1[1] * V2[0];
}

template <class T>
T operator%(const vect_component<T, 0>& V1,
            const vect_component<T, 0>& V2) noexcept {
  return T(0.0);
}

template <class T>
T operator%(const vect_component<T, 1>& V1,
            const vect_component<T, 1>& V2) noexcept {
  return T(0.0);
}

template <class T>
vect_component<T, 2> operator%(const vect_component<T, 0>& V1,
                               const vect_component<T, 1>& V2) noexcept {
  return vect_component<T, 2>(V1.q * V2.q);
}

template <class T>
vect_component<T, 2> operator%(const vect_component<T, 1>& V1,
                               const vect_component<T, 0>& V2) noexcept {
  return vect_component<T, 2>(-V1.q * V2.q);
}

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T, 2> operator%(const T& S, const vect<T, 2>& V) noexcept {
  vect<T, 2> result;
  result[0] = -V[1] * S;
  result[1] = V[0] * S;
  return result;
}

template <class T>
vect_component<T, 1> operator%(const T& S,
                               const vect_component<T, 0>& V) noexcept {
  return vect_component<T, 1>(V.q * S);
}

template <class T>
vect_component<T, 0> operator%(const T& S,
                               const vect_component<T, 1>& V) noexcept {
  return vect_component<T, 0>(-V.q * S);
}

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T, 2> operator%(const vect<T, 2>& V, const T& S) noexcept {
  vect<T, 2> result;
  result[0] = V[1] * S;
  result[1] = -V[0] * S;
  return result;
}

template <class T>
vect_component<T, 1> operator%(const vect_component<T, 0>& V,
                               const T& S) noexcept {
  return vect_component<T, 1>(-V.q * S);
}

template <class T>
vect_component<T, 0> operator%(const vect_component<T, 1>& V,
                               const T& S) noexcept {
  return vect_component<T, 0>(V.q * S);
}

/**
 * 3D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T, 3> operator%(const vect<T, 3>& V1, const vect<T, 3>& V2) noexcept {
  vect<T, 3> result;
  result[0] = V1[1] * V2[2] - V1[2] * V2[1];
  result[1] = V1[2] * V2[0] - V1[0] * V2[2];
  result[2] = V1[0] * V2[1] - V1[1] * V2[0];
  return result;
}

template <class T>
vect_component<T, 1> operator%(const vect_component<T, 0>& V1,
                               const vect_component<T, 2>& V2) noexcept {
  return vect_component<T, 1>(-V1.q * V2.q);
}

template <class T>
vect_component<T, 1> operator%(const vect_component<T, 2>& V1,
                               const vect_component<T, 0>& V2) noexcept {
  return vect_component<T, 1>(V1.q * V2.q);
}

template <class T>
vect_component<T, 0> operator%(const vect_component<T, 1>& V1,
                               const vect_component<T, 2>& V2) noexcept {
  return vect_component<T, 2>(V1.q * V2.q);
}

template <class T>
vect_component<T, 0> operator%(const vect_component<T, 2>& V1,
                               const vect_component<T, 1>& V2) noexcept {
  return vect_component<T, 2>(-V1.q * V2.q);
}

template <class T>
T operator%(const vect_component<T, 2>& V1,
            const vect_component<T, 2>& V2) noexcept {
  return T(0.0);
}

template <class T>
vect<T, 3> operator%(const vect<T, 3>& V1,
                     const vect_component<T, 0>& V2) noexcept {
  vect<T, 3> result;
  result[0] = T(0.0);
  result[1] = V1[2] * V2.q;
  result[2] = -V1[1] * V2.q;
  return result;
}

template <class T>
vect<T, 3> operator%(const vect<T, 3>& V1,
                     const vect_component<T, 1>& V2) noexcept {
  vect<T, 3> result;
  result[0] = -V1[2] * V2.q;
  result[1] = T(0.0);
  result[2] = V1[0] * V2.q;
  return result;
}

template <class T>
vect<T, 3> operator%(const vect<T, 3>& V1,
                     const vect_component<T, 2>& V2) noexcept {
  vect<T, 3> result;
  result[0] = V1[1] * V2.q;
  result[1] = -V1[0] * V2.q;
  result[2] = T(0.0);
  return result;
}

template <class T>
vect<T, 3> operator%(const vect_component<T, 0>& V1,
                     const vect<T, 3>& V2) noexcept {
  vect<T, 3> result;
  result[0] = T(0.0);
  result[1] = -V2[2] * V1.q;
  result[2] = V2[1] * V1.q;
  return result;
}

template <class T>
vect<T, 3> operator%(const vect_component<T, 1>& V1,
                     const vect<T, 3>& V2) noexcept {
  vect<T, 3> result;
  result[0] = V2[2] * V1.q;
  result[1] = T(0.0);
  result[2] = -V2[0] * V1.q;
  return result;
}

template <class T>
vect<T, 3> operator%(const vect_component<T, 2>& V1,
                     const vect<T, 3>& V2) noexcept {
  vect<T, 3> result;
  result[0] = -V2[1] * V1.q;
  result[1] = V2[0] * V1.q;
  result[2] = T(0.0);
  return result;
}

/**
 * This class implements a variable-size templated vector class which holds components of dimensional quantities.
 */
template <typename T, typename Allocator = std::allocator<T>>
class vect_n {
 public:
  using self = vect_n<T, Allocator>;

  using value_type = T;
  using reference = typename std::vector<T, Allocator>::reference;
  using const_reference = typename std::vector<T, Allocator>::const_reference;
  using pointer = typename std::vector<T, Allocator>::pointer;
  using const_pointer = typename std::vector<T, Allocator>::const_pointer;
  using allocator_type = Allocator;

  using iterator = typename std::vector<T, Allocator>::iterator;
  using const_iterator = typename std::vector<T, Allocator>::const_iterator;

  using size_type = typename std::vector<T, Allocator>::size_type;
  using difference_type = typename std::vector<T, Allocator>::difference_type;

  static constexpr std::size_t dimensions = 0;

  std::vector<T, Allocator> q; /**< Components of the vector. */

  /**
   * Returns the size of the vector.
   */
  size_type size() const noexcept { return q.size(); }
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const noexcept { return q.max_size(); }
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const noexcept { return q.capacity(); }
  /**
   * Resizes the vector.
   */
  void resize(size_type sz, T c = T()) { q.resize(sz, c); }
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return q.empty(); }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(size_type sz) { q.reserve(sz); }

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() noexcept { return q.begin(); }
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return q.begin(); }
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() noexcept { return q.end(); }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return q.end(); }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor.
   * \test PASSED
   */
  explicit vect_n(const allocator_type& aAlloc) : q(aAlloc) {}

  vect_n() : q(allocator_type()) {}

  /**
   * Constructor for a given size or dimension of vector.
   * \test PASSED
   */
  explicit vect_n(size_type aSize, const_reference Fill = value_type(),
                  const allocator_type& aAlloc = allocator_type())
      : q(aSize, Fill, aAlloc) {}

  /**
   * Constructor from an array of values of type "value_type".
   * \test PASSED
   */
  template <typename OtherAllocator>
  explicit vect_n(const std::vector<value_type, OtherAllocator>& Q,
                  const allocator_type& aAlloc = allocator_type())
      : q(Q.begin(), Q.end(), aAlloc) {}

  /**
   * Constructor from a forward iterator of values of type "value_type".
   * \test PASSED
   */
  template <typename InputIter>
  explicit vect_n(InputIter first, InputIter last,
                  const allocator_type& aAlloc = allocator_type())
      : q(first, last, aAlloc) {}

  /**
   * Standard Copy Constructor with standard semantics.
   * \test PASSED
   */
  vect_n(const self& V) : q(V.begin(), V.end(), V.get_allocator()) {}

  /**
   * Allocator-agnostic Copy Constructor with standard semantics.
   * \test PASSED
   */
  template <typename OtherAllocator>
  explicit vect_n(const vect_n<value_type, OtherAllocator>& V,
                  const allocator_type& aAlloc = allocator_type())
      : q(V.begin(), V.end(), aAlloc) {}

  /**
   * Constructor from a fixed-length vector.
   */
  template <unsigned int Size>
  explicit vect_n(const vect<value_type, Size>& V,
                  const allocator_type& aAlloc = allocator_type())
      : q(V.begin(), V.end(), aAlloc) {}

 private:
  static void set_value_impl(iterator& pval, const_reference a1) noexcept {
    *pval++ = a1;
  }

  template <typename... Args>
  static void set_value_impl(iterator& pval, const_reference a1,
                             const Args&... tail) noexcept {
    *pval++ = a1;
    set_value_impl(pval, tail...);
  }

 public:
  /**
   * Constructor for Size values.
   * \test PASSED
   */
  template <typename... Args>
  vect_n(const_reference a1, const_reference a2, const_reference a3,
         Args... args) noexcept
      : q(3 + sizeof...(Args)) {
    auto p_i = q.begin();
    set_value_impl(p_i, a1, a2, a3, args...);
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read/write. <
   * \test PASSED
   */
  reference operator[](size_type i) noexcept { return q[i]; }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](size_type i) const noexcept { return q[i]; }

  /**
   * Sub-vector operator, accessor for read/write.
   * \test PASSED
   */
  vect_ref_view<self> operator[](
      const std::pair<size_type, size_type>& r) noexcept {
    return sub(*this)[r];
  }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<size_type, size_type>& r) const noexcept {
    return sub(*this)[r];
  }

  /**
   * Array indexing operator, accessor for read/write. <
   * \test PASSED
   */
  reference operator()(size_type i) noexcept { return q[i]; }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()(size_type i) const noexcept { return q[i]; }

  /**
   * Returns the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return q.get_allocator(); }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  self& operator=(const self& V) {
    q.assign(V.q.begin(), V.q.end());
    return *this;
  }

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  template <typename Vector>
  self& operator=(const Vector& V) {
    static_assert(is_readable_vector_v<Vector>);
    q.resize(V.size());
    for (size_type i = 0; i < q.size(); ++i) {
      q[i] = V[i];
    }
    return *this;
  }
};

namespace rtti {

template <typename T, typename Allocator>
struct get_type_id<vect_n<T, Allocator>> {
  static constexpr unsigned int ID = 0x00000010;
  static constexpr auto type_name = std::string_view{"vect_n"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const vect_n<T, Allocator>&;
  using load_type = vect_n<T, Allocator>&;
};

template <typename T, typename Allocator, typename Tail>
struct get_type_info<vect_n<T, Allocator>, Tail> {
  using type =
      type_id<vect_n<T, Allocator>, typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<vect_n<T, Allocator>>::type_name,
                  lsl_left_bracket, get_type_id<T>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};
}  // namespace rtti

namespace serialization {
template <typename T, typename Allocator>
iarchive& operator>>(iarchive& in, vect_n<T, Allocator>& v) {
  in >> v.q;
  return in;
}

template <typename T, typename Allocator>
iarchive& operator&(iarchive& in,
                    const std::pair<std::string, vect_n<T, Allocator>&>& v) {
  in& RK_SERIAL_LOAD_WITH_ALIAS(v.first, v.second.q);
  return in;
}

template <typename T, typename Allocator>
oarchive& operator<<(oarchive& out, const vect_n<T, Allocator>& v) {
  out << v.q;
  return out;
}

template <typename T, typename Allocator>
oarchive& operator&(
    oarchive& out,
    const std::pair<std::string, const vect_n<T, Allocator>&>& v) {
  out& RK_SERIAL_SAVE_WITH_ALIAS(v.first, v.second.q);
  return out;
}
}  // namespace serialization

template <typename T, typename Allocator>
struct is_readable_vector<vect_n<T, Allocator>> {
  static constexpr bool value = true;
  using type = is_readable_vector<vect_n<T, Allocator>>;
};

template <typename T, typename Allocator>
struct is_writable_vector<vect_n<T, Allocator>> {
  static constexpr bool value = true;
  using type = is_writable_vector<vect_n<T, Allocator>>;
};

template <typename T, typename Allocator>
struct is_resizable_vector<vect_n<T, Allocator>> {
  static constexpr bool value = true;
  using type = is_resizable_vector<vect_n<T, Allocator>>;
};

template <typename T, typename Allocator>
struct has_allocator_vector<vect_n<T, Allocator>> {
  static constexpr bool value = true;
  using type = has_allocator_vector<vect_n<T, Allocator>>;
};

template <typename Vector>
struct vect_copy<vect_ref_view<Vector>> {
  using type = vect_n<vect_value_type_t<Vector>>;
};

template <typename Vector>
struct vect_copy<vect_const_ref_view<Vector>> {
  using type = vect_n<vect_value_type_t<Vector>>;
};

/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template <typename T, typename... Args>
vect_n<T> make_vect_n(const T& Q1, const T& Q2, const T& Q3,
                      const Args&... args) {
  return vect_n<T>(Q1, Q2, Q3, args...);
}

template <typename T, unsigned int Size>
template <typename U, typename Allocator>
vect<T, Size>::vect(const vect_n<U, Allocator>& V) noexcept {
  if (Size > V.size()) {
    std::copy(V.begin(), V.end(), q.begin());
    std::fill(q.begin() + V.size(), q.begin() + Size, T());
  } else {
    std::copy(V.begin(), V.begin() + Size, q.begin());
  }
}

/**
 * This class implements a variable-size templated vector class in which all vector-elements
 * have the same value.
 */
template <typename T, std::size_t Size = 0>
class vect_scalar {
 public:
  using self = vect_scalar<T>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using allocator_type = std::allocator<T>;

  using iterator = void;
  using const_iterator = vect_index_const_iter<self>;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t dimensions = Size;

 private:
  value_type q; /**< Components of the vector. */
 public:
  /**
   * Returns the size of the vector.
   */
  size_type size() const noexcept { return Size; }
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const noexcept { return Size; }
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const noexcept { return Size; }
  /**
   * Resizes the vector.
   */
  void resize(size_type sz, T c = T()) noexcept {}
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return true; }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(size_type sz) const noexcept {}

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return const_iterator(*this, Size); }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Constructor for a given value for the vector.
   * \test PASSED
   */
  explicit vect_scalar(const_reference aFill = value_type()) noexcept
      : q(aFill) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](size_type i) const noexcept { return q; }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<size_type, size_type>& r) const noexcept {
    return sub(*this)[r];
  }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()(size_type i) const noexcept { return q; }
};

template <typename T, std::size_t Size>
struct is_readable_vector<vect_scalar<T, Size>> {
  static constexpr bool value = true;
  using type = is_readable_vector<vect_scalar<T, Size>>;
};

template <typename T, std::size_t Size>
struct is_writable_vector<vect_scalar<T, Size>> {
  static constexpr bool value = false;
  using type = is_writable_vector<vect_scalar<T, Size>>;
};

template <typename T, std::size_t Size>
struct is_resizable_vector<vect_scalar<T, Size>> {
  static constexpr bool value = (Size == 0);
  using type = is_resizable_vector<vect_scalar<T, Size>>;
};

template <typename T, std::size_t Size>
struct has_allocator_vector<vect_scalar<T, Size>> {
  static constexpr bool value = false;
  using type = has_allocator_vector<vect_scalar<T, Size>>;
};

template <typename T, std::size_t Size>
struct vect_copy<vect_scalar<T, Size>> {
  using type = vect_n<T>;
};

/**
 * This class implements a variable-size templated vector class in which all vector-elements
 * have the same value.
 */
template <typename T>
class vect_scalar<T, 0> {
 public:
  using self = vect_scalar<T, 0>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using allocator_type = std::allocator<T>;

  using iterator = void;
  using const_iterator = vect_index_const_iter<self>;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t dimensions = 0;

 private:
  value_type q; /**< Components of the vector. */
  size_type count;

 public:
  /**
   * Returns the size of the vector.
   */
  size_type size() const noexcept { return count; }
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const noexcept {
    return std::numeric_limits<size_type>::max();
  }
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const noexcept {
    return std::numeric_limits<size_type>::max();
  }
  /**
   * Resizes the vector.
   */
  void resize(size_type sz, T c = T()) noexcept {
    count = sz;
    q = c;
  }
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return (count == 0); }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(size_type sz) const noexcept {}

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return const_iterator(*this, count); }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor.
   * \test PASSED
   */
  vect_scalar() noexcept : q(), count(0) {}

  /**
   * Constructor for a given size or dimension of vector.
   * \test PASSED
   */
  explicit vect_scalar(size_type aSize,
                       const_reference aFill = value_type()) noexcept
      : q(aFill), count(aSize) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](size_type i) const noexcept { return q; }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<size_type, size_type>& r) const noexcept {
    return sub(*this)[r];
  }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()(size_type i) const noexcept { return q; }
};

/*******************************************************************************
                         Basic Functions
*******************************************************************************/

/**
 * Square magnitude of the vector.
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_value_type_t<Vector>>
norm_2_sqr(const Vector& v) noexcept {
  vect_value_type_t<Vector> sum(0.0);
  for (int i = 0; i < v.size(); ++i) {
    sum += v[i] * v[i];
  }
  return sum;
}

/**
 * Magnitude of the vector.
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_value_type_t<Vector>>
norm_2(const Vector& v) noexcept {
  using std::sqrt;
  return sqrt(norm_2_sqr(v));
}

/**
 * Infinite norm of the vector.
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_value_type_t<Vector>>
norm_inf(const Vector& v) noexcept {
  using std::abs;
  vect_value_type_t<Vector> result(0.0);
  for (int i = 0; i < v.size(); ++i) {
    if (result < abs(v[i])) {
      result = abs(v[i]);
    }
  }
  return result;
}

/**
 * Square magnitude of the vector.
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_value_type_t<Vector>>
norm_1(const Vector& v) noexcept {
  using std::abs;
  vect_value_type_t<Vector> sum(0.0);
  for (int i = 0; i < v.size(); ++i) {
    sum += abs(v[i]);
  }
  return sum;
}

/**
 * Unit vector in the same direction.
 */
template <typename Vector>
auto unit(const Vector& v) {
  return v / norm_2(v);
}

/**
 * Checks if two vectors are colinear.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
bool colinear(const Vector1& v1, const Vector2& v2) noexcept {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;
  ValueType tmp_mag2 = norm_2(v2);
  ValueType tmp_mag1 = norm_2(v1);
  ValueType tmp_comb = norm_2(v1 + v2);
  return (((tmp_mag1 + tmp_mag2) *
               (ValueType(1.0) -
                ValueType(10.0) * std::numeric_limits<ValueType>::epsilon()) <
           tmp_comb) ||
          (abs(tmp_mag1 - tmp_mag2) *
               (ValueType(1.0) + std::numeric_limits<ValueType>::epsilon()) >
           tmp_comb));
}

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

/**
 * Standard add-and-store operator.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_writable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 Vector1&>
operator+=(Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] = v1[i] + v2[i];
  }
  return v1;
}

/**
 * Standard sub-and-store operator.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_writable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 Vector1&>
operator-=(Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] = v1[i] - v2[i];
  }
  return v1;
}

/**
 * Scalar multiply-and-store operator for gain.
 * \test PASSED
 */
template <typename T, typename Vector>
std::enable_if_t<is_writable_vector_v<Vector> && !is_readable_vector_v<T>,
                 Vector&>
operator*=(Vector& v, const T& S) noexcept {
  using ValueType = vect_value_type_t<Vector>;
  for (int i = 0; i < v.size(); ++i) {
    v[i] *= ValueType(S);
  }
  return v;
}

/**
 * Scalar divide-and-store operator for gain.
 * \test PASSED
 */
template <typename T, typename Vector>
std::enable_if_t<is_writable_vector_v<Vector> && !is_readable_vector_v<T>,
                 Vector&>
operator/=(Vector& v, const T& S) noexcept {
  using ValueType = vect_value_type_t<Vector>;
  for (int i = 0; i < v.size(); ++i) {
    v[i] /= ValueType(S);
  }
  return v;
}

/**
 * Add two vectors.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_writable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 vect_copy_t<Vector1>>
operator+(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_copy_t<Vector1> result;
  result = v1;
  result += v2;
  return result;
}

/**
 * Invert the vector.
 * \test PASSED
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_copy_t<Vector>> operator-(
    const Vector& v) {
  vect_copy_t<Vector> result;
  result = v;
  for (int i = 0; i < v.size(); ++i) {
    result[i] = -result[i];
  }
  return result;
}

/**
 * Sub two vectors.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 vect_copy_t<Vector1>>
operator-(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_copy_t<Vector1> result;
  result = v1;
  result -= v2;
  return result;
}

/**
 * Dot Product.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 vect_value_type_t<Vector1>>
operator*(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_value_type_t<Vector1> result(0);
  for (int i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

/**
 * Scalar-vector product.
 * \test PASSED
 */
template <typename T, typename Vector>
std::enable_if_t<is_readable_vector_v<Vector> && !is_readable_vector_v<T> &&
                     !is_readable_matrix_v<T>,
                 vect_copy_t<Vector>>
operator*(const Vector& v, const T& S) {
  vect_copy_t<Vector> result;
  result = v;
  result *= S;
  return result;
}

/**
 * Scalar-vector product.
 * \test PASSED
 */
template <typename T, typename Vector>
std::enable_if_t<is_readable_vector_v<Vector> && !is_readable_vector_v<T> &&
                     !is_readable_matrix_v<T>,
                 vect_copy_t<Vector>>
operator*(const T& S, const Vector& v) {
  vect_copy_t<Vector> result;
  result = v;
  result *= S;
  return result;
}

/**
 * Scalar-vector division.
 * \test PASSED
 */
template <typename T, typename Vector>
std::enable_if_t<is_readable_vector_v<Vector> && !is_readable_vector_v<T> &&
                     !is_readable_matrix_v<T>,
                 vect_copy_t<Vector>>
operator/(const Vector& v, const T& S) {
  vect_copy_t<Vector> result;
  result = v;
  result /= S;
  return result;
}

/**
 * Element-wise product of two vectors.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_writable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 vect_copy_t<Vector1>>
elem_product(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_copy_t<Vector1> result;
  result = v1;
  for (int i = 0; i < result.size(); ++i) {
    result[i] *= v2[i];
  }
  return result;
}

/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T, Allocator> diff(const vect_n<T, Allocator>& v1,
                          const vect_n<T, Allocator>& v2) {
  return v1 - v2;
}

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T, Allocator> add(const vect_n<T, Allocator>& v1,
                         const vect_n<T, Allocator>& v2) {
  return v1 + v2;
}

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator, component-wise.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator==(const Vector1& v1, const Vector2& v2) noexcept {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (int i = 0; i < v1.size(); ++i) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

/**
 * Inequality Comparison operator, component-wise.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator!=(const Vector1& v1, const Vector2& v2) noexcept {
  if (v1.size() != v2.size()) {
    return true;
  }
  for (int i = 0; i < v1.size(); ++i) {
    if (v1[i] != v2[i]) {
      return true;
    }
  }
  return false;
}

/**
 * Greater-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator>(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) > norm_2_sqr(v2));
}

/**
 * Smaller-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator<(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) < norm_2_sqr(v2));
}

/**
 * Greater-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator>=(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) >= norm_2_sqr(v2));
}

/**
 * Smaller-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename Vector1, typename Vector2>
std::enable_if_t<is_readable_vector_v<Vector1> && is_readable_vector_v<Vector2>,
                 bool>
operator<=(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) <= norm_2_sqr(v2));
}

/**
 * Prints a variable-size vector to an output stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, std::ostream&> operator<<(
    std::ostream& out_stream, const Vector& V) {
  out_stream << "(";
  if (V.size() > 0) {
    out_stream << V[0];
  }
  for (int i = 1; i < V.size(); ++i) {
    out_stream << "; " << V[i];
  }
  return out_stream << ")";
}

/**
 * Reads a variable-size vector to an input stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template <typename T>
std::istream& operator>>(std::istream& in_stream, vect_n<T>& V) {
  std::string tmp_str;
  std::getline(in_stream, tmp_str, '(');  // skip to opening bracket.
  std::getline(in_stream, tmp_str, ')');  // read to closing bracket.
  int sz = std::count(tmp_str.begin(), tmp_str.end(), ';') + 1;
  std::stringstream ss(tmp_str);
  V.resize(sz);
  std::string tmp2;
  for (int i = 0; ss >> V[i]; ++i) {
    std::getline(ss, tmp2, ';');
  }
  return in_stream;
}

/**
 * Reads a variable-size vector to an input stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template <typename T, unsigned int Size>
std::istream& operator>>(std::istream& in_stream, vect<T, Size>& V) {
  std::string tmp_str;
  std::getline(in_stream, tmp_str, '(');  // skip to opening bracket.
  std::getline(in_stream, tmp_str, ')');  // read to closing bracket.
  std::stringstream ss(tmp_str);
  std::string tmp2;
  for (int i = 0; (i < Size) && (ss >> V[i]); ++i) {
    std::getline(ss, tmp2, ';');
  }
  return in_stream;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_ALG_H_
