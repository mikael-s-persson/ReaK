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
  using size_type = int;

  value_type q;

  explicit vect_component(const_reference Q = 0) noexcept : q(Q) {}

  value_type operator[](int i) const noexcept {
    return (i == Index ? q : value_type());
  }

  value_type operator()(int i) const noexcept {
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

  template <typename OtherScalar, unsigned int OtherIndex>
  friend auto operator%(
      const self& V1,
      const vect_component<OtherScalar, OtherIndex>& V2) noexcept {
    if constexpr (Index == OtherIndex) {
      return value_type{0.0};
    } else if constexpr (Index == 0) {
      if constexpr (OtherIndex == 1) {
        auto val = V1.q * V2.q;
        return vect_component<decltype(val), 2>(val);
      } else if constexpr (OtherIndex == 2) {
        auto val = -V1.q * V2.q;
        return vect_component<decltype(val), 1>(val);
      } else {
        static_assert(
            OtherIndex <= 2,
            "Cross product not supported for greater than 3 dimensions!");
        return value_type{0.0};
      }
    } else if constexpr (Index == 1) {
      if constexpr (OtherIndex == 0) {
        auto val = -V1.q * V2.q;
        return vect_component<decltype(val), 2>(val);
      } else if constexpr (OtherIndex == 2) {
        auto val = V1.q * V2.q;
        return vect_component<decltype(val), 0>(val);
      } else {
        static_assert(
            OtherIndex <= 2,
            "Cross product not supported for greater than 3 dimensions!");
        return value_type{0.0};
      }
    } else if constexpr (Index == 2) {
      if constexpr (OtherIndex == 0) {
        auto val = V1.q * V2.q;
        return vect_component<decltype(val), 1>(val);
      } else if constexpr (OtherIndex == 1) {
        auto val = -V1.q * V2.q;
        return vect_component<decltype(val), 0>(val);
      } else {
        static_assert(
            OtherIndex <= 2,
            "Cross product not supported for greater than 3 dimensions!");
        return value_type{0.0};
      }
    } else {
      static_assert(
          Index <= 2,
          "Cross product not supported for greater than 3 dimensions!");
      return value_type{0.0};
    }
  }

  friend auto operator%(const_reference S, const self& V) noexcept {
    if constexpr (Index == 0) {
      return vect_component<T, 1>(V.q * S);
    } else if constexpr (Index == 1) {
      return vect_component<T, 0>(-V.q * S);
    } else {
      static_assert(Index <= 1,
                    "Cross product with scalar not supported for greater than "
                    "2 dimensions!");
      return value_type{0.0};
    }
  }

  friend auto operator%(const self& V, const_reference S) noexcept {
    if constexpr (Index == 0) {
      return vect_component<value_type, 1>(-V.q * S);
    } else if constexpr (Index == 1) {
      return vect_component<value_type, 0>(V.q * S);
    } else {
      static_assert(Index <= 1,
                    "Cross product with scalar not supported for greater than "
                    "2 dimensions!");
      return value_type{0.0};
    }
  }
};

static const vect_component<double, 0> vect_i = vect_component<double, 0>(1.0);
static const vect_component<double, 1> vect_j = vect_component<double, 1>(1.0);
static const vect_component<double, 2> vect_k = vect_component<double, 2>(1.0);

/**
 * This class implements a fixed-size templated vector class which holds components of primitive type.
 */
template <typename T, unsigned int Size>
class vect {
 public:
  using self = vect<T, Size>;

  static constexpr bool is_dynamic_size = (Size == 0);
  using container_type =
      std::conditional_t<is_dynamic_size, std::vector<T>, std::array<T, Size>>;

  using value_type = T;
  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;

  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  using size_type = int;
  using difference_type = int;

  static constexpr std::size_t dimensions = Size;

  /// Components of the vector.
  container_type q = {};

  int size() const noexcept { return q.size(); }
  int max_size() const noexcept {
    if constexpr (!is_dynamic_size) {
      return Size;
    } else {
      return q.max_size();
    }
  }
  int capacity() const noexcept {
    if constexpr (!is_dynamic_size) {
      return Size;
    } else {
      return q.capacity();
    }
  }
  void resize(int sz, T c = T()) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      q.resize(sz, c);
    }
  }
  bool empty() const noexcept { return q.empty(); }
  void reserve(int sz) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      q.reserve(sz);
    }
  }

  iterator begin() noexcept { return q.begin(); }
  const_iterator begin() const noexcept { return q.begin(); }
  iterator end() noexcept { return q.end(); }
  const_iterator end() const noexcept { return q.end(); }

  /// Default constructor: sets all to zero.
  vect() noexcept {
    for (auto& v : q) {
      v = value_type{};
    }
  }

  /// Constructor for a given size or dimension of vector, or single value,
  /// or pointer to statically sized array of values.
  template <typename U>
  explicit vect(const U& value_or_size) noexcept {
    if constexpr (is_dynamic_size) {
      if constexpr (std::is_integral_v<U>) {
        q.resize(value_or_size, value_type{});
      } else if constexpr (ReadableVector<U>) {
        q.resize(value_or_size.size());
        for (int i = 0; i < q.size(); ++i) {
          q[i] = value_or_size[i];
        }
      } else {
        q.resize(1);
        q[0] = value_or_size;
      }
    } else if constexpr (std::is_pointer_v<U>) {
      U p_val = value_or_size;
      for (auto& v : q) {
        v = *p_val++;
      }
    } else if constexpr (ReadableVector<U>) {
      if (Size != value_or_size.size()) {
        throw std::range_error("Vector size mismatch.");
      }
      for (int i = 0; i < q.size(); ++i) {
        q[i] = value_or_size[i];
      }
    } else {
      static_assert(Size == 1);
      q[0] = value_or_size;
    }
  }

  /// Constructor from an array of values of type "value_type".
  explicit vect(const std::vector<value_type>& Q) : q(Q.begin(), Q.end()) {}

  /// Constructor from a forward iterator of values of type "value_type".
  template <typename InputIter>
  explicit vect(
      InputIter first,
      std::enable_if_t<std::is_convertible_v<decltype(*first), value_type>,
                       InputIter>
          last)
      : q(first, last) {}

  vect(const self&) noexcept(!is_dynamic_size) = default;
  vect(self&&) noexcept = default;

  template <typename U, unsigned int OtherSize>
  vect(const vect<U, OtherSize>& rhs) noexcept {  // NOLINT
    if constexpr (is_dynamic_size) {
      q.resize(rhs.size());
      std::copy(rhs.begin(), rhs.end(), q.begin());
    } else {
      static_assert(Size >= OtherSize);
      std::copy(rhs.begin(), rhs.end(), q.begin());
      std::fill(q.begin() + rhs.size(), q.begin() + Size, value_type{});
    }
  }

  template <typename U, unsigned int Index>
  explicit vect(const vect_component<U, Index>& rhs) noexcept : vect() {
    if constexpr (is_dynamic_size) {
      q.resize(Index, value_type{});
    } else {
      static_assert(Size > Index);
    }
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
  /// Constructor for fixed number of values.
  template <typename Arg1, typename Arg2, typename... Args>
  vect(const Arg1& a1, const Arg2& a2,
       const Args&... args) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if constexpr (std::is_integral_v<Arg1> && (sizeof...(Args) == 0)) {
        q.resize(a1, a2);
        return;
      } else {
        q.resize(sizeof...(Args) + 2);
      }
    } else {
      static_assert(Size > sizeof...(Args) + 1);
    }
    auto it = q.begin();
    set_value_impl(it, a1, a2, args...);
    for (; it < q.end(); ++it) {
      *it = value_type{};
    }
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  reference operator[](int i) noexcept { return q[i]; }
  const_reference operator[](int i) const noexcept { return q[i]; }

  reference operator()(int i) noexcept { return q[i]; }
  const_reference operator()(int i) const noexcept { return q[i]; }

  /*******************************************************************************
                         Assignment Operators
  *******************************************************************************/

  self& operator=(const self& V) noexcept(!is_dynamic_size) = default;
  self& operator=(self&& V) noexcept = default;

  template <ReadableVector Vector>
  self& operator=(const Vector& rhs) {
    if constexpr (is_dynamic_size) {
      q.resize(rhs.size());
    } else {
      if (Size != rhs.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    }
    for (int i = 0; i < q.size(); ++i) {
      q[i] = rhs[i];
    }
    return *this;
  }

  template <typename U, unsigned int Index>
  self& operator=(const vect_component<U, Index>& rhs) noexcept(
      !is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      q.resize(std::max<int>(q.size(), Index));
    } else {
      static_assert(Size > Index);
    }
    std::fill(q.begin(), q.end(), value_type{});
    q[Index] = rhs.q;
  }

  template <typename U, unsigned int Index>
  friend self operator+(self lhs, const vect_component<U, Index>& rhs) noexcept(
      !is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (Index >= lhs.q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    lhs.q[Index] += rhs.q;
    return lhs;
  }

  template <typename U, unsigned int Index>
  self& operator+=(const vect_component<U, Index>& rhs) noexcept {
    if constexpr (is_dynamic_size) {
      if (Index >= q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    q[Index] += rhs.q;
    return *this;
  }

  template <typename U, unsigned int Index>
  friend self operator+(const vect_component<U, Index>& lhs,
                        self rhs) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (Index >= rhs.q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    rhs.q[Index] += lhs.q;
    return rhs;
  }

  template <typename U, unsigned int Index>
  friend self operator-(self lhs, const vect_component<U, Index>& rhs) noexcept(
      !is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (Index >= lhs.q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    lhs.q[Index] -= rhs.q;
    return lhs;
  }

  template <typename U, unsigned int Index>
  self& operator-=(const vect_component<U, Index>& rhs) noexcept {
    if constexpr (is_dynamic_size) {
      if (Index >= q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    q[Index] -= rhs.q;
    return *this;
  }

  template <typename U, unsigned int Index>
  friend self operator-(const vect_component<U, Index>& lhs,
                        const self& rhs) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (Index >= rhs.q.size()) {
        throw std::range_error("Vector size mismatch.");
      }
    } else {
      static_assert(Size > Index);
    }
    self result = -rhs;
    result.q[Index] += lhs.q;
    return result;
  }

  /*******************************************************************************
                           Special Vector Products / Operators
  *******************************************************************************/

  /// Cross-Product.
  template <class U, unsigned int OtherSize>
  friend auto operator%(
      const self& lhs,
      const vect<U, OtherSize>& rhs) noexcept(!is_dynamic_size) {
    if constexpr (Size == 2) {
      return lhs[0] * rhs[1] - lhs[1] * rhs[0];
    } else {
      if constexpr (is_dynamic_size) {
        if (lhs.size() != 3) {
          throw std::range_error(
              "Cross-product only supported for dynamic vector of size 3.");
        }
        if (rhs.size() != 3) {
          throw std::range_error(
              "Cross-product only supported for dynamic vector of size 3.");
        }
      } else {
        static_assert(
            Size == 3,
            "Cross-product only supported for static vector of size 2 or 3.");
        static_assert(
            OtherSize == 3,
            "Cross-product only supported for static vector of size 2 or 3.");
      }
      return self(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                  lhs[2] * rhs[0] - lhs[0] * rhs[2],
                  lhs[0] * rhs[1] - lhs[1] * rhs[0]);
    }
  }

  /// Cross-Product.
  friend self operator%(const_reference lhs,
                        const self& rhs) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (rhs.size() != 2) {
        throw std::range_error(
            "Cross product with scalar not supported for other than 2 "
            "dimensions!");
      }
    } else {
      static_assert(Size == 2,
                    "Cross product with scalar not supported for other than 2 "
                    "dimensions!");
    }
    return self(-rhs[1] * lhs, rhs[0] * lhs);
  }

  /// Cross-Product.
  friend self operator%(const self& lhs,
                        const_reference rhs) noexcept(!is_dynamic_size) {
    if constexpr (is_dynamic_size) {
      if (rhs.size() != 2) {
        throw std::range_error(
            "Cross product with scalar not supported for other than 2 "
            "dimensions!");
      }
    } else {
      static_assert(Size == 2,
                    "Cross product with scalar not supported for other than 2 "
                    "dimensions!");
    }
    return self(lhs[1] * rhs, -lhs[0] * rhs);
  }

  /// Cross-Product.
  template <typename U, unsigned int Index>
  friend auto operator%(
      const self& lhs,
      const vect_component<U, Index>& rhs) noexcept(!is_dynamic_size) {
    if constexpr (Size == 2) {
      if constexpr (Index == 0) {
        return -lhs[1] * rhs.q;
      } else if constexpr (Index == 1) {
        return lhs[0] * rhs.q;
      } else {
        static_assert(
            Index <= 1,
            "Cross-product only supported for vector components below 2.");
        return value_type{};
      }
    } else {
      if constexpr (is_dynamic_size) {
        if (lhs.size() != 3) {
          throw std::range_error(
              "Cross-product only supported for dynamic vector of size 3.");
        }
      } else {
        static_assert(
            Size == 3,
            "Cross-product only supported for static vector of size 2 or 3.");
      }
      if constexpr (Index == 0) {
        return self(0.0, lhs[2] * rhs.q, -lhs[1] * rhs.q);
      } else if constexpr (Index == 1) {
        return self(-lhs[2] * rhs.q, 0.0, lhs[0] * rhs.q);
      } else if constexpr (Index == 2) {
        return self(lhs[1] * rhs.q, -lhs[0] * rhs.q, 0.0);
      } else {
        static_assert(
            Index <= 2,
            "Cross-product only supported for vector components below 3.");
        return self{};
      }
    }
  }

  /// Cross-Product.
  template <typename U, unsigned int Index>
  friend auto operator%(const vect_component<U, Index>& lhs,
                        const self& rhs) noexcept(!is_dynamic_size) {
    if constexpr (Size == 2) {
      if constexpr (Index == 0) {
        return rhs[1] * lhs.q;
      } else if constexpr (Index == 1) {
        return -rhs[0] * lhs.q;
      } else {
        static_assert(
            Index <= 1,
            "Cross-product only supported for vector components below 2.");
        return value_type{};
      }
    } else {
      if constexpr (is_dynamic_size) {
        if (rhs.size() != 3) {
          throw std::range_error(
              "Cross-product only supported for dynamic vector of size 3.");
        }
      } else {
        static_assert(
            Size == 3,
            "Cross-product only supported for static vector of size 2 or 3.");
      }
      if constexpr (Index == 0) {
        return self(0.0, -rhs[2] * lhs.q, rhs[1] * lhs.q);
      } else if constexpr (Index == 1) {
        return self(rhs[2] * lhs.q, 0.0, -rhs[0] * lhs.q);
      } else if constexpr (Index == 2) {
        return self(-rhs[1] * lhs.q, rhs[0] * lhs.q, 0.0);
      } else {
        static_assert(
            Index <= 2,
            "Cross-product only supported for vector components below 3.");
        return self{};
      }
    }
  }

  friend auto diff(const self& v1, const self& v2) noexcept(!is_dynamic_size) {
    return v1 - v2;
  }

  friend auto add(const self& v1, const self& v2) noexcept(!is_dynamic_size) {
    return v1 + v2;
  }
};

namespace rtti {

template <typename T, unsigned int Size>
struct get_type_id<vect<T, Size>> {
  static constexpr unsigned int id = 0x00000011;
  static constexpr auto type_name = std::string_view{"vect"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const vect<T, Size>&;
  using load_type = vect<T, Size>&;
};

template <typename T, unsigned int Size, typename Tail>
struct get_type_info<vect<T, Size>, Tail> {
  using type = so_type_details::type_id<
      vect<T, Size>,
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
  if constexpr (Size == 0) {
    int sz = 0;
    in >> sz;
    v.resize(sz);
  }
  for (int i = 0; i != v.size(); ++i) {
    in >> v[i];
  }
  return in;
}

template <typename T, unsigned int Size>
iarchive& operator&(iarchive& in,
                    const std::pair<std::string, vect<T, Size>&>& v) {
  if constexpr (Size == 0) {
    int sz = 0;
    std::stringstream s_stream;
    s_stream << v.first << "_count";
    in& RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), sz);
    v.second.resize(sz);
  }
  for (int i = 0; i != v.second.size(); ++i) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    in& RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), v.second[i]);
  }
  return in;
}

template <typename T, unsigned int Size>
oarchive& operator<<(oarchive& out, const vect<T, Size>& v) {
  if constexpr (Size == 0) {
    out << v.size();
  }
  for (int i = 0; i != v.size(); ++i) {
    out << v[i];
  }
  return out;
}

template <typename T, unsigned int Size>
oarchive& operator&(oarchive& out,
                    const std::pair<std::string, const vect<T, Size>&>& v) {
  if constexpr (Size == 0) {
    std::stringstream s_stream;
    s_stream << v.first << "_count";
    out& RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), v.second.size());
  }
  for (int i = 0; i != v.second.size(); ++i) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    out& RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), v.second[i]);
  }
  return out;
}
}  // namespace serialization

/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template <typename T, typename... Args>
vect<T, sizeof...(Args) + 1> make_vect(const T& Q1,
                                       const Args&... args) noexcept {
  return {Q1, args...};
}

// This class implements a variable-size templated vector class which holds components of dimensional quantities.
template <typename T>
using vect_n = vect<T, 0>;

template <typename T, typename... Args>
vect<T, 0> make_vect_n(const T& Q1, const T& Q2, const T& Q3,
                       const Args&... args) {
  return vect<T, 0>(Q1, Q2, Q3, args...);
}

template <typename Vector>
struct vect_copy<vect_ref_view<Vector>> {
  using type = vect_n<vect_value_type_t<Vector>>;
};

template <typename Vector>
struct vect_copy<vect_const_ref_view<Vector>> {
  using type = vect_n<vect_value_type_t<Vector>>;
};

/**
 * This class implements a variable-size templated vector class in which all vector-elements
 * have the same value.
 */
template <typename T, int Size = 0>
class vect_scalar {
 public:
  using self = vect_scalar<T>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using iterator = void;
  using const_iterator = vect_index_const_iter<self>;

  using size_type = int;
  using difference_type = int;

  static constexpr std::size_t dimensions = Size;

  static constexpr bool is_dynamic_size = (Size == 0);

 private:
  struct data_alone {
    value_type q;
  };
  struct data_with_count {
    value_type q;
    int count;
  };
  using storage_type =
      std::conditional_t<is_dynamic_size, data_with_count, data_alone>;

  storage_type data;

 public:
  /// Returns the size of the vector.
  int size() const noexcept {
    if constexpr (is_dynamic_size) {
      return data.count;
    } else {
      return Size;
    }
  }
  /// Returns the max-size of the vector.
  int max_size() const noexcept {
    if constexpr (is_dynamic_size) {
      return std::numeric_limits<int>::max();
    } else {
      return Size;
    }
  }
  /// Returns the capacity of the vector.
  int capacity() const noexcept {
    if constexpr (is_dynamic_size) {
      return std::numeric_limits<int>::max();
    } else {
      return Size;
    }
  }
  /// Resizes the vector.
  void resize(int sz, T c = T()) noexcept {
    if constexpr (is_dynamic_size) {
      data.count = sz;
    }
  }
  /// Checks if the vector is empty.
  bool empty() const noexcept {
    if constexpr (is_dynamic_size) {
      return data.count == 0;
    }
    return false;
  }
  /// Reserve a capacity for the vector.
  void reserve(int sz) const noexcept {}

  /// Returns a const-iterator to the first element of the vector.
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// Returns a const-iterator to the one-past-last element of the vector.
  const_iterator end() const noexcept {
    if constexpr (is_dynamic_size) {
      return const_iterator(*this, data.count);
    } else {
      return const_iterator(*this, Size);
    }
  }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Constructor for a given value for the vector.
  explicit vect_scalar(const_reference aFill) noexcept { data.q = aFill; }

  explicit vect_scalar(int aSize, const_reference aFill =
                                      value_type()) noexcept(is_dynamic_size) {
    data.q = aFill;
    if constexpr (is_dynamic_size) {
      data.count = aSize;
    } else {
      if (aSize != Size) {
        throw std::range_error("Vector size mismatch!");
      }
    }
  }

  vect_scalar() : vect_scalar(value_type{}) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const noexcept { return data.q; }

  /// Array indexing operator, accessor for read only.
  const_reference operator()(int i) const noexcept { return data.q; }
};

template <typename T, int Size>
struct vect_copy<vect_scalar<T, Size>> {
  using type = vect_n<T>;
};

/*******************************************************************************
                         Basic Functions
*******************************************************************************/

/// Square magnitude of the vector.
template <ReadableVector Vector>
auto norm_2_sqr(const Vector& v) noexcept {
  vect_value_type_t<Vector> sum(0.0);
  for (int i = 0; i < v.size(); ++i) {
    sum += v[i] * v[i];
  }
  return sum;
}

/// Magnitude of the vector.
template <ReadableVector Vector>
auto norm_2(const Vector& v) noexcept {
  using std::sqrt;
  return sqrt(norm_2_sqr(v));
}

/// Magnitude of a scalar.
template <typename Scalar>
auto norm_2(const Scalar& v) noexcept {
  using std::abs;
  return abs(v);
}

/// Infinite norm of the vector.
template <ReadableVector Vector>
auto norm_inf(const Vector& v) noexcept {
  using std::abs;
  vect_value_type_t<Vector> result(0.0);
  for (int i = 0; i < v.size(); ++i) {
    if (result < abs(v[i])) {
      result = abs(v[i]);
    }
  }
  return result;
}

/// Square magnitude of the vector.
template <ReadableVector Vector>
auto norm_1(const Vector& v) noexcept {
  using std::abs;
  vect_value_type_t<Vector> sum(0.0);
  for (int i = 0; i < v.size(); ++i) {
    sum += abs(v[i]);
  }
  return sum;
}

/// Unit vector in the same direction.
template <typename Vector>
auto unit(const Vector& v) {
  return v / norm_2(v);
}

/// Checks if two vectors are colinear.
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

/// Standard add-and-store operator.
template <WritableVector Vector1, ReadableVector Vector2>
Vector1& operator+=(Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] = v1[i] + v2[i];
  }
  return v1;
}

/// Standard sub-and-store operator.
template <WritableVector Vector1, ReadableVector Vector2>
Vector1& operator-=(Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] = v1[i] - v2[i];
  }
  return v1;
}

/// Scalar multiply-and-store operator for gain.
template <WritableVector Vector,
          std::convertible_to<vect_value_type_t<Vector>> T>
Vector& operator*=(Vector& v, const T& S) noexcept {
  using ValueType = vect_value_type_t<Vector>;
  for (int i = 0; i < v.size(); ++i) {
    v[i] *= ValueType(S);
  }
  return v;
}

/// Scalar divide-and-store operator for gain.
template <WritableVector Vector,
          std::convertible_to<vect_value_type_t<Vector>> T>
Vector& operator/=(Vector& v, const T& S) noexcept {
  using ValueType = vect_value_type_t<Vector>;
  for (int i = 0; i < v.size(); ++i) {
    v[i] /= ValueType(S);
  }
  return v;
}

/// Add two vectors.
template <ReadableVector Vector1, ReadableVector Vector2>
vect_copy_t<Vector1> operator+(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_copy_t<Vector1> result;
  result = v1;
  result += v2;
  return result;
}

/// Invert the vector.
template <ReadableVector Vector>
vect_copy_t<Vector> operator-(const Vector& v) {
  vect_copy_t<Vector> result;
  result = v;
  for (int i = 0; i < v.size(); ++i) {
    result[i] = -result[i];
  }
  return result;
}

/// Sub two vectors.
template <ReadableVector Vector1, ReadableVector Vector2>
vect_copy_t<Vector1> operator-(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_copy_t<Vector1> result;
  result = v1;
  result -= v2;
  return result;
}

/// Dot Product.
template <ReadableVector Vector1, ReadableVector Vector2>
vect_value_type_t<Vector1> operator*(const Vector1& v1, const Vector2& v2) {
  if (v1.size() != v2.size()) {
    throw std::range_error("Vector size mismatch.");
  }
  vect_value_type_t<Vector1> result(0);
  for (int i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

/// Scalar-vector product.
template <ReadableVector Vector,
          std::convertible_to<vect_value_type_t<Vector>> T>
vect_copy_t<Vector> operator*(const Vector& v, const T& S) {
  vect_copy_t<Vector> result;
  result = v;
  result *= S;
  return result;
}

/// Scalar-vector product.
template <ReadableVector Vector,
          std::convertible_to<vect_value_type_t<Vector>> T>
vect_copy_t<Vector> operator*(const T& S, const Vector& v) {
  vect_copy_t<Vector> result;
  result = v;
  result *= S;
  return result;
}

/// Scalar-vector division.
template <ReadableVector Vector,
          std::convertible_to<vect_value_type_t<Vector>> T>
vect_copy_t<Vector> operator/(const Vector& v, const T& S) {
  vect_copy_t<Vector> result;
  result = v;
  result /= S;
  return result;
}

/// Element-wise product of two vectors.
template <ReadableVector Vector1, ReadableVector Vector2>
vect_copy_t<Vector1> elem_product(const Vector1& v1, const Vector2& v2) {
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

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/// Equality Comparison operator, component-wise.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator==(const Vector1& v1, const Vector2& v2) noexcept {
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

/// Inequality Comparison operator, component-wise.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator!=(const Vector1& v1, const Vector2& v2) noexcept {
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

/// Greater-than Comparison operator, Euclidean norm.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator>(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) > norm_2_sqr(v2));
}

/// Smaller-than Comparison operator, Euclidean norm.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator<(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) < norm_2_sqr(v2));
}

/// Greater-or-equal Comparison operator, Euclidean norm.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator>=(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) >= norm_2_sqr(v2));
}

/// Smaller-or-equal Comparison operator, Euclidean norm.
template <ReadableVector Vector1, ReadableVector Vector2>
bool operator<=(const Vector1& v1, const Vector2& v2) noexcept {
  return (norm_2_sqr(v1) <= norm_2_sqr(v2));
}

/// Prints a variable-size vector to an output stream as "(v1; v2; v3; ..; vN)".
template <typename T, unsigned int Size>
std::ostream& operator<<(std::ostream& out_stream, const vect<T, Size>& V) {
  out_stream << "(";
  if (V.size() > 0) {
    out_stream << V[0];
  }
  for (int i = 1; i < V.size(); ++i) {
    out_stream << "; " << V[i];
  }
  return out_stream << ")";
}

/// Reads a variable-size vector to an input stream as "(v1; v2; v3; ..; vN)".
template <typename T, unsigned int Size>
std::istream& operator>>(std::istream& in_stream, vect<T, Size>& V) {
  std::string tmp_str;
  std::getline(in_stream, tmp_str, '(');  // skip to opening bracket.
  std::getline(in_stream, tmp_str, ')');  // read to closing bracket.
  int sz = Size;
  if constexpr (Size == 0) {
    sz = std::count(tmp_str.begin(), tmp_str.end(), ';') + 1;
    V.resize(sz);
  }
  std::stringstream ss(tmp_str);
  std::string tmp2;
  for (int i = 0; (i < sz) && (ss >> V[i]); ++i) {
    std::getline(ss, tmp2, ';');
  }
  return in_stream;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_ALG_H_
