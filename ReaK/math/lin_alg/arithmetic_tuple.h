/**
 * \file arithmetic_tuple.h
 *
 * This library provides the arithmetic tuple class. This class is basically just a wrapper
 * of either the std::tuple class or the std::tuples::tuple class, either way, it provides
 * a meta-programming interface that is equivalent to the wrapped class and with the addition
 * of the support for all the basic arithmetic operators, requiring, of course, that these
 * arithmetic operators are also available on all the types contained in the tuple.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_MATH_LIN_ALG_ARITHMETIC_TUPLE_H_
#define REAK_MATH_LIN_ALG_ARITHMETIC_TUPLE_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/rtti/so_register_type.h"
#include "ReaK/core/rtti/typed_primitives.h"
#include "ReaK/core/serialization/archiver.h"

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace ReaK {

using std::get;

/// This meta-function computes a bool integral constant if the given type is an arithmetic-tuple.
/// \tparam Tuple The type to be tested as being an arithmetic-tuple or not.
template <typename Tuple>
struct is_arithmetic_tuple : std::false_type {};

template <typename Tuple>
static constexpr bool is_arithmetic_tuple_v = is_arithmetic_tuple<Tuple>::value;

/// This meta-function computes an integral constant describing the size (or length) of an arithmetic-tuple.
/// \tparam Tuple The arithmetic-tuple type.
template <typename Tuple>
struct arithmetic_tuple_size : std::integral_constant<std::size_t, 0> {};

template <typename Tuple>
struct arithmetic_tuple_size<const Tuple> : arithmetic_tuple_size<Tuple> {};

template <typename Tuple>
struct arithmetic_tuple_size<volatile Tuple> : arithmetic_tuple_size<Tuple> {};

template <typename Tuple>
struct arithmetic_tuple_size<const volatile Tuple>
    : arithmetic_tuple_size<Tuple> {};

template <typename Tuple>
static constexpr std::size_t arithmetic_tuple_size_v =
    arithmetic_tuple_size<Tuple>::value;

// Forward-declaration to specialize arithmetic_tuple_size and is_arithmetic_tuple.
template <typename... T>
class arithmetic_tuple;

/* Specialization, see general template docs. */
template <typename... T>
struct is_arithmetic_tuple<arithmetic_tuple<T...>> : std::true_type {};

/* Specialization, see general template docs. */
template <typename... T>
struct arithmetic_tuple_size<arithmetic_tuple<T...>>
    : std::integral_constant<std::size_t, sizeof...(T)> {};

// Implementation of tuple_for_each, needs to go before arithmetic_tuple,
// to implement its friend functions.
namespace {

template <typename Tuple, typename Func, std::size_t... Idx>
void tuple_for_each_impl(Tuple& arg1, Func f,
                         std::index_sequence<Idx...> /*unused*/) {
  // Using Louis Dionne's "swallow" trick:
  using swallow = int[];  // NOLINT
  (void)swallow{1, (f(get<Idx>(arg1)), void(), int())...};
}

template <typename Tuple1, typename Tuple2, typename Func, std::size_t... Idx>
void tuple_for_each_impl(Tuple1& arg1, Tuple2& arg2, Func f,
                         std::index_sequence<Idx...> /*unused*/) {
  // Using Louis Dionne's "swallow" trick:
  using swallow = int[];  // NOLINT
  (void)swallow{1, (f(get<Idx>(arg1), get<Idx>(arg2)), void(), int())...};
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Func,
          std::size_t... Idx>
void tuple_for_each_impl(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Func f,
                         std::index_sequence<Idx...> /*unused*/) {
  // Using Louis Dionne's "swallow" trick:
  using swallow = int[];  // NOLINT
  (void)swallow{
      1, (f(get<Idx>(arg1), get<Idx>(arg2), get<Idx>(arg3)), void(), int())...};
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Tuple4,
          typename Func, std::size_t... Idx>
void tuple_for_each_impl(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Tuple4& arg4,
                         Func f, std::index_sequence<Idx...> /*unused*/) {
  // Using Louis Dionne's "swallow" trick:
  using swallow = int[];  // NOLINT
  (void)swallow{
      1, (f(get<Idx>(arg1), get<Idx>(arg2), get<Idx>(arg3), get<Idx>(arg4)),
          void(), int())...};
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Tuple4,
          typename Tuple5, typename Func, std::size_t... Idx>
void tuple_for_each_impl(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Tuple4& arg4,
                         Tuple5& arg5, Func f,
                         std::index_sequence<Idx...> /*unused*/) {
  // Using Louis Dionne's "swallow" trick:
  using swallow = int[];  // NOLINT
  (void)swallow{1, (f(get<Idx>(arg1), get<Idx>(arg2), get<Idx>(arg3),
                      get<Idx>(arg4), get<Idx>(arg5)),
                    void(), int())...};
}

template <typename Tuple, typename Func>
void tuple_for_each(Tuple& arg1, Func f) {
  using idx_seq = std::make_index_sequence<arithmetic_tuple_size_v<Tuple>>;
  tuple_for_each_impl(arg1, f, idx_seq());
}

template <typename Tuple1, typename Tuple2, typename Func>
void tuple_for_each(Tuple1& arg1, Tuple2& arg2, Func f) {
  using idx_seq = std::make_index_sequence<arithmetic_tuple_size_v<Tuple1>>;
  tuple_for_each_impl(arg1, arg2, f, idx_seq());
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Func>
void tuple_for_each(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Func f) {
  using idx_seq = std::make_index_sequence<arithmetic_tuple_size_v<Tuple1>>;
  tuple_for_each_impl(arg1, arg2, arg3, f, idx_seq());
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Tuple4,
          typename Func>
void tuple_for_each(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Tuple4& arg4,
                    Func f) {
  using idx_seq = std::make_index_sequence<arithmetic_tuple_size_v<Tuple1>>;
  tuple_for_each_impl(arg1, arg2, arg3, arg4, f, idx_seq());
}

template <typename Tuple1, typename Tuple2, typename Tuple3, typename Tuple4,
          typename Tuple5, typename Func>
void tuple_for_each(Tuple1& arg1, Tuple2& arg2, Tuple3& arg3, Tuple4& arg4,
                    Tuple5& arg5, Func f) {
  using idx_seq = std::make_index_sequence<arithmetic_tuple_size_v<Tuple1>>;
  tuple_for_each_impl(arg1, arg2, arg3, arg4, arg5, f, idx_seq());
}

}  // namespace

/// This class template is a simple wrapper of a tuple with the addition of arithmetic operators.
/// This class is basically just a wrapper of the std::tuple class, and it provides
/// a meta-programming interface that is equivalent to std::tuple and with the addition
/// of the support for all the basic arithmetic operators, requiring, of course, that these
/// arithmetic operators are also available on all the types contained in the tuple.
/// \tparam T The types contained in the arithmetic-tuple.
template <typename... T>
class arithmetic_tuple : public std::tuple<T...> {
 public:
  using arithmetic_tuple_base_class = std::tuple<T...>;
  using self = arithmetic_tuple<T...>;

  template <typename... U>
  arithmetic_tuple(U&&... u)  // NOLINT
      : arithmetic_tuple_base_class(std::forward<U>(u)...) {}

  // Catch construction with other tuples, down-cast to std::tuple and let it
  // figure out if the param-pack matches <T...> or if it's just the first element.
  template <typename... U>
  arithmetic_tuple(arithmetic_tuple<U...>& rhs)  // NOLINT
      : arithmetic_tuple_base_class(static_cast<const std::tuple<U...>&>(rhs)) {
  }

  template <typename... U>
  arithmetic_tuple(const arithmetic_tuple<U...>& rhs)  // NOLINT
      : arithmetic_tuple_base_class(static_cast<const std::tuple<U...>&>(rhs)) {
  }

  template <typename... U>
  arithmetic_tuple(arithmetic_tuple<U...>&& rhs)  // NOLINT
      : arithmetic_tuple_base_class(static_cast<std::tuple<U...>&&>(rhs)) {}

  arithmetic_tuple(const self&) = default;
  arithmetic_tuple(self& rhs)
      : arithmetic_tuple(static_cast<const self&>(rhs)) {}
  arithmetic_tuple(self&&) noexcept = default;
  arithmetic_tuple& operator=(const self&) = default;
  arithmetic_tuple& operator=(self&&) noexcept = default;

  arithmetic_tuple_base_class& base() { return *this; }
  const arithmetic_tuple_base_class& base() const { return *this; }

  /// This function template is an overload of the addition operator on arithmetic-tuples.
  /// This function performs the addition of each element of the tuple, will only compile if
  /// all elements of the tuple support the addition operator.
  /// \param lhs Left-hand side of the addition.
  /// \param rhs Right-hand side of the addition.
  /// \return Added tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) + get<0>(rhs), get<1>(lhs) + get<1>(rhs), ...
  /// etc.).
  friend self operator+(const self& lhs, const self& rhs) {
    self result;
    tuple_for_each(
        result, lhs, rhs,
        [](auto& result_elem, const auto& lhs_elem, const auto& rhs_elem) {
          result_elem = lhs_elem + rhs_elem;
        });
    return result;
  }

  /// This function template is an overload of the subtraction operator on arithmetic-tuples.
  /// This function performs the subtraction of each element of the tuple, will only compile if
  /// all elements of the tuple support the subtraction operator.
  /// \param lhs Left-hand side of the subtraction.
  /// \param rhs Right-hand side of the subtraction.
  /// \return Subtracted tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) - get<0>(rhs), get<1>(lhs) - get<1>(rhs),
  /// ... etc.).
  friend self operator-(const self& lhs, const self& rhs) {
    self result;
    tuple_for_each(
        result, lhs, rhs,
        [](auto& result_elem, const auto& lhs_elem, const auto& rhs_elem) {
          result_elem = lhs_elem - rhs_elem;
        });
    return result;
  }

  /// This function template is an overload of the negation operator on arithmetic-tuples.
  /// This function performs the negation of each element of the tuple, will only compile if
  /// all elements of the tuple support the negation and assignment operators.
  /// \param lhs Left-hand side of the negation.
  /// \return Negated tuple, equivalent to make_arithmetic_tuple(-get<0>(lhs), -get<1>(lhs), ... etc.).
  friend self operator-(const self& lhs) {
    self result;
    tuple_for_each(result, lhs, [](auto& result_elem, const auto& lhs_elem) {
      result_elem = -lhs_elem;
    });
    return result;
  }

  /// This function template is an overload of the add-and-assign operator on arithmetic-tuples.
  /// This function performs the add-and-assign of each element of the tuple, will only compile if
  /// all elements of the tuple support the add-and-assign operator.
  /// \param lhs Left-hand side of the add-and-assign.
  /// \param rhs Right-hand side of the add-and-assign.
  /// \return Added-and-assigned tuple reference, equivalent to get<0>(lhs) += get<0>(rhs); get<1>(lhs) += get<1>(rhs); ...
  /// etc.
  friend self& operator+=(self& lhs, const self& rhs) {
    tuple_for_each(lhs, rhs, [](auto& lhs_elem, const auto& rhs_elem) {
      lhs_elem += rhs_elem;
    });
    return lhs;
  }

  /// This function template is an overload of the sub-and-assign operator on arithmetic-tuples.
  /// This function performs the sub-and-assign of each element of the tuple, will only compile if
  /// all elements of the tuple support the sub-and-assign operator.
  /// \param lhs Left-hand side of the sub-and-assign.
  /// \param rhs Right-hand side of the sub-and-assign.
  /// \return Subtracted-and-assigned tuple reference, equivalent to get<0>(lhs) -= get<0>(rhs); get<1>(lhs) -=
  /// get<1>(rhs); ... etc.
  friend self& operator-=(self& lhs, const self& rhs) {
    tuple_for_each(lhs, rhs, [](auto& lhs_elem, const auto& rhs_elem) {
      lhs_elem -= rhs_elem;
    });
    return lhs;
  }

  /// This function template is an overload of the multiplication operator on arithmetic-tuples.
  /// This function performs the multiplication of each element of the tuple, will only compile if
  /// all elements of the tuple support the multiplication operator.
  /// \param lhs Left-hand side of the multiplication.
  /// \param rhs Right-hand side of the multiplication.
  /// \return Multiplied tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) * get<0>(rhs), get<1>(lhs) * get<1>(rhs),
  /// ... etc.).
  friend self operator*(const self& lhs, const self& rhs) {
    self result;
    tuple_for_each(
        result, lhs, rhs,
        [](auto& result_elem, const auto& lhs_elem, const auto& rhs_elem) {
          result_elem = lhs_elem * rhs_elem;
        });
    return result;
  }

  /// This function template is an overload of the division operator on arithmetic-tuples.
  /// This function performs the division of each element of the tuple, will only compile if
  /// all elements of the tuple support the division operator.
  /// \param lhs Left-hand side of the division.
  /// \param rhs Right-hand side of the division.
  /// \return Divided tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) / get<0>(rhs), get<1>(lhs) / get<1>(rhs), ...
  /// etc.).
  friend self operator/(const self& lhs, const self& rhs) {
    self result;
    tuple_for_each(
        result, lhs, rhs,
        [](auto& result_elem, const auto& lhs_elem, const auto& rhs_elem) {
          result_elem = lhs_elem / rhs_elem;
        });
    return result;
  }

  /**
   * This function template is an overload of the scalar-multiplication operator on arithmetic-tuples.
   * This function performs the scalar-multiplication of each element of the tuple, will only compile if
   * all elements of the tuple support the scalar-multiplication operator.
   * \param lhs Left-hand side of the scalar-multiplication.
   * \param rhs Right-hand side of the scalar-multiplication (the scalar).
   * \return Multiplied tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) * rhs, get<1>(lhs) * rhs, ... etc.).
   * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
   * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same
   * type again).
   */
  template <typename Scalar>
  friend self operator*(const self& lhs, const Scalar& rhs) {
    self result;
    tuple_for_each(result, lhs, [&](auto& result_elem, const auto& lhs_elem) {
      result_elem = lhs_elem * rhs;
    });
    return result;
  }

  /**
   * This function template is an overload of the scalar-multiplication operator on arithmetic-tuples.
   * This function performs the scalar-multiplication of each element of the tuple, will only compile if
   * all elements of the tuple support the scalar-multiplication operator.
   * \param lhs Left-hand side of the scalar-multiplication (the scalar).
   * \param rhs Right-hand side of the scalar-multiplication.
   * \return Multiplied tuple, equivalent to make_arithmetic_tuple(lhs * get<0>(rhs), lhs * get<1>(rhs), ... etc.).
   * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
   * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same
   * type again).
   */
  template <typename Scalar>
  friend self operator*(const Scalar& lhs, const self& rhs) {
    self result;
    tuple_for_each(result, rhs, [&](auto& result_elem, const auto& rhs_elem) {
      result_elem = lhs * rhs_elem;
    });
    return result;
  }

  /**
   * This function template is an overload of the scalar-division operator on arithmetic-tuples.
   * This function performs the scalar-division of each element of the tuple, will only compile if
   * all elements of the tuple support the scalar-division operator.
   * \param lhs Left-hand side of the scalar-division.
   * \param rhs Right-hand side of the scalar-division (the scalar).
   * \return Divided tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) / rhs, get<1>(lhs) / rhs, ... etc.).
   * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
   * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same
   * type again).
   */
  template <typename Scalar>
  friend self operator/(const self& lhs, const Scalar& rhs) {
    self result;
    tuple_for_each(result, lhs, [&](auto& result_elem, const auto& lhs_elem) {
      result_elem = lhs_elem / rhs;
    });
    return result;
  }

  /// This function template is an overload of the multiply-and-assign operator on arithmetic-tuples.
  /// This function performs the multiply-and-assign of each element of the tuple, will only compile if
  /// all elements of the tuple support the multiply-and-assign operator.
  /// \param lhs Left-hand side of the multiply-and-assign.
  /// \param rhs Right-hand side of the multiply-and-assign.
  /// \return Multiplied tuple reference, equivalent to get<0>(lhs) *= get<0>(rhs); get<1>(lhs) *= get<1>(rhs); ... etc.
  /// \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
  friend self& operator*=(self& lhs, const self& rhs) {
    tuple_for_each(lhs, rhs, [](auto& lhs_elem, const auto& rhs_elem) {
      lhs_elem *= rhs_elem;
    });
    return lhs;
  }

  /// This function template is an overload of the divide-and-assign operator on arithmetic-tuples.
  /// This function performs the divide-and-assign of each element of the tuple, will only compile if
  /// all elements of the tuple support the divide-and-assign operator.
  /// \param lhs Left-hand side of the divide-and-assign.
  /// \param rhs Right-hand side of the divide-and-assign.
  /// \return Divided tuple reference, equivalent to get<0>(lhs) /= get<0>(rhs); get<1>(lhs) /= get<1>(rhs); ... etc.
  /// \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
  friend self& operator/=(self& lhs, const self& rhs) {
    tuple_for_each(lhs, rhs, [](auto& lhs_elem, const auto& rhs_elem) {
      lhs_elem /= rhs_elem;
    });
    return lhs;
  }

  /**
   * This function template is an overload of the scalar-multiply-and-assign operator on arithmetic-tuples.
   * This function performs the scalar-multiply-and-assign of each element of the tuple, will only compile if
   * all elements of the tuple support the scalar-multiply-and-assign operator.
   * \param lhs Left-hand side of the scalar-multiply-and-assign.
   * \param rhs Right-hand side of the scalar-multiply-and-assign (the scalar).
   * \return Multiplied tuple reference, equivalent to get<0>(lhs) *= rhs; get<1>(lhs) *= rhs; ... etc.
   * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
   * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same
   * type again).
   */
  template <typename Scalar>
  friend self& operator*=(self& lhs, const Scalar& rhs) {
    tuple_for_each(lhs, [&](auto& lhs_elem) { lhs_elem *= rhs; });
    return lhs;
  }

  /**
   * This function template is an overload of the scalar-divide-and-assign operator on arithmetic-tuples.
   * This function performs the scalar-divide-and-assign of each element of the tuple, will only compile if
   * all elements of the tuple support the scalar-divide-and-assign operator.
   * \param lhs Left-hand side of the scalar-divide-and-assign.
   * \param rhs Right-hand side of the scalar-divide-and-assign (the scalar).
   * \return Divided tuple reference, equivalent to get<0>(lhs) /= rhs; get<1>(lhs) /= rhs; ... etc.
   * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
   * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same
   * type again).
   */
  template <typename Scalar>
  friend self& operator/=(self& lhs, const Scalar& rhs) {
    tuple_for_each(lhs, [&](auto& lhs_elem) { lhs_elem /= rhs; });
    return lhs;
  }

  /// This function template is an overload of the unnamed archive-output operator for arithmetic-tuples.
  /// This function performs the unnamed archive-output of each element of the tuple, will only compile if
  /// all elements of the tuple support the unnamed archive-output operator.
  /// \param out The output archive.
  /// \param rhs The arithmetic-tuple object to output on the archive.
  /// \return The output archive.
  /// \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
  friend std::ostream& operator<<(std::ostream& out, const self& rhs) {
    out << "( ";
    bool is_first = true;
    tuple_for_each(rhs, [&](const auto& rhs_elem) {
      if (!is_first) {
        out << "; ";
      }
      out << rhs_elem;
      is_first = false;
    });
    out << ")";
    return out;
  }
};

template <typename... U>
arithmetic_tuple(U&&... u) -> arithmetic_tuple<std::decay_t<U>...>;

/// This function template can be used to create an arithmetic-tuple.
/// \tparam T The types contained in the arithmetic-tuple.
/// \param t The values that make up the arithmetic-tuple.
/// \return An arithmetic-tuple.
template <typename... T>
inline arithmetic_tuple<std::decay_t<T>...> make_arithmetic_tuple(T&&... t) {
  return arithmetic_tuple<std::decay_t<T>...>(std::forward<T>(t)...);
}

namespace detail {
template <int Idx, typename U, typename MaybeConstBase>
static U& get_by_type_impl(MaybeConstBase& value) {
  static_assert(Idx > 0, "The given type T was not found in the tuple!");
  if constexpr (std::is_same_v<
                    std::remove_const_t<U>,
                    std::tuple_element_t<
                        Idx - 1, std::remove_const_t<MaybeConstBase>>>) {
    using ReaK::get;
    return get<Idx - 1>(value);
  } else {
    return get_by_type_impl<Idx - 1, U, MaybeConstBase>(value);
  }
}
}  // namespace detail

template <typename U, typename... T>
U& get_by_type(arithmetic_tuple<T...>& value) {
  return detail::get_by_type_impl<sizeof...(T), U>(value.base());
}

template <typename U, typename... T>
const U& get_by_type(const arithmetic_tuple<T...>& value) {
  return detail::get_by_type_impl<sizeof...(T), const U>(value.base());
}

namespace detail {

template <typename Vector, typename RhsType>
void to_vect_impl(Vector& lhs, const RhsType& rhs) {
  if constexpr (is_readable_vector_v<RhsType>) {
    const auto i = lhs.size();
    lhs.resize(i + rhs.size());
    for (std::size_t j = 0; j < rhs.size(); ++j) {
      lhs[i + j] = rhs[j];
    }
  } else if constexpr (is_arithmetic_tuple_v<RhsType>) {
    tuple_for_each(rhs,
                   [&](const auto& rhs_elem) { to_vect_impl(lhs, rhs_elem); });
  } else {
    lhs.resize(lhs.size() + 1);
    lhs[lhs.size() - 1] = rhs;
  }
}

template <typename LhsType, typename Vector>
void from_vect_impl(LhsType& lhs, const Vector& rhs, std::size_t& i) {
  if constexpr (is_writable_vector_v<LhsType>) {
    for (std::size_t j = 0; j < rhs.size(); ++j, ++i) {
      lhs[j] = rhs[i];
    }
  } else if constexpr (is_arithmetic_tuple_v<LhsType>) {
    tuple_for_each(lhs,
                   [&](auto& lhs_elem) { from_vect_impl(lhs_elem, rhs, i); });
  } else {
    lhs = rhs[i++];
  }
}

}  // namespace detail

/// This function template converts anything to a vector, as long as the fundamental
/// value-type is compatible with the given output type.
/// \param v Something that can be converted to a vector.
/// \return A vector with the flattened content of the input.
template <typename ValueType, typename VectorType>
auto to_vect(const VectorType& v) {
  if constexpr (is_readable_vector_v<VectorType>) {
    return v;
  } else {
    using ReaK::detail::to_vect_impl;  // for ADL
    vect_n<ValueType> result_v;
    to_vect_impl(result_v, v);
    return result_v;
  }
}

/// This function template converts a vector into anything, as long as the fundamental
/// value-type is compatible with the vector type.
/// \param v A vector.
/// \return A vector with the flattened content of the input.
template <typename OutputType, typename VectorType>
OutputType from_vect(const VectorType& v) {
  if constexpr (is_writable_vector_v<OutputType>) {
    return OutputType(v);
  } else {
    using ReaK::detail::from_vect_impl;  // for ADL
    OutputType result_v;
    std::size_t i = 0;
    from_vect_impl(result_v, v, i);
    return result_v;
  }
}

namespace serialization {

/// This function template is an overload of the unnamed archive-output operator for arithmetic-tuples.
/// This function performs the unnamed archive-output of each element of the tuple, will only compile if
/// all elements of the tuple support the unnamed archive-output operator.
/// \param out The output archive.
/// \param rhs The arithmetic-tuple object to output on the archive.
/// \return The output archive.
template <typename... Args>
oarchive& operator<<(oarchive& out, const arithmetic_tuple<Args...>& rhs) {
  tuple_for_each(rhs, [&](const auto& rhs_elem) { out << rhs_elem; });
  return out;
}

/// This function template is an overload of the named archive-output operator for arithmetic-tuples.
/// This function performs the named archive-output of each element of the tuple, will only compile if
/// all elements of the tuple support the named archive-output operator.
/// \param out The output archive.
/// \param rhs The arithmetic-tuple object to output on the archive.
/// \return The output archive.
template <typename... Args>
oarchive& operator&(
    oarchive& out,
    const std::pair<std::string, const arithmetic_tuple<Args...>&>& rhs) {
  std::size_t cur_idx = 0;
  tuple_for_each(rhs.second, [&](const auto& rhs_elem) {
    std::stringstream ss(rhs.first);
    ss << "_q" << cur_idx;
    out& serialization::make_save_nvp(ss.str(), rhs_elem);
    ++cur_idx;
  });
  return out;
}

/// This function template is an overload of the unnamed archive-input operator for arithmetic-tuples.
/// This function performs the unnamed archive-input of each element of the tuple, will only compile if
/// all elements of the tuple support the unnamed archive-input operator.
/// \param in The input archive.
/// \param rhs The arithmetic-tuple object to input from the archive.
/// \return The input archive.
template <typename... Args>
iarchive& operator>>(iarchive& in, arithmetic_tuple<Args...>& rhs) {
  tuple_for_each(rhs, [&](auto& rhs_elem) { in >> rhs_elem; });
  return in;
}

/// This function template is an overload of the named archive-input operator for arithmetic-tuples.
/// This function performs the named archive-input of each element of the tuple, will only compile if
/// all elements of the tuple support the named archive-input operator.
/// \param in The input archive.
/// \param rhs The arithmetic-tuple object to input from the archive.
/// \return The input archive.
template <typename... Args>
iarchive& operator&(
    iarchive& in,
    const std::pair<std::string, arithmetic_tuple<Args...>&>& rhs) {
  std::size_t cur_idx = 0;
  tuple_for_each(rhs.second, [&](auto& rhs_elem) {
    std::stringstream ss(rhs.first);
    ss << "_q" << cur_idx;
    in& serialization::make_load_nvp(ss.str(), rhs_elem);
    ++cur_idx;
  });
  return in;
}

}  // namespace serialization

namespace rtti {

template <typename... T>
struct get_type_id<arithmetic_tuple<T...>> {
  static constexpr unsigned int ID = 0x0000002C;
  static constexpr auto type_name = std::string_view{"arithmetic_tuple"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const arithmetic_tuple<T...>&;
  using load_type = arithmetic_tuple<T...>&;
};

template <typename Tail, typename... T>
struct get_type_info<arithmetic_tuple<T...>, Tail> {
  using type = type_id<
      arithmetic_tuple<T...>,
      typename get_type_info_seq<T...>::template with_tail<Tail>::type::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<arithmetic_tuple<T...>>::type_name,
                  lsl_left_bracket, get_type_info_seq<T...>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};

}  // namespace rtti
}  // namespace ReaK

namespace ReaK {

template <int Idx, typename Tuple>
struct arithmetic_tuple_element {
  using type = std::tuple_element_t<Idx, Tuple>;
};

template <int Idx, typename... T>
struct arithmetic_tuple_element<Idx, arithmetic_tuple<T...>> {
  using type = std::tuple_element_t<
      Idx, typename arithmetic_tuple<T...>::arithmetic_tuple_base_class>;
};

template <int Idx, typename... T>
struct arithmetic_tuple_element<Idx, const arithmetic_tuple<T...>> {
  using type = std::tuple_element_t<
      Idx, const typename arithmetic_tuple<T...>::arithmetic_tuple_base_class>;
};

template <int Idx, typename... T>
struct arithmetic_tuple_element<Idx, volatile arithmetic_tuple<T...>> {
  using type = std::tuple_element_t<
      Idx,
      volatile typename arithmetic_tuple<T...>::arithmetic_tuple_base_class>;
};

template <int Idx, typename... T>
struct arithmetic_tuple_element<Idx, const volatile arithmetic_tuple<T...>> {
  using type =
      std::tuple_element_t<Idx, const volatile typename arithmetic_tuple<
                                    T...>::arithmetic_tuple_base_class>;
};

template <int Idx, typename Tuple>
using arithmetic_tuple_element_t =
    typename arithmetic_tuple_element<Idx, Tuple>::type;

}  // namespace ReaK

namespace ReaK {

namespace detail {

template <int Idx, typename T, typename Tuple>
struct arithmetic_tuple_index_of_impl;  // forward-decl

template <typename T, typename Tuple>
struct arithmetic_tuple_index_of_impl<0, T, Tuple> {
  static constexpr int value = std::numeric_limits<int>::max();
};

template <int Idx, typename T, typename Tuple>
struct arithmetic_tuple_index_of_impl {
  static constexpr int value =
      (std::is_same_v<T, arithmetic_tuple_element_t<Idx - 1, Tuple>>
           ? Idx - 1
           : arithmetic_tuple_index_of_impl<Idx - 1, T, Tuple>::value);
};

}  // namespace detail

template <typename T, typename Tuple>
static constexpr int arithmetic_tuple_index_of_v =
    detail::arithmetic_tuple_index_of_impl<arithmetic_tuple_size_v<Tuple>, T,
                                           Tuple>::value;

template <typename T, typename Tuple>
static constexpr bool arithmetic_tuple_has_type_v =
    (arithmetic_tuple_index_of_v<T, Tuple> != std::numeric_limits<int>::max());

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_ARITHMETIC_TUPLE_H_
