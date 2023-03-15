/**
 * \file vect_views.h
 *
 * This library provides a number of class templates to create vector views. A vector
 * view simply means that a sub-part of a vector is used as if it was a vector in its
 * own right. This can be very useful to set sub-parts to other values or to use a
 * sub-part in a vector expression (e.g. applying the operation on the entire vector
 * is not practical or efficient).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_MATH_LIN_ALG_VECT_VIEWS_H_
#define REAK_MATH_LIN_ALG_VECT_VIEWS_H_

#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_traits.h"

#include <stdexcept>
#include <type_traits>

namespace ReaK {

/**
 * This function can be used to generate a range of indices for the matrices.
 * A range is simple a pair of first and last indices (notice that the last index is
 * included in the range).
 */
inline std::pair<int, int> range(int aFirst, int aLast) {
  return {aFirst, aLast};
}

/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class takes a const reference to the given vector.
 * \tparam Vector A readable vector type.
 */
template <typename Vector>
class vect_const_ref_view {
 public:
  using self = vect_const_ref_view<Vector>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using iterator = typename vect_traits<Vector>::iterator;
  using const_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr std::size_t dimensions = vect_traits<Vector>::dimensions;

  BOOST_CONCEPT_ASSERT((ReadableVectorConcept<Vector>));

 private:
  const Vector* v;
  int offset;
  int count;

  self& operator=(const self&);
  explicit vect_const_ref_view(Vector&&);
  vect_const_ref_view(Vector&&, int, int aOffset = 0);

 public:
  /**
   * Constructs the sub-vector which represents the entire vector.
   * \param aV The vector from which the sub-part is taken.
   */
  explicit vect_const_ref_view(const Vector& aV) noexcept
      : v(&aV), offset(0), count(aV.size()) {}

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-part is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_const_ref_view(const Vector& aV, int aCount, int aOffset = 0) noexcept
      : v(&aV), offset(aOffset), count(aCount) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[](int i) const noexcept { return (*v)[offset + i]; }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<int, int>& r) const noexcept {
    return vect_const_ref_view<self>(*this, r.second - r.first, r.first);
  }

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()(int i) const noexcept { return (*v)[offset + i]; }

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  int size() const noexcept { return count; }
  /**
   * Returns the max-size of the vector.
   */
  int max_size() const noexcept { return count; }
  /**
   * Returns the capacity of the vector.
   */
  int capacity() const noexcept { return count; }
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return (count == 0); }

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return v->begin() + offset; }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return v->begin() + offset + count; }
};

template <typename Vector>
struct is_readable_vector<vect_const_ref_view<Vector>> {
  static constexpr bool value = is_readable_vector_v<Vector>;
  using type = is_readable_vector<Vector>;
};

template <typename Vector>
struct is_writable_vector<vect_const_ref_view<Vector>> {
  static constexpr bool value = false;
  using type = is_writable_vector<vect_const_ref_view<Vector>>;
};

template <typename Vector>
struct is_resizable_vector<vect_const_ref_view<Vector>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<vect_const_ref_view<Vector>>;
};

/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class takes a reference to the given vector.
 * \tparam Vector A readable vector type.
 */
template <typename Vector>
class vect_ref_view {
 public:
  using self = vect_ref_view<Vector>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using iterator = typename vect_traits<Vector>::iterator;
  using const_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr std::size_t dimensions = vect_traits<Vector>::dimensions;

  BOOST_CONCEPT_ASSERT((ReadableVectorConcept<Vector>));

 private:
  Vector* v;
  int offset;
  int count;

 public:
  /**
   * Constructs the sub-vector which represents the entire vector.
   * \param aV The vector from which the sub-part is taken.
   */
  explicit vect_ref_view(Vector& aV) noexcept
      : v(&aV), offset(0), count(aV.size()) {}

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-part is taken.
   * \param aCount The number of elements for the sub-part.
   * \param aOffset The offset from the start of the vector.
   */
  vect_ref_view(Vector& aV, int aCount, int aOffset = 0) noexcept
      : v(&aV), offset(aOffset), count(aCount) {}

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*v)[offset + i] = rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator+=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*v)[offset + i] += rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator-=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*v)[offset + i] -= rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Scalar>
  self& operator*=(const Scalar& rhs) noexcept {
    for (int i = 0; i < count; ++i) {
      (*v)[offset + i] *= rhs;
    }
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-write access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  reference operator[](int i) noexcept { return (*v)[offset + i]; }
  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[](int i) const noexcept { return (*v)[offset + i]; }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_ref_view<self> operator[](const std::pair<int, int>& r) noexcept {
    return vect_ref_view<self>(*this, r.second - r.first, r.first);
  }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<int, int>& r) const noexcept {
    return vect_const_ref_view<self>(*this, r.second - r.first, r.first);
  }

  /**
   * Vector indexing operator, accessor for read/write.
   * TEST PASSED
   */
  reference operator()(int i) noexcept { return (*v)[offset + i]; }

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()(int i) const noexcept { return (*v)[offset + i]; }

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  int size() const noexcept { return count; }
  /**
   * Returns the max-size of the vector.
   */
  int max_size() const noexcept { return count; }
  /**
   * Returns the capacity of the vector.
   */
  int capacity() const noexcept { return count; }
  /**
   * Resizes the vector.
   */
  void resize(int sz, value_type c = value_type()) const noexcept {}
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return (count == 0); }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(int sz) const noexcept {}

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() noexcept { return v->begin() + offset; }
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return v->begin() + offset; }
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() noexcept { return v->begin() + offset + count; }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return v->begin() + offset + count; }
};

template <typename Vector>
struct is_readable_vector<vect_ref_view<Vector>> {
  static constexpr bool value = is_readable_vector_v<Vector>;
  using type = is_readable_vector<Vector>;
};

template <typename Vector>
struct is_writable_vector<vect_ref_view<Vector>> {
  static constexpr bool value = is_writable_vector_v<Vector>;
  using type = is_writable_vector<Vector>;
};

template <typename Vector>
struct is_resizable_vector<vect_ref_view<Vector>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<vect_ref_view<Vector>>;
};

/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class makes a copy of the given vector (it is mainly meant to harmonize syntax
 * when rvalue vectors are involved, requires C++11).
 * \tparam Vector A readable vector type.
 */
template <typename Vector>
class vect_copy_view {
 public:
  using self = vect_copy_view<Vector>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using iterator = typename vect_traits<Vector>::iterator;
  using const_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr std::size_t dimensions = vect_traits<Vector>::dimensions;

  BOOST_CONCEPT_ASSERT((ReadableVectorConcept<Vector>));

 private:
  Vector v;
  int offset;
  int count;

 public:
  /**
   * Default constructor.
   */
  vect_copy_view() : v(), offset(0), count(0) {}

  /**
   * Constructs the sub-vector which represents the entire vector.
   */
  explicit vect_copy_view(const Vector& aV)
      : v(aV), offset(0), count(aV.size()) {}

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-block is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_copy_view(const Vector& aV, int aCount, int aOffset = 0)
      : v(aV), offset(aOffset), count(aCount) {}

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit vect_copy_view(Vector&& aV) : v(std::move(aV)), offset(0), count(0) {
    count = v.size();
  }

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-block is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_copy_view(Vector&& aV, int aCount, int aOffset = 0)
      : v(std::move(aV)), offset(aOffset), count(aCount) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.v, rhs.v);
    swap(lhs.offset, rhs.offset);
    swap(lhs.count, rhs.count);
  }

  vect_copy_view(const self& rhs) = default;
  vect_copy_view(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) noexcept = default;

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      v[offset + i] = rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator+=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      v[offset + i] += rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Vector2>
  self& operator-=(const Vector2& rhs) {
    static_assert(is_readable_vector_v<Vector2>);
    if (rhs.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      v[offset + i] -= rhs[i];
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Scalar>
  self& operator*=(const Scalar& rhs) noexcept {
    for (int i = 0; i < count; ++i) {
      v[offset + i] *= rhs;
    }
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-write access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  reference operator[](int i) noexcept { return v[offset + i]; }
  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[](int i) const noexcept { return v[offset + i]; }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_ref_view<self> operator[](const std::pair<int, int>& r) noexcept {
    return vect_ref_view<self>(*this, r.second - r.first, r.first);
  }

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view<self> operator[](
      const std::pair<int, int>& r) const noexcept {
    return vect_const_ref_view<self>(*this, r.second - r.first, r.first);
  }

  /**
   * Vector indexing operator, accessor for read/write.
   * TEST PASSED
   */
  reference operator()(int i) noexcept { return v[offset + i]; }

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()(int i) const noexcept { return v[offset + i]; }

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  int size() const noexcept { return count; }
  /**
   * Returns the max-size of the vector.
   */
  int max_size() const noexcept { return count; }
  /**
   * Returns the capacity of the vector.
   */
  int capacity() const noexcept { return count; }
  /**
   * Resizes the vector.
   */
  void resize(int sz, value_type c = value_type()) const noexcept {}
  /**
   * Checks if the vector is empty.
   */
  bool empty() const noexcept { return (count == 0); }
  /**
   * Reserve a capacity for the vector.
   */
  void reserve(int sz) const noexcept {}

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() noexcept { return v.begin() + offset; }
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const noexcept { return v.begin() + offset; }
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() noexcept { return v.begin() + offset + count; }
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const noexcept { return v.begin() + offset + count; }
};

template <typename Vector>
struct is_readable_vector<vect_copy_view<Vector>> {
  static constexpr bool value = is_readable_vector_v<Vector>;
  using type = is_readable_vector<Vector>;
};

template <typename Vector>
struct is_writable_vector<vect_copy_view<Vector>> {
  static constexpr bool value = is_writable_vector_v<Vector>;
  using type = is_writable_vector<Vector>;
};

template <typename Vector>
struct is_resizable_vector<vect_copy_view<Vector>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<vect_copy_view<Vector>>;
};

template <typename Vector>
struct vect_copy_view_factory {
  Vector v;
  explicit vect_copy_view_factory(const Vector& aV) : v(aV) {}
  explicit vect_copy_view_factory(Vector&& aV) : v(std::move(aV)) {}
  vect_copy_view<Vector> operator[](const std::pair<int, int>& indices) {
    return vect_copy_view<Vector>(std::move(v), indices.second - indices.first,
                                  indices.first);
  }
};

template <typename Vector>
struct vect_ref_view_factory {
  Vector& v;
  explicit vect_ref_view_factory(Vector& aV) noexcept : v(aV) {}
  vect_ref_view<Vector> operator[](
      const std::pair<int, int>& indices) noexcept {
    return vect_ref_view<Vector>(v, indices.second - indices.first,
                                 indices.first);
  }
};

template <typename Vector>
struct vect_const_ref_view_factory {
  const Vector& v;
  explicit vect_const_ref_view_factory(const Vector& aV) noexcept : v(aV) {}
  vect_const_ref_view<Vector> operator[](
      const std::pair<int, int>& indices) noexcept {
    return vect_const_ref_view<Vector>(v, indices.second - indices.first,
                                       indices.first);
  }
};

template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_ref_view_factory<Vector>>
sub(Vector& V) noexcept {
  return vect_ref_view_factory<Vector>(V);
}

template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>,
                 vect_const_ref_view_factory<Vector>>
sub(const Vector& V) noexcept {
  return vect_const_ref_view_factory<Vector>(V);
}

template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_copy_view_factory<Vector>>
sub_copy(const Vector& V) {
  return vect_copy_view_factory<Vector>(V);
}

template <typename Vector>
std::enable_if_t<is_readable_vector_v<Vector>, vect_copy_view_factory<Vector>>
sub(Vector&& V) {
  return vect_copy_view_factory<Vector>(std::move(V));
}
}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_VIEWS_H_
