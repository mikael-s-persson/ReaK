/**
 * \file vect_index_iterator.hpp
 *
 * This library provides classes to create an iterator into a vector where that
 * vector only provides indexing capabilities. This can be useful as an adaptor of
 * a vector class that cannot sensibly provide an iterator interface natively.
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

#ifndef REAK_VECT_INDEX_ITERATOR_HPP
#define REAK_VECT_INDEX_ITERATOR_HPP

#include "ReaK/math/lin_alg/vect_traits.hpp"

#include <iterator>
#include <stdexcept>

namespace ReaK {

/**
 * This class template implements an iterator via indexing into a vector, which it
 * takes by reference (internally held by pointer, to be copyable).
 * \tparam Vector A readable vector type.
 */
template <typename Vector>
class vect_index_iter {
 public:
  using value_type = vect_value_type_t<Vector>;
  using size_type = typename vect_traits<Vector>::size_type;
  using reference = typename vect_traits<Vector>::reference;
  using difference_type = typename vect_traits<Vector>::difference_type;
  using pointer = typename vect_traits<Vector>::pointer;
  using iterator_category = std::random_access_iterator_tag;
  using self = vect_index_iter<Vector>;

 private:
  Vector* v;    ///< Holds a reference to the vector.
  size_type i;  ///< Holds the current index into the vector.

 public:
  /**
   * Constructs an index-iterator from a vector reference, automatically starts at one-past-last element (it creates the
   * end-iterator).
   */
  explicit vect_index_iter(Vector& aV) : v(&aV), i(aV.size()) {}
  /**
   * Constructs an index-iterator from a vector reference and a starting index.
   */
  vect_index_iter(Vector& aV, size_type aI) : v(&aV), i(aI) {}

  // compiler-generated copy-constructor and assignment operator are correct.

  /**
   * Pre-increment operator.
   */
  self& operator++() {
    ++i;
    return *this;
  }
  /**
   * Post-increment operator.
   */
  self operator++(int) {
    self tmp = *this;
    ++i;
    return tmp;
  }
  /**
   * Add-and-store operator.
   */
  self& operator+=(difference_type aStep) {
    i += aStep;
    return *this;
  }
  /**
   * Pre-decrement operator.
   */
  self& operator--() {
    --i;
    return *this;
  }
  /**
   * Post-decrement operator.
   */
  self operator--(int) {
    self tmp = *this;
    --i;
    return tmp;
  }
  /**
   * Sub-and-store operator.
   */
  self& operator-=(difference_type aStep) {
    i -= aStep;
    return *this;
  }
  /**
   * Indexing operator.
   */
  reference operator[](difference_type aIdx) const { return (*v)[aIdx]; }
  /**
   * Dereference operator.
   */
  reference operator*() const { return (*v)[i]; }
  /**
   * Member-access operator.
   */
  pointer operator->() const { return &(*v)[i]; }

  /**
   * Addition operator.
   */
  friend self operator+(self it, difference_type n) { return it += n; }
  /**
   * Addition operator.
   */
  friend self operator+(difference_type n, self it) { return it += n; }

  /**
   * Subtraction operator.
   */
  friend difference_type operator-(const self& it1, const self& it2) {
    return it1.i - it2.i;
  }

  /**
   * Equality operator.
   */
  friend bool operator==(const self& it1, const self& it2) {
    return it1.i == it2.i;
  }
  /**
   * Inequality operator.
   */
  friend bool operator!=(const self& it1, const self& it2) {
    return it1.i != it2.i;
  }
  /**
   * Less-than operator.
   */
  friend bool operator<(const self& it1, const self& it2) {
    return it1.i < it2.i;
  }
  /**
   * Greater-than operator.
   */
  friend bool operator>(const self& it1, const self& it2) {
    return it1.i > it2.i;
  }
  /**
   * Less-or-equal-than operator.
   */
  friend bool operator<=(const self& it1, const self& it2) {
    return it1.i <= it2.i;
  }
  /**
   * Greater-or-equal-than operator.
   */
  friend bool operator>=(const self& it1, const self& it2) {
    return it1.i >= it2.i;
  }
};

/**
 * This class template implements an iterator via indexing into a vector, which it
 * takes by const-reference (internally held by const-pointer, to be copyable).
 * \tparam Vector A readable vector type.
 */
template <typename Vector>
class vect_index_const_iter {
 public:
  using value_type = vect_value_type_t<Vector>;
  using size_type = typename vect_traits<Vector>::size_type;
  using reference = typename vect_traits<Vector>::const_reference;
  using difference_type = typename vect_traits<Vector>::difference_type;
  using pointer = typename vect_traits<Vector>::const_pointer;
  using iterator_category = std::random_access_iterator_tag;
  using self = vect_index_const_iter<Vector>;

 private:
  const Vector* v;
  size_type i;

 public:
  /**
   * Constructs an index-iterator from a vector reference, automatically starts at one-past-last element (it creates the
   * end-iterator).
   */
  explicit vect_index_const_iter(const Vector& aV) : v(&aV), i(aV.size()){};
  /**
   * Constructs an index-iterator from a vector reference and a starting index.
   */
  vect_index_const_iter(const Vector& aV, size_type aI) : v(&aV), i(aI) {}

  /**
   * Pre-increment operator.
   */
  self& operator++() {
    ++i;
    return *this;
  }
  /**
   * Post-increment operator.
   */
  self operator++(int) {
    self tmp = *this;
    ++i;
    return tmp;
  }
  /**
   * Add-and-store operator.
   */
  self& operator+=(difference_type aStep) {
    i += aStep;
    return *this;
  }
  /**
   * Pre-decrement operator.
   */
  self& operator--() {
    --i;
    return *this;
  }
  /**
   * Post-decrement operator.
   */
  self operator--(int) {
    self tmp = *this;
    --i;
    return tmp;
  }
  /**
   * Sub-and-store operator.
   */
  self& operator-=(difference_type aStep) {
    i -= aStep;
    return *this;
  }
  /**
   * Indexing operator.
   */
  reference operator[](difference_type aIdx) const { return (*v)[aIdx]; }
  /**
   * Dereference operator.
   */
  reference operator*() const { return (*v)[i]; }
  /**
   * Member-access operator.
   */
  pointer operator->() const { return &(*v)[i]; }

  /**
   * Addition operator.
   */
  friend self operator+(self it, difference_type n) { return it += n; }
  /**
   * Addition operator.
   */
  friend self operator+(difference_type n, self it) { return it += n; }

  /**
   * Subtraction operator.
   */
  friend difference_type operator-(const self& it1, const self& it2) {
    return it1.i - it2.i;
  }

  /**
   * Equality operator.
   */
  friend bool operator==(const self& it1, const self& it2) {
    return it1.i == it2.i;
  }
  /**
   * Inequality operator.
   */
  friend bool operator!=(const self& it1, const self& it2) {
    return it1.i != it2.i;
  }
  /**
   * Less-than operator.
   */
  friend bool operator<(const self& it1, const self& it2) {
    return it1.i < it2.i;
  }
  /**
   * Greater-than operator.
   */
  friend bool operator>(const self& it1, const self& it2) {
    return it1.i > it2.i;
  }
  /**
   * Less-or-equal-than operator.
   */
  friend bool operator<=(const self& it1, const self& it2) {
    return it1.i <= it2.i;
  }
  /**
   * Greater-or-equal-than operator.
   */
  friend bool operator>=(const self& it1, const self& it2) {
    return it1.i >= it2.i;
  }
};
}  // namespace ReaK

#endif
