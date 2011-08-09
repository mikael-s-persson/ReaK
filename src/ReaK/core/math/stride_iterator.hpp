/**
 * \file stride_iterator.hpp
 * 
 * This library provides a class template which models the RandomAccessIterator concept of the 
 * STL and has the behaviour of adding an automatic stride to an underlying random-access iterator.
 * This class is used in ReaK to implement matrix iterators.
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

#ifndef STRIDE_ITERATOR_HPP
#define STRIDE_ITERATOR_HPP

#include <iterator>
#include <stdexcept>

namespace ReaK {

/**
 * This class template models a RandomAccessIterator concept (STL) and adds an
 * automatic stride to an underlying random-access iterator.
 * \tparam RAIter A random-access-iterator type.
 * \tparam Stride The stride (if known at compile-time), if not, the default is 0, signifying a dynamic stride value.
 */
template<class RAIter, unsigned int Stride = 0>
class stride_iterator {
  private:
    RAIter pos; ///< Holds the underlying iterator.
  public:
    typedef typename std::iterator_traits<RAIter>::value_type value_type;
    typedef typename std::iterator_traits<RAIter>::reference reference;
    typedef typename std::iterator_traits<RAIter>::difference_type difference_type;
    typedef typename std::iterator_traits<RAIter>::pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef stride_iterator<RAIter,Stride> self;

    /**
     * Default constructor.
     */
    stride_iterator( ) : pos(NULL) { };
    /**
     * Copy-constructor.
     */
    stride_iterator(const self& rhs) : pos(rhs.pos) { };
    /**
     * Constructs the stride iterator at a given position.
     * \param aPos the position of the iterator.
     */
    explicit stride_iterator(RAIter aPos) : pos(aPos) { };
    
    /**
     * Returns the underlying iterator.
     * \return the underlying iterator.
     */
    RAIter base() const { return pos; };

    /**
     * Standard assignment operator.
     */
    self& operator=(const self& rhs) { pos = rhs.pos; };
    /**
     * Pre-increment operator.
     */
    self& operator++() { pos += Stride; return *this; };
    /**
     * Post-increment operator.
     */
    self operator++(int) { self tmp = *this; pos += Stride; return tmp; };
    /**
     * Add-and-store operator.
     */
    self& operator+=(difference_type aStep) { pos += aStep * Stride; return *this; };
    /**
     * Pre-decrement operator.
     */
    self& operator--() { pos -= Stride; return *this; };
    /**
     * Post-decrement operator.
     */
    self operator--(int) { self tmp = *this; pos -= Stride; return tmp; };
    /**
     * Sub-and-store operator.
     */
    self& operator-=(difference_type aStep) { pos -= aStep * Stride; return *this; };
    /**
     * Indexing operator.
     */
    reference operator[](difference_type aIdx) { return pos[aIdx * Stride]; };
    /**
     * Dereference operator.
     */
    reference operator*() { return *pos; };
    /**
     * Member-access operator.
     */
    pointer operator->() { return pos.operator->(); };
    
    /**
     * Addition operator.
     */
    friend self operator+(self it, difference_type n) { return it += n; };
    /**
     * Addition operator.
     */
    friend self operator+(difference_type n, self it) { return it += n; };
    
    /**
     * Subtraction operator.
     */
    friend difference_type operator-(const self& it1, const self& it2) { return (it1.pos - it2.pos) / Stride; };

    /**
     * Equality operator.
     */
    friend bool operator==(const self& it1, const self& it2) { return it1.pos == it2.pos; };
    /**
     * Inequality operator.
     */
    friend bool operator!=(const self& it1, const self& it2) { return it1.pos != it2.pos; };
    /**
     * Less-than operator.
     */
    friend bool operator<(const self& it1, const self& it2) { return it1.pos < it2.pos; };
    /**
     * Greater-than operator.
     */
    friend bool operator>(const self& it1, const self& it2) { return it1.pos > it2.pos; };
    /**
     * Less-or-equal-than operator.
     */
    friend bool operator<=(const self& it1, const self& it2) { return it1.pos <= it2.pos; };
    /**
     * Greater-or-equal-than operator.
     */
    friend bool operator>=(const self& it1, const self& it2) { return it1.pos >= it2.pos; };
};

//Default template specialization, i.e., stride not given a compile time.
/**
 * This class template specialization models a RandomAccessIterator concept (STL) and adds an
 * automatic and dynamic stride to an underlying random-access iterator.
 * \tparam RAIter A random-access-iterator type.
 */
template<class RAIter>
class stride_iterator<RAIter,0> {
  public:
    typedef typename std::iterator_traits<RAIter>::value_type value_type;
    typedef typename std::iterator_traits<RAIter>::reference reference;
    typedef typename std::iterator_traits<RAIter>::difference_type difference_type;
    typedef typename std::iterator_traits<RAIter>::pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef stride_iterator<RAIter,0> self;
  private:
    RAIter pos; ///< Holds the underlying iterator.
    difference_type stride; ///< Holds the stride to apply to the iterator.
  public:
    /**
     * Default constructor.
     */
    stride_iterator<RAIter,0>( ) : pos(NULL), stride(0) { };
    /**
     * Copy-constructor.
     */
    stride_iterator<RAIter,0>(const self& rhs) : pos(rhs.pos), stride(rhs.stride) { };
    /**
     * Constructs the stride iterator at a given position.
     * \param aPos the position of the iterator.
     * \param aStride the stride to apply to the iterator.
     */
    stride_iterator<RAIter,0>(RAIter aPos, difference_type aStride) : pos(aPos), stride(aStride) { };

    /**
     * Returns the underlying iterator.
     * \return the underlying iterator.
     */
    RAIter base() const { return pos; };
    
    /**
     * Standard assignment operator.
     */
    self& operator=(const self& rhs) { pos = rhs.pos; stride = rhs.stride; };
    /**
     * Pre-increment operator.
     */
    self& operator++() { pos += stride; return *this; };
    /**
     * Post-increment operator.
     */
    self operator++(int) { self tmp = *this; pos += stride; return tmp; };
    /**
     * Add-and-store operator.
     */
    self& operator+=(difference_type aStep) { pos += aStep * stride; return *this; };
    /**
     * Pre-decrement operator.
     */
    self& operator--() { pos -= stride; return *this; };
    /**
     * Post-decrement operator.
     */
    self operator--(int) { self tmp = *this; pos -= stride; return tmp; };
    /**
     * Sub-and-store operator.
     */
    self& operator-=(difference_type aStep) { pos -= aStep * stride; return *this; };
    /**
     * Indexing operator.
     */
    reference operator[](difference_type aIdx) { return pos[aIdx * stride]; };
    /**
     * Dereference operator.
     */
    reference operator*() { return *pos; };
    /**
     * Member-access operator.
     */
    pointer operator->() { return pos.operator->(); };
    
    /**
     * Addition operator.
     */
    friend self operator+(self it1, difference_type n) { return it1 += n; };
    /**
     * Addition operator.
     */
    friend self operator+(difference_type n, self it) { return it += n; };
    
    /**
     * Subtraction operator.
     */
    friend difference_type operator-(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return (it1.pos - it2.pos) / it1.stride;
    };

    /**
     * Equality operator.
     */
    friend bool operator==(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos == it2.pos;
    };
    /**
     * Inequality operator.
     */
    friend bool operator!=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos != it2.pos;
    };
    /**
     * Less-than operator.
     */
    friend bool operator<(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos < it2.pos;
    };
    /**
     * Greater-than operator.
     */
    friend bool operator>(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos > it2.pos;
    };
    /**
     * Less-or-equal-than operator.
     */
    friend bool operator<=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos <= it2.pos;
    };
    /**
     * Greater-or-equal-than operator.
     */
    friend bool operator>=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos >= it2.pos;
    };
};

};


#endif
















