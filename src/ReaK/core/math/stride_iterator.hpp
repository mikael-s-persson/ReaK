
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

template<class RAIter, unsigned int Stride = 0>
class stride_iterator {
  private:
    RAIter pos;
  public:
    typedef typename std::iterator_traits<RAIter>::value_type value_type;
    typedef typename std::iterator_traits<RAIter>::reference reference;
    typedef typename std::iterator_traits<RAIter>::difference_type difference_type;
    typedef typename std::iterator_traits<RAIter>::pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef stride_iterator<RAIter,Stride> self;

    stride_iterator( ) : pos(NULL) { };
    stride_iterator(const self& rhs) : pos(rhs.pos) { };
    explicit stride_iterator(RAIter aPos) : pos(aPos) { };
    
    RAIter base() const { return pos; };

    self& operator=(const self& rhs) { pos = rhs.pos; };
    self& operator++() { pos += Stride; return *this; };
    self operator++(int) { self tmp = *this; pos += Stride; return tmp; };
    self& operator+=(difference_type aStep) { pos += aStep * Stride; return *this; };
    self& operator--() { pos -= Stride; return *this; };
    self operator--(int) { self tmp = *this; pos -= Stride; return tmp; };
    self& operator-=(difference_type aStep) { pos -= aStep * Stride; return *this; };
    reference operator[](difference_type aIdx) { return pos[aIdx * Stride]; };
    reference operator*() { return *pos; };
    pointer operator->() { return pos.operator->(); };
    
    friend self operator+(self it, difference_type n) { return it += n; };
    friend self operator+(difference_type n, self it) { return it += n; };
    
    friend difference_type operator-(const self& it1, const self& it2) { return (it1.pos - it2.pos) / Stride; };

    friend bool operator==(const self& it1, const self& it2) { return it1.pos == it2.pos; };
    friend bool operator!=(const self& it1, const self& it2) { return it1.pos != it2.pos; };
    friend bool operator<(const self& it1, const self& it2) { return it1.pos < it2.pos; };
    friend bool operator>(const self& it1, const self& it2) { return it1.pos > it2.pos; };
    friend bool operator<=(const self& it1, const self& it2) { return it1.pos <= it2.pos; };
    friend bool operator>=(const self& it1, const self& it2) { return it1.pos >= it2.pos; };
};

//Default template specialization, i.e., stride not given a compile time.
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
    RAIter pos;
    difference_type stride;
  public:
    stride_iterator<RAIter,0>( ) : pos(NULL), stride(0) { };
    stride_iterator<RAIter,0>(const self& rhs) : pos(rhs.pos), stride(rhs.stride) { };
    stride_iterator<RAIter,0>(RAIter aPos, difference_type aStride) : pos(aPos), stride(aStride) { };

    RAIter base() const { return pos; };
    
    self& operator=(const self& rhs) { pos = rhs.pos; stride = rhs.stride; };
    self& operator++() { pos += stride; return *this; };
    self operator++(int) { self tmp = *this; pos += stride; return tmp; };
    self& operator+=(difference_type aStep) { pos += aStep * stride; return *this; };
    self& operator--() { pos -= stride; return *this; };
    self operator--(int) { self tmp = *this; pos -= stride; return tmp; };
    self& operator-=(difference_type aStep) { pos -= aStep * stride; return *this; };
    reference operator[](difference_type aIdx) { return pos[aIdx * stride]; };
    reference operator*() { return *pos; };
    pointer operator->() { return pos.operator->(); };
    
    friend self operator+(self it1, difference_type n) { return it1 += n; };
    friend self operator+(difference_type n, self it) { return it += n; };
    
    friend difference_type operator-(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return (it1.pos - it2.pos) / it1.stride;
    };

    friend bool operator==(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos == it2.pos;
    };
    friend bool operator!=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos != it2.pos;
    };
    friend bool operator<(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos < it2.pos;
    };
    friend bool operator>(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos > it2.pos;
    };
    friend bool operator<=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos <= it2.pos;
    };
    friend bool operator>=(const self& it1, const self& it2) {
      if(it1.stride != it2.stride)
	throw std::range_error("Iterator stride mismatch.");
      return it1.pos >= it2.pos;
    };
};

};


#endif
















