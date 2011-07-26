
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

#ifndef VECT_INDEX_ITERATOR_HPP
#define VECT_INDEX_ITERATOR_HPP

#include "vect_traits.hpp"

#include <iterator>
#include <stdexcept>

namespace ReaK {


template <typename Vector>
class vect_index_iter {
  public:
    typedef typename vect_traits<Vector>::value_type value_type;
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::difference_type difference_type;
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef vect_index_iter<Vector> self;
    
  private:
    Vector* v;
    size_type i;
    
  public:
    vect_index_iter(Vector& aV) : v(&aV), i(aV.size()) { };
    vect_index_iter(Vector& aV, size_type aI) : v(&aV), i(aI) { };
    
    self& operator++() { ++i; return *this; };
    self operator++(int) { self tmp = *this; ++i; return tmp; };
    self& operator+=(difference_type aStep) { i += aStep; return *this; };
    self& operator--() { --i; return *this; };
    self operator--(int) { self tmp = *this; --i; return tmp; };
    self& operator-=(difference_type aStep) { i -= aStep; return *this; };
    reference operator[](difference_type aIdx) const { return (*v)[aIdx]; };
    reference operator*() const { return (*v)[i]; };
    pointer operator->() const { return &(*v)[i]; };
    
    friend self operator+(self it, difference_type n) { return it += n; };
    friend self operator+(difference_type n, self it) { return it += n; };
    
    friend difference_type operator-(const self& it1, const self& it2) { return it1.i - it2.i; };

    friend bool operator==(const self& it1, const self& it2) { return it1.i == it2.i; };
    friend bool operator!=(const self& it1, const self& it2) { return it1.i != it2.i; };
    friend bool operator< (const self& it1, const self& it2) { return it1.i <  it2.i; };
    friend bool operator> (const self& it1, const self& it2) { return it1.i >  it2.i; };
    friend bool operator<=(const self& it1, const self& it2) { return it1.i <= it2.i; };
    friend bool operator>=(const self& it1, const self& it2) { return it1.i >= it2.i; };

};



template <typename Vector>
class vect_index_const_iter {
  public:
    typedef typename vect_traits<Vector>::value_type value_type;
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::const_reference reference;
    typedef typename vect_traits<Vector>::difference_type difference_type;
    typedef typename vect_traits<Vector>::const_pointer pointer;
    typedef std::random_access_iterator_tag iterator_category;
    typedef vect_index_const_iter<Vector> self;
    
  private:
    const Vector* v;
    size_type i;
    
  public:
    vect_index_iter(Vector& aV) : v(&aV), i(aV.size()) { };
    vect_index_iter(Vector& aV, size_type aI) : v(&aV), i(aI) { };
    
    self& operator++() { ++i; return *this; };
    self operator++(int) { self tmp = *this; ++i; return tmp; };
    self& operator+=(difference_type aStep) { i += aStep; return *this; };
    self& operator--() { --i; return *this; };
    self operator--(int) { self tmp = *this; --i; return tmp; };
    self& operator-=(difference_type aStep) { i -= aStep; return *this; };
    reference operator[](difference_type aIdx) const { return (*v)[aIdx]; };
    reference operator*() const { return (*v)[i]; };
    pointer operator->() const { return &(*v)[i]; };
    
    friend self operator+(self it, difference_type n) { return it += n; };
    friend self operator+(difference_type n, self it) { return it += n; };
    
    friend difference_type operator-(const self& it1, const self& it2) { return it1.i - it2.i; };

    friend bool operator==(const self& it1, const self& it2) { return it1.i == it2.i; };
    friend bool operator!=(const self& it1, const self& it2) { return it1.i != it2.i; };
    friend bool operator< (const self& it1, const self& it2) { return it1.i <  it2.i; };
    friend bool operator> (const self& it1, const self& it2) { return it1.i >  it2.i; };
    friend bool operator<=(const self& it1, const self& it2) { return it1.i <= it2.i; };
    friend bool operator>=(const self& it1, const self& it2) { return it1.i >= it2.i; };

};






};

#endif











