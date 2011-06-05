
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

#ifndef VECT_TRAITS_HPP
#define VECT_TRAITS_HPP


#include <boost/config.hpp>

namespace ReaK {

template <typename Vector>
struct vect_traits {
  typedef typename Vector::value_type value_type;
  typedef typename Vector::reference reference;
  typedef typename Vector::const_reference const_reference;
  typedef typename Vector::pointer pointer;
  typedef typename Vector::const_pointer const_pointer;
  typedef typename Vector::allocator_type allocator_type;
  
  typedef typename Vector::iterator iterator;
  typedef typename Vector::const_iterator const_iterator;
  
  typedef typename Vector::size_type size_type;
  typedef typename Vector::difference_type difference_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = Vector::dimensions);
    
};





};


#endif









