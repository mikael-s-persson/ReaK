
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

#ifndef VECT_CONCEPTS_HPP
#define VECT_CONCEPTS_HPP

#include "vect_traits.hpp"
#include <vector>

#include <boost/concept_check.hpp>
#include <boost/config.hpp>

namespace ReaK {

template <typename Vector>
struct ReadableVectorConcept {
  Vector v;
  
  typename vect_traits<Vector>::size_type s;
  typename vect_traits<Vector>::value_type cr;
  
  typename vect_traits<Vector>::const_iterator it;
  
  void constraints() {
    cr = v[0]; //can be indexed and given an rvalue
    s = v.size();
    it = v.begin();
    ++it;
    it != v.end();
  };
  
};

template <typename Vector>
struct is_readable_vector {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_readable_vector<Vector> type;
};

template <typename T>
struct is_readable_vector< std::vector<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< std::vector<T> > type;
};



template <typename Vector>
struct WritableVectorConcept 
  : ReadableVectorConcept<Vector> { //must also be readable.
  
  typename vect_traits<Vector>::value_type r;
  
  void constraints() {
    this->v[0] = r; //can be indexed and given an lvalue
  };
  
};


template <typename Vector>
struct is_writable_vector {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector<Vector> type;
};

template <typename T>
struct is_writable_vector< std::vector<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< std::vector<T> > type;
};



template <typename Vector>
struct ResizableVectorConcept { 
  Vector v;
  
  typename vect_traits<Vector>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableVectorConcept) {
    v.resize(sz);
  };
  
};


template <typename Vector>
struct is_resizable_vector {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector<Vector> type;
};

template <typename T>
struct is_resizable_vector< std::vector<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_vector< std::vector<T> > type;
};




template <typename Vector>
struct DynAllocVectorConcept : ResizableVectorConcept<Vector> { 
  Vector v;
  
  typename vect_traits<Vector>::allocator_type al;
  
  BOOST_CONCEPT_USAGE(DynAllocVectorConcept) {
    al = v.get_allocator();
  };
  
};


template <typename Vector>
struct has_allocator_vector {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector<Vector> type;
};

template <typename T>
struct has_allocator_vector< std::vector<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_vector< std::vector<T> > type;
};


};



#endif


