
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

#ifndef STATE_VECTOR_CONCEPT_HPP
#define STATE_VECTOR_CONCEPT_HPP

#include "math/vect_concepts.hpp"
#include "math/vect_alg.hpp"
#include "math/mat_slices.hpp"

#include <boost/concept_check.hpp>


namespace ReaK {

namespace ctrl {


template <typename StateVector>
struct state_vector_traits {
  typedef typename StateVector::state_type state_type;
  typedef typename StateVector::state_difference_type state_difference_type;
  typedef typename StateVector::value_type value_type;
  typedef typename StateVector::size_type size_type;

  BOOST_STATIC_CONSTANT(std::size_t, dimensions = StateVector::dimensions);

};





template <typename StateVector>
struct StateVectorConcept {
  typename state_vector_traits<StateVector>::state_type s;
  typename state_vector_traits<StateVector>::state_difference_type ds;
  typename state_vector_traits<StateVector>::value_type v;
  typename state_vector_traits<StateVector>::size_type sz;
  
  void constraints() {
    boost::function_requires< ReadableVectorConcept< typename state_vector_traits<StateVector>::state_type > >();
    boost::function_requires< WritableVectorConcept< typename state_vector_traits<StateVector>::state_difference_type > >();

    ds = diff(s,s);
    s = add(s,ds);
    ds = v * ds;
    ds *= v;
    ds = ds + ds;
    ds += ds;
    ds = ds - ds;
    ds -= ds;
    ds = -ds;
    ds = unit(ds);
    v = norm(ds);
    sz = ds.size();
  };
  
};


template <typename StateVector>
struct is_state_vector {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_state_vector<StateVector> type;
};



template <typename T, typename Allocator>
struct state_vector_traits< vect_n<T,Allocator> > {
  typedef vect_n<T,Allocator> state_type;
  typedef vect_n<T,Allocator> state_difference_type;
  typedef typename vect_traits< vect_n<T,Allocator> >::value_type value_type;
  typedef typename vect_traits< vect_n<T,Allocator> >::size_type size_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = vect_traits< state_type >::dimensions);
  
};

template <typename T, typename Allocator>
struct is_state_vector< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_state_vector< vect_n<T,Allocator> > type;
};



template <typename T, unsigned int Size>
struct state_vector_traits< vect<T,Size> > {
  typedef vect<T,Size> state_type;
  typedef vect<T,Size> state_difference_type;
  typedef typename vect_traits< vect<T,Size> >::value_type value_type;
  typedef typename vect_traits< vect<T,Size> >::size_type size_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = vect_traits< state_type >::dimensions);
  
};

template <typename T, unsigned int Size>
struct is_state_vector< vect<T,Size> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_state_vector< vect<T,Size> > type;
};


};

};

#endif














