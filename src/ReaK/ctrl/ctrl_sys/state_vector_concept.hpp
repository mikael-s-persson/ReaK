/**
 * \file state_vector_concept.hpp
 * 
 * This library provides the traits class and the concept definition for a 
 * state-vector as used within the ReaK::ctrl namespace. A state-vector is 
 * an abstraction of the quantity that describes the state of a state-space
 * system (see SSSystemConcept).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_STATE_VECTOR_CONCEPT_HPP
#define REAK_STATE_VECTOR_CONCEPT_HPP

#include <ReaK/core/lin_alg/vect_concepts.hpp>
#include <ReaK/core/lin_alg/vect_alg.hpp>
#include <ReaK/core/lin_alg/mat_slices.hpp>

#include <boost/concept_check.hpp>


namespace ReaK {

namespace ctrl {

/**
 * This class template is the traits class that defines the traits that a 
 * state-vector should have.
 * \tparam StateVector The state-vector type for which the traits are sought.
 */
template <typename StateVector>
struct state_vector_traits {
  /** This is the type of the state-vector descriptor, usually the same as StateVector. */
  typedef typename StateVector::state_type state_type; 
  /** This is the type that describes the difference between two state-vectors. */
  typedef typename StateVector::state_difference_type state_difference_type;
  /** This is the value-type of the elements of the state-vector. */
  typedef typename StateVector::value_type value_type;
  /** This is the type that describes the size of the state-vector. */
  typedef typename StateVector::size_type size_type;

  /** This constant describes the dimension of the state-vector (0 if only known at run-time). */
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = StateVector::dimensions);

};




/**
 * This class template defines the concept that state-vectors should model.
 * 
 * Required concepts: 
 * 
 * the state-vector's state_difference_type should model ReadableVectorConcept.
 * 
 * Valid expressions (state_type s, state_difference_type ds, value_type v, size_type sz):
 * 
 * ds = diff(s,s);  The difference between state-vectors is obtained by the diff() function.
 * 
 * s = add(s,ds);  A state-difference can be added to a state-vector with the add() function.
 * 
 * ds = v * ds;  A state-difference is scalable.
 * 
 * ds = ds + ds;  State-differences can be added.
 * 
 * ds += ds;  State-differences can be added and stored.
 * 
 * ds = ds - ds;  State-differences can be subtracted.
 * 
 * ds -= ds;  State-differences can be subtracted and stored.
 * 
 * ds = -ds;  A state-difference can be negated.
 * 
 * ds = unit(ds);  A state-difference can be made into a unit-vector.
 * 
 * v = norm(ds);  A state-difference can be taken the norm of.
 * 
 * sz = ds.size();  A state-difference has a size.
 * 
 * \tparam StateVector The state-vector type to test for modeling the state-vector concept.
 */
template <typename StateVector>
struct StateVectorConcept {
  typename state_vector_traits<StateVector>::state_type s;
  typename state_vector_traits<StateVector>::state_difference_type ds;
  typename state_vector_traits<StateVector>::value_type v;
  typename state_vector_traits<StateVector>::size_type sz;
  
  BOOST_CONCEPT_ASSERT((ReadableVectorConcept< typename state_vector_traits<StateVector>::state_difference_type >));
  
  BOOST_CONCEPT_USAGE(StateVectorConcept)
  {
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
    v = norm_2(ds);
    sz = ds.size();
  };
  
};


/**
 * This meta-function is used to evaluate whether a type models the StateVectorConcept.
 * This does not attempt to instantiate the StateVectorConcept template because it would
 * break Sfinae rules. Instead, if one wants to make a new state-vector class, he should 
 * specialize this meta-function to evaluate to true.
 * \tparam StateVector The type which may or may not be a state-vector type.
 */
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














