/**
 * \file vect_concepts.hpp
 * 
 * This library defines the concepts related to generic vector types used in the ReaK library.
 * These concepts are based on the principle of minimum requirement, they were designed to 
 * require only the minimum set of valid expressions that will be used by the algorithms
 * that pertain to generic vector types.
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

#ifndef REAK_VECT_CONCEPTS_HPP
#define REAK_VECT_CONCEPTS_HPP

#include "vect_traits.hpp"
#include <vector>

#include <boost/concept_check.hpp>
#include <boost/config.hpp>

namespace ReaK {

/**
 * This concept class defines what makes a vector a readable vector, that is, 
 * a vector whose elements can be read.
 * 
 * Valid Expressions:
 * 
 * e = v[i];   can be indexed to be read.
 * 
 * s = v.size();   the size of the vector can be obtained.
 * 
 * cit = v.begin();   a const-iterator to the first vector element can be obtained.
 * 
 * ++cit;   the const-iterator can be incremented.
 * 
 * cit = v.end();   a const-iterator to the one-past-last vector element can be obtained.
 * 
 * \tparam Vector The vector type.
 */
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
    bool b = (it != v.end()); RK_UNUSED(b);
  };
  
};

/**
 * This meta-function evaluates whether a Vector class fulfills the ReadableVectorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
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



/**
 * This concept class defines what makes a vector a writable vector, that is, 
 * a vector whose elements can be written.
 * 
 * Valid Expressions (in addition to those of ReadableVectorConcept):
 * 
 * v[i] = e;   can be indexed to be written.
 * 
 * \tparam Vector The vector type.
 */
template <typename Vector>
struct WritableVectorConcept 
  : ReadableVectorConcept<Vector> { //must also be readable.
  
  typename vect_traits<Vector>::value_type r;
  
  void constraints() {
    this->v[0] = r; //can be indexed and given an lvalue
  };
  
};


/**
 * This meta-function evaluates whether a Vector class fulfills the WritableVectorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
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



/**
 * This concept class defines what makes a vector a resizable vector, that is, 
 * a vector whose size can be changed at run-time.
 * 
 * Valid Expressions:
 * 
 * v.resize(s);   can be resized.
 * 
 * \tparam Vector The vector type.
 */
template <typename Vector>
struct ResizableVectorConcept { 
  Vector v;
  
  typename vect_traits<Vector>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableVectorConcept) {
    v.resize(sz);
  };
  
};


/**
 * This meta-function evaluates whether a Vector class fulfills the ResizableVectorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
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




/**
 * This concept class defines what makes a vector a resizable vector, that is, 
 * a vector whose size can be changed at run-time.
 * 
 * Valid Expressions (in addition to those of ResizableVectorConcept):
 * 
 * al = v.get_allocator();   the allocator object can be obtained.
 * 
 * \tparam Vector The vector type.
 */
template <typename Vector>
struct DynAllocVectorConcept : ResizableVectorConcept<Vector> { 
  Vector v;
  
  typename vect_traits<Vector>::allocator_type al;
  
  BOOST_CONCEPT_USAGE(DynAllocVectorConcept) {
    al = v.get_allocator();
  };
  
};

/**
 * This meta-function evaluates whether a Vector class fulfills the DynAllocVectorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a vector class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new vector class.
 */
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


