/**
 * \file covariance_concept.hpp
 * 
 * This library provides a number of concepts and traits related to covariance matrix 
 * representations. A covariance matrix type is an abstract type that can provide a 
 * covariance matrix to described the uncertainty related to the components of a continuous 
 * vector. Since a covariance matrix can be represented in different ways, e.g., as a 
 * covariance matrix, as an information matrix, or decomposed, the algorithms in ReaK::ctrl
 * will not assume one representation, but rather relies on the concepts provided in this 
 * library to deal with covariance matrices.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_COVARIANCE_CONCEPT_HPP
#define REAK_COVARIANCE_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/core/lin_alg/mat_concepts.hpp>
#include <ReaK/core/lin_alg/arithmetic_tuple.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>

#include "state_vector_concept.hpp"

#include <boost/concept_check.hpp>


namespace ReaK {

namespace ctrl {
  
/** This namespace has a tag that tells the storage strategy for a covariance matrix. */
namespace covariance_storage {
  /** This tag tells the storage strategy for a covariance matrix. */
  enum tag {
    covariance = 1,
    information,
    decomposed,
    other
  };
};



/**
 * This traits class template defines the traits that characterize a covariance 
 * matrix (see CovarianceMatrixConcept).
 * \tparam CovarianceMatrix The covariance matrix type for which the traits are sought.
 */
template <typename CovarianceMatrix>
struct covariance_mat_traits {
  
  /** The type of the values of the components of the covariance matrix. */
  typedef typename CovarianceMatrix::value_type value_type;
  /** The type of the size of the covariance matrix. */
  typedef typename CovarianceMatrix::size_type size_type;
  
  /** The type of the actual covariance matrix that can be obtained with the required functions (see CovarianceMatrixConcept). */
  typedef typename CovarianceMatrix::matrix_type matrix_type;
  
  /** This constant tells the dimensions (or size) of the covariance matrix (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = CovarianceMatrix::dimensions);
  /** This constant tells the storage strategy of the covariance matrix (see covariance_storage::tag). */
  BOOST_STATIC_CONSTANT(covariance_storage::tag, storage = CovarianceMatrix::storage);
  
};


/** This namespace has a tag that gives the initial value of a covariance matrix. */
namespace covariance_initial_level {
  /** This tag gives the initial value of a covariance matrix. */
  enum tag {
    no_info = 0,
    full_info
  };
};


/**
 * This concept class template checks that a covariance matrix type models the concept 
 * of a covariance matrix as used in ReaK::ctrl. Basically, a covariance matrix can provide
 * a covariance matrix (actual matrix type modeling ReadableMatrixConcept) and an 
 * information matrix (actual matrix type modeling ReadableMatrixConcept) corresponding 
 * to the covariance matrix object. Also note that, in general, a regular matrix type will
 * not model this concept, see below.
 * 
 * Required concepts: 
 * 
 * the state-vector type should model the StateVectorConcept.
 * 
 * the matrix-type should model the ReadableMatrixConcept.
 * 
 * Valid expressions:
 * 
 * m = c.get_matrix();  The covariance matrix (actual matrix type) can be obtained.
 * 
 * m = c.get_inverse_matrix();  The information matrix (actual matrix type) can be obtained.
 *
 * s = dp * m * dp;  The matrix-type and state-difference type can be multiplied.
 * 
 * sz = c.size();  The covariance matrix type can provide its size.
 * 
 * \tparam CovarianceMatrix The covariance matrix type for which the traits are sought.
 * \tparam StateDiffType The state-difference type on which the covariance can apply.
 */
template <typename CovarianceMatrix, typename StateDiffType>
struct CovarianceMatrixConcept {
  
  typedef StateDiffType state_difference_type;
  
  BOOST_CONCEPT_ASSERT((ReadableMatrixConcept<typename covariance_mat_traits<CovarianceMatrix>::matrix_type>));
  
  CovarianceMatrix c;
  state_difference_type dp;
  typename covariance_mat_traits<CovarianceMatrix>::value_type s;
  typename covariance_mat_traits<CovarianceMatrix>::size_type sz;
  
  typename covariance_mat_traits<CovarianceMatrix>::matrix_type m;
  
  BOOST_CONCEPT_USAGE(CovarianceMatrixConcept)
  {
    typedef typename covariance_mat_traits<CovarianceMatrix>::value_type ValueType;
    using ReaK::to_vect;
    m = c.get_matrix();
    m = c.get_inverse_matrix();
    
    s = to_vect<ValueType>(dp) * m * to_vect<ValueType>(dp);
    sz = c.size();
  };
  
};



/**
 * This traits class template defines the traits that characterize a decomposed covariance 
 * matrix (see DecomposedCovarianceConcept).
 * \tparam CovarianceMatrix The covariance matrix type for which the traits are sought.
 */
template <typename CovarianceMatrix>
struct decomp_covariance_mat_traits {
  
  /** This type is the matrix type for the block components of the covariance matrix, should model ReadableMatrixConcept and WritableMatrixConcept. */
  typedef typename CovarianceMatrix::matrix_block_type matrix_block_type;
  
};


/**
 * This concept class template checks that a covariance matrix type models the concept 
 * of a decomposed covariance matrix as used in ReaK::ctrl. A common approach in estimation
 * is to decompose the covariance matrix P into a covarying block X and a informing-inverse 
 * block Y, such that the following relation applies: P = X * invert(Y). This is used, for 
 * example, to implement filters in header-file symplectic_kalman_filter.hpp 
 * or header-file aggregate_kalman_filter.hpp .
 * 
 * Required concepts: 
 * 
 * CovarianceMatrix should model the CovarianceMatrixConcept.
 * 
 * Valid expressions:
 * 
 * c = CovarianceMatrix(m,m);  A covariance matrix object can be created from the two blocks.
 * 
 * m = c.get_covarying_block();  The covarying block of the covariance matrix decomposition can be obtained.
 * 
 * m = c.get_informing_inv_block();  The informing-inverse block of the covariance matrix decomposition can be obtained.
 * 
 * \tparam CovarianceMatrix The covariance matrix type for which the traits are sought.
 * \tparam StateDiffType The state-difference type on which the covariance can apply.
 */
template <typename CovarianceMatrix, typename StateDiffType>
struct DecomposedCovarianceConcept : CovarianceMatrixConcept<CovarianceMatrix,StateDiffType> {
  
  typename decomp_covariance_mat_traits<CovarianceMatrix>::matrix_block_type mb;
  
  BOOST_CONCEPT_USAGE(DecomposedCovarianceConcept)
  {
    this->c = CovarianceMatrix(mb,mb);
    mb = this->c.get_covarying_block();
    mb = this->c.get_informing_inv_block();
  };
};




};

};

#endif













