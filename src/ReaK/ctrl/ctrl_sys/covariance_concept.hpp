
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

#ifndef COVARIANCE_CONCEPT_HPP
#define COVARIANCE_CONCEPT_HPP

#include "math/mat_concepts.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace ctrl {
  

namespace covariance_storage {
  enum tag {
    covariance = 1,
    information,
    decomposed,
    other
  };
};



template <typename CovarianceMatrix>
struct covariance_mat_traits {
  
  typedef CovarianceMatrix::point_type point_type;
  typedef CovarianceMatrix::point_difference_type point_difference_type;
  
  typedef CovarianceMatrix::value_type value_type;
  typedef CovarianceMatrix::size_type size_type;
  
  typedef CovarianceMatrix::matrix_type matrix_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = CovarianceMatrix::dimensions);
  BOOST_STATIC_CONSTANT(covariance_storage::tag, storage = CovarianceMatrix::storage);
  
};


enum covariance_initial_level {
  no_info = 0,
  full_info
};


template <typename CovarianceMatrix>
struct CovarianceMatrixConcept {
  CovarianceMatrix c;
  
  typename covariance_mat_traits<CovarianceMatrix>::point_type p;
  typename covariance_mat_traits<CovarianceMatrix>::point_difference_type dp;
  typename covariance_mat_traits<CovarianceMatrix>::value_type s;
  typename covariance_mat_traits<CovarianceMatrix>::size_type sz;
  
  typename covariance_mat_traits<CovarianceMatrix>::matrix_type m;
  
  void constraints() {
    boost::function_requires< ReadableMatrixConcept< typename covariance_mat_traits<CovarianceMatrix>::matrix_type > >();

    dp = p - p;
    dp = -dp;
    dp = unit(dp);
    
    m = c.get_matrix();
    m = c.get_inverse_matrix();
    
    s = dp * m * dp;
    sz = c.size();
  };
  
};



template <typename CovarianceMatrix>
struct decomp_covariance_mat_traits {
  
  typedef CovarianceMatrix::matrix_block_type matrix_block_type;
  
};


template <typename CovarianceMatrix>
struct DecomposedCovarianceConcept {
  CovarianceMatrix c;
  
  typename decomp_covariance_mat_traits<CovarianceMatrix>::matrix_block_type m;
  
  void constraints() {
    boost::function_requires< CovarianceMatrixConcept<CovarianceMatrix> >();
    c = CovarianceMatrix(m,m);
    m = c.get_covarying_block();
    m = c.get_informing_inv_block();
  };
};




};

};

#endif













