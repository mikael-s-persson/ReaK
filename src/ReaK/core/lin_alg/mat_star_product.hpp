/**
 * \file mat_star_product.hpp
 * 
 * This library provides a function to compute the Redeffer Star-product and the 
 * definition of the hamiltonian matrix type alias. The Redeffer star product can 
 * be applied to aggregate hamiltonian matrices together. This is used in hamiltonian
 * maps (i.e. diffusion operators).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_MAT_STAR_PRODUCT_HPP
#define REAK_MAT_STAR_PRODUCT_HPP

#include "mat_alg.hpp"
#include "mat_composite_adaptor.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_svd_method.hpp"

namespace ReaK {

/**
 * This class template is a type alias for what can represent a hamiltonian matrix (not strictly enforced).
 * 
 * Models: ReadableMatrixConcept and WritableMatrixConcept.
 * 
 * \tparam ValueType The value-type of the underlying matrices.
 */
template <typename ValueType>
struct hamiltonian_mat {
  typedef mat<ValueType, mat_structure::square> upper_left;
  typedef mat<ValueType, mat_structure::symmetric> upper_right;
  typedef mat<ValueType, mat_structure::symmetric> lower_left;
  typedef mat<ValueType, mat_structure::square> lower_right;
  
  typedef mat_horiz_cat< upper_left, upper_right > upper;
  typedef mat_horiz_cat< lower_left, lower_right > lower;
  
  typedef mat_vert_cat< upper, lower > type;
};


/**
 * This function template computes the Redeffer star-product of two matrices.
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \param M1 A square, hamiltonian matrix.
 * \param M2 A square, hamiltonian matrix.
 * \return A hamiltonian matrix which is the star product of M1 and M2.
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value,
typename hamiltonian_mat< typename mat_traits<Matrix1>::value_type >::type >::type 
 star_product(const Matrix1& M1, const Matrix2& M2) {
  typedef typename hamiltonian_mat<ValueType>::type result_type;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<result_type>::size_type SizeType;
  
  SizeType N = M1.get_row_count() / 2;
  if((M1.get_col_count() != 2 * N) || (M2.get_col_count() != 2 * N) || (M2.get_row_count() != 2 * N))
    throw std::range_error("Matrix dimensions mismatch.");
  
  mat<ValueType,mat_structure::square> Inv1(N);
  pseudoinvert_SVD(mat<ValueType,mat_structure::identity>(N)
                 - mat_const_sub_block< Matrix2 >(M2,N,N,N,0)
                 * mat_const_sub_block< Matrix1 >(M1,N,N,0,N),Inv1);
  
  mat<ValueType,mat_structure::square> DInv( mat_const_sub_block< Matrix1 >(M1,N,N,N,N) * Inv1 );
  mat<ValueType,mat_structure::square> WInvt( mat_const_sub_block< Matrix2 >(M2,N,N,0,0) * transpose_move(Inv1) );
  
  return result_type( 
    mat_horiz_cat< mat<ValueType,mat_structure::square>,
                   mat<ValueType,mat_structure::symmetric> >(
                     mat<ValueType,mat_structure::square>(WInvt * mat_const_sub_block<Matrix1>(M1,N,N,0,0)),
                     mat<ValueType,mat_structure::symmetric>(mat_const_sub_block<Matrix2>(M2,N,N,0,N)
                                                           + WInvt 
                                                           * mat_const_sub_block<Matrix1>(M1,N,N,0,N)
                                                           * mat_const_sub_block<Matrix2>(M2,N,N,N,N) ) ),
    mat_horiz_cat< mat<ValueType,mat_structure::symmetric>,
                   mat<ValueType,mat_structure::square> >(
                     mat<ValueType,mat_structure::symmetric>(mat_const_sub_block<Matrix1>(M1,N,N,N,0)
                                                           + DInv 
                                                           * mat_const_sub_block<Matrix2>(M2,N,N,N,0)
                                                           * mat_const_sub_block<Matrix1>(M1,N,N,0,0) ),
                     mat<ValueType,mat_structure::square>(DInv * mat_const_sub_block<Matrix2>(M2,N,N,N,N) ) ) );
};



};

#endif




