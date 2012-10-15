/**
 * \file mat_alg.hpp
 * 
 * This library declares (mostly by inclusions) all types of matrices currently available 
 * in the ReaK platform. These are all STL-vector based storage matrices or special matrices
 * such as nil or identity, or matrix adaptors, views, compositions and slices.
 * 
 * \todo Implement a suite of statically-sized matrix classes.
 * \todo Port the code related to the upper-triangular and lower-triangular matrices.
 * \todo Implement expression templates to optimize compound matrix expressions.
 * \todo Implement additional matrix views, for example, transposed view.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
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

#ifndef REAK_MAT_ALG_HPP
#define REAK_MAT_ALG_HPP

#include "mat_alg_general.hpp"
#include "mat_comparisons.hpp"
#include "mat_alg_rectangular.hpp"
#include "mat_alg_square.hpp"
#include "mat_alg_nil.hpp"
#include "mat_alg_identity.hpp"
#include "mat_alg_scalar.hpp"
#include "mat_alg_symmetric.hpp"
#include "mat_alg_skew_symmetric.hpp"
#include "mat_alg_diagonal.hpp"
#include "mat_alg_orthogonal.hpp"
#include "mat_alg_lower_triangular.hpp"
#include "mat_alg_upper_triangular.hpp"
#include "mat_alg_permutation.hpp"

#include "mat_operators.hpp"

#include "mat_vector_adaptor.hpp"
#include "mat_views.hpp"
#include "mat_transpose_view.hpp"
#include "mat_slices.hpp"
#include "mat_composite_adaptor.hpp"

/** Main namespace for ReaK */
namespace ReaK {
  


/****************************************************************************
                         Matrix Factory Functions
****************************************************************************/

/**
 * Builds a block diagonal matrix with two matrices of any dimension.
 * \param M1 first matrix (upper-left diagonal block).
 * \param M2 second matrix (lower-right diagonal block).
 * \return General block diagonal matrix.
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_writable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value,
Matrix1 >::type block_diag_mat(Matrix1 M1,const Matrix2& M2) {
  append_block_diag(M1,M2);
  return M1;
};

/**
 * Builds a block diagonal matrix with two matrices of any dimension.
 * \param M1 first matrix (upper-left diagonal block).
 * \param M2 second matrix (lower-right diagonal block).
 * \return General block diagonal matrix.
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_writable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
Matrix1 >::type block_diag_mat(Matrix1 M1,const Matrix2& M2,const Matrix3& M3) {
  append_block_diag(M1,M2);
  append_block_diag(M1,M3);
  return M1;
};

/**
 * Builds a block diagonal matrix with two matrices of any dimension.
 * \param M1 first matrix (upper-left diagonal block).
 * \param M2 second matrix (lower-right diagonal block).
 * \return General block diagonal matrix.
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_writable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value &&
                             is_readable_matrix<Matrix4>::value,
Matrix1 >::type block_diag_mat(Matrix1 M1,const Matrix2& M2,const Matrix3& M3,const Matrix4& M4) {
  append_block_diag(M1,M2);
  append_block_diag(M1,M3);
  append_block_diag(M1,M4);
  return M1;
};

/**
 * Builds a block diagonal matrix with two matrices of any dimension.
 * \param M1 first matrix (upper-left diagonal block).
 * \param M2 second matrix (lower-right diagonal block).
 * \return General block diagonal matrix.
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4, typename Matrix5>
typename boost::enable_if_c< is_writable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value &&
                             is_readable_matrix<Matrix4>::value &&
                             is_readable_matrix<Matrix5>::value,
Matrix1 >::type block_diag_mat(Matrix1 M1,const Matrix2& M2,const Matrix3& M3,const Matrix3& M4,const Matrix3& M5) {
  append_block_diag(M1,M2);
  append_block_diag(M1,M3);
  append_block_diag(M1,M4);
  append_block_diag(M1,M5);
  return M1;
};



/**
 * Builds a four-block matrix with four matrices of paired dimensions.
 * \param MUL first matrix (upper-left block).
 * \param MUR second matrix (upper-right block).
 * \param MLL third matrix (lower-left block).
 * \param MLR fourth matrix (lower-right block).
 * \return Compound four-block matrix.
 * \throw std::range_error if the dimensions of the four blocks don't allow proper juxtaposition.
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value &&
                             is_readable_matrix<Matrix4>::value,
Matrix1 >::type block_mat(Matrix1 MUL,const Matrix2& MUR,const Matrix3& MLL,const Matrix4& MLR) {
  if((MUL.get_row_count() != MUR.get_row_count()) || (MUL.get_col_count() != MLL.get_col_count()) || (MLL.get_row_count() != MLR.get_row_count()) || (MUR.get_col_count() != MLR.get_col_count()))
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType oldColCount = MUL.get_col_count();
  SizeType oldRowCount = MUL.get_row_count();
  append_block_diag(MUL,MLR);
  for(SizeType i = 0; i < MUR.get_row_count(); ++i)
    for(SizeType j = 0; j < MUR.get_col_count(); ++j)
      MUL(i,j + oldColCount) = MUR(i,j);
  for(SizeType i = 0; i < MLL.get_row_count(); ++i)
    for(SizeType j = 0; j < MLL.get_col_count(); ++j)
      MUL(i + oldRowCount,j) = MLL(i,j);
  return MUL;
};








};



#endif









