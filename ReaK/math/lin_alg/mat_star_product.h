/**
 * \file mat_star_product.h
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

#ifndef REAK_MATH_LIN_ALG_MAT_STAR_PRODUCT_H_
#define REAK_MATH_LIN_ALG_MAT_STAR_PRODUCT_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_composite_adaptor.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"

#include <type_traits>

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
  using upper_left = mat<ValueType, mat_structure::square>;
  using upper_right = mat<ValueType, mat_structure::symmetric>;
  using lower_left = mat<ValueType, mat_structure::symmetric>;
  using lower_right = mat<ValueType, mat_structure::square>;

  using upper = mat_horiz_cat<upper_left, upper_right>;
  using lower = mat_horiz_cat<lower_left, lower_right>;

  using type = mat_vert_cat<upper, lower>;
};

template <typename ValueType>
using hamiltonian_mat_t = typename hamiltonian_mat<ValueType>::type;

/**
 * This function template computes the Redeffer star-product of two matrices.
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \param M1 A square, hamiltonian matrix.
 * \param M2 A square, hamiltonian matrix.
 * \return A hamiltonian matrix which is the star product of M1 and M2.
 */
template <typename Matrix1, typename Matrix2>
auto star_product(const Matrix1& M1, const Matrix2& M2) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_readable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;

  int N = M1.get_row_count() / 2;
  if ((M1.get_col_count() != 2 * N) || (M2.get_col_count() != 2 * N) ||
      (M2.get_row_count() != 2 * N)) {
    throw std::range_error("Matrix dimensions mismatch.");
  }

  mat<ValueType, mat_structure::square> Inv1(N);
  pseudoinvert_SVD(mat<ValueType, mat_structure::identity>(N) -
                       mat_const_sub_block(M2, N, N, N, 0) *
                           mat_const_sub_block(M1, N, N, 0, N),
                   Inv1);

  mat<ValueType, mat_structure::square> DInv(
      mat_const_sub_block(M1, N, N, N, N) * Inv1);
  mat<ValueType, mat_structure::square> WInvt(
      mat_const_sub_block(M2, N, N, 0, 0) * transpose_move(Inv1));

  return hamiltonian_mat_t<ValueType>(
      mat_horiz_cat(mat<ValueType, mat_structure::square>(
                        WInvt * mat_const_sub_block(M1, N, N, 0, 0)),
                    mat<ValueType, mat_structure::symmetric>(
                        mat_const_sub_block(M2, N, N, 0, N) +
                        WInvt * mat_const_sub_block(M1, N, N, 0, N) *
                            mat_const_sub_block(M2, N, N, N, N))),
      mat_horiz_cat(mat<ValueType, mat_structure::symmetric>(
                        mat_const_sub_block(M1, N, N, N, 0) +
                        DInv * mat_const_sub_block(M2, N, N, N, 0) *
                            mat_const_sub_block(M1, N, N, 0, 0)),
                    mat<ValueType, mat_structure::square>(
                        DInv * mat_const_sub_block(M2, N, N, N, N))));
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_STAR_PRODUCT_H_
