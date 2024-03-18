/**
 * \file mat_ctrl_decomp.h
 *
 * This library provides function templates to perform various control-theoretic decompositions or
 * reductions, or state-space arrangements.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_MATH_LIN_ALG_MAT_CTRL_DECOMP_H_
#define REAK_MATH_LIN_ALG_MAT_CTRL_DECOMP_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/math/lin_alg/mat_householder.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/math/lin_alg/mat_norms.h"

#include <type_traits>

namespace ReaK {

namespace detail {}  // namespace detail

/**
 * This function will reduce a state-space representation of the state transition (A,B) into
 * the controllable and uncontrollable parts. The resulting reduction yields:
 * \n
 * $A = Q Ar Q^T$
 * \n
 * $B = Q Br Z^T$
 * \n
 * where Br has non-zero elements on the first r elements, and Ar has a N-r by N-r block of zeros at the bottom-left
 *corner.
 * The rank r is the order of the controllable space (invariant subspace) of the system (i.e., how many states are
 *controllable).
 * \n
 *
 * \param A square (n x n) matrix which represents state-to-state-derivative linear map (or state-to-state for a
 *discrete-time system).
 * \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map (or
 *input-to-state-increment for a discrete-time system).
 * \param Q accumulates a square (n x n) orthogonal transformation with the controllable / uncontrollable basis vectors.
 * \param Z accumulates a square (m x m) orthogonal transformation of the input vector (only for reordering
 *(permutations)).
 * \param NumTol tolerance for considering a value to be zero.
 * \return The numerical rank of the system, i.e., the number of controllable states.
 *
 * \throws std::range_error if the matrix dimensions are not consistent.
 *
 * \author Mikael Persson
 */
template <FullyWritableMatrix Matrix1, FullyWritableMatrix Matrix2, FullyWritableMatrix Matrix3,
          FullyWritableMatrix Matrix4>
int ctrl_reduction(Matrix1& A, Matrix2& B, Matrix3& Q, Matrix4& Z,
                   mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count())) {
    throw std::range_error(
        "The dimensions of the system matrices do not match! Should be A(n x "
        "n) and B(n x m).");
  }

  if ((A.get_row_count() == 0) || (B.get_col_count() == 0)) {
    return 0;
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  using std::swap;

  int N = A.get_row_count();
  int M = B.get_col_count();

  if ((N != Q.get_col_count()) || (Q.get_row_count() != Q.get_col_count())) {
    Q = mat<ValueType, mat_structure::identity>(N);
  }

  if ((M != Z.get_col_count()) || (Z.get_row_count() != Z.get_col_count())) {
    Z = mat<ValueType, mat_structure::identity>(N);
  }

  mat<ValueType, mat_structure::square> Tr(
      (mat<ValueType, mat_structure::identity>(N)));
  mat<ValueType, mat_structure::permutation> PCr(M);
  int r = detail::decompose_RRQR_impl(B, &Tr, PCr, NumTol);
  A = transpose_view(Tr) * A * Tr;
  Q = Q * Tr;
  Z = Z * PCr;
  if (r == N) {
    return N;
  }
  mat<ValueType, mat_structure::rectangular> AB_accum(B);
  for (int i = 1; i < N; ++i) {
    AB_accum =
        sub(A)(range(r, N), range(N - AB_accum.get_row_count(), N)) * AB_accum;
    mat<ValueType, mat_structure::square> Tr2(
        (mat<ValueType, mat_structure::identity>(N)));
    mat<ValueType, mat_structure::permutation> PCr2(M);
    mat_sub_block<mat<ValueType, mat_structure::square>> Tr2_sub(Tr2, N - r,
                                                                 N - r, r, r);
    r += detail::decompose_RRQR_impl(AB_accum, &Tr2_sub, PCr2, NumTol);
    A = transpose_view(Tr2) * A * Tr2;
    Q = Q * Tr2;
    if (r == N) {
      return N;
    }
  }
  return r;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_CTRL_DECOMP_H_
