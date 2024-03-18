/**
 * \file mat_hess_decomp.h
 *
 * This library provides methods to perform the Upper-Hessenberg decomposition of a matrix as
 * well as the Hessenberg-Triangular reduction of two matrices. These algorithms form the first
 * step in many algorithms such as the QR algorithm, the QZ algorithm, the Real-Schur decomposition,
 * etc. These algorithms were implemented as described in Golub and van Loan's classic book.
 *
 * This library simply provides various versions of the same algorithm (same underlying implementation,
 * with different interfaces). Versions differ based mostly on the type of matrices fed to the function
 * overloads and whether the orthogonal matrices that achieve the decompositions are required or not
 * (since explicitly forming those matrices can be expensive and should be avoided if it's not needed).
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

#ifndef REAK_MATH_LIN_ALG_MAT_HESS_DECOMP_H_
#define REAK_MATH_LIN_ALG_MAT_HESS_DECOMP_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/math/lin_alg/mat_givens_rot.h"
#include "ReaK/math/lin_alg/mat_householder.h"

#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/core/base/shared_object.h"

#include <type_traits>

namespace ReaK {

/*************************************************************************
                          Hessenberg Decompositions
*************************************************************************/

namespace detail {

template <typename Matrix1, typename Matrix2>
void decompose_Hess_impl(Matrix1& H, Matrix2* Q,
                         mat_value_type_t<Matrix1> absNumTol) {

  using ValueType = mat_value_type_t<Matrix1>;
  int N = H.get_row_count();
  householder_matrix<vect_n<ValueType>> hhm;

  for (int i = 0; i + 2 < N; ++i) {

    hhm.set(mat_row_slice<Matrix1>(H, i, i + 1, N - i - 1), absNumTol);

    mat_sub_block<Matrix1> subH1(H, N - i - 1, N - i, i + 1, i);
    householder_prod(hhm, subH1);  // P * H

    mat_sub_block<Matrix1> subH2(H, N, N - i - 1, 0, i + 1);
    householder_prod(subH2, hhm);  // H * P

    if (Q) {
      mat_sub_block<Matrix2> subQ(*Q, N, N - i - 1, 0, i + 1);
      householder_prod(subQ, hhm);  // Q_prev * P
    }
  }
}

/* tested and working. Golub and vanLoan Alg.-8.3.1 */
template <typename Matrix1, typename Matrix2>
void decompose_TriDiag_impl(Matrix1& T, Matrix2* Q,
                            mat_value_type_t<Matrix1> absNumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::sqrt;
  int N = T.get_row_count();
  householder_matrix<vect_n<ValueType>> hhm;

  vect_n<ValueType> p(N);
  vect_n<ValueType> w(N);

  for (int i = 0; i + 2 < N; ++i) {

    hhm.set(mat_row_slice<Matrix1>(T, i, i + 1, N - i - 1), absNumTol);

    mat_sub_block<Matrix1> subT1(T, N - i, N - i - 1, i, i + 1);
    householder_prod(subT1, hhm);  // Q_prev * P

    mat_sub_block<Matrix1> subT2(T, N - i - 1, N - i, i + 1, i);
    householder_prod(hhm, subT2);  // Q_prev * P

    if (Q) {
      // mat_sub_block<Matrix2> subQ(*Q,N - i - 1,N,i + 1,0);
      // householder_prod(hhm,subQ); // Q_prev * P
      mat_sub_block<Matrix2> subQ(*Q, N, N - i - 1, 0, i + 1);
      householder_prod(subQ, hhm);  // Q_prev * P
    }
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void reduce_HessTri_offset_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                                mat_size_type_t<Matrix1> row_offset,
                                mat_value_type_t<Matrix1> absNumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;

  int N = A.get_row_count() - row_offset;

  givens_rot_matrix<ValueType> G;

  householder_matrix<vect_n<ValueType>> hhm;

  for (int i = 0; i + 1 < N; ++i) {

    hhm.set(mat_row_slice<Matrix2>(B, i, row_offset + i, N - i), absNumTol);

    mat_sub_block<Matrix2> subB(B, N - i, B.get_col_count() - i, row_offset + i,
                                i);
    householder_prod(hhm, subB);  // P * R

    mat_sub_block<Matrix1> subA(A, N - i, A.get_col_count(), row_offset + i, 0);
    householder_prod(hhm, subA);  // P * R

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), N - i, 0, i);
      householder_prod(subQ, hhm);  // Q_prev * P
    }
  }

  for (int j = 0; j + 2 < N; ++j) {
    for (int i = N - 1; i > j + 1; --i) {
      if (abs(A(i, j)) < absNumTol) {
        continue;
      }

      G.set(A(row_offset + i - 1, j), A(row_offset + i, j));

      mat_sub_block<Matrix1> subA1(A, 2, A.get_col_count() - j,
                                   row_offset + i - 1, j);
      givens_rot_prod(G, subA1);  // G * A

      mat_sub_block<Matrix2> subB1(B, 2, B.get_col_count() - i + 1,
                                   row_offset + i - 1, i - 1);
      givens_rot_prod(G, subB1);  // G * B

      if (Q) {
        mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, i - 1);
        givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
      }

      G.set(-B(row_offset + i, i), B(row_offset + i, i - 1));
      G = transpose(G);

      mat_sub_block<Matrix2> subB2(B, i + 1, 2, row_offset, i - 1);
      givens_rot_prod(subB2, G);  // B * G^T

      mat_sub_block<Matrix1> subA2(A, N, 2, row_offset, i - 1);
      givens_rot_prod(subA2, G);  // A * G^T

      if (Z) {
        mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, i - 1);
        givens_rot_prod(subZ, G);  // Q_prev * G^T
      }
    }
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void reduce_HessTri_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                         mat_value_type_t<Matrix1> absNumTol) {
  reduce_HessTri_offset_impl(A, B, Q, Z, 0, absNumTol);
}
}  // namespace detail

/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \param A square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2, WritableMatrix Matrix3>
void decompose_Hess(const Matrix1& A, Matrix2& Q, Matrix3& H,
                    mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Upper-Hessenberg decomposition is only possible on a square matrix!");
  }

  auto Qtmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix2>) {
      Q = mat<mat_value_type_t<Matrix2>, mat_structure::identity>(
          A.get_row_count());
      return std::unique_ptr<Matrix2, null_deleter>(&Q, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix2>, mat_structure::square>>(
          mat<mat_value_type_t<Matrix2>, mat_structure::identity>(
              A.get_row_count()));
    }
  }();
  auto Htmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix3>) {
      H = A;
      return std::unique_ptr<Matrix3, null_deleter>(&H, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix3>, mat_structure::square>>(A);
    }
  }();

  detail::decompose_Hess_impl(*Htmp, Qtmp.get(), NumTol);

  if constexpr (!FullyWritableMatrix<Matrix2>) {
    Q = *Qtmp;
  }
  if constexpr (!FullyWritableMatrix<Matrix3>) {
    H = *Htmp;
  }
}

/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \param A square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2>
void decompose_Hess(const Matrix1& A, Matrix2& H,
                    mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Upper-Hessenberg decomposition is only possible on a square matrix!");
  }

  auto Htmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix2>) {
      H = A;
      return std::unique_ptr<Matrix2, null_deleter>(&H, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix2>, mat_structure::square>>(A);
    }
  }();

  detail::decompose_Hess_impl(*Htmp, static_cast<decltype(Htmp.get())>(nullptr),
                              NumTol);

  if constexpr (!FullyWritableMatrix<Matrix2>) {
    H = *Htmp;
  }
}

/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param Q holds as output, the unitary square matrix Q.
 * \param Z holds as output, the unitary square matrix Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, WritableMatrix Matrix3,
          WritableMatrix Matrix4, WritableMatrix Matrix5, WritableMatrix Matrix6>
void reduce_HessTri(const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R,
                    Matrix5& Q, Matrix6& Z,
                    mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Hessenberg-Triangular reduction is only possible on square matrices!");
  }

  auto Htmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix3>) {
      H = A;
      return std::unique_ptr<Matrix3, null_deleter>(&H, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix3>, mat_structure::square>>(A);
    }
  }();
  auto Rtmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix4>) {
      R = B;
      return std::unique_ptr<Matrix4, null_deleter>(&R, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix4>, mat_structure::square>>(B);
    }
  }();
  auto Qtmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix5>) {
      Q = mat<mat_value_type_t<Matrix5>, mat_structure::identity>(
          A.get_row_count());
      return std::unique_ptr<Matrix5, null_deleter>(&Q, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix5>, mat_structure::square>>(
          mat<mat_value_type_t<Matrix5>, mat_structure::identity>(
              A.get_row_count()));
    }
  }();
  auto Ztmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix6>) {
      Z = mat<mat_value_type_t<Matrix6>, mat_structure::identity>(
          A.get_row_count());
      return std::unique_ptr<Matrix6, null_deleter>(&Z, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix6>, mat_structure::square>>(
          mat<mat_value_type_t<Matrix6>, mat_structure::identity>(
              A.get_row_count()));
    }
  }();

  detail::reduce_HessTri_offset_impl(*Htmp, *Rtmp, Qtmp.get(), Ztmp.get(), 0,
                                     NumTol);

  if constexpr (!FullyWritableMatrix<Matrix3>) {
    H = *Htmp;
  }
  if constexpr (!FullyWritableMatrix<Matrix4>) {
    R = *Rtmp;
  }
  if constexpr (!FullyWritableMatrix<Matrix5>) {
    Q = *Qtmp;
  }
  if constexpr (!FullyWritableMatrix<Matrix6>) {
    Z = *Ztmp;
  }
}

/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \tparam Matrix3 A fully-writable matrix type.
 * \tparam Matrix4 A fully-writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, WritableMatrix Matrix3,
          WritableMatrix Matrix4>
void reduce_HessTri(const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R,
                    mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Hessenberg-Triangular reduction is only possible on square matrices!");
  }

  auto Htmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix3>) {
      H = A;
      return std::unique_ptr<Matrix3, null_deleter>(&H, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix3>, mat_structure::square>>(A);
    }
  }();
  auto Rtmp = [&]() {
    if constexpr (FullyWritableMatrix<Matrix4>) {
      R = B;
      return std::unique_ptr<Matrix4, null_deleter>(&R, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix4>, mat_structure::square>>(B);
    }
  }();

  detail::reduce_HessTri_offset_impl(
      *Htmp, *Rtmp, static_cast<decltype(Htmp.get())>(nullptr),
      static_cast<decltype(Rtmp.get())>(nullptr), 0, NumTol);

  if constexpr (!FullyWritableMatrix<Matrix3>) {
    H = *Htmp;
  }
  if constexpr (!FullyWritableMatrix<Matrix4>) {
    R = *Rtmp;
  }
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_HESS_DECOMP_H_
