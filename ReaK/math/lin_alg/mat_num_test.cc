
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

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_balance.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/mat_ctrl_decomp.h"
#include "ReaK/math/lin_alg/mat_gaussian_elim.h"
#include "ReaK/math/lin_alg/mat_jacobi_method.h"
#include "ReaK/math/lin_alg/mat_matchers.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/lin_alg/mat_schur_decomp.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"

#include <cstdio>
#include <fstream>
#include <iostream>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

using ::ReaK::testing::MatrixIsDiagonal;
using ::ReaK::testing::MatrixIsInverseOf;
using ::ReaK::testing::MatrixIsLowerTriangular;
using ::ReaK::testing::MatrixIsNear;
using ::ReaK::testing::MatrixIsNull;
using ::ReaK::testing::MatrixIsOrthogonal;
using ::ReaK::testing::MatrixIsTriDiagonal;
using ::ReaK::testing::MatrixIsUpperHessenberg;
using ::ReaK::testing::MatrixIsUpperTriangular;

auto OriginalMatrix() {
  return mat<double, mat_structure::symmetric>(2.0, -1.0, 0.0, 2.0, -1.0, 2.0);
}
auto InverseMatrix() {
  return mat<double, mat_structure::symmetric>(0.75, 0.5, 0.25, 1.0, 0.5, 0.75);
}

const double tolerance = std::numeric_limits<double>::epsilon();
const double desired_tolerance = std::numeric_limits<double>::epsilon();

TEST(MatInversionTests, FullRankGaussian) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_gauss_inv(3);
  EXPECT_NO_THROW(invert_gaussian(m_orig, m_gauss_inv, desired_tolerance));
  EXPECT_THAT(m_gauss_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_gauss_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankPLU) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_plu_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(invert_PLU(m_orig, m_plu_inv, desired_tolerance));
  EXPECT_THAT(m_plu_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_plu_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankGaussianSymmetric) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::symmetric> m_gauss_inv(3);
  EXPECT_NO_THROW(invert_gaussian(m_orig, m_gauss_inv, desired_tolerance));
  EXPECT_THAT(m_gauss_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_gauss_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankPLUSymmetric) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::symmetric> m_plu_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(invert_PLU(m_orig, m_plu_inv, desired_tolerance));
  EXPECT_THAT(m_plu_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_plu_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankCholeskySymmetric) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::symmetric> m_cholesky_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(invert_Cholesky(m_orig, m_cholesky_inv, desired_tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsNear(InverseMatrix(), 4.0 * tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsInverseOf(m_orig, 4.0 * tolerance));
}

TEST(MatInversionTests, FullRankCholeskySquare) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_cholesky_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(invert_Cholesky(m_orig, m_cholesky_inv, desired_tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsNear(InverseMatrix(), 4.0 * tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsInverseOf(m_orig, 4.0 * tolerance));
}

TEST(MatInversionTests, FullRankCholeskyFromDiagonal) {
  mat<double, mat_structure::diagonal> m_orig_diag(vect<double, 3>(2, 1, 0.5));
  mat<double, mat_structure::diagonal> m_trueinv_diag(
      vect<double, 3>(0.5, 1, 2.0));

  mat<double, mat_structure::square> m_cholesky_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(
      invert_Cholesky(m_orig_diag, m_cholesky_inv, desired_tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsNear(m_trueinv_diag, tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsInverseOf(m_orig_diag, tolerance));
}

TEST(MatInversionTests, FullRankCholeskyFromDiagonalToDiagonal) {
  mat<double, mat_structure::diagonal> m_orig_diag(vect<double, 3>(2, 1, 0.5));
  mat<double, mat_structure::diagonal> m_trueinv_diag(
      vect<double, 3>(0.5, 1, 2.0));
  mat<double, mat_structure::diagonal> m_cholesky_inv(
      mat<double, mat_structure::identity>(3));
  EXPECT_NO_THROW(
      invert_Cholesky(m_orig_diag, m_cholesky_inv, desired_tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsNear(m_trueinv_diag, tolerance));
  EXPECT_THAT(m_cholesky_inv, MatrixIsInverseOf(m_orig_diag, tolerance));
}

TEST(MatInversionTests, FullRankJacobiPseudoInverse) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::symmetric> m_jacobi_pinv(3);
  EXPECT_NO_THROW(
      pseudoinvert_Jacobi(m_orig, m_jacobi_pinv, desired_tolerance));
  EXPECT_THAT(m_jacobi_pinv, MatrixIsNear(InverseMatrix(), 4.0 * tolerance));
  EXPECT_THAT(m_jacobi_pinv, MatrixIsInverseOf(m_orig, 4.0 * tolerance));
}

TEST(MatInversionTests, FullRankQRPseudoInverse) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_qr_inv(3);
  EXPECT_NO_THROW(pseudoinvert_QR(m_orig, m_qr_inv, desired_tolerance));
  EXPECT_THAT(m_qr_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_qr_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankSVDPseudoInverse) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_svd_inv(3);
  EXPECT_NO_THROW(pseudoinvert_SVD(m_orig, m_svd_inv, desired_tolerance));
  EXPECT_THAT(m_svd_inv, MatrixIsNear(InverseMatrix(), 4.0 * tolerance));
  EXPECT_THAT(m_svd_inv, MatrixIsInverseOf(m_orig, 4.0 * tolerance));
}

TEST(MatInversionTests, FullRankLDL) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_ldl_inv(3);
  EXPECT_NO_THROW(invert_LDL(m_orig, m_ldl_inv, desired_tolerance));
  EXPECT_THAT(m_ldl_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_ldl_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankBandCholesky) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_bandchol_inv(3);
  EXPECT_NO_THROW(
      invert_BandCholesky(m_orig, m_bandchol_inv, 1, desired_tolerance));
  EXPECT_THAT(m_bandchol_inv, MatrixIsNear(InverseMatrix(), 4.0 * tolerance));
  EXPECT_THAT(m_bandchol_inv, MatrixIsInverseOf(m_orig, 4.0 * tolerance));
}

TEST(MatInversionTests, FullRankSymmetricQRPseudoInverse) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_symqr_pinv(3);
  EXPECT_NO_THROW(pseudoinvert_SymQR(m_orig, m_symqr_pinv, desired_tolerance));
  EXPECT_THAT(m_symqr_pinv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_symqr_pinv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatInversionTests, FullRankTriDiagonalLDL) {
  const auto m_orig = OriginalMatrix();
  mat<double, mat_structure::square> m_tridiagldl_inv(3);
  EXPECT_NO_THROW(
      invert_TriDiagLDL(m_orig, m_tridiagldl_inv, desired_tolerance));
  EXPECT_THAT(m_tridiagldl_inv, MatrixIsNear(InverseMatrix(), tolerance));
  EXPECT_THAT(m_tridiagldl_inv, MatrixIsInverseOf(m_orig, tolerance));
}

TEST(MatNum, MatLeftPseudoInversionTests) {

  mat<double, mat_structure::rectangular> m_test(3, 2);
  m_test(0, 0) = 1.0;
  m_test(0, 1) = 3.0;
  m_test(1, 0) = 4.0;
  m_test(1, 1) = 5.0;
  m_test(2, 0) = 7.0;
  m_test(2, 1) = 8.0;

  mat<double, mat_structure::rectangular> m_qr_inv(2, 3);
  EXPECT_NO_THROW(pseudoinvert_QR(m_test, m_qr_inv, desired_tolerance));
  EXPECT_THAT(m_qr_inv, MatrixIsInverseOf(m_test, 4.0 * tolerance));

  mat<double, mat_structure::rectangular> m_svd_inv(2, 3);
  EXPECT_NO_THROW(pseudoinvert_SVD(m_test, m_svd_inv, desired_tolerance));
  EXPECT_THAT(m_svd_inv, MatrixIsInverseOf(m_test, 10.0 * tolerance));

  mat<double, mat_structure::rectangular> m_rrqr_inv(2, 3);
  try {
    linlsq_RRQR(m_test, m_rrqr_inv, mat<double, mat_structure::identity>(3),
                desired_tolerance);
  } catch (std::exception& e) {
    std::cout << "Exception thrown: " << e.what() << std::endl;
  }
  EXPECT_THAT(m_rrqr_inv, MatrixIsInverseOf(m_test, 4.0 * tolerance));

  // last column is equal to first column + 0.5 * second column (i.e., it is not full column-rank)
  mat<double, mat_structure::rectangular> m_rankdef_test(4, 3);
  m_rankdef_test(0, 0) = 1.0;
  m_rankdef_test(0, 1) = 3.0;
  m_rankdef_test(0, 2) = 2.5;
  m_rankdef_test(1, 0) = 4.0;
  m_rankdef_test(1, 1) = 5.0;
  m_rankdef_test(1, 2) = 6.5;
  m_rankdef_test(2, 0) = 7.0;
  m_rankdef_test(2, 1) = 8.0;
  m_rankdef_test(2, 2) = 11.0;
  m_rankdef_test(3, 0) = 2.0;
  m_rankdef_test(3, 1) = 4.0;
  m_rankdef_test(3, 2) = 4.0;

  mat<double, mat_structure::rectangular> m_svd_rankdef_inv(3, 4);
  EXPECT_NO_THROW(pseudoinvert_SVD(m_rankdef_test, m_svd_rankdef_inv, 1e-10));

  // NOTE: I don't know what else to do to check that the solution is correct, besides comparing to SVD solution.

  mat<double, mat_structure::rectangular> m_rrqr_rankdef_inv(3, 4);
  EXPECT_NO_THROW(linlsq_RRQR(m_rankdef_test, m_rrqr_rankdef_inv,
                              mat<double, mat_structure::identity>(4), 1e-10));

  EXPECT_THAT(m_svd_rankdef_inv, MatrixIsNear(m_rrqr_rankdef_inv, tolerance));
}

TEST(MatNum, MatRightPseudoInversionTests) {

  mat<double, mat_structure::rectangular> m_test(2, 3);
  m_test(0, 0) = 1.0;
  m_test(0, 1) = 2.0;
  m_test(0, 2) = 3.0;
  m_test(1, 0) = 4.0;
  m_test(1, 1) = 5.0;
  m_test(1, 2) = 6.0;

  mat<double, mat_structure::rectangular> m_test_R = m_test;
  mat<double, mat_structure::square> m_test_Q(
      mat<double, mat_structure::identity>(2));
  EXPECT_NO_THROW(
      (detail::decompose_QR_impl<mat<double, mat_structure::rectangular>,
                                 mat<double, mat_structure::square>>(
          m_test_R, &m_test_Q, desired_tolerance)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(tolerance));
  EXPECT_THAT(m_test_R, MatrixIsUpperTriangular(tolerance));
  EXPECT_THAT(m_test_Q * m_test_R, MatrixIsNear(m_test, 10.0 * tolerance));

  mat<double, mat_structure::rectangular> m_test_minnorm_x(3, 2);
  mat<double, mat_structure::rectangular> m_test_minnorm_b(2, 2);
  m_test_minnorm_b(0, 0) = 4.0;
  m_test_minnorm_b(0, 1) = 5.0;
  m_test_minnorm_b(1, 0) = 7.0;
  m_test_minnorm_b(1, 1) = 8.0;
  EXPECT_NO_THROW(minnorm_QR(m_test, m_test_minnorm_x, m_test_minnorm_b,
                             desired_tolerance));
  EXPECT_THAT(m_test * m_test_minnorm_x,
              MatrixIsNear(m_test_minnorm_b, 10.0 * tolerance));

  mat<double, mat_structure::rectangular> m_test_minnormRRQR_x(3, 2);
  EXPECT_NO_THROW(minnorm_RRQR(m_test, m_test_minnormRRQR_x, m_test_minnorm_b,
                               desired_tolerance));
  EXPECT_THAT(m_test * m_test_minnormRRQR_x,
              MatrixIsNear(m_test_minnorm_b, 20.0 * tolerance));

  /* SVD pseudo-invert for a min-norm problem. */
  mat<double, mat_structure::rectangular> m_test_svd_pinv(3, 2);
  EXPECT_NO_THROW(pseudoinvert_SVD(m_test, m_test_svd_pinv, desired_tolerance));
  mat<double, mat_structure::rectangular> m_test_svd_x =
      m_test_svd_pinv * m_test_minnorm_b;
  EXPECT_THAT(m_test * m_test_svd_x,
              MatrixIsNear(m_test_minnorm_b, 40.0 * tolerance));

  // last row is equal to first column + 0.5 * second column (i.e., it is not full column-rank)
  mat<double, mat_structure::rectangular> m_rankdef_test(3, 4);
  m_rankdef_test(0, 0) = 1.0;
  m_rankdef_test(0, 1) = 4.0;
  m_rankdef_test(0, 2) = 7.0;
  m_rankdef_test(0, 3) = 2.0;
  m_rankdef_test(1, 0) = 3.0;
  m_rankdef_test(1, 1) = 5.0;
  m_rankdef_test(1, 2) = 8.0;
  m_rankdef_test(1, 3) = 4.0;
  m_rankdef_test(2, 0) = 2.5;
  m_rankdef_test(2, 1) = 6.5;
  m_rankdef_test(2, 2) = 11.0;
  m_rankdef_test(2, 3) = 4.0;

  mat<double, mat_structure::rectangular> m_rankdef_minnorm_b(3, 2);
  m_rankdef_minnorm_b(0, 0) = 4.0;
  m_rankdef_minnorm_b(0, 1) = 5.0;
  m_rankdef_minnorm_b(1, 0) = 7.0;
  m_rankdef_minnorm_b(1, 1) = 8.0;
  m_rankdef_minnorm_b(2, 0) = 4.0;
  m_rankdef_minnorm_b(2, 1) = 2.0;

  mat<double, mat_structure::rectangular> m_rankdef_minnormRRQR_x(4, 2);
  EXPECT_NO_THROW(minnorm_RRQR(m_rankdef_test, m_rankdef_minnormRRQR_x,
                               m_rankdef_minnorm_b, desired_tolerance));

  // NOTE: I don't know what else to do to check that the solution is correct, besides comparing to SVD solution.

  mat<double, mat_structure::rectangular> m_rankdef_svd_pinv(4, 3);
  EXPECT_NO_THROW(
      pseudoinvert_SVD(m_rankdef_test, m_rankdef_svd_pinv, desired_tolerance));
  mat<double, mat_structure::rectangular> m_rankdef_svd_x =
      m_rankdef_svd_pinv * m_rankdef_minnorm_b;

  EXPECT_THAT(m_rankdef_svd_x,
              MatrixIsNear(m_rankdef_minnormRRQR_x, 40.0 * tolerance));
}

TEST(MatNum, MatEigenValuesTests) {

  mat<double, mat_structure::symmetric> m_gauss(2.0, -1.0, 0.0, 2.0, -1.0, 2.0);

  vect_n<double> m_gauss_trueeig(1.0 / (1.0 + std::sqrt(0.5)), 2.0,
                                 2.0 + std::sqrt(2.0));

  mat<double, mat_structure::diagonal> m_jacobi_E(3);
  mat<double, mat_structure::square> m_jacobi_Q(3);
  EXPECT_NO_THROW(
      eigensolve_Jacobi(m_gauss, m_jacobi_E, m_jacobi_Q, desired_tolerance));
  vect_n<double> m_jacobi_Es(m_jacobi_E(0, 0), m_jacobi_E(1, 1),
                             m_jacobi_E(2, 2));
  std::sort(m_jacobi_Es.begin(), m_jacobi_Es.end());
  EXPECT_NEAR(norm_inf(m_jacobi_Es - m_gauss_trueeig), 0.0, tolerance);
  EXPECT_THAT(m_jacobi_Q * m_jacobi_E * transpose(m_jacobi_Q),
              MatrixIsNear(m_gauss, 10.0 * tolerance));

  mat<double, mat_structure::diagonal> m_qr_E(3, true);
  mat<double, mat_structure::square> m_qr_Q(3);
  EXPECT_NO_THROW(eigensolve_SymQR(m_gauss, m_qr_Q, m_qr_E, desired_tolerance));
  vect_n<double> m_qr_Es(m_qr_E(0, 0), m_qr_E(1, 1), m_qr_E(2, 2));
  std::sort(m_qr_Es.begin(), m_qr_Es.end());
  EXPECT_NEAR(norm_inf(m_qr_Es - m_gauss_trueeig), 0.0, 10.0 * tolerance);
  EXPECT_THAT(m_qr_Q * m_qr_E * transpose(m_qr_Q),
              MatrixIsNear(m_gauss, 10.0 * tolerance));

  mat<double, mat_structure::diagonal> m_svd_E(3);
  mat<double, mat_structure::square> m_svd_U(3);
  mat<double, mat_structure::square> m_svd_V(3);
  EXPECT_NO_THROW(
      decompose_SVD(m_gauss, m_svd_U, m_svd_E, m_svd_V, desired_tolerance));
  vect_n<double> m_svd_Es(m_svd_E(0, 0), m_svd_E(1, 1), m_svd_E(2, 2));
  std::sort(m_svd_Es.begin(), m_svd_Es.end());
  EXPECT_NEAR(norm_inf(m_svd_Es - m_gauss_trueeig), 0.0, 10.0 * tolerance);
  EXPECT_THAT(m_svd_U * m_svd_E * transpose(m_svd_V),
              MatrixIsNear(m_gauss, 10.0 * tolerance));
}

TEST(MatNum, MatDecompositionsTests) {

  mat<double, mat_structure::rectangular> m_test(3, 3);
  m_test(0, 0) = 1.0;
  m_test(0, 1) = 3.0;
  m_test(0, 2) = 0.0;
  m_test(1, 0) = 3.0;
  m_test(1, 1) = 5.0;
  m_test(1, 2) = 2.0;
  m_test(2, 0) = 0.0;
  m_test(2, 1) = 2.0;
  m_test(2, 2) = 4.0;

  mat<double, mat_structure::rectangular> m_test_Gram = m_test;
  EXPECT_NO_THROW(
      (orthogonalize_StableGramSchmidt(m_test_Gram, true, desired_tolerance)));
  EXPECT_THAT(m_test_Gram, MatrixIsOrthogonal(2.0 * tolerance));

  mat<double, mat_structure::rectangular> m_test_R = m_test;
  mat<double, mat_structure::square> m_test_Q(
      mat<double, mat_structure::identity>(3));

  EXPECT_NO_THROW((detail::symmetric_QR_step(m_test_R, &m_test_Q)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(2.0 * tolerance));
  EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q),
              MatrixIsNear(m_test, 40.0 * tolerance));

  EXPECT_NO_THROW((detail::symmetric_QR_step(m_test_R, &m_test_Q)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(2.0 * tolerance));
  EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q),
              MatrixIsNear(m_test, 40.0 * tolerance));

  EXPECT_NO_THROW((detail::symmetric_QR_step(m_test_R, &m_test_Q)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(2.0 * tolerance));
  EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q),
              MatrixIsNear(m_test, 40.0 * tolerance));

  m_test(0, 0) = 1.0;
  m_test(0, 1) = 3.0;
  m_test(0, 2) = 0.0;
  m_test(1, 0) = 3.0;
  m_test(1, 1) = 5.0;
  m_test(1, 2) = 2.0;
  m_test(2, 0) = 0.0;
  m_test(2, 1) = 2.0;
  m_test(2, 2) = 4.0;
  m_test_R = m_test;
  m_test_Q = mat<double, mat_structure::identity>(3);
  EXPECT_NO_THROW(
      (detail::symmetric_QRalg_impl(m_test_R, &m_test_Q, desired_tolerance)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(2.0 * tolerance));
  EXPECT_THAT(m_test_R, MatrixIsUpperTriangular(2.0 * tolerance));
  EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q),
              MatrixIsNear(m_test, 20.0 * tolerance));

  m_test(0, 0) = 1.0;
  m_test(0, 1) = 3.0;
  m_test(0, 2) = 2.0;
  m_test(1, 0) = 3.0;
  m_test(1, 1) = 5.0;
  m_test(1, 2) = 2.0;
  m_test(2, 0) = 2.0;
  m_test(2, 1) = 2.0;
  m_test(2, 2) = 4.0;
  m_test_R = m_test;
  m_test_Q = mat<double, mat_structure::identity>(3);
  EXPECT_NO_THROW(
      (detail::decompose_TriDiag_impl(m_test_R, &m_test_Q, desired_tolerance)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(2.0 * tolerance));
  EXPECT_THAT(m_test_R, MatrixIsTriDiagonal(2.0 * tolerance));
  EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q),
              MatrixIsNear(m_test, 20.0 * tolerance));

  m_test(0, 0) = 1.0;
  m_test(0, 1) = 3.0;
  m_test(0, 2) = 2.0;
  m_test(1, 0) = 3.0;
  m_test(1, 1) = 5.0;
  m_test(1, 2) = 2.0;
  m_test(2, 0) = 2.0;
  m_test(2, 1) = 2.0;
  m_test(2, 2) = 4.0;
  m_test_R = m_test;
  m_test_Q = mat<double, mat_structure::identity>(3);
  EXPECT_NO_THROW(
      (detail::symmetric_QRalg_impl(m_test_R, &m_test_Q, desired_tolerance)));
  EXPECT_THAT(m_test_Q, MatrixIsOrthogonal(4.0 * tolerance));
  EXPECT_THAT(m_test_R, MatrixIsUpperTriangular(4.0 * tolerance));
  // TODO This test is failing with too much error (~1e-11).
  // EXPECT_THAT(m_test_Q * m_test_R * transpose_view(m_test_Q), MatrixIsNear(m_test, 40.0 * tolerance));

  mat<double, mat_structure::rectangular> m_test4(4, 4);
  mat<double, mat_structure::square> m_test4_Q(
      mat<double, mat_structure::identity>(4));
  m_test4(0, 0) = 1.0;
  m_test4(0, 1) = 3.0;
  m_test4(0, 2) = 2.0;
  m_test4(0, 3) = -1.0;
  m_test4(1, 0) = 3.0;
  m_test4(1, 1) = 5.0;
  m_test4(1, 2) = 2.0;
  m_test4(1, 3) = 5.0;
  m_test4(2, 0) = 2.0;
  m_test4(2, 1) = 2.0;
  m_test4(2, 2) = 4.0;
  m_test4(2, 3) = 3.0;
  m_test4(3, 0) = -1.0;
  m_test4(3, 1) = 5.0;
  m_test4(3, 2) = 3.0;
  m_test4(3, 3) = -2.0;
  mat<double, mat_structure::rectangular> m_test4_R = m_test4;
  m_test4_Q = mat<double, mat_structure::identity>(4);
  EXPECT_NO_THROW(
      (detail::symmetric_QRalg_impl(m_test4_R, &m_test4_Q, desired_tolerance)));
  EXPECT_THAT(m_test4_Q, MatrixIsOrthogonal(4.0 * tolerance));
  EXPECT_THAT(m_test4_R, MatrixIsUpperTriangular(4.0 * tolerance));
  // TODO This test is failing with too much error (~1e-10).
  // EXPECT_THAT(m_test4_Q * m_test4_R * transpose_view(m_test4_Q), MatrixIsNear(m_test4, 40.0 * tolerance));

  mat<double, mat_structure::square> m_hess_Q(
      mat<double, mat_structure::identity>(4));
  mat<double, mat_structure::square> m_hess_H(m_test4);
  EXPECT_NO_THROW(
      (decompose_Hess(m_test4, m_hess_Q, m_hess_H, desired_tolerance)));
  EXPECT_THAT(m_hess_Q, MatrixIsOrthogonal(4.0 * tolerance));
  EXPECT_THAT(m_hess_H, MatrixIsUpperHessenberg(4.0 * tolerance));
  // A = Q T Q^T
  EXPECT_THAT(m_hess_Q * m_hess_H * transpose_view(m_hess_Q),
              MatrixIsNear(m_test4, 40.0 * tolerance));

  mat<double, mat_structure::square> m_realshur_Q(
      mat<double, mat_structure::identity>(4));
  mat<double, mat_structure::square> m_realshur_T(m_test4);
  EXPECT_NO_THROW((decompose_RealSchur(m_test4, m_realshur_Q, m_realshur_T,
                                       desired_tolerance)));
  EXPECT_THAT(m_realshur_Q, MatrixIsOrthogonal(8.0 * tolerance));
  EXPECT_THAT(m_realshur_T, MatrixIsUpperHessenberg(8.0 * tolerance));
  // A = Q T Q^T
  EXPECT_THAT(m_realshur_Q * m_realshur_T * transpose_view(m_realshur_Q),
              MatrixIsNear(m_test4, 80.0 * tolerance));

  mat<double, mat_structure::square> m_test_L(3, 0.0);
  mat<double, mat_structure::square> m_test_D(3, 0.0);
  mat<double, mat_structure::square> m_test_sqr(3, 0.0);
  mat<double, mat_structure::square> m_test_inv(3, 0.0);
  m_test_sqr(0, 0) = 6.0;
  m_test_sqr(0, 1) = 3.0;
  m_test_sqr(0, 2) = 2.0;
  m_test_sqr(1, 0) = 3.0;
  m_test_sqr(1, 1) = 5.0;
  m_test_sqr(1, 2) = 2.0;
  m_test_sqr(2, 0) = 2.0;
  m_test_sqr(2, 1) = 2.0;
  m_test_sqr(2, 2) = 4.0;
  EXPECT_NO_THROW(
      (decompose_LDL(m_test_sqr, m_test_L, m_test_D, desired_tolerance)));
  EXPECT_THAT(m_test_L, MatrixIsLowerTriangular(2.0 * tolerance));
  EXPECT_THAT(m_test_D, MatrixIsDiagonal(2.0 * tolerance));
  EXPECT_THAT(m_test_L * m_test_D * transpose_view(m_test_L),
              MatrixIsNear(m_test_sqr, 10.0 * tolerance));

  m_test_sqr(0, 0) = 6.0;
  m_test_sqr(0, 1) = 3.0;
  m_test_sqr(0, 2) = 0.0;
  m_test_sqr(1, 0) = 3.0;
  m_test_sqr(1, 1) = 5.0;
  m_test_sqr(1, 2) = 2.0;
  m_test_sqr(2, 0) = 0.0;
  m_test_sqr(2, 1) = 2.0;
  m_test_sqr(2, 2) = 4.0;
  EXPECT_NO_THROW(
      (decompose_BandCholesky(m_test_sqr, m_test_L, 1, desired_tolerance)));
  EXPECT_THAT(m_test_L, MatrixIsLowerTriangular(2.0 * tolerance));
  EXPECT_THAT(m_test_L * transpose_view(m_test_L),
              MatrixIsNear(m_test_sqr, 10.0 * tolerance));

  m_test_D = mat<double, mat_structure::identity>(3);
  m_test_sqr(0, 0) = 6.0;
  m_test_sqr(0, 1) = 3.0;
  m_test_sqr(0, 2) = 0.0;
  m_test_sqr(1, 0) = 3.0;
  m_test_sqr(1, 1) = 5.0;
  m_test_sqr(1, 2) = 2.0;
  m_test_sqr(2, 0) = 0.0;
  m_test_sqr(2, 1) = 2.0;
  m_test_sqr(2, 2) = 4.0;
  EXPECT_NO_THROW((
      decompose_TriDiagLDL(m_test_sqr, m_test_L, m_test_D, desired_tolerance)));
  EXPECT_THAT(m_test_L, MatrixIsLowerTriangular(2.0 * tolerance));
  EXPECT_THAT(m_test_D, MatrixIsDiagonal(2.0 * tolerance));
  EXPECT_THAT(m_test_L * m_test_D * transpose_view(m_test_L),
              MatrixIsNear(m_test_sqr, 10.0 * tolerance));
}

TEST(MatNum, MatCtrlReductionTests) {

  std::size_t N = 6;
  std::size_t M = 2;

  mat<double, mat_structure::rectangular> A(N, N);
  mat<double, mat_structure::rectangular> B(N, M);
  A(0, 0) = 0.0;
  A(0, 1) = 0.0;
  A(0, 2) = 0.0;
  A(0, 3) = 1.0;
  A(0, 4) = 0.0;
  A(0, 5) = 0.0;
  A(1, 0) = 0.0;
  A(1, 1) = 0.0;
  A(1, 2) = 0.0;
  A(1, 3) = 0.0;
  A(1, 4) = 1.0;
  A(1, 5) = 0.0;
  A(2, 0) = 0.0;
  A(2, 1) = 0.0;
  A(2, 2) = 0.0;
  A(2, 3) = 0.0;
  A(2, 4) = 0.0;
  A(2, 5) = 1.0;
  A(3, 0) = 0.0;
  A(3, 1) = 0.0;
  A(3, 2) = 0.0;
  A(3, 3) = 0.0;
  A(3, 4) = 0.0;
  A(3, 5) = 0.0;
  A(4, 0) = 0.0;
  A(4, 1) = 0.0;
  A(4, 2) = 0.0;
  A(4, 3) = 0.0;
  A(4, 4) = 0.0;
  A(4, 5) = 0.0;
  A(5, 0) = 0.0;
  A(5, 1) = 0.0;
  A(5, 2) = 0.0;
  A(5, 3) = 0.0;
  A(5, 4) = 0.0;
  A(5, 5) = 0.0;

  B(0, 0) = 0.0;
  B(0, 1) = 0.0;
  B(1, 0) = 0.0;
  B(1, 1) = 0.0;
  B(2, 0) = 0.0;
  B(2, 1) = 0.0;
  B(3, 0) = 0.0;
  B(3, 1) = 0.0;
  B(4, 0) = 1.0;
  B(4, 1) = 0.0;
  B(5, 0) = 0.0;
  B(5, 1) = 1.0;

  mat<double, mat_structure::rectangular> Ar = A;
  mat<double, mat_structure::rectangular> Br = B;

  mat<double, mat_structure::square> Q(
      (mat<double, mat_structure::identity>(N)));
  mat<double, mat_structure::square> Z(
      (mat<double, mat_structure::identity>(M)));
  std::size_t num_rank = N;
  EXPECT_NO_THROW(num_rank = ctrl_reduction(Ar, Br, Q, Z, desired_tolerance));

  // check basic structures, Q and Z being orthogonal, and equality is preserved (A = Q Ar Qt and B = Q Br Zt).
  EXPECT_THAT(Q, MatrixIsOrthogonal(4.0 * tolerance));
  EXPECT_THAT(Z, MatrixIsOrthogonal(4.0 * tolerance));
  EXPECT_THAT(Q * Ar * transpose(Q), MatrixIsNear(A, 10.0 * tolerance));
  EXPECT_THAT(Q * Br * transpose(Z), MatrixIsNear(B, 10.0 * tolerance));

  // Check the reduced controllability matrix:
  mat<double, mat_structure::rectangular> Cr(N, M * N);
  sub(Cr)(range(0, N), range(0, 2)) = Br;
  sub(Cr)(range(0, N), range(2, 4)) = Ar * Br;
  sub(Cr)(range(0, N), range(4, 6)) = Ar * Ar * Br;
  sub(Cr)(range(0, N), range(6, 8)) = Ar * Ar * Ar * Br;
  sub(Cr)(range(0, N), range(8, 10)) = Ar * Ar * Ar * Ar * Br;
  sub(Cr)(range(0, N), range(10, 12)) = Ar * Ar * Ar * Ar * Ar * Br;
  EXPECT_EQ(num_rank, 4);
  EXPECT_THAT(sub(Cr)(range(num_rank, 6), range(0, 12)),
              MatrixIsNull(10.0 * tolerance));

  /* Use row-rank SVD to check that the controllable parts are full row-rank. */
  mat<double, mat_structure::diagonal> Er(num_rank);
  mat<double, mat_structure::rectangular> Ur(num_rank, num_rank);
  mat<double, mat_structure::rectangular> Vr(M * N, num_rank);
  EXPECT_NO_THROW(decompose_SVD(sub(Cr)(range(0, num_rank), range(0, 12)), Ur,
                                Er, Vr, desired_tolerance));
  EXPECT_EQ((numrank_SVD(Er, desired_tolerance)), num_rank);
}

TEST(MatNum, MatBalancingTests) {

  mat<double, mat_structure::square> A(3);
  A(0, 0) = 1.0;
  A(0, 1) = 0.0;
  A(0, 2) = 1e-4;
  A(1, 0) = 1.0;
  A(1, 1) = 1.0;
  A(1, 2) = 1e-2;
  A(2, 0) = 1e4;
  A(2, 1) = 1e2;
  A(2, 2) = 1.0;

  mat<double, mat_structure::square> A_bal = A;
  vect_n<int> D_bal(3);
  EXPECT_NO_THROW(balance(A_bal, D_bal));

  mat<double, mat_structure::diagonal> A_E(3);
  mat<double, mat_structure::square> A_U(3);
  mat<double, mat_structure::square> A_V(3);
  EXPECT_NO_THROW(decompose_SVD(A, A_U, A_E, A_V, desired_tolerance));
  double cond_before = condition_number_SVD(A_E);

  mat<double, mat_structure::diagonal> A_bal_E(3);
  mat<double, mat_structure::square> A_bal_U(3);
  mat<double, mat_structure::square> A_bal_V(3);
  EXPECT_NO_THROW(
      decompose_SVD(A_bal, A_bal_U, A_bal_E, A_bal_V, desired_tolerance));
  double cond_after = condition_number_SVD(A_bal_E);

  // The condition number must have been reduced:
  EXPECT_LT(cond_after, cond_before);

  // Condition number must be with the N-th power of 2 (at most the difference between
  // any two eigen-value is a factor of 2, so, at most, the difference between the largest
  // and smallest eigen-value is 2^N. Here, N = 3.
  EXPECT_LT(cond_after, std::ldexp(1, 3));

  apply_left_bal_exp(D_bal, A_bal);
  apply_right_bal_inv_exp(A_bal, D_bal);

  // must be exact because exact floating-point arithmetic is used in balancing:
  EXPECT_THAT(A, MatrixIsNear(A_bal, 1e-200));

  mat<double, mat_structure::square> B(3);
  B(0, 0) = 1.0;
  B(0, 1) = 1e-2;
  B(0, 2) = 1e-4;
  B(1, 0) = 1e3;
  B(1, 1) = 1.0;
  B(1, 2) = 1.0;
  B(2, 0) = 1e4;
  B(2, 1) = 0.0;
  B(2, 2) = 2.0;

  mat<double, mat_structure::square> A_penbal = A;
  mat<double, mat_structure::square> B_penbal = B;
  vect_n<int> Dl_penbal(3);
  vect_n<int> Dr_penbal(3);
  EXPECT_NO_THROW(balance_pencil(A_penbal, B_penbal, Dl_penbal, Dr_penbal));

  mat<double, mat_structure::diagonal> B_E(3);
  mat<double, mat_structure::square> B_U(3);
  mat<double, mat_structure::square> B_V(3);
  EXPECT_NO_THROW(decompose_SVD(B, B_U, B_E, B_V, desired_tolerance));
  double cond_B_before = condition_number_SVD(B_E);

  mat<double, mat_structure::square> M_pen(3);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      M_pen(i, j) = A(i, j) * A(i, j) + B(i, j) * B(i, j);
    }
  }
  mat<double, mat_structure::diagonal> M_pen_E(3);
  mat<double, mat_structure::square> M_pen_U(3);
  mat<double, mat_structure::square> M_pen_V(3);
  EXPECT_NO_THROW(
      decompose_SVD(M_pen, M_pen_U, M_pen_E, M_pen_V, desired_tolerance));
  double cond_M_before = condition_number_SVD(M_pen_E);

  mat<double, mat_structure::diagonal> A_penbal_E(3);
  mat<double, mat_structure::square> A_penbal_U(3);
  mat<double, mat_structure::square> A_penbal_V(3);
  EXPECT_NO_THROW(decompose_SVD(A_penbal, A_penbal_U, A_penbal_E, A_penbal_V,
                                desired_tolerance));
  double cond_A_penbal = condition_number_SVD(A_penbal_E);

  mat<double, mat_structure::diagonal> B_penbal_E(3);
  mat<double, mat_structure::square> B_penbal_U(3);
  mat<double, mat_structure::square> B_penbal_V(3);
  EXPECT_NO_THROW(decompose_SVD(B_penbal, B_penbal_U, B_penbal_E, B_penbal_V,
                                desired_tolerance));
  double cond_B_penbal = condition_number_SVD(B_penbal_E);

  mat<double, mat_structure::square> M_penbal(3);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      M_penbal(i, j) =
          A_penbal(i, j) * A_penbal(i, j) + B_penbal(i, j) * B_penbal(i, j);
    }
  }
  mat<double, mat_structure::diagonal> M_penbal_E(3);
  mat<double, mat_structure::square> M_penbal_U(3);
  mat<double, mat_structure::square> M_penbal_V(3);
  EXPECT_NO_THROW(decompose_SVD(M_penbal, M_penbal_U, M_penbal_E, M_penbal_V,
                                desired_tolerance));
  double cond_M_after = condition_number_SVD(M_penbal_E);

  apply_left_bal_exp(Dl_penbal, A_penbal);
  apply_left_bal_exp(Dl_penbal, B_penbal);
  apply_right_bal_inv_exp(A_penbal, Dr_penbal);
  apply_right_bal_inv_exp(B_penbal, Dr_penbal);

  // The condition numbers must have been reduced:
  EXPECT_LT(cond_A_penbal, cond_before);
  EXPECT_LT(cond_B_penbal, cond_B_before);
  EXPECT_LT(cond_M_after, cond_M_before);

  EXPECT_LT(cond_A_penbal, 150.0);
  EXPECT_LT(cond_B_penbal, 10.0);
  EXPECT_LT(cond_M_after, 600.0);

  // must be exact because exact floating-point arithmetic is used in balancing:
  EXPECT_THAT(A, MatrixIsNear(A_penbal, 1e-200));
  EXPECT_THAT(B, MatrixIsNear(B_penbal, 1e-200));
}

}  // namespace
}  // namespace ReaK
