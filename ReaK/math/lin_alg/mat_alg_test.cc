
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
#include "ReaK/core/base/defs.h"

#include <cstdio>
#include <fstream>
#include <iostream>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

template <typename M>
M CreateM1234() {
  return M(1.0, 2.0, 3.0, 4.0);
}
template <typename M>
M CreateM1234Inv() {
  return M(-2.0, 1.0, 1.5, -0.5);
}
template <typename M>
M CreateM123Sym() {
  return M(1.0, 2.0, 3.0);
}
template <typename M>
M CreateM123SymInv() {
  return M(-3.0, 2.0, -1.0);
}

template <typename Scalar>
class SpecialMatrixTest : public ::testing::Test {
 protected:
  static constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
};

using SpecialMatrixTestTypes = ::testing::Types<double, float>;

TYPED_TEST_SUITE(SpecialMatrixTest, SpecialMatrixTestTypes);

TYPED_TEST(SpecialMatrixTest, IdentityAndNull) {
  using std::abs;
  using Scalar = TypeParam;
  const Scalar tolerance = this->tolerance;

  mat<Scalar, mat_structure::rectangular> m_ident(2, 2, true);
  EXPECT_TRUE(abs(m_ident(0, 0) - 1.0) < tolerance);
  EXPECT_TRUE(abs(m_ident(0, 1)) < tolerance);
  EXPECT_TRUE(abs(m_ident(1, 0)) < tolerance);
  EXPECT_TRUE(abs(m_ident(1, 1) - 1.0) < tolerance);

  EXPECT_TRUE(is_identity_mat(m_ident, tolerance));

  mat_identity_t<Scalar> m_ident2(2);
  EXPECT_TRUE(is_identity_mat(m_ident2, tolerance));

  mat_null_t<Scalar> m_zeroes(2, 2);
  EXPECT_TRUE(is_null_mat(m_zeroes, tolerance));
  EXPECT_TRUE(is_null_mat(m_zeroes * m_ident2, tolerance));
  EXPECT_TRUE(is_null_mat(m_ident2 * m_zeroes, tolerance));
}

TYPED_TEST(SpecialMatrixTest, IdentityVsNull) {
  using Scalar = TypeParam;
  const Scalar tolerance = this->tolerance;
  mat_identity_t<Scalar> m_ident2(2);
  mat_null_t<Scalar> m_zeroes(2, 2);
  EXPECT_TRUE(is_identity_mat((m_zeroes * m_ident2) + m_ident2, tolerance));
  EXPECT_TRUE(is_identity_mat((m_ident2 * m_zeroes) + m_ident2, tolerance));
  EXPECT_TRUE(is_identity_mat((m_zeroes + m_zeroes) + m_ident2, tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 - (m_zeroes + m_zeroes), tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 + (m_zeroes - m_zeroes), tolerance));
  EXPECT_TRUE(is_identity_mat(-((m_zeroes - m_zeroes) - m_ident2), tolerance));
}

TYPED_TEST(SpecialMatrixTest, IdentityVsIdentity) {
  using Scalar = TypeParam;
  const Scalar tolerance = this->tolerance;
  mat_identity_t<Scalar> m_ident2(2);
  EXPECT_TRUE(is_identity_mat((m_ident2 + m_ident2) * 0.5, tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m_ident2 + m_ident2), tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 * m_ident2, tolerance));
  EXPECT_TRUE(is_identity_mat((m_ident2 - m_ident2) + m_ident2, tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 + (m_ident2 - m_ident2), tolerance));
}

template <typename MatDense>
class DenseTest : public ::testing::Test {
 protected:
  using Scalar = mat_value_type_t<MatDense>;
  static constexpr mat_alignment::tag Alignment =
      mat_traits<MatDense>::alignment;
  static constexpr unsigned int RowCount =
      mat_traits<MatDense>::static_row_count;
  static constexpr unsigned int ColCount =
      mat_traits<MatDense>::static_col_count;
  DenseTest()
      : m1234(CreateM1234<MatDense>()), m1234_inv(CreateM1234Inv<MatDense>()) {}

  MatDense m1234;
  MatDense m1234_inv;
  static constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
};

using DenseTestTypes = ::testing::Types<
    mat<double, mat_structure::rectangular, mat_alignment::column_major>,
    mat<float, mat_structure::rectangular, mat_alignment::column_major>,
    mat<double, mat_structure::rectangular, mat_alignment::row_major>,
    mat<float, mat_structure::rectangular, mat_alignment::row_major>,
    mat<double, mat_structure::rectangular, mat_alignment::column_major, 2, 2>,
    mat<float, mat_structure::rectangular, mat_alignment::column_major, 2, 2>,
    mat<double, mat_structure::rectangular, mat_alignment::row_major, 2, 2>,
    mat<float, mat_structure::rectangular, mat_alignment::row_major, 2, 2>,
    mat<double, mat_structure::square, mat_alignment::column_major>,
    mat<float, mat_structure::square, mat_alignment::column_major>,
    mat<double, mat_structure::square, mat_alignment::row_major>,
    mat<float, mat_structure::square, mat_alignment::row_major>,
    mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>,
    mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>,
    mat<double, mat_structure::square, mat_alignment::row_major, 2, 2>,
    mat<float, mat_structure::square, mat_alignment::row_major, 2, 2>>;

TYPED_TEST_SUITE(DenseTest, DenseTestTypes);

TYPED_TEST(DenseTest, Transpose) {
  const auto tolerance = this->tolerance;
  auto& m1234 = this->m1234;
  using std::abs;
  EXPECT_TRUE(abs(m1234(0, 0) - 1.0) < tolerance);
  EXPECT_TRUE(abs(m1234(0, 1) - 2.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 0) - 3.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 1) - 4.0) < tolerance);

  m1234 = transpose(m1234);
  EXPECT_TRUE(abs(m1234(0, 0) - 1.0) < tolerance);
  EXPECT_TRUE(abs(m1234(0, 1) - 3.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 0) - 2.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 1) - 4.0) < tolerance);

  m1234 = transpose_move(m1234);
  EXPECT_TRUE(abs(m1234(0, 0) - 1.0) < tolerance);
  EXPECT_TRUE(abs(m1234(0, 1) - 2.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 0) - 3.0) < tolerance);
  EXPECT_TRUE(abs(m1234(1, 1) - 4.0) < tolerance);
}

TYPED_TEST(DenseTest, Operators) {
  const auto tolerance = this->tolerance;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  mat_identity_t<typename TestFixture::Scalar> m_ident2(2);
  EXPECT_TRUE(is_identity_mat(m1234_inv * m1234, tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m1234_inv * m1234 + m1234_inv * m1234),
                              tolerance));
  EXPECT_TRUE(is_identity_mat((m1234_inv * m1234 + m1234_inv * m1234) * 0.5,
                              tolerance));
  EXPECT_TRUE(is_identity_mat(
      (m1234_inv * m1234 - m1234_inv * m1234) + m_ident2, tolerance));
  EXPECT_TRUE(is_null_mat(m1234_inv * m1234 - m1234_inv * m1234, tolerance));
  EXPECT_TRUE(
      is_null_mat(m1234_inv * m1234 + (-(m1234_inv * m1234)), tolerance));
}

TYPED_TEST(DenseTest, VsDiagonal) {
  const auto tolerance = this->tolerance;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  using Scalar = typename TestFixture::Scalar;
  mat<Scalar, mat_structure::diagonal, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_ident_diag(2, Scalar{1.0});
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m_ident_diag * m1234), tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m1234 * m_ident_diag), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_diag), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m_ident_diag + m1234 * m1234_inv), tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident_diag + (m1234 * m1234_inv - m_ident_diag),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident_diag + (m_ident_diag - m1234 * m1234_inv),
                              tolerance));
}

TYPED_TEST(DenseTest, VsScalarMat) {
  const auto tolerance = this->tolerance;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  using Scalar = typename TestFixture::Scalar;
  mat<Scalar, mat_structure::scalar, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_ident_scalar(2, Scalar{1.0});
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m_ident_scalar * m1234), tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m1234 * m_ident_scalar), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_scalar), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m_ident_scalar + m1234 * m1234_inv), tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_scalar + (m1234 * m1234_inv - m_ident_scalar), tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_scalar + (m_ident_scalar - m1234 * m1234_inv), tolerance));
}

TYPED_TEST(DenseTest, VsIdentity) {
  const auto tolerance = this->tolerance;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  mat<typename TestFixture::Scalar, mat_structure::identity,
      TestFixture::Alignment, TestFixture::RowCount, TestFixture::RowCount>
      m_ident2(2);
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m_ident2 * m1234), tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m1234 * m_ident2), tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident2), tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m_ident2 + m1234 * m1234_inv), tolerance));
  EXPECT_TRUE(
      is_identity_mat(m_ident2 + (m1234 * m1234_inv - m_ident2), tolerance));
  EXPECT_TRUE(
      is_identity_mat(m_ident2 + (m_ident2 - m1234 * m1234_inv), tolerance));
}

TYPED_TEST(DenseTest, VsNull) {
  const auto tolerance = this->tolerance;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  mat<typename TestFixture::Scalar, mat_structure::nil, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_zeroes(2, 2);
  EXPECT_TRUE(is_null_mat(m1234 * m_zeroes, tolerance));
  EXPECT_TRUE(is_null_mat(m_zeroes * m1234, tolerance));
  EXPECT_TRUE(is_identity_mat(m_zeroes - m1234_inv * (-m1234), tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * m1234 - m_zeroes, tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * m1234 + m_zeroes, tolerance));
  EXPECT_TRUE(is_identity_mat(m_zeroes + m1234_inv * m1234, tolerance));
}

template <typename MatDenseAndDense>
class DenseDenseTest : public ::testing::Test {
 protected:
  using MatDense1 = std::tuple_element_t<0, MatDenseAndDense>;
  using MatDense2 = std::tuple_element_t<1, MatDenseAndDense>;
  using Scalar = mat_value_type_t<MatDense1>;
  DenseDenseTest()
      : m1234_1(CreateM1234<MatDense1>()),
        m1234_inv_1(CreateM1234Inv<MatDense1>()),
        m1234_2(CreateM1234<MatDense2>()),
        m1234_inv_2(CreateM1234Inv<MatDense2>()) {}

  MatDense1 m1234_1;
  MatDense1 m1234_inv_1;
  MatDense2 m1234_2;
  MatDense2 m1234_inv_2;
  static constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
};

using DenseDenseTestTypes = ::testing::Types<
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::column_major>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::column_major>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::row_major>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::row_major>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::column_major>,
        mat<double, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::column_major>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<mat<double, mat_structure::square, mat_alignment::column_major>,
               mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<float, mat_structure::square, mat_alignment::column_major>,
               mat<float, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<double, mat_structure::square, mat_alignment::row_major>,
               mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<float, mat_structure::square, mat_alignment::row_major>,
               mat<float, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<double, mat_structure::square, mat_alignment::column_major>,
               mat<double, mat_structure::rectangular, mat_alignment::row_major,
                   2, 2>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::column_major>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::row_major, 2, 2>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::row_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>,
        mat<double, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::column_major>,
        mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::column_major>,
        mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::row_major>,
        mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<float, mat_structure::rectangular, mat_alignment::row_major>,
               mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::column_major>,
        mat<double, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::column_major>,
        mat<float, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::column_major>,
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::column_major>,
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::row_major>,
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::row_major>,
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::column_major>,
        mat<double, mat_structure::square, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::column_major>,
        mat<float, mat_structure::square, mat_alignment::row_major, 2, 2>>,
    std::tuple<mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<float, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>,
               mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::rectangular, mat_alignment::row_major, 2, 2>,
        mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::rectangular, mat_alignment::row_major, 2, 2>,
        mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<mat<float, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>,
               mat<float, mat_structure::square, mat_alignment::row_major>>>;

TYPED_TEST_SUITE(DenseDenseTest, DenseDenseTestTypes);

TYPED_TEST(DenseDenseTest, Operators) {
  const auto tolerance = this->tolerance;
  const auto& m1234_1 = this->m1234_1;
  const auto& m1234_inv_1 = this->m1234_inv_1;
  const auto& m1234_2 = this->m1234_2;
  const auto& m1234_inv_2 = this->m1234_inv_2;
  EXPECT_TRUE(is_identity_mat(m1234_inv_1 * ((m1234_inv_2 * m1234_2) * m1234_1),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv_1 * (m1234_1 * (m1234_inv_2 * m1234_2)),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(
      0.5 * (m1234_1 * m1234_inv_1 + (m1234_inv_2 * m1234_2)), tolerance));
  EXPECT_TRUE(is_identity_mat(
      0.5 * ((m1234_inv_2 * m1234_2) + m1234_1 * m1234_inv_1), tolerance));
  EXPECT_TRUE(
      is_identity_mat((m1234_inv_2 * m1234_2) +
                          (m1234_1 * m1234_inv_1 - (m1234_inv_2 * m1234_2)),
                      tolerance));
  EXPECT_TRUE(
      is_identity_mat((m1234_inv_2 * m1234_2) +
                          ((m1234_inv_2 * m1234_2) - m1234_1 * m1234_inv_1),
                      tolerance));
}

template <typename MatSym>
class SymmetricTest : public ::testing::Test {
 protected:
  using Scalar = mat_value_type_t<MatSym>;
  static constexpr mat_alignment::tag Alignment = mat_traits<MatSym>::alignment;
  static constexpr unsigned int RowCount = mat_traits<MatSym>::static_row_count;
  SymmetricTest()
      : m123_sym(CreateM123Sym<MatSym>()),
        m123_inv_sym(CreateM123SymInv<MatSym>()) {}

  MatSym m123_sym;
  MatSym m123_inv_sym;
  static constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
};

using SymmetricTestTypes = ::testing::Types<
    mat<double, mat_structure::symmetric>, mat<float, mat_structure::symmetric>,
    mat<double, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
    mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>>;

TYPED_TEST_SUITE(SymmetricTest, SymmetricTestTypes);

TYPED_TEST(SymmetricTest, Operators) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  mat_identity_t<typename TestFixture::Scalar> m_ident2(2);
  EXPECT_TRUE(is_identity_mat(m123_inv_sym * m123_sym, tolerance));
  EXPECT_TRUE(is_identity_mat(
      0.5 * (m123_inv_sym * m123_sym + m123_inv_sym * m123_sym), tolerance));
  EXPECT_TRUE(is_identity_mat(
      (m123_inv_sym * m123_sym + m123_inv_sym * m123_sym) * 0.5, tolerance));
  EXPECT_TRUE(is_identity_mat(
      (m123_inv_sym * m123_sym - m123_inv_sym * m123_sym) + m_ident2,
      tolerance));
  EXPECT_TRUE(is_null_mat(m123_inv_sym * m123_sym - m123_inv_sym * m123_sym,
                          tolerance));
  EXPECT_TRUE(is_null_mat(
      m123_inv_sym * m123_sym + (-(m123_inv_sym * m123_sym)), tolerance));
}

TYPED_TEST(SymmetricTest, VsDiagonal) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  using Scalar = typename TestFixture::Scalar;
  mat<Scalar, mat_structure::diagonal, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_ident_diag(2, Scalar{1.0});
  EXPECT_TRUE(
      is_identity_mat(m123_inv_sym * (m_ident_diag * m123_sym), tolerance));
  EXPECT_TRUE(
      is_identity_mat(m123_inv_sym * (m123_sym * m_ident_diag), tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m123_sym * m123_inv_sym + m_ident_diag),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m_ident_diag + m123_sym * m123_inv_sym),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_diag + (m123_sym * m123_inv_sym - m_ident_diag), tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_diag + (m_ident_diag - m123_sym * m123_inv_sym), tolerance));
}

TYPED_TEST(SymmetricTest, VsScalarMat) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  using Scalar = typename TestFixture::Scalar;
  mat<Scalar, mat_structure::scalar, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_ident_scalar(2, Scalar{1.0});
  EXPECT_TRUE(
      is_identity_mat(m123_inv_sym * (m_ident_scalar * m123_sym), tolerance));
  EXPECT_TRUE(
      is_identity_mat(m123_inv_sym * (m123_sym * m_ident_scalar), tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m123_sym * m123_inv_sym + m_ident_scalar),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(0.5 * (m_ident_scalar + m123_sym * m123_inv_sym),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_scalar + (m123_sym * m123_inv_sym - m_ident_scalar), tolerance));
  EXPECT_TRUE(is_identity_mat(
      m_ident_scalar + (m_ident_scalar - m123_sym * m123_inv_sym), tolerance));
}

TYPED_TEST(SymmetricTest, VsIdentity) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  mat<typename TestFixture::Scalar, mat_structure::identity,
      TestFixture::Alignment, TestFixture::RowCount, TestFixture::RowCount>
      m_ident2(2);
  EXPECT_TRUE(is_identity_mat(m123_inv_sym * (m_ident2 * m123_sym), tolerance));
  EXPECT_TRUE(is_identity_mat(m123_inv_sym * (m123_sym * m_ident2), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m123_sym * m123_inv_sym + m_ident2), tolerance));
  EXPECT_TRUE(
      is_identity_mat(0.5 * (m_ident2 + m123_sym * m123_inv_sym), tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 + (m123_sym * m123_inv_sym - m_ident2),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(m_ident2 + (m_ident2 - m123_sym * m123_inv_sym),
                              tolerance));
}

TYPED_TEST(SymmetricTest, VsNull) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  mat<typename TestFixture::Scalar, mat_structure::nil, TestFixture::Alignment,
      TestFixture::RowCount, TestFixture::RowCount>
      m_zeroes(2, 2);
  EXPECT_TRUE(is_null_mat(m123_sym * m_zeroes, tolerance));
  EXPECT_TRUE(is_null_mat(m_zeroes * m123_sym, tolerance));
  EXPECT_TRUE(
      is_identity_mat(m_zeroes - m123_inv_sym * (-m123_sym), tolerance));
  EXPECT_TRUE(is_identity_mat(m123_inv_sym * m123_sym - m_zeroes, tolerance));
  EXPECT_TRUE(is_identity_mat(m123_inv_sym * m123_sym + m_zeroes, tolerance));
  EXPECT_TRUE(is_identity_mat(m_zeroes + m123_inv_sym * m123_sym, tolerance));
}

template <typename MatSymAndDense>
class SymmetricDenseTest : public ::testing::Test {
 protected:
  using MatSym = std::tuple_element_t<0, MatSymAndDense>;
  using MatDense = std::tuple_element_t<1, MatSymAndDense>;
  using Scalar = mat_value_type_t<MatSym>;
  SymmetricDenseTest()
      : m123_sym(CreateM123Sym<MatSym>()),
        m123_inv_sym(CreateM123SymInv<MatSym>()),
        m1234(CreateM1234<MatDense>()),
        m1234_inv(CreateM1234Inv<MatDense>()) {}

  MatSym m123_sym;
  MatSym m123_inv_sym;
  MatDense m1234;
  MatDense m1234_inv;
  static constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
};

using SymmetricDenseTestTypes = ::testing::Types<
    std::tuple<
        mat<double, mat_structure::symmetric>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::symmetric>,
        mat<double, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<mat<double, mat_structure::symmetric>,
               mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<float, mat_structure::symmetric>,
               mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<double, mat_structure::symmetric>,
               mat<double, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<mat<float, mat_structure::symmetric>,
               mat<float, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<mat<double, mat_structure::symmetric>,
               mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<float, mat_structure::symmetric>,
               mat<float, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<mat<double, mat_structure::symmetric>,
               mat<double, mat_structure::rectangular, mat_alignment::row_major,
                   2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::symmetric>,
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric>,
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::symmetric>,
        mat<double, mat_structure::square, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric>,
        mat<float, mat_structure::square, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::symmetric, mat_alignment::column_major, 2,
            2>,
        mat<double, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major>>,
    std::tuple<
        mat<double, mat_structure::symmetric, mat_alignment::column_major, 2,
            2>,
        mat<double, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major>>,
    std::tuple<mat<double, mat_structure::symmetric,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::square, mat_alignment::column_major>>,
    std::tuple<mat<double, mat_structure::symmetric,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::square, mat_alignment::row_major>>,
    std::tuple<mat<double, mat_structure::symmetric,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::rectangular,
                   mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::column_major, 2,
            2>>,
    std::tuple<mat<double, mat_structure::symmetric,
                   mat_alignment::column_major, 2, 2>,
               mat<double, mat_structure::rectangular, mat_alignment::row_major,
                   2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::rectangular, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::symmetric, mat_alignment::column_major, 2,
            2>,
        mat<double, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::square, mat_alignment::column_major, 2, 2>>,
    std::tuple<
        mat<double, mat_structure::symmetric, mat_alignment::column_major, 2,
            2>,
        mat<double, mat_structure::square, mat_alignment::row_major, 2, 2>>,
    std::tuple<
        mat<float, mat_structure::symmetric, mat_alignment::column_major, 2, 2>,
        mat<float, mat_structure::square, mat_alignment::row_major, 2, 2>>>;

TYPED_TEST_SUITE(SymmetricDenseTest, SymmetricDenseTestTypes);

TYPED_TEST(SymmetricDenseTest, Operators) {
  const auto tolerance = this->tolerance;
  const auto& m123_sym = this->m123_sym;
  const auto& m123_inv_sym = this->m123_inv_sym;
  const auto& m1234 = this->m1234;
  const auto& m1234_inv = this->m1234_inv;
  EXPECT_TRUE(is_identity_mat(m1234_inv * ((m123_sym * m123_inv_sym) * m1234),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(m1234_inv * (m1234 * (m123_sym * m123_inv_sym)),
                              tolerance));
  EXPECT_TRUE(is_identity_mat(
      0.5 * (m1234 * m1234_inv + (m123_sym * m123_inv_sym)), tolerance));
  EXPECT_TRUE(is_identity_mat(
      0.5 * ((m123_sym * m123_inv_sym) + m1234 * m1234_inv), tolerance));
  EXPECT_TRUE(
      is_identity_mat((m123_sym * m123_inv_sym) +
                          (m1234 * m1234_inv - (m123_sym * m123_inv_sym)),
                      tolerance));
  EXPECT_TRUE(
      is_identity_mat((m123_sym * m123_inv_sym) +
                          ((m123_sym * m123_inv_sym) - m1234 * m1234_inv),
                      tolerance));
}

TEST(MatAlg, MatOperatorTests) {
  constexpr double tolerance = std::numeric_limits<double>::epsilon();
  using std::abs;

  mat<double, mat_structure::rectangular> m1234(1.0, 2.0, 3.0, 4.0);
  mat<double, mat_structure::rectangular> m1234_inv(-2.0, 1.0, 1.5, -0.5);
  mat<double, mat_structure::square> m1234_sqr(1.0, 2.0, 3.0, 4.0);
  mat<double, mat_structure::square> m1234_inv_sqr(-2.0, 1.0, 1.5, -0.5);
  mat<double, mat_structure::symmetric> m123_sym(1.0, 2.0, 3.0);
  mat<double, mat_structure::symmetric> m123_inv_sym(-3.0, 2.0, -1.0);
  mat<double, mat_structure::rectangular> m_ident(2, 2, true);
  mat_identity_t<double> m_ident2(2);
  mat_null_t<double> m_zeroes(2, 2);

  mat<double, mat_structure::rectangular> m4321_orig(4.0, 2.0, 3.0, 1.0);
  std::array<double, 4> f4321 = {4.0, 3.0, 2.0, 1.0};
  std::vector<double> v4321(f4321.begin(), f4321.end());
  mat<double> m4321(v4321, 2, 2);
  EXPECT_TRUE((is_null_mat(m4321 - m4321_orig, tolerance)));

  mat<double> m4321_cpy(m4321);
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy = m4321;
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy += mat<double>(2, 2, true);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy += mat_identity_t<double>(2);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig - mat<double>(2, 2, true),
                               tolerance)));

  m4321_cpy -= mat<double>(2, 2, true);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy -= mat_identity_t<double>(2);
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy *= 2.0;
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig - m4321_orig, tolerance)));

  m4321_cpy = m4321_cpy * double(0.5);
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy = m4321_cpy + mat<double>(2, 2, true);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy = m4321_cpy + mat_identity_t<double>(2);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig - mat<double>(2, 2, true),
                               tolerance)));

  m4321_cpy = m4321_cpy - mat<double>(2, 2, true);
  EXPECT_TRUE((is_identity_mat(m4321_cpy - m4321_orig, tolerance)));

  m4321_cpy = m4321_cpy - mat_identity_t<double>(2);
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));
  EXPECT_TRUE((is_null_mat((-m4321_cpy) + m4321_orig, tolerance)));

  m4321_cpy = m4321_cpy * mat_identity_t<double>(2);
  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));

  EXPECT_TRUE((m4321_cpy == m4321));
  EXPECT_TRUE((m4321_cpy != mat<double>(2, 2, true)));
  EXPECT_TRUE((m4321_cpy != mat_identity_t<double>(2)));

  EXPECT_TRUE((is_null_mat(m4321_cpy - m4321_orig, tolerance)));
  vect_n<double> v85 = m4321_cpy * vect_n<double>(1.0, 1.0);
  EXPECT_NEAR(v85[0], 6.0, tolerance);
  EXPECT_NEAR(v85[1], 4.0, tolerance);
  vect_n<double> v104 = vect_n<double>(1.0, 1.0) * m4321_cpy;
  EXPECT_NEAR(v104[0], 7.0, tolerance);
  EXPECT_NEAR(v104[1], 3.0, tolerance);

  mat<double, mat_structure::rectangular> m_ident3(3, 3, true);
  set_block(m_ident3, m4321_cpy, 1, 1);
  EXPECT_TRUE(((abs(m_ident3(0, 0) - 1.0) < tolerance) &&
               (abs(m_ident3(0, 1) - 0.0) < tolerance) &&
               (abs(m_ident3(0, 2) - 0.0) < tolerance) &&
               (abs(m_ident3(1, 0) - 0.0) < tolerance) &&
               (abs(m_ident3(1, 1) - 4.0) < tolerance) &&
               (abs(m_ident3(1, 2) - 2.0) < tolerance) &&
               (abs(m_ident3(2, 0) - 0.0) < tolerance) &&
               (abs(m_ident3(2, 1) - 3.0) < tolerance) &&
               (abs(m_ident3(2, 2) - 1.0) < tolerance)));

  m4321_cpy = get_block(m_ident3, 1, 1, 2, 2);
  EXPECT_TRUE(((abs(m4321_cpy(0, 0) - 4.0) < tolerance) &&
               (abs(m4321_cpy(0, 1) - 2.0) < tolerance) &&
               (abs(m4321_cpy(1, 0) - 3.0) < tolerance) &&
               (abs(m4321_cpy(1, 1) - 1.0) < tolerance)));
  mat<double, mat_structure::rectangular> m4321_trans((transpose(m4321_cpy)));
  EXPECT_TRUE(((abs(m4321_trans(0, 0) - 4.0) < tolerance) &&
               (abs(m4321_trans(0, 1) - 3.0) < tolerance) &&
               (abs(m4321_trans(1, 0) - 2.0) < tolerance) &&
               (abs(m4321_trans(1, 1) - 1.0) < tolerance)));
  mat<double, mat_structure::symmetric> m4321_sym =
      mat<double, mat_structure::symmetric>(m4321_cpy);
  EXPECT_TRUE(((abs(m4321_sym(0, 0) - 4.0) < tolerance) &&
               (abs(m4321_sym(0, 1) - 2.5) < tolerance) &&
               (abs(m4321_sym(1, 0) - 2.5) < tolerance) &&
               (abs(m4321_sym(1, 1) - 1.0) < tolerance)));
  mat<double, mat_structure::skew_symmetric> m4321_skew =
      mat<double, mat_structure::skew_symmetric>(m4321_cpy);
  const mat<double, mat_structure::skew_symmetric>& m4321_skew_ref = m4321_skew;
  EXPECT_TRUE(((abs(m4321_skew_ref(0, 0) - 0.0) < tolerance) &&
               (abs(m4321_skew_ref(0, 1) + 0.5) < tolerance) &&
               (abs(m4321_skew_ref(1, 0) - 0.5) < tolerance) &&
               (abs(m4321_skew_ref(1, 1) - 0.0) < tolerance)));

  mat<double, mat_structure::diagonal> m123(vect_n<double>(1.0, 2.0, 3.0));
  const mat<double, mat_structure::diagonal>& m123_ref = m123;
  EXPECT_TRUE(
      ((abs(m123_ref(0, 0) - 1.0) < tolerance) &&
       (abs(m123_ref(1, 1) - 2.0) < tolerance) &&
       (abs(m123_ref(2, 2) - 3.0) < tolerance) &&
       (abs(m123_ref(0, 1)) < tolerance) && (abs(m123_ref(0, 2)) < tolerance) &&
       (abs(m123_ref(1, 0)) < tolerance) && (abs(m123_ref(1, 2)) < tolerance) &&
       (abs(m123_ref(2, 0)) < tolerance) && (abs(m123_ref(2, 1)) < tolerance)));

  mat<double, mat_structure::skew_symmetric> m123_skew =
      mat<double, mat_structure::skew_symmetric>(vect_n<double>(1.0, 2.0, 3.0));
  const mat<double, mat_structure::skew_symmetric>& m123_skew_ref = m123_skew;
  EXPECT_TRUE(((abs(m123_skew_ref(0, 0)) < tolerance) &&
               (abs(m123_skew_ref(1, 1)) < tolerance) &&
               (abs(m123_skew_ref(2, 2)) < tolerance) &&
               (abs(m123_skew_ref(0, 1) + 3.0) < tolerance) &&
               (abs(m123_skew_ref(0, 2) - 2.0) < tolerance) &&
               (abs(m123_skew_ref(1, 0) - 3.0) < tolerance) &&
               (abs(m123_skew_ref(1, 2) + 1.0) < tolerance) &&
               (abs(m123_skew_ref(2, 0) + 2.0) < tolerance) &&
               (abs(m123_skew_ref(2, 1) - 1.0) < tolerance)));
};

}  // namespace
}  // namespace ReaK
