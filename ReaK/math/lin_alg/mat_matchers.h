/// \file mat_matchers.h
///
/// This library provides several googletest matchers for the various matrix comparisons.
///
/// \author Sven Mikael Persson <mikael.s.persson@gmail.com>
/// \date March 2023

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

#ifndef REAK_MATH_LIN_ALG_MAT_MATCHERS_H_
#define REAK_MATH_LIN_ALG_MAT_MATCHERS_H_

#include "ReaK/math/lin_alg/mat_comparisons.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "gmock/gmock.h"

namespace ReaK::testing {

/// This function template computes the element-wise comparison of two matrices.
MATCHER_P2(MatrixIsNear, expected, tolerance, "") {
  if (!is_equal_mat(arg, expected, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << expected << "\nActual matrix =\n"
                     << arg << "\nDifference =\n"
                     << (arg - expected) << "\nTolerance = " << tolerance
                     << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is diagonal up to a tolerance.
MATCHER_P(MatrixIsDiagonal, tolerance, "") {
  if (!is_diagonal(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be diagonal, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is symmetric up to a tolerance.
MATCHER_P(MatrixIsSymmetric, tolerance, "") {
  if (!is_symmetric(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be symmetric, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is the identity up to a tolerance.
MATCHER_P(MatrixIsIdentity, tolerance, "") {
  if (!is_identity_mat(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be identity, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is null up to a tolerance.
MATCHER_P(MatrixIsNull, tolerance, "") {
  if (!is_null_mat(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be null, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is orthogonal up to a tolerance.
MATCHER_P(MatrixIsOrthogonal, tolerance, "") {
  if (!is_identity_mat(arg * transpose(arg), tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be orthogonal, but it is NOT!\n"
                     << "Matrix product with its transpose is =\n"
                     << (arg * transpose(arg)) << "\nTolerance = " << tolerance
                     << "\n";
    return false;
  }
  return true;
}

/// Verifies if two matrices are inverse of each other.
MATCHER_P2(MatrixIsInverseOf, arg_inv, tolerance, "") {
  if (!is_identity_mat(arg * arg_inv, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be the inverse of matrix = \n"
                     << arg_inv << ", but it is NOT!\n"
                     << "Matrix product with its expected inverse is =\n"
                     << (arg * arg_inv) << "\nTolerance = " << tolerance
                     << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is upper-triangular up to a tolerance.
MATCHER_P(MatrixIsUpperTriangular, tolerance, "") {
  if (!is_upper_triangular(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be upper triangular, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is lower-triangular up to a tolerance.
MATCHER_P(MatrixIsLowerTriangular, tolerance, "") {
  if (!is_lower_triangular(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be lower triangular, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is upper-hessenberg up to a tolerance.
MATCHER_P(MatrixIsUpperHessenberg, tolerance, "") {
  if (!is_upper_hessenberg(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be upper hessenberg, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is lower-hessenberg up to a tolerance.
MATCHER_P(MatrixIsLowerHessenberg, tolerance, "") {
  if (!is_lower_hessenberg(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be lower hessenberg, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

/// Verifies that the matrix A is tri-diagonal up to a tolerance.
MATCHER_P(MatrixIsTriDiagonal, tolerance, "") {
  if (!is_tri_diagonal(arg, tolerance)) {
    *result_listener << "\nExpected matrix =\n"
                     << arg << "\n  to be tri-diagonal, but it is NOT!"
                     << "\nTolerance = " << tolerance << "\n";
    return false;
  }
  return true;
}

}  // namespace ReaK::testing

#endif  // REAK_MATH_LIN_ALG_MAT_MATCHERS_H_
