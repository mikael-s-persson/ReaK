/**
 * \file mat_norms.h
 *
 * This library provides several functions to compute the various matrix norms.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_NORMS_H_
#define REAK_MATH_LIN_ALG_MAT_NORMS_H_

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_traits.h"

#include <type_traits>

namespace ReaK {

/**
 * This function template computes the 1-norm of a matrix.
 * \param M A matrix for which the 1-norm is sought.
 * \return the 1-norm of matrix M.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> norm_1(const Matrix& M) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  ValueType max = ValueType();
  for (int j = 0; j < M.get_col_count(); ++j) {
    ValueType sum = ValueType();
    for (int i = 0; i < M.get_row_count(); ++i) {
      sum += abs(M(i, j));
    }
    if (sum > max) {
      max = sum;
    }
  }
  return max;
}

/**
 * This function template computes the infinity-norm of a matrix.
 * \param M A matrix for which the infinity-norm is sought.
 * \return the infinity-norm of matrix M.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> norm_inf(const Matrix& M) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  ValueType max = ValueType();
  for (int i = 0; i < M.get_row_count(); ++i) {
    ValueType sum = ValueType();
    for (int j = 0; j < M.get_col_count(); ++j) {
      sum += abs(M(i, j));
    }
    if (sum > max) {
      max = sum;
    }
  }
  return max;
}

/**
 * This function template computes the element-wise 2-norm of a matrix.
 * \param M A matrix for which the element-wise 2-norm is sought.
 * \return the element-wise 2-norm of matrix M.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> elem_norm_2(const Matrix& M) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::sqrt;

  ValueType sum = ValueType();
  for (int i = 0; i < M.get_row_count(); ++i) {
    for (int j = 0; j < M.get_col_count(); ++j) {
      sum += M(i, j) * M(i, j);
    }
  }
  return sqrt(sum);
}

/**
 * This function template computes the Frobenius-norm of a matrix.
 * \param M A matrix for which the Frobenius-norm is sought.
 * \return the Frobenius-norm of matrix M.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> frobenius_norm(const Matrix& M) {
  return elem_norm_2(M);
}

/**
 * This function template computes the element-wise infinity-norm of a matrix.
 * \param M A matrix for which the element-wise infinity-norm is sought.
 * \return the element-wise infinity-norm of matrix M.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> elem_norm_max(const Matrix& M) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  ValueType max = ValueType();
  for (int i = 0; i < M.get_row_count(); ++i) {
    for (int j = 0; j < M.get_col_count(); ++j) {
      if (max < abs(M(i, j))) {
        max = abs(M(i, j));
      }
    }
  }
  return max;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_NORMS_H_
