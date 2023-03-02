/**
 * \file hessian_approx_update.hpp
 *
 * The following library provides methods to perform an update on an approximated Hessian or
 * inverse Hessian matrix (for quasi-Newton methods).
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_HESSIAN_APPROX_UPDATE_HPP
#define REAK_HESSIAN_APPROX_UPDATE_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include <type_traits>
#include <utility>

namespace ReaK::optim {

/**
 * This function updates an approximate Hessian matrix B using the symmetric rank-one update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param B The current approximation of the Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void sr1_hessian_update(Matrix& B, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector r = dy;
  r -= B * dx;
  ValueType denom = r * dx;

  if (abs(denom) < std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(r)) {
    throw singularity_error(
        "SR-1 Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < B.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      B(j, i) = (B(i, j) += r[i] * r[j] / denom);
    }
  }
}

struct hessian_update_sr1 {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& B, const Vector& dx, const Vector& dy) const {
    sr1_hessian_update(B, dx, dy);
  }
};

/**
 * This function updates an approximate inverse Hessian matrix H using the symmetric rank-one update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param H The current approximation of the inverse Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void sr1_inv_hessian_update(Matrix& H, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector r = dx;
  r -= H * dy;
  ValueType denom = r * dy;

  if (abs(denom) < std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(r)) {
    throw singularity_error(
        "SR-1 Inverse Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < H.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      H(j, i) = (H(i, j) += r[i] * r[j] / denom);
    }
  }
}

struct inv_hessian_update_sr1 {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& H, const Vector& dx, const Vector& dy) const {
    sr1_inv_hessian_update(H, dx, dy);
  }
};

/**
 * This function updates an approximate Hessian matrix B using the Davidon-Fletcher-Powell update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param B The current approximation of the Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void dfp_hessian_update(Matrix& B, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Bx = B * dx;
  ValueType denom = dy * dx;
  Vector xBxy = (dx * Bx) * dy;

  if (abs(denom) < std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dy)) {
    throw singularity_error(
        "DFP Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < B.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      B(j, i) = (B(i, j) += (dy[i] * xBxy[j] / denom - dy[i] * Bx[j] -
                             Bx[i] * dy[j] + dy[i] * dy[j]) /
                            denom);
    }
  }
}

struct hessian_update_dfp {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& B, const Vector& dx, const Vector& dy) const {
    dfp_hessian_update(B, dx, dy);
  }
};

/**
 * This function updates an approximate inverse Hessian matrix H using the Davidon-Fletcher-Powell update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param H The current approximation of the inverse Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void dfp_inv_hessian_update(Matrix& H, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Hy = H * dy;
  ValueType denom1 = dy * dx;
  ValueType denom2 = dy * Hy;

  if ((abs(denom1) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dx)) ||
      (abs(denom2) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(Hy))) {
    throw singularity_error(
        "DFP Inverse Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < H.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      H(j, i) = (H(i, j) += dx[i] * dx[j] / denom1 - Hy[i] * Hy[j] / denom2);
    }
  }
}

struct inv_hessian_update_dfp {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& H, const Vector& dx, const Vector& dy) const {
    dfp_inv_hessian_update(H, dx, dy);
  }
};

/**
 * This function updates an approximate Hessian matrix B using the
 * Broyden-Fletcher-Goldfarb-Shanno (BFGS) update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param B The current approximation of the Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void bfgs_hessian_update(Matrix& B, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Bx = B * dx;
  ValueType denom1 = dy * dx;
  ValueType denom2 = dx * Bx;

  if ((abs(denom1) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dy)) ||
      (abs(denom2) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(Bx))) {
    throw singularity_error(
        "BFGS Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < B.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      B(j, i) = (B(i, j) += dy[i] * dy[j] / denom1 - Bx[i] * Bx[j] / denom2);
    }
  }
}

struct hessian_update_bfgs {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& B, const Vector& dx, const Vector& dy) const {
    bfgs_hessian_update(B, dx, dy);
  }
};

/**
 * This function updates an approximate inverse Hessian matrix H using the
 * Broyden-Fletcher-Goldfarb-Shanno (BFGS) update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param H The current approximation of the inverse Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void bfgs_inv_hessian_update(Matrix& H, const Vector& dx, const Vector& dy) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Hy = H * dy;
  ValueType denom = dy * dx;
  Vector yHyx = (dy * Hy) * dx;

  if (abs(denom) < std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dx)) {
    throw singularity_error(
        "BFGS Inverse Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < H.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      H(j, i) = (H(i, j) += (dx[i] * yHyx[j] / denom - dx[i] * Hy[j] -
                             Hy[i] * dx[j] + dx[i] * dx[j]) /
                            denom);
    }
  }
}

struct inv_hessian_update_bfgs {
  template <typename Matrix, typename Vector>
  void operator()(Matrix& H, const Vector& dx, const Vector& dy) const {
    bfgs_inv_hessian_update(H, dx, dy);
  }
};

/**
 * This function updates an approximate Hessian matrix B using the Broyden-class update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param B The current approximation of the Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 */
template <typename Matrix, typename Vector>
void broyden_class_hessian_update(
    Matrix& B, const Vector& dx, const Vector& dy,
    vect_value_type_t<Vector> phi = vect_value_type_t<Vector>(0.5)) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Bx = B * dx;
  ValueType denom1 = dy * dx;
  ValueType denom2 = dx * Bx;
  Vector xBxy = denom2 * dy;

  if ((abs(denom1) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dy)) ||
      (abs(denom2) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(Bx))) {
    throw singularity_error(
        "Broyden-class Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < B.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      B(j, i) = (B(i, j) += (ValueType(1.0) - phi) * (dy[i] * dy[j] / denom1 -
                                                      Bx[i] * Bx[j] / denom2) +
                            phi *
                                (dy[i] * xBxy[j] / denom1 - dy[i] * Bx[j] -
                                 Bx[i] * dy[j] + dy[i] * dy[j]) /
                                denom1);
    }
  }
}

template <typename T>
struct hessian_update_broyden {
  T phi;
  explicit hessian_update_broyden(const T& aPhi = T(0.5)) : phi(aPhi) {}

  template <typename Matrix, typename Vector>
  void operator()(Matrix& B, const Vector& dx, const Vector& dy) const {
    broyden_class_hessian_update(B, dx, dy, phi);
  }
};

/**
 * This function updates an approximate inverse Hessian matrix H using the Broyden-class update.
 * TEST PASSED
 * \tparam Matrix A fully writable matrix type.
 * \tparam Vector A readable vector type.
 * \param H The current approximation of the inverse Hessian matrix.
 * \param dx The change in the independent vector.
 * \param dy The change in the gradient vector of the function whose Hessian is approximated.
 * \param phi The fraction defining the linear combination (from 0: BFGS update only, to 1: DFP update only).
 */
template <typename Matrix, typename Vector>
void broyden_class_inv_hessian_update(
    Matrix& H, const Vector& dx, const Vector& dy,
    vect_value_type_t<Vector> phi = vect_value_type_t<Vector>(0.5)) {
  static_assert(is_readable_vector_v<Vector>);
  static_assert(is_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using SizeType = mat_size_type_t<Matrix>;
  using std::abs;

  Vector Hy = H * dy;
  ValueType denom1 = dy * dx;
  ValueType denom2 = dy * Hy;
  Vector yHyx = denom2 * dx;

  if ((abs(denom1) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(dx)) ||
      (abs(denom2) <
       std::numeric_limits<ValueType>::epsilon() * norm_2_sqr(Hy))) {
    throw singularity_error(
        "Broyden-class Inverse Hessian update has detected a singular update!");
  }

  for (SizeType i = 0; i < H.get_row_count(); ++i) {
    for (SizeType j = 0; j <= i; ++j) {
      H(j, i) = (H(i, j) += (ValueType(1.0) - phi) * (dx[i] * dx[j] / denom1 -
                                                      Hy[i] * Hy[j] / denom2) +
                            phi *
                                (dx[i] * yHyx[j] / denom1 - dx[i] * Hy[j] -
                                 Hy[i] * dx[j] + dx[i] * dx[j]) /
                                denom1);
    }
  }
}

template <typename T>
struct inv_hessian_update_broyden {
  T phi;
  explicit inv_hessian_update_broyden(const T& aPhi = T(0.5)) : phi(aPhi) {}

  template <typename Matrix, typename Vector>
  void operator()(Matrix& H, const Vector& dx, const Vector& dy) const {
    broyden_class_inv_hessian_update(H, dx, dy, phi);
  }
};

template <typename HessianFunction>
struct hessian_update_dual_exact {
  HessianFunction fill_hessian;
  explicit hessian_update_dual_exact(HessianFunction aFillHessian)
      : fill_hessian(std::move(aFillHessian)) {}
  template <typename Matrix, typename T, typename Vector>
  void operator()(Matrix& B, const Vector& x, const T& x_value,
                  const Vector& x_grad, const Vector& /*unused*/,
                  const Vector& /*unused*/) const {
    fill_hessian(B, x, x_value, x_grad);
  }
  template <typename Matrix, typename T, typename Vector>
  void operator()(Matrix& B, const Vector& x, const T& x_value,
                  const Vector& x_grad) const {
    fill_hessian(B, x, x_value, x_grad);
  }
};

template <typename HessianUpdater>
struct hessian_update_dual_quasi {
  HessianUpdater update_hessian;
  explicit hessian_update_dual_quasi(HessianUpdater aUpdateHessian)
      : update_hessian(aUpdateHessian) {}
  template <typename Matrix, typename T, typename Vector>
  void operator()(Matrix& B, const Vector& /*unused*/, const T& /*unused*/,
                  const Vector& /*unused*/, const Vector& dx,
                  const Vector& dy) const {
    update_hessian(B, dx, dy);
  }
  template <typename Matrix, typename T, typename Vector>
  void operator()(Matrix& /*unused*/, const Vector& /*unused*/,
                  const T& /*unused*/, const Vector& /*unused*/) const {}
};

}  // namespace ReaK::optim

#endif
