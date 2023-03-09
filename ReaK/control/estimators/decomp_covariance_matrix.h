/**
 * \file decomp_covariance_matrix.h
 *
 * This library provides a class template to represent a covariance matrix as a
 * decomposition of the covariance matrix.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_CONTROL_CONTROLLERS_DECOMP_COVARIANCE_MATRIX_H_
#define REAK_CONTROL_CONTROLLERS_DECOMP_COVARIANCE_MATRIX_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_gaussian_elim.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/control/estimators/covariance_concept.h"

namespace ReaK::ctrl {

/**
 * This class template represents a covariance matrix type by holding
 * a decomposition of the covariance matrix. A common approach in estimation
 * is to decompose the covariance matrix P into a covarying block X and a informing-inverse
 * block Y, such that the following relation applies: P = X * invert(Y). This is used, for
 * example, to implement filters in header-file symplectic_kalman_filter.hpp
 * or header-file aggregate_kalman_filter.hpp .
 *
 * Models: CovarianceMatrixConcept and DecomposedCovarianceConcept.
 *
 * \tparam VectorType The state-vector type which the covariance matrix is the covariance of, should model
 *ReadableVectorConcept.
 */
template <typename VectorType>
class decomp_covariance_matrix : public named_object {
 public:
  BOOST_CONCEPT_ASSERT((ReadableVectorConcept<VectorType>));

  using self = decomp_covariance_matrix<VectorType>;

  using value_type = vect_value_type_t<VectorType>;
  using matrix_type = mat<value_type, mat_structure::symmetric>;
  using size_type = typename matrix_type::size_type;

  using matrix_block_type = mat<value_type, mat_structure::square>;

  static constexpr std::size_t dimensions = vect_traits<VectorType>::dimensions;
  static constexpr covariance_storage::tag storage = covariance_storage::other;

 private:
  matrix_block_type mat_X;
  matrix_block_type mat_Y;

 public:
  /**
   * Parametrized constructor.
   * \param aMatX The covarying matrix block.
   * \param aMatY The informing-inverse matrix block.
   */
  explicit decomp_covariance_matrix(
      const matrix_block_type& aMatX,
      const matrix_block_type& aMatY = matrix_block_type(),
      const std::string& aName = "")
      : mat_X(aMatX), mat_Y(aMatY) {
    setName(aName);
  }

  decomp_covariance_matrix() : decomp_covariance_matrix(matrix_block_type()) {}

  /**
   * Parametrized constructor.
   * \param aSize The size of the covariance matrix.
   * \param aLevel The information level to initialize this object with.
   */
  explicit decomp_covariance_matrix(
      size_type aSize,
      covariance_initial_level aLevel = covariance_initial_level::full_info,
      const std::string& aName = "")
      : mat_X(aSize, value_type(0)), mat_Y(aSize, value_type(0)) {
    setName(aName);
    if (aLevel == covariance_initial_level::full_info) {
      mat_Y = mat<value_type, mat_structure::identity>(aSize);
    } else {
      mat_X = mat<value_type, mat_structure::identity>(aSize);
    }
  }

  /**
   * Returns the covariance matrix (as a matrix object).
   * \return The covariance matrix (as a matrix object).
   */
  matrix_type get_matrix() const {
    try {
      mat<value_type, mat_structure::square> mY_inv;
      invert_PLU(mat_Y, mY_inv);
      return matrix_type(mat_X * mY_inv);
    } catch (singularity_error& e) {
      mat<value_type, mat_structure::square> mY_inv;
      pseudoinvert_QR(mat_Y, mY_inv);
      return matrix_type(mat_X * mY_inv);
    }
  }
  /**
   * Returns the inverse covariance matrix (information matrix) (as a matrix object).
   * \return The inverse covariance matrix (information matrix) (as a matrix object).
   */
  matrix_type get_inverse_matrix() const {
    try {
      mat<value_type, mat_structure::square> mX_inv;
      invert_PLU(mat_X, mX_inv);
      return matrix_type(mat_Y * mX_inv);
    } catch (singularity_error& e) {
      mat<value_type, mat_structure::square> mX_inv;
      pseudoinvert_QR(mat_X, mX_inv);
      return matrix_type(mat_Y * mX_inv);
    }
  }

  /**
   * Returns the covarying matrix block.
   * \return The covarying matrix block.
   */
  const matrix_block_type& get_covarying_block() const { return mat_X; }
  /**
   * Returns the informing-inverse matrix block.
   * \return The informing-inverse matrix block.
   */
  const matrix_block_type& get_informing_inv_block() const { return mat_Y; }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.mat_X, rhs.mat_X);
    swap(lhs.mat_Y, rhs.mat_Y);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(rhs, *this);
    return *this;
  }

  /**
   * Implicit conversion to a covariance matrix type.
   */
  operator matrix_type() const { return get_matrix(); }

  /**
   * Conversion to an information matrix type.
   */
  friend matrix_type invert(const self& aObj) {
    return aObj.get_inverse_matrix();
  }

  /**
   * Returns the size of the covariance matrix.
   * \return The size of the covariance matrix.
   */
  size_type size() const { return mat_X.get_row_count(); }

  void save(ReaK::serialization::oarchive& aA, unsigned int) const override {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(mat_X) & RK_SERIAL_SAVE_WITH_NAME(mat_Y);
  }
  void load(ReaK::serialization::iarchive& aA, unsigned int) override {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(mat_X) & RK_SERIAL_LOAD_WITH_NAME(mat_Y);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC230000A, 1, "decomp_covariance_matrix",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_DECOMP_COVARIANCE_MATRIX_H_
