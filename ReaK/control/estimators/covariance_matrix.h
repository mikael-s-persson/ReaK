/**
 * \file covariance_matrix.h
 *
 * This library provides a class template to represent a covariance matrix as a
 * covariance matrix.
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

#ifndef REAK_CONTROL_CONTROLLERS_COVARIANCE_MATRIX_H_
#define REAK_CONTROL_CONTROLLERS_COVARIANCE_MATRIX_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"

#include "ReaK/control/estimators/covariance_concept.h"

#include <type_traits>
#include <utility>

namespace ReaK::ctrl {

/**
 * This class template can represent a covariance matrix (as of the CovarianceMatrixConcept)
 * by containing the covariance matrix. Outputing
 * the covariance matrix induces no computation while outputing the information matrix induces
 * a matrix inversion.
 *
 * Models: CovarianceMatrixConcept
 *
 * \tparam VectorType The state-vector type which the covariance matrix is the covariance of, should model
 *ReadableVectorConcept.
 */
template <typename VectorType>
class covariance_matrix : public named_object {
 public:
  BOOST_CONCEPT_ASSERT((ReadableVectorConcept<VectorType>));

  using self = covariance_matrix<VectorType>;

  using value_type = vect_value_type_t<VectorType>;
  using matrix_type = mat<value_type, mat_structure::symmetric>;
  using size_type = typename matrix_type::size_type;

  static constexpr std::size_t dimensions = vect_traits<VectorType>::dimensions;
  static constexpr covariance_storage::tag storage =
      covariance_storage::covariance;

 private:
  matrix_type mat_cov;

 public:
  /**
   * Parametrized constructor.
   * \param aMat The covariance matrix to which to initialize this object.
   */
  template <typename Matrix>
  explicit covariance_matrix(Matrix aMat, const std::string& aName = "")
      : mat_cov(std::move(aMat)) {
    static_assert(is_readable_matrix_v<Matrix>);
    setName(aName);
  }

  /**
   * Parametrized constructor.
   * \param aSize The size of the covariance matrix.
   * \param aLevel The information level to initialize this object with.
   */
  explicit covariance_matrix(int aSize,
                             covariance_initial_level::tag aLevel =
                                 covariance_initial_level::full_info,
                             const std::string& aName = "")
      : mat_cov(
            aSize,
            value_type((aLevel == covariance_initial_level::full_info
                            ? 0
                            : std::numeric_limits<value_type>::infinity()))) {
    setName(aName);
  }

  covariance_matrix() : covariance_matrix(0) {}

  /**
   * Returns the covariance matrix (as a matrix object).
   * \return The covariance matrix (as a matrix object).
   */
  const matrix_type& get_matrix() const { return mat_cov; }
  /**
   * Returns the inverse covariance matrix (information matrix) (as a matrix object).
   * \return The inverse covariance matrix (information matrix) (as a matrix object).
   */
  matrix_type get_inverse_matrix() const {
    matrix_type m_inv;
    invert_Cholesky(mat_cov, m_inv, std::numeric_limits<value_type>::epsilon());
    return m_inv;
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.mat_cov, rhs.mat_cov);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(rhs, *this);
    return *this;
  }

  /**
   * Assignment to a readable matrix (covariance matrix).
   */
  template <typename Matrix>
  self& operator=(const Matrix& rhs) {
    static_assert(is_readable_matrix_v<Matrix>);
    mat_cov = rhs;
    return *this;
  }

  /**
   * Implicit conversion to a covariance matrix type.
   */
  explicit operator const matrix_type&() const { return mat_cov; }

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
  size_type size() const { return mat_cov.get_row_count(); }

  void save(ReaK::serialization::oarchive& aA,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(mat_cov);
  }
  void load(ReaK::serialization::iarchive& aA,
            unsigned int /*unused*/) override {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(mat_cov);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300008, 1, "covariance_matrix",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_COVARIANCE_MATRIX_H_
