/**
 * \file gaussian_belief_state.h
 *
 * This library provides a number of class templates that can be used to represent and
 * use a Gaussian belief-state (i.e. a mean-state and a covariance matrix). Those class
 * templates include the belief-state itself, of course, a sampler to generate random
 * states from the belief-state's probability distribution, and a PDF (probability
 * distribution function) to compute the probability of a given state vector.
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

#ifndef REAK_CONTROL_CONTROLLERS_GAUSSIAN_BELIEF_STATE_H_
#define REAK_CONTROL_CONTROLLERS_GAUSSIAN_BELIEF_STATE_H_

#include "ReaK/core/base/global_rng.h"
#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"

#include "ReaK/control/estimators/belief_state_concept.h"
#include "ReaK/control/estimators/covariance_concept.h"

#include <random>
#include <type_traits>
#include <utility>

namespace ReaK::ctrl {

/**
 * This class template is a functor that can compute the probability that a given state is part of
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * \tparam Storage The storage strategy of the Covariance matrix, this is used for specializing the PDF implementation
 * for the most efficient way to compute the probabilities.
 */
template <ContinuousBeliefState BState,
          covariance_storage::tag Storage =
              covariance_mat_traits<typename continuous_belief_state_traits<
                  BState>::covariance_type>::storage>
struct gaussian_pdf {
  using self = gaussian_pdf<BState, Storage>;

  using state_type =
      typename continuous_belief_state_traits<BState>::state_type;
  using covariance_type =
      typename continuous_belief_state_traits<BState>::covariance_type;

  using scalar_type =
      typename covariance_mat_traits<covariance_type>::value_type;
  using size_type = typename covariance_mat_traits<covariance_type>::size_type;
  using matrix_type =
      typename covariance_mat_traits<covariance_type>::matrix_type;

  state_type mean_state;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> L;
  scalar_type factor;

  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  explicit gaussian_pdf(const BState& aBelief)
      : mean_state(aBelief.get_mean_state()),
        L(aBelief.get_covariance().size()),
        factor(-1) {
    const matrix_type& E = aBelief.get_covariance().get_matrix();
    try {
      decompose_Cholesky(E, L);
    } catch (singularity_error&) {
      return;
    }
    factor = scalar_type(1);
    for (size_type i = 0; i < L.get_row_count(); ++i) {
      factor *= scalar_type(6.28318530718) * L(i, i);
    }
  }

  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs)
      : mean_state(rhs.mean_state), L(rhs.L), factor(rhs.factor) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.mean_state, rhs.mean_state);
    swap(lhs.L, rhs.L);
    swap(lhs.factor, rhs.factor);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <pp::Topology Space>
  scalar_type operator()(const state_type& v, const Space& space) const {
    using ReaK::to_vect;
    using std::exp;
    using std::sqrt;

    if (factor <= scalar_type(0)) {
      return scalar_type(0);
    }

    vect_n<scalar_type> d =
        to_vect<scalar_type>(space.difference(v, mean_state));
    mat<mat_value_type_t<matrix_type>, mat_structure::rectangular> b(d.size(),
                                                                     1);
    for (size_type i = 0; i < d.size(); ++i) {
      b(i, 0) = d[i];
    }
    ::ReaK::detail::backsub_Cholesky_impl(L, b);
    scalar_type sum = scalar_type(0);
    for (size_type i = 0; i < d.size(); ++i) {
      sum += d[i] * b(i, 0);
    }
    return exp(scalar_type(-0.5) * sum) / sqrt(factor);
  }

  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) *
           (log(P.factor) + scalar_type(P.L.get_row_count()));
  }
};

/**
 * This class template is a functor that can compute the probability that a given state is part of
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * This class template specialization uses the fact that the covariance is represented as an information
 * matrix in order to have a more efficient implementation.
 */
template <ContinuousBeliefState BState>
struct gaussian_pdf<BState, covariance_storage::information> {
  using self = gaussian_pdf<BState, covariance_storage::information>;

  using state_type =
      typename continuous_belief_state_traits<BState>::state_type;
  using covariance_type =
      typename continuous_belief_state_traits<BState>::covariance_type;

  using scalar_type =
      typename covariance_mat_traits<covariance_type>::value_type;
  using size_type = typename covariance_mat_traits<covariance_type>::size_type;
  using matrix_type =
      typename covariance_mat_traits<covariance_type>::matrix_type;

  state_type mean_state;
  matrix_type E_inv;
  scalar_type factor;

  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  explicit gaussian_pdf(const BState& aBelief)
      : mean_state(aBelief.get_mean_state()),
        E_inv(aBelief.get_covariance().get_inverse_matrix()),
        factor(-1) {
    factor = determinant_Cholesky(E_inv);
    if (abs(factor) < std::numeric_limits<scalar_type>::epsilon()) {
      factor = scalar_type(-1);
    } else {
      factor = scalar_type(1) / factor;
      for (size_type i = 0; i < E_inv.get_row_count(); ++i) {
        factor *= scalar_type(6.28318530718);
      }
    }
  }

  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs)
      : mean_state(rhs.mean_state), E_inv(rhs.E_inv), factor(rhs.factor) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.mean_state, rhs.mean_state);
    swap(lhs.E_inv, rhs.E_inv);
    swap(lhs.factor, rhs.factor);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <pp::Topology Space>
  scalar_type operator()(const state_type& v, const Space& space) const {
    using ReaK::to_vect;
    using std::abs;
    using std::exp;
    using std::sqrt;
    using state_difference_type = pp::topology_point_difference_type_t<Space>;
    static_assert(ReadableVector<state_difference_type>);
    static_assert(CovarianceMatrix<covariance_type, state_difference_type>);

    if (factor <= scalar_type(0)) {
      return scalar_type(0);
    }

    vect_n<scalar_type> d =
        to_vect<scalar_type>(space.difference(v, mean_state));
    return exp(scalar_type(-0.5) * (d * (E_inv * d))) / sqrt(factor);
  }

  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) *
           (log(P.factor) + scalar_type(P.E_inv.get_row_count()));
  }
};

/**
 * This class template is a functor that can compute the probability that a given state is part of
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * This class template specialization uses the fact that the covariance is represented as a decomposition of
 * the covariance matrix in order to have a more efficient implementation.
 */
template <ContinuousBeliefState BState>
struct gaussian_pdf<BState, covariance_storage::decomposed> {
  using self = gaussian_pdf<BState, covariance_storage::decomposed>;

  using state_type =
      typename continuous_belief_state_traits<BState>::state_type;
  using covariance_type =
      typename continuous_belief_state_traits<BState>::covariance_type;

  using scalar_type =
      typename covariance_mat_traits<covariance_type>::value_type;
  using size_type = typename covariance_mat_traits<covariance_type>::size_type;
  using matrix_type =
      typename covariance_mat_traits<covariance_type>::matrix_type;

  state_type mean_state;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> QX;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> RX;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> QY;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> RY;
  scalar_type factor;

  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  explicit gaussian_pdf(const BState& aBelief)
      : mean_state(aBelief.get_mean_state()), factor(-1) {
    decompose_QR(aBelief.get_covariance().get_covarying_block(), QX, RX);
    decompose_QR(aBelief.get_covariance().get_informing_inv_block(), QY, RY);

    factor = scalar_type(1);
    for (size_type i = 0; i < RX.get_row_count(); ++i) {
      factor *= scalar_type(6.28318530718) * RX(i, i) / RY(i, i);
    }
  }

  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs)
      : mean_state(rhs.mean_state),
        QX(rhs.QX),
        RX(rhs.RX),
        QY(rhs.QY),
        RY(rhs.RY),
        factor(rhs.factor) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.mean_state, rhs.mean_state);
    swap(lhs.QX, rhs.QX);
    swap(lhs.RX, rhs.RX);
    swap(lhs.QY, rhs.QY);
    swap(lhs.RY, rhs.RY);
    swap(lhs.factor, rhs.factor);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <pp::Topology Space>
  scalar_type operator()(const state_type& v, const Space& space) const {
    using ReaK::to_vect;
    using std::exp;
    using std::sqrt;
    using state_difference_type = pp::topology_point_difference_type_t<Space>;
    static_assert(WritableVector<state_difference_type>);
    static_assert(CovarianceMatrix<covariance_type, state_difference_type>);

    if (factor <= scalar_type(0)) {
      return scalar_type(0);
    }

    vect_n<scalar_type> d =
        to_vect<scalar_type>(space.difference(v, mean_state));
    vect_n<scalar_type> d_tmp = d * QX;  // QX^T d
    mat_vect_adaptor<state_difference_type> d_m(d_tmp);
    backsub_R(RX, d_m);
    return exp(scalar_type(-0.5) * (d * (QY * (RY * d_tmp)))) / sqrt(factor);
  }

  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) *
           (log(P.factor) + scalar_type(P.QX.get_row_count()));
  }
};

template <typename Vector, typename Matrix>
mat_value_type_t<Matrix> gaussian_pdf_at_diff(const Vector& dx,
                                              const Matrix& P) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::exp;
  using std::pow;
  using std::sqrt;

  ValueType det_sqrt = 1.0;
  Vector Pinv_dx = dx;
  try {
    mat<ValueType, mat_structure::square> L;
    decompose_Cholesky(P, L);
    for (int i = 0; i < L.get_col_count(); ++i) {
      det_sqrt *= L(i, i);
    }
    mat_vect_adaptor<Vector> Pinv_dx_mat(Pinv_dx);
    ::ReaK::detail::backsub_Cholesky_impl(L, Pinv_dx_mat);
  } catch (singularity_error&) {
    mat<double, mat_structure::diagonal> E(P.get_row_count());
    mat<double, mat_structure::square> U(P.get_row_count());
    mat<double, mat_structure::square> V(P.get_row_count());
    decompose_SVD(P, U, E, V);
    for (int i = 0; i < E.get_col_count(); ++i) {
      det_sqrt *= E(i, i);
    }
    det_sqrt = sqrt(det_sqrt);
    mat<ValueType, mat_structure::square> Pinv(P.get_row_count());
    pseudoinvert_SVD(U, E, V, Pinv, 1 - 8);
    Pinv_dx = Pinv * dx;
  }

  ValueType dist_x = -0.5 * (dx * Pinv_dx);
  ValueType factor = pow(2.0 * M_PI, -0.5 * P.get_row_count()) / det_sqrt;
  return factor * exp(dist_x);
}

template <typename Matrix>
mat_value_type_t<Matrix> gaussian_pdf_at_mean_value(const Matrix& P) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::pow;
  using std::sqrt;

  ValueType det_sqrt = 1.0;
  try {
    mat<ValueType, mat_structure::square> L;
    decompose_Cholesky(P, L);
    for (int i = 0; i < L.get_col_count(); ++i) {
      det_sqrt *= L(i, i);
    }
  } catch (singularity_error&) {
    mat<double, mat_structure::diagonal> E(P.get_row_count());
    mat<double, mat_structure::square> U(P.get_row_count());
    mat<double, mat_structure::square> V(P.get_row_count());
    decompose_SVD(P, U, E, V);
    for (int i = 0; i < E.get_col_count(); ++i) {
      det_sqrt *= E(i, i);
    }
    det_sqrt = sqrt(det_sqrt);
  }
  return pow(2.0 * M_PI, -0.5 * P.get_row_count()) / det_sqrt;
}

template <typename Vector, typename Matrix>
auto gaussian_likelihood_ratio_of_diff(const Vector& dx, const Matrix& P) {
  using ValueType = mat_value_type_t<Matrix>;
  Vector Pinv_dx = dx;
  try {
    mat<ValueType, mat_structure::square> L;
    decompose_Cholesky(P, L);
    mat_vect_adaptor<Vector> Pinv_dx_mat(Pinv_dx);
    ::ReaK::detail::backsub_Cholesky_impl(L, Pinv_dx_mat);
  } catch (singularity_error&) {
    mat<double, mat_structure::diagonal> E(P.get_row_count());
    mat<double, mat_structure::square> U(P.get_row_count());
    mat<double, mat_structure::square> V(P.get_row_count());
    decompose_SVD(P, U, E, V);
    mat<ValueType, mat_structure::square> Pinv(P.get_row_count());
    pseudoinvert_SVD(U, E, V, Pinv, 1 - 8);
    Pinv_dx = Pinv * dx;
  }

  ValueType dist_x = dx * Pinv_dx;
  return dist_x;
}

/**
 * This function template computes the symmetric KL-divergence between two Gaussian probability
 * distribution function objects.
 * \tparam BeliefState The belief-state type for the Gaussian PDFs.
 * \tparam Storage The storage strategy for the covariance matrices.
 * \tparam Space The topology type of the underlying state representations.
 * \param N0 The first PDF.
 * \param N1 The second PDF.
 * \param space The topology of the underlying state representations.
 * \return The symmetric KL-divergence between the two Gaussian PDFs.
 */
template <ContinuousBeliefState BState, covariance_storage::tag Storage,
          pp::Topology Space>
auto symKL_divergence(const gaussian_pdf<BState, Storage>& N0,
                      const gaussian_pdf<BState, Storage>& N1,
                      const Space& space) {
  using std::log;
  return -0.5 * log(N1(N0.mean_state, space) * N0(N1.mean_state, space)) -
         entropy(N0) - entropy(N1);
}

/**
 * This class template is a callable object (functor) which can generate random samples of
 * state-vectors taken from a gaussian belief-state.
 */
template <ContinuousBeliefState BState>
struct gaussian_sampler {
  using self = gaussian_sampler<BState>;

  using state_type =
      typename continuous_belief_state_traits<BState>::state_type;
  using covariance_type =
      typename continuous_belief_state_traits<BState>::covariance_type;

  using scalar_type =
      typename covariance_mat_traits<covariance_type>::value_type;
  using size_type = typename covariance_mat_traits<covariance_type>::size_type;
  using matrix_type =
      typename covariance_mat_traits<covariance_type>::matrix_type;

  state_type mean_state;
  mat<mat_value_type_t<matrix_type>, mat_structure::square> L;

  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  explicit gaussian_sampler(const BState& aBelief)
      : mean_state(aBelief.get_mean_state()),
        L(aBelief.get_covariance().size()) {
    using std::sqrt;
    const matrix_type& C = aBelief.get_covariance().get_matrix();
    try {
      decompose_Cholesky(C, L);
    } catch (singularity_error&) {
      mat<mat_value_type_t<matrix_type>, mat_structure::diagonal> E(
          aBelief.get_covariance().size());
      mat<mat_value_type_t<matrix_type>, mat_structure::square> U(
          aBelief.get_covariance().size());
      mat<mat_value_type_t<matrix_type>, mat_structure::square> V(
          aBelief.get_covariance().size());
      decompose_SVD(C, U, E, V);
      for (size_type i = 0; i < aBelief.get_covariance().size(); ++i) {
        E(i, i) = sqrt(E(i, i));
      }
      L = U * E;
    }
  }

  /**
   * Standard copy-constructor.
   */
  gaussian_sampler(const self& rhs) : mean_state(rhs.mean_state), L(rhs.L) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.mean_state, rhs.mean_state);
    swap(lhs.L, rhs.L);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * The call-operator which can be used to generate a random state-sample from the gaussian probability distribution.
   * \param space The state-space on which the state-vectors are distributed.
   * \return A random state-sample from the gaussian probability distribution
   */
  template <pp::Topology Space>
  state_type operator()(const Space& space) const {
    using ReaK::from_vect;
    using ReaK::to_vect;

    global_rng_type& rng = get_global_rng();
    std::normal_distribution<scalar_type> var_rnd;

    using state_difference_type = pp::topology_point_difference_type_t<Space>;
    static_assert(CovarianceMatrix<covariance_type, state_difference_type>);

    vect_n<scalar_type> z =
        to_vect<scalar_type>(space.difference(mean_state, mean_state));
    for (size_type i = 0; i < z.size(); ++i) {
      z[i] = var_rnd(rng);
    }

    return space.adjust(mean_state, from_vect<state_difference_type>(L * z));
  }
};

template <WritableVector Vector, ReadableMatrix Matrix>
Vector sample_gaussian_point(const Vector& mean, const Matrix& cov) {
  using ValueType = mat_value_type_t<Matrix>;
  using std::sqrt;

  mat<ValueType, mat_structure::square> L;
  try {
    decompose_Cholesky(cov, L);
  } catch (singularity_error&) {
    mat<ValueType, mat_structure::diagonal> E(cov.get_row_count());
    mat<ValueType, mat_structure::square> U(cov.get_row_count());
    mat<ValueType, mat_structure::square> V(cov.get_row_count());
    decompose_SVD(cov, U, E, V);
    for (int i = 0; i < cov.get_row_count(); ++i) {
      E(i, i) = sqrt(E(i, i));
    }
    L = U * E;
  }

  global_rng_type& rng = get_global_rng();
  std::normal_distribution<ValueType> var_rnd;

  Vector z = mean;
  for (int i = 0; i < z.size(); ++i) {
    z[i] = var_rnd(rng);
  }

  return mean + L * z;
}

template <WritableVector Vector, typename ValueType>
Vector sample_gaussian_point(
    Vector mean, const mat<ValueType, mat_structure::diagonal>& cov) {
  using std::sqrt;

  global_rng_type& rng = get_global_rng();
  std::normal_distribution<ValueType> var_rnd;

  for (int i = 0; i < mean.size(); ++i) {
    mean[i] += var_rnd(rng) * sqrt(cov(i, i));
  }

  return mean;
}

template <typename StateSpace, ReadableMatrix Matrix>
pp::topology_point_type_t<StateSpace> sample_gaussian_point(
    const StateSpace& space, const pp::topology_point_type_t<StateSpace>& mean,
    const Matrix& cov) {
  using DiffType = pp::topology_point_difference_type_t<StateSpace>;
  using ValueType = mat_value_type_t<Matrix>;
  using ReaK::from_vect;
  using ReaK::to_vect;
  using std::sqrt;

  mat<ValueType, mat_structure::square> L;
  try {
    decompose_Cholesky(cov, L);
  } catch (singularity_error&) {
    mat<ValueType, mat_structure::diagonal> E(cov.get_row_count());
    mat<ValueType, mat_structure::square> U(cov.get_row_count());
    mat<ValueType, mat_structure::square> V(cov.get_row_count());
    decompose_SVD(cov, U, E, V);
    for (int i = 0; i < cov.get_row_count(); ++i) {
      E(i, i) = sqrt(E(i, i));
    }
    L = U * E;
  }

  global_rng_type& rng = get_global_rng();
  std::normal_distribution<ValueType> var_rnd;

  vect_n<ValueType> z = to_vect<ValueType>(space.difference(mean, mean));
  for (int i = 0; i < z.size(); ++i) {
    z[i] = var_rnd(rng);
  }

  return space.adjust(mean, from_vect<DiffType>(L * z));
}

template <typename StateSpace, typename ValueType>
pp::topology_point_type_t<StateSpace> sample_gaussian_point(
    const StateSpace& space, const pp::topology_point_type_t<StateSpace>& mean,
    const mat<ValueType, mat_structure::diagonal>& cov) {
  using DiffType = pp::topology_point_difference_type_t<StateSpace>;
  using ReaK::from_vect;
  using ReaK::to_vect;
  using std::sqrt;

  global_rng_type& rng = get_global_rng();
  std::normal_distribution<ValueType> var_rnd;

  vect_n<ValueType> z = to_vect<ValueType>(space.difference(mean, mean));
  for (int i = 0; i < z.size(); ++i) {
    z[i] = var_rnd(rng) * sqrt(cov(i, i));
  }

  return space.adjust(mean, from_vect<DiffType>(z));
}

/**
 * This class template is used to represent a Gaussian belief-state, which is essentially a Gaussian
 * probability distribution which characterizes the estimation of a state-vector.
 *
 * \tparam StateType The state vector type which represents the mean-state of the belief, should model
 *         StateVectorConcept, by default it is the state-type associated with the covariance matrix type.
 * \tparam Covariance The covariance matrix type which represents the covariance of the state estimate.
 */
template <typename StateType, CovarianceMatrix<StateType> Covariance>
class gaussian_belief_state : public virtual shared_object {
 public:
  using self = gaussian_belief_state<StateType, Covariance>;

  using state_type = StateType;
  using state_difference_type = StateType;
  using covariance_type = Covariance;

  using scalar_type =
      typename covariance_mat_traits<covariance_type>::value_type;
  using size_type = typename covariance_mat_traits<covariance_type>::size_type;
  using matrix_type =
      typename covariance_mat_traits<covariance_type>::matrix_type;

  static constexpr belief_distribution::tag distribution =
      belief_distribution::unimodal;
  static constexpr belief_representation::tag representation =
      belief_representation::gaussian;

 private:
  state_type mean_state;
  covariance_type covar;

 public:
  /**
   * Returns the probability distribution functor associated with this belief-state's probability distribution.
   * \return The probability distribution functor associated with this belief-state's probability distribution.
   */
  auto get_pdf() const { return gaussian_pdf<self>(*this); }

  /**
   * Returns the most-likely state (i.e. the mean-state for a Gaussian distribution).
   * \return The most-likely state.
   */
  const state_type& get_most_likely_state() const { return mean_state; }

  /**
   * Returns the random sampler functor associated with this belief-state's probability distribution.
   * \return The random sampler functor associated with this belief-state's probability distribution.
   */
  auto get_random_sampler() const { return gaussian_sampler<self>(*this); }

  /**
   * Returns the mean-state state.
   * \return The mean-state.
   */
  const state_type& get_mean_state() const { return mean_state; }
  /**
   * Returns the covariance.
   * \return The covariance.
   */
  const covariance_type& get_covariance() const { return covar; }

  /**
   * Sets the mean-state.
   * \param aMeanState The new mean-state for this gaussian belief-state.
   */
  void set_mean_state(const state_type& aMeanState) { mean_state = aMeanState; }
  /**
   * Sets the covariance.
   * \param aCov The new covariance for this gaussian belief-state.
   */
  void set_covariance(const covariance_type& aCov) { covar = aCov; }

  /**
   * Returns the size of the covariance matrix of this gaussian belief-state.
   */
  size_type size() const { return covar.size(); }

  /**
   * Parametrized and default constructor.
   * \param aMeanState The mean-state of the gaussian belief-state.
   * \param aCov The covariance of the gaussian belief-state.
   */
  gaussian_belief_state(const state_type& aMeanState,  // NOLINT
                        const covariance_type& aCov)
      : mean_state(aMeanState), covar(aCov) {}
  explicit gaussian_belief_state(const state_type& aMeanState)  // NOLINT
      : gaussian_belief_state(aMeanState, covariance_type()) {}

  gaussian_belief_state()
      : gaussian_belief_state(state_type(), covariance_type()) {}

  /**
   * Standard copy-constructor.
   */
  gaussian_belief_state(const self& rhs)
      : mean_state(rhs.mean_state), covar(rhs.covar) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.mean_state, rhs.mean_state);
    swap(lhs.covar, rhs.covar);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  void save(ReaK::serialization::oarchive& aA,
            unsigned int /*Version*/) const override {
    ReaK::shared_object::save(
        aA, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(mean_state) & RK_SERIAL_SAVE_WITH_NAME(covar);
  }
  void load(ReaK::serialization::iarchive& aA,
            unsigned int /*Version*/) override {
    ReaK::shared_object::load(
        aA, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(mean_state) & RK_SERIAL_LOAD_WITH_NAME(covar);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300010, 1, "gaussian_belief_state",
                              shared_object)
};

/**
 * This function computes the symmetric KL-divergence between two belief-states.
 * \tparam StateType The state-type for the Gaussian belief-states.
 * \tparam Covariance The covariance matrix type for the Gaussian belief-states.
 * \tparam BSpace The belief-space on which the belief-states lie.
 * \param P The first Gaussian belief-state.
 * \param Q The second Gaussian belief-state.
 * \return The symmetric KL-divergence between the two Gaussian belief-states.
 */
template <typename StateType, CovarianceMatrix<StateType> Covariance,
          BeliefSpace BSpace>
auto symKL_divergence(const gaussian_belief_state<StateType, Covariance>& P,
                      const gaussian_belief_state<StateType, Covariance>& Q,
                      const BSpace& space) {
  return symKL_divergence(P.get_pdf(), Q.get_pdf(), space.get_state_topology());
}

/**
 * This function computes the symmetric KL-divergence between two belief-states.
 * \tparam StateType The state-type for the Gaussian belief-state.
 * \tparam Covariance The covariance matrix type for the Gaussian belief-state.
 * \param P The Gaussian belief-state.
 * \return The entropy of the Gaussian belief-state.
 */
template <typename StateType, CovarianceMatrix<StateType> Covariance>
auto entropy(const gaussian_belief_state<StateType, Covariance>& P) {
  return entropy(P.get_pdf());
}

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_GAUSSIAN_BELIEF_STATE_H_
