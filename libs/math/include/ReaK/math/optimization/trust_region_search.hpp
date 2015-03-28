/**
 * \file trust_region_search.hpp
 *
 * The following library is a collection of trust-region search algorithms.
 * For a quadratic model which is assumed to hold within a trust-region defined by a
 * radius on the norm of the step about the center-point, these methods find the minimum
 * of the quadratic model within that region (possibly, only an approximate solution).
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

#ifndef REAK_TRUST_REGION_SEARCH_HPP
#define REAK_TRUST_REGION_SEARCH_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>

#include "newton_search_directions.hpp"

#include <boost/bind.hpp>

namespace ReaK {


namespace optim {


namespace detail {


template < typename Vector, typename Matrix, typename T >
void compute_cauchy_point_impl( const Vector& g, const Matrix& B, Vector& p, T& norm_p, T radius, T abs_tol ) {
  T norm_g = norm_2( g );
  p = ( -radius / norm_g ) * g;
  norm_p = radius;
  T gBg = g * ( B * g );
  if( gBg > abs_tol ) {
    T tau = norm_g * norm_g * norm_g / ( radius * gBg );
    if( tau < T( 1.0 ) ) {
      p *= tau;
      norm_p *= tau;
    };
  };
};


template < typename Vector, typename Matrix, typename T, typename NewtonDirectioner >
void compute_dogleg_point_impl( const Vector& g, const Matrix& B, Vector& p, T& norm_p, T radius,
                                NewtonDirectioner get_direction, T abs_tol ) {
  using std::sqrt;
  T gg = g * g;
  T gBg = g * ( B * g );
  Vector pu = ( -gg / gBg ) * g;
  T norm_sqr_pu = ( gg * gg * gg ) / ( gBg * gBg );
  if( norm_sqr_pu > radius * radius ) {
    p = pu;
    p *= radius / sqrt( norm_sqr_pu );
    norm_p = radius;
    return;
  };

  Vector pb = g;
  try {
    get_direction( B, g, pb, abs_tol );
  } catch( singularity_error& ) {
    p = pu;
    norm_p = sqrt( norm_sqr_pu );
    return;
  };
  if( pb * ( B * pb ) < 0.0 ) {
    p = pu;
    norm_p = norm_2( p );
    return;
  };
  norm_p = norm_2( pb );
  if( norm_p < radius ) {
    p = pb;
    return;
  };

  Vector dp = pb;
  dp -= pu;
  T norm_sqr_dp = dp * dp;
  T cross_term = T( 2.0 ) * ( pu * dp );

  T temp = sqrt( cross_term * cross_term - T( 4.0 ) * norm_sqr_dp * ( norm_sqr_pu - radius * radius ) );
  T alpha1 = ( -cross_term + temp ) / ( T( 2.0 ) * norm_sqr_dp );
  T alpha2 = ( -cross_term - temp ) / ( T( 2.0 ) * norm_sqr_dp );
  if( ( alpha1 > T( 0.0 ) ) && ( alpha1 <= T( 1.0 ) ) ) {
    norm_p = radius;
    p = dp;
    p *= alpha1;
    p += pu;
  } else {
    norm_p = radius;
    p = dp;
    p *= alpha2;
    p += pu;
  };
};


template < typename Vector1, typename Matrix, typename Vector2, typename T, typename NewtonDirectioner >
void compute_right_pinv_dogleg_impl( const Vector1& c, const Matrix& A, Vector2& p, T& norm_p, T radius,
                                     NewtonDirectioner get_direction, T abs_tol ) {
  using std::sqrt;

  Vector2 Atc = c * A;
  Vector1 AAtc = A * Atc;
  T cAAc = Atc * Atc;
  T cAAAAc = AAtc * AAtc;
  Vector2 pu = ( -cAAc / cAAAAc ) * Atc;
  T norm_sqr_pu = ( cAAc * cAAc * cAAc ) / ( cAAAAc * cAAAAc );
  if( norm_sqr_pu > radius * radius ) {
    p = pu;
    p *= radius / sqrt( norm_sqr_pu );
    norm_p = radius;
    return;
  };

  Vector1 cb = c;
  mat< T, mat_structure::symmetric > AAt( A * transpose_view( A ) );
  try {
    get_direction( AAt, c, cb, abs_tol );
  } catch( singularity_error& ) {
    p = pu;
    norm_p = sqrt( norm_sqr_pu );
    return;
  };
  if( cb * ( AAt * cb ) < 0.0 ) {
    p = pu;
    norm_p = norm_2( p );
    return;
  };
  Vector2 pb = cb * A;
  norm_p = norm_2( pb );
  if( norm_p < radius ) {
    p = pb;
    return;
  };

  Vector2 dp = pb;
  dp -= pu;
  T norm_sqr_dp = dp * dp;
  T cross_term = T( 2.0 ) * ( pu * dp );

  T temp = sqrt( cross_term * cross_term - T( 4.0 ) * norm_sqr_dp * ( norm_sqr_pu - radius * radius ) );
  T alpha1 = ( -cross_term + temp ) / ( T( 2.0 ) * norm_sqr_dp );
  T alpha2 = ( -cross_term - temp ) / ( T( 2.0 ) * norm_sqr_dp );
  if( ( alpha1 > T( 0.0 ) ) && ( alpha1 <= T( 1.0 ) ) ) {
    norm_p = radius;
    p = dp;
    p *= alpha1;
    p += pu;
  } else {
    norm_p = radius;
    p = dp;
    p *= alpha2;
    p += pu;
  };
};
};


/**
 * This function computes the Cauchy point which is a steepest descent point that minimizes a
 * function within a trust-region of a given radius.
 * TEST PASSED
 * \tparam Vector A writable vector type.
 * \tparam Matrix A readable matrix type.
 * \param g The gradient vector of the function at the center of the trust-region.
 * \param B The Hessian (or approximate Hessian) of the function at the center of the trust-region.
 * \param p The resulting cauchy-point.
 * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
 * \param radius The radius of the trust-region.
 * \param abs_tol The tolerance at which to consider values to be zero.
 */
template < typename Vector, typename Matrix >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector >, is_readable_matrix< Matrix > >, void >::type
  compute_cauchy_point( const Vector& g, const Matrix& B, Vector& p, typename vect_traits< Vector >::value_type& norm_p,
                        typename vect_traits< Vector >::value_type radius,
                        typename vect_traits< Vector >::value_type abs_tol
                        = typename vect_traits< Vector >::value_type( 1e-6 ) ) {
  detail::compute_cauchy_point_impl( g, B, p, norm_p, radius, abs_tol );
};

/**
 * This functor class computes the Cauchy point which is a steepest descent point that minimizes a
 * function within a trust-region of a given radius.
 * TEST PASSED
 */
struct trust_region_solver_cauchy {
  /**
   * This function computes the Cauchy point which is a steepest descent point that minimizes a
   * function within a trust-region of a given radius.
   * \tparam Vector A writable vector type.
   * \tparam Matrix A readable matrix type.
   * \param g The gradient vector of the function at the center of the trust-region.
   * \param B The Hessian (or approximate Hessian) of the function at the center of the trust-region.
   * \param p The resulting cauchy-point.
   * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
   * \param radius The radius of the trust-region.
   * \param abs_tol The tolerance at which to consider values to be zero.
   */
  template < typename Vector, typename Matrix >
  void operator()( const Vector& g, const Matrix& B, Vector& p, typename vect_traits< Vector >::value_type& norm_p,
                   typename vect_traits< Vector >::value_type radius,
                   typename vect_traits< Vector >::value_type abs_tol
                   = typename vect_traits< Vector >::value_type( 1e-6 ) ) const {
    compute_cauchy_point( g, B, p, norm_p, radius, abs_tol );
  };
};


/**
 * This function computes the dogleg point which follows the steepest descent and then the Newton
 * direction in order to minimizes a quadratic function within a trust-region of a given radius.
 * TEST PASSED
 * \tparam Vector A writable vector type.
 * \tparam Matrix A readable matrix type.
 * \param g The gradient vector of the function at the center of the trust-region.
 * \param B The Hessian (or approximate Hessian) of the function at the center of the trust-region.
 * \param p The resulting cauchy-point.
 * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
 * \param radius The radius of the trust-region.
 * \param abs_tol The tolerance at which to consider values to be zero.
 */
template < typename Vector, typename Matrix >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector >, is_readable_matrix< Matrix > >, void >::type
  compute_dogleg_point( const Vector& g, const Matrix& B, Vector& p, typename vect_traits< Vector >::value_type& norm_p,
                        typename vect_traits< Vector >::value_type radius,
                        typename vect_traits< Vector >::value_type abs_tol
                        = typename vect_traits< Vector >::value_type( 1e-6 ) ) {
  detail::compute_dogleg_point_impl( g, B, p, norm_p, radius, newton_direction< Matrix, Vector >, abs_tol );
};

/**
 * This functor class computes the dogleg point which follows the steepest descent and then the Newton
 * direction in order to minimizes a quadratic function within a trust-region of a given radius.
 * TEST PASSED
 */
struct trust_region_solver_dogleg {
  /**
   * This function computes the dogleg point which follows the steepest descent and then the Newton
   * direction in order to minimizes a quadratic function within a trust-region of a given radius.
   * \tparam Vector A writable vector type.
   * \tparam Matrix A readable matrix type.
   * \param g The gradient vector of the function at the center of the trust-region.
   * \param B The Hessian (or approximate Hessian) of the function at the center of the trust-region.
   * \param p The resulting cauchy-point.
   * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
   * \param radius The radius of the trust-region.
   * \param abs_tol The tolerance at which to consider values to be zero.
   */
  template < typename Vector, typename Matrix >
  void operator()( const Vector& g, const Matrix& B, Vector& p, typename vect_traits< Vector >::value_type& norm_p,
                   typename vect_traits< Vector >::value_type radius,
                   typename vect_traits< Vector >::value_type abs_tol
                   = typename vect_traits< Vector >::value_type( 1e-6 ) ) const {
    compute_dogleg_point( g, B, p, norm_p, radius, abs_tol );
  };
};

/**
 * This functor class computes the dogleg point which follows the steepest descent and then the
 * regularized Newton direction in order to minimizes a quadratic function within a
 * trust-region of a given radius.
 * \test Must create a unit-test for this.
 */
template < typename T >
struct trust_region_solver_dogleg_reg {
  regularized_newton_directioner< T > get_reg_direction;

  /**
   * Parametrized Constructor.
   * \param aTau The initial relative damping factor for the damping the Hessian matrix (the actual damping factor is
   * relative to the trace of the Hessian).
   */
  trust_region_solver_dogleg_reg( T aTau = T( 1e-3 ) ) : get_reg_direction( aTau ){};
  /**
   * This function computes the dogleg point which follows the steepest descent and then the
   * regularized Newton direction in order to minimizes a quadratic function within a
   * trust-region of a given radius.
   * \tparam Vector A writable vector type.
   * \tparam Matrix A readable matrix type.
   * \param g The gradient vector of the function at the center of the trust-region.
   * \param B The Hessian (or approximate Hessian) of the function at the center of the trust-region.
   * \param p The resulting cauchy-point.
   * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
   * \param radius The radius of the trust-region.
   * \param abs_tol The tolerance at which to consider values to be zero.
   */
  template < typename Vector, typename Matrix >
  void operator()( const Vector& g, const Matrix& B, Vector& p, typename vect_traits< Vector >::value_type& norm_p,
                   typename vect_traits< Vector >::value_type radius,
                   typename vect_traits< Vector >::value_type abs_tol
                   = typename vect_traits< Vector >::value_type( 1e-6 ) ) const {
    detail::compute_dogleg_point_impl( g, B, p, norm_p, radius,
                                       boost::bind(
#ifdef _MSC_VER
                                         &(regularized_newton_directioner< T >::operator()< Matrix, Vector >),
#else
                                         &regularized_newton_directioner< T >::template operator()< Matrix, Vector >,
#endif
                                         &get_reg_direction, _1, _2, _3, _4 ),
                                       abs_tol );
  };
};


/**
 * This function computes the dogleg point which follows the steepest descent and then the minimum
 * norm pseudo-inverse direction in order to minimizes a quadratic function within a trust-region
 * of a given radius.
 * \test Must create a unit-test for this.
 * \tparam Vector1 A writable vector type.
 * \tparam Matrix A readable matrix type.
 * \tparam Vector2 A writable vector type.
 * \param g The gradient vector of the function at the center of the trust-region.
 * \param B The Jacobian matrix of the function at the center of the trust-region.
 * \param p The resulting point.
 * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
 * \param radius The radius of the trust-region.
 * \param abs_tol The tolerance at which to consider values to be zero.
 */
template < typename Vector1, typename Matrix, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector1 >, is_readable_matrix< Matrix >,
                                             is_writable_vector< Vector2 > >,
                           void >::type
  compute_right_pinv_dogleg( const Vector1& g, const Matrix& B, Vector2& p,
                             typename vect_traits< Vector2 >::value_type& norm_p,
                             typename vect_traits< Vector2 >::value_type radius,
                             typename vect_traits< Vector2 >::value_type abs_tol
                             = typename vect_traits< Vector2 >::value_type( 1e-6 ) ) {
  detail::compute_right_pinv_dogleg_impl(
    g, B, p, norm_p, radius,
    newton_direction< mat< typename vect_traits< Vector2 >::value_type, mat_structure::symmetric >, Vector1 >,
    abs_tol );
};

/**
 * This functor class computes the dogleg point which follows the steepest descent and then the minimum
 * norm pseudo-inverse direction in order to minimizes a quadratic function within a trust-region
 * of a given radius.
 * \test Must create a unit-test for this.
 */
struct tr_solver_right_pinv_dogleg {
  /**
   * This function computes the dogleg point which follows the steepest descent and then the minimum
   * norm pseudo-inverse direction in order to minimizes a quadratic function within a trust-region
   * of a given radius.
   * \tparam Vector1 A writable vector type.
   * \tparam Matrix A readable matrix type.
   * \tparam Vector2 A writable vector type.
   * \param g The gradient vector of the function at the center of the trust-region.
   * \param B The Jacobian matrix of the function at the center of the trust-region.
   * \param p The resulting point.
   * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
   * \param radius The radius of the trust-region.
   * \param abs_tol The tolerance at which to consider values to be zero.
   */
  template < typename Vector1, typename Matrix, typename Vector2 >
  void operator()( const Vector1& g, const Matrix& B, Vector2& p, typename vect_traits< Vector2 >::value_type& norm_p,
                   typename vect_traits< Vector2 >::value_type radius,
                   typename vect_traits< Vector2 >::value_type abs_tol
                   = typename vect_traits< Vector2 >::value_type( 1e-6 ) ) const {
    compute_right_pinv_dogleg( g, B, p, norm_p, radius, abs_tol );
  };
};

/**
 * This functor class computes the dogleg point which follows the steepest descent and then the
 * regularized minimum norm pseudo-inverse direction in order to minimizes a quadratic function
 * within a trust-region of a given radius.
 * \test Must create a unit-test for this.
 */
template < typename T >
struct tr_solver_right_pinv_dogleg_reg {
  regularized_newton_directioner< T > get_reg_direction;

  /**
   * Parametrized Constructor.
   * \param aTau The initial relative damping factor for the damping the Hessian matrix (the actual damping factor is
   * relative to the trace of the Hessian).
   */
  tr_solver_right_pinv_dogleg_reg( T aTau = T( 1e-3 ) ) : get_reg_direction( aTau ){};
  /**
   * This function computes the dogleg point which follows the steepest descent and then the
   * regularized minimum norm pseudo-inverse direction in order to minimizes a quadratic function
   * within a trust-region of a given radius.
   * \tparam Vector1 A writable vector type.
   * \tparam Matrix A readable matrix type.
   * \tparam Vector2 A writable vector type.
   * \param g The gradient vector of the function at the center of the trust-region.
   * \param B The Jacobian matrix of the function at the center of the trust-region.
   * \param p The resulting point.
   * \param norm_p The resulting norm of the cauchy-point (less-than or equal to the trust-region radius).
   * \param radius The radius of the trust-region.
   * \param abs_tol The tolerance at which to consider values to be zero.
   */
  template < typename Vector1, typename Matrix, typename Vector2 >
  void operator()( const Vector1& g, const Matrix& B, Vector2& p, typename vect_traits< Vector2 >::value_type& norm_p,
                   typename vect_traits< Vector2 >::value_type radius,
                   typename vect_traits< Vector2 >::value_type abs_tol
                   = typename vect_traits< Vector2 >::value_type( 1e-6 ) ) const {
    detail::compute_right_pinv_dogleg_impl(
      g, B, p, norm_p, radius,
      boost::bind(
#ifdef _MSC_VER
        &(regularized_newton_directioner< T >::
            operator()< mat< typename vect_traits< Vector2 >::value_type, mat_structure::symmetric >, Vector1 >),
#else
        &regularized_newton_directioner< T >::template
        operator()< mat< typename vect_traits< Vector2 >::value_type, mat_structure::symmetric >, Vector1 >,
#endif
        &get_reg_direction, _1, _2, _3, _4 ),
      abs_tol );
  };
};
};
};


#endif
