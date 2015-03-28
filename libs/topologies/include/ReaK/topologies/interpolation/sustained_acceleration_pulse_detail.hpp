/**
 * \file sustained_acceleration_pulse_detail.hpp
 *
 * This library contains the implementation details of the rate-limited sustained acceleration
 * pulse (SAP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_SUSTAINED_ACCELERATION_PULSE_DETAIL_HPP
#define REAK_SUSTAINED_ACCELERATION_PULSE_DETAIL_HPP

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/topologies/spaces/bounded_space_concept.hpp>
#include <ReaK/topologies/spaces/tangent_bundle_concept.hpp>

#include "sustained_velocity_pulse_detail.hpp"

#include <boost/bind.hpp>
#include <cmath>

namespace ReaK {

namespace pp {


namespace detail {

double sap_solve_for_min_dt_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                                  double dt_amax );

bool sap_min_dt_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol, double& result,
                           double dt_amax );

double sap_solve_for_no_slack_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                                    double delta_time, double dt_amax );

bool sap_no_slack_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol, double& slack,
                             double delta_time, double dt_amax );


template < typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 2 > >, void >::type
  sap_interpolate_HOT_impl( PointType& result, const PointDiff1& delta_second_order, const DiffSpace& space,
                            const TimeSpace& t_space, double t_factor ){/* nothing to do. */
  };

template < typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::equal_to< Idx, boost::mpl::size_t< 2 > >, void >::type
  sap_interpolate_HOT_impl( PointType& result, const PointDiff1& delta_second_order, const DiffSpace& space,
                            const TimeSpace& t_space, double t_factor ) {
  get< 2 >( result ) = lift_to_space< 2 >( delta_second_order, t_factor, space, t_space );
};


template < typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_constant_jerk_motion_HOT_impl( PointType& result, const PointDiff2& descended_jerk, const DiffSpace& space,
                                     const TimeSpace& t_space ){/* Nothing to do. */
  };

template < typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::equal_to< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_constant_jerk_motion_HOT_impl( PointType& result, const PointDiff2& descended_jerk, const DiffSpace& space,
                                     const TimeSpace& t_space ) {
  get< 3 >( result ) = lift_to_space< 3 >( descended_jerk, 1.0, space, t_space );
};

template < typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_constant_jerk_motion_HOT_impl( PointType& result, const PointDiff2& descended_jerk, const DiffSpace& space,
                                     const TimeSpace& t_space ) {
  sap_constant_jerk_motion_HOT_impl< typename boost::mpl::prior< Idx >::type >( result, descended_jerk, space,
                                                                                t_space );

  get< Idx::type::value >( result ) = get_space< Idx::type::value >( space, t_space ).origin();
};

template < typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace >
inline void sap_constant_jerk_motion_impl( PointType& result, const PointDiff2& descended_jerk, const DiffSpace& space,
                                           const TimeSpace& t_space, double dt ) {

  get< 0 >( result ) = get_space< 0 >( space, t_space ).adjust(
    get< 0 >( result ),
    descend_to_space< 0 >(
      get_space< 1 >( space, t_space )
        .adjust( get< 1 >( result ), descend_to_space< 1 >( get_space< 2 >( space, t_space ).adjust(
                                                              get< 2 >( result ), ( 1.0 / 6.0 ) * dt * descended_jerk ),
                                                            dt, space, t_space ) ),
      dt, space, t_space ) );


  get< 1 >( result ) = get_space< 1 >( space, t_space ).adjust(
    get< 1 >( result ),
    descend_to_space< 1 >( get_space< 2 >( space, t_space ).adjust( get< 2 >( result ), 0.5 * dt * descended_jerk ), dt,
                           space, t_space ) );

  get< 2 >( result ) = get_space< 2 >( space, t_space ).adjust( get< 2 >( result ), dt * descended_jerk );

  if( Idx::type::value > 2 )
    sap_constant_jerk_motion_HOT_impl< Idx >( result, descended_jerk, space, t_space );

  return;
};


template < typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace,
           typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 4 > >, void >::type
  sap_interpolate_impl( PointType& result, const PointType& start_point, const PointType& end_point,
                        const PointDiff0& delta_first_order, const PointType1& peak_velocity, const DiffSpace& space,
                        const TimeSpace& t_space, double dt, double dt_total ) {
  using std::sqrt;

  result = start_point;

  double dt_amax = get_space< 2 >( space, t_space ).get_radius();

  double dt_vp_1st = get( distance_metric, get_space< 1 >( space, t_space ) )( get< 1 >( start_point ), peak_velocity,
                                                                               get_space< 1 >( space, t_space ) );
  // we know that dt_vp_2nd = dt_vp_1st + dt_amax
  double dt_vp = dt_vp_1st - dt_amax;
  double dt_ap = dt_amax;
  if( dt_vp < 0.0 ) {
    // means that we don't have time to reach the maximum acceleration:
    dt_vp = 0.0;
    dt_ap = sqrt( dt_amax * dt_vp_1st );
  };

  double dt_vp2_1st = get( distance_metric, get_space< 1 >( space, t_space ) )( peak_velocity, get< 1 >( end_point ),
                                                                                get_space< 1 >( space, t_space ) );
  // we know that dt_vp_2nd = dt_vp_1st + dt_amax
  double dt_vp2 = dt_vp2_1st - dt_amax;
  double dt_ap2 = dt_amax;
  if( dt_vp2 < 0.0 ) {
    // means that we don't have time to reach the maximum acceleration:
    dt_vp2 = 0.0;
    dt_ap2 = sqrt( dt_amax * dt_vp2_1st );
  };
  dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp + 2.0 * dt_ap;

  if( dt_vp_1st > std::numeric_limits< double >::epsilon() ) {
    // Phase 1: in the jerk-up phase of velocity ramp-up.
    if( dt < dt_ap )
      dt_ap = dt;

    sap_constant_jerk_motion_impl< Idx >(
      result,
      get_space< 2 >( space, t_space ).difference(
        lift_to_space< 2 >( get_space< 1 >( space, t_space ).difference( peak_velocity, get< 1 >( start_point ) ),
                            dt_vp_1st * dt_amax, space, t_space ),
        get_space< 2 >( space, t_space ).origin() ),
      space, t_space, dt_ap );
    dt -= dt_ap;
    if( dt <= std::numeric_limits< double >::epsilon() )
      return;

    // Phase 2: in the constant accel phase of velocity ramp-up.
    if( dt_vp > std::numeric_limits< double >::epsilon() ) {
      if( dt < dt_vp )
        dt_vp = dt;

      svp_constant_accel_motion_impl< Idx >( result, descend_to_space< 1 >( get< 2 >( result ), 1.0, space, t_space ),
                                             space, t_space, dt_vp );
      dt -= dt_vp;
      if( dt <= std::numeric_limits< double >::epsilon() )
        return;
    };

    // Phase 3: in the jerk-down phase of velocity ramp-up.
    if( dt < dt_ap )
      dt_ap = dt;

    sap_constant_jerk_motion_impl< Idx >(
      result,
      get_space< 2 >( space, t_space ).difference(
        get_space< 2 >( space, t_space ).origin(),
        lift_to_space< 2 >( get_space< 1 >( space, t_space ).difference( peak_velocity, get< 1 >( start_point ) ),
                            dt_vp_1st * dt_amax, space, t_space ) ),
      space, t_space, dt_ap );
    dt -= dt_ap;
    if( dt <= std::numeric_limits< double >::epsilon() )
      return;
  };


  // Phase 4: in the cruise phase.
  if( dt < dt_total )
    dt_total = dt;

  svp_constant_vel_motion_impl< Idx >( result, descend_to_space< 0 >( get< 1 >( result ), 1.0, space, t_space ), space,
                                       t_space, dt_total );
  dt -= dt_total;
  if( dt <= std::numeric_limits< double >::epsilon() )
    return;

  if( dt_vp2_1st > std::numeric_limits< double >::epsilon() ) {
    // Phase 5: in the jerk-up phase of velocity ramp-down.
    if( dt < dt_ap2 )
      dt_ap2 = dt;

    sap_constant_jerk_motion_impl< Idx >(
      result, get_space< 2 >( space, t_space ).difference(
                lift_to_space< 2 >( get_space< 1 >( space, t_space ).difference( get< 1 >( end_point ), peak_velocity ),
                                    dt_vp2_1st * dt_amax, space, t_space ),
                get_space< 2 >( space, t_space ).origin() ),
      space, t_space, dt_ap2 );
    dt -= dt_ap2;
    if( dt <= std::numeric_limits< double >::epsilon() )
      return;

    // Phase 6: in the constant accel phase of velocity ramp-down.
    if( dt_vp2 > std::numeric_limits< double >::epsilon() ) {
      if( dt < dt_vp2 )
        dt_vp2 = dt;

      svp_constant_accel_motion_impl< Idx >( result, descend_to_space< 1 >( get< 2 >( result ), 1.0, space, t_space ),
                                             space, t_space, dt_vp2 );
      dt -= dt_vp2;
      if( dt <= std::numeric_limits< double >::epsilon() )
        return;
    };

    // Phase 7: in the jerk-down phase of velocity ramp-down.
    if( dt < dt_ap2 )
      dt_ap2 = dt;

    sap_constant_jerk_motion_impl< Idx >(
      result, get_space< 2 >( space, t_space ).difference(
                get_space< 2 >( space, t_space ).origin(),
                lift_to_space< 2 >( get_space< 1 >( space, t_space ).difference( get< 1 >( end_point ), peak_velocity ),
                                    dt_vp2_1st * dt_amax, space, t_space ) ),
      space, t_space, dt_ap2 );
    dt -= dt_ap2;
  };
};

template < typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace,
           typename TimeSpace >
inline typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_interpolate_impl( PointType& result, const PointType& start_point, const PointType& end_point,
                        const PointDiff0& delta_first_order, const PointType1& peak_velocity, const DiffSpace& space,
                        const TimeSpace& t_space, double dt, double dt_total ) {
  sap_interpolate_impl< typename boost::mpl::prior< Idx >::type >( result, start_point, end_point, delta_first_order,
                                                                   peak_velocity, space, t_space, dt, dt_total );

  get< Idx::type::value >( result ) = get_space< Idx::type::value >( space, t_space ).origin();
};

struct sap_update_delta_first_order {
  template < typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace >
  void operator()( const PointType& start_point, const PointType& end_point, PointDiff0& delta_first_order,
                   PointType1& peak_velocity, double& norm_delta, const double& beta, const DiffSpace& space,
                   const TimeSpace& t_space, double num_tol = 1E-6 ) {
    using std::fabs;
    using std::sqrt;

    double dt_amax = get_space< 2 >( space, t_space ).get_radius();

    double dt_vp1_1st = get( distance_metric, get_space< 1 >( space, t_space ) )(
      get< 1 >( start_point ), peak_velocity, get_space< 1 >( space, t_space ) );
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp1 = dt_vp1_1st - dt_amax;
    double dt_ap1 = dt_amax;
    if( dt_vp1 < 0.0 ) {
      // means that we don't have time to reach the maximum acceleration:
      dt_vp1 = 0.0;
      dt_ap1 = sqrt( dt_amax * dt_vp1_1st );
    };

    double dt_vp2_1st = get( distance_metric, get_space< 1 >( space, t_space ) )( peak_velocity, get< 1 >( end_point ),
                                                                                  get_space< 1 >( space, t_space ) );
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - dt_amax;
    double dt_ap2 = dt_amax;
    if( dt_vp2 < 0.0 ) {
      // means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt( dt_amax * dt_vp2_1st );
    };


    PointDiff0 start_to_peak
      = get_space< 0 >( space, t_space ).difference( get< 0 >( start_point ), get< 0 >( start_point ) );

    if( dt_vp1_1st > num_tol * get_space< 1 >( space, t_space ).get_radius() ) {
      start_to_peak
        = descend_to_space< 0 >( get< 1 >( start_point ), dt_vp1, space, t_space )
          + descend_to_space< 0 >(
              get_space< 1 >( space, t_space )
                .adjust( get< 1 >( start_point ),
                         ( ( 0.75 * dt_ap1 * ( dt_ap1 + dt_vp1 ) + 0.25 * dt_vp1 * dt_vp1 ) / ( dt_amax * dt_vp1_1st ) )
                         * get_space< 1 >( space, t_space ).difference( peak_velocity, get< 1 >( start_point ) ) ),
              2.0 * dt_ap1, space, t_space );
    };

    PointDiff0 peak_to_end
      = get_space< 0 >( space, t_space ).difference( get< 0 >( end_point ), get< 0 >( end_point ) );

    if( dt_vp2_1st > num_tol * get_space< 1 >( space, t_space ).get_radius() ) {
      peak_to_end
        = descend_to_space< 0 >( peak_velocity, dt_vp2, space, t_space )
          + descend_to_space< 0 >(
              get_space< 1 >( space, t_space )
                .adjust( peak_velocity,
                         ( ( 0.75 * dt_ap2 * ( dt_ap2 + dt_vp2 ) + 0.25 * dt_vp2 * dt_vp2 ) / ( dt_amax * dt_vp2_1st ) )
                         * get_space< 1 >( space, t_space ).difference( get< 1 >( end_point ), peak_velocity ) ),
              2.0 * dt_ap2, space, t_space );
    };

    delta_first_order = get_space< 0 >( space, t_space ).difference(
      get_space< 0 >( space, t_space ).adjust( get< 0 >( end_point ), -peak_to_end ),
      get_space< 0 >( space, t_space ).adjust( get< 0 >( start_point ), start_to_peak ) );

    norm_delta
      = get( distance_metric, get_space< 0 >( space, t_space ) )( delta_first_order, get_space< 0 >( space, t_space ) );
    PointDiff0 descended_peak_velocity = descend_to_space< 0 >( peak_velocity, 1.0, space, t_space );
    double normA = get( distance_metric, get_space< 0 >( space, t_space ) )( descended_peak_velocity,
                                                                             get_space< 0 >( space, t_space ) );
    double normC = get( distance_metric, get_space< 0 >( space, t_space ) )(
      descended_peak_velocity - delta_first_order, get_space< 0 >( space, t_space ) );
    if( normC * normC > normA * normA + norm_delta * norm_delta )
      norm_delta = -norm_delta;
    if( fabs( norm_delta ) > num_tol )
      peak_velocity = lift_to_space< 1 >( delta_first_order * beta, norm_delta, space, t_space );
  };
};


template < typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace >
double sap_compute_min_delta_time( const PointType& start_point, const PointType& end_point,
                                   PointDiff0& delta_first_order, PointType1& peak_velocity, double& norm_delta,
                                   double& beta, const DiffSpace& space, const TimeSpace& t_space,
                                   double num_tol = 1E-6, unsigned int max_iter = 20 ) {

  double result = 0.0;
  // RK_NOTICE(1,"*********  Starting Min-dt Iterations    *****************");

  svp_peak_velocity_iteration( start_point, end_point, delta_first_order, peak_velocity, norm_delta, beta, space,
                               t_space, sap_update_delta_first_order(),
                               boost::bind( sap_solve_for_min_dt_beta, _1, _2, _3, _4,
                                            boost::cref( get_space< 2 >( space, t_space ).get_radius() ) ),
                               boost::bind( sap_min_dt_predicate, _1, _2, _3, _4, boost::ref( result ),
                                            boost::cref( get_space< 2 >( space, t_space ).get_radius() ) ),
                               num_tol, max_iter );

  return result;
};

template < typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace >
double sap_compute_peak_velocity( const PointType& start_point, const PointType& end_point,
                                  PointDiff0& delta_first_order, PointType1& peak_velocity, double& norm_delta,
                                  double& beta, double delta_time, const DiffSpace& space, const TimeSpace& t_space,
                                  double num_tol = 1E-6, unsigned int max_iter = 20 ) {

  double slack = 0.0;
  // RK_NOTICE(1,"*********  Starting No-slack Iterations  *****************");

  svp_peak_velocity_iteration(
    start_point, end_point, delta_first_order, peak_velocity, norm_delta, beta, space, t_space,
    sap_update_delta_first_order(), boost::bind( sap_solve_for_no_slack_beta, _1, _2, _3, _4, boost::cref( delta_time ),
                                                 boost::cref( get_space< 2 >( space, t_space ).get_radius() ) ),
    boost::bind( sap_no_slack_predicate, _1, _2, _3, _4, boost::ref( slack ), boost::cref( delta_time ),
                 boost::cref( get_space< 2 >( space, t_space ).get_radius() ) ),
    num_tol, max_iter );

  return slack;
};

template < typename PointType, typename DiffSpace, typename TimeSpace >
double sap_compute_interpolation_data_impl(
  const PointType& start_point, const PointType& end_point,
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 0 >::type >::point_difference_type&
    delta_first_order,
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type&
    peak_velocity,
  const DiffSpace& space, const TimeSpace& t_space, double delta_time = 0.0,
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type
  * best_peak_velocity = NULL,
  double num_tol = 1e-6, unsigned int max_iter = 20 ) {

  delta_first_order = get_space< 0 >( space, t_space ).difference( get< 0 >( end_point ), get< 0 >( start_point ) );
  double norm_delta
    = get( distance_metric, get_space< 0 >( space, t_space ) )( delta_first_order, get_space< 0 >( space, t_space ) );
  double beta = 0.0;
  peak_velocity = get_space< 1 >( space, t_space ).origin();

  double min_delta_time = sap_compute_min_delta_time( start_point, end_point, delta_first_order, peak_velocity,
                                                      norm_delta, beta, space, t_space, num_tol, max_iter );

  if( best_peak_velocity )
    *best_peak_velocity = peak_velocity;

  if( min_delta_time > delta_time )
    return min_delta_time;

  beta = beta * min_delta_time / delta_time;
  peak_velocity = get_space< 1 >( space, t_space ).adjust(
    get_space< 1 >( space, t_space ).origin(),
    ( min_delta_time / delta_time )
    * get_space< 1 >( space, t_space ).difference( peak_velocity, get_space< 1 >( space, t_space ).origin() ) );

  sap_compute_peak_velocity( start_point, end_point, delta_first_order, peak_velocity, norm_delta, beta, delta_time,
                             space, t_space, num_tol, max_iter );

  return min_delta_time;
};
};
};
};

#endif
