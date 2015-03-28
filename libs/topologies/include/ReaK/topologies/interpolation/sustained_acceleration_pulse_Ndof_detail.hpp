/**
 * \file sustained_acceleration_pulse_Ndof_detail.hpp
 *
 * This library contains the implementation details of the rate-limited sustained acceleration
 * pulse (SAP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_SUSTAINED_ACCELERATION_PULSE_NDOF_DETAIL_HPP
#define REAK_SUSTAINED_ACCELERATION_PULSE_NDOF_DETAIL_HPP

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/optimization/optim_exceptions.hpp>

#include <ReaK/topologies/spaces/tangent_bundle_concept.hpp>

#include "sustained_velocity_pulse_Ndof_detail.hpp"

#include <cmath>

namespace ReaK {

namespace pp {


namespace detail {

void sap_Ndof_compute_interpolated_values( double start_position, double end_position, double start_velocity,
                                           double end_velocity, double peak_velocity, double max_velocity,
                                           double max_acceleration, double dt, double dt_total, double& result_pos,
                                           double& result_vel, double& result_acc, double& result_jerk );

double sap_Ndof_compute_min_delta_time( double start_position, double end_position, double start_velocity,
                                        double end_velocity, double& peak_velocity, double max_velocity,
                                        double max_acceleration );

void sap_Ndof_compute_peak_velocity( double start_position, double end_position, double start_velocity,
                                     double end_velocity, double& peak_velocity, double max_velocity,
                                     double max_acceleration, double delta_time );


template < typename Idx, typename PointType, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_Ndof_constant_jerk_motion_HOT_impl( PointType&, double, std::size_t, const DiffSpace&,
                                          const TimeSpace& ){/* Nothing to do. */
  };

template < typename Idx, typename PointType, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::equal_to< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_Ndof_constant_jerk_motion_HOT_impl( PointType& result, double descended_jerk, std::size_t i,
                                          const DiffSpace& space, const TimeSpace& t_space ) {
  typedef typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 2 >::type >::point_type
    PointType2;
  const PointType2& max_jerk = get_space< 3 >( space, t_space ).get_upper_corner();
  get< 3 >( result )[i] = descended_jerk * max_jerk[i];
};

template < typename Idx, typename PointType, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_Ndof_constant_jerk_motion_HOT_impl( PointType& result, double descended_jerk, std::size_t i,
                                          const DiffSpace& space, const TimeSpace& t_space ) {
  sap_Ndof_constant_jerk_motion_HOT_impl< typename boost::mpl::prior< Idx >::type >( result, descended_jerk, i, space,
                                                                                     t_space );

  get< Idx::type::value >( result )[i] = 0.0;
};


template < typename Idx, typename PointType, typename PointType1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 4 > >, void >::type
  sap_Ndof_interpolate_impl( PointType& result, const PointType& start_point, const PointType& end_point,
                             const PointType1& peak_velocity, const DiffSpace& space, const TimeSpace& t_space,
                             double dt, double dt_total ) {
  using std::sqrt;
  using std::fabs;

  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type max_velocity
    = get_space< 1 >( space, t_space ).get_upper_corner();
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 2 >::type >::point_type
    max_acceleration = get_space< 2 >( space, t_space ).get_upper_corner();

  for( std::size_t i = 0; i < peak_velocity.size(); ++i ) {

    double result_pos, result_vel, result_acc, result_desc_jerk;

    sap_Ndof_compute_interpolated_values( get< 0 >( start_point )[i], get< 0 >( end_point )[i],
                                          get< 1 >( start_point )[i], get< 1 >( end_point )[i], peak_velocity[i],
                                          max_velocity[i], max_acceleration[i], dt, dt_total, result_pos, result_vel,
                                          result_acc, result_desc_jerk );

    get< 0 >( result )[i] = result_pos;
    get< 1 >( result )[i] = result_vel;
    get< 2 >( result )[i] = result_acc;
    sap_Ndof_constant_jerk_motion_HOT_impl< Idx >( result, result_desc_jerk, i, space, t_space );
  };
};

template < typename Idx, typename PointType, typename PointType1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 3 > >, void >::type
  sap_Ndof_interpolate_impl( PointType& result, const PointType& start_point, const PointType& end_point,
                             const PointType1& peak_velocity, const DiffSpace& space, const TimeSpace& t_space,
                             double dt, double dt_total ) {
  sap_Ndof_interpolate_impl< typename boost::mpl::prior< Idx >::type >( result, start_point, end_point, peak_velocity,
                                                                        space, t_space, dt, dt_total );

  get< Idx::type::value >( result ) = get_space< Idx::type::value >( space, t_space ).origin();
};


template < typename PointType, typename DiffSpace, typename TimeSpace >
double sap_compute_Ndof_interpolation_data_impl(
  const PointType& start_point, const PointType& end_point,
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type&
    peak_velocity,
  const DiffSpace& space, const TimeSpace& t_space, double delta_time = 0.0,
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type
  * best_peak_velocity = NULL ) {
  using std::fabs;

  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1 >::type >::point_type max_velocity
    = get_space< 1 >( space, t_space ).get_upper_corner();
  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 2 >::type >::point_type
    max_acceleration = get_space< 2 >( space, t_space ).get_upper_corner();
  peak_velocity = max_velocity;
  double min_dt_final = 0.0;

  for( std::size_t i = 0; i < max_velocity.size(); ++i ) {
    double vp = 0.0;
    double min_delta_time = sap_Ndof_compute_min_delta_time( get< 0 >( start_point )[i], get< 0 >( end_point )[i],
                                                             get< 1 >( start_point )[i], get< 1 >( end_point )[i], vp,
                                                             max_velocity[i], max_acceleration[i] );
    peak_velocity[i] = vp;

    if( min_dt_final < min_delta_time )
      min_dt_final = min_delta_time;
  };

  if( best_peak_velocity )
    *best_peak_velocity = peak_velocity;

  if( min_dt_final > delta_time )
    delta_time = min_dt_final;

  for( std::size_t i = 0; i < peak_velocity.size(); ++i ) {
    double vp = 0.0;
    sap_Ndof_compute_peak_velocity( get< 0 >( start_point )[i], get< 0 >( end_point )[i], get< 1 >( start_point )[i],
                                    get< 1 >( end_point )[i], vp, max_velocity[i], max_acceleration[i], delta_time );
    peak_velocity[i] = vp;
  };

  return min_dt_final;
};
};
};
};

#endif
