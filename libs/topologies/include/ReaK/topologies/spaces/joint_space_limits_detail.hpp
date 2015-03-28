/**
 * \file joint_space_limits_detail.hpp
 *
 * This library contains the implementation details for classes in the joint_space_limits.hpp file.
 * \note This is not a stand-alone header file, it is simply to be included in the joint_space_limits.hpp header. This
 *"detail" header just serves to tuck away the details (TMP magic).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_JOINT_SPACE_LIMITS_DETAIL_HPP
#define REAK_JOINT_SPACE_LIMITS_DETAIL_HPP


#include "joint_space_topologies.hpp"
#include "se2_topologies.hpp"
#include "se3_topologies.hpp"
#include "Ndof_spaces.hpp"


#include <boost/mpl/and.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/prior.hpp>
#include <boost/mpl/less.hpp>

namespace ReaK {

namespace pp {

namespace detail {
namespace {


/*****************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACES
******************************************************************************************************/


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< OutSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 0 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace( arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type > >(
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 0 >( space_in ).getName() + "_rl",
      ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i],
      ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i] ) ) );
  ++gen_i;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< OutSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 1 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace(
    arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type >,
                      line_segment_topology< typename RateLimitMap::value_type > >(
      line_segment_topology< typename RateLimitMap::value_type >(
        get< 0 >( space_in ).getName() + "_rl",
        ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i],
        ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i] ),
      line_segment_topology< typename RateLimitMap::value_type >(
        get< 1 >( space_in ).getName() + "_rl",
        ( get< 1 >( space_in ).origin() - get< 1 >( space_in ).get_radius() ) / j_limits.gen_accel_limits[gen_i],
        ( get< 1 >( space_in ).origin() + get< 1 >( space_in ).get_radius() ) / j_limits.gen_accel_limits[gen_i] ) ),
    euclidean_tuple_distance(),
    differentiation_rule_array< 1, reach_time_differentiation >::type(
      reach_time_differentiation( j_limits.gen_speed_limits[gen_i] / j_limits.gen_accel_limits[gen_i] ) ) );
  ++gen_i;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< OutSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 2 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace(
    arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type >,
                      line_segment_topology< typename RateLimitMap::value_type >,
                      line_segment_topology< typename RateLimitMap::value_type > >(
      line_segment_topology< typename RateLimitMap::value_type >(
        get< 0 >( space_in ).getName() + "_rl",
        ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i],
        ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) / j_limits.gen_speed_limits[gen_i] ),
      line_segment_topology< typename RateLimitMap::value_type >(
        get< 1 >( space_in ).getName() + "_rl",
        ( get< 1 >( space_in ).origin() - get< 1 >( space_in ).get_radius() ) / j_limits.gen_accel_limits[gen_i],
        ( get< 1 >( space_in ).origin() + get< 1 >( space_in ).get_radius() ) / j_limits.gen_accel_limits[gen_i] ),
      line_segment_topology< typename RateLimitMap::value_type >(
        get< 2 >( space_in ).getName() + "_rl",
        ( get< 2 >( space_in ).origin() - get< 2 >( space_in ).get_radius() ) / j_limits.gen_jerk_limits[gen_i],
        ( get< 2 >( space_in ).origin() + get< 2 >( space_in ).get_radius() ) / j_limits.gen_jerk_limits[gen_i] ) ),
    euclidean_tuple_distance(),
    differentiation_rule_array< 2, reach_time_differentiation >::type(
      reach_time_differentiation( j_limits.gen_speed_limits[gen_i] / j_limits.gen_accel_limits[gen_i] ),
      reach_time_differentiation( j_limits.gen_accel_limits[gen_i] / j_limits.gen_jerk_limits[gen_i] ) ) );
  ++gen_i;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 0 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace( arithmetic_tuple< hyperbox_topology< VectorType, inf_norm_distance_metric > >(
    hyperbox_topology< VectorType, inf_norm_distance_metric >( get< 0 >( space_in ).getName() + "_rl", lower_bnd,
                                                               upper_bnd ) ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 1 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  VectorType speed_lim = get< 1 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    speed_lim[i] /= j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace( arithmetic_tuple< hyperbox_topology< VectorType, inf_norm_distance_metric >,
                                          hyperbox_topology< VectorType, inf_norm_distance_metric > >(
                          hyperbox_topology< VectorType, inf_norm_distance_metric >(
                            get< 0 >( space_in ).getName() + "_rl", lower_bnd, upper_bnd ),
                          hyperbox_topology< VectorType, inf_norm_distance_metric >(
                            get< 1 >( space_in ).getName() + "_rl", -speed_lim, speed_lim ) ),
                        manhattan_tuple_distance(), arithmetic_tuple< Ndof_reach_time_differentiation< VectorType > >(
                                                      Ndof_reach_time_differentiation< VectorType >( speed_lim ) ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 2 > > >,
                           void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  VectorType speed_lim = get< 1 >( space_in ).get_upper_corner();
  VectorType accel_lim = get< 2 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] /= j_limits.gen_speed_limits[gen_i];
    speed_lim[i] /= j_limits.gen_accel_limits[gen_i];
    accel_lim[i] /= j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace(
    arithmetic_tuple< hyperbox_topology< VectorType, inf_norm_distance_metric >,
                      hyperbox_topology< VectorType, inf_norm_distance_metric >,
                      hyperbox_topology< VectorType, inf_norm_distance_metric > >(
      hyperbox_topology< VectorType, inf_norm_distance_metric >( get< 0 >( space_in ).getName() + "_rl", lower_bnd,
                                                                 upper_bnd ),
      hyperbox_topology< VectorType, inf_norm_distance_metric >( get< 1 >( space_in ).getName() + "_rl", -speed_lim,
                                                                 speed_lim ),
      hyperbox_topology< VectorType, inf_norm_distance_metric >( get< 2 >( space_in ).getName() + "_rl", -accel_lim,
                                                                 accel_lim ) ),
    manhattan_tuple_distance(),
    arithmetic_tuple< Ndof_reach_time_differentiation< VectorType >, Ndof_reach_time_differentiation< VectorType > >(
      Ndof_reach_time_differentiation< VectorType >( speed_lim ),
      Ndof_reach_time_differentiation< VectorType >( accel_lim ) ) );
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 0 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef arithmetic_tuple< LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                                           get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                                           * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ),
                                           get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                                           * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ) ) ) ),
    RotTopoOutType( RotTopoTuple(
      LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl",
                   ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   / j_limits.frame2D_speed_limits[f2d_i + 1],
                   ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   / j_limits.frame2D_speed_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 1 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 2 > > BallTopo;
  typedef arithmetic_tuple< LineSegTopo, LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType(
      PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                             get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ),
                             get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ) ),
                    BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 1 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame2D_accel_limits[f2d_i] ),
                              get< 1 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame2D_accel_limits[f2d_i] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 1, reach_time_differentiation >::type(
        reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i] / j_limits.frame2D_accel_limits[f2d_i] ) ) ),
    RotTopoOutType(
      RotTopoTuple(
        LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl",
                     ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_speed_limits[f2d_i + 1],
                     ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_speed_limits[f2d_i + 1] ),
        LineSegTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_rl",
                     ( get< 1 >( get< 1 >( space_in ) ).origin() - get< 1 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_accel_limits[f2d_i + 1],
                     ( get< 1 >( get< 1 >( space_in ) ).origin() + get< 1 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_accel_limits[f2d_i + 1] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 1, reach_time_differentiation >::type( reach_time_differentiation(
        j_limits.frame2D_speed_limits[f2d_i + 1] / j_limits.frame2D_accel_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 2 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 2 > > BallTopo;
  typedef arithmetic_tuple< LineSegTopo, LineSegTopo, LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType(
      PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                             get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ),
                             get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame2D_speed_limits[f2d_i] ) ),
                    BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 1 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame2D_accel_limits[f2d_i] ),
                              get< 1 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame2D_accel_limits[f2d_i] ),
                    BallTopo( get< 2 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 2 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame2D_jerk_limits[f2d_i] ),
                              get< 2 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame2D_jerk_limits[f2d_i] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 2, reach_time_differentiation >::type(
        reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i] / j_limits.frame2D_accel_limits[f2d_i] ),
        reach_time_differentiation( j_limits.frame2D_accel_limits[f2d_i] / j_limits.frame2D_jerk_limits[f2d_i] ) ) ),
    RotTopoOutType(
      RotTopoTuple(
        LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl",
                     ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_speed_limits[f2d_i + 1],
                     ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_speed_limits[f2d_i + 1] ),
        LineSegTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_rl",
                     ( get< 1 >( get< 1 >( space_in ) ).origin() - get< 1 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_accel_limits[f2d_i + 1],
                     ( get< 1 >( get< 1 >( space_in ) ).origin() + get< 1 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_accel_limits[f2d_i + 1] ),
        LineSegTopo( get< 2 >( get< 1 >( space_in ) ).getName() + "_rl",
                     ( get< 2 >( get< 1 >( space_in ) ).origin() - get< 2 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_jerk_limits[f2d_i + 1],
                     ( get< 2 >( get< 1 >( space_in ) ).origin() + get< 2 >( get< 1 >( space_in ) ).get_radius() )
                     / j_limits.frame2D_jerk_limits[f2d_i + 1] ) ),
      euclidean_tuple_distance(), differentiation_rule_array< 2, reach_time_differentiation >::type(
                                    reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i + 1]
                                                                / j_limits.frame2D_accel_limits[f2d_i + 1] ),
                                    reach_time_differentiation( j_limits.frame2D_accel_limits[f2d_i + 1]
                                                                / j_limits.frame2D_jerk_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 0 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef rate_limited_quat_space< ValueType > QuatTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef arithmetic_tuple< QuatTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                                           get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                                           * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ),
                                           get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                                           * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ) ) ) ),
    RotTopoOutType( RotTopoTuple(
      QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl", j_limits.frame3D_speed_limits[f3d_i + 1] ) ) ) ) );
  f3d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 1 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef rate_limited_quat_space< ValueType > QuatTopo;
  typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 3 > > BallTopo;
  typedef arithmetic_tuple< QuatTopo, AngVelTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType(
      PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                             get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ),
                             get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ) ),
                    BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 1 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame3D_accel_limits[f3d_i] ),
                              get< 1 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame3D_accel_limits[f3d_i] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 1, reach_time_differentiation >::type(
        reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i] / j_limits.frame3D_accel_limits[f3d_i] ) ) ),
    RotTopoOutType(
      RotTopoTuple(
        QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl", j_limits.frame3D_speed_limits[f3d_i + 1] ),
        AngVelTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_rl",
                    get< 1 >( get< 1 >( space_in ) ).get_radius() / j_limits.frame3D_accel_limits[f3d_i + 1] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 1, reach_time_differentiation >::type( reach_time_differentiation(
        j_limits.frame3D_speed_limits[f3d_i + 1] / j_limits.frame3D_accel_limits[f3d_i + 1] ) ) ) ) );
  f3d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< OutSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 2 > > >,
             void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits, std::size_t&,
                              std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef rate_limited_quat_space< ValueType > QuatTopo;
  typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
  typedef ang_accel_3D_topology< ValueType > AngAccTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 3 > > BallTopo;
  typedef arithmetic_tuple< QuatTopo, AngVelTopo, AngAccTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType(
      PosTopoTuple( BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_rl",
                             get< 0 >( get< 0 >( space_in ) ).get_lower_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ),
                             get< 0 >( get< 0 >( space_in ) ).get_upper_corner()
                             * ( ValueType( 1.0 ) / j_limits.frame3D_speed_limits[f3d_i] ) ),
                    BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 1 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame3D_accel_limits[f3d_i] ),
                              get< 1 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame3D_accel_limits[f3d_i] ),
                    BallTopo( get< 2 >( get< 0 >( space_in ) ).getName() + "_rl",
                              get< 2 >( get< 0 >( space_in ) ).origin()
                              * ( ValueType( 1.0 ) / j_limits.frame3D_jerk_limits[f3d_i] ),
                              get< 2 >( get< 0 >( space_in ) ).get_radius() / j_limits.frame3D_jerk_limits[f3d_i] ) ),
      euclidean_tuple_distance(),
      differentiation_rule_array< 2, reach_time_differentiation >::type(
        reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i] / j_limits.frame3D_accel_limits[f3d_i] ),
        reach_time_differentiation( j_limits.frame3D_accel_limits[f3d_i] / j_limits.frame3D_jerk_limits[f3d_i] ) ) ),
    RotTopoOutType(
      RotTopoTuple(
        QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_rl", j_limits.frame3D_speed_limits[f3d_i + 1] ),
        AngVelTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_rl",
                    get< 1 >( get< 1 >( space_in ) ).get_radius() / j_limits.frame3D_accel_limits[f3d_i + 1] ),
        AngAccTopo( get< 2 >( get< 1 >( space_in ) ).getName() + "_rl",
                    get< 2 >( get< 1 >( space_in ) ).get_radius() / j_limits.frame3D_jerk_limits[f3d_i + 1] ) ),
      euclidean_tuple_distance(), differentiation_rule_array< 2, reach_time_differentiation >::type(
                                    reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i + 1]
                                                                / j_limits.frame3D_accel_limits[f3d_i + 1] ),
                                    reach_time_differentiation( j_limits.frame3D_accel_limits[f3d_i + 1]
                                                                / j_limits.frame3D_jerk_limits[f3d_i + 1] ) ) ) ) );
  f3d_i += 2;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< OutSpace >,
                                             is_rate_limited_se2_space< OutSpace >,
                                             is_rate_limited_se3_space< OutSpace >, is_Ndof_rl_space< OutSpace > >,
                            void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i );

template < typename RateLimitMap >
struct create_rl_joint_space_functor {
  std::size_t* p_gen_i;
  std::size_t* p_f2d_i;
  std::size_t* p_f3d_i;
  const RateLimitMap* p_j_limits;
  create_rl_joint_space_functor( std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                 const RateLimitMap& j_limits )
      : p_gen_i( &gen_i ), p_f2d_i( &f2d_i ), p_f3d_i( &f3d_i ), p_j_limits( &j_limits ){};
  template < typename OutSpace, typename InSpace >
  void operator()( OutSpace& space_out, const InSpace& space_in ) {
    create_rl_joint_space_impl( space_out, space_in, *p_j_limits, *p_gen_i, *p_f2d_i, *p_f3d_i );
  };
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< OutSpace >,
                                             is_rate_limited_se2_space< OutSpace >,
                                             is_rate_limited_se3_space< OutSpace >, is_Ndof_rl_space< OutSpace > >,
                            void >::type
  create_rl_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                              std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i ) {
  tuple_for_each( space_out, space_in, create_rl_joint_space_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< OutSpace >,
                                             is_rate_limited_se2_space< OutSpace >,
                                             is_rate_limited_se3_space< OutSpace >, is_Ndof_rl_space< OutSpace > >,
                            void >::type
  create_rl_joint_spaces_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  tuple_for_each( space_out, space_in, create_rl_joint_space_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::or_< is_rate_limited_joint_space< OutSpace >,
                                            is_rate_limited_se2_space< OutSpace >,
                                            is_rate_limited_se3_space< OutSpace >, is_Ndof_rl_space< OutSpace > >,
                           void >::type
  create_rl_joint_spaces_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_rl_joint_space_impl( space_out, space_in, j_limits, gen_i, f2d_i, f3d_i );
};


/***********************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACES
************************************************************************************************************/


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 0 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace( arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type > >(
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 0 >( space_in ).getName() + "_non_rl",
      ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i],
      ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i] ) ) );
  ++gen_i;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 1 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace( arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type >,
                                          line_segment_topology< typename RateLimitMap::value_type > >(
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 0 >( space_in ).getName() + "_non_rl",
      ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i],
      ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i] ),
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 1 >( space_in ).getName() + "_non_rl",
      ( get< 1 >( space_in ).origin() - get< 1 >( space_in ).get_radius() ) * j_limits.gen_accel_limits[gen_i],
      ( get< 1 >( space_in ).origin() + get< 1 >( space_in ).get_radius() ) * j_limits.gen_accel_limits[gen_i] ) ) );
  ++gen_i;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_rate_limited_joint_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 2 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {

  space_out = OutSpace( arithmetic_tuple< line_segment_topology< typename RateLimitMap::value_type >,
                                          line_segment_topology< typename RateLimitMap::value_type >,
                                          line_segment_topology< typename RateLimitMap::value_type > >(
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 0 >( space_in ).getName() + "_non_rl",
      ( get< 0 >( space_in ).origin() - get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i],
      ( get< 0 >( space_in ).origin() + get< 0 >( space_in ).get_radius() ) * j_limits.gen_speed_limits[gen_i] ),
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 1 >( space_in ).getName() + "_non_rl",
      ( get< 1 >( space_in ).origin() - get< 1 >( space_in ).get_radius() ) * j_limits.gen_accel_limits[gen_i],
      ( get< 1 >( space_in ).origin() + get< 1 >( space_in ).get_radius() ) * j_limits.gen_accel_limits[gen_i] ),
    line_segment_topology< typename RateLimitMap::value_type >(
      get< 2 >( space_in ).getName() + "_non_rl",
      ( get< 2 >( space_in ).origin() - get< 2 >( space_in ).get_radius() ) * j_limits.gen_jerk_limits[gen_i],
      ( get< 2 >( space_in ).origin() + get< 2 >( space_in ).get_radius() ) * j_limits.gen_jerk_limits[gen_i] ) ) );
  ++gen_i;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_rl_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 0 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace( arithmetic_tuple< hyperbox_topology< VectorType, manhattan_distance_metric > >(
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 0 >( space_in ).getName() + "_non_rl", lower_bnd,
                                                                upper_bnd ) ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_rl_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 1 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  VectorType speed_lim = get< 1 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    speed_lim[i] *= j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace( arithmetic_tuple< hyperbox_topology< VectorType, manhattan_distance_metric >,
                                          hyperbox_topology< VectorType, manhattan_distance_metric > >(
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 0 >( space_in ).getName() + "_non_rl", lower_bnd,
                                                                upper_bnd ),
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 1 >( space_in ).getName() + "_non_rl", -speed_lim,
                                                                speed_lim ) ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::and_< is_Ndof_rl_space< InSpace >,
                                             boost::mpl::equal_to< max_derivation_order< InSpace, time_topology >,
                                                                   boost::mpl::size_t< 2 > > >,
                           void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t&, std::size_t& ) {
  typedef typename derived_N_order_space< InSpace, time_topology, 0 >::type BoxTopo;
  typedef typename BoxTopo::point_type VectorType;
  VectorType lower_bnd = get< 0 >( space_in ).get_lower_corner();
  VectorType upper_bnd = get< 0 >( space_in ).get_upper_corner();
  VectorType speed_lim = get< 1 >( space_in ).get_upper_corner();
  VectorType accel_lim = get< 2 >( space_in ).get_upper_corner();
  for( std::size_t i = 0; i < lower_bnd.size(); ++i ) {
    lower_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    upper_bnd[i] *= j_limits.gen_speed_limits[gen_i];
    speed_lim[i] *= j_limits.gen_accel_limits[gen_i];
    accel_lim[i] *= j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
  space_out = OutSpace( arithmetic_tuple< hyperbox_topology< VectorType, manhattan_distance_metric >,
                                          hyperbox_topology< VectorType, manhattan_distance_metric >,
                                          hyperbox_topology< VectorType, manhattan_distance_metric > >(
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 0 >( space_in ).getName() + "_non_rl", lower_bnd,
                                                                upper_bnd ),
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 1 >( space_in ).getName() + "_non_rl", -speed_lim,
                                                                speed_lim ),
    hyperbox_topology< VectorType, manhattan_distance_metric >( get< 2 >( space_in ).getName() + "_non_rl", -accel_lim,
                                                                accel_lim ) ) );
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 0 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef arithmetic_tuple< LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i] ) ) ),
    RotTopoOutType( RotTopoTuple(
      LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1],
                   ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 1 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 2 > > BallTopo;
  typedef arithmetic_tuple< LineSegTopo, LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i] ),
      BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 1 >( get< 0 >( space_in ) ).origin() * j_limits.frame2D_accel_limits[f2d_i],
                get< 1 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame2D_accel_limits[f2d_i] ) ) ),
    RotTopoOutType( RotTopoTuple(
      LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1],
                   ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1] ),
      LineSegTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 1 >( get< 1 >( space_in ) ).origin() - get< 1 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_accel_limits[f2d_i + 1],
                   ( get< 1 >( get< 1 >( space_in ) ).origin() + get< 1 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_accel_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se2_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 2 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t& f2d_i, std::size_t& ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef line_segment_topology< ValueType > LineSegTopo;
  typedef hyperbox_topology< vect< ValueType, 2 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 2 > > BallTopo;
  typedef arithmetic_tuple< LineSegTopo, LineSegTopo, LineSegTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i] ),
      BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 1 >( get< 0 >( space_in ) ).origin() * j_limits.frame2D_accel_limits[f2d_i],
                get< 1 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame2D_accel_limits[f2d_i] ),
      BallTopo( get< 2 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 2 >( get< 0 >( space_in ) ).origin() * j_limits.frame2D_jerk_limits[f2d_i],
                get< 2 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame2D_jerk_limits[f2d_i] ) ) ),
    RotTopoOutType( RotTopoTuple(
      LineSegTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 0 >( get< 1 >( space_in ) ).origin() - get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1],
                   ( get< 0 >( get< 1 >( space_in ) ).origin() + get< 0 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_speed_limits[f2d_i + 1] ),
      LineSegTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 1 >( get< 1 >( space_in ) ).origin() - get< 1 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_accel_limits[f2d_i + 1],
                   ( get< 1 >( get< 1 >( space_in ) ).origin() + get< 1 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_accel_limits[f2d_i + 1] ),
      LineSegTopo( get< 2 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                   ( get< 2 >( get< 1 >( space_in ) ).origin() - get< 2 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_jerk_limits[f2d_i + 1],
                   ( get< 2 >( get< 1 >( space_in ) ).origin() + get< 2 >( get< 1 >( space_in ) ).get_radius() )
                   * j_limits.frame2D_jerk_limits[f2d_i + 1] ) ) ) ) );
  f2d_i += 2;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 0 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef quaternion_topology< ValueType > QuatTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef arithmetic_tuple< QuatTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i] ) ) ),
    RotTopoOutType( RotTopoTuple( QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl" ) ) ) ) );
  f3d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 1 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef quaternion_topology< ValueType > QuatTopo;
  typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 3 > > BallTopo;
  typedef arithmetic_tuple< QuatTopo, AngVelTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i] ),
      BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 1 >( get< 0 >( space_in ) ).origin() * j_limits.frame3D_accel_limits[f3d_i],
                get< 1 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame3D_accel_limits[f3d_i] ) ) ),
    RotTopoOutType( RotTopoTuple(
      QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl" ),
      AngVelTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                  get< 1 >( get< 1 >( space_in ) ).get_radius() * j_limits.frame3D_accel_limits[f3d_i + 1] ) ) ) ) );
  f3d_i += 2;
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::
  enable_if< boost::mpl::
               and_< is_rate_limited_se3_space< InSpace >,
                     boost::mpl::equal_to< max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type,
                                                                 time_topology >,
                                           boost::mpl::size_t< 2 > > >,
             void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t&, std::size_t&, std::size_t& f3d_i ) {

  typedef typename RateLimitMap::value_type ValueType;
  typedef quaternion_topology< ValueType > QuatTopo;
  typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
  typedef ang_accel_3D_topology< ValueType > AngAccTopo;
  typedef hyperbox_topology< vect< ValueType, 3 > > BoxTopo;
  typedef hyperball_topology< vect< ValueType, 3 > > BallTopo;
  typedef arithmetic_tuple< QuatTopo, AngVelTopo, AngAccTopo > RotTopoTuple;
  typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
  typedef typename arithmetic_tuple_element< 0, OutSpace >::type PosTopoOutType;
  typedef typename arithmetic_tuple_element< 1, OutSpace >::type RotTopoOutType;

  space_out = OutSpace( arithmetic_tuple< PosTopoOutType, RotTopoOutType >(
    PosTopoOutType( PosTopoTuple(
      BoxTopo( get< 0 >( get< 0 >( space_in ) ).getName() + "_non_rl",
               get< 0 >( get< 0 >( space_in ) ).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
               get< 0 >( get< 0 >( space_in ) ).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i] ),
      BallTopo( get< 1 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 1 >( get< 0 >( space_in ) ).origin() * j_limits.frame3D_accel_limits[f3d_i],
                get< 1 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame3D_accel_limits[f3d_i] ),
      BallTopo( get< 2 >( get< 0 >( space_in ) ).getName() + "_non_rl",
                get< 2 >( get< 0 >( space_in ) ).origin() * j_limits.frame3D_jerk_limits[f3d_i],
                get< 2 >( get< 0 >( space_in ) ).get_radius() * j_limits.frame3D_jerk_limits[f3d_i] ) ) ),
    RotTopoOutType( RotTopoTuple(
      QuatTopo( get< 0 >( get< 1 >( space_in ) ).getName() + "_non_rl" ),
      AngVelTopo( get< 1 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                  get< 1 >( get< 1 >( space_in ) ).get_radius() * j_limits.frame3D_accel_limits[f3d_i + 1] ),
      AngAccTopo( get< 2 >( get< 1 >( space_in ) ).getName() + "_non_rl",
                  get< 2 >( get< 1 >( space_in ) ).get_radius() * j_limits.frame3D_jerk_limits[f3d_i + 1] ) ) ) ) );
  f3d_i += 2;
};


template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< InSpace >,
                                             is_rate_limited_se2_space< InSpace >, is_rate_limited_se3_space< InSpace >,
                                             is_Ndof_rl_space< InSpace > >,
                            void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i );

template < typename RateLimitMap >
struct create_normal_joint_space_functor {
  std::size_t* p_gen_i;
  std::size_t* p_f2d_i;
  std::size_t* p_f3d_i;
  const RateLimitMap* p_j_limits;
  create_normal_joint_space_functor( std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                     const RateLimitMap& j_limits )
      : p_gen_i( &gen_i ), p_f2d_i( &f2d_i ), p_f3d_i( &f3d_i ), p_j_limits( &j_limits ){};
  template < typename OutSpace, typename InSpace >
  void operator()( OutSpace& space_out, const InSpace& space_in ) {
    create_normal_joint_space_impl( space_out, space_in, *p_j_limits, *p_gen_i, *p_f2d_i, *p_f3d_i );
  };
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< InSpace >,
                                             is_rate_limited_se2_space< InSpace >, is_rate_limited_se3_space< InSpace >,
                                             is_Ndof_rl_space< InSpace > >,
                            void >::type
  create_normal_joint_space_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits,
                                  std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i ) {
  tuple_for_each( space_out, space_in,
                  create_normal_joint_space_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::disable_if< boost::mpl::or_< is_rate_limited_joint_space< InSpace >,
                                             is_rate_limited_se2_space< InSpace >, is_rate_limited_se3_space< InSpace >,
                                             is_Ndof_rl_space< InSpace > >,
                            void >::type
  create_normal_joint_spaces_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  tuple_for_each( space_out, space_in,
                  create_normal_joint_space_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutSpace, typename InSpace, typename RateLimitMap >
typename boost::enable_if< boost::mpl::or_< is_rate_limited_joint_space< InSpace >,
                                            is_rate_limited_se2_space< InSpace >, is_rate_limited_se3_space< InSpace >,
                                            is_Ndof_rl_space< InSpace > >,
                           void >::type
  create_normal_joint_spaces_impl( OutSpace& space_out, const InSpace& space_in, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_normal_joint_space_impl( space_out, space_in, j_limits, gen_i, f2d_i, f3d_i );
};


/**************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACE VECTORS
***************************************************************************************************/


template < typename RateLimitMap >
void create_rl_joint_vector_impl( arithmetic_tuple< typename RateLimitMap::value_type >& result,
                                  const arithmetic_tuple< typename RateLimitMap::value_type >& pt,
                                  const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) / j_limits.gen_speed_limits[gen_i];
  ++gen_i;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type >& result,
  const arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type >& pt,
  const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) / j_limits.gen_speed_limits[gen_i];
  get< 1 >( result ) = get< 1 >( pt ) / j_limits.gen_accel_limits[gen_i];
  ++gen_i;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                    typename RateLimitMap::value_type >& result,
  const arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                          typename RateLimitMap::value_type >& pt,
  const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) / j_limits.gen_speed_limits[gen_i];
  get< 1 >( result ) = get< 1 >( pt ) / j_limits.gen_accel_limits[gen_i];
  get< 2 >( result ) = get< 2 >( pt ) / j_limits.gen_jerk_limits[gen_i];
  ++gen_i;
};


template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_rl_joint_vectors_impl( arithmetic_tuple< Vector >& result, const arithmetic_tuple< Vector >& pt,
                                const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] / j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };
};

template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_rl_joint_vectors_impl( arithmetic_tuple< Vector, Vector >& result,
                                const arithmetic_tuple< Vector, Vector >& pt, const RateLimitMap& j_limits,
                                std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] / j_limits.gen_speed_limits[gen_i];
    get< 1 >( result )[i] = get< 1 >( pt )[i] / j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };
};

template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_rl_joint_vectors_impl( arithmetic_tuple< Vector, Vector, Vector >& result,
                                const arithmetic_tuple< Vector, Vector, Vector >& pt, const RateLimitMap& j_limits,
                                std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] / j_limits.gen_speed_limits[gen_i];
    get< 1 >( result )[i] = get< 1 >( pt )[i] / j_limits.gen_accel_limits[gen_i];
    get< 2 >( result )[i] = get< 2 >( pt )[i] / j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
};


template < typename RateLimitMap >
void
  create_rl_joint_vector_impl( arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 > >,
                                                 arithmetic_tuple< typename RateLimitMap::value_type > >& result,
                               const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 > >,
                                                       arithmetic_tuple< typename RateLimitMap::value_type > >& pt,
                               const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_speed_limits[f2d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) / j_limits.frame2D_speed_limits[f2d_i + 1];
  f2d_i += 2;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 > >,
                    arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 > >,
                          arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type > >&
    pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_speed_limits[f2d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) / j_limits.frame2D_speed_limits[f2d_i + 1];
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_accel_limits[f2d_i] );
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) / j_limits.frame2D_accel_limits[f2d_i + 1];
  f2d_i += 2;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 > >,
                    arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                                      typename RateLimitMap::value_type > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 > >,
                          arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                                            typename RateLimitMap::value_type > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_speed_limits[f2d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) / j_limits.frame2D_speed_limits[f2d_i + 1];
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_accel_limits[f2d_i] );
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) / j_limits.frame2D_accel_limits[f2d_i + 1];
  get< 2 >( get< 0 >( result ) ) = get< 2 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame2D_jerk_limits[f2d_i] );
  get< 2 >( get< 1 >( result ) ) = get< 2 >( get< 1 >( pt ) ) / j_limits.frame2D_jerk_limits[f2d_i + 1];
  f2d_i += 2;
};


template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_speed_limits[f3d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  f3d_i += 2;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                      vect< typename RateLimitMap::value_type, 3 > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                            vect< typename RateLimitMap::value_type, 3 > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_speed_limits[f3d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_accel_limits[f3d_i] );
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * ( 1.0 / j_limits.frame3D_accel_limits[f3d_i + 1] );
  f3d_i += 2;
};

template < typename RateLimitMap >
void create_rl_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                      vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                            vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_speed_limits[f3d_i] );
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_accel_limits[f3d_i] );
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * ( 1.0 / j_limits.frame3D_accel_limits[f3d_i + 1] );
  get< 2 >( get< 0 >( result ) ) = get< 2 >( get< 0 >( pt ) ) * ( 1.0 / j_limits.frame3D_jerk_limits[f3d_i] );
  get< 2 >( get< 1 >( result ) ) = get< 2 >( get< 1 >( pt ) ) * ( 1.0 / j_limits.frame3D_jerk_limits[f3d_i + 1] );
  f3d_i += 2;
};


template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_rl_joint_vector_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits, std::size_t& gen_i,
                                  std::size_t& f2d_i, std::size_t& f3d_i );

template < typename RateLimitMap >
struct create_rl_joint_vector_functor {
  std::size_t* p_gen_i;
  std::size_t* p_f2d_i;
  std::size_t* p_f3d_i;
  const RateLimitMap* p_j_limits;
  create_rl_joint_vector_functor( std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                  const RateLimitMap& j_limits )
      : p_gen_i( &gen_i ), p_f2d_i( &f2d_i ), p_f3d_i( &f3d_i ), p_j_limits( &j_limits ){};
  template < typename OutPoint, typename InPoint >
  void operator()( OutPoint& result, const InPoint& pt ) {
    create_rl_joint_vector_impl( result, pt, *p_j_limits, *p_gen_i, *p_f2d_i, *p_f3d_i );
  };
};

template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_rl_joint_vector_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits, std::size_t& gen_i,
                                  std::size_t& f2d_i, std::size_t& f3d_i ) {
  tuple_for_each( result, pt, create_rl_joint_vector_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_rl_joint_vectors_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  tuple_for_each( result, pt, create_rl_joint_vector_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};


/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACE VECTORS
*******************************************************************************************************************/


template < typename RateLimitMap >
void create_normal_joint_vector_impl( arithmetic_tuple< typename RateLimitMap::value_type >& result,
                                      const arithmetic_tuple< typename RateLimitMap::value_type >& pt,
                                      const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) * j_limits.gen_speed_limits[gen_i];
  ++gen_i;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type >& result,
  const arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type >& pt,
  const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) * j_limits.gen_speed_limits[gen_i];
  get< 1 >( result ) = get< 1 >( pt ) * j_limits.gen_accel_limits[gen_i];
  ++gen_i;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                    typename RateLimitMap::value_type >& result,
  const arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                          typename RateLimitMap::value_type >& pt,
  const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  get< 0 >( result ) = get< 0 >( pt ) * j_limits.gen_speed_limits[gen_i];
  get< 1 >( result ) = get< 1 >( pt ) * j_limits.gen_accel_limits[gen_i];
  get< 2 >( result ) = get< 2 >( pt ) * j_limits.gen_jerk_limits[gen_i];
  ++gen_i;
};


template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_normal_joint_vectors_impl( arithmetic_tuple< Vector >& result, const arithmetic_tuple< Vector >& pt,
                                    const RateLimitMap& j_limits, std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] * j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };
};

template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_normal_joint_vectors_impl( arithmetic_tuple< Vector, Vector >& result,
                                    const arithmetic_tuple< Vector, Vector >& pt, const RateLimitMap& j_limits,
                                    std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] * j_limits.gen_speed_limits[gen_i];
    get< 1 >( result )[i] = get< 1 >( pt )[i] * j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };
};

template < typename Idx, typename Vector, typename RateLimitMap >
typename boost::enable_if< is_writable_vector< Vector >, void >::type
  create_normal_joint_vectors_impl( arithmetic_tuple< Vector, Vector, Vector >& result,
                                    const arithmetic_tuple< Vector, Vector, Vector >& pt, const RateLimitMap& j_limits,
                                    std::size_t& gen_i, std::size_t&, std::size_t& ) {
  for( std::size_t i = 0; i < get< 0 >( pt ).size(); ++i ) {
    get< 0 >( result )[i] = get< 0 >( pt )[i] * j_limits.gen_speed_limits[gen_i];
    get< 1 >( result )[i] = get< 1 >( pt )[i] * j_limits.gen_accel_limits[gen_i];
    get< 2 >( result )[i] = get< 2 >( pt )[i] * j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
};


template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 > >,
                    arithmetic_tuple< typename RateLimitMap::value_type > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 > >,
                          arithmetic_tuple< typename RateLimitMap::value_type > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i + 1];
  f2d_i += 2;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 > >,
                    arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 > >,
                          arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type > >&
    pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i + 1];
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * j_limits.frame2D_accel_limits[f2d_i];
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * j_limits.frame2D_accel_limits[f2d_i + 1];
  f2d_i += 2;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 >,
                                      vect< typename RateLimitMap::value_type, 2 > >,
                    arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                                      typename RateLimitMap::value_type > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 >,
                                            vect< typename RateLimitMap::value_type, 2 > >,
                          arithmetic_tuple< typename RateLimitMap::value_type, typename RateLimitMap::value_type,
                                            typename RateLimitMap::value_type > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t& f2d_i, std::size_t& ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) ) * j_limits.frame2D_speed_limits[f2d_i + 1];
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * j_limits.frame2D_accel_limits[f2d_i];
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * j_limits.frame2D_accel_limits[f2d_i + 1];
  get< 2 >( get< 0 >( result ) ) = get< 2 >( get< 0 >( pt ) ) * j_limits.frame2D_jerk_limits[f2d_i];
  get< 2 >( get< 1 >( result ) ) = get< 2 >( get< 1 >( pt ) ) * j_limits.frame2D_jerk_limits[f2d_i + 1];
  f2d_i += 2;
};


template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame3D_speed_limits[f3d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  f3d_i += 2;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                      vect< typename RateLimitMap::value_type, 3 > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                            vect< typename RateLimitMap::value_type, 3 > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame3D_speed_limits[f3d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * j_limits.frame3D_accel_limits[f3d_i];
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * j_limits.frame3D_accel_limits[f3d_i + 1];
  f3d_i += 2;
};

template < typename RateLimitMap >
void create_normal_joint_vector_impl(
  arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > >,
                    arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                      vect< typename RateLimitMap::value_type, 3 >,
                                      vect< typename RateLimitMap::value_type, 3 > > >& result,
  const arithmetic_tuple< arithmetic_tuple< vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > >,
                          arithmetic_tuple< unit_quat< typename RateLimitMap::value_type >,
                                            vect< typename RateLimitMap::value_type, 3 >,
                                            vect< typename RateLimitMap::value_type, 3 > > >& pt,
  const RateLimitMap& j_limits, std::size_t&, std::size_t&, std::size_t& f3d_i ) {
  get< 0 >( get< 0 >( result ) ) = get< 0 >( get< 0 >( pt ) ) * j_limits.frame3D_speed_limits[f3d_i];
  get< 0 >( get< 1 >( result ) ) = get< 0 >( get< 1 >( pt ) );
  get< 1 >( get< 0 >( result ) ) = get< 1 >( get< 0 >( pt ) ) * j_limits.frame3D_accel_limits[f3d_i];
  get< 1 >( get< 1 >( result ) ) = get< 1 >( get< 1 >( pt ) ) * j_limits.frame3D_accel_limits[f3d_i + 1];
  get< 2 >( get< 0 >( result ) ) = get< 2 >( get< 0 >( pt ) ) * j_limits.frame3D_jerk_limits[f3d_i];
  get< 2 >( get< 1 >( result ) ) = get< 2 >( get< 1 >( pt ) ) * j_limits.frame3D_jerk_limits[f3d_i + 1];
  f3d_i += 2;
};


template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_normal_joint_vector_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits,
                                      std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i );

template < typename RateLimitMap >
struct create_normal_joint_vector_functor {
  std::size_t* p_gen_i;
  std::size_t* p_f2d_i;
  std::size_t* p_f3d_i;
  const RateLimitMap* p_j_limits;
  create_normal_joint_vector_functor( std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                      const RateLimitMap& j_limits )
      : p_gen_i( &gen_i ), p_f2d_i( &f2d_i ), p_f3d_i( &f3d_i ), p_j_limits( &j_limits ){};
  template < typename OutPoint, typename InPoint >
  void operator()( OutPoint& result, const InPoint& pt ) {
    create_normal_joint_vector_impl( result, pt, *p_j_limits, *p_gen_i, *p_f2d_i, *p_f3d_i );
  };
};

template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_normal_joint_vector_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits,
                                      std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i ) {
  tuple_for_each( result, pt, create_normal_joint_vector_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};

template < typename OutPoint, typename InPoint, typename RateLimitMap >
void create_normal_joint_vectors_impl( OutPoint& result, const InPoint& pt, const RateLimitMap& j_limits ) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  tuple_for_each( result, pt, create_normal_joint_vector_functor< RateLimitMap >( gen_i, f2d_i, f3d_i, j_limits ) );
};
};
};
};
};

#endif
