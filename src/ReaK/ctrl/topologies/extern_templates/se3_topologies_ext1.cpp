
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

#include "base/defs.hpp"

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "topologies/se3_topologies.hpp"

#include "topologies/joint_space_limits.tpp"

namespace ReaK {

namespace pp {
  
// se3_0th_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

template se3_0th_order_topology<double>::type make_se3_space(const std::string& aName, 
                                                             const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner);
// se3_1st_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double>, ang_velocity_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

template se3_1st_order_topology<double>::type make_se3_space(const std::string& aName,
                                                             const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner,
                                                             const double& aMaxSpeed, const double& aMaxAngularSpeed);

// se3_2nd_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double>, ang_velocity_3D_topology<double>, ang_accel_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

template se3_2nd_order_topology<double>::type make_se3_space(const std::string& aName,
                                                             const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner,
                                                             const double& aMaxSpeed, const double& aMaxAngularSpeed,
                                                             const double& aMaxAcceleration, const double& aMaxAngularAccel);

// se3_0th_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

template se3_0th_order_rl_topology<double>::type make_rl_se3_space(const std::string& aName,
                                                                   const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner, 
                                                                   const double& aMaxSpeed, const double& aMaxAngularSpeed);

// se3_1st_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

template se3_1st_order_rl_topology<double>::type make_rl_se3_space(const std::string& aName,
                                                                   const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner, 
                                                                   const double& aMaxSpeed, const double& aMaxAngularSpeed,
                                                                   const double& aMaxAcceleration, const double& aMaxAngularAccel);

// se3_2nd_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double>, ang_accel_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >; 

template se3_2nd_order_rl_topology<double>::type make_rl_se3_space(const std::string& aName,
                                                                   const vect<double,3>& aMinCorner, const vect<double,3>& aMaxCorner,
                                                                   const double& aMaxSpeed, const double& aMaxAngularSpeed,
                                                                   const double& aMaxAcceleration, const double& aMaxAngularAccel,
                                                                   const double& aMaxJerk, const double& aMaxAngularJerk);


// se3_0th_order_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se3_1st_order_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double>, ang_velocity_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se3_2nd_order_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology<double>, ang_velocity_3D_topology<double>, ang_accel_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;


// se3_0th_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se3_1st_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se3_2nd_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double>, ang_accel_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>; 


// se3_0th_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;

// se3_1st_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;

// se3_2nd_order_rl_topology
template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,3> >, hyperball_topology< vect<double,3> >, hyperball_topology< vect<double,3> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space<double>, ang_velocity_3D_topology<double>, ang_accel_3D_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>; 





template metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_0th_order_topology<double>::type, 1>::type&) const;
template metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_1st_order_topology<double>::type, 1>::type&) const;
template metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type&) const;

template metric_space_array< se3_0th_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type&) const;
template metric_space_array< se3_1st_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type&) const;
template metric_space_array< se3_2nd_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type&) const;

template 
topology_traits< metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_0th_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_0th_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type& ) const;
template 
topology_traits< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_1st_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_1st_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type& ) const;
template 
topology_traits< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_2nd_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type& ) const;

template 
topology_traits< metric_space_array< se3_0th_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_0th_order_topology<double>::type, 1>::type& ) const;
template 
topology_traits< metric_space_array< se3_1st_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_1st_order_topology<double>::type, 1>::type& ) const;
template 
topology_traits< metric_space_array< se3_2nd_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type& ) const;



};



template frame_3D<double> get_frame_3D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template frame_3D<double> get_frame_3D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt);

template frame_3D<double> get_frame_3D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                          arithmetic_tuple< unit_quat<double> > >& pt);

template void set_frame_3D(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const frame_3D<double>& p);

template void set_frame_3D(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt,
  const frame_3D<double>& p);

template void set_frame_3D(
  arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                    arithmetic_tuple< unit_quat<double> > >& pt,
  const frame_3D<double>& p);

template pose_3D<double> get_pose_3D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                          arithmetic_tuple< unit_quat<double> > >& pt);

template void set_pose_3D(
  arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                    arithmetic_tuple< unit_quat<double> > >& pt,
  const pose_3D<double>& p);

template const unit_quat<double>& get_quaternion(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template const unit_quat<double>& get_quaternion(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt);

template const unit_quat<double>& get_quaternion(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                          arithmetic_tuple< unit_quat<double> > >& pt);

template void set_quaternion(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const unit_quat<double>& q);

template void set_quaternion(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt,
  const unit_quat<double>& q);

template void set_quaternion(
  arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                    arithmetic_tuple< unit_quat<double> > >& pt,
  const unit_quat<double>& q);

template const vect<double,3>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template const vect<double,3>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt);

template const vect<double,3>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                          arithmetic_tuple< unit_quat<double> > >& pt);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const vect<double,3>& p);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt,
  const vect<double,3>& p);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,3> >,
                    arithmetic_tuple< unit_quat<double> > >& pt,
  const vect<double,3>& p);

template const vect<double,3>& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template const vect<double,3>& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt);

template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const vect<double,3>& p);

template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt,
  const vect<double,3>& p);

template const vect<double,3>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template const vect<double,3>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt);

template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const vect<double,3>& p);

template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3> > >& pt,
  const vect<double,3>& p);


template const vect<double,3>& get_ang_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template void set_ang_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const vect<double,3>& p);



template const vect<double,3>& get_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                          arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt);

template void set_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,3>,    vect<double,3>, vect<double,3> >,
                    arithmetic_tuple< unit_quat<double>, vect<double,3>, vect<double,3> > >& pt,
  const vect<double,3>& p);


};

#else

namespace ReaK {

namespace pp {

void dummy_se3_topologies_externs_1_symbol() { };

};

};

#endif














