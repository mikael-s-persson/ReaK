/**
 * \file manip_free_workspace.tpp
 * 
 * This library defines a class for path-planning for a manipulator moving inside an environment 
 * with obstacles and physical limits (joint limits). This class is essentially just an assembly 
 * of many of the building blocks in the ReaK path-planning library. This class can also draw the 
 * elements of the motion graph (edges) as end-effector trajectories in a Coin3D scene-graph.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_MANIP_FREE_WORKSPACE_TPP
#define REAK_MANIP_FREE_WORKSPACE_TPP

#include "manip_free_workspace.hpp"


#include "path_planning/random_sampler_concept.hpp"
#include "path_planning/metric_space_concept.hpp"
#include "path_planning/spatial_trajectory_concept.hpp"  // for SpatialTrajectoryConcept

#include "basic_distance_metrics.hpp"
#include "default_random_sampler.hpp"

#include "rate_limited_spaces.hpp"

#include "proximity/proxy_query_model.hpp"  // for proxy-query class

#include "direct_kinematics_topomap.hpp"       // for write_joint_coordinates_impl
#include "joint_space_limits.tpp"              // for joint_limits_collection and create_normal_joint_vectors_impl.  (needed for manip_dynamic_env and manip_quasi_static_env)
#include "kte_models/direct_kinematics_model.hpp"

#include "interpolation/generic_interpolator_factory.hpp"

namespace ReaK {

namespace pp {


namespace detail {


template <typename PointType, typename RateLimitedJointSpace>
bool manip_dk_proxy_env_impl::is_free(const PointType& pt, const RateLimitedJointSpace& space) const {
  typedef typename get_rate_illimited_space< RateLimitedJointSpace >::type NormalJointSpace;
  NormalJointSpace normal_j_space; // dummy
  typename topology_traits< NormalJointSpace >::point_type pt_inter;
  detail::create_normal_joint_vectors_impl(pt_inter, pt, *m_joint_limits_map);
  detail::write_joint_coordinates_impl(pt_inter, normal_j_space, m_model);
  // update the kinematics model with the given joint states.
  m_model->doDirectMotion();
  
  // NOTE: it is assumed that the proxy environments are properly linked with the manipulator kinematic model.
  
  for( std::vector< shared_ptr< geom::proxy_query_pair_2D > >::const_iterator it = m_proxy_env_2D.begin(); it != m_proxy_env_2D.end(); ++it) {
    shared_ptr< geom::proximity_finder_2D > tmp = (*it)->findMinimumDistance();
    if((tmp) && (tmp->getLastResult().mDistance < 0.0))
      return false;
  };
  for( std::vector< shared_ptr< geom::proxy_query_pair_3D > >::const_iterator it = m_proxy_env_3D.begin(); it != m_proxy_env_3D.end(); ++it) {
    shared_ptr< geom::proximity_finder_3D > tmp = (*it)->findMinimumDistance();
    if((tmp) && (tmp->getLastResult().mDistance < 0.0))
      return false;
  };
  
  return true;
};


};


template <typename RateLimitedJointSpace, typename InterpMethodTag>
typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::move_position_toward(const typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p1, double fraction, const typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p2) const {
  typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::type InterpType;
  typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::pseudo_factory_type InterpFactoryType;
  typedef typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type PointType;
  
  InterpType interp;
  double dt_min = m_distance(p1, p2, m_space);
  interp.initialize(p1, p2, dt_min, m_space, time_topology(), InterpFactoryType());
  double dt = dt_min * fraction;
  dt = (dt < max_edge_length ? dt : max_edge_length);
  double d = min_interval;
  PointType result = p1;
  PointType last_result = p1;
  while(d < dt) {
    interp.compute_point(result, p1, p2, m_space, time_topology(), d, dt_min, InterpFactoryType());
    if(!is_free(result))
      return last_result;
    d += min_interval;
    last_result = result;
  };
  if((fraction == 1.0) && (dt_min < max_edge_length)) //these equal comparison are used for when exact end fractions are used.
    return p2;
  else if(fraction == 0.0)
    return p1;
  else {
    interp.compute_point(result, p1, p2, m_space, time_topology(), dt, dt_min, InterpFactoryType());
    return result;
  };
};


template <typename RateLimitedJointSpace, typename InterpMethodTag>
std::pair<typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type, bool> manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::random_walk(const typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p_u) const {
  typedef typename manip_quasi_static_env<RateLimitedJointSpace, InterpMethodTag>::point_type PointType;
  
  PointType p_rnd, p_v;
  unsigned int i = 0;
  do {
    p_rnd = m_rand_sampler(m_space);
    double dist = m_distance(p_u, p_rnd, m_space);
    p_v = move_position_toward(p_u, max_edge_length / dist, p_rnd);
    ++i;
  } while((m_distance(p_u, p_v, m_space) < min_interval) && (i <= 10));
  if(i > 10) {
    //could not expand vertex u, then just output a random C-free point.
    return std::make_pair(p_v, false);
  };
  return std::make_pair(p_v, true);
};


};

};



#endif

