/**
 * \file manip_free_dynamic_workspace.tpp
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

#ifndef REAK_MANIP_FREE_DYNAMIC_WORKSPACE_TPP
#define REAK_MANIP_FREE_DYNAMIC_WORKSPACE_TPP

#include "manip_free_dynamic_workspace.hpp"
#include "manip_free_workspace.tpp"


namespace ReaK {

namespace pp {

  
template <typename RateLimitedJointSpace, typename InterpMethodTag>
typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::move_position_toward(const typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p1, double fraction, const typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p2) const {
  if(p1.time > p2.time) // Am I trying to go backwards in time (impossible)?
    return p1; //p2 is not reachable from p1.
  
  typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::type InterpType;
  typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::pseudo_factory_type InterpFactoryType;
  typedef typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type PointType;
  
  InterpType interp;
  double reach_time = m_distance(p1.pt, p2.pt, m_space.get_space_topology());
  double dt_total = (p2.time - p1.time);  // the free time that I have along the path.
  if(dt_total < reach_time) // There is not enough time to reach the end-point.
    return p1;
  interp.initialize(p1.pt, p2.pt, (p2.time - p1.time), m_space.get_space_topology(), m_space.get_time_topology(), InterpFactoryType());
  double dt = dt_total * fraction;
  dt = (dt < max_edge_length ? dt : max_edge_length);
  double d = min_interval;
  PointType result = p1;
  PointType last_result = p1;
  while(d < dt) {
    interp.compute_point(result.pt, p1.pt, p2.pt, m_space.get_space_topology(), m_space.get_time_topology(), d, dt_total, InterpFactoryType());
    result.time = p1.time + d;
    if(!is_free(result))
      return last_result;
    d += min_interval;
    last_result = result;
  };
  if((fraction == 1.0) && (dt_total < max_edge_length)) //these equal comparison are used for when exact end fractions are used.
    return p2;
  else if(fraction == 0.0)
    return p1;
  else {
    interp.compute_point(result.pt, p1.pt, p2.pt, m_space.get_space_topology(), m_space.get_time_topology(), dt, dt_total, InterpFactoryType());
    result.time = p1.time + dt;
    return result;
  };
};

template <typename RateLimitedJointSpace, typename InterpMethodTag>
std::pair<typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type, bool> manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::random_walk(const typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type& p_u) const {
  typedef typename manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag>::point_type PointType;
  
  PointType p_rnd, p_v;
  unsigned int i = 0;
  do {
    p_rnd = point_type(0.0, m_rand_sampler(m_space.get_space_topology()));
    double reach_time = m_distance(p_u.pt, p_rnd.pt, m_space.get_space_topology());
    p_rnd.time = m_space.get_time_topology().random_point() + reach_time + p_u.time;
    p_v = move_position_toward(p_u, max_edge_length / reach_time, p_rnd);
    ++i;
  } while((p_v.time - p_u.time < min_interval) && (i <= 20));
  if(i > 20) {
    //could not expand vertex u, then just output a arbitrary C-free point.
    return std::make_pair(p_v, false);
  };
  return std::make_pair(p_v, true);
};


};

};



#endif

