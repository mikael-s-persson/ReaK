
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

#ifndef SPATIAL_TRAJECTORY_CONCEPT_HPP
#define SPATIAL_TRAJECTORY_CONCEPT_HPP

#include "spatial_path_concept.hpp"
#include "temporal_space.hpp"

#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {

template <typename SpatialTrajectory>
struct spatial_trajectory_traits {
  typedef typename SpatialTrajectory::point_type point_type;
  typedef typename SpatialTrajectory::point_difference_type point_difference_type;
  
  typedef typename SpatialTrajectory::waypoint_descriptor waypoint_descriptor;
  typedef typename SpatialTrajectory::const_waypoint_descriptor const_waypoint_descriptor;
  
  typedef typename SpatialTrajectory::topology topology;
  typedef typename temporal_topology_traits<topology>::time_topology time_topology;
  typedef typename temporal_topology_traits<topology>::space_topology space_topology;
  typedef typename SpatialTrajectory::distance_metric distance_metric;
  
  
};


template <typename SpatialTrajectory, typename Topology>
struct SpatialTrajectoryConcept {
  SpatialTrajectory p;
  typename temporal_topology_traits<Topology>::point_type pt;
  std::pair< typename spatial_path_traits<SpatialTrajectory>::const_waypoint_descriptor, 
             typename temporal_topology_traits<Topology>::point_type> w_p;
  typedef typename temporal_topology_traits<Topology>::time_topology time_topology;
  typename time_topology::point_difference_type dt;
  typename time_topology::point_type t;
  double d;
  void constraints() {
    boost::function_requires< TemporalSpaceConcept<Topology> >();
    boost::function_requires< TemporalDistMetricConcept< typename spatial_trajectory_traits<SpatialTrajectory>::distance_metric, Topology > >();
    pt = p.move_away_from(pt, dt);
    d = p.travel_distance(pt, pt);
    pt = p.get_point(t);
    w_p = p.move_away_from(w_p, dt);
    d = p.travel_distance(w_p, w_p);
    w_p = p.get_point(t);
  };
  
};


};

};

#endif














