/**
 * \file predicted_trajectory_concept.hpp
 * 
 * This library defines the concept that represents a predicted trajectory. This concept 
 * supplements the SpatialTrajectoryConcept with a few additional requirements which 
 * characterize a predicted trajectory (see PredictedTrajectoryConcept).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef PREDICTED_TRAJECTORY_CONCEPT_HPP
#define PREDICTED_TRAJECTORY_CONCEPT_HPP

#include "spatial_path_concept.hpp"
#include "temporal_space.hpp"
#include "spatial_trajectory_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {


/**
 * This concept represents a predicted trajectory, as used in ReaK::pp. This concept 
 * supplements the SpatialTrajectoryConcept with a few additional requirements which 
 * characterize a predicted trajectory (see PredictedTrajectoryConcept).
 * 
 * Required concepts:
 * 
 * The topology should model the TemporalSpaceConcept.
 * 
 * The trajectory should model the SpatialTrajectoryConcept on the given topology.
 * 
 * Valid expressions:
 * 
 * std::pair< const_waypoint_descriptor, point_type> w_p;  A waypoint is a pair of a const-waypoint descriptor and a point on the topology.
 * 
 * p.set_initial_point(pt, t);  The initial point (pt) at time (t) can be set to seed the predicted trajectory (p).
 * 
 * p.set_initial_point(w_p, t);  The initial waypoint (w_p) at time (t) can be set to seed the predicted trajectory (p).
 * 
 * \tparam PredictedTrajectory The trajectory type to be checked for compliance to this concept.
 * \tparam Topology The temporal-topology type that can contain the trajectory.
 */
template <typename PredictedTrajectory, typename Topology>
struct PredictedTrajectoryConcept {
  PredictedTrajectory p;
  typename temporal_topology_traits<Topology>::point_type pt;
  std::pair< typename spatial_path_traits<PredictedTrajectory>::const_waypoint_descriptor, 
             typename temporal_topology_traits<Topology>::point_type> w_p;
  typedef typename temporal_topology_traits<Topology>::time_topology time_topology;
  time_topology::point_type t;
  void constraints() {
    boost::function_requires< SpatialTrajectoryConcept<PredictedTrajectory,Topology> >();
    p.set_initial_point(pt, t);
    p.set_initial_point(w_p, t);
  };
  
};




};

}; 



#endif







