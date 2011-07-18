
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







