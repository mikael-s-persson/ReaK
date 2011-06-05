
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

#ifndef REACHABILITY_SPACE_CONCEPT_HPP
#define REACHABILITY_SPACE_CONCEPT_HPP

#include "temporal_space_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {


template <typename ReachabilityTopology>
struct reachability_topology_traits {
  typedef ReachabilityTopology::point_type point_type;
  typedef ReachabilityTopology::point_difference_type point_difference_type;
  
  typedef ReachabilityTopology::temporal_space_type temporal_space_type;
  typedef ReachabilityTopology::time_topology time_topology;
  typedef ReachabilityTopology::space_topology space_topology;
  typedef ReachabilityTopology::distance_metric distance_metric;
  
  BOOST_STATIC_CONSTANT(std::size_t, time_dimensions = time_topology::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, space_dimensions = space_topology::point_type::dimensions);
  
};


template <typename Topology>
struct ReachabilitySpaceConcept {
  typename reachability_topology_traits<Topology>::point_type p1;
  typename reachability_topology_traits<Topology>::point_difference_type pd;
  Topology reachable_space;
  double d;
  void constraints() {
    boost::function_requires< TemporalSpaceConcept< typename reachability_topology_traits<Topology>::temporal_space_type > >();
    d = reachable_space.forward_reach(p1);
    d = reachable_space.backward_reach(p1);
    d = reachable_space.forward_norm(pd);
    d = reachable_space.backward_norm(pd);
  };
  
};


};

};

#endif



