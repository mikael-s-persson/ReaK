
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

#ifndef SPATIAL_PATH_CONCEPT_HPP
#define SPATIAL_PATH_CONCEPT_HPP


#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

#include <exception>
#include <string>

namespace ReaK {

namespace pp {

  
  
template <typename SpatialPath>
struct spatial_path_traits {
  typedef SpatialPath::point_type point_type;
  typedef SpatialPath::point_difference_type point_difference_type;
  
  typedef SpatialPath::waypoint_descriptor waypoint_descriptor;
  typedef SpatialPath::const_waypoint_descriptor const_waypoint_descriptor;
  
  typedef SpatialPath::topology topology;
  typedef SpatialPath::distance_metric distance_metric;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology::point_type::dimensions);
  
};



template <typename SpatialPath, typename Topology>
struct SpatialPathConcept {
  SpatialPath p;
  Topology::point_type pt;
  std::pair< typename spatial_path_traits<SpatialPath>::const_waypoint_descriptor, Topology::point_type> w_p;
  double d;
  void constraints() {
    boost::function_requires< MetricSpaceConcept<Topology> >();
    boost::function_requires< DistanceMetricConcept< typename spatial_path_traits<SpatialPath>::distance_metric, Topology > >();
    pt  = p.move_away_from(pt, d);
    d   = p.travel_distance(pt, pt);
    w_p = p.move_away_from(w_p, d);
    d   = p.travel_distance(w_p,w_p);
  };
  
};



struct invalid_path : public std::exception {
  std::string message;
  invalid_path(const std::string& aSender) : message(std::string("Invalid path reported! Originating from ") + aSender) { };
  const char* what() const throw() {
    return message.c_str();
  };
};






};

};


#endif









