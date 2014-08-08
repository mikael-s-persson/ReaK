/**
 * \file spatial_path_concept.hpp
 * 
 * This library defines the traits and concepts related to a spatial path. A 
 * path is simply a continuous curve in a topology (or space).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_SPATIAL_PATH_CONCEPT_HPP
#define REAK_SPATIAL_PATH_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <cmath>
#include <boost/concept_check.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>

#include <exception>
#include <string>

namespace ReaK {

namespace pp {

  
/**
 * This traits class defines the traits that characterize a spatial path within a 
 * topology.
 * \tparam SpatialPath The spatial path type for which the traits are sought.
 */
template <typename SpatialPath>
struct spatial_path_traits {
  /** This type describes a point in the space or topology. */
  typedef typename SpatialPath::point_type point_type;
  /** This type describes the difference between two points in the space or topology. */
  typedef typename SpatialPath::point_difference_type point_difference_type;
  
  /** This type describes waypoints used by the path to quickly access local parameters of the path. */
  typedef typename SpatialPath::waypoint_descriptor waypoint_descriptor;
  /** This type describes const-waypoints used by the path to quickly access local parameters of the path. */
  typedef typename SpatialPath::const_waypoint_descriptor const_waypoint_descriptor;
  
  /** This type is the topology type in which the path exists. */
  typedef typename SpatialPath::topology topology;
  /** This type is the distance metric type used on the topology and defining the travel distances along the path. */
  typedef typename SpatialPath::distance_metric distance_metric;
  
};


/**
 * This concept class defines the requirements for a type to model a spatial-path as used in 
 * ReaK::pp. A spatial path is a continuous curve within a topology.
 * 
 * Required concepts:
 * 
 * The topology should model the MetricSpaceConcept.
 * 
 * The distance-metric should model the DistanceMetricConcept.
 * 
 * Valid expressions:
 * 
 * std::pair< const_waypoint_descriptor, point_type > w_p;  A waypoint-point pair is used to associate a point and a waypoint.
 * 
 * pt = p.move_away_from(pt, d);  A point along the path (p), at a distance (d) from a point (pt), can be obtained.
 * 
 * d = p.travel_distance(pt,pt);  The travel distance, along the path (p), between two points (pt,pt), can be obtained.
 * 
 * w_p = p.move_away_from(w_p, d);  A waypoint along the path (p), at a distance (d) from a waypoint (w_p), can be obtained.
 * 
 * d = p.travel_distance(w_p,w_p);  The travel distance, along the path (p), between two waypoints (w_p,w_p), can be obtained.
 * 
 * pt = p.get_start_point();  The starting point of the path can be obtained.
 * 
 * w_p = p.get_start_waypoint();  The starting waypoint-pair of the path can be obtained.
 * 
 * pt = p.get_end_point();  The end point of the path can be obtained.
 * 
 * w_p = p.get_end_waypoint();  The end waypoint-pair of the path can be obtained.
 * 
 * \tparam SpatialPath The type to be checked for the requirements of this concept.
 * \tparam Topology The topology in which the spatial-path should reside.
 */
template <typename SpatialPath, typename Topology>
struct SpatialPathConcept {
  
  BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));
  BOOST_CONCEPT_ASSERT((DistanceMetricConcept< typename spatial_path_traits<SpatialPath>::distance_metric, Topology >));
  
  SpatialPath p;
  typename topology_traits<Topology>::point_type pt;
  std::pair< typename spatial_path_traits<SpatialPath>::const_waypoint_descriptor, typename topology_traits<Topology>::point_type> w_p;
  double d;
  BOOST_CONCEPT_USAGE(SpatialPathConcept)
  {
    pt  = p.move_away_from(pt, d);
    d   = p.travel_distance(pt, pt);
    w_p = p.move_away_from(w_p, d);
    d   = p.travel_distance(w_p,w_p);
    
    pt = p.get_start_point();
    w_p = p.get_start_waypoint();
    pt = p.get_end_point();
    w_p = p.get_end_waypoint();
  };
  
};


/**
 * This exception class represents an exception occurring when a path is invalid.
 */
struct invalid_path : public std::exception {
  std::string message;
  invalid_path(const std::string& aSender) : message(std::string("Invalid path reported! Originating from ") + aSender) { };
  virtual ~invalid_path() throw() { };
  const char* what() const throw() {
    return message.c_str();
  };
};






};

};


#endif









