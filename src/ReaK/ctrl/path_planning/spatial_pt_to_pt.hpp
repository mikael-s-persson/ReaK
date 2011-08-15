/**
 * \file spatial_pt_to_pt.hpp
 * 
 * This library provides an implementation of a simple path within a topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation.
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

#ifndef SPATIAL_PATH_HPP
#define SPATIAL_PATH_HPP

#include "spatial_path_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>

namespace ReaK {

namespace pp {

  
/**
 * This class implementats a simple path within a topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 * \tparam DistanceMetric The distance metric type used over the topology, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
struct point_to_point_path {
  public:
    typedef point_to_point_path<Topology,DistanceMetric> self;
    typedef Topology topology;
    typedef Topology::point_type point_type;
    typedef Topology::point_difference_type point_difference_type;
    
    typedef typename std::list<point_type>::iterator waypoint_descriptor;
    typedef typename std::list<point_type>::const_iterator const_waypoint_descriptor;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
    typedef DistanceMetric distance_metric;
    
  private:
    
    const topology& space;
    point_type start_point;
    point_type end_point;
    
    std::list<point_type> waypoints;
    
    typedef std::pair<const_waypoint_descriptor, const_waypoint_descriptor> const_waypoint_bounds;
    
    point_to_point_path(const point_to_point_path<Topology,DistanceMetric>&); //non-copyable.
    point_to_point_path<Topology,DistanceMetric>& operator=(const point_to_point_path<Topology,DistanceMetric>&);
    
    void constraints() {
      boost::function_requires< SpatialPathConcept<point_to_point_path<Topology,DistanceMetric>, Topology> >();
    };
    
    const_waypoint_bounds get_waypoint_bounds(const point_type& p, const_waypoint_descriptor start) const {
      const_waypoint_descriptor it_end = waypoints.end();
      if(start == it_end)
	throw invalid_path("Point-to-point path (exhausted waypoints during waypoint query)");
      double d1 = space.distance(*start,p);
      const_waypoint_descriptor it2 = start; ++it2;
      if(it2 == it_end)
	return std::make_pair(start,start);
      double d2 = space.distance(*it2,p);
      const_waypoint_descriptor it3 = it2; ++it3;
      while(it3 != it_end) {
	double d3 = space.distance(*it3,p);
	if(d3 > d2) {
	  if(d3 < d1) {
	    start = it2;
	    it2 = it3;
	    break;
	  } else
	    break;
	} else {
	  start = it2; d1 = d2;
	  it2 = it3; d2 = d3;
	  ++it3;
	};
      };
      return std::make_pair(start,it2);
    };
    
    double travel_distance_impl(const point_type& a, const const_waypoint_bounds& wpb_a, 
				const point_type& b, const const_waypoint_bounds& wpb_b) const {
      if((wpb_a.first == wpb_a.second) || (wpb_a.first == wpb_b.first))
	return space.distance(a,b); //this means that a is at the end of the path, thus, the "path" goes directly towards b.
      double sum = space.distance(a, *(wpb_a.second));
      const_waypoint_descriptor it = wpb_a.second;
      while(it != wpb_b.first) {
	const_waypoint_descriptor it_prev = it;
	sum += space.distance(*it_prev, *(++it));
      };
      sum += space.distance(*(wpb_b.first), b);
      return sum;
    };
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     */
    explicit point_to_point_path(const topology& aSpace) : 
                                 space(aSpace), 
                                 start_point(aSpace.origin()), 
                                 end_point(aSpace.origin()),
                                 waypoints() { 
      waypoints.push_back(start_point);
      waypoints.push_back(end_point);
    };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     */
    point_to_point_path(const topology& aSpace, const point_type& aStart, const point_type& aEnd) :
                        space(aSpace), start_point(aEnd), end_point(aEnd), waypoints() {
      waypoints.push_back(start_point);
      waypoints.push_back(end_point);
    };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     */
    template <typename ForwardIter>
    point_to_point_path(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace) : 
                        space(aSpace), waypoints(aBegin, aEnd) {
      if(waypoints.size() > 0) {
	start_point = waypoints.front();
	end_point = waypoints.back();
      } else 
	throw invalid_path("Point-to-point path (empty list of waypoints)");
    };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& getSpace() const throw() { return space; };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      std::swap(lhs.start_point, rhs.start_point);
      std::swap(lhs.end_point, rhs.end_point);
      lhs.waypoints.swap(rhs.waypoints);
    };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(const point_type& a, const point_type& b) const {
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_a = get_waypoint_bounds(a, waypoints.begin());
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_b = get_waypoint_bounds(b, wpb_a.first);
      return travel_distance_impl(a, wpb_a, b, wpb_b);
    };
    
    /**
     * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
     * \param a The first waypoint-point-pair.
     * \param b The second waypoint-point-pair.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_a = get_waypoint_bounds(a.second, a.first);
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_b = get_waypoint_bounds(b.second, b.first);
      a.first = wpb_a.first; b.first = wpb_b.first;
      return travel_distance_impl(a.second, wpb_a, b.second, wpb_b);
    };
    
    /**
     * Computes the point that is a distance away from a point on the path.
     * \param a The point on the path.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    point_type move_away_from(const point_type& a, double d) const {
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_a = get_waypoint_bounds(a, waypoints.begin());
      point_type prev = a;
      const_waypoint_descriptor it = wpb_a.second;
      while((it != waypoints.end()) && (d > std::numeric_limits<double>::epsilon())) {
	double d1 = space.distance(prev, *it);
	if(d1 > d)
	  return space.move_position_toward(prev, d / d1, *(wpb_a.second));
	d -= d1; prev = *it; ++it;
      };
      return prev;
    };
    
    /**
     * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the path.
     * \param a The waypoint-point-pair on the path.
     * \param d The distance to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a distance away from the given waypoint-point-pair.
     */
    waypoint_pair move_away_from(const waypoint_pair& a, double d) const {
      std::pair<const_waypoint_descriptor, const_waypoint_descriptor> wpb_a = get_waypoint_bounds(a.second, a.first);
      const_waypoint_descriptor it_prev = wpb_a.first;
      point_type prev = a.second;
      const_waypoint_descriptor it = wpb_a.second;
      while((it != waypoints.end()) && (d > std::numeric_limits<double>::epsilon())) {
	double d1 = space.distance(prev, *it);
	if(d1 > d)
	  return std::make_pair(it_prev, space.move_position_toward(prev, d / d1, *(wpb_a.second)));
	d -= d1; prev = *it; it_prev = it; ++it;
      };
      return std::make_pair(it_prev,prev);
    };
    
};




};

};

#endif









