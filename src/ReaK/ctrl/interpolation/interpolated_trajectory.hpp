/**
 * \file interpolated_trajectory.hpp
 * 
 * This library provides an implementation of an interpolated trajectory within a temporal topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with an interpolation functor provided as a template argument.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_INTERPOLATED_TRAJECTORY_HPP
#define REAK_INTERPOLATED_TRAJECTORY_HPP

#include "path_planning/spatial_trajectory_concept.hpp"

#include "waypoint_container.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include <lin_alg/mat_num_exceptions.hpp>

namespace ReaK {

namespace pp {
  

/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept.
 * \tparam Interpolator The interpolation functor type on which the points and the path can reside, should be callable as p_int = interpolator(p1, p2, time, temporal_space).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename Interpolator, typename DistanceMetric = default_distance_metric>
class interpolated_trajectory : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef interpolated_trajectory<Topology,Interpolator,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef Interpolator interpolator_type;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  private:
    
    interpolator_type interp;
    
    interpolated_trajectory(const cubic_hermite_interp<Topology,Interpolator,DistanceMetric>&); //non-copyable.
    interpolated_trajectory<Topology,Interpolator,DistanceMetric>& operator=(const interpolated_trajectory<Topology,Interpolator,DistanceMetric>&);
    
    double travel_distance_impl(const point_type& a, const const_waypoint_bounds& wpb_a, 
				const point_type& b, const const_waypoint_bounds& wpb_b) const {
      if((wpb_a.first == wpb_a.second) || (wpb_a.first == wpb_b.first))
	return this->dist(a,b,this->space); //this means that a is at the end of the path, thus, the "path" goes directly towards b.
      double sum = this->dist(a, *(wpb_a.second), this->space);
      const_waypoint_descriptor it = wpb_a.second;
      while(it != wpb_b.first) {
	const_waypoint_descriptor it_prev = it;
	sum += this->dist(*it_prev, *(++it),this->space);
      };
      sum += this->dist(*(wpb_b.first), b, this->space);
      return sum;
    };
    
  public:
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the trajectory is.
     * \param aDist The distance metric functor that the trajectory should use.
     * \param aInterp The interpolator functor that the trajectory should use.
     */
    explicit cubic_hermite_interp(const topology& aSpace, const distance_metric& aDist = distance_metric(), const interpolator_type& aInterp = interpolator_type()) : 
                                  base_class_type(aSpace, aDist), interp(aInterp) { };
    
    /**
     * Constructs the trajectory from a space, the start and end points.
     * \param aSpace The space on which the trajectory is.
     * \param aStart The start point of the trajectory.
     * \param aEnd The end-point of the trajectory.
     * \param aDist The distance metric functor that the trajectory should use.
     * \param aInterp The interpolator functor that the trajectory should use.
     */
    cubic_hermite_interp(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric(), const interpolator_type& aInterp = interpolator_type()) :
                         base_class_type(aSpace, aStart, aEnd, aDist), interp(aInterp) { };
			
    /**
     * Constructs the trajectory from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the trajectory with.
     * \param aBegin An iterator to the first point of the trajectory.
     * \param aEnd An iterator to the on-past-last point of the trajectory.
     * \param aSpace The space on which the trajectory is.
     * \param aDist The distance metric functor that the trajectory should use.
     * \param aInterp The interpolator functor that the trajectory should use.
     */
    template <typename ForwardIter>
    cubic_hermite_interp(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric(), const interpolator_type& aInterp = interpolator_type()) : 
                         base_class_type(aBegin, aEnd, aSpace, aDist), interp(aInterp) { };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(static_cast<base_class_type&>(lhs),static_cast<base_class_type&>(rhs));
    };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(const point_type& a, const point_type& b) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b, wpb_a.first);
      return this->travel_distance_impl(a, wpb_a, b, wpb_b);
    };
    
    /**
     * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
     * \param a The first waypoint-point-pair.
     * \param b The second waypoint-point-pair.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b.second, b.first);
      a.first = wpb_a.first; b.first = wpb_b.first;
      return this->travel_distance_impl(a.second, wpb_a, b.second, wpb_b);
    };
    
    
    /**
     * Computes the point that is a time-difference away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param dt The time to move away from the point.
     * \return The point that is a time away from the given point.
     */
    point_type move_time_diff_from(const point_type& a, double dt) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      const point_type* prev = &a;
      const_waypoint_descriptor it = wpb_a.second;
      while((it != this->waypoints.end()) && (dt > std::numeric_limits<double>::epsilon())) {
	double d1 = it->time - prev->time;
	if(d1 > dt)
	  return interp(*prev, *(wpb_a.second), prev->time + dt, this->space);
	dt -= d1; prev = &(*it); ++it;
      };
      return *prev;
    };
    
    /**
     * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
     * \param a The waypoint-point-pair on the trajectory.
     * \param dt The time to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_time_diff_from(const waypoint_pair& a, double dt) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      const_waypoint_descriptor it_prev = wpb_a.first;
      const point_type* prev = &(a.second);
      const_waypoint_descriptor it = wpb_a.second;
      while((it != this->waypoints.end()) && (dt > std::numeric_limits<double>::epsilon())) {
	double d1 = it->time - prev->time;
	if(d1 > dt)
	  return std::make_pair(it_prev, interp(*prev, *(wpb_a.second), prev->time + dt, this->space));
	d -= d1; prev = &(*it); it_prev = it; ++it;
      };
      return std::make_pair(it_prev,*prev);
    };
       
    /**
     * Computes the point that is on the trajectory at the given time.
     * \param t The time at which the point is sought.
     * \return The point that is on the trajectory at the given time.
     */
    point_type get_point_at_time(double t) const {
      const_waypoint_descriptor start = this->waypoints.begin();
      point_type p = *start;
      p.time = t;
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(p, start);
      if(wpb_p.first == wpb_p.second)
	return *(wpb_p.first);
      return interp(*(wpb_p.first), *(wpb_p.second), t, this->space);
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given time.
     * \param t The time at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given time.
     */
    waypoint_pair get_waypoint_at_time(double t) const {
      const_waypoint_descriptor start = this->waypoints.begin();
      point_type p = *start;
      p.time = t;
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(p, start);
      if(wpb_p.first == wpb_p.second)
	return std::make_pair(wpb_p.first, *(wpb_p.first));
      return std::make_pair(wpb_p.first, interp(*(wpb_p.first), *(wpb_p.second), t, this->space));
    };
    
    
};



};

};

#endif









