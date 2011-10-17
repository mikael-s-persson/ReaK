/**
 * \file cubic_hermite_interp.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a cubic Hermite interpolation (cubic hermite spline, or cspline).
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

#ifndef REAK_CUBIC_HERMITE_INTERP_HPP
#define REAK_CUBIC_HERMITE_INTERP_HPP

#include "path_planning/spatial_trajectory_concept.hpp"

#include "path_planning/differentiable_space_concept.hpp"

#include "waypoint_container.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

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
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class cubic_hermite_interp : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 1, typename temporal_topology_traits<Topology>::time_topology >));
    
    typedef cubic_hermite_interp<Topology,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  private:
    
    cubic_hermite_interp(const cubic_hermite_interp<Topology,DistanceMetric>&); //non-copyable.
    cubic_hermite_interp<Topology,DistanceMetric>& operator=(const cubic_hermite_interp<Topology,DistanceMetric>&);
    
    double travel_distance_impl(const point_type& a, const const_waypoint_bounds& wpb_a, 
				const point_type& b, const const_waypoint_bounds& wpb_b) const {
      if((wpb_a.first == wpb_a.second) || (wpb_a.first == wpb_b.first))
	return dist(a,b,space); //this means that a is at the end of the path, thus, the "path" goes directly towards b.
      double sum = dist(a, *(wpb_a.second), space);
      const_waypoint_descriptor it = wpb_a.second;
      while(it != wpb_b.first) {
	const_waypoint_descriptor it_prev = it;
	sum += dist(*it_prev, *(++it),space);
      };
      sum += dist(*(wpb_b.first), b, space);
      return sum;
    };
    
    point_type compute_interpolation_impl(const point_type& a, const point_type& b, double t) const {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::get;
#endif
      
      double t_factor = b.time - a.time;
      if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
	throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
      double t_normal = (t - a.time) / (b.time - a.time);
      
      point_type result = this->space.move_position_toward(a,t_normal,b); //NOTE this is just to get other higher-order values linearly interpolated.
      
      typedef typename temporal_topology_traits<topology>::space_topology SpaceType;
      const SpaceType& pt_space = this->space.get_space_topology(); 
      const typename temporal_topology_traits<topology>::time_topology& t_space = this->space.get_time_topology(); 
      get<0>(result.pt) = pt_space.get_space<0>(t_space).adjust(get<0>(a.pt), 
        ((3.0 - 2.0 * t_normal) * t_normal * t_normal) * pt_space.get_space<0>(t_space).difference(get<0>(b.pt),get<0>(a.pt))
        + (t_normal * (1.0 - t_normal * (2.0 + t_normal))) * pt_space.descend_to_space<0>(get<1>(a.pt),t_factor)
        + (t_normal * t_normal * (t_normal - 1.0)) * pt_space.descend_to_space<0>(get<1>(b.pt),t_factor) );
      
      // TODO: This is a bit weird, I have to think about it some more.
      get<1>(result.pt) = pt_space.get_space<1>(t_space).adjust(
        pt_space.lift_to_space<1>(((6.0 * t_normal - 6.0) * t_normal) * pt_space.get_space<0>(t_space).difference(get<0>(b.pt),get<0>(a.pt)), t_factor),
        (t_normal * (1.0 - t_normal * (2.0 + t_normal))) * pt_space.descend_to_space<0>(get<1>(a.pt),t_factor)
        + (t_normal * t_normal * (t_normal - 1.0)) * pt_space.descend_to_space<0>(get<1>(b.pt),t_factor) );
      
    };
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit cubic_hermite_interp(const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                  base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    cubic_hermite_interp(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                         base_class_type(aSpace, aStart, aEnd, aDist) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    cubic_hermite_interp(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                         base_class_type(aBegin, aEnd, aSpace, aDist) { };
    
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
      BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      const point_type* prev = &a;
      const_waypoint_descriptor it = wpb_a.second;
      while((it != this->waypoints.end()) && (dt > std::numeric_limits<double>::epsilon())) {
	double d1 = it->time - prev->time;
	if(d1 > dt)
	  return this->space.move_position_toward(*prev, dt / d1, *(wpb_a.second));
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
      BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      const_waypoint_descriptor it_prev = wpb_a.first;
      const point_type* prev = &(a.second);
      const_waypoint_descriptor it = wpb_a.second;
      while((it != this->waypoints.end()) && (dt > std::numeric_limits<double>::epsilon())) {
	double d1 = it->time - prev->time;
	if(d1 > dt)
	  return std::make_pair(it_prev, this->space.move_position_toward(*prev, dt / d1, *(wpb_a.second)));
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
      BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
      const_waypoint_descriptor start = this->waypoints.begin();
      point_type p = *start;
      p.time = t;
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(p, start);
      if(wpb_p.first == wpb_p.second)
	return *(wpb_p.first);
      return this->space.move_position_toward(*(wpb_p.first),(t - wpb_p.first->time) / (wpb_p.second->time - wpb_p.first->time),*(wpb_p.second));
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given time.
     * \param t The time at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given time.
     */
    waypoint_pair get_waypoint_at_time(double t) const {
      BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
      const_waypoint_descriptor start = this->waypoints.begin();
      point_type p = *start;
      p.time = t;
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(p, start);
      if(wpb_p.first == wpb_p.second)
	return std::make_pair(wpb_p.first, *(wpb_p.first));
      return std::make_pair(wpb_p.first, this->space.move_position_toward(*(wpb_p.first),(t - wpb_p.first->time) / (wpb_p.second->time - wpb_p.first->time),*(wpb_p.second)));
    };
    
    
};



};

};

#endif









