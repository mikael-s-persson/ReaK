/**
 * \file point_to_point_path.hpp
 * 
 * This library provides an implementation of a linear-interpolated path within a topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation based on the distance.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_POINT_TO_POINT_PATH_HPP
#define REAK_POINT_TO_POINT_PATH_HPP

#include "base/defs.hpp"

#include "path_planning/spatial_path_concept.hpp"

#include "waypoint_container.hpp"

#include "topologies/basic_distance_metrics.hpp"

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
 * This class implements a linear-interpolated path within a topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation based on the distance. 
 * This class models the SpatialPathConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the TemporalDistMetricConcept.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class point_to_point_path : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef point_to_point_path<Topology,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  private:
    
    double travel_distance_impl(const point_type& a, const const_waypoint_bounds& wpb_a, 
                                const point_type& b, const const_waypoint_bounds& wpb_b) const {
      if(std::distance(wpb_a.first, this->waypoints.end()) < 
         std::distance(wpb_b.first, this->waypoints.end())) 
        return travel_distance_impl(b,wpb_b,a,wpb_a);
      
      double sum = 0;
      if(((wpb_a.first == wpb_b.first) && (wpb_a.second == wpb_b.second)) ||
         (wpb_a.second == wpb_b.first)) {
        //this means that a and b are in the same segment.
        return this->dist(a,b,*(this->space));
      };
      
      sum += this->dist(a, *(wpb_a.second), *(this->space));
      
      const_waypoint_descriptor it = wpb_a.second;
      const_waypoint_descriptor it_prev = it;
      while(++it != wpb_b.first)
        sum += this->dist(*it_prev++, *it, *(this->space));
      
      sum += this->dist(*it_prev, *(wpb_b.first), *(this->space));
      
      sum += this->dist(*(wpb_b.first), b, *(this->space));
      
      return sum;
    };
    
    waypoint_pair get_point_at_distance_impl(double d, const const_waypoint_bounds& wpb_a) const {
      const_waypoint_descriptor it_prev = wpb_a.first;
      const_waypoint_descriptor it = it_prev; ++it;
      
      if(it == this->waypoints.end())
        return waypoint_pair(it_prev,*it_prev);
      
      double total_d = this->dist(*it_prev, *it, *(this->space));
      
      return waypoint_pair(it_prev, this->space->move_position_toward(*it_prev, d / total_d, *it));
    };
    
    waypoint_pair move_away_from_impl(point_type a, const_waypoint_bounds wpb_a, double d) const {
      while(true) {
        double d_to_second = this->dist(a, *(wpb_a.second), *(this->space));
        if(d_to_second >= d)
          return waypoint_pair(wpb_a.first, this->space->move_position_toward(a, d / d_to_second, *(wpb_a.second)));
        const_waypoint_descriptor it = wpb_a.second; ++it;
        if(it == this->waypoints.end())
          return waypoint_pair(wpb_a.first, *(wpb_a.second));
        d -= d_to_second;
        a = *(wpb_a.second);
        wpb_a.first = wpb_a.second;
        wpb_a.second = it;
      };
    };
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit point_to_point_path(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), 
                                 const distance_metric& aDist = distance_metric()) : 
                                 base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    point_to_point_path(const shared_ptr<topology>& aSpace, 
                        const point_type& aStart, 
                        const point_type& aEnd, 
                        const distance_metric& aDist = distance_metric()) :
                        base_class_type(aSpace, aStart, aEnd, aDist) { };
    
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the on-past-last point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    point_to_point_path(ForwardIter aBegin, ForwardIter aEnd, 
                        const shared_ptr<topology>& aSpace, 
                        const distance_metric& aDist = distance_metric()) : 
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
     * Computes the point that is a distance away from a point on the path.
     * \param a The point on the path.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    point_type move_away_from(const point_type& a, double d) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      return move_away_from_impl(a,wpb_a,d).second;
    };
    
    /**
     * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the path.
     * \param a The waypoint-point-pair on the path.
     * \param dt The distance to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_away_from(const waypoint_pair& a, double d) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      return move_away_from_impl(a.second,wpb_a,d);
    };
       
    /**
     * Computes the point that is on the path at the given distance from the start.
     * \param d The distance (from start) at which the point is sought.
     * \return The point that is on the path at the given distance from the start.
     */
    point_type get_point_at_distance(double d) const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return move_away_from_impl(*start, const_waypoint_bounds(start,start), d).second;
    };
    
    /**
     * Computes the waypoint-point pair that is on the path at the given distance from the start.
     * \param t The distance from the start at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given distance from the start.
     */
    waypoint_pair get_waypoint_at_distance(double d) const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return move_away_from_impl(*start, const_waypoint_bounds(start,start), d);
    };
    
    /**
     * Returns the total length of the path.
     * \return The total length of the path.
     */
    double get_total_length() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return travel_distance_impl(*start, const_waypoint_bounds(start,start),
                                  *end, const_waypoint_bounds(end,end));
    };
    
    /**
     * Returns the starting point of the path.
     * \return The starting point of the path.
     */
    point_type get_start_point() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return *start;
    };
    
    /**
     * Returns the starting waypoint-point-pair of the path.
     * \return The starting waypoint-point-pair of the path.
     */
    waypoint_pair get_start_waypoint() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return waypoint_pair(start,*start);
    };
    
    /**
     * Returns the end point of the path.
     * \return The end point of the path.
     */
    point_type get_end_point() const {
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return *end;
    };
    
    /**
     * Returns the end waypoint-point-pair of the path.
     * \return The end waypoint-point-pair of the path.
     */
    waypoint_pair get_end_waypoint() const {
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return waypoint_pair(end,*end);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC244000B,1,"point_to_point_path",base_class_type)
    
};



};

};

#endif









