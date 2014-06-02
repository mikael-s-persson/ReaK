/**
 * \file discrete_point_path.hpp
 * 
 * This library provides an implementation of a path represented by discrete points within a topology.
 * The path is represented by a set of waypoints (presumably close to each other) and traveling along the 
 * path is restricted to hopping between discrete waypoints.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_DISCRETE_POINT_PATH_HPP
#define REAK_DISCRETE_POINT_PATH_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/ctrl/path_planning/spatial_path_concept.hpp>

#include "waypoint_container.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include <utility>

namespace ReaK {

namespace pp {
  

/**
 * This class implements a point-to-point path within a topology (interpolated by whatever method is used in 
 * the topology's move function). The path is represented by a set of waypoints and all intermediate points 
 * are computed with the topology's move_position_toward function based on the distance or fraction of travel. 
 * This class models the SpatialPathConcept and SequentialPathConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class discrete_point_path : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef discrete_point_path<Topology,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  private:
    
    typedef typename base_class_type::container_type container_type;
    
    const container_type& get_waypoints() const { return this->waypoints; };
    const distance_metric& get_dist() const { return this->dist; };
    const topology& get_space() const { return *(this->space); };
    
    
    double travel_distance_impl(const_waypoint_descriptor a, 
                                const_waypoint_descriptor b) const {
      // first assume that a is before b:
      double sum = 0.0;
      const_waypoint_descriptor a_cpy = a;
      while(a_cpy != b) {
        const_waypoint_descriptor a_next = a_cpy; ++a_next;
        if(a_next == this->waypoints.end()) {
          if(b == this->waypoints.end()) 
            return sum; // we're done.
          else
            return travel_distance_impl(b,a); // this means that b must be before a, so, reverse the computation.
        };
        sum += get_dist()(*a_cpy, *a_next, get_space());
        a_cpy = a_next;
      };
      return sum;
    };
    
    const_waypoint_descriptor move_away_from_impl(const_waypoint_descriptor a, double d) const {
      bool go_backwards = false;
      if(d < 0.0) {
        go_backwards = true;
        d = -d;
      };
      
      double sum = 0.0;
      while(sum < d) {
        const_waypoint_descriptor a_next = a; 
        if(go_backwards) {
          if(a == this->waypoints.begin())
            return a;
          --a_next;
          sum += get_dist()(*a_next,*a,get_space());
        } else {
          ++a_next;
          if(a_next == this->waypoints.end())
            return a;
          sum += get_dist()(*a,*a_next,get_space());
        };
        a = a_next;
      };
      return a;
    };
    
    const_waypoint_descriptor move_fraction_away_from_impl(const_waypoint_descriptor a, double d) const {
      bool go_backwards = false;
      if(d < 0.0) {
        go_backwards = true;
        d = -d;
      };
      
      double sum = 0.0;
      while(sum < d) {
        const_waypoint_descriptor a_next = a; 
        if(go_backwards) {
          if(a == this->waypoints.begin())
            return a;
          --a_next;
        } else {
          ++a_next;
          if(a_next == this->waypoints.end())
            return a;
        };
        sum += 1.0;
        a = a_next;
      };
      return a;
    };
    
  public:
    
    
    struct point_distance_iterator {
      const discrete_point_path* parent;
      const_waypoint_descriptor current_wpt;
      
      point_distance_iterator(const discrete_point_path* aParent, 
                              const const_waypoint_descriptor& aWPt) :
                              parent(aParent), current_wpt(aWPt) { };
      
      point_distance_iterator& operator+=(double rhs) {
        current_wpt = parent->move_away_from_impl(current_wpt, rhs);
        return *this;
      };
      
      friend 
      point_distance_iterator operator+(point_distance_iterator lhs, double rhs) {
        return (lhs += rhs);
      };
      
      friend 
      point_distance_iterator operator+(double lhs, point_distance_iterator rhs) {
        return (rhs += lhs);
      };
      
      friend 
      point_distance_iterator operator-(point_distance_iterator lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend
      point_distance_iterator& operator-=(point_distance_iterator& lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend 
      bool operator==(const point_distance_iterator& lhs, 
                      const point_distance_iterator& rhs) {
        return ((lhs.parent == rhs.parent) && (lhs.current_wpt == rhs.current_wpt));
      };
      
      friend 
      bool operator!=(const point_distance_iterator& lhs, 
                      const point_distance_iterator& rhs) {
        return !(lhs == rhs);
      };
      
      const point_type& operator*() const {
        return *current_wpt;
      };
      
    };
    
    
    struct point_fraction_iterator {
      const discrete_point_path* parent;
      const_waypoint_descriptor current_wpt;
      
      point_fraction_iterator(const discrete_point_path* aParent, 
                              const const_waypoint_descriptor& aWPt) :
                              parent(aParent), current_wpt(aWPt) { };
      
      point_fraction_iterator& operator+=(double rhs) {
        current_wpt = parent->move_fraction_away_from_impl(current_wpt, rhs);
        return *this;
      };
      
      friend 
      point_fraction_iterator operator+(point_fraction_iterator lhs, double rhs) {
        return (lhs += rhs);
      };
      
      friend 
      point_fraction_iterator operator+(double lhs, point_fraction_iterator rhs) {
        return (rhs += lhs);
      };
      
      friend 
      point_fraction_iterator operator-(point_fraction_iterator lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend
      point_fraction_iterator& operator-=(point_fraction_iterator& lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend 
      bool operator==(const point_fraction_iterator& lhs, 
                      const point_fraction_iterator& rhs) {
        return ((lhs.parent == rhs.parent) && (lhs.current_wpt == rhs.current_wpt));
      };
      
      friend 
      bool operator!=(const point_fraction_iterator& lhs, 
                      const point_fraction_iterator& rhs) {
        return !(lhs == rhs);
      };
      
      const point_type& operator*() const {
        return *current_wpt;
      };
      
    };
    
    
    
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit discrete_point_path(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), 
                                 const distance_metric& aDist = distance_metric()) : 
                                 base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    discrete_point_path(const shared_ptr<topology>& aSpace, 
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
    discrete_point_path(ForwardIter aBegin, ForwardIter aEnd, 
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
     * Returns the starting distance-iterator along the path.
     * \return The starting distance-iterator along the path.
     */
    point_distance_iterator begin_distance_travel() const {
      return point_distance_iterator(this, this->waypoints.begin());
    };
    
    /**
     * Returns the end distance-iterator along the path.
     * \return The end distance-iterator along the path.
     */
    point_distance_iterator end_distance_travel() const {
      const_waypoint_descriptor it = this->waypoints.end(); --it;
      return point_distance_iterator(this, it);
    };
    
    /**
     * Returns the starting fraction-iterator along the path.
     * \return The starting fraction-iterator along the path.
     */
    point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(this, this->waypoints.begin());
    };
    
    /**
     * Returns the end fraction-iterator along the path.
     * \return The end fraction-iterator along the path.
     */
    point_fraction_iterator end_fraction_travel() const {
      const_waypoint_descriptor it = this->waypoints.end(); --it;
      return point_fraction_iterator(this, it);
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
      return this->travel_distance_impl(wpb_a.first, wpb_b.first);
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
      return this->travel_distance_impl(wpb_a.first, wpb_b.first);
    };
    
    
    /**
     * Computes the point that is a distance away from a point on the path.
     * \param a The point on the path.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    point_type move_away_from(const point_type& a, double d) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      return *(move_away_from_impl(wpb_a.first, d));
    };
    
    /**
     * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the path.
     * \param a The waypoint-point-pair on the path.
     * \param dt The distance to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_away_from(const waypoint_pair& a, double d) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      const_waypoint_descriptor res = move_away_from_impl(wpb_a.first, d);
      return waypoint_pair(res, *res);
    };
       
    /**
     * Computes the point that is on the path at the given distance from the start.
     * \param d The distance (from start) at which the point is sought.
     * \return The point that is on the path at the given distance from the start.
     */
    point_type get_point_at_distance(double d) const {
      return *(move_away_from_impl(this->waypoints.begin(), d));
    };
    
    /**
     * Computes the waypoint-point pair that is on the path at the given distance from the start.
     * \param t The distance from the start at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given distance from the start.
     */
    waypoint_pair get_waypoint_at_distance(double d) const {
      const_waypoint_descriptor start = move_away_from_impl(this->waypoints.begin(), d);
      return waypoint_pair(start, *start);
    };
    
    /**
     * Returns the total length of the path.
     * \return The total length of the path.
     */
    double get_total_length() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return travel_distance_impl(start, end);
    };
    
    /**
     * Returns the starting point of the path.
     * \return The starting point of the path.
     */
    point_type get_start_point() const {
      return *(this->waypoints.begin());
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440013,1,"discrete_point_path",base_class_type)
    
};



};

};

#endif









