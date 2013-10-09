/**
 * \file point_to_point_trajectory.hpp
 * 
 * This library provides an implementation of a point-to-point trajectory within a temporal topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with the topology's functions.
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

#ifndef REAK_POINT_TO_POINT_TRAJECTORY_HPP
#define REAK_POINT_TO_POINT_TRAJECTORY_HPP

#include "base/defs.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"

#include "waypoint_container.hpp"

#include "topologies/basic_distance_metrics.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include <utility>
#include <limits>
#include <cmath>

namespace ReaK {

namespace pp {
  

/**
 * This class implements a point-to-point trajectory within a temporal topology (interpolated by whatever method is used in 
 * the topology's move function). The path is represented by a set of waypoints and all intermediate points 
 * are computed with the topology's move_position_toward function based on the distance or fraction of travel. 
 * This class models the SpatialTrajectoryConcept and SequentialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class point_to_point_trajectory : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef point_to_point_trajectory<Topology,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef typename base_class_type::waypoint_pair waypoint_pair;
    
    
    
  private:
    
    typedef typename base_class_type::container_type container_type;
    
    const container_type& get_waypoints() const { return this->waypoints; };
    const distance_metric& get_dist() const { return this->dist; };
    const topology& get_space() const { return *(this->space); };
    
    
    double travel_distance_impl(const point_type& a, const const_waypoint_bounds& wpb_a, 
                                const point_type& b, const const_waypoint_bounds& wpb_b) const {
      if(a.time > b.time)
        return travel_distance_impl(b,wpb_b,a,wpb_a);
      
      double sum = 0;
      if((wpb_a.first == wpb_b.first) && (wpb_a.first == wpb_b.first)) {
        //this means that a and b are in the same segment.
        return this->dist(a, b, *(this->space));
      };
      
      sum += this->dist(a, wpb_a.second->second, *(this->space));
      
      const_waypoint_descriptor it = wpb_a.second;
      const_waypoint_descriptor it_prev = it;
      while(it++ != wpb_b.first)
        sum += this->dist((it_prev++)->second, it->second, *(this->space));
      
      sum += this->dist(wpb_b.first->second, b, *(this->space));
      
      return sum;
    };
    
    waypoint_pair get_point_at_time_impl(double t, const const_waypoint_bounds& wpb_a) const {
      using std::fabs;
      
      const_waypoint_descriptor it_prev = wpb_a.first;
      const_waypoint_descriptor it = it_prev; ++it;
      
      if(it == this->waypoints.end()) {
        if(it_prev == this->waypoints.begin()) {
          waypoint_pair result(it_prev, it_prev->second);
          result.second.time = t;
          return result;
        };
        it = it_prev; --it_prev;
      };
      
      double t_total = it->first - it_prev->first;
      double dt = t - it->first;
      if( t_total < 1e-6 * fabs(it->first) ) {
        waypoint_pair result(it_prev, it_prev->second);
        result.second.time = t;
        return result;
      };
      
      return waypoint_pair(it_prev, this->space->move_position_toward(it_prev->second, dt / t_total, it->second));
    };
    
    
    waypoint_pair move_time_diff_from_impl(const point_type& a, const const_waypoint_bounds& wpb_a, double dt) const {
      if((dt > 0.0) && (wpb_a.second->first > a.time + dt))
        return waypoint_pair(wpb_a.first, this->space->move_position_toward(a, dt / (wpb_a.second->first - a.time), wpb_a.second->second));
      else if((dt <= 0.0) && (wpb_a.first->first < a.time + dt))
        return waypoint_pair(wpb_a.first, this->space->move_position_toward(wpb_a.first->second, -dt / (a.time - wpb_a.first->first), a));
      
      point_type result = a;
      result.time += dt;
      return get_point_at_time_impl(result.time,this->get_waypoint_bounds(result, wpb_a.first));
    };
    
  public:
    
    
    struct point_time_iterator {
      const point_to_point_trajectory* parent;
      const_waypoint_bounds current_wpbound;
      point_type current_pt;
      double current_dt;
      
      point_time_iterator(const point_to_point_trajectory* aParent,
                          const const_waypoint_bounds& aWPB, 
                          const point_type& aCurrentPt) :
                          parent(aParent),
                          current_wpbound(aWPB), 
                          current_pt(aCurrentPt),
                          current_dt(current_pt.time - current_wpbound.first->first) {
        if(current_wpbound.first == current_wpbound.second)
          current_dt = 0.0;
      };
      
      point_time_iterator(const point_to_point_trajectory* aParent,
                          const const_waypoint_bounds& aWPB) :
                          parent(aParent),
                          current_wpbound(aWPB), 
                          current_pt(current_wpbound.first->second), 
                          current_dt(0.0) { };
                              
#ifdef RK_ENABLE_CXX11_FEATURES
      point_time_iterator(const point_to_point_trajectory* aParent,
                          const_waypoint_bounds&& aWPB, 
                          point_type&& aCurrentPt) :
                          parent(aParent),
                          current_wpbound(std::move(aWPB)), 
                          current_pt(std::move(aCurrentPt)), 
                          current_dt(current_pt.time - current_wpbound.first->first) {
        if(current_wpbound.first == current_wpbound.second)
          current_dt = 0.0;
      };
      
      point_time_iterator(const point_to_point_trajectory* aParent,
                          const_waypoint_bounds&& aWPB) :
                          parent(aParent),
                          current_wpbound(std::move(aWPB)), 
                          current_pt(current_wpbound.first->second), 
                          current_dt(0.0) { };
#endif
      
      point_time_iterator& operator+=(double rhs) {
        if(current_wpbound.first == current_wpbound.second) {
          // then we are at some edge of the set of waypoints.
          if(current_wpbound.first == parent->get_waypoints().begin()) {
            // we are at the beginning.
            if(rhs < 0.0)
              return *this;  // cannot move further back.
            ++(current_wpbound.second);
            if(current_wpbound.second == parent->get_waypoints().end()) {
              current_wpbound.second = current_wpbound.first;
              return *this;  // nothing beyond the first point.
            };
          } else {
            // we are at the end.
            if(rhs > 0.0) 
              return *this;  // cannot move further forward.
            --(current_wpbound.first);
            current_dt = current_wpbound.second->first - current_wpbound.first->first;
          };
        };
        
        if( rhs < 0.0) {
          while(current_dt + rhs < 0.0) {  // while spilling over before the current bounds.
            if(current_wpbound.first == parent->get_waypoints().begin()) {
              current_pt = current_wpbound.first->second;
              current_wpbound.second = current_wpbound.first;
              current_dt = 0.0;
              return *this;
            };
            --(current_wpbound.first);
            --(current_wpbound.second);
            rhs += current_dt;
            current_dt = current_wpbound.second->first - current_wpbound.first->first;
          };
        } else {
          while(current_wpbound.second->first - current_wpbound.first->first <= current_dt + rhs) {  // while spilling over beyond the current bounds.
            rhs -= current_wpbound.second->first - current_wpbound.first->first - current_dt;
            ++(current_wpbound.first);
            ++(current_wpbound.second);
            current_dt = 0.0;
            if(current_wpbound.second == parent->get_waypoints().end()) {
              current_pt = current_wpbound.first->second;
              current_wpbound.second = current_wpbound.first;
              current_dt = 0.0;
              return *this;
            };
          };
        };
        
        double frac = (current_dt + rhs) / (current_wpbound.second->first - current_wpbound.first->first);
        current_pt = parent->get_space().move_position_toward(current_wpbound.first->second, frac, current_wpbound.second->second);
        current_dt += rhs;
        return *this;
      };
      
      friend 
      point_time_iterator operator+(point_time_iterator lhs, double rhs) {
        return (lhs += rhs);
      };
      
      friend 
      point_time_iterator operator+(double lhs, point_time_iterator rhs) {
        return (rhs += lhs);
      };
      
      friend 
      point_time_iterator operator-(point_time_iterator lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend
      point_time_iterator& operator-=(point_time_iterator& lhs, double rhs) {
        return (lhs += -rhs);
      };
      
      friend 
      bool operator==(const point_time_iterator& lhs, 
                      const point_time_iterator& rhs) {
        using std::fabs;
        if( ( lhs.parent == rhs.parent ) && 
            ( lhs.current_wpbound.first == rhs.current_wpbound.first ) &&
            ( lhs.current_wpbound.second == rhs.current_wpbound.second ) &&
            ( fabs(lhs.current_dt - rhs.current_dt) <= 1e-6 * (lhs.current_wpbound.second->first - lhs.current_wpbound.first->first) ) )
          return true;
        else 
          return false;
      };
      
      friend 
      bool operator!=(const point_time_iterator& lhs, 
                      const point_time_iterator& rhs) {
        return !(lhs == rhs);
      };
      
      const point_type& operator*() const {
        return current_pt;
      };
      
    };
    
    
    struct point_fraction_iterator {
      const point_to_point_trajectory* parent;
      const_waypoint_bounds current_wpbound;
      point_type current_pt;
      double current_fraction;
      
      point_fraction_iterator(const point_to_point_trajectory* aParent,
                              const const_waypoint_bounds& aWPB, 
                              const point_type& aCurrentPt) :
                              parent(aParent),
                              current_wpbound(aWPB), 
                              current_pt(aCurrentPt),
                              current_fraction((current_pt.time - current_wpbound.first->first) / (current_wpbound.second->first - current_wpbound.first->first)) {
        if(current_wpbound.first == current_wpbound.second)
          current_fraction = 0.0;
      };
      
      point_fraction_iterator(const point_to_point_trajectory* aParent,
                              const const_waypoint_bounds& aWPB) :
                              parent(aParent),
                              current_wpbound(aWPB), 
                              current_pt(current_wpbound.first->second),
                              current_fraction(0.0) { };
                              
#ifdef RK_ENABLE_CXX11_FEATURES
      point_fraction_iterator(const point_to_point_trajectory* aParent,
                              const_waypoint_bounds&& aWPB, 
                              point_type&& aCurrentPt) :
                              parent(aParent),
                              current_wpbound(std::move(aWPB)), 
                              current_pt(std::move(aCurrentPt)),
                              current_fraction((current_pt.time - current_wpbound.first->first) / (current_wpbound.second->first - current_wpbound.first->first)) {
        if(current_wpbound.first == current_wpbound.second)
          current_fraction = 0.0;
      };
      
      point_fraction_iterator(const point_to_point_trajectory* aParent,
                              const_waypoint_bounds&& aWPB) :
                              parent(aParent),
                              current_wpbound(std::move(aWPB)), 
                              current_pt(current_wpbound.first->second),
                              current_fraction(0.0) { };
#endif
      
      point_fraction_iterator& operator+=(double rhs) {
        if(current_wpbound.first == current_wpbound.second) {
          // then we are at some edge of the set of waypoints.
          if(current_wpbound.first == parent->get_waypoints().begin()) {
            // we are at the beginning.
            if(rhs < 0.0)
              return *this;  // cannot move further back.
            ++(current_wpbound.second);
            if(current_wpbound.second == parent->get_waypoints().end()) {
              current_wpbound.second = current_wpbound.first;
              return *this;  // nothing beyond the first point.
            };
          } else {
            // we are at the end.
            if(rhs > 0.0) 
              return *this;  // cannot move further forward.
            --(current_wpbound.first);
            current_fraction = 1.0;
          };
        };
        
        if(rhs < 0.0) {
          while(current_fraction + rhs < 0.0) {  // while spilling over before the current bounds.
            if(current_wpbound.first == parent->get_waypoints().begin()) {
              current_pt = current_wpbound.first->second;
              current_wpbound.second = current_wpbound.first;
              current_fraction = 0.0;
              return *this;
            };
            --(current_wpbound.first);
            --(current_wpbound.second);
            rhs += current_fraction;
            current_fraction = 1.0;
          };
        } else {
          while(1.0 <= current_fraction + rhs) {  // while spilling over beyond the current bounds.
            ++(current_wpbound.first);
            ++(current_wpbound.second);
            rhs -= 1.0 - current_fraction;
            current_fraction = 0.0;
            if(current_wpbound.second == parent->get_waypoints().end()) {
              current_pt = current_wpbound.first->second;
              current_wpbound.second = current_wpbound.first;
              current_fraction = 0.0;
              return *this;
            };
          };
        };
        
        current_pt = parent->get_space().move_position_toward(
          current_wpbound.first->second, current_fraction + rhs, current_wpbound.second->second);
        current_fraction += rhs;
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
        using std::fabs;
        if((lhs.parent == rhs.parent) && 
           (lhs.current_wpbound.first == rhs.current_wpbound.first) &&
           (lhs.current_wpbound.second == rhs.current_wpbound.second) &&
           (fabs(lhs.current_fraction - rhs.current_fraction) <= 1e-6))
          return true;
        else 
          return false;
      };
      
      friend 
      bool operator!=(const point_fraction_iterator& lhs, 
                      const point_fraction_iterator& rhs) {
        return !(lhs == rhs);
      };
      
      const point_type& operator*() const {
        return current_pt;
      };
      
    };
    
    
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the trajectory is.
     * \param aDist The distance metric functor that the trajectory should use.
     */
    explicit point_to_point_trajectory(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), 
                                       const distance_metric& aDist = distance_metric()) : 
                                       base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the trajectory from a space, the start and end points.
     * \param aSpace The space on which the trajectory is.
     * \param aStart The start point of the trajectory.
     * \param aEnd The end-point of the trajectory.
     * \param aDist The distance metric functor that the trajectory should use.
     */
    point_to_point_trajectory(const shared_ptr<topology>& aSpace, 
                              const point_type& aStart, 
                              const point_type& aEnd, 
                              const distance_metric& aDist = distance_metric()) :
                              base_class_type(aSpace, aStart, aEnd, aDist) { };
    
    /**
     * Constructs the trajectory from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the trajectory with.
     * \param aBegin An iterator to the first point of the trajectory.
     * \param aEnd An iterator to the on-past-last point of the trajectory.
     * \param aSpace The space on which the trajectory is.
     * \param aDist The distance metric functor that the trajectory should use.
     */
    template <typename ForwardIter>
    point_to_point_trajectory(ForwardIter aBegin, ForwardIter aEnd, 
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
     * Returns the starting time-iterator along the trajectory.
     * \return The starting time-iterator along the trajectory.
     */
    point_time_iterator begin_time_travel() const {
      return point_time_iterator(this, const_waypoint_bounds(this->waypoints.begin(),this->waypoints.begin()));
    };
    
    /**
     * Returns the end time-iterator along the trajectory.
     * \return The end time-iterator along the trajectory.
     */
    point_time_iterator end_time_travel() const {
      const_waypoint_descriptor it = this->waypoints.end(); --it;
      return point_time_iterator(this, const_waypoint_bounds(it,it));
    };
    
    /**
     * Returns the starting fraction-iterator along the path.
     * \return The starting fraction-iterator along the path.
     */
    point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(this, const_waypoint_bounds(this->waypoints.begin(),this->waypoints.begin()));
    };
    
    /**
     * Returns the end fraction-iterator along the path.
     * \return The end fraction-iterator along the path.
     */
    point_fraction_iterator end_fraction_travel() const {
      const_waypoint_descriptor it = this->waypoints.end(); --it;
      return point_fraction_iterator(this, const_waypoint_bounds(it,it));
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
     * Returns the total length of the path.
     * \return The total length of the path.
     */
    double get_total_length() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return travel_distance_impl(start->second, const_waypoint_bounds(start,start),
                                  end->second, const_waypoint_bounds(end,end));
    };
    
    
    /**
     * Computes the point that is a distance away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    point_type move_time_diff_from(point_type a, double d) const {
      a.time += d;
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a, this->waypoints.begin());
      return move_time_diff_from_impl(wpb_a.first, wpb_a, a.time - wpb_a.first->first)->second;
    };
    
    /**
     * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the trajectory.
     * \param a The waypoint-point-pair on the trajectory.
     * \param dt The distance to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_time_diff_from(waypoint_pair a, double d) const {
      a.second.time += d;
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
      const_waypoint_descriptor res = move_time_diff_from_impl(wpb_a.first, wpb_a, a.time - wpb_a.first->first);
      return waypoint_pair(res, res->second);
    };
    
    /**
     * Computes the point that is on the trajectory at the given time from the start.
     * \param d The time (from start) at which the point is sought.
     * \return The point that is on the trajectory at the given time from the start.
     */
    point_type get_point_at_time(double d) const {
      return move_time_diff_from(waypoint_pair(this->waypoints.begin(), this->waypoints.begin()->second), d).second;
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given distance from the start.
     * \param t The distance from the start at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given distance from the start.
     */
    waypoint_pair get_waypoint_at_time(double d) const {
      return move_time_diff_from(waypoint_pair(this->waypoints.begin(), this->waypoints.begin()->second), d);
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440017,1,"point_to_point_trajectory",base_class_type)
    
};



};

};

#endif









