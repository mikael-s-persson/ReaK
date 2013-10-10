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

#include "base/defs.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"
#include "path_planning/interpolator_concept.hpp"

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
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept.
 * \tparam InterpolatorFactory The interpolation factory type which can create interpolators for on the given topology, should model the InterpolatorFactoryConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the TemporalDistMetricConcept.
 */
template <typename Topology, typename InterpolatorFactory, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class interpolated_trajectory : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((InterpolatorFactoryConcept<InterpolatorFactory,Topology,DistanceMetric>));
    
    typedef interpolated_trajectory<Topology,InterpolatorFactory,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef InterpolatorFactory interpolator_factory_type;
    typedef typename interpolator_factory_traits<interpolator_factory_type>::interpolator_type interpolator_type;
    
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef typename base_class_type::waypoint_pair waypoint_pair;
    
    typedef std::map< const point_type*, interpolator_type > interpolator_map_type;
    
  private:
    
    interpolator_factory_type interp_fact;
    mutable interpolator_map_type interp_segments; //this is mutable for JIT construction of it.
    
    virtual double travel_distance_impl(const point_type& a, const point_type& b) const {
      if(a.time > b.time)
        return travel_distance_impl(b, a);
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.time);
      const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b.time);
      
      double sum = 0;
      if((wpb_a.first == wpb_b.first) && (wpb_a.first == wpb_b.first)) {
        //this means that a and b are in the same segment.
        interpolator_type seg = interp_fact.create_interpolator(&a,&b);
        return seg.travel_distance_to(b, this->dist);
      };
      if(wpb_a.first == wpb_a.second) {
        //this means that a is before the first point in the trajectory.
        interpolator_type seg = interp_fact.create_interpolator(&a,&(wpb_a.second->second));
        sum += seg.travel_distance_from(a, this->dist);
      } else {
        typename interpolator_map_type::iterator it_int = interp_segments.find(&(wpb_a.first->second));
        if(it_int == interp_segments.end()) {
          interp_segments[&(wpb_a.first->second)] = interp_fact.create_interpolator(&(wpb_a.first->second),&(wpb_a.second->second));
          it_int = interp_segments.find(&(wpb_a.first->second));
        } else if( it_int->second.get_end_point() != &(wpb_a.second->second) ) {
          it_int->second.set_segment(&(wpb_a.first->second),&(wpb_a.second->second));
        };
        sum += it_int->second.travel_distance_from(a, this->dist);
      };
      const_waypoint_descriptor it = wpb_a.second;
      const_waypoint_descriptor it_prev = it;
      while(it++ != wpb_b.first) {
        typename interpolator_map_type::iterator it_int = interp_segments.find(&(it_prev->second));
        if(it_int == interp_segments.end()) {
          interp_segments[&(it_prev->second)] = interp_fact.create_interpolator(&(it_prev->second),&(it->second));
          it_int = interp_segments.find(&(it_prev->second));
        } else if( it_int->second.get_end_point() != &(it->second) ) {
          it_int->second.set_segment(&(it_prev->second),&(it->second));
        };
        sum += it_int->second.travel_distance_from(it_prev->second, this->dist);
        ++it_prev;
      };
      {
        typename interpolator_map_type::iterator it_int = interp_segments.find(&(wpb_b.first->second));
        if(it_int == interp_segments.end()) {
          if(wpb_b.first == wpb_b.second) {
            interpolator_type seg = interp_fact.create_interpolator(&(wpb_b.first->second),&b);
            sum += seg.travel_distance_to(b, this->dist);
          } else {
            interp_segments[&(wpb_b.first->second)] = interp_fact.create_interpolator(&(wpb_b.first->second),&(wpb_b.second->second));
            sum += interp_segments[&(wpb_b.first->second)].travel_distance_to(b, this->dist);
          };
        } else if( it_int->second.get_end_point() != &(wpb_b.second->second) ) {
          it_int->second.set_segment(&(wpb_b.first->second),&(wpb_b.second->second));
          sum += it_int->second.travel_distance_to(b, this->dist);
        } else {
          sum += it_int->second.travel_distance_to(b, this->dist);
        };
      };
      return sum;
    };
    
    waypoint_pair get_point_at_time_impl(double t, const const_waypoint_bounds& wpb_a) const {
      
      if( wpb_a.first == wpb_a.second ) {
        // one way or another, the point is at the boundary:
        waypoint_pair result(wpb_a.first, wpb_a.first->second);
        result.second.time = t;
        return result;
      };
      
      typename interpolator_map_type::iterator it_int = interp_segments.find(&(wpb_a.first->second));
      if(it_int == interp_segments.end()) {
        return waypoint_pair( wpb_a.first,
          (interp_segments[&(wpb_a.first->second)] = interp_fact.create_interpolator(&(wpb_a.first->second),&(wpb_a.second->second)))
          .get_point_at_time(t));
      } else if(it_int->second.get_end_point() != &(wpb_a.second->second)) {
        it_int->second.set_segment(&(wpb_a.first->second),&(wpb_a.second->second));
      };
      return waypoint_pair(wpb_a.first, it_int->second.get_point_at_time(t));
    };
    
    virtual waypoint_pair move_time_diff_from_impl(const point_type& a, const const_waypoint_bounds& wpb_a, double dt) const {
      if( ( a.time + dt >= wpb_a.first->first ) && ( a.time + dt <= wpb_a.second->first ) )
        return get_point_at_time_impl(a.time + dt, wpb_a);
      else
        return get_point_at_time_impl(a.time + dt, this->get_waypoint_bounds(a.time + dt));
    };
    
  public:
    
    
    
    struct point_time_iterator {
      const interpolated_trajectory* parent;
      const_waypoint_bounds current_wpbound;
      point_type current_pt;
      double current_dt;
      
      point_time_iterator(const interpolated_trajectory* aParent,
                          const const_waypoint_bounds& aWPB, 
                          const point_type& aCurrentPt) :
                          parent(aParent),
                          current_wpbound(aWPB), 
                          current_pt(aCurrentPt),
                          current_dt(current_pt.time - current_wpbound.first->first) {
        if(current_wpbound.first == current_wpbound.second)
          current_dt = 0.0;
      };
      
      point_time_iterator(const interpolated_trajectory* aParent,
                          const const_waypoint_bounds& aWPB) :
                          parent(aParent),
                          current_wpbound(aWPB), 
                          current_pt(current_wpbound.first->second), 
                          current_dt(0.0) { };
                              
#ifdef RK_ENABLE_CXX11_FEATURES
      point_time_iterator(const interpolated_trajectory* aParent,
                          const_waypoint_bounds&& aWPB, 
                          point_type&& aCurrentPt) :
                          parent(aParent),
                          current_wpbound(std::move(aWPB)), 
                          current_pt(std::move(aCurrentPt)), 
                          current_dt(current_pt.time - current_wpbound.first->first) {
        if(current_wpbound.first == current_wpbound.second)
          current_dt = 0.0;
      };
      
      point_time_iterator(const interpolated_trajectory* aParent,
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
        
        current_pt = parent->move_time_diff_from_impl(current_pt, current_wpbound, rhs).second;
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
      const interpolated_trajectory* parent;
      const_waypoint_bounds current_wpbound;
      point_type current_pt;
      double current_fraction;
      
      point_fraction_iterator(const interpolated_trajectory* aParent,
                              const const_waypoint_bounds& aWPB, 
                              const point_type& aCurrentPt) :
                              parent(aParent),
                              current_wpbound(aWPB), 
                              current_pt(aCurrentPt),
                              current_fraction((current_pt.time - current_wpbound.first->first) / (current_wpbound.second->first - current_wpbound.first->first)) {
        if(current_wpbound.first == current_wpbound.second)
          current_fraction = 0.0;
      };
      
      point_fraction_iterator(const interpolated_trajectory* aParent,
                              const const_waypoint_bounds& aWPB) :
                              parent(aParent),
                              current_wpbound(aWPB), 
                              current_pt(current_wpbound.first->second),
                              current_fraction(0.0) { };
                              
#ifdef RK_ENABLE_CXX11_FEATURES
      point_fraction_iterator(const interpolated_trajectory* aParent,
                              const_waypoint_bounds&& aWPB, 
                              point_type&& aCurrentPt) :
                              parent(aParent),
                              current_wpbound(std::move(aWPB)), 
                              current_pt(std::move(aCurrentPt)),
                              current_fraction((current_pt.time - current_wpbound.first->first) / (current_wpbound.second->first - current_wpbound.first->first)) {
        if(current_wpbound.first == current_wpbound.second)
          current_fraction = 0.0;
      };
      
      point_fraction_iterator(const interpolated_trajectory* aParent,
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
        
        current_pt = parent->move_time_diff_from_impl(current_pt, current_wpbound, 
                                                      rhs * (current_wpbound.second->second - current_wpbound.first->second)).second;
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
     * \param aInterp The interpolator functor that the trajectory should use.
     */
    explicit interpolated_trajectory(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), const distance_metric& aDist = distance_metric(), const interpolator_factory_type& aInterpFactory = interpolator_factory_type()) : 
                                     base_class_type(aSpace, aDist), interp_fact(aInterpFactory), interp_segments() { 
      interp_fact.set_temporal_space(this->space);
    };
    
    /**
     * Constructs the trajectory from a space, the start and end points.
     * \param aSpace The space on which the trajectory is.
     * \param aStart The start point of the trajectory.
     * \param aEnd The end-point of the trajectory.
     * \param aDist The distance metric functor that the trajectory should use.
     * \param aInterp The interpolator functor that the trajectory should use.
     */
    interpolated_trajectory(const shared_ptr<topology>& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric(), const interpolator_factory_type& aInterpFactory = interpolator_factory_type()) :
                            base_class_type(aSpace, aStart, aEnd, aDist), interp_fact(aInterpFactory), interp_segments() { 
      interp_fact.set_temporal_space(this->space);
    };
			
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
    interpolated_trajectory(ForwardIter aBegin, ForwardIter aEnd, const shared_ptr<topology>& aSpace, const distance_metric& aDist = distance_metric(), const interpolator_factory_type& aInterpFactory = interpolator_factory_type()) : 
                            base_class_type(aBegin, aEnd, aSpace, aDist), interp_fact(aInterpFactory), interp_segments() { 
      interp_fact.set_temporal_space(this->space);
    };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(static_cast<base_class_type&>(lhs),static_cast<base_class_type&>(rhs));
      swap(lhs.interp_fact, rhs.interp_fact);
      swap(lhs.interp_segments, rhs.interp_segments);
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
     * Computes the point that is a time-difference away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param dt The time to move away from the point.
     * \return The point that is a time away from the given point.
     */
    point_type move_time_diff_from(const point_type& a, double dt) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.time + dt);
      return this->move_time_diff_from_impl(wpb_a.first->second, wpb_a, a.time + dt - wpb_a.first->first).second;
    };
    
    /**
     * Computes the waypoint-point-pair that is a time away from a waypoint-point-pair on the trajectory.
     * \param a The waypoint-point-pair on the trajectory.
     * \param dt The time to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_time_diff_from(const waypoint_pair& a, double dt) const {
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second.time + dt);
      return this->move_time_diff_from_impl(wpb_a.first->second, wpb_a, a.second.time + dt - wpb_a.first->first);
    };
       
    /**
     * Computes the point that is on the trajectory at the given time.
     * \param t The time at which the point is sought.
     * \return The point that is on the trajectory at the given time.
     */
    point_type get_point_at_time(double t) const {
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(t);
      return this->move_time_diff_from_impl(wpb_p.first->second, wpb_p, t - wpb_p.first->first).second;
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given time.
     * \param t The time at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given time.
     */
    waypoint_pair get_waypoint_at_time(double t) const {
      const_waypoint_bounds wpb_p = this->get_waypoint_bounds(t);
      return this->move_time_diff_from_impl(wpb_p.first->second, wpb_p, t - wpb_p.first->first);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(interp_fact);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(interp_fact);
      interp_segments.clear();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440002,1,"interpolated_trajectory",base_class_type)
    
};



};

};

#endif









