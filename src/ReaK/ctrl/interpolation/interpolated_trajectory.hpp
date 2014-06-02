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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/lin_alg/mat_num_exceptions.hpp>

#include <ReaK/ctrl/path_planning/spatial_trajectory_concept.hpp>
#include <ReaK/ctrl/path_planning/interpolator_concept.hpp>

#include "waypoint_container.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>

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
    
  protected:
    
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









