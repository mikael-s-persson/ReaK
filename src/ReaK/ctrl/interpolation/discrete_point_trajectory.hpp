/**
 * \file discrete_point_trajectory.hpp
 * 
 * This library provides an implementation of a trajectory represented by discrete points within a temporal topology.
 * The trajectory is represented by a set of waypoints (presumably close to each other) and traveling along the 
 * trajectory is restricted to hopping between discrete waypoints.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_DISCRETE_POINT_TRAJECTORY_HPP
#define REAK_DISCRETE_POINT_TRAJECTORY_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/ctrl/path_planning/spatial_trajectory_concept.hpp>

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
class discrete_point_trajectory : public waypoint_container<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef discrete_point_trajectory<Topology,DistanceMetric> self;
    typedef waypoint_container<Topology,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef typename base_class_type::waypoint_pair waypoint_pair;
    
  protected:
    
    typedef typename base_class_type::container_type container_type;
    
    
    virtual double travel_distance_impl(const point_type& a, const point_type& b) const {
      if(a.time > b.time)
        return travel_distance_impl(b, a);
      const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.time);
      const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b.time);
      
      // first assume that a is before b:
      double sum = 0.0;
      const_waypoint_descriptor a_cpy = wpb_a.first;
      while(a_cpy != wpb_b.first) {
        const_waypoint_descriptor a_next = a_cpy; ++a_next;
        if(a_next == this->waypoints.end())
          break; // we're done.
        sum += this->getDistanceMetric()(a_cpy->second, a_next->second, this->getSpace());
        a_cpy = a_next;
      };
      return sum;
    };
    
    virtual waypoint_pair move_time_diff_from_impl(const point_type& a, const const_waypoint_bounds& wpb_a, double dt) const {
      const_waypoint_bounds wpb_p = wpb_a;
      if( ( a.time + dt < wpb_p.first->first ) || ( a.time + dt > wpb_p.second->first ) ) {
        wpb_p = this->get_waypoint_bounds(a.time + dt);
        dt = a.time + dt - wpb_p.first->first;
      };
      
      if(dt < 0.0) // is at start
        return waypoint_pair(wpb_p.first, wpb_p.first->second);
      
      const_waypoint_descriptor a_next = wpb_p.first; ++a_next;
      if( a_next == this->waypoints.end() ) // is at end
        return waypoint_pair(wpb_p.first, wpb_p.first->second);
      
      if( dt < 0.5 * (a_next->first - wpb_p.first->first) ) // round up or down
        return waypoint_pair(wpb_p.first, wpb_p.first->second);
      else
        return waypoint_pair(a_next, a_next->second);
    };
    
  public:
    
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the trajectory is.
     * \param aDist The distance metric functor that the trajectory should use.
     */
    explicit discrete_point_trajectory(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), 
                                       const distance_metric& aDist = distance_metric()) : 
                                       base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the trajectory from a space, the start and end points.
     * \param aSpace The space on which the trajectory is.
     * \param aStart The start point of the trajectory.
     * \param aEnd The end-point of the trajectory.
     * \param aDist The distance metric functor that the trajectory should use.
     */
    discrete_point_trajectory(const shared_ptr<topology>& aSpace, 
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
    discrete_point_trajectory(ForwardIter aBegin, ForwardIter aEnd, 
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
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440016,1,"discrete_point_trajectory",base_class_type)
    
};



};

};

#endif









