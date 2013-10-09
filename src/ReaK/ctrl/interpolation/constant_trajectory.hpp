/**
 * \file constant_trajectory.hpp
 * 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
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

#ifndef REAK_CONSTANT_TRAJECTORY_HPP
#define REAK_CONSTANT_TRAJECTORY_HPP

#include "path_planning/spatial_path_concept.hpp"
#include "path_planning/temporal_space_concept.hpp"

#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

#include "topologies/time_topology.hpp"
#include "topologies/temporal_space.hpp"

namespace ReaK {

namespace pp {

  
template <typename SpaceType>
class constant_trajectory : public shared_object {
  public:
    typedef constant_trajectory<SpaceType> self;
    
    struct waypoint_descriptor { };
    struct const_waypoint_descriptor { };
    
    typedef temporal_space< SpaceType, ::ReaK::pp::time_topology > topology;
    typedef ::ReaK::pp::time_topology time_topology;
    typedef SpaceType space_topology;
    typedef spatial_distance_only distance_metric;
    
    typedef typename topology_traits< topology >::point_type point_type;
    typedef typename topology_traits< topology >::point_difference_type point_difference_type;
    
    typedef typename topology_traits< SpaceType >::point_type space_point_type;
    
    typedef std::pair< const_waypoint_descriptor, point_type> waypoint_pair_type;
    
  private:
    point_type start_point;
    
  public:
    
    
    struct point_time_iterator {
      point_type current_pt;
      
      explicit point_time_iterator(const point_type& aCurrentPt) : current_pt(aCurrentPt) { };
                              
#ifdef RK_ENABLE_CXX11_FEATURES
      explicit point_time_iterator(point_type&& aCurrentPt) : current_pt(std::move(aCurrentPt)) { };
#endif
      
      friend point_time_iterator operator+(const point_time_iterator& lhs, double) { return lhs; };
      friend point_time_iterator operator+(double, const point_time_iterator& rhs) { return rhs; };
      friend point_time_iterator operator-(const point_time_iterator& lhs, double) { return lhs; };
      
      friend point_time_iterator& operator+=(point_time_iterator& lhs, double) { return lhs; };
      friend point_time_iterator& operator-=(point_time_iterator& lhs, double) { return lhs; };
      
      friend bool operator==(const point_time_iterator& lhs, const point_time_iterator& rhs) { return true; };
      friend bool operator!=(const point_time_iterator& lhs, const point_time_iterator& rhs) { return false; };
      
      const point_type& operator*() const { return current_pt; };
      
    };
    
    typedef point_time_iterator point_fraction_iterator;
    
    
    /**
     * Returns the starting time-iterator along the trajectory.
     * \return The starting time-iterator along the trajectory.
     */
    point_time_iterator begin_time_travel() const { return point_time_iterator(start_point); };
    
    /**
     * Returns the end time-iterator along the trajectory.
     * \return The end time-iterator along the trajectory.
     */
    point_time_iterator end_time_travel() const { return point_time_iterator(start_point); };
    
    /**
     * Returns the starting fraction-iterator along the trajectory.
     * \return The starting fraction-iterator along the trajectory.
     */
    point_fraction_iterator begin_fraction_travel() const { return point_fraction_iterator(start_point); };
    
    /**
     * Returns the end fraction-iterator along the trajectory.
     * \return The end fraction-iterator along the trajectory.
     */
    point_fraction_iterator end_fraction_travel() const { return point_fraction_iterator(start_point); };
    
    
    
    
    constant_trajectory(const space_point_type& aStartPoint = space_point_type()) : start_point(0.0, aStartPoint) { };
    
    const point_type* get_start_point() const { return &start_point; };
    void set_start_point(const point_type* aStartPoint) { 
      if(aStartPoint)
        start_point = *aStartPoint; 
    };
    
    point_type get_point_at_time(double t) const { return point_type(t, start_point.pt); };
    
    point_type move_time_diff_from(const point_type& p, double dt) const { return point_type(p.time + dt, start_point.pt); };
    
    double travel_distance(const point_type&, const point_type&) const { return 0.0; };
    
    waypoint_pair_type get_waypoint_at_time(double t) const { 
      return waypoint_pair_type(const_waypoint_descriptor(), point_type(t, start_point.pt)); 
    };
    
    waypoint_pair_type move_time_diff_from(const waypoint_pair_type& w_p, double dt) const { 
      return waypoint_pair_type(const_waypoint_descriptor(), point_type(w_p.second.time + dt, start_point.pt)); 
    };
    
    double travel_distance(const waypoint_pair_type&, const waypoint_pair_type&) const { return 0.0; };
    
    double get_start_time() const { return -std::numeric_limits<double>::infinity(); };
    double get_end_time() const { return std::numeric_limits<double>::infinity(); }; 
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(start_point);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      ReaK::shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(start_point);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440010,1,"constant_trajectory",shared_object)
};



};

};

#endif














