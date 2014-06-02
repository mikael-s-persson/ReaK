/**
 * \file path_wrapper.hpp
 * 
 * This library provides a path-wrapper class template which makes an OOP-compatible path class 
 * for a given path (in the generic sense).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PATH_WRAPPER_HPP
#define REAK_PATH_WRAPPER_HPP

#include <ReaK/core/base/defs.hpp>

#include "path_base.hpp"
#include "spatial_path_concept.hpp"

#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  


/**
 * This class wraps a generic spatial path class into an OOP interface.
 * It, itself, also models the generic SpatialPathConcept, so this wrapper can 
 * be used for both purposes.
 * \tparam SpatialPath The path type to be wrapped.
 */
template <typename SpatialPath>
class path_wrapper : public path_base< typename spatial_path_traits<SpatialPath>::topology > {
  public:
    typedef path_base< typename spatial_path_traits<SpatialPath>::topology > base_type;
    typedef path_wrapper<SpatialPath> self;
    
    typedef typename base_type::topology topology;
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    BOOST_CONCEPT_ASSERT((SpatialPathConcept<SpatialPath,topology>));
    
    typedef typename spatial_path_traits<SpatialPath>::waypoint_descriptor waypoint_descriptor;
    typedef typename spatial_path_traits<SpatialPath>::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename spatial_path_traits<SpatialPath>::distance_metric distance_metric;
    
    typedef std::pair< const_waypoint_descriptor, point_type> waypoint_pair;
    
  protected:
    SpatialPath m_traj;
    mutable waypoint_pair m_last_waypoint;
    
  public:
    
    SpatialPath& get_underlying_path() { return m_traj; };
    const SpatialPath& get_underlying_path() const { return m_traj; };
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     * \param aTraj The wrapped path object to use.
     */
    explicit path_wrapper(const std::string& aName = "",
                          const SpatialPath& aTraj = SpatialPath()) : 
                          base_type(aName),
                          m_traj(aTraj), m_last_waypoint() { };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    virtual double travel_distance(const point_type& a, const point_type& b) const {
      if(m_last_waypoint.first == const_waypoint_descriptor())
        m_last_waypoint = m_traj.get_start_waypoint();
      waypoint_pair next_waypoint = m_last_waypoint;
      m_last_waypoint.second = a;
      next_waypoint.second = b;
      return m_traj.travel_distance(m_last_waypoint, next_waypoint);
    };
    
    /**
     * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
     * \param a The first waypoint-point-pair.
     * \param b The second waypoint-point-pair.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
      m_last_waypoint = a;
      return m_traj.travel_distance(m_last_waypoint, b);
    };
    
    /**
     * Computes the point that is a distance away from a point on the path.
     * \param a The point on the path.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    virtual point_type move_away_from(const point_type& a, double d) const {
      if(m_last_waypoint.first == const_waypoint_descriptor())
        m_last_waypoint = m_traj.get_start_waypoint();
      m_last_waypoint.second = a;
      m_last_waypoint = m_traj.move_away_from(m_last_waypoint, d);
      return m_last_waypoint.second;
    };
    
    /**
     * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the path.
     * \param a The waypoint-point-pair on the path.
     * \param dt The distance to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a distance away from the given waypoint-point-pair.
     */
    waypoint_pair move_away_from(const waypoint_pair& a, double d) const {
      m_last_waypoint = m_traj.move_away_from(a, d);
      return m_last_waypoint;
    };
    
    /**
     * Returns the starting point of the path.
     * \return The starting point of the path.
     */
    virtual point_type get_start_point() const {
      return m_traj.get_start_point();
    };
    
    /**
     * Returns the end point of the path.
     * \return The end point of the path.
     */
    virtual point_type get_end_point() const {
      return m_traj.get_end_point();
    };
    
    /**
     * Returns the starting waypoint-point-pair of the path.
     * \return The starting waypoint-point-pair of the path.
     */
    waypoint_pair get_start_waypoint() const {
      return m_traj.get_start_waypoint();
    };
    
    /**
     * Returns the end waypoint-point-pair of the path.
     * \return The end waypoint-point-pair of the path.
     */
    waypoint_pair get_end_waypoint() const {
      return m_traj.get_end_waypoint();
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_traj);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_traj);
      m_last_waypoint = waypoint_pair();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC244000D,1,"path_wrapper",base_type)
  
};


};

};

#endif









