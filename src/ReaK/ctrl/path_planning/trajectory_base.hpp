/**
 * \file trajectory_base.hpp
 * 
 * This library provides the base-class for trajectories within a temporal topology.
 * This is a base-class that stems the object-oriented compatibility of other temporal
 * trajectory classes. Then, this library provides a trajectory-wrapper class template 
 * which makes an OOP-compatible trajectory class for a given trajectory.
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

#ifndef REAK_TRAJECTORY_BASE_HPP
#define REAK_TRAJECTORY_BASE_HPP

#include "base/defs.hpp"

#include "base/named_object.hpp"

#include "spatial_trajectory_concept.hpp"
#include "sequential_trajectory_concept.hpp"

#include "seq_trajectory_base.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include <lin_alg/mat_num_exceptions.hpp>

namespace ReaK {

namespace pp {
  

/**
 * This class defines the OOP interface for a trajectory in a temporal topology.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept.
 */
template <typename Topology>
class trajectory_base : public seq_trajectory_base<Topology> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef Topology topology;
    typedef typename topology_traits<topology>::point_type point_type;
    typedef typename topology_traits<topology>::point_difference_type point_difference_type;
    typedef trajectory_base<Topology> self;
    typedef seq_trajectory_base<Topology> base_type;
    
    
  public:
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     */
    explicit trajectory_base(const std::string& aName) : 
                             base_type(aName) {};
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    virtual double travel_distance(const point_type& a, const point_type& b) const = 0;
    
    /**
     * Computes the point that is a time-difference away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param dt The time to move away from the point.
     * \return The point that is a time away from the given point.
     */
    virtual point_type move_time_diff_from(const point_type& a, double dt) const = 0;
    
    /**
     * Computes the point that is on the trajectory at the given time.
     * \param t The time at which the point is sought.
     * \return The point that is on the trajectory at the given time.
     */
    virtual point_type get_point_at_time(double t) const = 0;
    
    /**
     * Returns the starting time of the trajectory.
     * \return The starting time of the trajectory.
     */
    virtual double get_start_time() const = 0;
    
    /**
     * Returns the end time of the trajectory.
     * \return The end time of the trajectory.
     */
    virtual double get_end_time() const = 0;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2440009,1,"trajectory_base",base_type)
};


/**
 * This class wraps a generic spatial trajectory class into an OOP interface.
 * It, itself, also models the generic SpatialTrajectoryConcept and SequentialTrajectoryConcept, 
 * so this wrapper can be used for both purposes (i.e. this is a type-erasure class).
 * \tparam SpatialTrajectory The trajectory type to be wrapped.
 */
template <typename SpatialTrajectory>
class trajectory_wrapper : public trajectory_base< typename spatial_trajectory_traits<SpatialTrajectory>::topology > {
  public:
    typedef trajectory_base< typename spatial_trajectory_traits<SpatialTrajectory>::topology > base_type;
    typedef typename base_type::base_type seq_base_type;
    typedef trajectory_wrapper<SpatialTrajectory> self;
    
    typedef typename base_type::topology topology;
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    BOOST_CONCEPT_ASSERT((SequentialTrajectoryConcept<SpatialTrajectory,topology>));
    
    typedef typename spatial_trajectory_traits<SpatialTrajectory>::waypoint_descriptor waypoint_descriptor;
    typedef typename spatial_trajectory_traits<SpatialTrajectory>::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename temporal_space_traits<topology>::time_topology time_topology;
    typedef typename temporal_space_traits<topology>::space_topology space_topology;
    typedef typename spatial_trajectory_traits<SpatialTrajectory>::distance_metric distance_metric;
    
    typedef std::pair< const_waypoint_descriptor, point_type> waypoint_pair;
    
    typedef SpatialTrajectory wrapped_type;
    
  protected:
    SpatialTrajectory m_traj;
    mutable waypoint_pair m_last_waypoint;
    
    typedef typename seq_base_type::point_time_iterator_impl base_pt_time_iterator_impl;
    typedef typename sequential_trajectory_traits<SpatialTrajectory>::point_time_iterator gen_pt_time_iterator;
    
    struct point_time_iterator_impl : public base_pt_time_iterator_impl {
      
      gen_pt_time_iterator base_it;
      
      point_time_iterator_impl(gen_pt_time_iterator aBaseIt) : base_it(aBaseIt) { };
      
      virtual ~point_time_iterator_impl() { };
      
      virtual void move_by_time(double d) {
        base_it += d;
      };
      
      virtual bool is_equal_to(const base_pt_time_iterator_impl* rhs) const {
        return (base_it == static_cast<const point_time_iterator_impl*>(rhs)->base_it);
      };
      
      virtual point_type get_point() const {
        return *base_it;
      };
      
      virtual base_pt_time_iterator_impl* clone() const {
        return new point_time_iterator_impl(base_it);
      };
      
    };
    
    typedef typename seq_base_type::point_fraction_iterator_impl base_pt_frac_iterator_impl;
    typedef typename sequential_trajectory_traits<SpatialTrajectory>::point_fraction_iterator gen_pt_frac_iterator;
    
    struct point_fraction_iterator_impl : public base_pt_frac_iterator_impl {
      
      gen_pt_frac_iterator base_it;
      
      point_fraction_iterator_impl(gen_pt_frac_iterator aBaseIt) : base_it(aBaseIt) { };
      
      virtual ~point_fraction_iterator_impl() { };
      
      virtual void move_by_fraction(double f) {
        base_it += f;
      };
      
      virtual bool is_equal_to(const base_pt_frac_iterator_impl* rhs) const {
        return (base_it == static_cast<const point_fraction_iterator_impl*>(rhs)->base_it);
      };
      
      virtual point_type get_point() const {
        return *base_it;
      };
      
      virtual base_pt_frac_iterator_impl* clone() const {
        return new point_fraction_iterator_impl(base_it);
      };
      
    };
    
  public:
    
    typedef typename base_type::point_time_iterator point_time_iterator;
    typedef typename base_type::point_fraction_iterator point_fraction_iterator;
    
    wrapped_type& get_underlying_trajectory() { return m_traj; };
    const wrapped_type& get_underlying_trajectory() const { return m_traj; };
    
    wrapped_type& get_wrapped_object() { return m_traj; };
    const wrapped_type& get_wrapped_object() const { return m_traj; };
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     * \param aTraj The wrapped trajectory object to use.
     */
    explicit trajectory_wrapper(const std::string& aName = "",
                                const SpatialTrajectory& aTraj = SpatialTrajectory()) : 
                                base_type(aName),
                                m_traj(aTraj), m_last_waypoint() { };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    virtual double travel_distance(const point_type& a, const point_type& b) const {
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
     * Computes the point that is a time-difference away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param dt The time to move away from the point.
     * \return The point that is a time away from the given point.
     */
    virtual point_type move_time_diff_from(const point_type& a, double dt) const {
      m_last_waypoint.second = a;
      m_last_waypoint = m_traj.move_time_diff_from(m_last_waypoint, dt);
      return m_last_waypoint.second;
    };
    
    /**
     * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
     * \param a The waypoint-point-pair on the trajectory.
     * \param dt The time to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_time_diff_from(const waypoint_pair& a, double dt) const {
      m_last_waypoint = m_traj.move_time_diff_from(a, dt);
      return m_last_waypoint;
    };
    
    /**
     * Computes the point that is on the trajectory at the given time.
     * \param t The time at which the point is sought.
     * \return The point that is on the trajectory at the given time.
     */
    virtual point_type get_point_at_time(double t) const {
      m_last_waypoint = m_traj.get_waypoint_at_time(t);
      return m_last_waypoint.second;
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given time.
     * \param t The time at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given time.
     */
    waypoint_pair get_waypoint_at_time(double t) const {
      return (m_last_waypoint = m_traj.get_waypoint_at_time(t));
    };
    
    /**
     * Returns the starting time of the trajectory.
     * \return The starting time of the trajectory.
     */
    virtual double get_start_time() const {
      return m_traj.get_start_time();
    };
    
    /**
     * Returns the end time of the trajectory.
     * \return The end time of the trajectory.
     */
    virtual double get_end_time() const {
      return m_traj.get_end_time();
    };
    
    /**
     * Returns the starting time-iterator along the trajectory.
     * \return The starting time-iterator along the trajectory.
     */
    virtual point_time_iterator begin_time_travel() const {
      return point_time_iterator(new point_time_iterator_impl(m_traj.begin_time_travel()));
    };
    
    /**
     * Returns the end time-iterator along the trajectory.
     * \return The end time-iterator along the trajectory.
     */
    virtual point_time_iterator end_time_travel() const {
      return point_time_iterator(new point_time_iterator_impl(m_traj.end_time_travel()));
    };
    
    /**
     * Returns the starting fraction-iterator along the trajectory.
     * \return The starting fraction-iterator along the trajectory.
     */
    virtual point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.begin_fraction_travel()));
    };
    
    /**
     * Returns the end fraction-iterator along the trajectory.
     * \return The end fraction-iterator along the trajectory.
     */
    virtual point_fraction_iterator end_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.end_fraction_travel()));
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC244000A,1,"trajectory_wrapper",base_type)
  
};


};

};

#endif









