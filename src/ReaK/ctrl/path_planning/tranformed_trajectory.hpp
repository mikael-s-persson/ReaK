/**
 * \file transformed_trajectory.hpp
 * 
 * This library provides an implementation of a trajectory which is a map between an underlying
 * trajectory (on its topology) to another topology using a homeomorphism across the topologies.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_TRANSFORMED_TRAJECTORY_HPP
#define REAK_TRANSFORMED_TRAJECTORY_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include "spatial_trajectory_concept.hpp"
#include "topological_map_concepts.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include <map>

namespace ReaK {

namespace pp {
  

/**
 * This class implements a trajectory which is a map between an underlying trajectory (on its 
 * own topology) to another topology using a homeomorphism across the topologies. 
 * This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, 
 *                  should model the TemporalSpaceConcept.
 * \tparam InputTrajectory The underlying trajectory, should model the SpatialTrajectoryConcept.
 * \tparam Mapping The homeomorphic mapping between the spatial topology of the underlying 
 *                 topology (of the InputTrajectory) and the given topology.
 */
template <typename Topology, typename InputTrajectory, typename Mapping>
class transformed_trajectory : public shared_object {
  public:
    typedef typename spatial_trajectory_traits<InputTrajectory>::topology input_topology;
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SpatialTrajectoryConcept<InputTrajectory,input_topology>));
    BOOST_CONCEPT_ASSERT((HomeomorphismConcept<Mapping,input_topology,Topology>));
    
    typedef transformed_trajectory<Topology,InputTrajectory,Mapping> self;
    
    typedef typename spatial_trajectory_traits<InputTrajectory>::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename topology_traits<Topology>::point_type point_type;
    typedef Topology topology;
    typedef typename spatial_trajectory_traits<InputTrajectory>::distance_metric distance_metric;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  private:
    
    typename shared_pointer<const Topology>::type space;
    typename shared_pointer<const InputTrajectory>::type traject;
    Mapping map;
    
  public:
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the trajectory is.
     * \param aTrajectory The underlying trajectory to use.
     * \param aMap The homeomorphic mapping object to use.
     */
    explicit transformed_trajectory(const typename shared_pointer<const Topology>::type& aSpace = typename shared_pointer<const Topology>::type(new Topology()), 
				    const typename shared_pointer<const InputTrajectory>::type& aTrajectory = typename shared_pointer<const InputTrajectory>::type(new InputTrajectory()), 
				    const Mapping& aMap = Mapping()) : 
                                    space(aSpace), traject(aTrajectory), map(aMap) { };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(const point_type& a, const point_type& b) const {
      return traject->travel_distance(map.map_to_space(a,*space,traject->get_temporal_space()),
			 	      map.map_to_space(b,*space,traject->get_temporal_space()));
    };
    
    /**
     * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
     * \param a The first waypoint-point-pair.
     * \param b The second waypoint-point-pair.
     * \return The travel distance between two points if traveling along the path.
     */
    double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
      return traject->travel_distance(std::make_pair(a.first, map.map_to_space(a.second,*space,traject->get_temporal_space())),
				      std::make_pair(b.first, map.map_to_space(b.second,*space,traject->get_temporal_space())));
    };
    
    
    /**
     * Computes the point that is a time-difference away from a point on the trajectory.
     * \param a The point on the trajectory.
     * \param dt The time to move away from the point.
     * \return The point that is a time away from the given point.
     */
    point_type move_time_diff_from(const point_type& a, double dt) const {
      return map.map_to_space(
	traject->move_time_diff_from( map.map_to_space(a,*space,traject->get_temporal_space()),dt),
	traject->get_temporal_space(), *space);
    };
    
    /**
     * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
     * \param a The waypoint-point-pair on the trajectory.
     * \param dt The time to move away from the waypoint-point-pair.
     * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
     */
    waypoint_pair move_time_diff_from(const waypoint_pair& a, double dt) const {
      std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits<InputTrajectory>::point_type> result = 
        traject->move_time_diff_from( std::make_pair(a.second,map.map_to_space(a.second,*space,traject->get_temporal_space())), dt);
      return std::make_pair( result.first, map.map_to_space(result.second, traject->get_temporal_space(), *space));
    };
       
    /**
     * Computes the point that is on the trajectory at the given time.
     * \param t The time at which the point is sought.
     * \return The point that is on the trajectory at the given time.
     */
    point_type get_point_at_time(double t) const {
      return map.map_to_space(traject->get_point_at_time(t), traject->get_temporal_space(), *space);
    };
    
    /**
     * Computes the waypoint-point pair that is on the trajectory at the given time.
     * \param t The time at which the waypoint-point pair is sought.
     * \return The waypoint-point pair that is on the trajectory at the given time.
     */
    waypoint_pair get_waypoint_at_time(double t) const {
      std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits<InputTrajectory>::point_type> result = 
        traject->get_waypoint_at_time(t);
      return std::make_pair( result.first, map.map_to_space(result.second, traject->get_temporal_space(), *space));
    };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& get_temporal_space() const throw() { return *space; };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(space)
        & RK_SERIAL_SAVE_WITH_NAME(traject)
	& RK_SERIAL_SAVE_WITH_NAME(map);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(space)
        & RK_SERIAL_LOAD_WITH_NAME(traject)
	& RK_SERIAL_LOAD_WITH_NAME(map);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440008,1,"transformed_trajectory",shared_object)
    
};



};

};

#endif









