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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include "spatial_trajectory_concept.hpp"
#include <ReaK/topologies/spaces/topological_map_concepts.hpp>
#include "sequential_trajectory_concept.hpp"

#include <boost/concept_check.hpp>

#include <map>

namespace ReaK {

namespace pp {


namespace detail {


template < typename ParentType, typename InputIterator >
struct transformed_point_iterator {
  typedef transformed_point_iterator< ParentType, InputIterator > self;

  const ParentType* p_parent;
  InputIterator base_it;

  transformed_point_iterator( const ParentType* aParent, InputIterator aBaseIt )
      : p_parent( aParent ), base_it( aBaseIt ){};

  friend self& operator+=( self& lhs, double rhs ) {
    lhs.base_it += rhs;
    return lhs;
  };

  friend self& operator-=( self& lhs, double rhs ) {
    lhs.base_it -= rhs;
    return lhs;
  };

  friend self operator+( self lhs, double rhs ) { return ( lhs += rhs ); };
  friend self operator+( double lhs, self rhs ) { return ( rhs += lhs ); };
  friend self operator-( self lhs, double rhs ) { return ( lhs -= rhs ); };

  friend bool operator==( const self& lhs, const self& rhs ) { return ( lhs.base_it == rhs.base_it ); };

  friend bool operator!=( const self& lhs, const self& rhs ) { return !( lhs == rhs ); };

  typename ParentType::point_type operator*() const { return p_parent->map_point_forward( *base_it ); };
};
};


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
template < typename Topology, typename InputTrajectory, typename Mapping >
class transformed_trajectory : public shared_object {
public:
  typedef typename spatial_trajectory_traits< InputTrajectory >::topology input_topology;

  BOOST_CONCEPT_ASSERT( ( TemporalSpaceConcept< Topology > ) );
  BOOST_CONCEPT_ASSERT( ( SpatialTrajectoryConcept< InputTrajectory, input_topology > ) );
  BOOST_CONCEPT_ASSERT( ( BijectionConcept< Mapping, input_topology, Topology > ) );

  typedef transformed_trajectory< Topology, InputTrajectory, Mapping > self;

  typedef typename spatial_trajectory_traits< InputTrajectory >::waypoint_descriptor waypoint_descriptor;
  typedef typename spatial_trajectory_traits< InputTrajectory >::const_waypoint_descriptor const_waypoint_descriptor;
  typedef typename topology_traits< Topology >::point_type point_type;
  typedef typename topology_traits< Topology >::point_difference_type point_difference_type;
  typedef Topology topology;
  typedef typename spatial_trajectory_traits< InputTrajectory >::distance_metric distance_metric;

  typedef std::pair< const_waypoint_descriptor, point_type > waypoint_pair;

  typedef detail::transformed_point_iterator< self, typename sequential_trajectory_traits< InputTrajectory >::
                                                      point_time_iterator > point_time_iterator;

  typedef detail::transformed_point_iterator< self, typename sequential_trajectory_traits< InputTrajectory >::
                                                      point_fraction_iterator > point_fraction_iterator;

  friend struct detail::
    transformed_point_iterator< self, typename sequential_trajectory_traits< InputTrajectory >::point_time_iterator >;

  friend struct detail::transformed_point_iterator< self, typename sequential_trajectory_traits< InputTrajectory >::
                                                            point_fraction_iterator >;

private:
  shared_ptr< Topology > space;
  shared_ptr< InputTrajectory > traject;
  Mapping map;


  point_type map_point_forward( const typename spatial_trajectory_traits< InputTrajectory >::point_type& pt_in ) const {
    return map.map_to_space( pt_in, traject->get_temporal_space(), *space );
  };

  typename spatial_trajectory_traits< InputTrajectory >::point_type
    map_point_backward( const point_type& pt_in ) const {
    return map.map_to_space( pt_in, *space, traject->get_temporal_space() );
  };

public:
  /**
   * Returns the space on which the path resides.
   * \return The space on which the path resides.
   */
  const topology& get_temporal_space() const throw() { return *space; };


  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the trajectory is.
   * \param aTrajectory The underlying trajectory to use.
   * \param aMap The homeomorphic mapping object to use.
   */
  explicit transformed_trajectory( const shared_ptr< Topology >& aSpace = shared_ptr< Topology >(),
                                   const shared_ptr< InputTrajectory >& aTrajectory = shared_ptr< InputTrajectory >(),
                                   const Mapping& aMap = Mapping() )
      : space( aSpace ), traject( aTrajectory ), map( aMap ){};

  /**
   * Computes the travel distance between two points, if traveling along the path.
   * \param a The first point.
   * \param b The second point.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance( const point_type& a, const point_type& b ) const {
    return get( ReaK::pp::distance_metric, *space )( a, b, *space );
    // NOTE: requiring a bijection is too much:
    //       return traject->travel_distance(map_point_backward(a), map_point_backward(b));
  };

  /**
   * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
   * \param a The first waypoint-point-pair.
   * \param b The second waypoint-point-pair.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance( waypoint_pair& a, waypoint_pair& b ) const {
    return get( ReaK::pp::distance_metric, *space )( a.second, b.second, *space );
    // NOTE: requiring a bijection is too much:
    //       return traject->travel_distance(std::make_pair(a.first, map_point_backward(a.second)),
    //                                       std::make_pair(b.first, map_point_backward(b.second)));
  };

  /**
   * Returns the total travel-distance of the trajectory.
   * \return The total travel-distance of the trajectory.
   */
  double get_total_length() const { return traject->get_total_length(); };


  /**
   * Computes the point that is a time-difference away from a point on the trajectory.
   * \param a The point on the trajectory.
   * \param dt The time to move away from the point.
   * \return The point that is a time away from the given point.
   */
  point_type move_time_diff_from( const point_type& a, double dt ) const {
    return get_point_at_time( a.time + dt );
    // NOTE: requiring a bijection is too much!
    //       return map_point_forward(traject->move_time_diff_from(map_point_backward(a),dt));
  };

  /**
   * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
   * \param a The waypoint-point-pair on the trajectory.
   * \param dt The time to move away from the waypoint-point-pair.
   * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
   */
  waypoint_pair move_time_diff_from( const waypoint_pair& a, double dt ) const {
    return get_waypoint_at_time( a.second.time + dt );
    // NOTE: requiring a bijection is too much!
    //       std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits<InputTrajectory>::point_type>
    //       result =
    //         traject->move_time_diff_from( std::make_pair(a.first,map_point_backward(a.second)), dt);
    //       return std::make_pair( result.first, map_point_forward(result.second));
  };

  /**
   * Computes the point that is on the trajectory at the given time.
   * \param t The time at which the point is sought.
   * \return The point that is on the trajectory at the given time.
   */
  point_type get_point_at_time( double t ) const { return map_point_forward( traject->get_point_at_time( t ) ); };

  /**
   * Computes the waypoint-point pair that is on the trajectory at the given time.
   * \param t The time at which the waypoint-point pair is sought.
   * \return The waypoint-point pair that is on the trajectory at the given time.
   */
  waypoint_pair get_waypoint_at_time( double t ) const {
    std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits< InputTrajectory >::point_type > result
      = traject->get_waypoint_at_time( t );
    return std::make_pair( result.first, map_point_forward( result.second ) );
  };

  /**
   * Returns the starting time of the trajectory.
   * \return The starting time of the trajectory.
   */
  double get_start_time() const { return traject->get_start_time(); };

  /**
   * Returns the end time of the trajectory.
   * \return The end time of the trajectory.
   */
  double get_end_time() const { return traject->get_end_time(); };


  /**
   * Returns the starting point of the waypoints.
   * \return The starting point of the waypoints.
   */
  const point_type& get_start_point() const { return map_point_forward( traject->get_start_point() ); };

  /**
   * Returns the starting waypoint-point-pair of the waypoints.
   * \return The starting waypoint-point-pair of the waypoints.
   */
  waypoint_pair get_start_waypoint() const {
    std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits< InputTrajectory >::point_type > result
      = traject->get_start_waypoint();
    return std::make_pair( result.first, map_point_forward( result.second ) );
  };

  /**
   * Returns the end point of the waypoints.
   * \return The end point of the waypoints.
   */
  const point_type& get_end_point() const { return map_point_forward( traject->get_end_point() ); };

  /**
   * Returns the end waypoint-point-pair of the waypoints.
   * \return The end waypoint-point-pair of the waypoints.
   */
  waypoint_pair get_end_waypoint() const {
    std::pair< const_waypoint_descriptor, typename spatial_trajectory_traits< InputTrajectory >::point_type > result
      = traject->get_end_waypoint();
    return std::make_pair( result.first, map_point_forward( result.second ) );
  };


  /**
   * Returns the starting time-iterator along the trajectory.
   * \return The starting time-iterator along the trajectory.
   */
  point_time_iterator begin_time_travel() const { return point_time_iterator( this, traject->begin_time_travel() ); };

  /**
   * Returns the end time-iterator along the trajectory.
   * \return The end time-iterator along the trajectory.
   */
  point_time_iterator end_time_travel() const { return point_time_iterator( this, traject->end_time_travel() ); };

  /**
   * Returns the starting fraction-iterator along the trajectory.
   * \return The starting fraction-iterator along the trajectory.
   */
  point_fraction_iterator begin_fraction_travel() const {
    return point_fraction_iterator( this, traject->begin_fraction_travel() );
  };

  /**
   * Returns the end fraction-iterator along the trajectory.
   * \return The end fraction-iterator along the trajectory.
   */
  point_fraction_iterator end_fraction_travel() const {
    return point_fraction_iterator( this, traject->end_fraction_travel() );
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    shared_object::save( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( space ) & RK_SERIAL_SAVE_WITH_NAME( traject ) & RK_SERIAL_SAVE_WITH_NAME( map );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    shared_object::load( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( space ) & RK_SERIAL_LOAD_WITH_NAME( traject ) & RK_SERIAL_LOAD_WITH_NAME( map );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2440008, 1, "transformed_trajectory", shared_object )
};
};
};

#endif
