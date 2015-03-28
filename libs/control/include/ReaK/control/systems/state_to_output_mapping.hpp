/**
 * \file state_to_output_mapping.hpp
 *
 * This library provides a class template which can map state-vectors to elements
 * of an output topology.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2013
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

#ifndef REAK_STATE_TO_OUTPUT_MAPPING_HPP
#define REAK_STATE_TO_OUTPUT_MAPPING_HPP

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include "state_space_sys_concept.hpp"

#include <boost/utility.hpp>

namespace ReaK {

namespace ctrl {


/**
 * This class is a bijection mapping from a belief-space to a state-space (topology)
 * by making the assumption of maximum likelihood. In other words, this mapping reduces
 * belief-states into their most likely value only.
 *
 * Models: BijectionConcept between a state-space topology and a compatible output topology.
 */
template < typename StateSpaceSystem >
struct state_to_output_map : public named_object {

  typedef state_to_output_map< StateSpaceSystem > self;

  typedef typename pp::topology_traits< BeliefSpace >::point_type belief_state_type;
  typedef typename belief_state_traits< belief_state_type >::state_type state_type;

  explicit state_to_output_map( const shared_ptr< StateSpaceSystem >& aSys = shared_ptr< StateSpaceSystem >() )
      : named_object(), ss_system( aSys ) {
    setName( "state_to_output_map" );
  };

  shared_ptr< StateSpaceSystem > ss_system;

  /**
   * This function computes the output corresponding to a given state-point.
   * \tparam StateSpace A state-space topology whose point-types are compatible with the state type of the state-space
   * system used by this object.
   * \tparam OutputSpace A topology whose point-types are compatible with the output-type of the state-space system used
   * by this object.
   * \param p The state-point for which the output is sought.
   * \param ss_space The state-space topology object.
   * \return The output corresponding to the given state-point.
   * \note This specialization is used for temporal state-spaces.
   */
  template < typename StateSpace, typename OutputSpace >
  typename boost::enable_if< pp::is_temporal_space< StateSpace >, pp::topology_traits< OutputSpace > >::type::point_type
    map_to_space( const typename pp::topology_traits< StateSpace >::point_type& p, const StateSpace& ss_space,
                  const OutputSpace& ) const {
    typedef typename ss_system_traits< StateSpaceSystem >::input_type Input;
    typedef typename pp::topology_traits< OutputSpace >::point_type OutputTopoPoint;
    return OutputTopoPoint( p.time, ss_system->get_output( ss_space.get_space_topology(), p.pt, Input(), p.time ) );
  };

  /**
   * This function computes the output corresponding to a given state-point.
   * \tparam StateSpace A state-space topology whose point-types are compatible with the state type of the state-space
   * system used by this object.
   * \tparam OutputSpace A topology whose point-types are compatible with the output-type of the state-space system used
   * by this object.
   * \param p The state-point for which the output is sought.
   * \param ss_space The state-space topology object.
   * \return The output corresponding to the given state-point.
   * \note This specialization is used for non-temporal state-spaces.
   */
  template < typename StateSpace, typename OutputSpace >
  typename boost::disable_if< pp::is_temporal_space< StateSpace >,
                              pp::topology_traits< OutputSpace > >::type::point_type
    map_to_space( const typename pp::topology_traits< StateSpace >::point_type& p, const StateSpace& ss_space,
                  const OutputSpace& ) const {
    typedef typename ss_system_traits< StateSpaceSystem >::input_type Input;
    typedef typename pp::topology_traits< OutputSpace >::point_type OutputTopoPoint;
    return OutputTopoPoint( ss_system->get_output( ss_space, p, Input(), 0.0 ) );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( ss_system );
  };
  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    named_object::load( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( ss_system );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2300016, 1, "state_to_output_map", named_object )
};
};
};

#endif
