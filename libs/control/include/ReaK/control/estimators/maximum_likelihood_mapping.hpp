/**
 * \file maximum_likelihood_mapping.hpp
 * 
 * This library provides a class template which can map belief-state points into points of a 
 * compatibly state-space topology. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_MAXIMUM_LIKELIHOOD_MAPPING_HPP
#define REAK_MAXIMUM_LIKELIHOOD_MAPPING_HPP

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include "belief_state_concept.hpp"

#include <boost/utility.hpp>

namespace ReaK {

namespace ctrl {



/**
 * This class is a bijection mapping from a belief-space to a state-space (topology) 
 * by making the assumption of maximum likelihood. In other words, this mapping reduces
 * belief-states into their most likely value only.
 * 
 * Models: BijectionConcept between a belief-state topology and a compatible state topology.
 */
struct maximum_likelihood_map : public named_object {
  
  typedef maximum_likelihood_map self;
  
  maximum_likelihood_map() : named_object() { setName("maximum_likelihood_map"); };
  
  /**
   * This function extracts the most-likely-value from a belief-state.
   * \tparam BeliefSpace A belief-state topoloogy whose state-space is compatible with the given state-space topology.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with the state type that the belief-state type would produce.
   * \param b The belief-state from which the maximum likelihood value is sought.
   * \return The maximum likelihood value of the belief-state.
   */
  template <typename BeliefSpace, typename StateSpaceOut>
  typename boost::enable_if<
    pp::is_temporal_space<BeliefSpace>,
    pp::topology_traits<StateSpaceOut> >::type::point_type
      map_to_space(const typename pp::topology_traits<BeliefSpace>::point_type& b, const BeliefSpace&, const StateSpaceOut&) const {
    typedef typename pp::topology_traits<StateSpaceOut>::point_type OutPointType;
    return OutPointType(b.time, b.pt.get_most_likely_state());
  };
  
  /**
   * This function extracts the most-likely-value from a belief-state.
   * \tparam BeliefSpace A belief-state topoloogy whose state-space is compatible with the given state-space topology.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with the state type that the belief-state type would produce.
   * \param b The belief-state from which the maximum likelihood value is sought.
   * \return The maximum likelihood value of the belief-state.
   */
  template <typename BeliefSpace, typename StateSpaceOut>
  typename boost::disable_if<
    pp::is_temporal_space<BeliefSpace>,
    pp::topology_traits<StateSpaceOut> >::type::point_type
      map_to_space(const typename pp::topology_traits<BeliefSpace>::point_type& b, const BeliefSpace&, const StateSpaceOut&) const {
    typedef typename pp::topology_traits<StateSpaceOut>::point_type OutPointType;
    return OutPointType(b.get_most_likely_state());
  };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300014,1,"maximum_likelihood_map",named_object)
    
};



};

};

#endif











