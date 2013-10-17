/**
 * \file maximum_likelihood_mapping.hpp
 * 
 * This library provides a class template which can join together a state-vector topology
 * and a covariance matrix topology (see covar_topology.hpp) to create a Gaussian belief-state
 * topology. The distance pseudo-metric used is the symmetric KL-divergence between two 
 * belief-states.
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

#ifndef REAK_GAUSSIAN_BELIEF_SPACE_HPP
#define REAK_GAUSSIAN_BELIEF_SPACE_HPP

#include "base/named_object.hpp"

#include "belief_state_concept.hpp"

namespace ReaK {

namespace ctrl {



/**
 * This class is a bijection mapping from a belief-space to a state-space (topology) 
 * by making the assumption of maximum likelihood. In other words, this mapping reduces
 * belief-states into their most likely value only.
 * 
 * Models: BijectionConcept between a BeliefSpace class and a compatible state topology.
 */
template <typename BeliefSpace>
struct maximum_likelihood_map : public named_object {
  
  typedef maximum_likelihood_map<BeliefSpace> self;
  
  typedef typename pp::topology_traits< BeliefSpace >::point_type belief_state_type;
  typedef belief_state_traits< belief_state_type >::state_type state_type;
  
  maximum_likelihood_map() : named_object() { setName("maximum_likelihood_map"); };
  
  /**
   * This function extracts the most-likely-value from a belief-state.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with the state type that the belief-state type would produce.
   * \param b The belief-state from which the maximum likelihood value is sought.
   * \return The maximum likelihood value of the belief-state.
   */
  template <typename StateSpaceOut>
  state_type map_to_space(const belief_state_type& b, const BeliefSpace&, const StateSpaceOut&) const {
    return state_type(b.get_most_likely_state());
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











