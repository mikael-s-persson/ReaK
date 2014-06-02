/**
 * \file se3_random_samplers.hpp
 * 
 * This library defines SE(3) random-sampler classes that work on points of an SE(3) topology.
 * All the classes satisfy the RandomSamplerConcept.
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

#ifndef REAK_SE3_RANDOM_SAMPLERS_HPP
#define REAK_SE3_RANDOM_SAMPLERS_HPP

#include <ReaK/core/base/serializable.hpp>

#include <ReaK/ctrl/path_planning/random_sampler_concept.hpp>

#include "se3_topologies.hpp"

namespace ReaK {

namespace pp {
  
/**
 */
struct position_only_sampler : public serialization::serializable {
  
  position_only_sampler() { };
  
  template <typename Topology>
  position_only_sampler(const Topology&) { };
  
  /** 
   * This function returns a random sample point on a topology.
   * \tparam SE3Topology The SE(3) topology.
   * \param s The topology or space on which the points lie.
   * \return A random sample point on the given topology.
   */
  template <typename SE3Topology>
  typename topology_traits< SE3Topology >::point_type operator()(const SE3Topology& s) const {
    typedef typename topology_traits< SE3Topology >::point_type point_type;
    point_type result = s.origin();
    point_type rnd = s.random_point();
    set_position(result, get_position(rnd));
    return result;
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(position_only_sampler,0xC2450006,1,"position_only_sampler",serialization::serializable)
};

};

};


#endif


