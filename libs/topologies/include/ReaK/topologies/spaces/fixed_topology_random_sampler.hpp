/**
 * \file fixed_topology_random_sampler.hpp
 *
 * This library defines the default random-sampler class that work on points of a fixed topology.
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

#ifndef REAK_FIXED_TOPOLOGY_RANDOM_SAMPLER_HPP
#define REAK_FIXED_TOPOLOGY_RANDOM_SAMPLER_HPP

#include <ReaK/core/base/serializable.hpp>

#include "random_sampler_concept.hpp"

namespace ReaK {

namespace pp {

/**
 * This class is the default random-sampler functor which models the RandomSamplerConcept.
 * This class will simply rely on the random_point() function included in the
 * given topology.
 * \note Do not use this random-sampler to define a topology, because it will be cyclic (infinite recursion).
 */
template < typename Topology >
struct fixed_topology_random_sampler : public serializable {
  typedef typename topology_traits< Topology >::point_type PointType;

  const Topology* m_space;

  fixed_topology_random_sampler( const Topology* pSpace = nullptr ) : m_space( pSpace ){};

  template < typename Topology2 >
  PointType operator()( const Topology2& s ) const {
    if( m_space )
      return m_space->random_point();
    else
      return s.random_point();
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( fixed_topology_random_sampler, 0xC2450005, 1, "fixed_topology_random_sampler",
                               serializable )
};
};
};


#endif
