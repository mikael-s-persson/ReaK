/**
 * \file topological_map_concepts.hpp
 * 
 * This library defines the traits and concepts that pertain to mappings between metric spaces
 * or topologies in general. This includes the basic concept of a bijection, a homeomorphism,
 * and a diffeomorphism. Other related mappings cannot be enforced by concepts (would require axioms).
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

#ifndef REAK_TOPOLOGICAL_MAP_CONCEPTS_HPP
#define REAK_TOPOLOGICAL_MAP_CONCEPTS_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"


#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"
#include "tangent_bundle_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {
  
  
/**
 * This concept defines the requirements to fulfill in order to model a bijection between 
 * two metric-spaces as used in ReaK::pp.
 * 
 * Required Concepts:
 * 
 * Both spaces should model the TopologyConcept.
 * 
 * Valid expressions (m of type Mapping):
 * 
 * p_out = m.map_to_space(p_in, space_in, space_out);  The point (p_out) in the output space (space_out) can be obtained from a point (p_in) in the input space (space_in).
 * 
 * \tparam Mapping The mapping type to be checked for this concept.
 * \tparam InSpace The input space type of the mapping.
 * \tparam OutSpace The output space type of the mapping.
 */
template <typename Mapping, typename InSpace, typename OutSpace>
struct BijectionConcept :
  public MetricSpaceConcept< InSpace >,
  public MetricSpaceConcept< OutSpace > {
  
  InSpace space_in;
  OutSpace space_out;
  typename topology_traits<InSpace>::point_type p_in;
  typename topology_traits<OutSpace>::point_type p_out;
  Mapping m;
  
  BOOST_CONCEPT_ASSERT((TopologyConcept<InSpace>));
  BOOST_CONCEPT_ASSERT((TopologyConcept<OutSpace>));
  
  BOOST_CONCEPT_USAGE(BijectionConcept) 
  {
    p_out = m.map_to_space(p_in, space_in, space_out);
  };
  
};


template <typename OuterBijection, typename InnerBijection, typename MiddleSpace>
struct bijection_cascade : public shared_object {
  typedef bijection_cascade<OuterBijection, InnerBijection, MiddleSpace> self;
  
  OuterBijection map_outer;
  InnerBijection map_inner;
  shared_ptr< MiddleSpace > mid_space;
  
  bijection_cascade(const OuterBijection& aMapOuter, 
		    const InnerBijection& aMapInner, 
		    const shared_ptr< MiddleSpace >& aMidSpace) : 
		    map_outer(aMapOuter), 
		    map_inner(aMapInner),
		    mid_space(aMidSpace) { };
  
  template <typename PointType, typename SpaceIn, typename SpaceOut>
  typename topology_traits<SpaceOut>::point_type map_to_space(const PointType& p_in, const SpaceIn& s_in, const SpaceOut& s_out) const {
    return map_outer.map_to_space( map_inner.map_to_space(p_in, s_in, *mid_space), *mid_space, s_out);
  };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_SAVE_WITH_NAME(map_outer)
      & RK_SERIAL_SAVE_WITH_NAME(map_inner)
      & RK_SERIAL_SAVE_WITH_NAME(mid_space);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_LOAD_WITH_NAME(map_outer)
      & RK_SERIAL_LOAD_WITH_NAME(map_inner)
      & RK_SERIAL_LOAD_WITH_NAME(mid_space);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240002B,1,"bijection_cascade",shared_object)
    
};


struct identity_topo_map : public shared_object {
  typedef identity_topo_map self;
  
  identity_topo_map() { };
  
  template <typename PointType, typename SpaceIn, typename SpaceOut>
  PointType map_to_space(const PointType& p_in, const SpaceIn&, const SpaceOut&) const {
    return p_in;
  };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240002B,1,"identity_topo_map",shared_object)
    
};


/**
 * This concept defines the requirements to fulfill in order to model a homeomorphism between 
 * two metric-spaces as used in ReaK::pp.
 * 
 * Required Concepts:
 * 
 * The mapping is a bijection between both spaces in both directions.
 * 
 * \tparam Mapping The mapping type to be checked for this concept.
 * \tparam InSpace The input space type of the mapping.
 * \tparam OutSpace The output space type of the mapping.
 */
template <typename Mapping, typename InSpace, typename OutSpace>
struct HomeomorphismConcept :
  public BijectionConcept< Mapping, InSpace, OutSpace >,
  public BijectionConcept< Mapping, OutSpace, InSpace > {
  
  BOOST_CONCEPT_USAGE(HomeomorphismConcept) 
  { };
  
};



/**
 * This concept defines the requirements to fulfill in order to model a diffeomorphism between 
 * two differentiable-spaces as used in ReaK::pp.
 * 
 * Required Concepts:
 * 
 * The mapping is a bijection between both spaces in both directions.
 *
 * Both spaces are differentiable by the same independent space and to the given order.
 * 
 * \tparam Mapping The mapping type to be checked for this concept.
 * \tparam InSpace The input space type of the mapping.
 * \tparam OutSpace The output space type of the mapping.
 * \tparam Order The order of differentiation required by the diffeomorphism.
 * \tparam IndependentSpace The independent space against which the differentiation is taken on the spaces.
 */
template <typename Mapping, typename InSpace, typename OutSpace,
          unsigned int Order, typename IndependentSpace>
struct DiffeomorphismConcept :
  public BijectionConcept< Mapping, InSpace, OutSpace >,
  public BijectionConcept< Mapping, OutSpace, InSpace >,
  public TangentBundleConcept< InSpace, Order, IndependentSpace>,
  public TangentBundleConcept< OutSpace, Order, IndependentSpace> {
  
  BOOST_CONCEPT_USAGE(DiffeomorphismConcept) 
  { };
  
};


};

};


#endif


