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

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/shared_object.hpp"

#include <cmath>
#include "boost/concept_check.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/tangent_bundle_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK::pp {

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
 * p_out = m.map_to_space(p_in, space_in, space_out);  The point (p_out) in the output space (space_out) can be obtained
 *from a point (p_in) in the input space (space_in).
 *
 * \tparam Mapping The mapping type to be checked for this concept.
 * \tparam InSpace The input space type of the mapping.
 * \tparam OutSpace The output space type of the mapping.
 */
template <typename Mapping, typename InSpace, typename OutSpace>
struct BijectionConcept {

  InSpace space_in;
  OutSpace space_out;
  topology_point_type_t<InSpace> p_in;
  topology_point_type_t<OutSpace> p_out;
  Mapping m;

  BOOST_CONCEPT_ASSERT((TopologyConcept<InSpace>));
  BOOST_CONCEPT_ASSERT((TopologyConcept<OutSpace>));

  BOOST_CONCEPT_USAGE(BijectionConcept) {
    p_out = m.map_to_space(p_in, space_in, space_out);
  }
};

/**
 * This class is a simple composition of two topological maps. Provided an output map,
 * an intermediate space, and an input map, this topological map will first use the
 * input map to map a given source point to the intermediate space, and then use
 * the output map to make the final mapping to the given output space.
 * \note This class template is limited to implementing a BijectionConcept, i.e., a directional mapping, since it cannot
 * perform the inverse mapping.
 * \tparam OuterBijection The type of the map used to produce the output of the bijection.
 * \tparam InnerBijection The type of the map used to map the given point to the intermediate space of the composition.
 * \tparam MiddleSpace The intermediate space type between the two bijections.
 */
template <typename OuterBijection, typename InnerBijection,
          typename MiddleSpace>
struct bijection_cascade : public shared_object {
  using self = bijection_cascade<OuterBijection, InnerBijection, MiddleSpace>;

  OuterBijection map_outer;
  InnerBijection map_inner;
  std::shared_ptr<MiddleSpace> mid_space;

  bijection_cascade(const OuterBijection& aMapOuter,
                    const InnerBijection& aMapInner,
                    const std::shared_ptr<MiddleSpace>& aMidSpace)
      : map_outer(aMapOuter), map_inner(aMapInner), mid_space(aMidSpace) {}

  /**
   * This function template maps a given point in an input space
   * to an output point in the output space. The mapping simply forwards the point.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space.
   * \tparam OutSpace The type of the output space.
   * \param p_in The point in the input space.
   * \param s_in The input space.
   * \param s_out The output space.
   * \return A point in the output space, identical in value and type to the input point.
   */
  template <typename PointType, typename SpaceIn, typename SpaceOut>
  topology_point_type_t<SpaceOut> map_to_space(const PointType& p_in,
                                               const SpaceIn& s_in,
                                               const SpaceOut& s_out) const {
    return map_outer.map_to_space(
        map_inner.map_to_space(p_in, s_in, *mid_space), *mid_space, s_out);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(map_outer) &
        RK_SERIAL_SAVE_WITH_NAME(map_inner) &
        RK_SERIAL_SAVE_WITH_NAME(mid_space);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(map_outer) &
        RK_SERIAL_LOAD_WITH_NAME(map_inner) &
        RK_SERIAL_LOAD_WITH_NAME(mid_space);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002B, 1, "bijection_cascade",
                              shared_object)
};

/**
 * This class is a simple identity topological map, that is, a map that simply
 * forwards the point, regardless of the input versus output spaces provided.
 * \note This class will only be a valid map between spaces with compatible point-types.
 */
struct identity_topo_map : public shared_object {
  using self = identity_topo_map;

  identity_topo_map() = default;

  /**
   * This function template maps a given point in an input space
   * to an output point in the output space. The mapping simply forwards the point.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space.
   * \tparam OutSpace The type of the output space.
   * \param p_in The point in the input space.
   * \return A point in the output space, identical in value and type to the input point.
   */
  template <typename PointType, typename SpaceIn, typename SpaceOut>
  PointType map_to_space(const PointType& p_in, const SpaceIn& /*unused*/,
                         const SpaceOut& /*unused*/) const {
    return p_in;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002D, 1, "identity_topo_map",
                              shared_object)
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
struct HomeomorphismConcept
    : public BijectionConcept<Mapping, InSpace, OutSpace>,
      public BijectionConcept<Mapping, OutSpace, InSpace> {

  BOOST_CONCEPT_USAGE(HomeomorphismConcept) {}
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
struct DiffeomorphismConcept
    : public BijectionConcept<Mapping, InSpace, OutSpace>,
      public BijectionConcept<Mapping, OutSpace, InSpace>,
      public TangentBundleConcept<InSpace, Order, IndependentSpace>,
      public TangentBundleConcept<OutSpace, Order, IndependentSpace> {

  BOOST_CONCEPT_USAGE(DiffeomorphismConcept) {}
};

}  // namespace ReaK::pp

#endif
