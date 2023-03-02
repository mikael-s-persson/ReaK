/**
 * \file proxy_model_updater.hpp
 *
 * This library defines an abstract base class for to "update" a proximity-query model for a given time (for
 * dynamic models). A typical derived class implementation would be to query a spatial-trajectory (see
 *SpatialTrajectoryConcept)
 * for the state at the given time and then apply that state to the geometry used by the proximity-query method.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_PROXY_MODEL_UPDATER_HPP
#define REAK_PROXY_MODEL_UPDATER_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include "metric_space_concept.hpp"
#include "topological_map_concepts.hpp"

namespace ReaK::pp {

/**
 * This base-class is used to update dynamic proximity-query models to a given time value.
 */
class proxy_model_updater : public shared_object {
 public:
  using self = proxy_model_updater;

  /**
   * Default constructor.
   * \note This is an abstract base-class, so, it can't be constructed.
   */
  proxy_model_updater() = default;

  /**
   * This virtual function is meant to update the configuration of some (proximity-query) model
   * such that it matches its configuration for the given time (e.g., along its trajectory).
   * \param t The time to which the model should be synchronized.
   */
  virtual void synchronize_proxy_model(double t) const = 0;

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2400029, 1, "proxy_model_updater",
                              shared_object)
};

/**
 * This base-class is used to apply a configuration to a (proximity-query) model.
 */
template <typename JointSpace>
class proxy_model_applicator : public shared_object {
 public:
  using self = proxy_model_applicator<JointSpace>;
  using joint_space = JointSpace;
  using point_type = topology_point_type_t<JointSpace>;

  proxy_model_applicator() = default;

  /**
   * This function applies the given configuration / joint-state onto some underlying (proximity-query) model.
   * \param pt The point in the joint-space, i.e. the joint coordinates / configuration.
   * \param jt_space The joint-space.
   */
  virtual void apply_to_model(const point_type& pt,
                              const joint_space& jt_space) const = 0;

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC240003B, 1, "proxy_model_applicator",
                              shared_object)
};

/**
 * This type-erasure derived-class is used to apply a configuration to a (proximity-query) model.
 * \tparam JointSpace The joint-space of the configurations given to the applicator.
 * \tparam DKTopoMap The generic direct-kinematics topological-map type that is used to apply the configuration to the
 * model.
 * \tparam JointSpaceMapping The map type that can map points between joint-space and intermediate-space.
 * \tparam IntermediateSpace The type of an intermediate space required by the direct-kinematics map.
 */
template <typename JointSpace, typename DKTopoMap,
          typename JointSpaceMapping = identity_topo_map,
          typename IntermediateSpace = JointSpace>
class any_model_applicator : public proxy_model_applicator<JointSpace> {
 public:
  using base_type = proxy_model_applicator<JointSpace>;
  using self = any_model_applicator<JointSpace, DKTopoMap, JointSpaceMapping,
                                    IntermediateSpace>;
  using joint_space = JointSpace;
  using point_type = topology_point_type_t<JointSpace>;

  DKTopoMap dk_topomap;
  JointSpaceMapping map_to_jt_space;
  std::shared_ptr<IntermediateSpace> inter_space;

  /**
   * Parametrized / default constructor.
   * \param aDKTopoMap The generic direct-kinematics topological-map object that is used to apply the configuration to
   * the model.
   * \param aMapToJtSpace The map that can map the points of the joint-space topology into points of the intermediate
   * space.
   * \param aInterSpace The shared-pointer to the intermediate space.
   */
  explicit any_model_applicator(
      const DKTopoMap& aDKTopoMap = {},
      const JointSpaceMapping& aMapToJtSpace = {},
      const std::shared_ptr<IntermediateSpace>& aInterSpace = {})
      : dk_topomap(aDKTopoMap),
        map_to_jt_space(aMapToJtSpace),
        inter_space(aInterSpace) {}

  /**
   * This function applies the given configuration / joint-state onto some underlying (proximity-query) model.
   * \param pt The point in the joint-space, i.e. the joint coordinates / configuration.
   * \param jt_space The joint-space.
   */
  void apply_to_model(const point_type& pt,
                      const joint_space& jt_space) const override {
    if (!inter_space) {
      return;
    }
    dk_topomap.apply_to_model(
        map_to_jt_space.map_to_space(pt, jt_space, *inter_space), *inter_space);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(dk_topomap) &
        RK_SERIAL_SAVE_WITH_NAME(map_to_jt_space) &
        RK_SERIAL_SAVE_WITH_NAME(inter_space);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(dk_topomap) &
        RK_SERIAL_LOAD_WITH_NAME(map_to_jt_space) &
        RK_SERIAL_LOAD_WITH_NAME(inter_space);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240003C, 1, "any_model_applicator",
                              base_type)
};

/**
 * This type-erasure derived-class is used to apply a configuration to a (proximity-query) model.
 * \tparam JointSpace The joint-space of the configurations given to the applicator.
 * \tparam DKTopoMap The generic direct-kinematics topological-map type that is used to apply the configuration to the
 * model.
 */
template <typename JointSpace, typename DKTopoMap>
class any_model_applicator<JointSpace, DKTopoMap, identity_topo_map, JointSpace>
    : public proxy_model_applicator<JointSpace> {
 public:
  using base_type = proxy_model_applicator<JointSpace>;
  using self = any_model_applicator<JointSpace, DKTopoMap>;
  using joint_space = JointSpace;
  using point_type = topology_point_type_t<JointSpace>;

  DKTopoMap dk_topomap;

  /**
   * Parametrized / default constructor.
   * \param aDKTopoMap The generic direct-kinematics topological-map object that is used to apply the configuration to
   * the model.
   */
  explicit any_model_applicator(const DKTopoMap& aDKTopoMap = {})
      : dk_topomap(aDKTopoMap) {}

  /**
   * This function applies the given configuration / joint-state onto some underlying (proximity-query) model.
   * \param pt The point in the joint-space, i.e. the joint coordinates / configuration.
   * \param jt_space The joint-space.
   */
  void apply_to_model(const point_type& pt,
                      const joint_space& jt_space) const override {
    dk_topomap.apply_to_model(pt, jt_space);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(dk_topomap);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(dk_topomap);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240003C, 1, "any_model_applicator",
                              base_type)
};

/**
 * This function template is used to create a type-erasure object that can be used
 * to apply a configuration to a (proximity-query) model.
 * \tparam JointSpace The joint-space of the configurations given to the applicator.
 * \tparam DKTopoMap The generic direct-kinematics topological-map type that is used to apply the configuration to the
 * model.
 * \param aDKMap The generic direct-kinematics topological-map object that is used to apply the configuration to the
 * model.
 * \return A shared-pointer to a type-erased object that can be used to apply a configuration to a (proximity-query)
 * model.
 */
template <typename JointSpace, typename DKTopoMap>
std::shared_ptr<any_model_applicator<JointSpace, DKTopoMap>>
make_any_model_applicator(const DKTopoMap& aDKMap) {
  return std::make_shared<any_model_applicator<JointSpace, DKTopoMap>>(aDKMap);
}

/**
 * This function template is used to create a type-erasure object that can be used
 * to apply a configuration to a (proximity-query) model.
 * \tparam JointSpace The joint-space of the configurations given to the applicator.
 * \tparam DKTopoMap The generic direct-kinematics topological-map type that is used to apply the configuration to the
 * model.
 * \tparam JointSpaceMapping The map type that can map points between joint-space and intermediate-space.
 * \tparam IntermediateSpace The type of an intermediate space required by the direct-kinematics map.
 * \param aDKMap The generic direct-kinematics topological-map object that is used to apply the configuration to the
 * model.
 * \param aMapToJtSpace The map that can map the points of the joint-space topology into points of the intermediate
 * space.
 * \param aInterSpace The shared-pointer to the intermediate space.
 * \return A shared-pointer to a type-erased object that can be used to apply a configuration to a (proximity-query)
 * model.
 */
template <typename JointSpace, typename DKTopoMap, typename JointSpaceMapping,
          typename IntermediateSpace>
std::shared_ptr<any_model_applicator<JointSpace, DKTopoMap, JointSpaceMapping,
                                     IntermediateSpace>>
make_any_model_applicator(
    const DKTopoMap& aDKMap, const JointSpaceMapping& aMapToJtSpace,
    const std::shared_ptr<IntermediateSpace>& aInterSpace) {
  return std::make_shared<any_model_applicator<
      JointSpace, DKTopoMap, JointSpaceMapping, IntermediateSpace>>(
      aDKMap, aMapToJtSpace, aInterSpace);
}

}  // namespace ReaK::pp

#endif
