/**
 * \file any_sbmp_reporter.h
 *
 * This library defines a type-erasure base-class for sampling-based motion/path planning reporters.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_PLANNING_PATH_PLANNING_ANY_SBMP_REPORTER_H_
#define REAK_PLANNING_PATH_PLANNING_ANY_SBMP_REPORTER_H_

#include "ReaK/core/base/shared_object.h"
#include "ReaK/planning/path_planning/any_motion_graphs.h"
#include "ReaK/planning/path_planning/planning_queries.h"
#include "ReaK/topologies/interpolation/seq_path_base.h"
#include "ReaK/topologies/interpolation/seq_trajectory_base.h"
#include "ReaK/topologies/spaces/subspace_concept.h"
#include "any_motion_graphs.h"

#include <functional>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class can be used as the base for a dynamically polymorphic SBMP/SBPP Reporter
 * (SBMPReporterConcept and SBPPReporterConcept). This operates on type-erasure via the ReaK::graph::any_graph class.
 * \tparam FreeSpaceType The C-free topology type.
 */
template <typename FreeSpaceType>
class any_sbmp_reporter : public shared_object {
 public:
  using self = any_sbmp_reporter<FreeSpaceType>;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  using wrapped = std::reference_wrapper<const self>;

  using solution_record_ptr =
      typename planning_query<FreeSpaceType>::solution_record_ptr;

  virtual void draw_any_motion_graph(
      const FreeSpaceType& /*unused*/,
      const bagl::dynamic_graph_observer& /*unused*/) const {}
  virtual void draw_any_solution(const FreeSpaceType& /*unused*/,
                                 const solution_record_ptr& /*unused*/) const {}

  ~any_sbmp_reporter() override = default;

  virtual void reset_internal_state() {}

  /**
   * Draws the entire motion-graph.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   */
  template <typename MotionGraph, typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& space, const MotionGraph& g,
                         PositionMap /*unused*/) const {
    bagl::dynamic_properties dp;
    add_all_motion_graph_property_maps<FreeSpaceType>(g, dp);
    bagl::dynamic_graph_observer_wrapper<const MotionGraph> mg(g, dp);
    this->draw_any_motion_graph(space, mg);
  }

  /**
   * Draws the solution trajectory.
   */
  void draw_solution(const FreeSpaceType& space,
                     const solution_record_ptr& path) const {
    this->draw_any_solution(space, path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::get_static_object_type()->version());
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460013, 1, "any_sbmp_reporter",
                              shared_object)
};

namespace detail {
namespace {

template <bool IsSteerable, typename FreeSpaceType>
struct get_sbmp_reporter_any_property_type {
  using type = topology_point_type_t<FreeSpaceType>;
  static std::string name() { return "vertex_position"; }
};

template <typename FreeSpaceType>
struct get_sbmp_reporter_any_property_type<true, FreeSpaceType> {
  using type = steerable_space_steer_record_t<FreeSpaceType>;
  static std::string name() { return "edge_steer_rec"; }
};

}  // namespace
}  // namespace detail

/**
 * This class can be used to wrap a SBMP Reporter into a dynamically polymorphic SBMP/SBPP Reporter
 * (SBMPReporterConcept and SBPPReporterConcept). This operates on type-erasure via the ReaK::graph::any_graph class.
 * \tparam FreeSpaceType The C-free topology type.
 * \tparam Reporter The reporter object type to be encapsulated by this type-erasure class.
 */
template <typename FreeSpaceType, typename Reporter>
class type_erased_sbmp_reporter : public any_sbmp_reporter<FreeSpaceType> {
 public:
  using base_type = any_sbmp_reporter<FreeSpaceType>;
  using self = type_erased_sbmp_reporter<FreeSpaceType, Reporter>;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  using solution_record_ptr = typename base_type::solution_record_ptr;

 protected:
  Reporter reporter;

 public:
  void draw_any_motion_graph(
      const FreeSpaceType& space,
      const bagl::dynamic_graph_observer& g) const override {
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      reporter.draw_motion_graph(
          space, g,
          bagl::get_dynamic_property_map<
              const steerable_space_steer_record_t<FreeSpaceType>&>(
              "edge_steer_rec", g.get_properties()));
    } else {
      reporter.draw_motion_graph(
          space, g,
          bagl::get_dynamic_property_map<
              const topology_point_type_t<FreeSpaceType>&>("vertex_position",
                                                           g.get_properties()));
    }
  }
  void draw_any_solution(const FreeSpaceType& space,
                         const solution_record_ptr& traj) const override {
    reporter.draw_solution(space, traj);
  }

  void reset_internal_state() override { reporter.reset_internal_state(); }

  explicit type_erased_sbmp_reporter(Reporter aReporter)
      : reporter(aReporter) {}

  type_erased_sbmp_reporter() : type_erased_sbmp_reporter(Reporter()) {}

  ~type_erased_sbmp_reporter() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(reporter);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(reporter);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460014, 1, "type_erased_sbmp_reporter",
                              base_type)
};

template <typename FreeSpaceType, typename Reporter>
class type_erased_sbmp_reporter<FreeSpaceType, Reporter*>
    : public any_sbmp_reporter<FreeSpaceType> {
 public:
  using base_type = any_sbmp_reporter<FreeSpaceType>;
  using self = type_erased_sbmp_reporter<FreeSpaceType, Reporter*>;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  using solution_record_ptr = typename base_type::solution_record_ptr;

 protected:
  Reporter* reporter;

 public:
  void draw_any_motion_graph(
      const FreeSpaceType& space,
      const bagl::dynamic_graph_observer& g) const override {
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      reporter->draw_motion_graph(
          space, g,
          bagl::get_dynamic_property_map<
              const steerable_space_steer_record_t<FreeSpaceType>&>(
              "edge_steer_rec", g.get_properties()));
    } else {
      reporter->draw_motion_graph(
          space, g,
          bagl::get_dynamic_property_map<
              const topology_point_type_t<FreeSpaceType>&>("vertex_position",
                                                           g.get_properties()));
    }
  }
  void draw_any_solution(const FreeSpaceType& space,
                         const solution_record_ptr& traj) const override {
    reporter->draw_solution(space, traj);
  }

  void reset_internal_state() override { reporter->reset_internal_state(); }

  explicit type_erased_sbmp_reporter(Reporter* aReporter)
      : reporter(aReporter) {}

  ~type_erased_sbmp_reporter() override = default;
};

/**
 * This class can be used as the base for a dynamically polymorphic SBMP/SBPP Reporter
 * (SBMPReporterConcept and SBPPReporterConcept). This operates on type-erasure via the ReaK::graph::any_graph class.
 * \tparam FreeSpaceType The C-free topology type.
 */
template <typename FreeSpaceType>
class any_sbmp_reporter_chain : public shared_object {
 public:
  using self = any_sbmp_reporter_chain<FreeSpaceType>;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  using solution_record_ptr =
      typename planning_query<FreeSpaceType>::solution_record_ptr;

 private:
  std::vector<std::shared_ptr<any_sbmp_reporter<FreeSpaceType>>> reporters;

 public:
  /**
   * Add a reporter to this collection of dynamically-dispatched (type-erased) reporters.
   * \tparam Reporter The reporter type to use to report the progress of the path-planning, should model
   * SBMPReporterConcept and SBPPReporterConcept.
   * \param aReporter The reporter object to use to report the progress of the path-planner.
   */
  template <typename Reporter>
  void add_reporter(Reporter aReporter) {
    reporters.emplace_back(
        std::make_shared<type_erased_sbmp_reporter<FreeSpaceType, Reporter>>(
            aReporter));
  }

  template <typename Reporter>
  void add_reporter(const std::reference_wrapper<Reporter>& aReporter) {
    reporters.emplace_back(
        std::make_shared<type_erased_sbmp_reporter<FreeSpaceType, Reporter*>>(
            &aReporter.get()));
  }

  /**
   * Draws the entire motion-graph.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   */
  template <typename MotionGraph, typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& space, const MotionGraph& g,
                         PositionMap /*unused*/) const {
    bagl::dynamic_properties dp;
    add_all_motion_graph_property_maps<FreeSpaceType>(g, dp);
    bagl::dynamic_graph_observer_wrapper<const MotionGraph> mg(g, dp);
    for (auto& r : reporters) {
      r->draw_any_motion_graph(space, mg);
    }
  }

  /**
   * Draws the solution trajectory.
   */
  void draw_solution(const FreeSpaceType& space,
                     const solution_record_ptr& path) const {
    for (auto& r : reporters) {
      r->draw_any_solution(space, path);
    }
  }

  void reset_internal_state() {
    for (auto& r : reporters) {
      r->reset_internal_state();
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(reporters);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(reporters);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460015, 1, "any_sbmp_reporter_chain",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_ANY_SBMP_REPORTER_H_
