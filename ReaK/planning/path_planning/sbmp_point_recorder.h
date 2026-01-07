/**
 * \file sbmp_point_recorder.h
 *
 * This library defines a sampling-based motion/path planning reporter for recording the
 * points along a motion graph and solution paths.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_PLANNING_PATH_PLANNING_SBMP_POINT_RECORDER_H_
#define REAK_PLANNING_PATH_PLANNING_SBMP_POINT_RECORDER_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/shared_object.h"
#include "ReaK/core/recorders/ascii_recorder.h"
#include "ReaK/core/recorders/data_record.h"
#include "ReaK/planning/path_planning/basic_sbmp_reporters.h"
#include "ReaK/topologies/interpolation/seq_trajectory_base.h"
#include "ReaK/topologies/spaces/topological_map_concepts.h"

#include <tuple>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and uses the recorder library to print out the points of a motion-graph and solution paths
 * as space-separated files (i.e., loadable into Matlab / Octave).
 * \tparam JointStateSpace A joint-state space type.
 * \tparam JointStateMapping A map type that can map the points of the path-planning topology into points of the
 * joint-state space.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template <typename JointStateSpace,
          typename JointStateMapping = identity_topo_map,
          typename NextReporter = no_sbmp_report>
class sbmp_point_recorder : public shared_object {
 public:
  using self =
      sbmp_point_recorder<JointStateSpace, JointStateMapping, NextReporter>;

  /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
  NextReporter next_reporter;

 protected:
  /// A shared-pointer to the joint-state space.
  std::shared_ptr<JointStateSpace> jt_space;
  /// A map that can map the points of the path-planning topology into points of the joint-state space.
  JointStateMapping map_to_jt_space;
  /// Holds the interval-size between output points of the solution trajectory/path.
  double interval_size;
  /// Holds the file-path where to output the reports.
  std::string file_path;
  mutable std::size_t solution_count;

 public:
  /**
   * Parametrized constructor.
   * \param aJointSpace The shared-pointer to the joint-state space.
   * \param aMapToJtSpace The map that can map the points of the path-planning topology into points of the joint-state
   * space.
   * \param aFilePath The path where to create the output files.
   * \param aIntervalSize The interval-size between output points of the solution trajectory/path.
   * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
   */
  explicit sbmp_point_recorder(
      const std::shared_ptr<JointStateSpace>& aJointSpace,
      const JointStateMapping& aMapToJtSpace = JointStateMapping(),
      const std::string& aFilePath = "", double aIntervalSize = 0.1,
      NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        jt_space(aJointSpace),
        map_to_jt_space(aMapToJtSpace),
        interval_size(aIntervalSize),
        file_path(aFilePath),
        solution_count(0) {}

  explicit sbmp_point_recorder(const std::string& aFilePath,
                               double aIntervalSize = 0.1,
                               NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        jt_space(),
        map_to_jt_space(),
        interval_size(aIntervalSize),
        file_path(aFilePath),
        solution_count(0) {}

  sbmp_point_recorder()
      : sbmp_point_recorder(std::shared_ptr<JointStateSpace>()) {}

  void reset_internal_state() { solution_count = 0; }

  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   * \param free_space The C-free topology.
   * \param g The motion-graph.
   * \param pos The position-map to obtain positions of the motion-graph vertices.
   */
  template <typename FreeSpaceType, typename MotionGraph,
            typename SteerRecOrPosMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g,
                         SteerRecOrPosMap steer_rec_or_pos) const {
    using ReaK::to_vect;

    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
    recorder::ascii_recorder rec_out(file_path + "progress_" + ss.str());
    bool not_initialized_yet = true;

    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      for (auto [vi, vi_end] = vertices(g); vi != vi_end; ++vi) {
        for (auto [ei, ei_end] = out_edges(*vi, g); ei != ei_end; ++ei) {
          const auto& st_rec = get(steer_rec_or_pos, *ei);
          for (auto it = st_rec.begin_fraction_travel();
               it != st_rec.end_fraction_travel(); it += 0.1) {
            auto s_new = map_to_jt_space.map_to_space(
                *it, free_space.get_super_space(), *jt_space);
            vect_n<double> v_s_new = to_vect<double>(s_new);

            if (not_initialized_yet) {
              for (std::size_t i = 0; i < v_s_new.size(); ++i) {
                std::stringstream ss2;
                ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
                rec_out << ss2.str();
              }
              rec_out << recorder::data_recorder::end_name_row;
              not_initialized_yet = false;
            }

            for (double i : v_s_new) {
              rec_out << i;
            }
            rec_out << recorder::data_recorder::end_value_row;
          }
        }
      }
    } else {
      for (auto [vi, vi_end] = vertices(g); vi != vi_end; ++vi) {
        auto p_u = get(steer_rec_or_pos, *vi);
        auto s_u = map_to_jt_space.map_to_space(
            p_u, free_space.get_super_space(), *jt_space);
        vect_n<double> v_s_u = to_vect<double>(s_u);

        if (not_initialized_yet) {
          for (int i = 0; i < v_s_u.size(); ++i) {
            std::stringstream ss2;
            ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
            rec_out << ss2.str();
          }
          rec_out << recorder::data_recorder::end_name_row;
          not_initialized_yet = false;
        }
        for (double i : v_s_u) {
          rec_out << i;
        }
        rec_out << recorder::data_recorder::end_value_row;

        for (auto [ei, ei_end] = out_edges(*vi, g); ei != ei_end; ++ei) {
          auto p_v = get(steer_rec_or_pos, target(*ei, g));
          auto s_v = map_to_jt_space.map_to_space(
              p_v, free_space.get_super_space(), *jt_space);
          v_s_u = to_vect<double>(s_u);

          for (double i : v_s_u) {
            rec_out << i;
          }
          rec_out << recorder::data_recorder::end_value_row;

          for (double j = 0.1; j <= 1.01; j += 0.1) {
            auto p_new =
                free_space.get_super_space().move_position_toward(p_u, j, p_v);
            auto s_new = map_to_jt_space.map_to_space(
                p_new, free_space.get_super_space(), *jt_space);
            vect_n<double> v_s_new = to_vect<double>(s_new);
            for (double i : v_s_new) {
              rec_out << i;
            }
            rec_out << recorder::data_recorder::end_value_row;
          }
        }
      }
    }
    rec_out << recorder::data_recorder::flush;
    rec_out << recorder::data_recorder::close;

    next_reporter.draw_motion_graph(free_space, g, steer_rec_or_pos);
  }

  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  auto draw_solution(const FreeSpaceType& free_space,
                     const TrajOrPathPtr& traj_or_path) const {
    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << (solution_count++) << ".ssv";
    recorder::ascii_recorder rec_out(file_path + "solution_" + ss.str());

    auto get_start_it = [&]() {
      if constexpr (is_temporal_space_v<FreeSpaceType>) {
        return traj_or_path->get_start_time();
      } else {
        return traj_or_path->begin_fraction_travel();
      }
    };
    auto deref_it = [&](auto& it) {
      if constexpr (is_temporal_space_v<FreeSpaceType>) {
        return traj_or_path->get_point_at_time(it);
      } else {
        return *it;
      }
    };
    auto get_end_it = [&]() {
      if constexpr (is_temporal_space_v<FreeSpaceType>) {
        return traj_or_path->get_end_time();
      } else {
        return traj_or_path->end_fraction_travel();
      }
    };

    auto it = get_start_it();
    auto u_pt = deref_it(it);
    auto s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(),
                                            *jt_space);
    vect_n<double> v_s_u = to_vect<double>(s_u);

    for (std::size_t i = 0; i < v_s_u.size(); ++i) {
      std::stringstream ss2;
      ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
      rec_out << ss2.str();
    }
    rec_out << recorder::data_recorder::end_name_row;
    for (double i : v_s_u) {
      rec_out << i;
    }
    rec_out << recorder::data_recorder::end_value_row;

    while (it != get_end_it()) {
      it += interval_size;
      u_pt = deref_it(it);
      s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(),
                                         *jt_space);
      v_s_u = to_vect<double>(s_u);
      for (double i : v_s_u) {
        rec_out << i;
      }
      rec_out << recorder::data_recorder::end_value_row;
    }
    rec_out << recorder::data_recorder::flush;
    rec_out << recorder::data_recorder::close;

    next_reporter.draw_solution(free_space, traj_or_path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter);
    if constexpr (!std::is_same_v<JointStateMapping, identity_topo_map>) {
      A& RK_SERIAL_SAVE_WITH_NAME(jt_space) &
          RK_SERIAL_SAVE_WITH_NAME(map_to_jt_space);
    }
    A& RK_SERIAL_SAVE_WITH_NAME(interval_size) &
        RK_SERIAL_SAVE_WITH_NAME(file_path);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter);
    if constexpr (!std::is_same_v<JointStateMapping, identity_topo_map>) {
      A& RK_SERIAL_LOAD_WITH_NAME(jt_space) &
          RK_SERIAL_LOAD_WITH_NAME(map_to_jt_space);
    }
    A& RK_SERIAL_LOAD_WITH_NAME(interval_size) &
        RK_SERIAL_LOAD_WITH_NAME(file_path);
    solution_count = 0;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC246000F, 1, "sbmp_point_recorder",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_SBMP_POINT_RECORDER_H_
