/**
 * \file basic_sbmp_reporters.h
 *
 * This library defines simple sampling-based motion/path planning reporters.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_PLANNING_PATH_PLANNING_SIMPLE_SBMP_REPORTERS_H_
#define REAK_PLANNING_PATH_PLANNING_SIMPLE_SBMP_REPORTERS_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/shared_object.h"
#include "ReaK/topologies/interpolation/seq_path_base.h"
#include "ReaK/topologies/interpolation/seq_trajectory_base.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include <chrono>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and reports nothing.
 */
struct no_sbmp_report : public shared_object {

  void reset_internal_state() {}

  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   */
  template <typename FreeSpaceType, typename MotionGraph,
            typename SteerRecOrPosMap>
  void draw_motion_graph(const FreeSpaceType& /*unused*/,
                         const MotionGraph& /*unused*/,
                         SteerRecOrPosMap /*unused*/) const {}

  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  void draw_solution(const FreeSpaceType& /*unused*/,
                     const TrajOrPathPtr& /*unused*/) const {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(no_sbmp_report, 0xC2460002, 1, "no_sbmp_report",
                              shared_object)
};

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and uses the underlying C-free (free_space) to draw individual edges of the motion graph or
 * solution trajectory. The underlying space should have the functions:
 *
 * void reset_output() const;
 *
 * void draw_edge(const point_type& a, const point_type& b, bool is_solution_path) const;
 *
 * void save_output(const std::string& filename) const;
 */
template <typename NextReporter = no_sbmp_report>
struct differ_sbmp_report_to_space : public shared_object {
  using self = differ_sbmp_report_to_space<NextReporter>;

  NextReporter next_reporter;
  /// Holds the interval-size between output points of the solution trajectory/path.
  double interval_size;
  /// Holds the file-path where to output the reports.
  std::string file_path;

  mutable std::size_t progress_count;
  mutable std::size_t solution_count;

  explicit differ_sbmp_report_to_space(
      const std::string& aFilePath, double aIntervalSize = 0.1,
      NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        interval_size(aIntervalSize),
        file_path(aFilePath),
        progress_count(0),
        solution_count(0) {}

  differ_sbmp_report_to_space() : differ_sbmp_report_to_space("") {}

  void reset_internal_state() {
    progress_count = 0;
    solution_count = 0;

    next_reporter.reset_internal_state();
  }

  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam SteerRecOrPosMap The property-map type that can map motion-graph edge descriptors into steer-records or position.
   * \param free_space The C-free topology.
   * \param g The motion-graph.
   * \param steer_rec_or_pos The steer-record-map to obtain steer-records of the motion-graph edges,
   *                         or the position-map to obtain positions of the motion-graph vertices.
   */
  template <typename FreeSpaceType, typename MotionGraph,
            typename SteerRecOrPosMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g,
                         SteerRecOrPosMap steer_rec_or_pos) const {
    free_space.reset_output();

    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      using PointType = topology_point_type_t<FreeSpaceType>;
      for (auto v : vertices(g)) {
        for (auto e : out_edges(v, g)) {
          const auto& st_rec = get(steer_rec_or_pos, e);
          auto it = st_rec.begin_fraction_travel();
          auto prev_it = it;
          it += 0.1;
          for (; prev_it != st_rec.end_fraction_travel(); it += 0.1) {
            free_space.draw_edge(PointType(*prev_it), PointType(*it), false);
            prev_it = it;
          }
        }
      }
    } else {
      for (auto v : vertices(g)) {
        for (auto e : out_edges(v, g)) {
          free_space.draw_edge(get(steer_rec_or_pos, v),
                               get(steer_rec_or_pos, target(e, g)), false);
        }
      }
    }

    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << (progress_count++);
    free_space.save_output(file_path + "progress_" + ss.str());

    next_reporter.draw_motion_graph(free_space, g, steer_rec_or_pos);
  }

  /**
   * Draws the solution.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj_or_path The solution.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  void draw_solution(const FreeSpaceType& free_space,
                     const TrajOrPathPtr& traj_or_path) const {
    free_space.reset_output();

    double total_dist = 0.0;

    if constexpr (is_temporal_space_v<FreeSpaceType>) {
      double t = traj_or_path->get_start_time();
      auto u_pt = traj_or_path->get_point_at_time(t);
      auto v_pt = u_pt;
      while (t < traj_or_path->get_end_time()) {
        t += interval_size;
        v_pt = traj_or_path->get_point_at_time(t);
        free_space.draw_edge(u_pt, v_pt, true);
        u_pt = v_pt;
      }
      total_dist =
          (traj_or_path->get_end_time() - traj_or_path->get_start_time());
    } else {
      auto u_it = traj_or_path->begin_distance_travel();
      auto v_it = u_it;
      v_it += interval_size;
      while (u_it != traj_or_path->end_distance_travel()) {
        total_dist += interval_size;
        free_space.draw_edge(*u_it, *v_it, true);
        u_it = v_it;
        v_it += interval_size;
      }
    }

    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << (solution_count++) << "_"
       << total_dist;
    free_space.save_output(file_path + "solution_" + ss.str());

    next_reporter.draw_solution(free_space, traj_or_path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter) &
        RK_SERIAL_SAVE_WITH_NAME(interval_size) &
        RK_SERIAL_SAVE_WITH_NAME(file_path);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter) &
        RK_SERIAL_LOAD_WITH_NAME(interval_size) &
        RK_SERIAL_LOAD_WITH_NAME(file_path);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460003, 1,
                              "differ_sbmp_report_to_space", shared_object)
};

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and simply times the progress of the planner, which it outputs to a given output stream.
 */
template <typename NextReporter = no_sbmp_report>
struct timing_sbmp_report : public shared_object {
  using self = timing_sbmp_report<NextReporter>;

  NextReporter next_reporter;
  /// Holds the interval-size between output points of the solution trajectory/path.
  mutable std::chrono::high_resolution_clock::time_point last_time;
  std::ostream* p_out;

  explicit timing_sbmp_report(std::ostream& aOutStream,
                              NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        last_time(std::chrono::high_resolution_clock::now()),
        p_out(&aOutStream) {}

  timing_sbmp_report() : timing_sbmp_report(std::cout) {}

  void reset_internal_state() {
    last_time = std::chrono::high_resolution_clock::now();

    next_reporter.reset_internal_state();
  }

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
    std::chrono::high_resolution_clock::duration dt =
        std::chrono::high_resolution_clock::now() - last_time;
    (*p_out)
        << num_vertices(g) << " "
        << std::chrono::duration_cast<std::chrono::microseconds>(dt).count()
        << std::endl;

    next_reporter.draw_motion_graph(free_space, g, steer_rec_or_pos);

    last_time = std::chrono::high_resolution_clock::now();
  }

  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  void draw_solution(const FreeSpaceType& free_space,
                     const TrajOrPathPtr& traj_or_path) const {
    next_reporter.draw_solution(free_space, traj_or_path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter);
    last_time = std::chrono::high_resolution_clock::now();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460004, 1, "timing_sbmp_report",
                              shared_object)
};

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and simply prints the progress of the planner.
 */
template <typename NextReporter = no_sbmp_report>
struct print_sbmp_progress : public shared_object {
  using self = print_sbmp_progress<NextReporter>;

  NextReporter next_reporter;

  explicit print_sbmp_progress(NextReporter aNextReporter)
      : next_reporter(aNextReporter) {}

  print_sbmp_progress() : print_sbmp_progress(NextReporter()) {}

  void reset_internal_state() { next_reporter.reset_internal_state(); }

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
    std::cout << "\r" << std::setw(15) << num_vertices(g);
    std::cout.flush();

    next_reporter.draw_motion_graph(free_space, g, steer_rec_or_pos);
  }

  /**
   * Draws the solution.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj_or_path The solution.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  void draw_solution(const FreeSpaceType& free_space,
                     const TrajOrPathPtr& traj_or_path) const {
    std::cout << "Solution Found!" << std::endl;

    next_reporter.draw_solution(free_space, traj_or_path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460005, 1, "print_sbmp_progress",
                              shared_object)
};

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and records the current best solution for each progress interval, which it outputs to a given output stream.
 */
template <typename NextReporter = no_sbmp_report>
struct least_cost_sbmp_report : public shared_object {
  using self = least_cost_sbmp_report<NextReporter>;

  NextReporter next_reporter;
  std::ostream* p_out;
  std::ostream* p_sol;
  mutable double current_best;
  mutable std::size_t last_node_count;

  explicit least_cost_sbmp_report(std::ostream& aOutStream,
                                  std::ostream* aPSolutionOutput = nullptr,
                                  NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        p_out(&aOutStream),
        p_sol(aPSolutionOutput),
        current_best(1e10),
        last_node_count(0) {}

  least_cost_sbmp_report() : least_cost_sbmp_report(std::cout) {}

  void reset_internal_state() {
    current_best = 1e10;
    last_node_count = 0;

    next_reporter.reset_internal_state();
  }

  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam SteerRecOrPosMap The property-map type.
   * \param free_space The C-free topology.
   * \param g The motion-graph.
   * \param pos The map to obtain position or steer record.
   */
  template <typename FreeSpaceType, typename MotionGraph,
            typename SteerRecOrPosMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g,
                         SteerRecOrPosMap steer_rec_or_pos) const {
    last_node_count = num_vertices(g);
    (*p_out) << num_vertices(g) << " " << current_best << std::endl;

    next_reporter.draw_motion_graph(free_space, g, steer_rec_or_pos);
  }

  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType, typename TrajOrPathPtr>
  void draw_solution(const FreeSpaceType& free_space,
                     const TrajOrPathPtr& traj_or_path) const {
    double total_cost = 0.0;
    if constexpr (is_temporal_space_v<FreeSpaceType>) {
      total_cost =
          traj_or_path->travel_distance(*(traj_or_path->begin_time_travel()),
                                        *(traj_or_path->end_time_travel()));
    } else {
      for (auto it = traj_or_path->begin_fraction_travel();
           it != traj_or_path->end_fraction_travel();) {
        auto it_next = it;
        it_next += 1.0;
        total_cost += get(distance_metric, free_space.get_super_space())(
            *it, *it_next, free_space.get_super_space());
        it = it_next;
      }
    }
    if (total_cost < current_best) {
      current_best = total_cost;
    }

    if (p_sol) {
      (*p_sol) << last_node_count << " " << current_best << std::endl;
    }

    next_reporter.draw_solution(free_space, traj_or_path);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter);
    current_best = 1e10;
    last_node_count = 0;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460006, 1, "least_cost_sbmp_report",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_SIMPLE_SBMP_REPORTERS_H_
