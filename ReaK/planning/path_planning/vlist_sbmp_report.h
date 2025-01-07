/**
 * \file vlist_sbmp_report.h
 *
 * This library defines a sampling-based motion/path planning reporter that prints a list of
 * vertex and their properties.
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

#ifndef REAK_PLANNING_PATH_PLANNING_VLIST_SBMP_REPORT_H_
#define REAK_PLANNING_PATH_PLANNING_VLIST_SBMP_REPORT_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/shared_object.h"
#include "ReaK/topologies/interpolation/seq_path_base.h"
#include "ReaK/topologies/interpolation/seq_trajectory_base.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include <fstream>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and prints all the vertices (the properties) of the motion graph into files.
 * \tparam VertexPrinter A functor type that can print, in a single line, all the properties of a vertex of the
 * motion-graph.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template <typename VertexPrinter, typename NextReporter = no_sbmp_report>
struct vlist_sbmp_report : public shared_object {
  using self = vlist_sbmp_report<VertexPrinter, NextReporter>;

  /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
  NextReporter next_reporter;
  /// Holds the functor that can print, in a single line, all the properties of a vertex of the motion-graph.
  VertexPrinter print_to_stream;
  /// Holds the file-path where to output the reports.
  std::string file_path;

  /**
   * Parametrized constructor.
   * \param aFilePath The path where to create the output files.
   * \param aPrinter The functor that can print, in a single line, all the properties of a vertex of the motion-graph.
   * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
   */
  explicit vlist_sbmp_report(const std::string& aFilePath,
                             VertexPrinter aPrinter = VertexPrinter(),
                             NextReporter aNextReporter = NextReporter())
      : next_reporter(aNextReporter),
        print_to_stream(aPrinter),
        file_path(aFilePath) {}

  vlist_sbmp_report() : vlist_sbmp_report("") {}

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
  template <typename FreeSpaceType, typename MotionGraph, typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g,
                         PositionMap pos) const {
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << num_vertices(g);
    std::ofstream file_out(file_path + "vlist_" + ss.str());

    for (auto v : vertices(g)) {
      print_to_stream(file_out, v, g);
    }

    next_reporter.draw_motion_graph(free_space, g, pos);
  }

  /**
   * Draws the solution trajectory or path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj_or_path The solution trajectory or path.
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
    A& RK_SERIAL_SAVE_WITH_NAME(next_reporter) &
        RK_SERIAL_SAVE_WITH_NAME(print_to_stream) &
        RK_SERIAL_SAVE_WITH_NAME(file_path);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(next_reporter) &
        RK_SERIAL_LOAD_WITH_NAME(print_to_stream) &
        RK_SERIAL_LOAD_WITH_NAME(file_path);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460010, 1, "vlist_sbmp_report",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_VLIST_SBMP_REPORT_H_
