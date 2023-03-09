
/*
 *    Copyright 2023 Sven Mikael Persson
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
#include "ReaK/planning/path_planning/planning_space_options_po.h"
#include "ReaK/planning/path_planning/planning_space_options.h"

#include "ReaK/core/serialization/archiver_factory.h"

#include "absl/flags/flag.h"

// Planning-space options
ABSL_FLAG(
    std::string, pp_space_definition, "",
    "Specify the file containing the space settings (order, interp, etc.).");
ABSL_FLAG(int, pp_space_order, -1,
          "Differential order of the planning space (0,1,2,...).");
ABSL_FLAG(std::string, pp_interpolation_method, "",
          "The interpolation method to use in the planning (can be: linear, "
          "cubic, quintic, svp, sap).");
ABSL_FLAG(bool, pp_use_temporal_space, false,
          "Specify that the planning space is temporal (space-time).");
ABSL_FLAG(bool, pp_use_rate_limited_space, false,
          "Specify that the planning space is rate-limited (reach-time "
          "normalization).");
ABSL_FLAG(double, pp_min_travel_dist, -1.0,
          "Minimum distance to travel between 'distinct' points in the "
          "planning-space.");
ABSL_FLAG(double, pp_max_travel_dist, -1.0,
          "Maximum distance to travel between points in the planning-space "
          "(max reachable).");
ABSL_FLAG(int, pp_output_space_order, -1,
          "Differential order of the output-space (0,1,2,...).");
ABSL_FLAG(std::string, pp_output_interp_method, "",
          "The interpolation method to use in the output (can be: linear, "
          "cubic, quintic, svp, sap).");
ABSL_FLAG(bool, pp_use_rate_limited_output_space, false,
          "Specify that the output-space is rate-limited (reach-time "
          "normalization).");

namespace ReaK::pp {

planning_space_options get_planning_space_options_from_flags() {
  planning_space_options space_options;

  if (!absl::GetFlag(FLAGS_pp_space_definition).empty()) {
    try {
      (*serialization::open_iarchive(
          absl::GetFlag(FLAGS_pp_space_definition))) >>
          space_options;
    } catch (std::exception& e) {
      RK_UNUSED(e);
    }
  }

  if (absl::GetFlag(FLAGS_pp_space_order) >= 0) {
    space_options.set_space_order(absl::GetFlag(FLAGS_pp_space_order));
  }

  if (!absl::GetFlag(FLAGS_pp_interpolation_method).empty()) {
    if (absl::GetFlag(FLAGS_pp_interpolation_method) == "linear") {
      space_options.set_interp_id(0);
    } else if (absl::GetFlag(FLAGS_pp_interpolation_method) == "cubic") {
      space_options.set_interp_id(1);
    } else if (absl::GetFlag(FLAGS_pp_interpolation_method) == "quintic") {
      space_options.set_interp_id(2);
    } else if (absl::GetFlag(FLAGS_pp_interpolation_method) == "p") {
      space_options.set_interp_id(3);
    } else if (absl::GetFlag(FLAGS_pp_interpolation_method) == "sap") {
      space_options.set_interp_id(4);
    }
  }

  space_options.set_temporal_space(absl::GetFlag(FLAGS_pp_use_temporal_space));
  space_options.set_temporal_output_space(
      absl::GetFlag(FLAGS_pp_use_temporal_space));

  space_options.set_rate_limited(
      absl::GetFlag(FLAGS_pp_use_rate_limited_space));

  if (absl::GetFlag(FLAGS_pp_min_travel_dist) > 0.0) {
    space_options.min_travel = absl::GetFlag(FLAGS_pp_min_travel_dist);
  }
  if (absl::GetFlag(FLAGS_pp_max_travel_dist) > 0.0) {
    space_options.max_travel = absl::GetFlag(FLAGS_pp_max_travel_dist);
  }

  if (absl::GetFlag(FLAGS_pp_output_space_order) >= 0) {
    space_options.set_output_space_order(
        absl::GetFlag(FLAGS_pp_output_space_order));
  }

  if (!absl::GetFlag(FLAGS_pp_output_interp_method).empty()) {
    if (absl::GetFlag(FLAGS_pp_output_interp_method) == "linear") {
      space_options.set_output_interp_id(0);
    } else if (absl::GetFlag(FLAGS_pp_output_interp_method) == "cubic") {
      space_options.set_output_interp_id(1);
    } else if (absl::GetFlag(FLAGS_pp_output_interp_method) == "quintic") {
      space_options.set_output_interp_id(2);
    } else if (absl::GetFlag(FLAGS_pp_output_interp_method) == "p") {
      space_options.set_output_interp_id(3);
    } else if (absl::GetFlag(FLAGS_pp_output_interp_method) == "sap") {
      space_options.set_output_interp_id(4);
    }
  }

  space_options.set_output_rate_limited(
      absl::GetFlag(FLAGS_pp_use_rate_limited_output_space));

  return space_options;
}

}  // namespace ReaK::pp
