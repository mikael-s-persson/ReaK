/**
 * \file planning_space_options_po.hpp
 *
 * This library defines functions to load planning-space options from Boost.Program-Options.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PLANNING_SPACE_OPTIONS_PO_HPP
#define REAK_PLANNING_SPACE_OPTIONS_PO_HPP

#include "planning_space_options.hpp"

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <boost/program_options.hpp>

namespace ReaK::pp {

boost::program_options::options_description
get_planning_space_options_po_desc() {
  namespace po = boost::program_options;
  po::options_description planning_space_options("Planning-space options");
  planning_space_options.add_options()(
      "space-definition", po::value<std::string>(),
      "specify the file containing the space settings (order, interp, etc.).")

      ("space-order", po::value<std::size_t>(),
       "differential order of the planning space (0,1,2,...).")(
          "interpolation-method", po::value<std::string>(),
          "the interpolation method to use in the planning (can be: linear, "
          "cubic, quintic, svp, sap).")

          ("use-temporal-space",
           "specify that the planning space is temporal (space-time).")(
              "use-rate-limited-space",
              "specify that the planning space is rate-limited (reach-time "
              "normalization).")

              ("min-travel-dist", po::value<double>(),
               "minimum distance to travel between 'distinct' points in the "
               "planning-space.")("max-travel-dist", po::value<double>(),
                                  "maximum distance to travel between points "
                                  "in the planning-space (max reachable).")

                  ("output-space-order", po::value<std::size_t>(),
                   "differential order of the output-space (0,1,2,...).")(
                      "output-interp-method", po::value<std::string>(),
                      "the interpolation method to use in the output (can be: "
                      "linear, cubic, quintic, svp, sap).")(
                      "use-rate-limited-output-space",
                      "specify that the output-space is rate-limited "
                      "(reach-time normalization).");
  return planning_space_options;
}

planning_space_options get_planning_space_options_from_po(
    boost::program_options::variables_map& vm) {
  planning_space_options space_options;

  if (vm.count("space-definition")) {
    try {
      (*serialization::open_iarchive(
          vm["space-definition"].as<std::string>())) >>
          space_options;
    } catch (std::exception& e) {
      RK_UNUSED(e);
    }
  }

  if (vm.count("space-order")) {
    space_options.set_space_order(vm["space-order"].as<std::size_t>());
  }

  if (vm.count("interpolation-method")) {
    if (vm["interpolation-method"].as<std::string>() == "linear") {
      space_options.set_interp_id(0);
    } else if (vm["interpolation-method"].as<std::string>() == "cubic") {
      space_options.set_interp_id(1);
    } else if (vm["interpolation-method"].as<std::string>() == "quintic") {
      space_options.set_interp_id(2);
    } else if (vm["interpolation-method"].as<std::string>() == "p") {
      space_options.set_interp_id(3);
    } else if (vm["interpolation-method"].as<std::string>() == "sap") {
      space_options.set_interp_id(4);
    }
  }

  if (vm.count("use-temporal-space")) {
    space_options.set_temporal_space(true);
    space_options.set_temporal_output_space(true);
  }
  if (vm.count("use-rate-limited-space")) {
    space_options.set_rate_limited(true);
  }

  if (vm.count("min-travel-dist")) {
    space_options.min_travel = vm["min-travel-dist"].as<double>();
  }
  if (vm.count("max-travel-dist")) {
    space_options.max_travel = vm["max-travel-dist"].as<double>();
  }

  if (vm.count("output-space-order")) {
    space_options.set_output_space_order(
        vm["output-space-order"].as<std::size_t>());
  }

  if (vm.count("output-interp-method")) {
    if (vm["output-interp-method"].as<std::string>() == "linear") {
      space_options.set_output_interp_id(0);
    } else if (vm["output-interp-method"].as<std::string>() == "cubic") {
      space_options.set_output_interp_id(1);
    } else if (vm["output-interp-method"].as<std::string>() == "quintic") {
      space_options.set_output_interp_id(2);
    } else if (vm["output-interp-method"].as<std::string>() == "p") {
      space_options.set_output_interp_id(3);
    } else if (vm["output-interp-method"].as<std::string>() == "sap") {
      space_options.set_output_interp_id(4);
    }
  }

  if (vm.count("use-rate-limited-output-space")) {
    space_options.set_output_rate_limited(true);
  }

  return space_options;
}

}  // namespace ReaK::pp

#endif
