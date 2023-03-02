/**
 * \file chaser_target_model_data_po.hpp
 *
 * This library defines functions that provide Boost.Program-Options support for the chaser-target
 * model data (chaser_target_data class).
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

#ifndef REAK_CHASER_TARGET_MODEL_DATA_PO_HPP
#define REAK_CHASER_TARGET_MODEL_DATA_PO_HPP

#include <ReaK/core/serialization/archiver_factory.hpp>
#include "chaser_target_model_data.hpp"

#include <boost/program_options.hpp>

namespace ReaK::kte {

boost::program_options::options_description get_chaser_target_data_po_desc() {
  namespace po = boost::program_options;
  po::options_description planner_alg_options("Chaser-target scenario models");
  planner_alg_options.add_options()(
      "chaser-target-env", po::value<std::string>(),
      "specify the file containing the chaser-target-env models.")

      ("chaser-model-file", po::value<std::string>(),
       "specify the file containing the chaser model.")(
          "target-model-file", po::value<std::string>(),
          "specify the file containing the target model.")(
          "environment-models",
          po::value<std::vector<std::string>>()->multitoken(),
          "specify the file(s) containing the environment's geometric models.");
  return planner_alg_options;
}

chaser_target_data get_chaser_target_data_from_po(
    boost::program_options::variables_map& vm) {
  chaser_target_data scene_data;

  if (vm.count("chaser-target-env")) {
    try {
      (*serialization::open_iarchive(
          vm["chaser-target-env"].as<std::string>())) >>
          scene_data;
    } catch (std::exception& e) {
      RK_UNUSED(e);
    }
  }

  if (vm.count("chaser-model-file")) {
    scene_data.load_chaser(vm["chaser-model-file"].as<std::string>());
  }

  if (vm.count("target-model-file")) {
    scene_data.load_target(vm["target-model-file"].as<std::string>());
  }

  if (vm.count("environment-models")) {
    const std::vector<std::string>& vf =
        vm["environment-models"].as<std::vector<std::string>>();
    for (std::vector<std::string>::const_iterator it = vf.begin();
         it != vf.end(); ++it) {
      scene_data.load_environment(*it);
    }
  }

  return scene_data;
}
}  // namespace ReaK::kte

#endif
