
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
#include "ReaK/mbd/models/chaser_target_model_data_po.h"
#include "ReaK/mbd/models/chaser_target_model_data.h"

#include <filesystem>

#include "absl/flags/flag.h"

// Chaser-target scenario models
ABSL_FLAG(std::string, chaser_target_env, "",
          "Specify the file containing the chaser-target-env models.");
ABSL_FLAG(std::string, chaser_model_file, "",
          "Specify the file containing the chaser model.");
ABSL_FLAG(std::string, target_model_file, "",
          "Specify the file containing the target model.");
ABSL_FLAG(std::vector<std::string>, environment_models,
          std::vector<std::string>{},
          "Specify the file(s) containing the environment's geometric models.");

namespace ReaK::kte {

chaser_target_data get_chaser_target_data_from_flags() {
  chaser_target_data scene_data;

  if (!absl::GetFlag(FLAGS_chaser_target_env).empty()) {
    try {
      (*serialization::open_iarchive(absl::GetFlag(FLAGS_chaser_target_env))) >>
          scene_data;
    } catch ([[maybe_unused]] std::exception& e) {}
  }

  if (!absl::GetFlag(FLAGS_chaser_model_file).empty()) {
    scene_data.load_chaser(absl::GetFlag(FLAGS_chaser_model_file));
  }

  if (!absl::GetFlag(FLAGS_target_model_file).empty()) {
    scene_data.load_target(absl::GetFlag(FLAGS_target_model_file));
  }

  if (!absl::GetFlag(FLAGS_environment_models).empty()) {
    const std::vector<std::string>& vf =
        absl::GetFlag(FLAGS_environment_models);
    for (const auto& f_name : vf) {
      scene_data.load_environment(f_name);
    }
  }

  return scene_data;
}

}  // namespace ReaK::kte
