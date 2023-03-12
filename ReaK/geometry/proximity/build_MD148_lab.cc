/**
 * \file build_MD148_lab.cpp
 *
 * This application constructs a geometric model of the relevant objects in the MD148 laboratory
 * environment for the experimental tests with the CRS robot and the airship.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2013
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

#include "ReaK/geometry/shapes/box.h"
#include "ReaK/geometry/shapes/capped_cylinder.h"
#include "ReaK/geometry/shapes/colored_model.h"
#include "ReaK/geometry/shapes/coord_arrows_3D.h"
#include "ReaK/geometry/shapes/plane.h"

#include "ReaK/geometry/proximity/proxy_query_model.h"

#include "ReaK/core/serialization/archiver_factory.h"

#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

namespace fs = std::filesystem;

// I/O options
ABSL_FLAG(std::string, output_path, "models",
          "Specify the output path (default is 'models').");
ABSL_FLAG(std::string, output_name, "MD148_lab",
          "Specify the output base-name (default is 'MD148_lab').");
ABSL_FLAG(std::string, format, "xml",
          "Specify the format that should be outputted (default is 'xml', but "
          "can also be 'bin' or 'protobuf').");

int main(int argc, char** argv) {

  absl::ParseCommandLine(argc, argv);

  std::string output_base_name = absl::GetFlag(FLAGS_output_name);

  std::string output_path_name = absl::GetFlag(FLAGS_output_path);
  while (output_path_name[output_path_name.length() - 1] == '/') {
    output_path_name.erase(output_path_name.length() - 1, 1);
  }

  fs::create_directory(output_path_name.c_str());

  std::string output_extension = ".rkx";
  if (absl::GetFlag(FLAGS_format) == "bin") {
    output_extension = ".rkb";
  } else if (absl::GetFlag(FLAGS_format) == "protobuf") {
    output_extension = ".pb";
  }

  using namespace ReaK;
  using namespace geom;

  auto MD148_basic_lab =
      std::make_shared<colored_model_3D>("MD148_basic_lab_render");
  auto MD148_lab_proxy =
      std::make_shared<proxy_query_model_3D>("MD148_basic_lab_proxy");

  auto lab_floor = std::make_shared<plane>(
      "MD148_floor", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(-0.8, -1.0, 0.0), quaternion<double>()),
      vect<double, 2>(4.0, 6.0));

  auto lab_n_wall = std::make_shared<plane>(
      "MD148_north_wall", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(1.2, -1.0, 1.5),
          axis_angle<double>(M_PI * 0.5, vect<double, 3>(0.0, -1.0, 0.0))
              .getQuaternion()),
      vect<double, 2>(3.0, 6.0));

  auto lab_w_wall = std::make_shared<plane>(
      "MD148_west_wall", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(-0.8, 2.0, 1.5),
          axis_angle<double>(M_PI * 0.5, vect<double, 3>(1.0, 0.0, 0.0))
              .getQuaternion()),
      vect<double, 2>(4.0, 3.0));

  auto lab_robot_track = std::make_shared<box>(
      "MD148_robot_track", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(0.0, -1.71, 0.15), quaternion<double>()),
      vect<double, 3>(0.4, 3.42, 0.3));

  auto lab_operator_wall = std::make_shared<plane>(
      "MD148_operator_wall", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(-0.8, -3.5, 1.5),
          axis_angle<double>(M_PI * 0.5, vect<double, 3>(-1.0, 0.0, 0.0))
              .getQuaternion()),
      vect<double, 2>(4.0, 3.0));

  auto lab_robot_track_left = std::make_shared<capped_cylinder>(
      "MD148_robot_track_left", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(0.1, -1.71, 0.15),
          axis_angle<double>(M_PI * 0.5, vect<double, 3>(1.0, 0.0, 0.0))
              .getQuaternion()),
      3.42, 0.18);

  auto lab_robot_track_right = std::make_shared<capped_cylinder>(
      "MD148_robot_track_right", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(-0.1, -1.71, 0.15),
          axis_angle<double>(M_PI * 0.5, vect<double, 3>(1.0, 0.0, 0.0))
              .getQuaternion()),
      3.42, 0.18);

  auto lab_global_arrows = std::make_shared<coord_arrows_3D>(
      "global_frame_arrows", std::shared_ptr<pose_3D<double>>(),
      pose_3D<double>(), 1.0);

  (*MD148_basic_lab)
      .addElement(color(0, 0, 0), lab_global_arrows)
      .addElement(color(0.3, 0.3, 0.3), lab_floor)
      .addElement(color(0.35, 0.35, 0.35), lab_robot_track);

  (*MD148_lab_proxy)
      .addShape(lab_floor)
      .addShape(lab_n_wall)
      .addShape(lab_w_wall)
      .addShape(lab_operator_wall)
      .addShape(lab_robot_track_left)
      .addShape(lab_robot_track_right);

  (*serialization::open_oarchive(output_path_name + "/" + output_base_name +
                                 ".geom" + output_extension))
      << MD148_basic_lab << MD148_lab_proxy;

  return 0;
}
