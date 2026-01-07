
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

#include "ReaK/mbd/kte/driving_actuator.h"
#include "ReaK/mbd/kte/free_joints.h"
#include "ReaK/mbd/kte/inertia.h"
#include "ReaK/mbd/kte/kte_map_chain.h"
#include "ReaK/mbd/kte/state_measures.h"

#include "ReaK/mbd/models/manip_dynamics_model.h"

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/geometry/shapes/box.h"
#include "ReaK/geometry/shapes/colored_model.h"
#include "ReaK/geometry/shapes/coord_arrows_3D.h"
#include "ReaK/geometry/shapes/sphere.h"

#include "ReaK/core/serialization/archiver_factory.h"

#include <filesystem>
#include <numbers>
#include <memory>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

namespace fs = std::filesystem;

// I/O options
ABSL_FLAG(std::string, output_path, "models",
          "Specify the output path (default is 'models').");
ABSL_FLAG(std::string, output_name, "airship3D",
          "Specify the output base-name (default is 'airship3D').");
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
  using namespace kte;

  auto global_base = std::make_shared<frame_3D<double>>();

  //  Output from calibration program:
  global_base->Position = vect<double, 3>(0.176438, -0.356764, -0.0521845);
  global_base->Quat = quaternion<double>(
      vect<double, 4>(0.999986, 0.00301881, -0.00257751, 0.00343278));

  auto airship3D_frame = std::make_shared<frame_3D<double>>();

  auto airship3D_output_frame = std::make_shared<frame_3D<double>>();

  auto airship3D_joint_jac = std::make_shared<jacobian_3D_3D<double>>();

  auto airship3D_joint = std::make_shared<free_joint_3D>(
      "airship3D_joint", airship3D_frame, global_base, airship3D_output_frame,
      airship3D_joint_jac);

  auto airship3D_dep_frame =
      std::make_shared<joint_dependent_frame_3D>(airship3D_output_frame);
  airship3D_dep_frame->add_joint(airship3D_frame, airship3D_joint_jac);

  auto airship3D_inertia = std::make_shared<inertia_3D>(
      "airship3D_inertia", airship3D_dep_frame, 1.0,
      mat<double, mat_structure::symmetric>(
          mat<double, mat_structure::identity>(3)));

  auto airship3D_actuator = std::make_shared<driving_actuator_3D>(
      "airship3D_actuator", airship3D_frame, airship3D_joint);

  auto airship3D_position = std::make_shared<position_measure_3D>(
      "airship3D_position", airship3D_frame);

  auto airship3D_rotation = std::make_shared<rotation_measure_3D>(
      "airship3D_rotation", airship3D_frame);

  auto airship3D_model = std::make_shared<kte_map_chain>("airship3D_model");

  (*airship3D_model) << airship3D_position << airship3D_rotation
                     << airship3D_actuator << airship3D_joint
                     << airship3D_inertia;

  auto airship3D_dyn_model =
      std::make_shared<manipulator_dynamics_model>("airship3D_dyn_model");
  airship3D_dyn_model->setModel(airship3D_model);
  (*airship3D_dyn_model) << airship3D_frame;
  (*airship3D_dyn_model)
      << airship3D_inertia;  // also adds the dependent frame 'airship3D_dep_frame'

  (*airship3D_dyn_model) << airship3D_actuator;
  (*airship3D_dyn_model) << airship3D_position;
  (*airship3D_dyn_model) << airship3D_rotation;

  // New position (lower, during bottom-heavy tests)
  double gr_radius = 0.93;
  double gr_arc_length = 0.584;

  double gr_beta = gr_arc_length / gr_radius;

  auto airship3D_grasp_frame = std::make_shared<frame_3D<double>>(
      airship3D_output_frame,
      vect<double, 3>(gr_radius * std::cos(gr_beta), 0.0,
                      gr_radius * std::sin(gr_beta)),
      quaternion<double>::yrot(-0.5 * std::numbers::pi - gr_beta).getQuaternion(),
      vect<double, 3>(0.0, 0.0, 0.0), vect<double, 3>(0.0, 0.0, 0.0),
      vect<double, 3>(0.0, 0.0, 0.0), vect<double, 3>(0.0, 0.0, 0.0),
      vect<double, 3>(0.0, 0.0, 0.0), vect<double, 3>(0.0, 0.0, 0.0));
  airship3D_grasp_frame->Position +=
      airship3D_grasp_frame->Quat * (-0.3 * vect_k);

  auto airship3D_dep_grasp_frame =
      std::make_shared<joint_dependent_frame_3D>(airship3D_grasp_frame);
  airship3D_dep_grasp_frame->add_joint(airship3D_frame, airship3D_joint_jac);

  auto hull = std::make_shared<sphere>("airship3D_hull", airship3D_output_frame,
                                       pose_3D<double>(), 0.93);

  auto grapple = std::make_shared<box>(
      "airship3D_grapple", airship3D_output_frame,
      pose_3D<double>(
          std::shared_ptr<pose_3D<double>>(),
          vect<double, 3>((gr_radius + 0.05) * std::cos(gr_beta), 0.0,
                          (gr_radius + 0.05) * std::sin(gr_beta)),
          quaternion<double>::yrot(-0.5 * std::numbers::pi - gr_beta).getQuaternion()),
      vect<double, 3>(0.08, 0.003, 0.1));

  auto grapple_arrows = std::make_shared<coord_arrows_3D>(
      "airship3D_grapple_arrows", airship3D_grasp_frame, pose_3D<double>(),
      0.3);

  auto geom_mdl = std::make_shared<colored_model_3D>("airship3D_geom_render");
  (*geom_mdl)
      .addAnchor(airship3D_output_frame)
      .addElement(color(0, 0, 0), grapple_arrows)
      .addElement(color(1, 1, 1), hull)
      .addElement(color(0.2, 0.2, 0.2), grapple);

  auto proxy_mdl =
      std::make_shared<proxy_query_model_3D>("airship3D_geom_render");
  (*proxy_mdl).addShape(hull);

  auto airship3D_kin_model =
      std::make_shared<manipulator_kinematics_model>("airship3D_kin_model");
  airship3D_kin_model->setModel(airship3D_model);
  (*airship3D_kin_model) << airship3D_frame;
  (*airship3D_kin_model) << airship3D_dep_grasp_frame;

  (*serialization::open_oarchive(output_path_name + "/" + output_base_name +
                                 ".model" + output_extension))
      << global_base << airship3D_kin_model << airship3D_grasp_frame << geom_mdl
      << proxy_mdl;
}
