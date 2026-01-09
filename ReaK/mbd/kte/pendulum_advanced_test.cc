
/*
 *    Copyright 2011 Sven Mikael Persson
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

#include "ReaK/mbd/kte/kte_map_chain.h"

#include "ReaK/mbd/kte/inertia.h"
#include "ReaK/mbd/kte/revolute_joint.h"
#include "ReaK/mbd/kte/rigid_link.h"

#include "ReaK/mbd/kte/force_actuator.h"
#include "ReaK/mbd/kte/joint_friction.h"

#include "ReaK/core/recorders/ascii_recorder.h"

#include <fstream>
#include <memory>
#include "ReaK/core/serialization/bin_archiver.h"
#include "ReaK/core/serialization/protobuf_archiver.h"
#include "ReaK/core/serialization/xml_archiver.h"

using namespace ReaK;

using namespace ReaK::kte;

using namespace serialization;

void simulate_system(std::shared_ptr<gen_coord<double>> joint_coord,
                     kte_map_chain& adv_pendulum) {
  recorder::ascii_recorder output_rec("adv_pendulum_results.ssv");
  output_rec << "time"
             << "q"
             << "qd"
             << "qdd"
             << "f" << recorder::data_recorder::end_name_row;

  double sim_time = 0.0;
  double last_sim_time = -0.01;

  for (; sim_time < 5; sim_time += 0.00001) {

    joint_coord->q_ddot = 0.0;
    adv_pendulum.doMotion();
    adv_pendulum.clearForce();
    adv_pendulum.doForce();
    double f_nl = joint_coord->f;

    joint_coord->q_ddot = 1.0;
    adv_pendulum.doMotion();
    adv_pendulum.clearForce();
    adv_pendulum.doForce();
    double f_nl_1 = joint_coord->f;

    joint_coord->q_ddot = f_nl / (f_nl - f_nl_1);

    if (sim_time >= last_sim_time + 0.01) {
      last_sim_time = sim_time;
      std::cout << "\r" << sim_time;
      std::cout.flush();
      output_rec << sim_time << joint_coord->q << joint_coord->q_dot
                 << joint_coord->q_ddot << f_nl
                 << recorder::data_recorder::end_value_row;
    };

    joint_coord->q += joint_coord->q_dot * 0.00001;
    joint_coord->q_dot += joint_coord->q_ddot * 0.00001;
  };

  output_rec << recorder::data_recorder::close;
};

int main() {

#if 1
  std::shared_ptr<frame_2D<double>> base_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  std::shared_ptr<frame_2D<double>> joint_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  std::shared_ptr<frame_2D<double>> end_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  std::shared_ptr<gen_coord<double>> joint_coord =
      rtti::rk_dynamic_ptr_cast<gen_coord<double>>(gen_coord<double>::Create());

  base_frame->Acceleration = vect<double, 2>(0, 9.81);  // add gravity

  // create motor inertia
  auto motor_inertia = std::make_shared<inertia_gen>(
      "motor_inertia", std::make_shared<joint_dependent_gen_coord>(joint_coord),
      5);
  // create friction
  auto friction = std::make_shared<joint_dry_microslip_gen>(
      "friction", joint_coord, 1E-6, 2E-6, 1, 0.9);
  // create revolute joint
  auto rev_joint = std::make_shared<revolute_joint_2D>("joint1", joint_coord,
                                                       base_frame, joint_frame);
  // create actuator
  auto actuator =
      std::make_shared<force_actuator_gen>("actuator", joint_coord, rev_joint);
  // create link of lenght 0.5 meters
  auto link1 = std::make_shared<rigid_link_2D>(
      "link1", joint_frame, end_frame,
      pose_2D<double>(std::weak_ptr<pose_2D<double>>(),
                      vect<double, 2>(0.5, 0.0), rot_mat_2D<double>(0.0)));
  // create end mass of 1.0 kg (point mass only)
  auto mass1 = std::make_shared<inertia_2D>(
      "mass1", std::make_shared<joint_dependent_frame_2D>(end_frame), 1.0, 0.0);

  kte_map_chain adv_pendulum("adv_pendulum");

  adv_pendulum << motor_inertia << friction << rev_joint << link1 << mass1;

  {
    xml_oarchive adv_pendulum_arc("models/adv_pendulum.rkx");
    adv_pendulum_arc << joint_coord << adv_pendulum << end_frame;
  };

  {
    bin_oarchive adv_pendulum_arc("models/adv_pendulum.rkb");
    adv_pendulum_arc << joint_coord << adv_pendulum << end_frame;
  };

  {
    protobuf_oarchive adv_pendulum_arc("models/adv_pendulum.rkp");
    adv_pendulum_arc << joint_coord << adv_pendulum << end_frame;
  };

  {
    protobuf_schemer adv_pendulum_sch;
    adv_pendulum_sch << joint_coord << adv_pendulum << end_frame;
    std::ofstream out_file("models/adv_pendulum.proto");
    adv_pendulum_sch.print_schemes(out_file);
  };

  kte_map_chain adv_motorized_pendulum("models/adv_motorized_pendulum");

  adv_motorized_pendulum << actuator << motor_inertia << friction << rev_joint
                         << link1 << mass1;

  {
    xml_oarchive adv_motorized_pendulum_arc(
        "models/adv_motorized_pendulum.rkx");
    adv_motorized_pendulum_arc << joint_coord << adv_pendulum << end_frame;
  };

  std::cout << "Pendulum model created.. starting simulation!" << std::endl;

  simulate_system(joint_coord, adv_pendulum);

  std::cout << "Done!" << std::endl;

#else
  std::shared_ptr<gen_coord<double>> joint_coord;
  std::shared_ptr<frame_2D<double>> end_frame;

#if 1
  kte_map_chain adv_pendulum;

  {
    xml_iarchive adv_pendulum_arc("models/adv_pendulum.rkx");
    adv_pendulum_arc >> joint_coord >> adv_pendulum >> end_frame;
  };

  std::cout << "Pendulum model loaded from xml file.. starting simulation!"
            << std::endl;

  simulate_system(joint_coord, adv_pendulum);

  std::cout << "Done!" << std::endl;

  {
    bin_iarchive adv_pendulum_arc("models/adv_pendulum.rkb");
    adv_pendulum_arc >> joint_coord >> adv_pendulum >> end_frame;
  };

  std::cout << "Pendulum model loaded from binary file.. starting simulation!"
            << std::endl;

  simulate_system(joint_coord, adv_pendulum);

  std::cout << "Done!" << std::endl;

  {
    protobuf_iarchive adv_pendulum_arc("models/adv_pendulum.rkp");
    adv_pendulum_arc >> joint_coord >> adv_pendulum >> end_frame;
  };

  std::cout << "Pendulum model loaded from protobuf file.. starting simulation!"
            << std::endl;

  simulate_system(joint_coord, adv_pendulum);

  std::cout << "Done!" << std::endl;

#else
  kte_map_chain adv_motorized_pendulum;

  {
    xml_iarchive adv_motorized_pendulum_arc(
        "models/adv_motorized_pendulum.rkx");
    iarchive& arc_ref = adv_motorized_pendulum_arc;
    arc_ref >> joint_coord >> adv_motorized_pendulum >> end_frame;
  };

  std::cout << "Motorized pendulum model loaded.. starting simulation!"
            << std::endl;

  simulate_system(joint_coord, adv_motorized_pendulum);

  std::cout << "Done!" << std::endl;

#endif

#endif

  return 0;
};
