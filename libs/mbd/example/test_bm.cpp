
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

#include <ReaK/mbd/kte/kte_map_chain.hpp>

#include <ReaK/mbd/kte/inertia.hpp>
#include <ReaK/mbd/kte/jacobian_joint_map.hpp>
#include <ReaK/mbd/kte/mass_matrix_calculator.hpp>
#include <ReaK/mbd/kte/revolute_joint.hpp>
#include <ReaK/mbd/kte/rigid_link.hpp>

#include <ReaK/core/recorders/ascii_recorder.hpp>

#include <ReaK/core/serialization/xml_archiver.hpp>

using namespace ReaK;

using namespace ReaK::kte;

using namespace serialization;

int main() {

#if 1
  auto base_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  auto joint_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  auto joint_jacobian = rtti::rk_dynamic_ptr_cast<jacobian_gen_2D<double>>(
      jacobian_gen_2D<double>::Create());
  auto end_frame =
      rtti::rk_dynamic_ptr_cast<frame_2D<double>>(frame_2D<double>::Create());
  auto joint_coord =
      rtti::rk_dynamic_ptr_cast<gen_coord<double>>(gen_coord<double>::Create());

  base_frame->Acceleration = vect<double, 2>(0, 9.81);  // add gravity

  // create revolute joint
  auto rev_joint = std::make_shared<revolute_joint_2D>(
      "joint1", joint_coord, base_frame, joint_frame, joint_jacobian);
  // create link of lenght 0.5 meters
  auto link1 = std::make_shared<rigid_link_2D>(
      "link1", joint_frame, end_frame,
      pose_2D<double>(std::weak_ptr<pose_2D<double>>(),
                      vect<double, 2>(0.5, 0.0), rot_mat_2D<double>(0.0)));

  jacobian_joint_map_2D joint_map1;
  joint_map1[joint_coord] = joint_jacobian;
  // create end mass of 1.0 kg (point mass only)
  auto mass1 = std::make_shared<inertia_2D>(
      "mass1",
      std::make_shared<joint_dependent_frame_2D>(end_frame, joint_map1), 1.0,
      0.0);

  auto mass_mat1 = std::make_shared<mass_matrix_calc>("mass_mat1");

  *mass_mat1 << mass1;
  *mass_mat1 << joint_coord;

  kte_map_chain pendulum("pendulum");
  pendulum << rev_joint << link1 << mass1;
  {
    xml_oarchive pendulum_arc("pendulum.xml");
    oarchive& arc_ref = pendulum_arc;
    arc_ref << joint_coord << pendulum << end_frame << mass_mat1;
  };
#else
  std::shared_ptr<gen_coord<double>> joint_coord;
  std::shared_ptr<frame_2D<double>> end_frame;
  std::shared_ptr<mass_matrix_calc> mass_mat1;

  kte_map_chain pendulum;

  {
    xml_iarchive pendulum_arc("pendulum.xml");
    iarchive& arc_ref = pendulum_arc;
    arc_ref >> joint_coord >> pendulum >> end_frame >> mass_mat1;
  };

  std::cout << "Pendulum model loaded.. staring simulation!" << std::endl;

#endif

  recorder::ascii_recorder output_rec("pendulum_results.ssv");
  output_rec << "time"
             << "q"
             << "qd"
             << "qdd"
             << "f" << recorder::data_recorder::end_name_row;

  double sim_time = 0.0;

  for (; sim_time < 50.0; sim_time += 0.001) {

    joint_coord->q_ddot = 0.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl = joint_coord->f;

    mat<double, mat_structure::symmetric> M;
    mass_mat1->getMassMatrix(M);

    joint_coord->q_ddot = 1.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl_in = joint_coord->f;

    // std::cout << (f_nl - f_nl_in) << " " << M(0,0) << std::endl;

    joint_coord->q_ddot = f_nl / (f_nl - f_nl_in);

    output_rec << sim_time << joint_coord->q << joint_coord->q_dot
               << joint_coord->q_ddot << f_nl
               << recorder::data_recorder::end_value_row;

    joint_coord->q += joint_coord->q_dot * 0.001;
    joint_coord->q_dot += joint_coord->q_ddot * 0.001;
  };

  output_rec << recorder::data_recorder::close;
};
