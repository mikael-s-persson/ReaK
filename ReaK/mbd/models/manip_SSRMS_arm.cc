
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

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/optimization/optim_exceptions.h"

#include "ReaK/mbd/kte/revolute_joint.h"
#include "ReaK/mbd/kte/rigid_link.h"
#include "ReaK/mbd/models/manip_SSRMS_arm.h"

#include <cmath>
#include <utility>

namespace ReaK::kte {

//     std::shared_ptr< frame_3D<double> > m_base_frame;
//     std::vector< std::shared_ptr< gen_coord<double> > > m_joints;
//     std::shared_ptr< joint_dependent_frame_3D > m_EE;
//     vect_n<double> link_lengths;
//     vect_n<double> joint_offsets;
//     vect_n<double> joint_lower_bounds;
//     vect_n<double> joint_upper_bounds;

manip_SSRMS_kinematics::manip_SSRMS_kinematics(
    const std::string& aName, std::shared_ptr<frame_3D<double>> aBaseFrame,
    const vect_n<double>& aLinkLengths, const vect_n<double>& aJointOffsets,
    const vect_n<double>& aJointLowerBounds,
    const vect_n<double>& aJointUpperBounds)
    : inverse_kinematics_model(aName),
      m_base_frame(std::move(aBaseFrame)),
      link_lengths(aLinkLengths),
      joint_offsets(aJointOffsets),
      joint_lower_bounds(aJointLowerBounds),
      joint_upper_bounds(aJointUpperBounds) {
  if (!m_base_frame) {
    m_base_frame = std::make_shared<frame_3D<double>>();
  }

  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());

  // declare all the intermediate frames.
  auto joint_1_base = std::make_shared<frame_3D<double>>();
  auto joint_1_end = std::make_shared<frame_3D<double>>();
  auto joint_2_base = std::make_shared<frame_3D<double>>();
  auto joint_2_end = std::make_shared<frame_3D<double>>();
  auto joint_3_base = std::make_shared<frame_3D<double>>();
  auto joint_3_end = std::make_shared<frame_3D<double>>();
  auto joint_4_base = std::make_shared<frame_3D<double>>();
  auto joint_4_end = std::make_shared<frame_3D<double>>();
  auto joint_5_base = std::make_shared<frame_3D<double>>();
  auto joint_5_end = std::make_shared<frame_3D<double>>();
  auto joint_6_base = std::make_shared<frame_3D<double>>();
  auto joint_6_end = std::make_shared<frame_3D<double>>();
  auto joint_7_base = std::make_shared<frame_3D<double>>();
  auto joint_7_end = std::make_shared<frame_3D<double>>();
  auto arm_EE = std::make_shared<frame_3D<double>>();

  // declare all the joint jacobians.
  auto joint_1_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_2_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_3_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_4_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_5_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_6_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_7_jacobian = std::make_shared<jacobian_gen_3D<double>>();

  // create link
  auto link_0 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_0", m_base_frame, joint_1_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(0.0, 0.0, joint_offsets[0]),
                      quaternion<double>()));

  // create revolute joint
  auto joint_1 = std::make_shared<revolute_joint_3D>(
      "manip_SSRMS_joint_1", m_joints[0], vect<double, 3>(0.0, 0.0, 1.0),
      joint_1_base, joint_1_end, joint_1_jacobian);

  // create link
  auto link_1 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_1", joint_1_end, joint_2_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[0], -joint_offsets[1], 0.0),
                      quaternion<double>()));

  // create revolute joint
  auto joint_2 = std::make_shared<revolute_joint_3D>(
      "manip_SSRMS_joint_2", m_joints[1], vect<double, 3>(0.0, -1.0, 0.0),
      joint_2_base, joint_2_end, joint_2_jacobian);

  // create link
  auto link_2 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_2", joint_2_end, joint_3_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[1], 0.0, joint_offsets[2]),
                      quaternion<double>()));

  // create revolute joint
  auto joint_3 = std::make_shared<revolute_joint_3D>(
      "manip_SSRMS_joint_3", m_joints[2], vect<double, 3>(0.0, 0.0, 1.0),
      joint_3_base, joint_3_end, joint_3_jacobian);

  // create link
  auto link_3 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_3", joint_3_end, joint_4_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[2], 0.0, joint_offsets[3]),
                      quaternion<double>()));

  // create revolute joint
  auto joint_4 = std::make_shared<revolute_joint_3D>(
      "manip_SSRMS_joint_4", m_joints[3], vect<double, 3>(0.0, 0.0, 1.0),
      joint_4_base, joint_4_end, joint_4_jacobian);

  // create link
  auto link_4 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_4", joint_4_end, joint_5_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[3], 0.0, joint_offsets[4]),
                      quaternion<double>()));

  // create revolute joint
  auto joint_5 = std::make_shared<revolute_joint_3D>(
      "manip_SSRMS_joint_5", m_joints[4], vect<double, 3>(0.0, 0.0, 1.0),
      joint_5_base, joint_5_end, joint_5_jacobian);

  // create link
  auto link_5 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_5", joint_5_end, joint_6_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[4], 0.0, 0.0),
                      quaternion<double>()));

  // create revolute joint
  auto joint_6 = std::make_shared<kte::revolute_joint_3D>(
      "manip_SSRMS_joint_6", m_joints[5], vect<double, 3>(0.0, -1.0, 0.0),
      joint_6_base, joint_6_end, joint_6_jacobian);

  // create link
  auto link_6 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_6", joint_6_end, joint_7_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(link_lengths[5], -joint_offsets[5], 0.0),
                      quaternion<double>::xrot(M_PI).getQuaternion()));

  // create revolute joint
  auto joint_7 = std::make_shared<kte::revolute_joint_3D>(
      "manip_SSRMS_joint_7", m_joints[5], vect<double, 3>(0.0, 0.0, 1.0),
      joint_7_base, joint_7_end, joint_7_jacobian);

  // create link
  auto link_7 = std::make_shared<rigid_link_3D>(
      "manip_SSRMS_link_7", joint_7_end, arm_EE,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(0.0, 0.0, joint_offsets[6]),
                      quaternion<double>()));

  // create inertia
  m_EE = std::make_shared<joint_dependent_frame_3D>(arm_EE);

  m_EE->add_joint(m_joints[0], joint_1_jacobian);
  m_EE->add_joint(m_joints[1], joint_2_jacobian);
  m_EE->add_joint(m_joints[2], joint_3_jacobian);
  m_EE->add_joint(m_joints[3], joint_4_jacobian);
  m_EE->add_joint(m_joints[4], joint_5_jacobian);
  m_EE->add_joint(m_joints[5], joint_6_jacobian);
  m_EE->add_joint(m_joints[6], joint_7_jacobian);

  m_chain = std::make_shared<kte_map_chain>("manip_SSRMS_kin_model");

  *m_chain << link_0 << joint_1 << link_1 << joint_2 << link_2 << joint_3
           << link_3 << joint_4 << link_4 << joint_5 << link_5 << joint_6
           << link_6 << joint_7 << link_7;
}

void manip_SSRMS_kinematics::doDirectMotion() {
  m_chain->doMotion();
}

static double clamp_to_pi_range(double a) {
  return (a > M_PI ? a - 2.0 * M_PI : (a < -M_PI ? a + 2.0 * M_PI : a));
}

void manip_SSRMS_kinematics::doInverseMotion() {
  using std::abs;
  using std::acos;
  using std::asin;
  using std::atan2;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  const double wrist_to_wrist_dist = link_lengths[2] + link_lengths[3];
  const double plane_offset =
      joint_offsets[2] + joint_offsets[3] + joint_offsets[4];

  const double pos_epsilon = 1e-4 * wrist_to_wrist_dist;
  const double extend_epsilon = .01 * wrist_to_wrist_dist;

  frame_3D<double> EE_fr = m_EE->mFrame->getFrameRelativeTo(m_base_frame);

  vect<double, 3> wrist_to_wrist = EE_fr.Position;
  wrist_to_wrist[2] -= joint_offsets[0];
  vect<double, 3> EE_x_axis = EE_fr.Quat * vect_i;
  vect<double, 3> EE_y_axis = EE_fr.Quat * vect_j;
  vect<double, 3> EE_z_axis = EE_fr.Quat * vect_k;
  wrist_to_wrist -= joint_offsets[6] * EE_z_axis;

  vect<double, 3> w1_y_axis(0.0, -1.0, 0.0);
  vect<double, 3> w2_y_axis(0.0, -1.0, 0.0);
  vect<double, 3> w1_z_axis(0.0, 0.0, 1.0);

  //   std::cout << "wrist to wrist = " << wrist_to_wrist << std::endl;

  if (1.0 - abs(EE_z_axis[2]) > pos_epsilon) {
    double wrist_plane_to_base_offset = wrist_to_wrist * EE_z_axis;
    if (abs(wrist_plane_to_base_offset) < plane_offset) {

      //       std::cout << "Warning: EE plane not supported!" << std::endl;

      // Fix the y1 axis:
      if (wrist_plane_to_base_offset > 0.0) {
        w1_y_axis = unit(vect<double, 3>(EE_z_axis[0], EE_z_axis[1], 0.0));
      } else {
        w1_y_axis = unit(vect<double, 3>(-EE_z_axis[0], -EE_z_axis[1], 0.0));
      }

      // Solve for the y2 axis and cross vector.
      vect<double, 3> wrist_to_wrist_planar =
          wrist_to_wrist - abs(wrist_plane_to_base_offset) * w1_y_axis;
      // now the planar w2w vector should be in the w2 plane.
      vect<double, 3> y1_cross_w2w = unit(w1_y_axis % wrist_to_wrist_planar);
      double w2w_planar_mag = norm_2(wrist_to_wrist_planar);
      double w2_w2w_ang = asin(plane_offset / w2w_planar_mag);
      w2_y_axis = wrist_to_wrist_planar -
                  plane_offset * (cos(w2_w2w_ang) * y1_cross_w2w +
                                  (sin(w2_w2w_ang) / w2w_planar_mag) *
                                      wrist_to_wrist_planar);
      w2_y_axis = unit(w2_y_axis);
      w1_z_axis = unit(w1_y_axis % w2_y_axis);
    } else if (abs(wrist_to_wrist[2]) < plane_offset) {

      //       std::cout << "Warning: horizontal plane not supported!" << std::endl;

      // Fix the y2 axis:
      if (wrist_to_wrist[2] > 0.0) {
        w2_y_axis =
            unit(vect<double, 3>(0.0, 0.0, 1.0) - EE_z_axis[2] * EE_z_axis);
      } else {
        w2_y_axis =
            unit(vect<double, 3>(0.0, 0.0, -1.0) + EE_z_axis[2] * EE_z_axis);
      }

      // Solve for the y1 axis and cross vector.
      vect<double, 3> wrist_to_wrist_planar =
          wrist_to_wrist - abs(wrist_to_wrist[2]) * w2_y_axis;
      // now the planar w2w vector should be in the w2 plane.
      vect<double, 3> y2_cross_w2w = unit(wrist_to_wrist_planar % w2_y_axis);
      double w2w_planar_mag = norm_2(wrist_to_wrist_planar);
      double w1_w2w_ang = asin(plane_offset / w2w_planar_mag);
      w1_y_axis = wrist_to_wrist_planar -
                  plane_offset * (cos(w1_w2w_ang) * y2_cross_w2w -
                                  (sin(w1_w2w_ang) / w2w_planar_mag) *
                                      wrist_to_wrist_planar);
      w1_y_axis = unit(w1_y_axis);
      w1_z_axis = unit(w1_y_axis % w2_y_axis);
    } else {

      //       std::cout << "Warning: general solution being computed!" << std::endl;

      vect<double, 3> y1y2_cross(0.0, 0.0, 0.0);
      vect<double, 3> wrist_to_wrist_planar = wrist_to_wrist;
      double denom = EE_z_axis[0] * EE_z_axis[0] + EE_z_axis[1] * EE_z_axis[1];
      for (std::size_t i = 0; i < 20; ++i) {
        double d = wrist_to_wrist_planar * EE_z_axis;
        w1_y_axis[0] = EE_z_axis[0] * d / denom;
        w1_y_axis[1] = EE_z_axis[1] * d / denom;
        w1_y_axis[2] = 0.0;
        w2_y_axis = wrist_to_wrist_planar - w1_y_axis;
        //         std::cout << "w1 cross w2w = " << (w1_y_axis % wrist_to_wrist_planar) << std::endl;
        //         std::cout << "w2w cross w2 = " << (wrist_to_wrist_planar % w2_y_axis) << std::endl;
        vect<double, 3> y1y2_new = unit(w1_y_axis % wrist_to_wrist_planar +
                                        wrist_to_wrist_planar % w2_y_axis);
        if (norm_2(y1y2_new - y1y2_cross) * plane_offset > pos_epsilon) {
          y1y2_cross = unit((y1y2_new + y1y2_cross) * 0.5);
        } else {
          break;
        }
        wrist_to_wrist_planar = wrist_to_wrist - plane_offset * y1y2_cross;
      }
      // normalize the y vectors:
      w1_y_axis *= (1.0 / norm_2(w1_y_axis));
      w2_y_axis *= (1.0 / norm_2(w2_y_axis));
      w1_z_axis = unit(w1_y_axis % w2_y_axis);
    }
  } else {

    //     std::cout << "Warning: aligned EE not supported!" << std::endl;

    // project the w2w vector onto the base plane:
    vect<double, 3> wrist_to_wrist_planar(wrist_to_wrist[0], wrist_to_wrist[1],
                                          0.0);
    double w2w_planar_mag = norm_2(wrist_to_wrist_planar);
    if (w2w_planar_mag < plane_offset) {
      throw optim::infeasible_problem(
          "Inverse kinematics problem is infeasible! End-effector pose is "
          "out-of-reach! "
          "Desired wrist position is too close to the center axis of the SSRMS "
          "arm.");
    }
    // now the planar w2w vector should be in the w2 plane.
    vect<double, 3> y_cross_w2w(-wrist_to_wrist[1] / w2w_planar_mag,
                                wrist_to_wrist[0] / w2w_planar_mag, 0.0);
    double w_w2w_ang = asin(plane_offset / w2w_planar_mag);

    w1_z_axis = cos(w_w2w_ang) * y_cross_w2w +
                (sin(w_w2w_ang) / w2w_planar_mag) * wrist_to_wrist_planar;
    w1_y_axis = wrist_to_wrist_planar - plane_offset * w1_z_axis;
    w1_y_axis = unit(w1_y_axis);
    w2_y_axis = w1_y_axis;
  }

  vect<double, 3> w1_x_axis = w1_z_axis % w1_y_axis;
  //   std::cout << "w1 x-axis = " << w1_x_axis << std::endl;
  //   std::cout << "w1 y-axis = " << w1_y_axis << std::endl;
  //   std::cout << "w1 z-axis = " << w1_z_axis << std::endl;
  wrist_to_wrist -= plane_offset * w1_z_axis;

  double c1 = -w1_y_axis[1];
  double s1 = w1_y_axis[0];
  m_joints[0]->q = atan2(s1, c1);

  double c2 = -w1_x_axis[0] * w1_y_axis[1] + w1_x_axis[1] * w1_y_axis[0];
  double s2 = w1_x_axis[2];
  m_joints[1]->q = atan2(s2, c2);

  vect<double, 3> w2_z_axis = w1_z_axis;
  vect<double, 3> w2_x_axis = w2_z_axis % w2_y_axis;
  //   std::cout << "w2 x-axis = " << w2_x_axis << std::endl;
  //   std::cout << "w2 y-axis = " << w2_y_axis << std::endl;
  //   std::cout << "w2 z-axis = " << w2_z_axis << std::endl;

  double c7 = -w2_y_axis * EE_y_axis;
  double s7 = w2_y_axis * EE_x_axis;
  m_joints[6]->q = atan2(s7, c7);

  double c6 = -w2_x_axis * (w2_y_axis % EE_z_axis);
  double s6 = -w2_x_axis * EE_z_axis;
  m_joints[5]->q = atan2(s6, c6);

  vect<double, 3> shoulder_to_shoulder = wrist_to_wrist -
                                         link_lengths[1] * w1_x_axis -
                                         link_lengths[4] * w2_x_axis;

  shoulder_to_shoulder -= joint_offsets[1] * w1_y_axis;
  shoulder_to_shoulder -= joint_offsets[5] * w2_y_axis;
  double a3_offset = atan2(-shoulder_to_shoulder * w1_y_axis,
                           shoulder_to_shoulder * w1_x_axis);
  double a5_offset =
      atan2(shoulder_to_shoulder * w2_y_axis, shoulder_to_shoulder * w2_x_axis);

  //   std::cout << "shouder-to-shoulder = " << shoulder_to_shoulder << std::endl;

  double s2s_dist_sqr = shoulder_to_shoulder * shoulder_to_shoulder;
  if (s2s_dist_sqr > (wrist_to_wrist_dist - extend_epsilon) *
                         (wrist_to_wrist_dist - extend_epsilon)) {
    throw optim::infeasible_problem(
        "Inverse kinematics problem is infeasible! End-effector pose is "
        "out-of-reach! "
        "Desired wrist position is outside the spherical workspace envelope "
        "(fully-extended arm).");
  }

  double s2s_dist = sqrt(s2s_dist_sqr);
  double l2_sqr = link_lengths[2] * link_lengths[2];
  double l3_sqr = link_lengths[3] * link_lengths[3];
  double c3 =
      (s2s_dist_sqr + l2_sqr - l3_sqr) / (2.0 * link_lengths[2] * s2s_dist);
  double c4 = (s2s_dist_sqr - l2_sqr - l3_sqr) /
              (2.0 * link_lengths[2] * link_lengths[3]);
  double c5 =
      (s2s_dist_sqr - l2_sqr + l3_sqr) / (2.0 * s2s_dist * link_lengths[3]);

  m_joints[2]->q = clamp_to_pi_range(acos(c3) + a3_offset);
  m_joints[3]->q = -acos(c4);
  m_joints[4]->q = clamp_to_pi_range(acos(c5) + a5_offset);

  if ((abs(clamp_to_pi_range(-acos(c3) + a3_offset)) < abs(m_joints[2]->q)) ||
      (abs(clamp_to_pi_range(-acos(c5) + a5_offset)) < abs(m_joints[4]->q))) {
    m_joints[2]->q = clamp_to_pi_range(-acos(c3) + a3_offset);
    m_joints[3]->q = acos(c4);
    m_joints[4]->q = clamp_to_pi_range(-acos(c5) + a5_offset);
  }

  // then, use the solution to compute the jacobian matrix:
  mat<double, mat_structure::rectangular> jac(6, 7);
  getJacobianMatrix(jac);

  // finally, use the jacobian to find a solution for the joint velocities:
  mat<double, mat_structure::rectangular> x(7, 1);
  mat<double, mat_structure::rectangular> b(6, 1);
  b(0, 0) = EE_fr.Velocity[0];
  b(1, 0) = EE_fr.Velocity[1];
  b(2, 0) = EE_fr.Velocity[2];
  b(3, 0) = EE_fr.AngVelocity[0];
  b(4, 0) = EE_fr.AngVelocity[1];
  b(5, 0) = EE_fr.AngVelocity[2];
  minnorm_RRQR(jac, x, b, pos_epsilon);

  m_joints[0]->q_dot = x(0, 0);
  m_joints[1]->q_dot = x(1, 0);
  m_joints[2]->q_dot = x(2, 0);
  m_joints[3]->q_dot = x(3, 0);
  m_joints[4]->q_dot = x(4, 0);
  m_joints[5]->q_dot = x(5, 0);
  m_joints[6]->q_dot = x(6, 0);
  // acceleration is irrelevant (not part of start variables).
  m_joints[0]->q_ddot = 0.0;
  m_joints[1]->q_ddot = 0.0;
  m_joints[2]->q_ddot = 0.0;
  m_joints[3]->q_ddot = 0.0;
  m_joints[4]->q_ddot = 0.0;
  m_joints[5]->q_ddot = 0.0;
  m_joints[6]->q_ddot = 0.0;

  m_chain->doMotion();
}

void manip_SSRMS_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  /* calculate individual rotations */
  quaternion<double>::zrot q1(m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::zrot q3(m_joints[2]->q);
  quaternion<double>::zrot q4(m_joints[3]->q);
  quaternion<double>::zrot q5(m_joints[4]->q);
  quaternion<double>::yrot q6(-m_joints[5]->q);
  quaternion<double>::xrot q67(M_PI);
  quaternion<double>::zrot q7(m_joints[6]->q);

  vect<double, 3> a0(0.0, 0.0, joint_offsets[0]);
  quaternion<double> q_accum = q1.getQuaternion();
  vect<double, 3> e1(0.0, 0.0, 1.0);
  vect<double, 3> a1 =
      q_accum * vect<double, 3>(link_lengths[0], -joint_offsets[1], 0.0);
  vect<double, 3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double, 3> a2 =
      q_accum * vect<double, 3>(link_lengths[1], 0.0, joint_offsets[2]);
  vect<double, 3> e3 = q_accum * vect_k;
  q_accum *= q3;
  vect<double, 3> a3 =
      q_accum * vect<double, 3>(link_lengths[2], 0.0, joint_offsets[3]);
  vect<double, 3> e4 = q_accum * vect_k;
  q_accum *= q4;
  vect<double, 3> a4 =
      q_accum * vect<double, 3>(link_lengths[3], 0.0, joint_offsets[4]);
  vect<double, 3> e5 = q_accum * vect_k;
  q_accum *= q5;
  vect<double, 3> a5 = q_accum * vect<double, 3>(link_lengths[4], 0.0, 0.0);
  vect<double, 3> e6 = -(q_accum * vect_j);
  q_accum *= q6;
  vect<double, 3> a6 =
      q_accum * vect<double, 3>(link_lengths[5], -joint_offsets[5], 0.0);
  q_accum *= q67;
  vect<double, 3> e7 = q_accum * vect_k;
  q_accum *= q7;
  vect<double, 3> a7 = q_accum * vect<double, 3>(0.0, 0.0, joint_offsets[6]);

  vect<double, 3> a67 = a6 + a7;
  vect<double, 3> a57 = a5 + a67;
  vect<double, 3> a47 = a4 + a57;
  vect<double, 3> a37 = a3 + a47;
  vect<double, 3> a27 = a2 + a37;
  vect<double, 3> a17 = a1 + a27;

  Jac.resize(std::make_pair(6, 7));
  vect<double, 3> v1 = vect_k % a17;
  Jac(0, 0) = v1[0];
  Jac(1, 0) = v1[1];
  Jac(2, 0) = v1[2];
  Jac(3, 0) = 0.0;
  Jac(4, 0) = 0.0;
  Jac(5, 0) = 1.0;
  vect<double, 3> v2 = e2 % a27;
  Jac(0, 1) = v2[0];
  Jac(1, 1) = v2[1];
  Jac(2, 1) = v2[2];
  Jac(3, 1) = e2[0];
  Jac(4, 1) = e2[1];
  Jac(5, 1) = e2[2];
  vect<double, 3> v3 = e3 % a37;
  Jac(0, 2) = v3[0];
  Jac(1, 2) = v3[1];
  Jac(2, 2) = v3[2];
  Jac(3, 2) = e3[0];
  Jac(4, 2) = e3[1];
  Jac(5, 2) = e3[2];
  vect<double, 3> v4 = e4 % a47;
  Jac(0, 3) = v4[0];
  Jac(1, 3) = v4[1];
  Jac(2, 3) = v4[2];
  Jac(3, 3) = e4[0];
  Jac(4, 3) = e4[1];
  Jac(5, 3) = e4[2];
  vect<double, 3> v5 = e5 % a57;
  Jac(0, 4) = v5[0];
  Jac(1, 4) = v5[1];
  Jac(2, 4) = v5[2];
  Jac(3, 4) = e5[0];
  Jac(4, 4) = e5[1];
  Jac(5, 4) = e5[2];
  vect<double, 3> v6 = e6 % a67;
  Jac(0, 5) = v6[0];
  Jac(1, 5) = v6[1];
  Jac(2, 5) = v6[2];
  Jac(3, 5) = e6[0];
  Jac(4, 5) = e6[1];
  Jac(5, 5) = e6[2];
  Jac(0, 6) = 0.0;
  Jac(1, 6) = 0.0;
  Jac(2, 6) = 0.0;
  Jac(3, 6) = e7[0];
  Jac(4, 6) = e7[1];
  Jac(5, 6) = e7[2];
}

void manip_SSRMS_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  /* calculate individual rotations */
  vect<double, 3> ey(0.0, -1.0, 0.0);
  vect<double, 3> ez(0.0, 0.0, 1.0);

  quaternion<double>::zrot q1(m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::zrot q3(m_joints[2]->q);
  quaternion<double>::zrot q4(m_joints[3]->q);
  quaternion<double>::zrot q5(m_joints[4]->q);
  quaternion<double>::yrot q6(-m_joints[5]->q);
  quaternion<double>::xrot q67(M_PI);
  quaternion<double>::zrot q7(m_joints[6]->q);

  vect<double, 3> a1_p(link_lengths[0], -joint_offsets[1], 0.0);
  vect<double, 3> a2_p(link_lengths[1], 0.0, joint_offsets[2]);
  vect<double, 3> a3_p(link_lengths[2], 0.0, joint_offsets[3]);
  vect<double, 3> a4_p(link_lengths[3], 0.0, joint_offsets[4]);
  vect<double, 3> a5_p(link_lengths[4], 0.0, 0.0);
  vect<double, 3> a6_p(link_lengths[5], -joint_offsets[5], 0.0);
  vect<double, 3> a7_p(0.0, 0.0, joint_offsets[6]);

  quaternion<double> q_accum = q1.getQuaternion();
  vect<double, 3> a1 = q_accum * a1_p;
  vect<double, 3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double, 3> a2 = q_accum * a2_p;
  vect<double, 3> e3 = q_accum * vect_k;
  q_accum *= q3;
  vect<double, 3> a3 = q_accum * a3_p;
  vect<double, 3> e4 = q_accum * vect_k;
  q_accum *= q4;
  vect<double, 3> a4 = q_accum * a4_p;
  vect<double, 3> e5 = q_accum * vect_k;
  q_accum *= q5;
  vect<double, 3> a5 = q_accum * a5_p;
  vect<double, 3> e6 = -(q_accum * vect_j);
  q_accum *= q6;
  vect<double, 3> a6 = q_accum * a6_p;
  q_accum *= q67;
  vect<double, 3> e7 = q_accum * vect_k;
  q_accum *= q7;
  vect<double, 3> a7 = q_accum * a7_p;

  vect<double, 3> a67 = a6 + a7;
  vect<double, 3> a57 = a5 + a67;
  vect<double, 3> a47 = a4 + a57;
  vect<double, 3> a37 = a3 + a47;
  vect<double, 3> a27 = a2 + a37;
  vect<double, 3> a17 = a1 + a27;

  Jac.resize(std::make_pair(6, 7));
  vect<double, 3> v1 = vect_k % a17;
  Jac(0, 0) = v1[0];
  Jac(1, 0) = v1[1];
  Jac(2, 0) = v1[2];
  Jac(3, 0) = 0.0;
  Jac(4, 0) = 0.0;
  Jac(5, 0) = 1.0;
  vect<double, 3> v2 = e2 % a27;
  Jac(0, 1) = v2[0];
  Jac(1, 1) = v2[1];
  Jac(2, 1) = v2[2];
  Jac(3, 1) = e2[0];
  Jac(4, 1) = e2[1];
  Jac(5, 1) = e2[2];
  vect<double, 3> v3 = e3 % a37;
  Jac(0, 2) = v3[0];
  Jac(1, 2) = v3[1];
  Jac(2, 2) = v3[2];
  Jac(3, 2) = e3[0];
  Jac(4, 2) = e3[1];
  Jac(5, 2) = e3[2];
  vect<double, 3> v4 = e4 % a47;
  Jac(0, 3) = v4[0];
  Jac(1, 3) = v4[1];
  Jac(2, 3) = v4[2];
  Jac(3, 3) = e4[0];
  Jac(4, 3) = e4[1];
  Jac(5, 3) = e4[2];
  vect<double, 3> v5 = e5 % a57;
  Jac(0, 4) = v5[0];
  Jac(1, 4) = v5[1];
  Jac(2, 4) = v5[2];
  Jac(3, 4) = e5[0];
  Jac(4, 4) = e5[1];
  Jac(5, 4) = e5[2];
  vect<double, 3> v6 = e6 % a67;
  Jac(0, 5) = v6[0];
  Jac(1, 5) = v6[1];
  Jac(2, 5) = v6[2];
  Jac(3, 5) = e6[0];
  Jac(4, 5) = e6[1];
  Jac(5, 5) = e6[2];
  Jac(0, 6) = 0.0;
  Jac(1, 6) = 0.0;
  Jac(2, 6) = 0.0;
  Jac(3, 6) = e7[0];
  Jac(4, 6) = e7[1];
  Jac(5, 6) = e7[2];

  // vect<double,3> q1_dot = m_joints[0]->q_dot * ez;
  // vect<double,3> q2_dot = m_joints[1]->q_dot * ey;
  // vect<double,3> q3_dot = m_joints[2]->q_dot * ez;
  // vect<double,3> q4_dot = m_joints[3]->q_dot * ez;
  // vect<double,3> q5_dot = m_joints[4]->q_dot * ez;
  // vect<double,3> q6_dot = m_joints[5]->q_dot * ey;
  // vect<double,3> q7_dot = m_joints[6]->q_dot * ez;
  vect<double, 3> a1_dot = q1 * (m_joints[0]->q_dot * ez % a1_p);

  vect<double, 3> a2_dot = q2 * (m_joints[1]->q_dot * ey % a2_p);
  a2_p = q2 * a2_p;
  a2_dot = q1 * (m_joints[0]->q_dot * ez % a2_p + a2_dot);

  vect<double, 3> a3_dot = q3 * (m_joints[2]->q_dot * ez % a3_p);
  a3_p = q3 * a3_p;
  a3_dot = q2 * (m_joints[1]->q_dot * ey % a3_p + a3_dot);
  a3_p = q2 * a3_p;
  a3_dot = q1 * (m_joints[0]->q_dot * ez % a3_p + a3_dot);

  vect<double, 3> a4_dot = q4 * (m_joints[3]->q_dot * ez % a4_p);
  a4_p = q4 * a4_p;
  a4_dot = q3 * (m_joints[2]->q_dot * ez % a4_p + a4_dot);
  a4_p = q3 * a4_p;
  a4_dot = q2 * (m_joints[1]->q_dot * ey % a4_p + a4_dot);
  a4_p = q2 * a4_p;
  a4_dot = q1 * (m_joints[0]->q_dot * ez % a4_p + a4_dot);

  vect<double, 3> a5_dot = q5 * (m_joints[4]->q_dot * ez % a5_p);
  a5_p = q5 * a5_p;
  a5_dot = q4 * (m_joints[3]->q_dot * ez % a5_p + a5_dot);
  a5_p = q4 * a5_p;
  a5_dot = q3 * (m_joints[2]->q_dot * ez % a5_p + a5_dot);
  a5_p = q3 * a5_p;
  a5_dot = q2 * (m_joints[1]->q_dot * ey % a5_p + a5_dot);
  a5_p = q2 * a5_p;
  a5_dot = q1 * (m_joints[0]->q_dot * ez % a5_p + a5_dot);

  vect<double, 3> a6_dot = q6 * (m_joints[5]->q_dot * ey % a6_p);
  a6_p = q6 * a6_p;
  a6_dot = q5 * (m_joints[4]->q_dot * ez % a6_p + a6_dot);
  a6_p = q5 * a6_p;
  a6_dot = q4 * (m_joints[3]->q_dot * ez % a6_p + a6_dot);
  a6_p = q4 * a6_p;
  a6_dot = q3 * (m_joints[2]->q_dot * ez % a6_p + a6_dot);
  a6_p = q3 * a6_p;
  a6_dot = q2 * (m_joints[1]->q_dot * ey % a6_p + a6_dot);
  a6_p = q2 * a6_p;
  a6_dot = q1 * (m_joints[0]->q_dot * ez % a6_p + a6_dot);

  vect<double, 3> a7_dot = q67 * q7 * (m_joints[6]->q_dot * ez % a7_p);
  a7_p = q67 * q7 * a7_p;
  a7_dot = q6 * (m_joints[5]->q_dot * ey % a7_p + a7_dot);
  a7_p = q6 * a7_p;
  a7_dot = q5 * (m_joints[4]->q_dot * ez % a7_p + a7_dot);
  a7_p = q5 * a7_p;
  a7_dot = q4 * (m_joints[3]->q_dot * ez % a7_p + a7_dot);
  a7_p = q4 * a7_p;
  a7_dot = q3 * (m_joints[2]->q_dot * ez % a7_p + a7_dot);
  a7_p = q3 * a7_p;
  a7_dot = q2 * (m_joints[1]->q_dot * ey % a7_p + a7_dot);
  a7_p = q2 * a7_p;
  a7_dot = q1 * (m_joints[0]->q_dot * ez % a7_p + a7_dot);

  vect<double, 3> ey_p;
  vect<double, 3> ez_p;

  vect<double, 3> e2_dot = q1 * (m_joints[0]->q_dot * ez % ey);

  vect<double, 3> e3_dot = q2 * (m_joints[1]->q_dot * ey % ez);
  ez_p = q2 * ez;
  e3_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e3_dot);

  ez_p = q3 * ez;
  vect<double, 3> e4_dot = q2 * (m_joints[1]->q_dot * ey % ez_p);
  ez_p = q2 * ez_p;
  e4_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e4_dot);

  ez_p = q4 * ez;
  vect<double, 3> e5_dot = q3 * (m_joints[2]->q_dot * ez % ez_p);
  ez_p = q3 * ez_p;
  e5_dot = q2 * (m_joints[1]->q_dot * ey % ez_p + e5_dot);
  ez_p = q2 * ez_p;
  e5_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e5_dot);

  vect<double, 3> e6_dot = q5 * (m_joints[4]->q_dot * ez % ey);
  ey_p = q5 * ey;
  e6_dot = q4 * (m_joints[3]->q_dot * ez % ey_p + e6_dot);
  ey_p = q4 * ey_p;
  e6_dot = q3 * (m_joints[2]->q_dot * ez % ey_p + e6_dot);
  ey_p = q3 * ey_p;
  e6_dot = q2 * (m_joints[1]->q_dot * ey % ey_p + e6_dot);
  ey_p = q2 * ey_p;
  e6_dot = q1 * (m_joints[0]->q_dot * ez % ey_p + e6_dot);

  ez_p = q67 * ez;
  vect<double, 3> e7_dot = q6 * (m_joints[5]->q_dot * ey % ez_p);
  ez_p = q6 * ez_p;
  e7_dot = q5 * (m_joints[4]->q_dot * ez % ez_p + e7_dot);
  ez_p = q5 * ez_p;
  e7_dot = q4 * (m_joints[3]->q_dot * ez % ez_p + e7_dot);
  ez_p = q4 * ez_p;
  e7_dot = q3 * (m_joints[2]->q_dot * ez % ez_p + e7_dot);
  ez_p = q3 * ez_p;
  e7_dot = q2 * (m_joints[1]->q_dot * ey % ez_p + e7_dot);
  ez_p = q2 * ez_p;
  e7_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e7_dot);

  vect<double, 3> a67_dot = a6_dot + a7_dot;
  vect<double, 3> a57_dot = a5_dot + a67_dot;
  vect<double, 3> a47_dot = a4_dot + a57_dot;
  vect<double, 3> a37_dot = a3_dot + a47_dot;
  vect<double, 3> a27_dot = a2_dot + a37_dot;
  vect<double, 3> a17_dot = a1_dot + a27_dot;

  JacDot.resize(std::make_pair(6, 7));
  vect<double, 3> v1_dot = vect_k % a17_dot;
  JacDot(0, 0) = v1_dot[0];
  JacDot(1, 0) = v1_dot[1];
  JacDot(2, 0) = v1_dot[2];
  JacDot(3, 0) = 0.0;
  JacDot(4, 0) = 0.0;
  JacDot(5, 0) = 0.0;
  vect<double, 3> v2_dot = e2_dot % a27 + e2 % a27_dot;
  JacDot(0, 1) = v2_dot[0];
  JacDot(1, 1) = v2_dot[1];
  JacDot(2, 1) = v2_dot[2];
  JacDot(3, 1) = e2_dot[0];
  JacDot(4, 1) = e2_dot[1];
  JacDot(5, 1) = e2_dot[2];
  vect<double, 3> v3_dot = e3_dot % a37 + e3 % a37_dot;
  JacDot(0, 2) = v3_dot[0];
  JacDot(1, 2) = v3_dot[1];
  JacDot(2, 2) = v3_dot[2];
  JacDot(3, 2) = e3_dot[0];
  JacDot(4, 2) = e3_dot[1];
  JacDot(5, 2) = e3_dot[2];
  vect<double, 3> v4_dot = e4_dot % a47 + e4 % a47_dot;
  JacDot(0, 3) = v4_dot[0];
  JacDot(1, 3) = v4_dot[1];
  JacDot(2, 3) = v4_dot[2];
  JacDot(3, 3) = e4_dot[0];
  JacDot(4, 3) = e4_dot[1];
  JacDot(5, 3) = e4_dot[2];
  vect<double, 3> v5_dot = e5_dot % a57 + e5 % a57_dot;
  JacDot(0, 4) = v5_dot[0];
  JacDot(1, 4) = v5_dot[1];
  JacDot(2, 4) = v5_dot[2];
  JacDot(3, 4) = e5_dot[0];
  JacDot(4, 4) = e5_dot[1];
  JacDot(5, 4) = e5_dot[2];
  vect<double, 3> v6_dot = e6_dot % a67 + e6 % a67_dot;
  JacDot(0, 5) = v6_dot[0];
  JacDot(1, 5) = v6_dot[1];
  JacDot(2, 5) = v6_dot[2];
  JacDot(3, 5) = e6_dot[0];
  JacDot(4, 5) = e6_dot[1];
  JacDot(5, 5) = e6_dot[2];
  JacDot(0, 6) = 0.0;
  JacDot(1, 6) = 0.0;
  JacDot(2, 6) = 0.0;
  JacDot(3, 6) = e7_dot[0];
  JacDot(4, 6) = e7_dot[1];
  JacDot(5, 6) = e7_dot[2];
}

vect_n<double> manip_SSRMS_kinematics::getJointPositions() const {
  return {m_joints[0]->q, m_joints[1]->q, m_joints[2]->q, m_joints[3]->q,
          m_joints[4]->q, m_joints[5]->q, m_joints[6]->q};
}

void manip_SSRMS_kinematics::setJointPositions(
    const vect_n<double>& aJointPositions) {
  m_joints[0]->q = aJointPositions[0];
  m_joints[1]->q = aJointPositions[1];
  m_joints[2]->q = aJointPositions[2];
  m_joints[3]->q = aJointPositions[3];
  m_joints[4]->q = aJointPositions[4];
  m_joints[5]->q = aJointPositions[5];
  m_joints[6]->q = aJointPositions[6];
}

vect_n<double> manip_SSRMS_kinematics::getJointVelocities() const {
  return {m_joints[0]->q_dot, m_joints[1]->q_dot, m_joints[2]->q_dot,
          m_joints[3]->q_dot, m_joints[4]->q_dot, m_joints[5]->q_dot,
          m_joints[6]->q_dot};
}

void manip_SSRMS_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_joints[0]->q_dot = aJointVelocities[0];
  m_joints[1]->q_dot = aJointVelocities[1];
  m_joints[2]->q_dot = aJointVelocities[2];
  m_joints[3]->q_dot = aJointVelocities[3];
  m_joints[4]->q_dot = aJointVelocities[4];
  m_joints[5]->q_dot = aJointVelocities[5];
  m_joints[6]->q_dot = aJointVelocities[6];
}

vect_n<double> manip_SSRMS_kinematics::getJointAccelerations() const {
  return {m_joints[0]->q_ddot, m_joints[1]->q_ddot, m_joints[2]->q_ddot,
          m_joints[3]->q_ddot, m_joints[4]->q_ddot, m_joints[5]->q_ddot,
          m_joints[6]->q_ddot};
}

void manip_SSRMS_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_joints[0]->q_ddot = aJointAccelerations[0];
  m_joints[1]->q_ddot = aJointAccelerations[1];
  m_joints[2]->q_ddot = aJointAccelerations[2];
  m_joints[3]->q_ddot = aJointAccelerations[3];
  m_joints[4]->q_ddot = aJointAccelerations[4];
  m_joints[5]->q_ddot = aJointAccelerations[5];
  m_joints[6]->q_ddot = aJointAccelerations[6];
}

vect_n<double> manip_SSRMS_kinematics::getDependentPositions() const {
  return {m_EE->mFrame->Position[0], m_EE->mFrame->Position[1],
          m_EE->mFrame->Position[2], m_EE->mFrame->Quat[0],
          m_EE->mFrame->Quat[1],     m_EE->mFrame->Quat[2],
          m_EE->mFrame->Quat[3]};
}

void manip_SSRMS_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  m_EE->mFrame->Position[0] = aDepPositions[0];
  m_EE->mFrame->Position[1] = aDepPositions[1];
  m_EE->mFrame->Position[2] = aDepPositions[2];
  m_EE->mFrame->Quat = quaternion<double>(vect<double, 4>(
      aDepPositions[3], aDepPositions[4], aDepPositions[5], aDepPositions[6]));
}

vect_n<double> manip_SSRMS_kinematics::getDependentVelocities() const {
  return {m_EE->mFrame->Velocity[0],    m_EE->mFrame->Velocity[1],
          m_EE->mFrame->Velocity[2],    m_EE->mFrame->AngVelocity[0],
          m_EE->mFrame->AngVelocity[1], m_EE->mFrame->AngVelocity[2]};
}

void manip_SSRMS_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  m_EE->mFrame->Velocity[0] = aDepVelocities[0];
  m_EE->mFrame->Velocity[1] = aDepVelocities[1];
  m_EE->mFrame->Velocity[2] = aDepVelocities[2];
  m_EE->mFrame->AngVelocity[0] = aDepVelocities[3];
  m_EE->mFrame->AngVelocity[1] = aDepVelocities[4];
  m_EE->mFrame->AngVelocity[2] = aDepVelocities[5];
}

vect_n<double> manip_SSRMS_kinematics::getDependentAccelerations() const {
  return {m_EE->mFrame->Acceleration[0],    m_EE->mFrame->Acceleration[1],
          m_EE->mFrame->Acceleration[2],    m_EE->mFrame->AngAcceleration[0],
          m_EE->mFrame->AngAcceleration[1], m_EE->mFrame->AngAcceleration[2]};
}

void manip_SSRMS_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  m_EE->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_EE->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_EE->mFrame->Acceleration[2] = aDepAccelerations[2];
  m_EE->mFrame->AngAcceleration[0] = aDepAccelerations[3];
  m_EE->mFrame->AngAcceleration[1] = aDepAccelerations[4];
  m_EE->mFrame->AngAcceleration[2] = aDepAccelerations[5];
}

void manip_SSRMS_kinematics::save(serialization::oarchive& A,
                                  unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_joints) & RK_SERIAL_SAVE_WITH_NAME(m_EE) &
      RK_SERIAL_SAVE_WITH_NAME(link_lengths) &
      RK_SERIAL_SAVE_WITH_NAME(joint_offsets) &
      RK_SERIAL_SAVE_WITH_NAME(joint_lower_bounds) &
      RK_SERIAL_SAVE_WITH_NAME(joint_upper_bounds) &
      RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void manip_SSRMS_kinematics::load(serialization::iarchive& A,
                                  unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_joints) & RK_SERIAL_LOAD_WITH_NAME(m_EE) &
      RK_SERIAL_LOAD_WITH_NAME(link_lengths) &
      RK_SERIAL_LOAD_WITH_NAME(joint_offsets) &
      RK_SERIAL_LOAD_WITH_NAME(joint_lower_bounds) &
      RK_SERIAL_LOAD_WITH_NAME(joint_upper_bounds) &
      RK_SERIAL_LOAD_WITH_NAME(m_chain);
}
}  // namespace ReaK::kte
