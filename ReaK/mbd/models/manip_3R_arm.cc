
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
#include "ReaK/mbd/models/manip_3R_arm.h"

#include <cmath>
#include <utility>

namespace ReaK::kte {

//     std::shared_ptr< frame_2D<double> > m_base_frame;
//     std::vector< std::shared_ptr< gen_coord<double> > > m_joints;
//     std::shared_ptr< joint_dependent_frame_2D > m_EE;
//     double m_link1_length;
//     double m_link2_length;
//     double m_link3_length;

manip_3R_2D_kinematics::manip_3R_2D_kinematics(
    const std::string& aName, std::shared_ptr<frame_2D<double>> aBaseFrame,
    double aLink1Length, double aLink2Length, double aLink3Length)
    : inverse_kinematics_model(aName),
      m_base_frame(std::move(aBaseFrame)),
      m_link1_length(aLink1Length),
      m_link2_length(aLink2Length),
      m_link3_length(aLink3Length) {
  if (!m_base_frame) {
    m_base_frame = std::make_shared<frame_2D<double>>();
  }

  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());

  // declare all the intermediate frames.
  auto joint_1_end = std::make_shared<frame_2D<double>>();
  auto joint_2_base = std::make_shared<frame_2D<double>>();
  auto joint_2_end = std::make_shared<frame_2D<double>>();
  auto joint_3_base = std::make_shared<frame_2D<double>>();
  auto joint_3_end = std::make_shared<frame_2D<double>>();
  auto arm_EE = std::make_shared<frame_2D<double>>();

  // declare all the joint jacobians.
  auto joint_1_jacobian = std::make_shared<jacobian_gen_2D<double>>();
  auto joint_2_jacobian = std::make_shared<jacobian_gen_2D<double>>();
  auto joint_3_jacobian = std::make_shared<jacobian_gen_2D<double>>();

  // create revolute joint
  auto joint_1 = std::make_shared<revolute_joint_2D>(
      "manip_3R_joint_1", m_joints[0], m_base_frame, joint_1_end,
      joint_1_jacobian);

  // create link from F to CM (note that this is very approximate!!!)
  auto link_1 = std::make_shared<rigid_link_2D>(
      "manip_3R_link_1", joint_1_end, joint_2_base,
      pose_2D<double>(std::weak_ptr<pose_2D<double>>(),
                      vect<double, 2>(m_link1_length, 0.0),
                      rot_mat_2D<double>()));

  // create revolute joint
  auto joint_2 = std::make_shared<revolute_joint_2D>(
      "manip_3R_joint_2", m_joints[1], joint_2_base, joint_2_end,
      joint_2_jacobian);

  // create link
  auto link_2 = std::make_shared<rigid_link_2D>(
      "manip_3R_link_2", joint_2_end, joint_3_base,
      pose_2D<double>(std::weak_ptr<pose_2D<double>>(),
                      vect<double, 2>(m_link2_length, 0.0),
                      rot_mat_2D<double>()));

  // create revolute joint
  auto joint_3 = std::make_shared<revolute_joint_2D>(
      "manip_3R_joint_3", m_joints[2], joint_3_base, joint_3_end,
      joint_3_jacobian);

  // create link
  auto link_3 = std::make_shared<rigid_link_2D>(
      "manip_3R_link_3", joint_3_end, arm_EE,
      pose_2D<double>(std::weak_ptr<pose_2D<double>>(),
                      vect<double, 2>(m_link3_length, 0.0),
                      rot_mat_2D<double>()));

  // create inertia
  m_EE = std::make_shared<joint_dependent_frame_2D>(arm_EE);

  m_EE->add_joint(m_joints[0], joint_1_jacobian);
  m_EE->add_joint(m_joints[1], joint_2_jacobian);
  m_EE->add_joint(m_joints[2], joint_3_jacobian);

  m_chain = std::make_shared<kte_map_chain>("manip_3R_kin_model");

  *m_chain << joint_1 << link_1 << joint_2 << link_2 << joint_3 << link_3;
}

void manip_3R_2D_kinematics::doDirectMotion() {
  m_chain->doMotion();
}

void manip_3R_2D_kinematics::doInverseMotion() {
  using std::abs;
  using std::acos;
  using std::atan2;
  using std::sin;

  vect<double, 2> w_pos =
      m_EE->mFrame->transformToGlobal(vect<double, 2>(-m_link3_length, 0.0));
  vect<double, 2> w_v = m_EE->mFrame->rotateToGlobal(vect<double, 2>(1.0, 0.0));
  frame_2D<double> EE_fr = m_EE->mFrame->getFrameRelativeTo(m_base_frame);

  if (m_base_frame) {
    w_pos = m_base_frame->transformFromGlobal(w_pos);
    w_v = m_base_frame->rotateFromGlobal(w_v);
  }

  double w_dist = norm_2(w_pos);
  if (w_dist > m_link1_length + m_link2_length) {
    throw optim::infeasible_problem(
        "Inverse kinematics problem is infeasible! Desired end-effector pose "
        "of 2D 3R "
        "manipulator is outside the reachable workspace!");
  }
  if (w_dist < 1.000001 * abs(m_link1_length - m_link2_length)) {
    throw optim::infeasible_problem(
        "Inverse kinematics problem is infeasible! Desired end-effector pose "
        "of 2D 3R "
        "manipulator is too close to the center (inner work-space limit)!");
  }

  double a2 = acos(-(m_link1_length * m_link1_length +
                     m_link2_length * m_link2_length - w_dist * w_dist) /
                   (2.0 * m_link1_length * m_link2_length));
  if (abs(a2 - m_joints[0]->q) > abs(a2 + m_joints[0]->q)) {
    a2 = -a2;
  }

  double c1_p = (m_link1_length * m_link1_length + w_dist * w_dist -
                 m_link2_length * m_link2_length) /
                (2.0 * m_link1_length * w_dist);
  double s1_p = m_link2_length * sin(a2) / w_dist;
  double c1_c = w_pos[0] / w_dist;
  double s1_c = w_pos[1] / w_dist;
  double c1 = c1_c * c1_p - s1_c * s1_p;
  double s1 = s1_c * c1_p + c1_c * s1_p;
  double a1 = atan2(s1, c1);

  auto R1 = rot_mat_2D<double>(a1);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(a1 + a2);
  vect<double, 2> w_v_local = w_v * R12;
  double a3 = atan2(w_v_local[1], w_v_local[0]);

  vect<double, 2> cl1 =
      double(1.0) % (R1 * vect<double, 2>(m_link1_length, 0.0));
  vect<double, 2> cl2 =
      double(1.0) % (R12 * vect<double, 2>(m_link2_length, 0.0));
  vect<double, 2> cl3 = double(1.0) % (R12 * rot_mat_2D<double>(a3) *
                                       vect<double, 2>(m_link3_length, 0.0));

  mat<double, mat_structure::rectangular> A(3, 3);
  A(0, 0) = cl1[0];
  A(1, 0) = cl1[1];
  A(2, 0) = 1.0;
  A(0, 1) = cl2[0];
  A(1, 1) = cl2[1];
  A(2, 1) = 1.0;
  A(0, 2) = cl3[0];
  A(1, 2) = cl3[1];
  A(2, 2) = 1.0;
  mat<double, mat_structure::rectangular> b(3, 1);
  b(0, 0) = EE_fr.Velocity[0];
  b(1, 0) = EE_fr.Velocity[1];
  b(2, 0) = EE_fr.AngVelocity;
  mat<double, mat_structure::rectangular> jt_vel(3, 1);
  linlsq_RRQR(A, jt_vel, b);  // solve for the joint velocities.

  m_joints[0]->q = a1;
  m_joints[0]->q_dot = jt_vel(0, 0);
  m_joints[0]->q_ddot = 0.0;
  m_joints[1]->q = a2;
  m_joints[1]->q_dot = jt_vel(1, 0);
  m_joints[1]->q_ddot = 0.0;
  m_joints[2]->q = a3;
  m_joints[2]->q_dot = jt_vel(2, 0);
  m_joints[2]->q_ddot = 0.0;

  m_chain->doMotion();
}

void manip_3R_2D_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  auto R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  rot_mat_2D<double> R123 =
      rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q + m_joints[2]->q);

  vect<double, 2> cl1 =
      double(1.0) % (R1 * vect<double, 2>(m_link1_length, 0.0));
  vect<double, 2> cl2 =
      double(1.0) % (R12 * vect<double, 2>(m_link2_length, 0.0));
  vect<double, 2> cl3 =
      double(1.0) % (R123 * vect<double, 2>(m_link3_length, 0.0));

  Jac.resize(std::make_pair(3, 3));
  Jac(0, 0) = cl1[0];
  Jac(1, 0) = cl1[1];
  Jac(2, 0) = 1.0;
  Jac(0, 1) = cl2[0];
  Jac(1, 1) = cl2[1];
  Jac(2, 1) = 1.0;
  Jac(0, 2) = cl3[0];
  Jac(1, 2) = cl3[1];
  Jac(2, 2) = 1.0;
}

void manip_3R_2D_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  auto R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  rot_mat_2D<double> R123 =
      rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q + m_joints[2]->q);

  vect<double, 2> l1 = R1 * vect<double, 2>(m_link1_length, 0.0);
  vect<double, 2> l2 = R12 * vect<double, 2>(m_link2_length, 0.0);
  vect<double, 2> l3 = R123 * vect<double, 2>(m_link3_length, 0.0);

  Jac.resize(std::make_pair(3, 3));
  Jac(0, 0) = -l1[1];
  Jac(1, 0) = l1[0];
  Jac(2, 0) = 1.0;
  Jac(0, 1) = -l2[1];
  Jac(1, 1) = l2[0];
  Jac(2, 1) = 1.0;
  Jac(0, 2) = -l3[1];
  Jac(1, 2) = l3[0];
  Jac(2, 2) = 1.0;

  JacDot.resize(std::make_pair(3, 3));
  JacDot(0, 0) = -m_joints[0]->q_dot * l1[0];
  JacDot(1, 0) = -m_joints[0]->q_dot * l1[1];
  JacDot(2, 0) = 0.0;
  JacDot(0, 1) = -m_joints[1]->q_dot * l2[0];
  JacDot(1, 1) = -m_joints[1]->q_dot * l2[1];
  JacDot(2, 1) = 0.0;
  JacDot(0, 2) = -m_joints[2]->q_dot * l3[0];
  JacDot(1, 2) = -m_joints[2]->q_dot * l3[1];
  JacDot(2, 2) = 0.0;
}

vect_n<double> manip_3R_2D_kinematics::getJointPositions() const {
  return {m_joints[0]->q, m_joints[1]->q, m_joints[2]->q};
}

void manip_3R_2D_kinematics::setJointPositions(
    const vect_n<double>& aJointPositions) {
  m_joints[0]->q = aJointPositions[0];
  m_joints[1]->q = aJointPositions[1];
  m_joints[2]->q = aJointPositions[2];
}

vect_n<double> manip_3R_2D_kinematics::getJointVelocities() const {
  return {m_joints[0]->q_dot, m_joints[1]->q_dot, m_joints[2]->q_dot};
}

void manip_3R_2D_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_joints[0]->q_dot = aJointVelocities[0];
  m_joints[1]->q_dot = aJointVelocities[1];
  m_joints[2]->q_dot = aJointVelocities[2];
}

vect_n<double> manip_3R_2D_kinematics::getJointAccelerations() const {
  return {m_joints[0]->q_ddot, m_joints[1]->q_ddot, m_joints[2]->q_ddot};
}

void manip_3R_2D_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_joints[0]->q_ddot = aJointAccelerations[0];
  m_joints[1]->q_ddot = aJointAccelerations[1];
  m_joints[2]->q_ddot = aJointAccelerations[2];
}

vect_n<double> manip_3R_2D_kinematics::getDependentPositions() const {
  return {m_EE->mFrame->Position[0], m_EE->mFrame->Position[1],
          m_EE->mFrame->Rotation[0], m_EE->mFrame->Rotation[1]};
}

void manip_3R_2D_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  m_EE->mFrame->Position[0] = aDepPositions[0];
  m_EE->mFrame->Position[1] = aDepPositions[1];
  m_EE->mFrame->Rotation =
      rot_mat_2D<double>(vect<double, 2>(aDepPositions[2], aDepPositions[3]));
}

vect_n<double> manip_3R_2D_kinematics::getDependentVelocities() const {
  return {m_EE->mFrame->Velocity[0], m_EE->mFrame->Velocity[1],
          m_EE->mFrame->AngVelocity};
}

void manip_3R_2D_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  m_EE->mFrame->Velocity[0] = aDepVelocities[0];
  m_EE->mFrame->Velocity[1] = aDepVelocities[1];
  m_EE->mFrame->AngVelocity = aDepVelocities[2];
}

vect_n<double> manip_3R_2D_kinematics::getDependentAccelerations() const {
  return {m_EE->mFrame->Acceleration[0], m_EE->mFrame->Acceleration[1],
          m_EE->mFrame->AngAcceleration};
}

void manip_3R_2D_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  m_EE->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_EE->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_EE->mFrame->AngAcceleration = aDepAccelerations[2];
}

void manip_3R_2D_kinematics::save(serialization::oarchive& A,
                                  unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_joints) & RK_SERIAL_SAVE_WITH_NAME(m_EE) &
      RK_SERIAL_SAVE_WITH_NAME(m_link1_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_link2_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_link3_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void manip_3R_2D_kinematics::load(serialization::iarchive& A,
                                  unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_joints) & RK_SERIAL_LOAD_WITH_NAME(m_EE) &
      RK_SERIAL_LOAD_WITH_NAME(m_link1_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_link2_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_link3_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_chain);
}

//     std::shared_ptr< frame_3D<double> > m_base_frame;
//     std::vector< std::shared_ptr< gen_coord<double> > > m_joints;
//     std::shared_ptr< joint_dependent_frame_3D > m_EE;
//     double m_link1_length;
//     double m_link1_dz;
//     double m_link2_length;
//     double m_link2_dz;
//     double m_link3_length;
//     double m_link3_dz;

manip_3R_3D_kinematics::manip_3R_3D_kinematics(
    const std::string& aName, std::shared_ptr<frame_3D<double>> aBaseFrame,
    double aLink1Length, double aLink1DZ, double aLink2Length, double aLink2DZ,
    double aLink3Length, double aLink3DZ)
    : inverse_kinematics_model(aName),
      m_base_frame(std::move(aBaseFrame)),
      m_link1_length(aLink1Length),
      m_link1_dz(aLink1DZ),
      m_link2_length(aLink2Length),
      m_link2_dz(aLink2DZ),
      m_link3_length(aLink3Length),
      m_link3_dz(aLink3DZ) {
  if (!m_base_frame) {
    m_base_frame = std::make_shared<frame_3D<double>>();
  }

  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());
  m_joints.push_back(std::make_shared<gen_coord<double>>());

  // declare all the intermediate frames.
  auto joint_1_end = std::make_shared<frame_3D<double>>();
  auto joint_2_base = std::make_shared<frame_3D<double>>();
  auto joint_2_end = std::make_shared<frame_3D<double>>();
  auto joint_3_base = std::make_shared<frame_3D<double>>();
  auto joint_3_end = std::make_shared<frame_3D<double>>();
  auto arm_EE = std::make_shared<frame_3D<double>>();

  // declare all the joint jacobians.
  auto joint_1_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_2_jacobian = std::make_shared<jacobian_gen_3D<double>>();
  auto joint_3_jacobian = std::make_shared<jacobian_gen_3D<double>>();

  // create revolute joint
  auto joint_1 = std::make_shared<kte::revolute_joint_3D>(
      "manip_3R_joint_1", m_joints[0], vect<double, 3>(0.0, 0.0, 1.0),
      m_base_frame, joint_1_end, joint_1_jacobian);

  // create link from F to CM (note that this is very approximate!!!)
  auto link_1 = std::make_shared<kte::rigid_link_3D>(
      "manip_3R_link_1", joint_1_end, joint_2_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(m_link1_length, 0.0, m_link1_dz),
                      quaternion<double>()));

  // create revolute joint
  auto joint_2 = std::make_shared<kte::revolute_joint_3D>(
      "manip_3R_joint_2", m_joints[1], vect<double, 3>(0.0, 0.0, 1.0),
      joint_2_base, joint_2_end, joint_2_jacobian);

  // create link
  auto link_2 = std::make_shared<kte::rigid_link_3D>(
      "manip_3R_link_2", joint_2_end, joint_3_base,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(m_link2_length, 0.0, m_link2_dz),
                      quaternion<double>()));

  // create revolute joint
  auto joint_3 = std::make_shared<kte::revolute_joint_3D>(
      "manip_3R_joint_3", m_joints[2], vect<double, 3>(0.0, 0.0, 1.0),
      joint_3_base, joint_3_end, joint_3_jacobian);

  // create link
  auto link_3 = std::make_shared<kte::rigid_link_3D>(
      "manip_3R_link_3", joint_3_end, arm_EE,
      pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                      vect<double, 3>(m_link3_length, 0.0, m_link3_dz),
                      quaternion<double>()));

  // create inertia
  m_EE = std::make_shared<joint_dependent_frame_3D>(arm_EE);

  m_EE->add_joint(m_joints[0], joint_1_jacobian);
  m_EE->add_joint(m_joints[1], joint_2_jacobian);
  m_EE->add_joint(m_joints[2], joint_3_jacobian);

  m_chain = std::make_shared<kte_map_chain>("manip_3R_kin_model");

  *m_chain << joint_1 << link_1 << joint_2 << link_2 << joint_3 << link_3;
}

void manip_3R_3D_kinematics::doDirectMotion() {
  m_chain->doMotion();
}

void manip_3R_3D_kinematics::doInverseMotion() {
  using std::abs;
  using std::acos;
  using std::atan2;
  using std::sin;
  using std::sqrt;

  vect<double, 3> w_pos = m_EE->mFrame->transformToGlobal(
      vect<double, 3>(-m_link3_length, 0.0, -m_link3_dz));
  vect<double, 3> w_v =
      m_EE->mFrame->rotateToGlobal(vect<double, 3>(1.0, 0.0, 0.0));
  frame_3D<double> EE_fr = m_EE->mFrame->getFrameRelativeTo(m_base_frame);

  if (m_base_frame) {
    w_pos = m_base_frame->transformFromGlobal(w_pos);
    w_v = m_base_frame->rotateFromGlobal(w_v);
  }

  double w_dist = sqrt(w_pos[0] * w_pos[0] + w_pos[1] * w_pos[1]);
  if (w_dist > m_link1_length + m_link2_length) {
    throw optim::infeasible_problem(
        "Inverse kinematics problem is infeasible! Desired end-effector pose "
        "of 3D 3R "
        "manipulator is outside the reachable workspace!");
  }
  if (w_dist < 1.000001 * abs(m_link1_length - m_link2_length)) {
    throw optim::infeasible_problem(
        "Inverse kinematics problem is infeasible! Desired end-effector pose "
        "of 3D 3R "
        "manipulator is too close to the center (inner work-space limit)!");
  }

  double a2 = acos(-(m_link1_length * m_link1_length +
                     m_link2_length * m_link2_length - w_dist * w_dist) /
                   (2.0 * m_link1_length * m_link2_length));
  if (abs(a2 - m_joints[0]->q) > abs(a2 + m_joints[0]->q)) {
    a2 = -a2;
  }

  double c1_p = (m_link1_length * m_link1_length + w_dist * w_dist -
                 m_link2_length * m_link2_length) /
                (2.0 * m_link1_length * w_dist);
  double s1_p = m_link2_length * sin(a2) / w_dist;
  double c1_c = w_pos[0] / w_dist;
  double s1_c = w_pos[1] / w_dist;
  double c1 = c1_c * c1_p - s1_c * s1_p;
  double s1 = s1_c * c1_p + c1_c * s1_p;
  double a1 = atan2(s1, c1);

  auto R1 = rot_mat_2D<double>(a1);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(a1 + a2);
  vect<double, 2> w_v_local = vect<double, 2>(w_v[0], w_v[1]) * R12;
  double a3 = atan2(w_v_local[1], w_v_local[0]);
  rot_mat_2D<double> R123 = rot_mat_2D<double>(a1 + a2 + a3);

  vect<double, 2> cl1 =
      double(1.0) % (R1 * vect<double, 2>(m_link1_length, 0.0));
  vect<double, 2> cl2 =
      double(1.0) % (R12 * vect<double, 2>(m_link2_length, 0.0));
  vect<double, 2> cl3 =
      double(1.0) % (R123 * vect<double, 2>(m_link3_length, 0.0));

  mat<double, mat_structure::rectangular> A(3, 3);
  A(0, 0) = cl1[0];
  A(1, 0) = cl1[1];
  A(2, 0) = 1.0;
  A(0, 1) = cl2[0];
  A(1, 1) = cl2[1];
  A(2, 1) = 1.0;
  A(0, 2) = cl3[0];
  A(1, 2) = cl3[1];
  A(2, 2) = 1.0;
  mat<double, mat_structure::rectangular> b(3, 1);
  b(0, 0) = EE_fr.Velocity[0];
  b(1, 0) = EE_fr.Velocity[1];
  b(2, 0) = EE_fr.AngVelocity[2];
  mat<double, mat_structure::rectangular> jt_vel(3, 1);
  linlsq_RRQR(A, jt_vel, b);  // solve for the joint velocities.

  m_joints[0]->q = a1;
  m_joints[0]->q_dot = jt_vel(0, 0);
  m_joints[0]->q_ddot = 0.0;
  m_joints[1]->q = a2;
  m_joints[1]->q_dot = jt_vel(1, 0);
  m_joints[1]->q_ddot = 0.0;
  m_joints[2]->q = a3;
  m_joints[2]->q_dot = jt_vel(2, 0);
  m_joints[2]->q_ddot = 0.0;

  m_chain->doMotion();
}

void manip_3R_3D_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  auto R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  rot_mat_2D<double> R123 =
      rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q + m_joints[2]->q);

  vect<double, 2> cl1 =
      double(1.0) % (R1 * vect<double, 2>(m_link1_length, 0.0));
  vect<double, 2> cl2 =
      double(1.0) % (R12 * vect<double, 2>(m_link2_length, 0.0));
  vect<double, 2> cl3 =
      double(1.0) % (R123 * vect<double, 2>(m_link3_length, 0.0));

  Jac.resize(std::make_pair(6, 3));
  Jac(0, 0) = cl1[0];
  Jac(1, 0) = cl1[1];
  Jac(2, 0) = 0.0;
  Jac(3, 0) = 0.0;
  Jac(4, 0) = 0.0;
  Jac(5, 0) = 1.0;
  Jac(0, 1) = cl2[0];
  Jac(1, 1) = cl2[1];
  Jac(2, 1) = 0.0;
  Jac(3, 1) = 0.0;
  Jac(4, 1) = 0.0;
  Jac(5, 1) = 1.0;
  Jac(0, 2) = cl3[0];
  Jac(1, 2) = cl3[1];
  Jac(2, 2) = 0.0;
  Jac(3, 2) = 0.0;
  Jac(4, 2) = 0.0;
  Jac(5, 2) = 1.0;
}

void manip_3R_3D_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  auto R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  rot_mat_2D<double> R123 =
      rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q + m_joints[2]->q);

  vect<double, 2> l1 = R1 * vect<double, 2>(m_link1_length, 0.0);
  vect<double, 2> l2 = R12 * vect<double, 2>(m_link2_length, 0.0);
  vect<double, 2> l3 = R123 * vect<double, 2>(m_link3_length, 0.0);

  Jac.resize(std::make_pair(6, 3));
  Jac(0, 0) = -l1[1];
  Jac(1, 0) = l1[0];
  Jac(2, 0) = 0.0;
  Jac(3, 0) = 0.0;
  Jac(4, 0) = 0.0;
  Jac(5, 0) = 1.0;
  Jac(0, 1) = -l2[1];
  Jac(1, 1) = l2[0];
  Jac(2, 1) = 0.0;
  Jac(3, 1) = 0.0;
  Jac(4, 1) = 0.0;
  Jac(5, 1) = 1.0;
  Jac(0, 2) = -l3[1];
  Jac(1, 2) = l3[0];
  Jac(2, 2) = 0.0;
  Jac(3, 2) = 0.0;
  Jac(4, 2) = 0.0;
  Jac(5, 2) = 1.0;

  JacDot.resize(std::make_pair(6, 3));
  JacDot(0, 0) = -m_joints[0]->q_dot * l1[0];
  JacDot(1, 0) = -m_joints[0]->q_dot * l1[1];
  JacDot(2, 0) = 0.0;
  JacDot(3, 0) = 0.0;
  JacDot(4, 0) = 0.0;
  JacDot(5, 0) = 0.0;
  JacDot(0, 1) = -m_joints[1]->q_dot * l2[0];
  JacDot(1, 1) = -m_joints[1]->q_dot * l2[1];
  JacDot(2, 1) = 0.0;
  JacDot(3, 1) = 0.0;
  JacDot(4, 1) = 0.0;
  JacDot(5, 1) = 0.0;
  JacDot(0, 2) = -m_joints[2]->q_dot * l3[0];
  JacDot(1, 2) = -m_joints[2]->q_dot * l3[1];
  JacDot(2, 2) = 0.0;
  JacDot(3, 2) = 0.0;
  JacDot(4, 2) = 0.0;
  JacDot(5, 2) = 0.0;
}

vect_n<double> manip_3R_3D_kinematics::getJointPositions() const {
  return {m_joints[0]->q, m_joints[1]->q, m_joints[2]->q};
}

void manip_3R_3D_kinematics::setJointPositions(
    const vect_n<double>& aJointPositions) {
  m_joints[0]->q = aJointPositions[0];
  m_joints[1]->q = aJointPositions[1];
  m_joints[2]->q = aJointPositions[2];
}

vect_n<double> manip_3R_3D_kinematics::getJointVelocities() const {
  return {m_joints[0]->q_dot, m_joints[1]->q_dot, m_joints[2]->q_dot};
}

void manip_3R_3D_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_joints[0]->q_dot = aJointVelocities[0];
  m_joints[1]->q_dot = aJointVelocities[1];
  m_joints[2]->q_dot = aJointVelocities[2];
}

vect_n<double> manip_3R_3D_kinematics::getJointAccelerations() const {
  return {m_joints[0]->q_ddot, m_joints[1]->q_ddot, m_joints[2]->q_ddot};
}

void manip_3R_3D_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_joints[0]->q_ddot = aJointAccelerations[0];
  m_joints[1]->q_ddot = aJointAccelerations[1];
  m_joints[2]->q_ddot = aJointAccelerations[2];
}

vect_n<double> manip_3R_3D_kinematics::getDependentPositions() const {
  return {m_EE->mFrame->Position[0], m_EE->mFrame->Position[1],
          m_EE->mFrame->Position[2], m_EE->mFrame->Quat[0],
          m_EE->mFrame->Quat[1],     m_EE->mFrame->Quat[2],
          m_EE->mFrame->Quat[3]};
}

void manip_3R_3D_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  m_EE->mFrame->Position[0] = aDepPositions[0];
  m_EE->mFrame->Position[1] = aDepPositions[1];
  m_EE->mFrame->Position[2] = aDepPositions[2];
  m_EE->mFrame->Quat = quaternion<double>(vect<double, 4>(
      aDepPositions[3], aDepPositions[4], aDepPositions[5], aDepPositions[6]));
}

vect_n<double> manip_3R_3D_kinematics::getDependentVelocities() const {
  return {m_EE->mFrame->Velocity[0],    m_EE->mFrame->Velocity[1],
          m_EE->mFrame->Velocity[2],    m_EE->mFrame->AngVelocity[0],
          m_EE->mFrame->AngVelocity[1], m_EE->mFrame->AngVelocity[2]};
}

void manip_3R_3D_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  m_EE->mFrame->Velocity[0] = aDepVelocities[0];
  m_EE->mFrame->Velocity[1] = aDepVelocities[1];
  m_EE->mFrame->Velocity[2] = aDepVelocities[2];
  m_EE->mFrame->AngVelocity[0] = aDepVelocities[3];
  m_EE->mFrame->AngVelocity[1] = aDepVelocities[4];
  m_EE->mFrame->AngVelocity[2] = aDepVelocities[5];
}

vect_n<double> manip_3R_3D_kinematics::getDependentAccelerations() const {
  return {m_EE->mFrame->Acceleration[0],    m_EE->mFrame->Acceleration[1],
          m_EE->mFrame->Acceleration[2],    m_EE->mFrame->AngAcceleration[0],
          m_EE->mFrame->AngAcceleration[1], m_EE->mFrame->AngAcceleration[2]};
}

void manip_3R_3D_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  m_EE->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_EE->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_EE->mFrame->Acceleration[2] = aDepAccelerations[2];
  m_EE->mFrame->AngAcceleration[0] = aDepAccelerations[3];
  m_EE->mFrame->AngAcceleration[1] = aDepAccelerations[4];
  m_EE->mFrame->AngAcceleration[2] = aDepAccelerations[5];
}

void manip_3R_3D_kinematics::save(serialization::oarchive& A,
                                  unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_joints) & RK_SERIAL_SAVE_WITH_NAME(m_EE) &
      RK_SERIAL_SAVE_WITH_NAME(m_link1_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_link1_dz) &
      RK_SERIAL_SAVE_WITH_NAME(m_link2_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_link2_dz) &
      RK_SERIAL_SAVE_WITH_NAME(m_link3_length) &
      RK_SERIAL_SAVE_WITH_NAME(m_link3_dz) & RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void manip_3R_3D_kinematics::load(serialization::iarchive& A,
                                  unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_joints) & RK_SERIAL_LOAD_WITH_NAME(m_EE) &
      RK_SERIAL_LOAD_WITH_NAME(m_link1_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_link1_dz) &
      RK_SERIAL_LOAD_WITH_NAME(m_link2_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_link2_dz) &
      RK_SERIAL_LOAD_WITH_NAME(m_link3_length) &
      RK_SERIAL_LOAD_WITH_NAME(m_link3_dz) & RK_SERIAL_LOAD_WITH_NAME(m_chain);
}
}  // namespace ReaK::kte
