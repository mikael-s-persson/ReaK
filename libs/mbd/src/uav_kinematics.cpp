
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

#include "ReaK/math/lin_alg/mat_alg.hpp"
#include "ReaK/math/lin_alg/mat_qr_decomp.hpp"
#include "ReaK/math/optimization/optim_exceptions.hpp"

#include "ReaK/mbd/kte/free_joints.hpp"
#include "ReaK/mbd/models/uav_kinematics.hpp"

#include <cmath>
#include <utility>

namespace ReaK::kte {

UAV_kinematics::UAV_kinematics(const std::string& aName,
                               std::shared_ptr<frame_3D<double>> aBaseFrame)
    : inverse_kinematics_model(aName), m_base_frame(std::move(aBaseFrame)) {
  m_motion_frame = std::make_shared<frame_3D<double>>();

  auto EE_frame = std::make_shared<frame_3D<double>>();
  auto joint_jacobian = std::make_shared<jacobian_3D_3D<double>>();

  auto joint_1 = std::make_shared<free_joint_3D>(
      "UAV_free_joint", m_motion_frame, m_base_frame, EE_frame, joint_jacobian);

  m_output_frame = std::make_shared<joint_dependent_frame_3D>(EE_frame);
  m_output_frame->add_joint(m_motion_frame, joint_jacobian);

  m_chain = std::make_shared<kte_map_chain>("UAV_quadrotor_kin_model");

  *m_chain << joint_1;
}

void UAV_kinematics::doDirectMotion() {
  m_chain->doMotion();
}

void UAV_kinematics::doInverseMotion() {
  frame_3D<double> f = m_output_frame->mFrame->getFrameRelativeTo(m_base_frame);
  m_motion_frame->Position = f.Position;
  m_motion_frame->Quat = f.Quat;
  m_motion_frame->Velocity = f.Velocity;
  m_motion_frame->AngVelocity = f.AngVelocity;
  m_motion_frame->Acceleration = f.Acceleration;
  m_motion_frame->AngAcceleration = f.AngAcceleration;
  m_chain->doMotion();
}

void UAV_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  Jac = mat<double, mat_structure::identity>(6);
}

void UAV_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {

  Jac = mat<double, mat_structure::identity>(6);
  JacDot = mat<double, mat_structure::nil>(6, 6);
}

vect_n<double> UAV_kinematics::getJointPositions() const {
  return {m_motion_frame->Position[0], m_motion_frame->Position[1],
          m_motion_frame->Position[2], m_motion_frame->Quat[0],
          m_motion_frame->Quat[1],     m_motion_frame->Quat[2],
          m_motion_frame->Quat[3]};
}

void UAV_kinematics::setJointPositions(const vect_n<double>& aJointPositions) {
  m_motion_frame->Position[0] = aJointPositions[0];
  m_motion_frame->Position[1] = aJointPositions[1];
  m_motion_frame->Position[2] = aJointPositions[2];
  m_motion_frame->Quat = quaternion<double>(
      vect<double, 4>(aJointPositions[3], aJointPositions[4],
                      aJointPositions[5], aJointPositions[6]));
}

vect_n<double> UAV_kinematics::getJointVelocities() const {
  return {m_motion_frame->Velocity[0],    m_motion_frame->Velocity[1],
          m_motion_frame->Velocity[2],    m_motion_frame->AngVelocity[0],
          m_motion_frame->AngVelocity[1], m_motion_frame->AngVelocity[2]};
}

void UAV_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_motion_frame->Velocity[0] = aJointVelocities[0];
  m_motion_frame->Velocity[1] = aJointVelocities[1];
  m_motion_frame->Velocity[2] = aJointVelocities[2];
  m_motion_frame->AngVelocity[0] = aJointVelocities[3];
  m_motion_frame->AngVelocity[1] = aJointVelocities[4];
  m_motion_frame->AngVelocity[2] = aJointVelocities[5];
}

vect_n<double> UAV_kinematics::getJointAccelerations() const {
  return {
      m_motion_frame->Acceleration[0],    m_motion_frame->Acceleration[1],
      m_motion_frame->Acceleration[2],    m_motion_frame->AngAcceleration[0],
      m_motion_frame->AngAcceleration[1], m_motion_frame->AngAcceleration[2]};
}

void UAV_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_motion_frame->Acceleration[0] = aJointAccelerations[0];
  m_motion_frame->Acceleration[1] = aJointAccelerations[1];
  m_motion_frame->Acceleration[2] = aJointAccelerations[2];
  m_motion_frame->AngAcceleration[0] = aJointAccelerations[3];
  m_motion_frame->AngAcceleration[1] = aJointAccelerations[4];
  m_motion_frame->AngAcceleration[2] = aJointAccelerations[5];
}

vect_n<double> UAV_kinematics::getDependentPositions() const {
  return {
      m_output_frame->mFrame->Position[0], m_output_frame->mFrame->Position[1],
      m_output_frame->mFrame->Position[2], m_output_frame->mFrame->Quat[0],
      m_output_frame->mFrame->Quat[1],     m_output_frame->mFrame->Quat[2],
      m_output_frame->mFrame->Quat[3]};
}

void UAV_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  m_output_frame->mFrame->Position[0] = aDepPositions[0];
  m_output_frame->mFrame->Position[1] = aDepPositions[1];
  m_output_frame->mFrame->Position[2] = aDepPositions[2];
  m_output_frame->mFrame->Quat = quaternion<double>(vect<double, 4>(
      aDepPositions[3], aDepPositions[4], aDepPositions[5], aDepPositions[6]));
}

vect_n<double> UAV_kinematics::getDependentVelocities() const {
  return {m_output_frame->mFrame->Velocity[0],
          m_output_frame->mFrame->Velocity[1],
          m_output_frame->mFrame->Velocity[2],
          m_output_frame->mFrame->AngVelocity[0],
          m_output_frame->mFrame->AngVelocity[1],
          m_output_frame->mFrame->AngVelocity[2]};
}

void UAV_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  m_output_frame->mFrame->Velocity[0] = aDepVelocities[0];
  m_output_frame->mFrame->Velocity[1] = aDepVelocities[1];
  m_output_frame->mFrame->Velocity[2] = aDepVelocities[2];
  m_output_frame->mFrame->AngVelocity[0] = aDepVelocities[3];
  m_output_frame->mFrame->AngVelocity[1] = aDepVelocities[4];
  m_output_frame->mFrame->AngVelocity[2] = aDepVelocities[5];
}

vect_n<double> UAV_kinematics::getDependentAccelerations() const {
  return {m_output_frame->mFrame->Acceleration[0],
          m_output_frame->mFrame->Acceleration[1],
          m_output_frame->mFrame->Acceleration[2],
          m_output_frame->mFrame->AngAcceleration[0],
          m_output_frame->mFrame->AngAcceleration[1],
          m_output_frame->mFrame->AngAcceleration[2]};
}

void UAV_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  m_output_frame->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_output_frame->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_output_frame->mFrame->Acceleration[2] = aDepAccelerations[2];
  m_output_frame->mFrame->AngAcceleration[0] = aDepAccelerations[3];
  m_output_frame->mFrame->AngAcceleration[1] = aDepAccelerations[4];
  m_output_frame->mFrame->AngAcceleration[2] = aDepAccelerations[5];
}

void UAV_kinematics::save(serialization::oarchive& A,
                          unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_motion_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_output_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void UAV_kinematics::load(serialization::iarchive& A, unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_motion_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_output_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_chain);
}
}  // namespace ReaK::kte
