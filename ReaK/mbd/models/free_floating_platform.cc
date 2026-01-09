
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

#include "ReaK/mbd/models/free_floating_platform.h"

#include "ReaK/mbd/kte/free_joints.h"
#include "ReaK/mbd/kte/rigid_link.h"

#include <cmath>
#include <utility>
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/optimization/optim_exceptions.h"

namespace ReaK::kte {

//     std::shared_ptr< frame_2D<double> > m_base_frame;
//     std::shared_ptr< frame_2D<double> > m_state_frame;
//     std::shared_ptr< jacobian_2D_2D<double> > m_state_jacobian;
//     std::shared_ptr< frame_2D<double> > m_output_frame;
//     std::vector< std::shared_ptr< joint_dependent_frame_2D > > m_EEs;
//     std::vector< pose_2D<double> > m_EEposes;  // relative to output-frame
//
//     std::shared_ptr< kte_map_chain > m_chain;

free_floater_2D_kinematics::free_floater_2D_kinematics(
    const std::string& aName, std::shared_ptr<frame_2D<double>> aBaseFrame)
    : inverse_kinematics_model(aName), m_base_frame(std::move(aBaseFrame)) {
  if (!m_base_frame) {
    m_base_frame = std::make_shared<frame_2D<double>>();
  }

  m_state_frame = std::make_shared<frame_2D<double>>();
  m_output_frame = std::make_shared<frame_2D<double>>();
  m_output_frame->Parent = m_state_frame;
  m_state_jacobian = std::make_shared<jacobian_2D_2D<double>>();
  m_state_jacobian->Parent = m_state_frame;

  // create free-joint
  auto joint_1 = std::make_shared<free_joint_2D>(
      "free_floater_joint_1", m_state_frame, m_base_frame, m_output_frame,
      m_state_jacobian);

  m_chain = std::make_shared<kte_map_chain>("free_floater_kin_model");

  *m_chain << joint_1;
}

void free_floater_2D_kinematics::resyncEndEffectors() const {
  if (m_EEposes.size() < m_EEs.size()) {
    m_EEs.erase(m_EEs.begin() + m_EEposes.size(), m_EEs.end());
  }
  while (m_EEposes.size() > m_EEs.size()) {
    auto tmp = std::make_shared<frame_2D<double>>();
    tmp->Parent = m_output_frame;
    m_EEs.push_back(std::make_shared<joint_dependent_frame_2D>(tmp));
    m_EEs.back()->add_joint(m_state_frame, m_state_jacobian);
  }
}

void free_floater_2D_kinematics::doDirectMotion() {
  m_chain->doMotion();

  resyncEndEffectors();
  for (std::size_t i = 0; i < m_EEs.size(); ++i) {
    m_EEs[i]->mFrame->Parent = m_output_frame;
    m_EEs[i]->mFrame->Position = m_EEposes[i].Position;
    m_EEs[i]->mFrame->Rotation = m_EEposes[i].Rotation;
    m_EEs[i]->mFrame->Velocity = vect<double, 2>(0.0, 0.0);
    m_EEs[i]->mFrame->AngVelocity = 0.0;
    m_EEs[i]->mFrame->Acceleration = vect<double, 2>(0.0, 0.0);
    m_EEs[i]->mFrame->AngAcceleration = 0.0;
  }
}

void free_floater_2D_kinematics::doInverseMotion() {
  resyncEndEffectors();
  if (m_EEs.empty()) {
    return;
  }

  // in theory, all EE_outputs should be equal.
  std::vector<frame_2D<double>> EE_outputs(m_EEs.size());
  for (std::size_t i = 0; i < m_EEs.size(); ++i) {
    EE_outputs[i] = m_EEs[i]->mFrame->getFrameRelativeTo(m_base_frame);
    EE_outputs[i] *= ~(m_EEposes[i]);
  }

  // average out the EE_outputs:
  frame_2D<double> EE_fr = EE_outputs[0];
  for (std::size_t i = 1; i < m_EEs.size(); ++i) {
    EE_fr.Position = EE_fr.Position * double(i) + EE_outputs[i].Position;
    EE_fr.Position /= double(i + 1);

    EE_fr.Rotation *= rot_mat_2D<double>(
        (invert(EE_fr.Rotation) * EE_outputs[i].Rotation).getAngle() /
        double(i + 1));

    EE_fr.Velocity = EE_fr.Velocity * double(i) + EE_outputs[i].Velocity;
    EE_fr.Velocity /= double(i + 1);

    EE_fr.AngVelocity =
        (EE_fr.AngVelocity * double(i) + EE_outputs[i].AngVelocity) /
        double(i + 1);

    EE_fr.Acceleration =
        EE_fr.Acceleration * double(i) + EE_outputs[i].Acceleration;
    EE_fr.Acceleration /= double(i + 1);

    EE_fr.AngAcceleration =
        (EE_fr.AngAcceleration * double(i) + EE_outputs[i].AngAcceleration) /
        double(i + 1);
  }

  // now, EE_fr is the desired output-frame, relative to the base-frame. That's the state-frame.
  (*m_output_frame) = EE_fr;
  (*m_state_frame) = EE_fr;
  m_state_frame->Parent = std::weak_ptr<
      pose_2D<double>>();  // state-frame should not be relative to anything.

  m_chain->doMotion();
}

void free_floater_2D_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  getJacobianMatrixAndDerivativeImpl(Jac, nullptr);
}

void free_floater_2D_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  getJacobianMatrixAndDerivativeImpl(Jac, &JacDot);
}

void free_floater_2D_kinematics::getJacobianMatrixAndDerivativeImpl(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>* JacDot) const {
  using Matrix = mat<double, mat_structure::rectangular>;

  std::size_t m = getDependentVelocitiesCount();
  std::size_t n = getJointVelocitiesCount();
  Jac = mat<double, mat_structure::nil>(m, n);
  if (JacDot != nullptr) {
    *JacDot = mat<double, mat_structure::nil>(m, n);
  }

  std::size_t RowInd = 0;
  for (auto& m_EE : m_EEs) {
    if (m_EE->mUpStream2DJoints.find(m_state_frame) !=
        m_EE->mUpStream2DJoints.end()) {
      mat_sub_block<Matrix> subJac =
          sub(Jac)(range(RowInd, RowInd + 3), range(0, 3));
      if (JacDot != nullptr) {
        mat_sub_block<Matrix> subJacDot =
            sub(*JacDot)(range(RowInd, RowInd + 3), range(0, 3));
        m_EE->mUpStream2DJoints[m_state_frame]
            ->get_jac_relative_to(m_EE->mFrame)
            .write_to_matrices(subJac, subJacDot);
      } else {
        m_EE->mUpStream2DJoints[m_state_frame]
            ->get_jac_relative_to(m_EE->mFrame)
            .write_to_matrices(subJac);
      }
    }
    RowInd += 3;
  }
}

vect_n<double> free_floater_2D_kinematics::getJointPositions() const {
  return {m_state_frame->Position[0], m_state_frame->Position[1],
          m_state_frame->Rotation[0], m_state_frame->Rotation[1]};
}

void free_floater_2D_kinematics::setJointPositions(
    const vect_n<double>& aJointPositions) {
  m_state_frame->Position[0] = aJointPositions[0];
  m_state_frame->Position[1] = aJointPositions[1];
  m_state_frame->Rotation = rot_mat_2D<double>(
      vect<double, 2>(aJointPositions[2], aJointPositions[3]));
}

vect_n<double> free_floater_2D_kinematics::getJointVelocities() const {
  return {m_state_frame->Velocity[0], m_state_frame->Velocity[1],
          m_state_frame->AngVelocity};
}

void free_floater_2D_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_state_frame->Velocity[0] = aJointVelocities[0];
  m_state_frame->Velocity[1] = aJointVelocities[1];
  m_state_frame->AngVelocity = aJointVelocities[2];
}

vect_n<double> free_floater_2D_kinematics::getJointAccelerations() const {
  return {m_state_frame->Acceleration[0], m_state_frame->Acceleration[1],
          m_state_frame->AngAcceleration};
}

void free_floater_2D_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_state_frame->Acceleration[0] = aJointAccelerations[0];
  m_state_frame->Acceleration[1] = aJointAccelerations[1];
  m_state_frame->AngAcceleration = aJointAccelerations[2];
}

vect_n<double> free_floater_2D_kinematics::getDependentPositions() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentPositionsCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    pose_2D<double> p_gbl = m_EE->mFrame->getGlobalPose();
    result[j] = p_gbl.Position[0];
    ++j;
    result[j] = p_gbl.Position[1];
    ++j;
    result[j] = p_gbl.Rotation[0];
    ++j;
    result[j] = p_gbl.Rotation[1];
    ++j;
  }
  return result;
}

void free_floater_2D_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    pose_2D<double> p_gbl = m_EE->mFrame->getGlobalPose();
    pose_2D<double> p_inc(
        std::weak_ptr<pose_2D<double>>(),
        vect<double, 2>(aDepPositions[j], aDepPositions[j + 1]) -
            p_gbl.Position,
        invert(p_gbl.Rotation) *
            rot_mat_2D<double>(
                vect<double, 2>(aDepPositions[j + 2], aDepPositions[j + 3])));
    (*m_EE->mFrame) *= p_inc;
    j += 4;
  }
}

vect_n<double> free_floater_2D_kinematics::getDependentVelocities() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentVelocitiesCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_2D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    result[j] = p_gbl.Velocity[0];
    ++j;
    result[j] = p_gbl.Velocity[1];
    ++j;
    result[j] = p_gbl.AngVelocity;
    ++j;
  }
  return result;
}

void free_floater_2D_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_2D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    m_EE->mFrame->Velocity += p_gbl.rotateFromParent(
        vect<double, 2>(aDepVelocities[j], aDepVelocities[j + 1]) -
        p_gbl.Velocity);
    m_EE->mFrame->AngVelocity += aDepVelocities[j + 2] - p_gbl.AngVelocity;
    j += 3;
  }
}

vect_n<double> free_floater_2D_kinematics::getDependentAccelerations() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentAccelerationsCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_2D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    result[j] = p_gbl.Acceleration[0];
    ++j;
    result[j] = p_gbl.Acceleration[1];
    ++j;
    result[j] = p_gbl.AngAcceleration;
    ++j;
  }
  return result;
}

void free_floater_2D_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_2D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    m_EE->mFrame->Acceleration += p_gbl.rotateFromParent(
        vect<double, 2>(aDepAccelerations[j], aDepAccelerations[j + 1]) -
        p_gbl.Acceleration);
    m_EE->mFrame->AngAcceleration +=
        aDepAccelerations[j + 2] - p_gbl.AngAcceleration;
    j += 3;
  }
}

void free_floater_2D_kinematics::save(serialization::oarchive& A,
                                      unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::get_static_object_type()->version());
  resyncEndEffectors();
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_state_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_state_jacobian) &
      RK_SERIAL_SAVE_WITH_NAME(m_output_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_EEs) & RK_SERIAL_SAVE_WITH_NAME(m_EEposes) &
      RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void free_floater_2D_kinematics::load(serialization::iarchive& A,
                                      unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::get_static_object_type()->version());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_state_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_state_jacobian) &
      RK_SERIAL_LOAD_WITH_NAME(m_output_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_EEs) & RK_SERIAL_LOAD_WITH_NAME(m_EEposes) &
      RK_SERIAL_LOAD_WITH_NAME(m_chain);
  resyncEndEffectors();
}

//     std::shared_ptr< frame_3D<double> > m_base_frame;
//     std::shared_ptr< frame_3D<double> > m_state_frame;
//     std::shared_ptr< jacobian_2D_2D<double> > m_state_jacobian;
//     std::shared_ptr< frame_3D<double> > m_output_frame;
//     std::vector< std::shared_ptr< joint_dependent_frame_3D > > m_EEs;
//     std::vector< pose_3D<double> > m_EEposes;  // relative to output-frame
//
//     std::shared_ptr< kte_map_chain > m_chain;

free_floater_3D_kinematics::free_floater_3D_kinematics(
    const std::string& aName, std::shared_ptr<frame_3D<double>> aBaseFrame)
    : inverse_kinematics_model(aName), m_base_frame(std::move(aBaseFrame)) {
  if (!m_base_frame) {
    m_base_frame = std::make_shared<frame_3D<double>>();
  }

  m_state_frame = std::make_shared<frame_3D<double>>();
  m_output_frame = std::make_shared<frame_3D<double>>();
  m_output_frame->Parent = m_state_frame;
  m_state_jacobian = std::make_shared<jacobian_3D_3D<double>>();
  m_state_jacobian->Parent = m_state_frame;

  // create free-joint
  auto joint_1 = std::make_shared<free_joint_3D>(
      "free_floater_joint_1", m_state_frame, m_base_frame, m_output_frame,
      m_state_jacobian);

  m_chain = std::make_shared<kte_map_chain>("free_floater_kin_model");

  *m_chain << joint_1;
}

void free_floater_3D_kinematics::resyncEndEffectors() const {
  if (m_EEposes.size() < m_EEs.size()) {
    m_EEs.erase(m_EEs.begin() + m_EEposes.size(), m_EEs.end());
  }
  while (m_EEposes.size() > m_EEs.size()) {
    auto tmp = std::make_shared<frame_3D<double>>();
    tmp->Parent = m_output_frame;
    m_EEs.push_back(std::make_shared<joint_dependent_frame_3D>(tmp));
    m_EEs.back()->add_joint(m_state_frame, m_state_jacobian);
  }
}

void free_floater_3D_kinematics::doDirectMotion() {
  m_chain->doMotion();

  resyncEndEffectors();
  for (std::size_t i = 0; i < m_EEs.size(); ++i) {
    m_EEs[i]->mFrame->Parent = m_output_frame;
    m_EEs[i]->mFrame->Position = m_EEposes[i].Position;
    m_EEs[i]->mFrame->Quat = m_EEposes[i].Quat;
    m_EEs[i]->mFrame->Velocity = vect<double, 3>(0.0, 0.0, 0.0);
    m_EEs[i]->mFrame->AngVelocity = vect<double, 3>(0.0, 0.0, 0.0);
    m_EEs[i]->mFrame->Acceleration = vect<double, 3>(0.0, 0.0, 0.0);
    m_EEs[i]->mFrame->AngAcceleration = vect<double, 3>(0.0, 0.0, 0.0);
  }
}

void free_floater_3D_kinematics::doInverseMotion() {
  resyncEndEffectors();
  if (m_EEs.empty()) {
    return;
  }

  // in theory, all EE_outputs should be equal.
  std::vector<frame_3D<double>> EE_outputs(m_EEs.size());
  for (std::size_t i = 0; i < m_EEs.size(); ++i) {
    EE_outputs[i] = m_EEs[i]->mFrame->getFrameRelativeTo(m_base_frame);
    EE_outputs[i] *= ~(m_EEposes[i]);
  }

  // average out the EE_outputs:
  frame_3D<double> EE_fr = EE_outputs[0];
  std::shared_ptr<frame_3D<double>> EE_fr_ptr(&EE_fr, null_deleter());
  for (std::size_t i = 1; i < m_EEs.size(); ++i) {
    frame_3D<double> EE_diff = EE_outputs[i].getFrameRelativeTo(EE_fr_ptr);
    EE_diff.Position /= double(i + 1);
    EE_diff.Velocity /= double(i + 1);
    EE_diff.Acceleration /= double(i + 1);
    auto aa = axis_angle<double>(EE_diff.Quat);
    aa.angle() /= double(i + 1);
    EE_diff.Quat = aa;
    EE_diff.AngVelocity /= double(i + 1);
    EE_diff.AngAcceleration /= double(i + 1);
    EE_fr *= EE_diff;
  }

  // now, EE_fr is the desired output-frame, relative to the base-frame. That's the state-frame.
  (*m_output_frame) = EE_fr;
  (*m_state_frame) = EE_fr;
  m_state_frame->Parent = std::weak_ptr<
      pose_3D<double>>();  // state-frame should not be relative to anything.

  m_chain->doMotion();
}

void free_floater_3D_kinematics::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  getJacobianMatrixAndDerivativeImpl(Jac, nullptr);
}

void free_floater_3D_kinematics::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  getJacobianMatrixAndDerivativeImpl(Jac, &JacDot);
}

void free_floater_3D_kinematics::getJacobianMatrixAndDerivativeImpl(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>* JacDot) const {
  using Matrix = mat<double, mat_structure::rectangular>;

  std::size_t m = getDependentVelocitiesCount();
  std::size_t n = getJointVelocitiesCount();
  Jac = mat<double, mat_structure::nil>(m, n);
  if (JacDot != nullptr) {
    *JacDot = mat<double, mat_structure::nil>(m, n);
  }

  std::size_t RowInd = 0;
  for (auto& m_EE : m_EEs) {
    if (m_EE->mUpStream3DJoints.find(m_state_frame) !=
        m_EE->mUpStream3DJoints.end()) {
      mat_sub_block<Matrix> subJac =
          sub(Jac)(range(RowInd, RowInd + 6), range(0, 6));
      if (JacDot != nullptr) {
        mat_sub_block<Matrix> subJacDot =
            sub(*JacDot)(range(RowInd, RowInd + 6), range(0, 6));
        m_EE->mUpStream3DJoints[m_state_frame]
            ->get_jac_relative_to(m_EE->mFrame)
            .write_to_matrices(subJac, subJacDot);
      } else {
        m_EE->mUpStream3DJoints[m_state_frame]
            ->get_jac_relative_to(m_EE->mFrame)
            .write_to_matrices(subJac);
      }
    }
    RowInd += 6;
  }
}

vect_n<double> free_floater_3D_kinematics::getJointPositions() const {
  return {m_state_frame->Position[0], m_state_frame->Position[1],
          m_state_frame->Position[2], m_state_frame->Quat[0],
          m_state_frame->Quat[1],     m_state_frame->Quat[2],
          m_state_frame->Quat[3]};
}

void free_floater_3D_kinematics::setJointPositions(
    const vect_n<double>& aJointPositions) {
  m_state_frame->Position[0] = aJointPositions[0];
  m_state_frame->Position[1] = aJointPositions[1];
  m_state_frame->Position[2] = aJointPositions[2];
  m_state_frame->Quat = quaternion<double>(
      vect<double, 4>(aJointPositions[3], aJointPositions[4],
                      aJointPositions[5], aJointPositions[6]));
}

vect_n<double> free_floater_3D_kinematics::getJointVelocities() const {
  return {m_state_frame->Velocity[0],    m_state_frame->Velocity[1],
          m_state_frame->Velocity[2],    m_state_frame->AngVelocity[0],
          m_state_frame->AngVelocity[1], m_state_frame->AngVelocity[2]};
}

void free_floater_3D_kinematics::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  m_state_frame->Velocity[0] = aJointVelocities[0];
  m_state_frame->Velocity[1] = aJointVelocities[1];
  m_state_frame->Velocity[2] = aJointVelocities[2];
  m_state_frame->AngVelocity[0] = aJointVelocities[3];
  m_state_frame->AngVelocity[1] = aJointVelocities[4];
  m_state_frame->AngVelocity[2] = aJointVelocities[5];
}

vect_n<double> free_floater_3D_kinematics::getJointAccelerations() const {
  return {m_state_frame->Acceleration[0],    m_state_frame->Acceleration[1],
          m_state_frame->Acceleration[2],    m_state_frame->AngAcceleration[0],
          m_state_frame->AngAcceleration[1], m_state_frame->AngAcceleration[2]};
}

void free_floater_3D_kinematics::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  m_state_frame->Acceleration[0] = aJointAccelerations[0];
  m_state_frame->Acceleration[1] = aJointAccelerations[1];
  m_state_frame->Acceleration[2] = aJointAccelerations[2];
  m_state_frame->AngAcceleration[0] = aJointAccelerations[3];
  m_state_frame->AngAcceleration[1] = aJointAccelerations[4];
  m_state_frame->AngAcceleration[2] = aJointAccelerations[5];
}

vect_n<double> free_floater_3D_kinematics::getDependentPositions() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentPositionsCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    pose_3D<double> p_gbl = m_EE->mFrame->getGlobalPose();
    result[j] = p_gbl.Position[0];
    ++j;
    result[j] = p_gbl.Position[1];
    ++j;
    result[j] = p_gbl.Position[2];
    ++j;
    result[j] = p_gbl.Quat[0];
    ++j;
    result[j] = p_gbl.Quat[1];
    ++j;
    result[j] = p_gbl.Quat[2];
    ++j;
    result[j] = p_gbl.Quat[3];
    ++j;
  }
  return result;
}

void free_floater_3D_kinematics::setDependentPositions(
    const vect_n<double>& aDepPositions) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    pose_3D<double> p_gbl = m_EE->mFrame->getGlobalPose();
    pose_3D<double> p_inc(
        std::weak_ptr<pose_3D<double>>(),
        vect<double, 3>(aDepPositions[j], aDepPositions[j + 1],
                        aDepPositions[j + 2]) -
            p_gbl.Position,
        invert(p_gbl.Quat) * quaternion<double>(vect<double, 4>(
                                 aDepPositions[j + 3], aDepPositions[j + 4],
                                 aDepPositions[j + 5], aDepPositions[j + 6])));
    (*m_EE->mFrame) *= p_inc;
    j += 7;
  }
}

vect_n<double> free_floater_3D_kinematics::getDependentVelocities() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentVelocitiesCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_3D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    result[j] = p_gbl.Velocity[0];
    ++j;
    result[j] = p_gbl.Velocity[1];
    ++j;
    result[j] = p_gbl.Velocity[2];
    ++j;
    result[j] = p_gbl.AngVelocity[0];
    ++j;
    result[j] = p_gbl.AngVelocity[1];
    ++j;
    result[j] = p_gbl.AngVelocity[2];
    ++j;
  }
  return result;
}

void free_floater_3D_kinematics::setDependentVelocities(
    const vect_n<double>& aDepVelocities) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_3D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    m_EE->mFrame->Velocity += p_gbl.rotateFromParent(
        vect<double, 3>(aDepVelocities[j], aDepVelocities[j + 1],
                        aDepVelocities[j + 2]) -
        p_gbl.Velocity);
    m_EE->mFrame->AngVelocity +=
        vect<double, 3>(aDepVelocities[j + 3], aDepVelocities[j + 4],
                        aDepVelocities[j + 5]) -
        p_gbl.AngVelocity;
    j += 6;
  }
}

vect_n<double> free_floater_3D_kinematics::getDependentAccelerations() const {
  resyncEndEffectors();

  vect_n<double> result(getDependentAccelerationsCount());
  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_3D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    result[j] = p_gbl.Acceleration[0];
    ++j;
    result[j] = p_gbl.Acceleration[1];
    ++j;
    result[j] = p_gbl.Acceleration[2];
    ++j;
    result[j] = p_gbl.AngAcceleration[0];
    ++j;
    result[j] = p_gbl.AngAcceleration[1];
    ++j;
    result[j] = p_gbl.AngAcceleration[2];
    ++j;
  }
  return result;
}

void free_floater_3D_kinematics::setDependentAccelerations(
    const vect_n<double>& aDepAccelerations) {
  resyncEndEffectors();

  std::size_t j = 0;
  for (auto& m_EE : m_EEs) {
    frame_3D<double> p_gbl = m_EE->mFrame->getGlobalFrame();
    m_EE->mFrame->Acceleration += p_gbl.rotateFromParent(
        vect<double, 3>(aDepAccelerations[j], aDepAccelerations[j + 1],
                        aDepAccelerations[j + 2]) -
        p_gbl.Acceleration);
    m_EE->mFrame->AngAcceleration +=
        vect<double, 3>(aDepAccelerations[j + 3], aDepAccelerations[j + 4],
                        aDepAccelerations[j + 5]) -
        p_gbl.AngAcceleration;
    j += 6;
  }
}

void free_floater_3D_kinematics::save(serialization::oarchive& A,
                                      unsigned int /*unused*/) const {
  inverse_kinematics_model::save(
      A, inverse_kinematics_model::get_static_object_type()->version());
  resyncEndEffectors();
  A& RK_SERIAL_SAVE_WITH_NAME(m_base_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_state_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_state_jacobian) &
      RK_SERIAL_SAVE_WITH_NAME(m_output_frame) &
      RK_SERIAL_SAVE_WITH_NAME(m_EEs) & RK_SERIAL_SAVE_WITH_NAME(m_EEposes) &
      RK_SERIAL_SAVE_WITH_NAME(m_chain);
}

void free_floater_3D_kinematics::load(serialization::iarchive& A,
                                      unsigned int /*unused*/) {
  inverse_kinematics_model::load(
      A, inverse_kinematics_model::get_static_object_type()->version());
  A& RK_SERIAL_LOAD_WITH_NAME(m_base_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_state_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_state_jacobian) &
      RK_SERIAL_LOAD_WITH_NAME(m_output_frame) &
      RK_SERIAL_LOAD_WITH_NAME(m_EEs) & RK_SERIAL_LOAD_WITH_NAME(m_EEposes) &
      RK_SERIAL_LOAD_WITH_NAME(m_chain);
  resyncEndEffectors();
}
}  // namespace ReaK::kte
