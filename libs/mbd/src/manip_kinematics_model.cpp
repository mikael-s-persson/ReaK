
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

#include "ReaK/mbd/models/manip_kinematics_model.hpp"

#include "ReaK/mbd/models/manip_kinematics_helper.hpp"

namespace ReaK::kte {

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<gen_coord<double>>& aCoord) {
  if (aCoord) {
    mCoords.push_back(aCoord);
  }
  return *this;
}

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<frame_2D<double>>& aFrame2D) {
  if (aFrame2D) {
    mFrames2D.push_back(aFrame2D);
  }
  return *this;
}

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<frame_3D<double>>& aFrame3D) {
  if (aFrame3D) {
    mFrames3D.push_back(aFrame3D);
  }
  return *this;
}

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<joint_dependent_gen_coord>& aDependentGenCoord) {
  if (aDependentGenCoord) {
    mDependentGenCoords.push_back(aDependentGenCoord);
  }
  return *this;
}

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<joint_dependent_frame_2D>& aDependent2DFrame) {
  if (aDependent2DFrame) {
    mDependent2DFrames.push_back(aDependent2DFrame);
  }
  return *this;
}

manipulator_kinematics_model& manipulator_kinematics_model::operator<<(
    const std::shared_ptr<joint_dependent_frame_3D>& aDependent3DFrame) {
  if (aDependent3DFrame) {
    mDependent3DFrames.push_back(aDependent3DFrame);
  }
  return *this;
}

void manipulator_kinematics_model::getJacobianMatrix(
    mat<double, mat_structure::rectangular>& Jac) const {
  manip_kin_mdl_jac_calculator(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getJacobianMatrix(Jac);
}

void manipulator_kinematics_model::getJacobianMatrixAndDerivative(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>& JacDot) const {
  manip_kin_mdl_jac_calculator(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getJacobianMatrixAndDerivative(Jac, JacDot);
}

vect_n<double> manipulator_kinematics_model::getJointPositions() const {
  vect_n<double> result(getJointPositionsCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getJointPositions(&result[0]);

  return result;
}

void manipulator_kinematics_model::setJointPositions(
    const vect_n<double>& aJointPositions) {
  if (aJointPositions.size() != getJointPositionsCount()) {
    throw std::range_error("Joint-position vector has incorrect dimensions!");
  }

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .setJointPositions(&aJointPositions[0]);
}

vect_n<double> manipulator_kinematics_model::getJointVelocities() const {
  vect_n<double> result(getJointVelocitiesCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getJointVelocities(&result[0]);

  return result;
}

void manipulator_kinematics_model::setJointVelocities(
    const vect_n<double>& aJointVelocities) {
  if (aJointVelocities.size() != getJointVelocitiesCount()) {
    throw std::range_error("Joint-velocity vector has incorrect dimensions!");
  }

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .setJointVelocities(&aJointVelocities[0]);
}

vect_n<double> manipulator_kinematics_model::getJointAccelerations() const {
  vect_n<double> result(getJointAccelerationsCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getJointAccelerations(&result[0]);

  return result;
}

void manipulator_kinematics_model::setJointAccelerations(
    const vect_n<double>& aJointAccelerations) {
  if (aJointAccelerations.size() != getJointAccelerationsCount()) {
    throw std::range_error(
        "Joint-acceleration vector has incorrect dimensions!");
  }

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .setJointAccelerations(&aJointAccelerations[0]);
}

vect_n<double> manipulator_kinematics_model::getDependentPositions() const {
  vect_n<double> result(getDependentPositionsCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getDependentPositions(&result[0]);

  return result;
}

vect_n<double> manipulator_kinematics_model::getDependentVelocities() const {
  vect_n<double> result(getDependentVelocitiesCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getDependentVelocities(&result[0]);

  return result;
}

vect_n<double> manipulator_kinematics_model::getDependentAccelerations() const {
  vect_n<double> result(getDependentAccelerationsCount());

  manip_kin_mdl_joint_io(
      std::shared_ptr<const direct_kinematics_model>(this, null_deleter()))
      .getDependentAccelerations(&result[0]);

  return result;
}
}  // namespace ReaK::kte
