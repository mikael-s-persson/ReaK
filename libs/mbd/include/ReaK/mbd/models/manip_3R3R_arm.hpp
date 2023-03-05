/**
 * \file manip_3R3R_arm.hpp
 *
 * This library declares a class to represent a kte-based model of a 3R-3R manipulator in 3D, i.e.,
 * a 3R-3R manipulator refers to a 6-dof manipulator in a decoupled architecture (2-dof shoulder + elbow + 3-dof wrist).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date November 2012
 */

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

#ifndef REAK_MANIP_3R3R_ARM_HPP
#define REAK_MANIP_3R3R_ARM_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/mbd/kte/kte_map_chain.hpp"
#include "ReaK/mbd/models/inverse_kinematics_model.hpp"

namespace ReaK::kte {

/**
 * This class that models a 3D manipulator with 6 revolute joints in a decoupled architecture
 * with a 2-dof shoulder, an elbow, and a 3-dof wrist. This class is only a kinematics model.
 * \note In the zero-configuration (all joints at zero), the arm is pointing straight up (z-axis),
 * and, at that configuration, joints 1, 4 and 6 turn about the positive z-axis and joints 2, 3, and 5
 * turn about the negative y-axis, leading to the end-effector (flange) to have its local z-axis pointing
 * outwards from the flange and its local y-axis pointing to the side.
 */
class manip_3R3R_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::vector<std::shared_ptr<gen_coord<double>>> m_joints;
  std::shared_ptr<joint_dependent_frame_3D> m_EE;
  double base_to_shoulder;
  double shoulder_to_elbow;
  double elbow_to_joint_4;
  double joint_4_to_wrist;
  double wrist_to_flange;
  std::shared_ptr<kte_map_chain> m_chain;

 public:
  static constexpr std::size_t degrees_of_freedom = 6;

  vect_n<double> joint_lower_bounds;
  vect_n<double> joint_upper_bounds;

  double getBaseToShoulder() const { return base_to_shoulder; }
  double getShoulderToElbow() const { return shoulder_to_elbow; }
  double getElbowToJoint4() const { return elbow_to_joint_4; }
  double getJoint4ToWrist() const { return joint_4_to_wrist; }
  double getWristToFlange() const { return wrist_to_flange; }

  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  /**
   * Default constructor.
   */
  explicit manip_3R3R_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_3D<double>> aBaseFrame =
          std::shared_ptr<frame_3D<double>>(),
      double aBaseToShoulder = 0.0, double aShoulderToElbow = 1.0,
      double aElbowToJoint4 = 0.5, double aJoint4ToWrist = 0.5,
      double aWristToFlange = 0.2,
      const vect_n<double>& aJointLowerBounds = vect_n<double>(-M_PI, -M_PI,
                                                               -M_PI, -M_PI,
                                                               -M_PI, -M_PI),
      const vect_n<double>& aJointUpperBounds = vect_n<double>(M_PI, M_PI, M_PI,
                                                               M_PI, M_PI,
                                                               M_PI));

  ~manip_3R3R_kinematics() override = default;

  std::size_t getJointPositionsCount() const override { return 6; }

  std::size_t getJointVelocitiesCount() const override { return 6; }

  std::size_t getJointAccelerationsCount() const override { return 6; }

  std::size_t getDependentPositionsCount() const override { return 7; }

  std::size_t getDependentVelocitiesCount() const override { return 6; }

  std::size_t getDependentAccelerationsCount() const override { return 6; }

  std::size_t getCoordsCount() const override { return 6; }

  std::shared_ptr<gen_coord<double>> getCoord(std::size_t i) const override {
    return m_joints[i];
  }

  std::size_t getDependentFrames3DCount() const override { return 1; }

  std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const override {
    return m_EE;
  }

  void doDirectMotion() override;

  void doInverseMotion() override;

  void getJacobianMatrix(
      mat<double, mat_structure::rectangular>& Jac) const override;

  void getJacobianMatrixAndDerivative(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>& JacDot) const override;

  vect_n<double> getJointPositionLowerBounds() const override {
    return joint_lower_bounds;
  }

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override {
    joint_lower_bounds = aJointLowerBounds;
  }

  vect_n<double> getJointPositionUpperBounds() const override {
    return joint_upper_bounds;
  }

  void setJointPositionUpperBounds(
      const vect_n<double>& aJointUpperBounds) override {
    joint_upper_bounds = aJointUpperBounds;
  }

  vect_n<double> getJointPositions() const override;

  void setJointPositions(const vect_n<double>& aJointPositions) override;

  vect_n<double> getJointVelocities() const override;

  void setJointVelocities(const vect_n<double>& aJointVelocities) override;

  vect_n<double> getJointAccelerations() const override;

  void setJointAccelerations(
      const vect_n<double>& aJointAccelerations) override;

  vect_n<double> getDependentPositions() const override;

  vect_n<double> getDependentVelocities() const override;

  vect_n<double> getDependentAccelerations() const override;

  void setDependentPositions(const vect_n<double>& aDepPositions) override;

  void setDependentVelocities(const vect_n<double>& aDepVelocities) override;

  void setDependentAccelerations(
      const vect_n<double>& aDepAccelerations) override;

  void save(serialization::oarchive& A, unsigned int /*unused*/) const override;

  void load(serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_3R3R_kinematics, 0xC2100055, 1,
                              "manip_3R3R_kinematics", inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif
