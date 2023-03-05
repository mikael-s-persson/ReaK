/**
 * \file manip_P3R3R_arm.hpp
 *
 * This library declares a class to represent a kte-based model of a P-3R-3R manipulator in 3D, i.e.,
 * a P-3R-3R manipulator refers to a 6-dof manipulator in a decoupled architecture (2-dof shoulder + elbow + 3-dof
 *wrist)
 * mounted on a linear track (1-dof).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2013
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

#ifndef REAK_MANIP_P3R3R_ARM_HPP
#define REAK_MANIP_P3R3R_ARM_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/mbd/kte/kte_map_chain.hpp"
#include "ReaK/mbd/kte/prismatic_joint.hpp"

#include "ReaK/mbd/models/inverse_kinematics_model.hpp"
#include "ReaK/mbd/models/manip_3R3R_arm.hpp"

namespace ReaK::kte {

/**
 * This class that models a 3D manipulator with 6 revolute joints in a decoupled architecture
 * with a 2-dof shoulder, an elbow, and a 3-dof wrist, with the base mounted on a linear track.
 * This class is only a kinematics model.
 * \note In the zero-configuration (all joints at zero), the arm is pointing straight up (z-axis),
 * and, at that configuration, joints 1, 4 and 6 turn about the positive z-axis and joints 2, 3, and 5
 * turn about the negative y-axis, leading to the end-effector (flange) to have its local z-axis pointing
 * outwards from the flange and its local y-axis pointing to the side. The track is aligned with the
 * positive x-axis of the base-frame.
 */
class manip_P3R3R_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::shared_ptr<gen_coord<double>> m_track_coord;
  std::shared_ptr<frame_3D<double>> m_output_frame;
  std::shared_ptr<prismatic_joint_3D> m_track_joint;
  double m_track_lower_bound;
  double m_track_upper_bound;

  manip_3R3R_kinematics m_arm_model;
  std::shared_ptr<kte_map_chain> m_chain;

 public:
  static constexpr std::size_t degrees_of_freedom = 7;

  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  /**
   * Default constructor.
   */
  explicit manip_P3R3R_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_3D<double>> aBaseFrame =
          std::shared_ptr<frame_3D<double>>(),
      double aBaseToShoulder = 0.0, double aShoulderToElbow = 1.0,
      double aElbowToJoint4 = 0.5, double aJoint4ToWrist = 0.5,
      double aWristToFlange = 0.2,
      const vect_n<double>& aJointLowerBounds =
          vect_n<double>(0.0, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI),
      const vect_n<double>& aJointUpperBounds = vect_n<double>(1.0, M_PI, M_PI,
                                                               M_PI, M_PI, M_PI,
                                                               M_PI));

  ~manip_P3R3R_kinematics() override = default;

  std::size_t getJointPositionsCount() const override {
    return 1 + m_arm_model.getJointPositionsCount();
  }

  std::size_t getJointVelocitiesCount() const override {
    return 1 + m_arm_model.getJointVelocitiesCount();
  }

  std::size_t getJointAccelerationsCount() const override {
    return 1 + m_arm_model.getJointAccelerationsCount();
  }

  std::size_t getDependentPositionsCount() const override {
    return m_arm_model.getDependentPositionsCount();
  }

  std::size_t getDependentVelocitiesCount() const override {
    return m_arm_model.getDependentVelocitiesCount();
  }

  std::size_t getDependentAccelerationsCount() const override {
    return m_arm_model.getDependentAccelerationsCount();
  }

  std::size_t getCoordsCount() const override {
    return 1 + m_arm_model.getCoordsCount();
  }

  std::shared_ptr<gen_coord<double>> getCoord(std::size_t i) const override {
    if (i == 0) {
      return m_track_coord;
    }
    return m_arm_model.getCoord(i - 1);
  }

  std::size_t getDependentFrames3DCount() const override {
    return m_arm_model.getDependentFrames3DCount();
  }

  std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const override {
    return m_arm_model.getDependentFrame3D(i);
  }

  void doDirectMotion() override;

  void doInverseMotion() override;

  void getJacobianMatrix(
      mat<double, mat_structure::rectangular>& Jac) const override;

  void getJacobianMatrixAndDerivative(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>& JacDot) const override;

  vect_n<double> getJointPositionLowerBounds() const override;

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override;

  vect_n<double> getJointPositionUpperBounds() const override;

  void setJointPositionUpperBounds(
      const vect_n<double>& aJointUpperBounds) override;

  vect_n<double> getJointPositions() const override;

  void setJointPositions(const vect_n<double>& aJointPositions) override;

  vect_n<double> getJointVelocities() const override;

  void setJointVelocities(const vect_n<double>& aJointVelocities) override;

  vect_n<double> getJointAccelerations() const override;

  void setJointAccelerations(
      const vect_n<double>& aJointAccelerations) override;

  vect_n<double> getDependentPositions() const override {
    return m_arm_model.getDependentPositions();
  }

  vect_n<double> getDependentVelocities() const override {
    return m_arm_model.getDependentVelocities();
  }

  vect_n<double> getDependentAccelerations() const override {
    return m_arm_model.getDependentAccelerations();
  }

  void setDependentPositions(const vect_n<double>& aDepPositions) override {
    m_arm_model.setDependentPositions(aDepPositions);
  }

  void setDependentVelocities(const vect_n<double>& aDepVelocities) override {
    m_arm_model.setDependentVelocities(aDepVelocities);
  }

  void setDependentAccelerations(
      const vect_n<double>& aDepAccelerations) override {
    m_arm_model.setDependentAccelerations(aDepAccelerations);
  }

  void save(serialization::oarchive& A, unsigned int /*unused*/) const override;

  void load(serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_P3R3R_kinematics, 0xC2100059, 1,
                              "manip_P3R3R_kinematics",
                              inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif
