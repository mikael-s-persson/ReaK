/**
 * \file manip_SSRMS_arm.h
 *
 * This library declares a class to represent a kte-based model of a SSRMS manipulator in 3D, i.e.,
 * the Canadarm-2 manipulator.
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

#ifndef REAK_MBD_MODELS_MANIP_SSRMS_ARM_H_
#define REAK_MBD_MODELS_MANIP_SSRMS_ARM_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/mbd/kte/kte_map_chain.h"

#include "ReaK/mbd/models/inverse_kinematics_model.h"

namespace ReaK::kte {

/**
 * This class to represent a kte-based model of a SSRMS manipulator in 3D, i.e.,
 * the Canadarm-2 manipulator.
 */
class manip_SSRMS_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::vector<std::shared_ptr<gen_coord<double>>> m_joints;
  std::shared_ptr<joint_dependent_frame_3D> m_EE;
  vect_n<double> link_lengths;
  vect_n<double> joint_offsets;
  std::shared_ptr<kte_map_chain> m_chain;

 public:
  static constexpr std::size_t degrees_of_freedom = 7;

  vect_n<double> joint_lower_bounds;
  vect_n<double> joint_upper_bounds;

  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  /**
   * Default constructor.
   */
  explicit manip_SSRMS_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_3D<double>> aBaseFrame =
          std::shared_ptr<frame_3D<double>>(),
      const vect_n<double>& aLinkLengths = vect_n<double>(0.0, 0.380, 6.850,
                                                          6.850, 0.380, 0.0),
      const vect_n<double>& aJointOffsets = vect_n<double>(0.0, 0.635, 0.504,
                                                           0.504, 0.504, 0.635,
                                                           0.0),
      const vect_n<double>& aJointLowerBounds =
          vect_n<double>(-1.5 * M_PI, -1.5 * M_PI, -1.5 * M_PI, -1.5 * M_PI,
                         -1.5 * M_PI, -1.5 * M_PI, -1.5 * M_PI),
      const vect_n<double>& aJointUpperBounds =
          vect_n<double>(1.5 * M_PI, 1.5 * M_PI, 1.5 * M_PI, 1.5 * M_PI,
                         1.5 * M_PI, 1.5 * M_PI, 1.5 * M_PI));

  ~manip_SSRMS_kinematics() override = default;
  ;

  std::size_t getJointPositionsCount() const override { return 7; }

  std::size_t getJointVelocitiesCount() const override { return 7; }

  std::size_t getJointAccelerationsCount() const override { return 7; }

  std::size_t getDependentPositionsCount() const override { return 7; }

  std::size_t getDependentVelocitiesCount() const override { return 6; }

  std::size_t getDependentAccelerationsCount() const override { return 6; }

  std::size_t getCoordsCount() const override { return 7; }

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

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_SSRMS_kinematics, 0xC2100056, 1,
                              "manip_SSRMS_kinematics",
                              inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_MODELS_MANIP_SSRMS_ARM_H_
