/**
 * \file manip_ERA_arm.hpp
 *
 * This library declares a class to represent a kte-based model of a ERA manipulator, i.e.,
 * the European Robotic Arm manipulator.
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

#ifndef REAK_MANIP_ERA_ARM_HPP
#define REAK_MANIP_ERA_ARM_HPP

#include <memory>
#include "ReaK/core/base/defs.hpp"
#include "ReaK/mbd/kte/kte_map_chain.hpp"
#include "ReaK/mbd/models/inverse_kinematics_model.hpp"

namespace ReaK::kte {

/**
 * This class to represent a kte-based model of a ERA manipulator, i.e.,
 * the European Robotic Arm manipulator.
 */
class manip_ERA_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::vector<std::shared_ptr<gen_coord<double>>> m_joints;
  std::shared_ptr<joint_dependent_frame_3D> m_EE;
  vect_n<double> link_lengths;
  std::shared_ptr<kte_map_chain> m_chain;

 public:
  static constexpr std::size_t degrees_of_freedom = 7;

  vect<double, 3> preferred_elbow_dir;
  vect_n<double> joint_lower_bounds;
  vect_n<double> joint_upper_bounds;

  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  /**
   * Default constructor.
   */
  explicit manip_ERA_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_3D<double>> aBaseFrame =
          std::make_shared<frame_3D<double>>(),
      const vect_n<double>& aLinkLengths = vect_n<double>(1.228, 0.340, 4.072,
                                                          4.072, 0.340, 1.228),
      const vect<double, 3>& aPreferredElbowDir = (vect<double, 3>(0.0, 1.0,
                                                                   0.0)),
      const vect_n<double>& aJointLowerBounds =
          vect_n<double>(-M_PI, -0.75 * M_PI, -0.75 * M_PI, -0.75 * M_PI,
                         -0.75 * M_PI, -0.75 * M_PI, -M_PI),
      const vect_n<double>& aJointUpperBounds =
          vect_n<double>(M_PI, 0.75 * M_PI, 0.75 * M_PI, 0.75 * M_PI,
                         0.75 * M_PI, 0.75 * M_PI, M_PI));

  ~manip_ERA_kinematics() override = default;
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

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_ERA_kinematics, 0xC2100057, 1,
                              "manip_ERA_kinematics", inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif
