/**
 * \file manip_3R_arm.hpp
 *
 * This library declares a class to represent a kte-based model of a 3R manipulator in 2D or 3D.
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

#ifndef REAK_MANIP_3R_ARM_HPP
#define REAK_MANIP_3R_ARM_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/mbd/kte/kte_map_chain.hpp"
#include "ReaK/mbd/models/inverse_kinematics_model.hpp"

namespace ReaK::kte {

/**
 * This class that models a 2D manipulator with 3 revolute joints (i.e., shoulder-elbow-wrist joints).
 * This class is only a kinematics model.
 */
class manip_3R_2D_kinematics : public inverse_kinematics_model {
 protected:
  std::shared_ptr<frame_2D<double>> m_base_frame;
  std::vector<std::shared_ptr<gen_coord<double>>> m_joints;
  std::shared_ptr<joint_dependent_frame_2D> m_EE;
  double m_link1_length;
  double m_link2_length;
  double m_link3_length;

  std::shared_ptr<kte_map_chain> m_chain;

 public:
  static constexpr std::size_t degrees_of_freedom = 3;

  /**
   * Default constructor.
   */
  explicit manip_3R_2D_kinematics(const std::string& aName = "",
                                  std::shared_ptr<frame_2D<double>> aBaseFrame =
                                      std::shared_ptr<frame_2D<double>>(),
                                  double aLink1Length = 1.0,
                                  double aLink2Length = 1.0,
                                  double aLink3Length = 0.0);

  ~manip_3R_2D_kinematics() override = default;
  ;

  std::size_t getJointPositionsCount() const override { return 3; }

  std::size_t getJointVelocitiesCount() const override { return 3; }

  std::size_t getJointAccelerationsCount() const override { return 3; }

  std::size_t getDependentPositionsCount() const override { return 4; }

  std::size_t getDependentVelocitiesCount() const override { return 3; }

  std::size_t getDependentAccelerationsCount() const override { return 3; }

  std::size_t getCoordsCount() const override { return 3; }

  std::shared_ptr<gen_coord<double>> getCoord(std::size_t i) const override {
    return m_joints[i];
  }

  std::size_t getDependentFrames2DCount() const override { return 1; }

  std::shared_ptr<joint_dependent_frame_2D> getDependentFrame2D(
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
    return {-std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
  }

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override {}

  vect_n<double> getJointPositionUpperBounds() const override {
    return {std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
  }

  void setJointPositionUpperBounds(
      const vect_n<double>& aJointUpperBounds) override {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_3R_2D_kinematics, 0xC2100052, 1,
                              "manip_3R_2D_kinematics",
                              inverse_kinematics_model)
};

/**
 * This class that models a 3D manipulator with 3 revolute joints (i.e., shoulder-elbow-wrist joints
 * all aligned along the z-axis). This class is only a kinematics model.
 */
class manip_3R_3D_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::vector<std::shared_ptr<gen_coord<double>>> m_joints;
  std::shared_ptr<joint_dependent_frame_3D> m_EE;
  double m_link1_length;
  double m_link1_dz;
  double m_link2_length;
  double m_link2_dz;
  double m_link3_length;
  double m_link3_dz;

  std::shared_ptr<kte_map_chain> m_chain;

 public:
  /**
   * Default constructor.
   */
  explicit manip_3R_3D_kinematics(const std::string& aName = "",
                                  std::shared_ptr<frame_3D<double>> aBaseFrame =
                                      std::shared_ptr<frame_3D<double>>(),
                                  double aLink1Length = 1.0,
                                  double aLink1DZ = 0.0,
                                  double aLink2Length = 1.0,
                                  double aLink2DZ = 0.0,
                                  double aLink3Length = 0.0,
                                  double aLink3DZ = 0.0);

  ~manip_3R_3D_kinematics() override = default;
  ;

  std::size_t getJointPositionsCount() const override { return 3; }

  std::size_t getJointVelocitiesCount() const override { return 3; }

  std::size_t getJointAccelerationsCount() const override { return 3; }

  std::size_t getDependentPositionsCount() const override { return 7; }

  std::size_t getDependentVelocitiesCount() const override { return 6; }

  std::size_t getDependentAccelerationsCount() const override { return 6; }

  std::size_t getCoordsCount() const override { return 3; }

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
    return {-std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
  }

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override {}

  vect_n<double> getJointPositionUpperBounds() const override {
    return {std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
  }

  void setJointPositionUpperBounds(
      const vect_n<double>& aJointUpperBounds) override {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(manip_3R_3D_kinematics, 0xC2100053, 1,
                              "manip_3R_3D_kinematics",
                              inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif
