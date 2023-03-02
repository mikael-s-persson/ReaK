/**
 * \file free_floating_platform.hpp
 *
 * This library declares a class to represent the kinematics model of a free-floating platform in 2D or 3D.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date December 2013
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

#ifndef REAK_FREE_FLOATING_PLATFORM_HPP
#define REAK_FREE_FLOATING_PLATFORM_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/mbd/kte/kte_map_chain.hpp>
#include "inverse_kinematics_model.hpp"

namespace ReaK::kte {

/**
 * This class that models a 2D free-floating platform.
 * This class is only a kinematics model (which is rather trivial).
 */
class free_floater_2D_kinematics : public inverse_kinematics_model {
 protected:
  std::shared_ptr<frame_2D<double>> m_base_frame;
  std::shared_ptr<frame_2D<double>> m_state_frame;
  std::shared_ptr<jacobian_2D_2D<double>> m_state_jacobian;
  std::shared_ptr<frame_2D<double>> m_output_frame;
  mutable std::vector<std::shared_ptr<joint_dependent_frame_2D>> m_EEs;
  std::vector<pose_2D<double>> m_EEposes;  // relative to output-frame

  std::shared_ptr<kte_map_chain> m_chain;

  void resyncEndEffectors() const;
  void getJacobianMatrixAndDerivativeImpl(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>* JacDot) const;

 public:
  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  static constexpr std::size_t degrees_of_freedom = 3;

  std::shared_ptr<frame_2D<double>> getBaseFrame() const {
    return m_base_frame;
  }
  void setBaseFrame(const std::shared_ptr<frame_2D<double>>& aBaseFrame) {
    m_base_frame = aBaseFrame;
    m_output_frame->Parent = aBaseFrame;
  }

  std::shared_ptr<frame_2D<double>> getStateFrame() const {
    return m_state_frame;
  }
  void setStateFrame(const std::shared_ptr<frame_2D<double>>& aStateFrame) {
    for (auto& m_EE : m_EEs) {
      m_EE->mUpStream2DJoints.erase(m_state_frame);
      m_EE->mUpStream2DJoints[aStateFrame] = m_state_jacobian;
    }
    m_state_frame = aStateFrame;
  }

  std::shared_ptr<frame_2D<double>> getOutputFrame() const {
    return m_output_frame;
  }

  const std::vector<pose_2D<double>>& getEEPoses() const { return m_EEposes; }
  std::vector<pose_2D<double>>& getEEPoses() { return m_EEposes; }

  /**
   * Default constructor.
   */
  explicit free_floater_2D_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_2D<double>> aBaseFrame =
          std::shared_ptr<frame_2D<double>>());

  ~free_floater_2D_kinematics() override = default;

  std::size_t getJointPositionsCount() const override { return 4; }

  std::size_t getJointVelocitiesCount() const override { return 3; }

  std::size_t getJointAccelerationsCount() const override { return 3; }

  std::size_t getDependentPositionsCount() const override {
    return 4 * m_EEposes.size();
  }

  std::size_t getDependentVelocitiesCount() const override {
    return 3 * m_EEposes.size();
  }

  std::size_t getDependentAccelerationsCount() const override {
    return 3 * m_EEposes.size();
  }

  std::size_t getFrames2DCount() const override { return 1; }

  std::shared_ptr<frame_2D<double>> getFrame2D(std::size_t i) const override {
    return m_state_frame;
  }

  std::size_t getDependentFrames2DCount() const override {
    return m_EEposes.size();
  }

  std::shared_ptr<joint_dependent_frame_2D> getDependentFrame2D(
      std::size_t i) const override {
    resyncEndEffectors();
    return m_EEs[i];
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
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
  }

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override {}

  vect_n<double> getJointPositionUpperBounds() const override {
    return {std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
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

  RK_RTTI_MAKE_CONCRETE_1BASE(free_floater_2D_kinematics, 0xC210005B, 1,
                              "free_floater_2D_kinematics",
                              inverse_kinematics_model)
};

/**
 * This class that models a 3D manipulator with 3 revolute joints (i.e., shoulder-elbow-wrist joints
 * all aligned along the z-axis). This class is only a kinematics model.
 */
class free_floater_3D_kinematics : public inverse_kinematics_model {
 protected:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::shared_ptr<frame_3D<double>> m_state_frame;
  std::shared_ptr<jacobian_3D_3D<double>> m_state_jacobian;
  std::shared_ptr<frame_3D<double>> m_output_frame;
  mutable std::vector<std::shared_ptr<joint_dependent_frame_3D>> m_EEs;
  std::vector<pose_3D<double>> m_EEposes;  // relative to output-frame

  std::shared_ptr<kte_map_chain> m_chain;

  void resyncEndEffectors() const;
  void getJacobianMatrixAndDerivativeImpl(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>* JacDot) const;

 public:
  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  static constexpr std::size_t degrees_of_freedom = 6;

  std::shared_ptr<frame_3D<double>> getBaseFrame() const {
    return m_base_frame;
  }
  void setBaseFrame(const std::shared_ptr<frame_3D<double>>& aBaseFrame) {
    m_base_frame = aBaseFrame;
    m_output_frame->Parent = aBaseFrame;
  }

  std::shared_ptr<frame_3D<double>> getStateFrame() const {
    return m_state_frame;
  }
  void setStateFrame(const std::shared_ptr<frame_3D<double>>& aStateFrame) {
    for (auto& m_EE : m_EEs) {
      m_EE->mUpStream3DJoints.erase(m_state_frame);
      m_EE->mUpStream3DJoints[aStateFrame] = m_state_jacobian;
    }
    m_state_frame = aStateFrame;
  }

  std::shared_ptr<frame_3D<double>> getOutputFrame() const {
    return m_output_frame;
  }

  const std::vector<pose_3D<double>>& getEEPoses() const { return m_EEposes; }
  std::vector<pose_3D<double>>& getEEPoses() { return m_EEposes; }

  /**
   * Default constructor.
   */
  explicit free_floater_3D_kinematics(
      const std::string& aName = "",
      std::shared_ptr<frame_3D<double>> aBaseFrame =
          std::shared_ptr<frame_3D<double>>());

  ~free_floater_3D_kinematics() override = default;
  ;

  std::size_t getJointPositionsCount() const override { return 7; }

  std::size_t getJointVelocitiesCount() const override { return 6; }

  std::size_t getJointAccelerationsCount() const override { return 6; }

  std::size_t getDependentPositionsCount() const override {
    return 7 * m_EEposes.size();
  }

  std::size_t getDependentVelocitiesCount() const override {
    return 6 * m_EEposes.size();
  }

  std::size_t getDependentAccelerationsCount() const override {
    return 6 * m_EEposes.size();
  }

  std::size_t getFrames3DCount() const override { return 1; }

  std::shared_ptr<frame_3D<double>> getFrame3D(std::size_t i) const override {
    return m_state_frame;
  }

  std::size_t getDependentFrames3DCount() const override {
    return m_EEposes.size();
  }

  std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const override {
    resyncEndEffectors();
    return m_EEs[i];
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
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
  }

  void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) override {}

  vect_n<double> getJointPositionUpperBounds() const override {
    return {std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
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

  RK_RTTI_MAKE_CONCRETE_1BASE(free_floater_3D_kinematics, 0xC210005C, 1,
                              "free_floater_3D_kinematics",
                              inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif
