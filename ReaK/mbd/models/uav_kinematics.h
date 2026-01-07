/**
 * \file uav_kinematics.h
 *
 * This library declares a class to represent a kte-based model of the kinematics a UAV.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
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

#ifndef REAK_MBD_MODELS_UAV_KINEMATICS_H_
#define REAK_MBD_MODELS_UAV_KINEMATICS_H_

#include <memory>
#include "ReaK/mbd/kte/kte_map_chain.h"
#include "ReaK/mbd/models/inverse_kinematics_model.h"

namespace ReaK::kte {

/**
 * This class to represent a kte-based model of a UAV.
 */
class UAV_kinematics : public inverse_kinematics_model {
 private:
  std::shared_ptr<frame_3D<double>> m_base_frame;
  std::shared_ptr<frame_3D<double>> m_motion_frame;
  std::shared_ptr<joint_dependent_frame_3D> m_output_frame;
  std::shared_ptr<kte_map_chain> m_chain;

 public:
  std::shared_ptr<kte_map_chain> getKTEChain() const override {
    return m_chain;
  }

  /**
   * Default constructor.
   */
  explicit UAV_kinematics(const std::string& aName = "",
                          std::shared_ptr<frame_3D<double>> aBaseFrame =
                              std::make_shared<frame_3D<double>>());

  ~UAV_kinematics() override = default;
  ;

  std::size_t getJointPositionsCount() const override { return 7; }

  std::size_t getJointVelocitiesCount() const override { return 6; }

  std::size_t getJointAccelerationsCount() const override { return 6; }

  std::size_t getDependentPositionsCount() const override { return 7; }

  std::size_t getDependentVelocitiesCount() const override { return 6; }

  std::size_t getDependentAccelerationsCount() const override { return 6; }

  std::size_t getFrames3DCount() const override { return 1; }

  std::shared_ptr<frame_3D<double>> getFrame3D(
      std::size_t /*i*/) const override {
    return m_motion_frame;
  }

  std::size_t getDependentFrames3DCount() const override { return 1; }

  std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const override {
    return m_output_frame;
  }

  void doDirectMotion() override;

  void doInverseMotion() override;

  void getJacobianMatrix(
      mat<double, mat_structure::rectangular>& Jac) const override;

  void getJacobianMatrixAndDerivative(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>& JacDot) const override;

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

  RK_RTTI_MAKE_CONCRETE_1BASE(UAV_kinematics, 0xC2100058, 1, "UAV_kinematics",
                              inverse_kinematics_model)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_MODELS_UAV_KINEMATICS_H_
