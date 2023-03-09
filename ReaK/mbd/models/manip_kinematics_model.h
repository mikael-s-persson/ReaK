/**
 * \file manip_kinematics_model.h
 *
 * This library declares classes to represent manipulator kinematic systems. Essentially,
 * the model of the manipulator is only a KTE chain provided by the user, but these
 * manipulator-model classes take care of grouping the joints, their limits, and their
 * jacobian matrices. The jacobian matrices are computed from jacobian mappings of up-stream
 * joints, for each output-frame to determine the twist-shaping matrices.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date September 2010
 */

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

#ifndef REAK_MBD_MODELS_MANIP_KINEMATICS_MODEL_H_
#define REAK_MBD_MODELS_MANIP_KINEMATICS_MODEL_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map_chain.h"
#include "ReaK/mbd/models/direct_kinematics_model.h"

#include <vector>

namespace ReaK::kte {

/**
 * This class stores the required information to represent the kinematic model of a manipulator.
 * Here, a manipulator is defined as a kinematic chain with "input" coordinates (or frames) and
 * "output" coordinates (or frames). For example, a typical serial manipulator
 * could have a set of generalized coordinates (joint coordinates) as well as one or more frames
 * for the end-effector(s) (or additional link motions). This class is basically used to
 * regroup all that information and provides a certain number of functions related to the
 * use of a manipulator model (like computing jacobians).
 */
class manipulator_kinematics_model : public direct_kinematics_model {
 protected:
  std::vector<std::shared_ptr<gen_coord<double>>>
      mCoords;  ///< Holds the list of generalized coordinates in the system.
  std::vector<std::shared_ptr<frame_2D<double>>>
      mFrames2D;  ///< Holds the list of 2D coordinates frame in the system.
  std::vector<std::shared_ptr<frame_3D<double>>>
      mFrames3D;  ///< Holds the list of 3D coordinates frame in the system.

  std::vector<std::shared_ptr<joint_dependent_gen_coord>>
      mDependentGenCoords;  ///< Holds the list of dependent generalized coordinates.
  std::vector<std::shared_ptr<joint_dependent_frame_2D>>
      mDependent2DFrames;  ///< Holds the list of dependent 2D frames.
  std::vector<std::shared_ptr<joint_dependent_frame_3D>>
      mDependent3DFrames;  ///< Holds the list of dependent 3D frames.

  std::shared_ptr<kte_map_chain>
      mModel;  ///< Holds the model of the manipulator as a kte-chain.

 public:
  /**
   * Default constructor.
   */
  explicit manipulator_kinematics_model(const std::string& aName = "")
      : direct_kinematics_model(aName) {}

  /**
   * Default destructor.
   */
  ~manipulator_kinematics_model() override = default;
  ;

  /**
   * Sets the manipulator KTE model to use in this object.
   * \param aModel The manipulator KTE model to use in this object.
   */
  virtual void setModel(const std::shared_ptr<kte_map_chain>& aModel) {
    mModel = aModel;
  }

  /**
   * Gets the manipulator KTE model used by this object.
   * \return The manipulator KTE model used by this object.
   */
  std::shared_ptr<kte_map_chain> getModel() const { return mModel; }

  /**
   * Add a dependent generalized coordinate to the jacobian calculation.
   * \param aDependentGenCoord a dependent generalized coordinate to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<joint_dependent_gen_coord>& aDependentGenCoord);

  /**
   * Add a dependent 2D frame to the jacobian calculation.
   * \param aDependent2DFrame a dependent 2D frame  to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<joint_dependent_frame_2D>& aDependent2DFrame);

  /**
   * Add a dependent 3D frame to the jacobian calculation.
   * \param aDependent3DFrame a dependent 3D frame to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<joint_dependent_frame_3D>& aDependent3DFrame);

  /**
   * Add a system generalized coordinate.
   * \param aCoord a system generalized coordinate to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<gen_coord<double>>& aCoord);

  /**
   * Add a system 2D frame.
   * \param aFrame2D a system 2D frame to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<frame_2D<double>>& aFrame2D);

  /**
   * Add a system 3D frame.
   * \param aFrame3D a system 3D frame to add.
   * \return reference to this.
   */
  virtual manipulator_kinematics_model& operator<<(
      const std::shared_ptr<frame_3D<double>>& aFrame3D);

  /******************************************************************************************
   *  direct_kinematics_model: joint and dependent positions, velocities and accelerations  *
   ******************************************************************************************/

  std::size_t getJointPositionsCount() const override {
    return mCoords.size() + 4 * mFrames2D.size() + 7 * mFrames3D.size();
  }

  std::size_t getJointVelocitiesCount() const override {
    return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
  }

  std::size_t getJointAccelerationsCount() const override {
    return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
  }

  std::size_t getDependentPositionsCount() const override {
    return mDependentGenCoords.size() + 4 * mDependent2DFrames.size() +
           7 * mDependent3DFrames.size();
  }

  std::size_t getDependentVelocitiesCount() const override {
    return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() +
           6 * mDependent3DFrames.size();
  }

  std::size_t getDependentAccelerationsCount() const override {
    return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() +
           6 * mDependent3DFrames.size();
  }

  /*************************************************************************
   *  direct_kinematics_model: joint and dependent coordinates and frames  *
   *************************************************************************/

  std::size_t getCoordsCount() const override { return mCoords.size(); }

  std::shared_ptr<gen_coord<double>> getCoord(std::size_t i) const override {
    return mCoords[i];
  }

  std::size_t getFrames2DCount() const override { return mFrames2D.size(); }

  std::shared_ptr<frame_2D<double>> getFrame2D(std::size_t i) const override {
    return mFrames2D[i];
  }

  std::size_t getFrames3DCount() const override { return mFrames3D.size(); }

  std::shared_ptr<frame_3D<double>> getFrame3D(std::size_t i) const override {
    return mFrames3D[i];
  }

  std::size_t getDependentCoordsCount() const override {
    return mDependentGenCoords.size();
  }

  std::shared_ptr<joint_dependent_gen_coord> getDependentCoord(
      std::size_t i) const override {
    return mDependentGenCoords[i];
  }

  std::size_t getDependentFrames2DCount() const override {
    return mDependent2DFrames.size();
  }

  std::shared_ptr<joint_dependent_frame_2D> getDependentFrame2D(
      std::size_t i) const override {
    return mDependent2DFrames[i];
  }

  std::size_t getDependentFrames3DCount() const override {
    return mDependent3DFrames.size();
  }

  std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const override {
    return mDependent3DFrames[i];
  }

  /************************************************************************
   *  direct_kinematics_model: kinematic functions (motion and Jacobian)  *
   ************************************************************************/

  void doDirectMotion() override {
    if (mModel) {
      mModel->doMotion();
    }
  }

  void getJacobianMatrix(
      mat<double, mat_structure::rectangular>& Jac) const override;

  void getJacobianMatrixAndDerivative(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>& JacDot) const override;

  /***********************************************************************************
   *  direct_kinematics_model: contatenated positions, velocities and accelerations  *
   ***********************************************************************************/

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

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    direct_kinematics_model::save(
        A, direct_kinematics_model::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mCoords) & RK_SERIAL_SAVE_WITH_NAME(mFrames2D) &
        RK_SERIAL_SAVE_WITH_NAME(mFrames3D) &
        RK_SERIAL_SAVE_WITH_NAME(mDependentGenCoords) &
        RK_SERIAL_SAVE_WITH_NAME(mDependent2DFrames) &
        RK_SERIAL_SAVE_WITH_NAME(mDependent3DFrames) &
        RK_SERIAL_SAVE_WITH_NAME(mModel);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    direct_kinematics_model::load(
        A, direct_kinematics_model::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mCoords) & RK_SERIAL_LOAD_WITH_NAME(mFrames2D) &
        RK_SERIAL_LOAD_WITH_NAME(mFrames3D) &
        RK_SERIAL_LOAD_WITH_NAME(mDependentGenCoords) &
        RK_SERIAL_LOAD_WITH_NAME(mDependent2DFrames) &
        RK_SERIAL_LOAD_WITH_NAME(mDependent3DFrames) &
        RK_SERIAL_LOAD_WITH_NAME(mModel);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(manipulator_kinematics_model, 0xC210004D, 1,
                              "manipulator_kinematics_model",
                              direct_kinematics_model)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_MODELS_MANIP_KINEMATICS_MODEL_H_
