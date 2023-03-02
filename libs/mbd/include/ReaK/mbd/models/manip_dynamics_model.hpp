/**
 * \file manip_dynamics_model.hpp
 *
 *
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date September 2010
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

#ifndef REAK_MANIP_DYNAMICS_MODEL_HPP
#define REAK_MANIP_DYNAMICS_MODEL_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/integrators/integrator.hpp>
#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <ReaK/mbd/kte/kte_system_input.hpp>
#include <ReaK/mbd/kte/kte_system_output.hpp>
#include <ReaK/mbd/kte/mass_matrix_calculator.hpp>

#include "inverse_dynamics_model.hpp"
#include "manip_kinematics_model.hpp"

#include <vector>

namespace ReaK::kte {

/**
 * This class stores the required information to represent the dynamics model of a manipulator.
 * Here, a manipulator is defined as a kinematic chain with "joint" coordinates (or frames) and
 * "dependent" coordinates (or frames). For example, a typical serial manipulator
 * could have a set of generalized coordinates (joint coordinates) as well as one or more frames
 * for the end-effector(s) (or additional link motions). This class is basically used to
 * regroup all that information and provides a certain number of functions related to the
 * use of a manipulator model (like computing mass-matrix). Additionally, system inputs and outputs
 * can also be registered (such as joint driver inputs, or state-measurement outputs).
 */
class manipulator_dynamics_model : public manipulator_kinematics_model,
                                   public inverse_dynamics_model,
                                   public state_rate_function_with_io<double> {
 protected:
  mass_matrix_calc mMassCalc;  ///< Holds the model's mass-matrix calculator.

  std::vector<std::shared_ptr<system_input>>
      mInputs;  ///< Holds the list of system input objects that are part of the KTE model.
  std::vector<std::shared_ptr<system_output>>
      mOutputs;  ///< Holds the list of system output objects that are part of the KTE model.

 public:
  using manipulator_kinematics_model::operator<<;

  /**
   * Default constructor.
   */
  explicit manipulator_dynamics_model(const std::string& aName = "")
      : manipulator_kinematics_model(aName),
        inverse_dynamics_model(aName),
        mMassCalc(aName + "_mass_calc") {}

  /**
   * Default destructor.
   */
  ~manipulator_dynamics_model() override = default;

  /**
   * Gets the manipulator KTE mass-matrix calculator used by this object.
   * \return The manipulator KTE mass-matrix calculator used by this object.
   */
  const mass_matrix_calc& getMassCalc() const { return mMassCalc; }

  /**
   * Add a system generalized coordinate.
   * \param aCoord a system generalized coordinate to add.
   * \return reference to this.
   */
  manipulator_kinematics_model& operator<<(
      const std::shared_ptr<gen_coord<double>>& aCoord) override;

  /**
   * Add a system 2D frame.
   * \param aFrame2D a system 2D frame to add.
   * \return reference to this.
   */
  manipulator_kinematics_model& operator<<(
      const std::shared_ptr<frame_2D<double>>& aFrame2D) override;

  /**
   * Add a system 3D frame.
   * \param aFrame3D a system 3D frame to add.
   * \return reference to this.
   */
  manipulator_kinematics_model& operator<<(
      const std::shared_ptr<frame_3D<double>>& aFrame3D) override;

  /**
   * Add a generalized inertial element.
   * \param aInertiaGen a generalized inertial element to add.
   * \return reference to this.
   */
  virtual manipulator_dynamics_model& operator<<(
      const std::shared_ptr<inertia_gen>& aInertiaGen);

  /**
   * Add a 2D inertial element.
   * \param aInertia2D a 2D inertial element to add.
   * \return reference to this.
   */
  virtual manipulator_dynamics_model& operator<<(
      const std::shared_ptr<inertia_2D>& aInertia2D);

  /**
   * Add a 3D inertial element.
   * \param aInertia3D a 3D inertial element to add.
   * \return reference to this.
   */
  virtual manipulator_dynamics_model& operator<<(
      const std::shared_ptr<inertia_3D>& aInertia3D);

  /**
   * Add a system input.
   * \param aInput A KTE system input to add.
   * \return reference to this.
   */
  virtual manipulator_dynamics_model& operator<<(
      const std::shared_ptr<system_input>& aInput);

  /**
   * Add a system output.
   * \param aOutput A KTE system output to add.
   * \return reference to this.
   */
  virtual manipulator_dynamics_model& operator<<(
      const std::shared_ptr<system_output>& aOutput);

  /**
   * Get the total number of state values for all the joint frames concatenated.
   * \return The total number of state values for all the joint frames concatenated.
   */
  std::size_t getJointStatesCount() const override {
    return getJointPositionsCount() + getJointVelocitiesCount();
  }

  /**
   * Get the total number of state values for all the dependent frames concatenated.
   * \return The total number of state values for all the dependent frames concatenated.
   */
  std::size_t getDependentStatesCount() const {
    return getDependentPositionsCount() + getDependentVelocitiesCount();
  }

  /**
   * Get the total number of KTE system inputs, if all concatenated in one vector.
   * \return The total number of KTE system inputs, if all concatenated in one vector.
   */
  std::size_t getInputsCount() const {
    std::size_t result = 0;
    for (const auto& mInput : mInputs) {
      result += mInput->getInputCount();
    }
    return result;
  }

  /**
   * Get the total number of KTE system outputs, if all concatenated in one vector.
   * \return The total number of KTE system outputs, if all concatenated in one vector.
   */
  std::size_t getOutputsCount() const {
    std::size_t result = 0;
    for (const auto& mOutput : mOutputs) {
      result += mOutput->getOutputCount();
    }
    return result;
  }

  /**
   * Obtain a vector that contains all the joint states concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized States, 2D States, 3D States ),
   * where the joints are sorted in the same order as in the container returned by Coords(),
   * Frames2D() and Frames3D(), respectively. In other words, the first 2 * Coords().size() elements
   * are the joint states of generalized coordinate joints, the next 7 * Frames2D().size()
   * elements are the states (and angular states) of the 2D frame joints, and the final 13 * Frames3D().size()
   * are the states (and angular states) of the 3D frame joints.
   * \return All the joint states concatenated into one vector.
   */
  vect_n<double> getJointStates() const override;

  /**
   * Set all the joint states of the manipulator to a vector of concatenated joint-states.
   * The ordering in the vector is as follows: ( Generalized States, 2D States, 3D States ),
   * where the joints are sorted in the same order as in the container returned by Coords(),
   * Frames2D() and Frames3D(), respectively. In other words, the first 2 * Coords().size() elements
   * are the joint states of generalized coordinate joints, the next 7 * Frames2D().size()
   * elements are the states (and angular states) of the 2D frame joints, and the final 13 * Frames3D().size()
   * are the states (and angular states) of the 3D frame joints.
   * \param aJointAccelerations All the joint states concatenated into one vector.
   */
  void setJointStates(const vect_n<double>& aJointStates) override;

  /**
   * This function computes the output-vector corresponding to a state vector.
   *
   * \param aTime current integration time
   * \param aState current state vector
   * \param aOutput holds, as output, the output-vector
   */
  void computeOutput(double aTime, const ReaK::vect_n<double>& aState,
                     ReaK::vect_n<double>& aOutput) override;

  /**
   * This function sets the input-vector.
   *
   * \param aInput current input-vector
   */
  void setInput(const ReaK::vect_n<double>& aInput) override;

  /**
   * This function gets the currently set input-vector for the system.
   *
   * \return The currently set input-vector for the system.
   */
  vect_n<double> getInput() const;

  /**
   * Computes the time-derivative of the state-vector of all the joints concatenated into one vector.
   * The vector of state-derivatives corresponds to the vector obtained from getJointStates().
   * \param aTime current integration time
   * \param aState current state vector
   * \param aStateRate holds, as output, the time-derivative of the state vector
   */
  void computeStateRate(double aTime, const ReaK::vect_n<double>& aState,
                        ReaK::vect_n<double>& aStateRate) override;

  /**
   * Obtain a vector that contains all the dependent states concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependent frames are sorted in the same order as in the container returned by DependentCoords(),
   * DependentFrames2D() and DependentFrames3D(), respectively. In other words, the first 2 * DependentCoords().size()
   * elements
   * are the states of dependent generalized coordinates, the next 7 * DependentFrames2D().size()
   * elements are the states of the 2D dependent frames, and the final 13 * DependentFrames3D().size()
   * are the states of the 3D dependent frames.
   * \return All the dependent states concatenated into one vector.
   */
  vect_n<double> getDependentStates() const;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mModel) {
      mModel->doMotion(aFlag, aStorage);
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {
    if (mModel) {
      mModel->doForce(aFlag, aStorage);
    }
  }

  void clearForce() override {
    if (mModel) {
      mModel->clearForce();
    }
  }

  /**
   * Get the mass matrix for the system.
   * \param M stores, as output, the calculated system's mass-matrix.
   */
  void getMassMatrix(mat<double, mat_structure::symmetric>& M) override;

  /**
   * Get the mass matrix for the system and its time-derivative.
   * \param M stores, as output, the calculated system's mass-matrix.
   * \param M_dot stores, as output, the calculated time-derivative of the system's mass matrix.
   */
  void getMassMatrixAndDerivative(
      mat<double, mat_structure::symmetric>& M,
      mat<double, mat_structure::square>& M_dot) override;

  /**
   * Get the twist-shaping matrix, the block-diagonal, link mass-matrix, and the time-derivative of the twist-shaping
   * matrix.
   * \param Tcm stores, as output, the calculated system's twist-shaping matrix.
   * \param Mcm stores, as output, the calculated block-diagonal, link mass matrix.
   * \param Tcm_dot storse, as output, the calculated time-derivative of the system's twist-shaping matrix.
   */
  void get_TMT_TdMT(mat<double, mat_structure::rectangular>& Tcm,
                    mat<double, mat_structure::symmetric>& Mcm,
                    mat<double, mat_structure::rectangular>& Tcm_dot) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    manipulator_kinematics_model::save(
        A, manipulator_kinematics_model::getStaticObjectType()->TypeVersion());
    state_rate_function<double>::save(
        A, state_rate_function<double>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mMassCalc);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    manipulator_kinematics_model::load(
        A, manipulator_kinematics_model::getStaticObjectType()->TypeVersion());
    state_rate_function<double>::load(
        A, state_rate_function<double>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mMassCalc);
  }

  RK_RTTI_MAKE_CONCRETE_3BASE(manipulator_dynamics_model, 0xC210004E, 1,
                              "manipulator_dynamics_model",
                              manipulator_kinematics_model,
                              inverse_dynamics_model,
                              state_rate_function<double>)
};

}  // namespace ReaK::kte

#endif
