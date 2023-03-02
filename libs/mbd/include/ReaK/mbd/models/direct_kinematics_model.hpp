/**
 * \file direct_kinematics_model.hpp
 *
 * This library declares a base-class to represent classes that can be used to compute the
 * direct kinematics of some (kte-based) model.
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

#ifndef REAK_DIRECT_KINEMATICS_MODEL_HPP
#define REAK_DIRECT_KINEMATICS_MODEL_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>
#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/mbd/kte/jacobian_joint_map.hpp>
#include <ReaK/mbd/kte/kte_map_chain.hpp>

namespace ReaK::kte {

/**
 * This base-class represents a class that can be used to compute the direct kinematics
 * of some (kte-based) model. This class makes a distinction between joint coordinates
 * and dependent coordinates, where the former are generally coordinates (generalized,
 * and 2D/3D frames) determining the state of the joints of a kinematic chain, while the
 * latter are generally coordinates (generalized, and 2D/3D frames) whose state are dependent
 * on the joint states of a kinematic chain.
 */
class direct_kinematics_model : public virtual named_object {
 public:
  /**
   * Returns a pointer to the KTE-chain that represents this direct kinematics model.
   * \return A pointer to the KTE-chain that represents this direct kinematics model.
   */
  virtual std::shared_ptr<kte_map_chain> getKTEChain() const { return {}; }

  /**
   * Default constructor.
   */
  explicit direct_kinematics_model(const std::string& aName = "") {
    setName(aName);
  }

  /**
   * Default destructor.
   */
  ~direct_kinematics_model() override = default;

  /**
   * Get the total number of position values for all the joint frames concatenated.
   * \return The total number of position values for all the joint frames concatenated.
   */
  virtual std::size_t getJointPositionsCount() const { return 0; }

  /**
   * Get the total number of velocity values for all the joint frames concatenated.
   * \return The total number of velocity values for all the joint frames concatenated.
   */
  virtual std::size_t getJointVelocitiesCount() const { return 0; }

  /**
   * Get the total number of acceleration values for all the joint frames concatenated.
   * \return The total number of acceleration values for all the joint frames concatenated.
   */
  virtual std::size_t getJointAccelerationsCount() const { return 0; }

  /**
   * Get the total number of position values for all the dependent frames concatenated.
   * \return The total number of position values for all the dependent frames concatenated.
   */
  virtual std::size_t getDependentPositionsCount() const { return 0; }

  /**
   * Get the total number of velocity values for all the dependent frames concatenated.
   * \return The total number of velocity values for all the dependent frames concatenated.
   */
  virtual std::size_t getDependentVelocitiesCount() const { return 0; }

  /**
   * Get the total number of acceleration values for all the dependent frames concatenated.
   * \return The total number of acceleration values for all the dependent frames concatenated.
   */
  virtual std::size_t getDependentAccelerationsCount() const { return 0; }

  /**
   * Get the number of generalized joint coordinates.
   * \return The number of generalized joint coordinates.
   */
  virtual std::size_t getCoordsCount() const { return 0; }

  /**
   * Get access to one of the generalized joint coordinates.
   * \param i The index of the generalized joint coordinate.
   * \return Access to one of the generalized joint coordinates.
   */
  virtual std::shared_ptr<gen_coord<double>> getCoord(std::size_t i) const {
    return {};
  }

  /**
   * Get the number of joint 2D coordinate frames.
   * \return The number of joint 2D coordinate frames.
   */
  virtual std::size_t getFrames2DCount() const { return 0; }

  /**
   * Get access to one of the joint 2D coordinate frames.
   * \param i The index of the joint 2D coordinate frame.
   * \return Access to one of the joint 2D coordinate frames.
   */
  virtual std::shared_ptr<frame_2D<double>> getFrame2D(std::size_t i) const {
    return {};
  }

  /**
   * Get the number of joint 3D coordinate frames.
   * \return The number of joint 3D coordinate frames.
   */
  virtual std::size_t getFrames3DCount() const { return 0; }

  /**
   * Get access to one of the joint 3D coordinate frames.
   * \param i The index of the joint 3D coordinate frame.
   * \return Access to one of the joint 3D coordinate frames.
   */
  virtual std::shared_ptr<frame_3D<double>> getFrame3D(std::size_t i) const {
    return {};
  }

  /**
   * Get the number of generalized dependent coordinates.
   * \return The number of generalized dependent coordinates.
   */
  virtual std::size_t getDependentCoordsCount() const { return 0; }

  /**
   * Get access to one of the generalized dependent coordinates.
   * \param i The index of the generalized dependent coordinate.
   * \return Access to one of the generalized dependent coordinates.
   */
  virtual std::shared_ptr<joint_dependent_gen_coord> getDependentCoord(
      std::size_t i) const {
    return {};
  }

  /**
   * Get the number of dependent 2D coordinate frames.
   * \return The number of dependent 2D coordinate frames.
   */
  virtual std::size_t getDependentFrames2DCount() const { return 0; }

  /**
   * Get access to one of the dependent 2D coordinate frames.
   * \param i The index of the dependent 2D coordinate frame.
   * \return Access to one of the dependent 2D coordinate frames.
   */
  virtual std::shared_ptr<joint_dependent_frame_2D> getDependentFrame2D(
      std::size_t i) const {
    return {};
  }

  /**
   * Get the number of dependent 3D coordinate frames.
   * \return The number of dependent 3D coordinate frames.
   */
  virtual std::size_t getDependentFrames3DCount() const { return 0; }

  /**
   * Get access to one of the dependent 3D coordinate frames.
   * \param i The index of the dependent 3D coordinate frame.
   * \return Access to one of the dependent 3D coordinate frames.
   */
  virtual std::shared_ptr<joint_dependent_frame_3D> getDependentFrame3D(
      std::size_t i) const {
    return {};
  }

  /**
   * This function performs the direct-kinematics computation. In other words, prior to
   * calling this function, the joint coordinates have been set to a value.
   * Then, this function is called and will fill the dependent coordinates with the
   * values obtained by the motion on the joint coordinates (e.g., end-effector).
   */
  virtual void doDirectMotion() {}

  /**
   * Get the Jacobian matrix for the system (or twist-shaping matrix). The Jacobian takes the velocity
   * information of the system coordinates and frames, and maps them to velocity information
   * of the system's dependent coordinates and frames.
   * \param Jac stores, as output, the calculated system's Jacobian matrix.
   */
  virtual void getJacobianMatrix(
      mat<double, mat_structure::rectangular>& Jac) const {}

  /**
   * Get the Jacobian matrix for the system (or twist-shaping matrix), and its time-derivative.
   * The Jacobian takes the velocity information of the system coordinates and frames, and maps
   * them to velocity information of the system's dependent coordinates and frames. The time-derivative
   * of the Jacobian matrix will map the velocity information of the system coordinates and frames
   * to the acceleration information of the system's dependent coordinates and frames.
   * \param Jac stores, as output, the calculated system's Jacobian matrix.
   * \param JacDot stores, as output, the calculated time-derivative of the system's Jacobian matrix.
   */
  virtual void getJacobianMatrixAndDerivative(
      mat<double, mat_structure::rectangular>& Jac,
      mat<double, mat_structure::rectangular>& JacDot) const {}

  /**
   * Obtain a vector that contains all the joint position lower-bounds concatenated into one vector.
   * \return All the joint position lower-bounds concatenated into one vector.
   */
  virtual vect_n<double> getJointPositionLowerBounds() const { return {}; }

  /**
   * Sets all the joint position lower-bounds concatenated into one vector.
   * \param aJointLowerBounds All the joint position lower-bounds concatenated into one vector.
   */
  virtual void setJointPositionLowerBounds(
      const vect_n<double>& aJointLowerBounds) {}

  /**
   * Obtain a vector that contains all the joint position upper-bounds concatenated into one vector.
   * \return All the joint position upper-bounds concatenated into one vector.
   */
  virtual vect_n<double> getJointPositionUpperBounds() const { return {}; }

  /**
   * Sets all the joint position lower-bounds concatenated into one vector.
   * \param aJointUpperBounds All the joint position lower-bounds concatenated into one vector.
   */
  virtual void setJointPositionUpperBounds(
      const vect_n<double>& aJointUpperBounds) {}

  /**
   * Obtain a vector that contains all the joint positions concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint positions of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the position (and orientation) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the position (and orientation) of the 3D frame joints.
   * \return All the joint positions concatenated into one vector.
   */
  virtual vect_n<double> getJointPositions() const { return {}; }

  /**
   * Set all the joint positions of the manipulator to a vector of concatenated joint-positions.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint positions of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the position (and orientation) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the position (and orientation) of the 3D frame joints.
   * \param aJointPositions All the joint positions concatenated into one vector.
   */
  virtual void setJointPositions(const vect_n<double>& aJointPositions) {}

  /**
   * Obtain a vector that contains all the joint velocities concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint velocities of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the velocity (and ang. vel.) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the velocity (and ang. vel.) of the 3D frame joints.
   * \note The vector obtained by this function can be multiplied by the Jacobian matrix to obtain
   *       the velocities of the output frames. Similarly, it can be multiplied by the time-derivative
   *       of the Jacobian matrix to obtain the joint-velocity contribution to the acceleration on the output frames.
   * \return All the joint velocities concatenated into one vector.
   */
  virtual vect_n<double> getJointVelocities() const { return {}; }

  /**
   * Set all the joint velocities of the manipulator to a vector of concatenated joint-velocities.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint velocities of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the velocity (and ang. vel.) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the velocity (and ang. vel.) of the 3D frame joints.
   * \param aJointVelocities All the joint velocities concatenated into one vector.
   */
  virtual void setJointVelocities(const vect_n<double>& aJointVelocities) {}

  /**
   * Obtain a vector that contains all the joint accelerations concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint accelerations of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the acceleration (and ang. acc.) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the acceleration (and ang. acc.) of the 3D frame joints.
   * \note The vector obtained by this function can be multiplied by the Jacobian matrix to obtain
   *       the joint-acceleration contribution to the acceleration on the output frames.
   * \return All the joint accelerations concatenated into one vector.
   */
  virtual vect_n<double> getJointAccelerations() const { return {}; }

  /**
   * Set all the joint accelerations of the manipulator to a vector of concatenated joint-accelerations.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the joints are sorted in the same order as in the container returned by getCoord(),
   * getFrame2D() and getFrame3D(), respectively. In other words, the first getCoordsCount()
   * elements are the joint accelerations of generalized coordinate joints, the next 4 * getFrames2DCount()
   * elements are the acceleration (and ang. acc.) of the 2D frame joints, and the
   * final 7 * getFrames3DCount() are the acceleration (and ang. acc.) of the 3D frame joints.
   * \param aJointAccelerations All the joint accelerations concatenated into one vector.
   */
  virtual void setJointAccelerations(
      const vect_n<double>& aJointAccelerations) {}

  /**
   * Obtain a vector that contains all the dependent positions concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependent frames are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In
   * other words, the first getDependentCoordsCount() elements are the positions of
   * dependent generalized coordinates, the next 4 * getDependentFrames2DCount() elements are
   * the position (and orientation) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount()
   * are the position (and orientation) of the 3D dependent frames.
   * \return All the dependent positions concatenated into one vector.
   */
  virtual vect_n<double> getDependentPositions() const { return {}; }

  /**
   * Obtain a vector that contains all the dependent velocities concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependent frames are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In
   * other words, the first getDependentCoordsCount() elements are the velocities of
   * dependent generalized coordinates, the next 4 * getDependentFrames2DCount() elements are
   * the velocity (and ang. vel.) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount()
   * are the velocity (and ang. vel.) of the 3D dependent frames.
   * \return All the dependent velocities concatenated into one vector.
   */
  virtual vect_n<double> getDependentVelocities() const { return {}; }

  /**
   * Obtain a vector that contains all the dependent accelerations concatenated into one vector.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependent frames are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In
   * other words, the first getDependentCoordsCount() elements are the accelerations of
   * dependent generalized coordinates, the next 4 * getDependentFrames2DCount() elements are
   * the acceleration (and ang. acc.) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount()
   * are the acceleration (and ang. acc.) of the 3D dependent frames.
   * \return All the dependent accelerations concatenated into one vector.
   */
  virtual vect_n<double> getDependentAccelerations() const { return {}; }

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(direct_kinematics_model, 0xC210004F, 1,
                              "direct_kinematics_model", named_object)
};
}  // namespace ReaK::kte

#endif
