/**
 * \file inverse_kinematics_model.hpp
 *
 * This library declares a base-class to represent classes that can be used to compute the
 * inverse kinematics of some (kte-based) model.
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

#ifndef REAK_INVERSE_KINEMATICS_MODEL_HPP
#define REAK_INVERSE_KINEMATICS_MODEL_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/math/kinetostatics/kinetostatics.hpp"
#include "ReaK/mbd/kte/jacobian_joint_map.hpp"
#include "ReaK/mbd/models/direct_kinematics_model.hpp"

namespace ReaK::kte {

/**
 * This base-class represents a class that can be used to compute the inverse kinematics
 * of some (kte-based) model. This class makes a distinction between joint coordinates
 * and dependent coordinates, where the former are generally coordinates (generalized,
 * and 2D/3D frames) determining the state of the joints of a kinematic chain, while the
 * latter are generally coordinates (generalized, and 2D/3D frames) whose state are dependent
 * on the joint states of a kinematic chain.
 */
class inverse_kinematics_model : public direct_kinematics_model {
 public:
  /**
   * Default constructor.
   */
  explicit inverse_kinematics_model(const std::string& aName = "")
      : direct_kinematics_model(aName) {}

  /**
   * Default destructor.
   */
  ~inverse_kinematics_model() override = default;

  /**
   * This function performs the inverse-kinematics computation. In other words, prior to
   * calling this function, the dependent coordinates have been set to a desired value.
   * Then, this function is called and will fill the joint coordinates with the necessary
   * values to obtained the desired motion on the dependent coordinates (e.g., end-effector).
   */
  virtual void doInverseMotion() {}

  /**
   * Set all the dependent positions of the manipulator to a vector of concatenated dependent positions.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependents are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In other
   * words, the first getDependentCoordsCount() elements are the joint positions of generalized
   * coordinate joints, the next 4 * getDependentFrames2DCount() elements are the position (and
   * orientation) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount() are
   * the position (and orientation) of the 3D dependent frames.
   * \param aDepPositions All the dependent positions concatenated into one vector.
   */
  virtual void setDependentPositions(const vect_n<double>& aDepPositions) {}

  /**
   * Set all the dependent velocities of the manipulator to a vector of concatenated dependent velocities.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependents are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In other
   * words, the first getDependentCoordsCount() elements are the dependent velocities of generalized
   * coordinates, the next 4 * getDependentFrames2DCount() elements are the velocity (and
   * ang. vel.) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount() are
   * the velocity (and ang. vel.) of the 3D dependent frames.
   * \param aDepVelocities All the dependent velocities concatenated into one vector.
   */
  virtual void setDependentVelocities(const vect_n<double>& aDepVelocities) {}

  /**
   * Set all the dependent accelerations of the manipulator to a vector of concatenated dependent accelerations.
   * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ),
   * where the dependents are sorted in the same order as in the container returned by
   * getDependentCoord(), getDependentFrame2D() and getDependentFrame3D(), respectively. In other
   * words, the first getDependentCoordsCount() elements are the dependent accelerations of generalized
   * coordinates, the next 4 * getDependentFrames2DCount() elements are the acceleration (and
   * ang. acc.) of the 2D dependent frames, and the final 7 * getDependentFrames3DCount() are
   * the acceleration (and ang. acc.) of the 3D dependent frames.
   * \param aDepAccelerations All the dependent accelerations concatenated into one vector.
   */
  virtual void setDependentAccelerations(
      const vect_n<double>& aDepAccelerations) {}

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    direct_kinematics_model::save(
        A, direct_kinematics_model::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    direct_kinematics_model::load(
        A, direct_kinematics_model::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(inverse_kinematics_model, 0xC2100050, 1,
                              "inverse_kinematics_model",
                              direct_kinematics_model)
};

}  // namespace ReaK::kte

#endif
