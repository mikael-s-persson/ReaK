/**
 * \file inverse_dynamics_model.h
 *
 * This library declares a base-class to represent classes that can be used to compute the
 * inverse dynamics of some (kte-based) model.
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

#ifndef REAK_MBD_MODELS_INVERSE_DYNAMICS_MODEL_H_
#define REAK_MBD_MODELS_INVERSE_DYNAMICS_MODEL_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * This base-class represents a class that can be used to compute the inverse dynamics
 * of some (kte-based) model. This class makes a distinction between joint coordinates
 * and dependent coordinates, where the former are generally coordinates (generalized,
 * and 2D/3D frames) determining the state of the joints of a kinematic chain, while the
 * latter are generally coordinates (generalized, and 2D/3D frames) whose state are dependent
 * on the joint states of a kinematic chain.
 */
class inverse_dynamics_model : public kte_map {
 public:
  /**
   * Default constructor.
   */
  explicit inverse_dynamics_model(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Default destructor.
   */
  ~inverse_dynamics_model() override = default;
  ;

  /**
   * Get the total number of position values for all the joint frames concatenated.
   * \return The total number of position values for all the joint frames concatenated.
   */
  virtual std::size_t getJointStatesCount() const { return 0; }

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
  virtual vect_n<double> getJointStates() const { return {}; }

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
  virtual void setJointStates(const vect_n<double>& aJointStates) {}

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {}

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  /**
   * Get the mass matrix for the system.
   * \param M stores, as output, the calculated system's mass-matrix.
   */
  virtual void getMassMatrix(mat<double, mat_structure::symmetric>& M) {}

  /**
   * Get the mass matrix for the system and its time-derivative.
   * \param M stores, as output, the calculated system's mass-matrix.
   * \param M_dot stores, as output, the calculated time-derivative of the system's mass matrix.
   */
  virtual void getMassMatrixAndDerivative(
      mat<double, mat_structure::symmetric>& M,
      mat<double, mat_structure::square>& M_dot) {}

  /**
   * Get the twist-shaping matrix, the block-diagonal, link mass-matrix, and the time-derivative of the twist-shaping
   * matrix.
   * \param Tcm stores, as output, the calculated system's twist-shaping matrix.
   * \param Mcm stores, as output, the calculated block-diagonal, link mass matrix.
   * \param Tcm_dot storse, as output, the calculated time-derivative of the system's twist-shaping matrix.
   */
  virtual void get_TMT_TdMT(mat<double, mat_structure::rectangular>& Tcm,
                            mat<double, mat_structure::symmetric>& Mcm,
                            mat<double, mat_structure::rectangular>& Tcm_dot) {}

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(inverse_dynamics_model, 0xC210004F, 1,
                              "inverse_dynamics_model", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_MODELS_INVERSE_DYNAMICS_MODEL_H_
