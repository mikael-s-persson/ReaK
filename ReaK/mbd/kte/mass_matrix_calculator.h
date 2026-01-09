/**
 * \file mass_matrix_calculator.h
 *
 * This library declares classes to compute the system's mass matrix in a single pass.
 * The mass matrix calculators use the jacobian mappings of up-stream joints, for each
 * inertial element to determine the twist-shaping matrices and the block-diagonal, link mass matrix.
 * The similarity transform of the latter by the former gives the system's mass matrix, and
 * similarily, the time-derivative of the mass matrix.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_MBD_KTE_MASS_MATRIX_CALCULATOR_H_
#define REAK_MBD_KTE_MASS_MATRIX_CALCULATOR_H_

#include "ReaK/mbd/kte/inertia.h"

#include <vector>

namespace ReaK::kte {

/**
 * This class is a mass matrix calculator for a system. It holds lists of all
 * inertial elements (which contain up-stream joint Jacobians) as well as a list of coordinates
 * which it uses to compute the mass matrix and its time-derivative using the twist-shaping
 * matrix formulation.
 */
class mass_matrix_calc : public named_object {
 private:
  std::vector<std::shared_ptr<inertia_gen>>
      mGenInertias;  ///< Holds the list of generalized coordinate inertial elements.
  std::vector<std::shared_ptr<inertia_2D>>
      m2DInertias;  ///< Holds the list of 2D inertial elements.
  std::vector<std::shared_ptr<inertia_3D>>
      m3DInertias;  ///< Holds the list of 3D inertial elements.

  std::vector<std::shared_ptr<gen_coord<double>>>
      mCoords;  ///< Holds the list of generalized coordinates in the system.
  std::vector<std::shared_ptr<frame_2D<double>>>
      mFrames2D;  ///< Holds the list of 2D coordinates frame in the system.
  std::vector<std::shared_ptr<frame_3D<double>>>
      mFrames3D;  ///< Holds the list of 3D coordinates frame in the system.

 public:
  /**
   * Default constructor.
   */
  explicit mass_matrix_calc(const std::string& aName = "") {
    this->set_name(aName);
  }

  /**
   * Add a generalized coordinate inertial element to the mass matrix calculation.
   * \param aGenInertia a generalized coordinate inertial element to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(const std::shared_ptr<inertia_gen>& aGenInertia);

  /**
   * Add a 2D inertial element to the mass matrix calculation.
   * \param a2DInertia a 2D inertial element to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(const std::shared_ptr<inertia_2D>& a2DInertia);

  /**
   * Add a 3D inertial element to the mass matrix calculation.
   * \param a3DInertia a 3D inertial element to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(const std::shared_ptr<inertia_3D>& a3DInertia);

  /**
   * Add a system generalized coordinate to the mass matrix calculation.
   * \param aCoord a system generalized coordinate to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(
      const std::shared_ptr<gen_coord<double>>& aCoord);

  /**
   * Add a system generalized coordinate to the mass matrix calculation.
   * \param aFrame2D a system 2D frame to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(
      const std::shared_ptr<frame_2D<double>>& aFrame2D);

  /**
   * Add a system generalized coordinate to the mass matrix calculation.
   * \param aFrame3D a system 3D frame to add.
   * \return reference to this.
   */
  mass_matrix_calc& operator<<(
      const std::shared_ptr<frame_3D<double>>& aFrame3D);

  /** Get read-only access to the list of generalized coordinates. */
  const std::vector<std::shared_ptr<gen_coord<double>>>& Coords() const {
    return mCoords;
  }

  /** Get read-only access to the list of 2D coordinate frames. */
  const std::vector<std::shared_ptr<frame_2D<double>>>& Frames2D() const {
    return mFrames2D;
  }

  /** Get read-only access to the list of 3D coordinate frames. */
  const std::vector<std::shared_ptr<frame_3D<double>>>& Frames3D() const {
    return mFrames3D;
  }

  /**
   * Get the mass matrix for the system.
   * \param M stores, as output, the calculated system's mass-matrix.
   */
  void getMassMatrix(mat<double, mat_structure::symmetric>& M);

  /**
   * Get the mass matrix for the system and its time-derivative.
   * \param M stores, as output, the calculated system's mass-matrix.
   * \param M_dot stores, as output, the calculated time-derivative of the system's mass matrix.
   */
  void getMassMatrixAndDerivative(mat<double, mat_structure::symmetric>& M,
                                  mat<double, mat_structure::square>& M_dot);

  /**
   * Get the twist-shaping matrix, the block-diagonal, link mass-matrix, and the time-derivative of the twist-shaping
   * matrix.
   * \param Tcm stores, as output, the calculated system's twist-shaping matrix.
   * \param Mcm stores, as output, the calculated block-diagonal, link mass matrix.
   * \param Tcm_dot storse, as output, the calculated time-derivative of the system's twist-shaping matrix.
   */
  void get_TMT_TdMT(mat<double, mat_structure::rectangular>& Tcm,
                    mat<double, mat_structure::symmetric>& Mcm,
                    mat<double, mat_structure::rectangular>& Tcm_dot);

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mGenInertias) &
        RK_SERIAL_SAVE_WITH_NAME(m2DInertias) &
        RK_SERIAL_SAVE_WITH_NAME(m3DInertias) &
        RK_SERIAL_SAVE_WITH_NAME(mCoords) &
        RK_SERIAL_SAVE_WITH_NAME(mFrames2D) &
        RK_SERIAL_SAVE_WITH_NAME(mFrames3D);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    named_object::load(A, named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mGenInertias) &
        RK_SERIAL_LOAD_WITH_NAME(m2DInertias) &
        RK_SERIAL_LOAD_WITH_NAME(m3DInertias) &
        RK_SERIAL_LOAD_WITH_NAME(mCoords) &
        RK_SERIAL_LOAD_WITH_NAME(mFrames2D) &
        RK_SERIAL_LOAD_WITH_NAME(mFrames3D);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(mass_matrix_calc, 0xC2000001, 1,
                              "mass_matrix_calc", named_object)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_MASS_MATRIX_CALCULATOR_H_
