/**
 * \file manipulator_model.hpp
 *
 * This library declares classes to represent manipulator systems, both kinematic only or 
 * dynamic as well. Essentially, the model of the manipulator is only a KTE chain provided 
 * by the user, but these manipulator-model classes take care of grouping the joints, their
 * limits, and their jacobian matrices.
 * The jacobian matrices are computed from jacobian mappings of up-stream joints, for each
 * output-frame to determine the twist-shaping matrices.
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

#ifndef REAK_MANIPULATOR_MODEL_HPP
#define REAK_MANIPULATOR_MODEL_HPP

#include "kinetostatics/kinetostatics.hpp"
#include "kte_map_chain.hpp"
#include "mass_matrix_calculator.hpp"

#include <vector>

namespace ReaK {

namespace kte {


/**
 * This class is a 
 */
class manipulator_kinematics_model : public named_object {
  private:
    std::vector< shared_pointer< gen_coord<double> >::type > mCoords; ///< Holds the list of generalized coordinates in the system.
    std::vector< shared_pointer< frame_2D<double> >::type > mFrames2D; ///< Holds the list of 2D coordinates frame in the system.
    std::vector< shared_pointer< frame_3D<double> >::type > mFrames3D; ///< Holds the list of 3D coordinates frame in the system.

    std::vector< shared_pointer< inertia_gen >::type > mGenInertias; ///< Holds the list of generalized coordinate inertial elements.
    std::vector< shared_pointer< inertia_2D >::type > m2DInertias; ///< Holds the list of 2D inertial elements.
    std::vector< shared_pointer< inertia_3D >::type > m3DInertias; ///< Holds the list of 3D inertial elements.

    
  public:

    /**
     * Default constructor.
     */
    mass_matrix_calc(const std::string& aName = "") : named_object(),
                                                      mGenInertias(),
                                                      m2DInertias(),
                                                      m3DInertias(),
                                                      mCoords(),
                                                      mFrames2D(),
                                                      mFrames3D() {
      this->setName(aName);
    };

    /**
     * Add a generalized coordinate inertial element to the mass matrix calculation.
     * \param aGenInertia a generalized coordinate inertial element to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< inertia_gen >::type& aGenInertia);

    /**
     * Add a 2D inertial element to the mass matrix calculation.
     * \param a2DInertia a 2D inertial element to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< inertia_2D >::type& a2DInertia);

    /**
     * Add a 3D inertial element to the mass matrix calculation.
     * \param a3DInertia a 3D inertial element to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< inertia_3D >::type& a3DInertia);

    /**
     * Add a system generalized coordinate to the mass matrix calculation.
     * \param aCoord a system generalized coordinate to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< gen_coord<double> >::type& aCoord);

    /**
     * Add a system generalized coordinate to the mass matrix calculation.
     * \param aFrame2D a system 2D frame to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< frame_2D<double> >::type& aFrame2D);

    /**
     * Add a system generalized coordinate to the mass matrix calculation.
     * \param aFrame3D a system 3D frame to add.
     * \return reference to this.
     */
    mass_matrix_calc& operator <<(const shared_pointer< frame_3D<double> >::type& aFrame3D);

    /** Get read-only access to the list of generalized coordinates. */
    const std::vector< shared_pointer< gen_coord<double> >::type >& Coords() const { return mCoords; };

    /** Get read-only access to the list of 2D coordinate frames. */
    const std::vector< shared_pointer< frame_2D<double> >::type >& Frames2D() const { return mFrames2D; };

    /** Get read-only access to the list of 3D coordinate frames. */
    const std::vector< shared_pointer< frame_3D<double> >::type >& Frames3D() const { return mFrames3D; };

    /**
     * Get the mass matrix for the system.
     * \param M stores, as output, the calculated system's mass-matrix.
     */
    void getMassMatrix(mat<double,mat_structure::symmetric>& M);

    /**
     * Get the mass matrix for the system and its time-derivative.
     * \param M stores, as output, the calculated system's mass-matrix.
     * \param M_dot stores, as output, the calculated time-derivative of the system's mass matrix.
     */
    void getMassMatrixAndDerivative(mat<double,mat_structure::symmetric>& M, mat<double,mat_structure::square>& M_dot);

    /**
     * Get the twist-shaping matrix, the block-diagonal, link mass-matrix, and the time-derivative of the twist-shaping matrix.
     * \param Tcm stores, as output, the calculated system's twist-shaping matrix.
     * \param Mcm stores, as output, the calculated block-diagonal, link mass matrix.
     * \param Tcm_dot storse, as output, the calculated time-derivative of the system's twist-shaping matrix.
     */
    void get_TMT_TdMT(mat<double,mat_structure::rectangular>& Tcm, mat<double,mat_structure::symmetric>& Mcm, mat<double,mat_structure::rectangular>& Tcm_dot);

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mGenInertias)
        & RK_SERIAL_SAVE_WITH_NAME(m2DInertias)
        & RK_SERIAL_SAVE_WITH_NAME(m3DInertias)
        & RK_SERIAL_SAVE_WITH_NAME(mCoords)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames2D)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames3D);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mGenInertias)
        & RK_SERIAL_LOAD_WITH_NAME(m2DInertias)
        & RK_SERIAL_LOAD_WITH_NAME(m3DInertias)
        & RK_SERIAL_LOAD_WITH_NAME(mCoords)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames2D)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames3D);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(mass_matrix_calc,0xC2000001,1,"mass_matrix_calc",named_object)

};

};

};

#endif











