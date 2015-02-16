/**
 * \file manip_kinematics_model.hpp
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

#ifndef REAK_MANIP_KINEMATICS_MODEL_HPP
#define REAK_MANIP_KINEMATICS_MODEL_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <ReaK/mbd/kte/kte_map_chain.hpp>
#include "direct_kinematics_model.hpp"

#include <vector>

namespace ReaK {

namespace kte {


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
    std::vector< shared_ptr< gen_coord<double> > > mCoords; ///< Holds the list of generalized coordinates in the system.
    std::vector< shared_ptr< frame_2D<double> > > mFrames2D; ///< Holds the list of 2D coordinates frame in the system.
    std::vector< shared_ptr< frame_3D<double> > > mFrames3D; ///< Holds the list of 3D coordinates frame in the system.

    std::vector< shared_ptr< joint_dependent_gen_coord > > mDependentGenCoords; ///< Holds the list of dependent generalized coordinates.
    std::vector< shared_ptr< joint_dependent_frame_2D > > mDependent2DFrames; ///< Holds the list of dependent 2D frames.
    std::vector< shared_ptr< joint_dependent_frame_3D > > mDependent3DFrames; ///< Holds the list of dependent 3D frames.

    shared_ptr< kte_map_chain > mModel; ///< Holds the model of the manipulator as a kte-chain.
    
  public:
    
    /**
     * Default constructor.
     */
    manipulator_kinematics_model(const std::string& aName = "") : direct_kinematics_model(aName),
                                                                  mCoords(),
                                                                  mFrames2D(),
                                                                  mFrames3D(),
                                                                  mDependentGenCoords(),
                                                                  mDependent2DFrames(),
                                                                  mDependent3DFrames(),
                                                                  mModel() { };
    
    /**
     * Default destructor.
     */
    virtual ~manipulator_kinematics_model() { };
    
    /**
     * Sets the manipulator KTE model to use in this object.
     * \param aModel The manipulator KTE model to use in this object.
     */
    virtual void setModel(const shared_ptr< kte_map_chain >& aModel) { mModel = aModel; };
    
    /**
     * Gets the manipulator KTE model used by this object.
     * \return The manipulator KTE model used by this object.
     */
    shared_ptr< kte_map_chain > getModel() const { return mModel; };

    /**
     * Add a dependent generalized coordinate to the jacobian calculation.
     * \param aDependentGenCoord a dependent generalized coordinate to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< joint_dependent_gen_coord >& aDependentGenCoord);

    /**
     * Add a dependent 2D frame to the jacobian calculation.
     * \param aDependent2DFrame a dependent 2D frame  to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< joint_dependent_frame_2D >& aDependent2DFrame);

    /**
     * Add a dependent 3D frame to the jacobian calculation.
     * \param aDependent3DFrame a dependent 3D frame to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< joint_dependent_frame_3D >& aDependent3DFrame);

    /**
     * Add a system generalized coordinate.
     * \param aCoord a system generalized coordinate to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< gen_coord<double> >& aCoord);

    /**
     * Add a system 2D frame.
     * \param aFrame2D a system 2D frame to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< frame_2D<double> >& aFrame2D);

    /**
     * Add a system 3D frame.
     * \param aFrame3D a system 3D frame to add.
     * \return reference to this.
     */
    virtual manipulator_kinematics_model& operator <<(const shared_ptr< frame_3D<double> >& aFrame3D);
    
    
    /******************************************************************************************
     *  direct_kinematics_model: joint and dependent positions, velocities and accelerations  *
     ******************************************************************************************/
    
    virtual std::size_t getJointPositionsCount() const {
      return mCoords.size() + 4 * mFrames2D.size() + 7 * mFrames3D.size();
    };
    
    virtual std::size_t getJointVelocitiesCount() const {
      return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
    };
    
    virtual std::size_t getJointAccelerationsCount() const {
      return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
    };
    
    virtual std::size_t getDependentPositionsCount() const {
      return mDependentGenCoords.size() + 4 * mDependent2DFrames.size() + 7 * mDependent3DFrames.size();
    };
    
    virtual std::size_t getDependentVelocitiesCount() const {
      return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() + 6 * mDependent3DFrames.size();
    };
    
    virtual std::size_t getDependentAccelerationsCount() const {
      return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() + 6 * mDependent3DFrames.size();
    };
    
    
    
    /*************************************************************************
     *  direct_kinematics_model: joint and dependent coordinates and frames  *
     *************************************************************************/
    
    virtual std::size_t getCoordsCount() const { return mCoords.size(); };
    
    virtual shared_ptr< gen_coord<double> > getCoord(std::size_t i) const { 
      return mCoords[i];
    };
    
    virtual std::size_t getFrames2DCount() const { return mFrames2D.size(); };
    
    virtual shared_ptr< frame_2D<double> > getFrame2D(std::size_t i) const { 
      return mFrames2D[i];
    };
    
    virtual std::size_t getFrames3DCount() const { return mFrames3D.size(); };
    
    virtual shared_ptr< frame_3D<double> > getFrame3D(std::size_t i) const { 
      return mFrames3D[i];
    };
    
    virtual std::size_t getDependentCoordsCount() const { return mDependentGenCoords.size(); };
    
    virtual shared_ptr< joint_dependent_gen_coord > getDependentCoord(std::size_t i) const { 
      return mDependentGenCoords[i];
    };
    
    virtual std::size_t getDependentFrames2DCount() const { return mDependent2DFrames.size(); };
    
    virtual shared_ptr< joint_dependent_frame_2D > getDependentFrame2D(std::size_t i) const { 
      return mDependent2DFrames[i];
    };
    
    virtual std::size_t getDependentFrames3DCount() const { return mDependent3DFrames.size(); };
    
    virtual shared_ptr< joint_dependent_frame_3D > getDependentFrame3D(std::size_t i) const { 
      return mDependent3DFrames[i];
    };
    
    
    /************************************************************************
     *  direct_kinematics_model: kinematic functions (motion and Jacobian)  *
     ************************************************************************/
    
    virtual void doDirectMotion() {
      if(mModel)
        mModel->doMotion();
    };
    
    virtual void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const;
    
    virtual void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const;
    
    
    /***********************************************************************************
     *  direct_kinematics_model: contatenated positions, velocities and accelerations  *
     ***********************************************************************************/
    
    virtual vect_n<double> getJointPositions() const;
    
    virtual void setJointPositions(const vect_n<double>& aJointPositions);
    
    virtual vect_n<double> getJointVelocities() const;
    
    virtual void setJointVelocities(const vect_n<double>& aJointVelocities);
    
    virtual vect_n<double> getJointAccelerations() const;
    
    virtual void setJointAccelerations(const vect_n<double>& aJointAccelerations);
    
    virtual vect_n<double> getDependentPositions() const;
    
    virtual vect_n<double> getDependentVelocities() const;
    
    virtual vect_n<double> getDependentAccelerations() const;
    
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      direct_kinematics_model::save(A,direct_kinematics_model::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCoords)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames2D)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames3D)
        & RK_SERIAL_SAVE_WITH_NAME(mDependentGenCoords)
        & RK_SERIAL_SAVE_WITH_NAME(mDependent2DFrames)
        & RK_SERIAL_SAVE_WITH_NAME(mDependent3DFrames)
        & RK_SERIAL_SAVE_WITH_NAME(mModel);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      direct_kinematics_model::load(A,direct_kinematics_model::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCoords)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames2D)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames3D)
        & RK_SERIAL_LOAD_WITH_NAME(mDependentGenCoords)
        & RK_SERIAL_LOAD_WITH_NAME(mDependent2DFrames)
        & RK_SERIAL_LOAD_WITH_NAME(mDependent3DFrames)
        & RK_SERIAL_LOAD_WITH_NAME(mModel);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(manipulator_kinematics_model,0xC210004D,1,"manipulator_kinematics_model",direct_kinematics_model)

};




};

};

#endif











