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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/kinetostatics/kinetostatics.hpp>
#include "kte_map_chain.hpp"
#include "mass_matrix_calculator.hpp"
#include "kte_system_input.hpp"
#include "kte_system_output.hpp"

#include <ReaK/core/integrators/integrator.hpp>

#include <vector>

namespace ReaK {

namespace kte {


class manip_kin_mdl_joint_io; //forward-declaration.
class manip_kin_mdl_jac_calculator; // forward-declaration.
  

/**
 * This class stores the required information to represent the kinematic model of a manipulator.
 * Here, a manipulator is defined as a kinematic chain with "input" coordinates (or frames) and 
 * "output" coordinates (or frames). For example, a typical serial manipulator 
 * could have a set of generalized coordinates (joint coordinates) as well as one or more frames 
 * for the end-effector(s) (or additional link motions). This class is basically used to 
 * regroup all that information and provides a certain number of functions related to the 
 * use of a manipulator model (like computing jacobians).
 */
class manipulator_kinematics_model : public kte_map {
  protected:
    std::vector< shared_ptr< gen_coord<double> > > mCoords; ///< Holds the list of generalized coordinates in the system.
    std::vector< shared_ptr< frame_2D<double> > > mFrames2D; ///< Holds the list of 2D coordinates frame in the system.
    std::vector< shared_ptr< frame_3D<double> > > mFrames3D; ///< Holds the list of 3D coordinates frame in the system.

    std::vector< shared_ptr< joint_dependent_gen_coord > > mDependentGenCoords; ///< Holds the list of dependent generalized coordinates.
    std::vector< shared_ptr< joint_dependent_frame_2D > > mDependent2DFrames; ///< Holds the list of dependent 2D frames.
    std::vector< shared_ptr< joint_dependent_frame_3D > > mDependent3DFrames; ///< Holds the list of dependent 3D frames.

    shared_ptr< kte_map_chain > mModel; ///< Holds the model of the manipulator as a kte-chain.
    
    friend class manip_kin_mdl_joint_io;
    friend class manip_kin_mdl_jac_calculator;
    
  public:
    
    /**
     * Default constructor.
     */
    manipulator_kinematics_model(const std::string& aName = "") : kte_map(aName),
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

    /**
     * Get the total number of position values for all the joint frames concatenated.
     * \return The total number of position values for all the joint frames concatenated.
     */
    unsigned int getJointPositionsCount() const {
      return mCoords.size() + 4 * mFrames2D.size() + 7 * mFrames3D.size();
    };
    
    /**
     * Get the total number of velocity values for all the joint frames concatenated.
     * \return The total number of velocity values for all the joint frames concatenated.
     */
    unsigned int getJointVelocitiesCount() const {
      return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
    };
    
    /**
     * Get the total number of acceleration values for all the joint frames concatenated.
     * \return The total number of acceleration values for all the joint frames concatenated.
     */
    unsigned int getJointAccelerationsCount() const {
      return mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
    };
    
    /**
     * Get the total number of position values for all the dependent frames concatenated.
     * \return The total number of position values for all the dependent frames concatenated.
     */
    unsigned int getDependentPositionsCount() const {
      return mDependentGenCoords.size() + 4 * mDependent2DFrames.size() + 7 * mDependent3DFrames.size();
    };
    
    /**
     * Get the total number of velocity values for all the dependent frames concatenated.
     * \return The total number of velocity values for all the dependent frames concatenated.
     */
    unsigned int getDependentVelocitiesCount() const {
      return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() + 6 * mDependent3DFrames.size();
    };
    
    /**
     * Get the total number of acceleration values for all the dependent frames concatenated.
     * \return The total number of acceleration values for all the dependent frames concatenated.
     */
    unsigned int getDependentAccelerationsCount() const {
      return mDependentGenCoords.size() + 3 * mDependent2DFrames.size() + 6 * mDependent3DFrames.size();
    };
    
    /** Get read-only access to the list of generalized coordinates. */
    const std::vector< shared_ptr< gen_coord<double> > >& Coords() const { return mCoords; };

    /** Get read-only access to the list of 2D coordinate frames. */
    const std::vector< shared_ptr< frame_2D<double> > >& Frames2D() const { return mFrames2D; };

    /** Get read-only access to the list of 3D coordinate frames. */
    const std::vector< shared_ptr< frame_3D<double> > >& Frames3D() const { return mFrames3D; };

    /** Get read-only access to the list of generalized coordinates. */
    const std::vector< shared_ptr< joint_dependent_gen_coord > >& DependentCoords() const { return mDependentGenCoords; };

    /** Get read-only access to the list of 2D coordinate frames. */
    const std::vector< shared_ptr< joint_dependent_frame_2D > >& DependentFrames2D() const { return mDependent2DFrames; };

    /** Get read-only access to the list of 3D coordinate frames. */
    const std::vector< shared_ptr< joint_dependent_frame_3D > >& DependentFrames3D() const { return mDependent3DFrames; };

    /**
     * Get the Jacobian matrix for the system (or twist-shaping matrix). The Jacobian takes the velocity 
     * information of the system coordinates and frames, and maps them to velocity information 
     * of the system's dependent coordinates and frames.
     * \param Jac stores, as output, the calculated system's Jacobian matrix.
     */
    virtual void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const;

    /**
     * Get the Jacobian matrix for the system (or twist-shaping matrix), and its time-derivative. 
     * The Jacobian takes the velocity information of the system coordinates and frames, and maps 
     * them to velocity information of the system's dependent coordinates and frames. The time-derivative
     * of the Jacobian matrix will map the velocity information of the system coordinates and frames 
     * to the acceleration information of the system's dependent coordinates and frames.
     * \param Jac stores, as output, the calculated system's Jacobian matrix.
     * \param JacDot stores, as output, the calculated time-derivative of the system's Jacobian matrix.
     */
    virtual void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const;

    
    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>()) {
      if(mModel)
        mModel->doMotion(aFlag,aStorage);
    };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>()) {
      if(mModel)
        mModel->doForce(aFlag,aStorage);
    };

    virtual void clearForce() { 
      if(mModel)
        mModel->clearForce();
    };
    
    /**
     * Obtain a vector that contains all the joint positions concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint positions of generalized coordinate joints, the next 4 * Frames2D().size() 
     * elements are the position (and orientation) of the 2D frame joints, and the final 7 * Frames3D().size()
     * are the position (and orientation) of the 3D frame joints.
     * \return All the joint positions concatenated into one vector.
     */
    vect_n<double> getJointPositions() const;
    
    /**
     * Set all the joint positions of the manipulator to a vector of concatenated joint-positions.
     * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint positions of generalized coordinate joints, the next 4 * Frames2D().size() 
     * elements are the position (and orientation) of the 2D frame joints, and the final 7 * Frames3D().size()
     * are the position (and orientation) of the 3D frame joints.
     * \param aJointPositions All the joint positions concatenated into one vector.
     */
    void setJointPositions(const vect_n<double>& aJointPositions);
    
    /**
     * Obtain a vector that contains all the joint velocities concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Velocities, 2D Velocities, 3D Velocities ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint velocities of generalized coordinate joints, the next 3 * Frames2D().size() 
     * elements are the velocities (and angular vel.) of the 2D frame joints, and the final 6 * Frames3D().size()
     * are the velocities (and angular vel.) of the 3D frame joints.
     * \note The vector obtained by this function can be multiplied by the Jacobian matrix to obtain
     *       the velocities of the output frames. Similarly, it can be multiplied by the time-derivative
     *       of the Jacobian matrix to obtain the joint-velocity contribution to the acceleration on the output frames.
     * \return All the joint velocities concatenated into one vector.
     */
    vect_n<double> getJointVelocities() const;
    
    /**
     * Set all the joint velocities of the manipulator to a vector of concatenated joint-velocities.
     * The ordering in the vector is as follows: ( Generalized Velocities, 2D Velocities, 3D Velocities ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint velocities of generalized coordinate joints, the next 3 * Frames2D().size() 
     * elements are the velocities (and angular vel.) of the 2D frame joints, and the final 6 * Frames3D().size()
     * are the velocities (and angular vel.) of the 3D frame joints.
     * \param aJointVelocities All the joint velocities concatenated into one vector.
     */
    void setJointVelocities(const vect_n<double>& aJointVelocities);
    
    /**
     * Obtain a vector that contains all the joint accelerations concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Accelerations, 2D Accelerations, 3D Accelerations ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint accelerations of generalized coordinate joints, the next 3 * Frames2D().size() 
     * elements are the accelerations (and angular acc.) of the 2D frame joints, and the final 6 * Frames3D().size()
     * are the accelerations (and angular acc.) of the 3D frame joints.
     * \note The vector obtained by this function can be multiplied by the Jacobian matrix to obtain
     *       the joint-acceleration contribution to the acceleration on the output frames.
     * \return All the joint accelerations concatenated into one vector.
     */
    vect_n<double> getJointAccelerations() const;
    
    /**
     * Set all the joint accelerations of the manipulator to a vector of concatenated joint-accelerations.
     * The ordering in the vector is as follows: ( Generalized Accelerations, 2D Accelerations, 3D Accelerations ), 
     * where the joints are sorted in the same order as in the container returned by Coords(), 
     * Frames2D() and Frames3D(), respectively. In other words, the first Coords().size() elements
     * are the joint accelerations of generalized coordinate joints, the next 3 * Frames2D().size() 
     * elements are the accelerations (and angular acc.) of the 2D frame joints, and the final 6 * Frames3D().size()
     * are the accelerations (and angular acc.) of the 3D frame joints.
     * \param aJointAccelerations All the joint accelerations concatenated into one vector.
     */
    void setJointAccelerations(const vect_n<double>& aJointAccelerations);
    
    
    /**
     * Obtain a vector that contains all the dependent positions concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ), 
     * where the dependent frames are sorted in the same order as in the container returned by DependentCoords(), 
     * DependentFrames2D() and DependentFrames3D(), respectively. In other words, the first DependentCoords().size() elements
     * are the positions of dependent generalized coordinates, the next 4 * DependentFrames2D().size() 
     * elements are the position (and orientation) of the 2D dependent frames, and the final 7 * DependentFrames3D().size()
     * are the position (and orientation) of the 3D dependent frames.
     * \return All the dependent positions concatenated into one vector.
     */
    vect_n<double> getDependentPositions() const;
    
    /**
     * Obtain a vector that contains all the dependent velocities concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Velocities, 2D Velocities, 3D Velocities ), 
     * where the dependent frames are sorted in the same order as in the container returned by DependentCoords(), 
     * DependentFrames2D() and DependentFrames3D(), respectively. In other words, the first DependentCoords().size() elements
     * are the velocities of dependent generalized coordinates, the next 3 * DependentFrames2D().size() 
     * elements are the velocities (and angular vel.) of the 2D dependent frames, and the final 6 * DependentFrames3D().size()
     * are the velocities (and angular vel.) of the 3D dependent frames.
     * \return All the dependent velocities concatenated into one vector.
     */
    vect_n<double> getDependentVelocities() const;
    
    /**
     * Obtain a vector that contains all the dependent accelerations concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Accelerations, 2D Accelerations, 3D Accelerations ), 
     * where the dependent frames are sorted in the same order as in the container returned by DependentCoords(), 
     * DependentFrames2D() and DependentFrames3D(), respectively. In other words, the first DependentCoords().size() elements
     * are the joint accelerations of dependent generalized coordinates, the next 3 * DependentFrames2D().size() 
     * elements are the accelerations (and angular acc.) of the 2D dependent frames, and the final 6 * DependentFrames3D().size()
     * are the accelerations (and angular acc.) of the 3D dependent frames.
     * \return All the dependent accelerations concatenated into one vector.
     */
    vect_n<double> getDependentAccelerations() const;
    
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCoords)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames2D)
        & RK_SERIAL_SAVE_WITH_NAME(mFrames3D)
        & RK_SERIAL_SAVE_WITH_NAME(mDependentGenCoords)
        & RK_SERIAL_SAVE_WITH_NAME(mDependent2DFrames)
        & RK_SERIAL_SAVE_WITH_NAME(mDependent3DFrames)
        & RK_SERIAL_SAVE_WITH_NAME(mModel);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCoords)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames2D)
        & RK_SERIAL_LOAD_WITH_NAME(mFrames3D)
        & RK_SERIAL_LOAD_WITH_NAME(mDependentGenCoords)
        & RK_SERIAL_LOAD_WITH_NAME(mDependent2DFrames)
        & RK_SERIAL_LOAD_WITH_NAME(mDependent3DFrames)
        & RK_SERIAL_LOAD_WITH_NAME(mModel);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(manipulator_kinematics_model,0xC210FFFE,1,"manipulator_kinematics_model",kte_map)

};





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
class manipulator_dynamics_model : public manipulator_kinematics_model, public state_rate_function_with_io<double> {
  protected:
    mass_matrix_calc mMassCalc; ///< Holds the model's mass-matrix calculator.
    
    std::vector< shared_ptr< system_input > > mInputs; ///< Holds the list of system input objects that are part of the KTE model.
    std::vector< shared_ptr< system_output > > mOutputs; ///< Holds the list of system output objects that are part of the KTE model.
    
    
  public:
    
    using manipulator_kinematics_model::operator<<;

    /**
     * Default constructor.
     */
    manipulator_dynamics_model(const std::string& aName = "") : manipulator_kinematics_model(aName),
                                                                mMassCalc(aName + "_mass_calc") { };
    
    /**
     * Default destructor.
     */
    virtual ~manipulator_dynamics_model() { };
    
    /**
     * Gets the manipulator KTE mass-matrix calculator used by this object.
     * \return The manipulator KTE mass-matrix calculator used by this object.
     */
    const mass_matrix_calc& getMassCalc() const { return mMassCalc; };

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
    
    /**
     * Add a generalized inertial element.
     * \param aInertiaGen a generalized inertial element to add.
     * \return reference to this.
     */
    virtual manipulator_dynamics_model& operator <<(const shared_ptr< inertia_gen >& aInertiaGen);

    /**
     * Add a 2D inertial element.
     * \param aInertia2D a 2D inertial element to add.
     * \return reference to this.
     */
    virtual manipulator_dynamics_model& operator <<(const shared_ptr< inertia_2D >& aInertia2D);

    /**
     * Add a 3D inertial element.
     * \param aInertia3D a 3D inertial element to add.
     * \return reference to this.
     */
    virtual manipulator_dynamics_model& operator <<(const shared_ptr< inertia_3D >& aInertia3D);
    
    /**
     * Add a system input.
     * \param aInput A KTE system input to add.
     * \return reference to this.
     */
    virtual manipulator_dynamics_model& operator <<(const shared_ptr< system_input >& aInput);

    /**
     * Add a system output.
     * \param aOutput A KTE system output to add.
     * \return reference to this.
     */
    virtual manipulator_dynamics_model& operator <<(const shared_ptr< system_output >& aOutput);
    
    /**
     * Get the total number of state values for all the joint frames concatenated.
     * \return The total number of state values for all the joint frames concatenated.
     */
    unsigned int getJointStatesCount() const {
      return getJointPositionsCount() + getJointVelocitiesCount();
    };
    
    /**
     * Get the total number of state values for all the dependent frames concatenated.
     * \return The total number of state values for all the dependent frames concatenated.
     */
    unsigned int getDependentStatesCount() const {
      return getDependentPositionsCount() + getDependentVelocitiesCount();
    };
    
    /**
     * Get the total number of KTE system inputs, if all concatenated in one vector.
     * \return The total number of KTE system inputs, if all concatenated in one vector.
     */
    unsigned int getInputsCount() const {
      unsigned int result = 0;
      for(std::vector< shared_ptr< system_input > >::const_iterator it = mInputs.begin(); it != mInputs.end(); ++it)
        result += (*it)->getInputCount();
      return result;
    };
    
    /**
     * Get the total number of KTE system outputs, if all concatenated in one vector.
     * \return The total number of KTE system outputs, if all concatenated in one vector.
     */
    unsigned int getOutputsCount() const {
      unsigned int result = 0;
      for(std::vector< shared_ptr< system_output > >::const_iterator it = mOutputs.begin(); it != mOutputs.end(); ++it)
        result += (*it)->getOutputCount();
      return result;
    };
    
    
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
    vect_n<double> getJointStates() const;
    
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
    void setJointStates(const vect_n<double>& aJointStates);
    
    /**
     * This function computes the output-vector corresponding to a state vector.
     * 
     * \param aTime current integration time
     * \param aState current state vector
     * \param aOutput holds, as output, the output-vector
     */
    virtual void RK_CALL computeOutput(double aTime,const ReaK::vect_n<double>& aState, ReaK::vect_n<double>& aOutput);
    
    /**
     * This function sets the input-vector.
     *
     * \param aInput current input-vector
     */
    virtual void RK_CALL setInput(const ReaK::vect_n<double>& aInput);
    
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
    virtual void RK_CALL computeStateRate(double aTime,const ReaK::vect_n<double>& aState, ReaK::vect_n<double>& aStateRate);
    
    /**
     * Obtain a vector that contains all the dependent states concatenated into one vector.
     * The ordering in the vector is as follows: ( Generalized Coordinates, 2D Poses, 3D Poses ), 
     * where the dependent frames are sorted in the same order as in the container returned by DependentCoords(), 
     * DependentFrames2D() and DependentFrames3D(), respectively. In other words, the first 2 * DependentCoords().size() elements
     * are the states of dependent generalized coordinates, the next 7 * DependentFrames2D().size() 
     * elements are the states of the 2D dependent frames, and the final 13 * DependentFrames3D().size()
     * are the states of the 3D dependent frames.
     * \return All the dependent states concatenated into one vector.
     */
    vect_n<double> getDependentStates() const;
    
    
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
      manipulator_kinematics_model::save(A,manipulator_kinematics_model::getStaticObjectType()->TypeVersion());
      state_rate_function<double>::save(A,state_rate_function<double>::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mMassCalc);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      manipulator_kinematics_model::load(A,manipulator_kinematics_model::getStaticObjectType()->TypeVersion());
      state_rate_function<double>::load(A,state_rate_function<double>::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mMassCalc);
    };

    RK_RTTI_MAKE_CONCRETE_2BASE(manipulator_dynamics_model,0xC210FFFF,1,"manipulator_dynamics_model",manipulator_kinematics_model,state_rate_function<double>)
    
};





};

};

#endif











