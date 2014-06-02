/**
 * \file free_floating_platform.hpp
 * 
 * This library declares a class to represent the kinematics model of a free-floating platform in 2D or 3D.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date December 2013
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

#ifndef REAK_FREE_FLOATING_PLATFORM_HPP
#define REAK_FREE_FLOATING_PLATFORM_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/ctrl/mbd_kte/kte_map_chain.hpp>
#include "inverse_kinematics_model.hpp"

namespace ReaK {

namespace kte {


/**
 * This class that models a 2D free-floating platform.
 * This class is only a kinematics model (which is rather trivial).
 */
class free_floater_2D_kinematics : public inverse_kinematics_model {
  protected:
    shared_ptr< frame_2D<double> > m_base_frame;
    shared_ptr< frame_2D<double> > m_state_frame;
    shared_ptr< jacobian_2D_2D<double> > m_state_jacobian;
    shared_ptr< frame_2D<double> > m_output_frame;
    mutable std::vector< shared_ptr< joint_dependent_frame_2D > > m_EEs;
    std::vector< pose_2D<double> > m_EEposes;  // relative to output-frame
    
    shared_ptr< kte_map_chain > m_chain;
    
    void resyncEndEffectors() const;
    void getJacobianMatrixAndDerivativeImpl(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>* JacDot) const;
    
  public:
    
    virtual shared_ptr< kte_map_chain > getKTEChain() const { return m_chain; };
    
    BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = 3);
    
    shared_ptr< frame_2D<double> > getBaseFrame() const { return m_base_frame; };
    void setBaseFrame(const shared_ptr< frame_2D<double> >& aBaseFrame) { 
      m_base_frame = aBaseFrame; 
      m_output_frame->Parent = aBaseFrame;
    };
    
    shared_ptr< frame_2D<double> > getStateFrame() const { return m_state_frame; };
    void setStateFrame(const shared_ptr< frame_2D<double> >& aStateFrame) { 
      for(std::size_t i = 0; i < m_EEs.size(); ++i) {
        m_EEs[i]->mUpStream2DJoints.erase(m_state_frame);
        m_EEs[i]->mUpStream2DJoints[aStateFrame] = m_state_jacobian;
      };
      m_state_frame = aStateFrame; 
    };
    
    shared_ptr< frame_2D<double> > getOutputFrame() const { return m_output_frame; };
    
    const std::vector< pose_2D<double> >& getEEPoses() const { return m_EEposes; };
    std::vector< pose_2D<double> >& getEEPoses() { return m_EEposes; };
    
    /**
     * Default constructor.
     */
    free_floater_2D_kinematics(const std::string& aName = "",
                               const shared_ptr< frame_2D<double> >& aBaseFrame = shared_ptr< frame_2D<double> >());
    
    virtual ~free_floater_2D_kinematics() { };
    
    virtual std::size_t getJointPositionsCount() const { return 4; };
    
    virtual std::size_t getJointVelocitiesCount() const { return 3; };
    
    virtual std::size_t getJointAccelerationsCount() const { return 3; };
    
    virtual std::size_t getDependentPositionsCount() const { return 4 * m_EEposes.size(); };
    
    virtual std::size_t getDependentVelocitiesCount() const { return 3 * m_EEposes.size(); };
    
    virtual std::size_t getDependentAccelerationsCount() const { return 3 * m_EEposes.size(); };
    
    virtual std::size_t getFrames2DCount() const { return 1; };
    
    virtual shared_ptr< frame_2D<double> > getFrame2D(std::size_t i) const { 
      return m_state_frame;
    };
    
    virtual std::size_t getDependentFrames2DCount() const { return m_EEposes.size(); };
    
    virtual shared_ptr< joint_dependent_frame_2D > getDependentFrame2D(std::size_t i) const { 
      resyncEndEffectors();
      return m_EEs[i];
    };
    
    virtual void doDirectMotion();
    
    virtual void doInverseMotion();
    
    virtual void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const;
    
    virtual void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const;
    
    virtual vect_n<double> getJointPositionLowerBounds() const { return vect_n<double>(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()); };
    
    virtual void setJointPositionLowerBounds(const vect_n<double>& aJointLowerBounds) { };
    
    virtual vect_n<double> getJointPositionUpperBounds() const { return vect_n<double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()); };
    
    virtual void setJointPositionUpperBounds(const vect_n<double>& aJointUpperBounds) { };
    
    virtual vect_n<double> getJointPositions() const;
    
    virtual void setJointPositions(const vect_n<double>& aJointPositions);
    
    virtual vect_n<double> getJointVelocities() const;
    
    virtual void setJointVelocities(const vect_n<double>& aJointVelocities);
    
    virtual vect_n<double> getJointAccelerations() const;
    
    virtual void setJointAccelerations(const vect_n<double>& aJointAccelerations);
    
    virtual vect_n<double> getDependentPositions() const;
    
    virtual vect_n<double> getDependentVelocities() const;
    
    virtual vect_n<double> getDependentAccelerations() const;
    
    virtual void setDependentPositions(const vect_n<double>& aDepPositions);
    
    virtual void setDependentVelocities(const vect_n<double>& aDepVelocities);
    
    virtual void setDependentAccelerations(const vect_n<double>& aDepAccelerations);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(free_floater_2D_kinematics,0xC210005B,1,"free_floater_2D_kinematics",inverse_kinematics_model)
    
};



/**
 * This class that models a 3D manipulator with 3 revolute joints (i.e., shoulder-elbow-wrist joints
 * all aligned along the z-axis). This class is only a kinematics model.
 */
class free_floater_3D_kinematics : public inverse_kinematics_model {
  protected:
    shared_ptr< frame_3D<double> > m_base_frame;
    shared_ptr< frame_3D<double> > m_state_frame;
    shared_ptr< jacobian_3D_3D<double> > m_state_jacobian;
    shared_ptr< frame_3D<double> > m_output_frame;
    mutable std::vector< shared_ptr< joint_dependent_frame_3D > > m_EEs;
    std::vector< pose_3D<double> > m_EEposes;  // relative to output-frame
    
    shared_ptr< kte_map_chain > m_chain;
    
    void resyncEndEffectors() const;
    void getJacobianMatrixAndDerivativeImpl(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>* JacDot) const;
    
  public:
    
    virtual shared_ptr< kte_map_chain > getKTEChain() const { return m_chain; };
    
    BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = 6);
    
    shared_ptr< frame_3D<double> > getBaseFrame() const { return m_base_frame; };
    void setBaseFrame(const shared_ptr< frame_3D<double> >& aBaseFrame) { 
      m_base_frame = aBaseFrame; 
      m_output_frame->Parent = aBaseFrame;
    };
    
    shared_ptr< frame_3D<double> > getStateFrame() const { return m_state_frame; };
    void setStateFrame(const shared_ptr< frame_3D<double> >& aStateFrame) { 
      for(std::size_t i = 0; i < m_EEs.size(); ++i) {
        m_EEs[i]->mUpStream3DJoints.erase(m_state_frame);
        m_EEs[i]->mUpStream3DJoints[aStateFrame] = m_state_jacobian;
      };
      m_state_frame = aStateFrame; 
    };
    
    shared_ptr< frame_3D<double> > getOutputFrame() const { return m_output_frame; };
    
    const std::vector< pose_3D<double> >& getEEPoses() const { return m_EEposes; };
    std::vector< pose_3D<double> >& getEEPoses() { return m_EEposes; };
    
    /**
     * Default constructor.
     */
    free_floater_3D_kinematics(const std::string& aName = "",
                               const shared_ptr< frame_3D<double> >& aBaseFrame = shared_ptr< frame_3D<double> >());
    
    virtual ~free_floater_3D_kinematics() { };
    
    virtual std::size_t getJointPositionsCount() const { return 7; };
    
    virtual std::size_t getJointVelocitiesCount() const { return 6; };
    
    virtual std::size_t getJointAccelerationsCount() const { return 6; };
    
    virtual std::size_t getDependentPositionsCount() const { return 7 * m_EEposes.size(); };
    
    virtual std::size_t getDependentVelocitiesCount() const { return 6 * m_EEposes.size(); };
    
    virtual std::size_t getDependentAccelerationsCount() const { return 6 * m_EEposes.size(); };
    
    virtual std::size_t getFrames3DCount() const { return 1; };
    
    virtual shared_ptr< frame_3D<double> > getFrame3D(std::size_t i) const { 
      return m_state_frame;
    };
    
    virtual std::size_t getDependentFrames3DCount() const { return m_EEposes.size(); };
    
    virtual shared_ptr< joint_dependent_frame_3D > getDependentFrame3D(std::size_t i) const { 
      resyncEndEffectors();
      return m_EEs[i];
    };
    
    virtual void doDirectMotion();
    
    virtual void doInverseMotion();
    
    virtual void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const;
    
    virtual void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const;
    
    virtual vect_n<double> getJointPositionLowerBounds() const { return vect_n<double>(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()); };
    
    virtual void setJointPositionLowerBounds(const vect_n<double>& aJointLowerBounds) { };
    
    virtual vect_n<double> getJointPositionUpperBounds() const { return vect_n<double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()); };
    
    virtual void setJointPositionUpperBounds(const vect_n<double>& aJointUpperBounds) { };
    
    virtual vect_n<double> getJointPositions() const;
    
    virtual void setJointPositions(const vect_n<double>& aJointPositions);
    
    virtual vect_n<double> getJointVelocities() const;
    
    virtual void setJointVelocities(const vect_n<double>& aJointVelocities);
    
    virtual vect_n<double> getJointAccelerations() const;
    
    virtual void setJointAccelerations(const vect_n<double>& aJointAccelerations);
    
    virtual vect_n<double> getDependentPositions() const;
    
    virtual vect_n<double> getDependentVelocities() const;
    
    virtual vect_n<double> getDependentAccelerations() const;
    
    virtual void setDependentPositions(const vect_n<double>& aDepPositions);
    
    virtual void setDependentVelocities(const vect_n<double>& aDepVelocities);
    
    virtual void setDependentAccelerations(const vect_n<double>& aDepAccelerations);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(free_floater_3D_kinematics,0xC210005C,1,"free_floater_3D_kinematics",inverse_kinematics_model)
    
};



};

};

#endif











