/**
 * \file manip_ERA_arm.hpp
 * 
 * This library declares a class to represent a kte-based model of a ERA manipulator, i.e., 
 * the European Robotic Arm manipulator.
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

#ifndef REAK_MANIP_ERA_ARM_HPP
#define REAK_MANIP_ERA_ARM_HPP

#include <ReaK/core/base/defs.hpp>
#include "inverse_kinematics_model.hpp"
#include <ReaK/mbd/kte/kte_map_chain.hpp>

namespace ReaK {

namespace kte {



/**
 * This class to represent a kte-based model of a ERA manipulator, i.e., 
 * the European Robotic Arm manipulator.
 */
class manip_ERA_kinematics : public inverse_kinematics_model {
  private:
    shared_ptr< frame_3D<double> > m_base_frame;
    std::vector< shared_ptr< gen_coord<double> > > m_joints;
    shared_ptr< joint_dependent_frame_3D > m_EE;
    vect_n<double> link_lengths;
    shared_ptr< kte_map_chain > m_chain;
    
  public:
    
    BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = 7);
    
    vect<double,3> preferred_elbow_dir;
    vect_n<double> joint_lower_bounds;
    vect_n<double> joint_upper_bounds;
    
    shared_ptr< kte_map_chain > getKTEChain() const { return m_chain; };
    
    /**
     * Default constructor.
     */
    manip_ERA_kinematics(const std::string& aName = "",
                         const shared_ptr< frame_3D<double> >& aBaseFrame = shared_ptr< frame_3D<double> >(new frame_3D<double>()),
                         const vect_n<double>& aLinkLengths = vect_n<double>(1.228, 0.340, 4.072, 4.072, 0.340, 1.228), 
                         const vect<double,3>& aPreferredElbowDir = (vect<double,3>(0.0, 1.0, 0.0)), 
                         const vect_n<double>& aJointLowerBounds = vect_n<double>(-M_PI, -0.75 * M_PI, -0.75 * M_PI, -0.75 * M_PI, -0.75 * M_PI, -0.75 * M_PI, -M_PI),
                         const vect_n<double>& aJointUpperBounds = vect_n<double>( M_PI,  0.75 * M_PI,  0.75 * M_PI,  0.75 * M_PI,  0.75 * M_PI,  0.75 * M_PI,  M_PI));
    
    virtual ~manip_ERA_kinematics() { };
    
    virtual std::size_t getJointPositionsCount() const { return 7; };
    
    virtual std::size_t getJointVelocitiesCount() const { return 7; };
    
    virtual std::size_t getJointAccelerationsCount() const { return 7; };
    
    virtual std::size_t getDependentPositionsCount() const { return 7; };
    
    virtual std::size_t getDependentVelocitiesCount() const { return 6; };
    
    virtual std::size_t getDependentAccelerationsCount() const { return 6; };
    
    virtual std::size_t getCoordsCount() const { return 7; };
    
    virtual shared_ptr< gen_coord<double> > getCoord(std::size_t i) const { 
      return m_joints[i];
    };
    
    virtual std::size_t getDependentFrames3DCount() const { return 1; };
    
    virtual shared_ptr< joint_dependent_frame_3D > getDependentFrame3D(std::size_t i) const { 
      return m_EE;
    };
    
    virtual void doDirectMotion();
    
    virtual void doInverseMotion();
    
    virtual void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const;
    
    virtual void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const;
    
    virtual vect_n<double> getJointPositionLowerBounds() const { return joint_lower_bounds; };
    
    virtual void setJointPositionLowerBounds(const vect_n<double>& aJointLowerBounds) { joint_lower_bounds = aJointLowerBounds; };
    
    virtual vect_n<double> getJointPositionUpperBounds() const { return joint_upper_bounds; };
    
    virtual void setJointPositionUpperBounds(const vect_n<double>& aJointUpperBounds) { joint_upper_bounds = aJointUpperBounds; };
    
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
    
    RK_RTTI_MAKE_CONCRETE_1BASE(manip_ERA_kinematics,0xC2100057,1,"manip_ERA_kinematics",inverse_kinematics_model)
    
};



};

};

#endif











