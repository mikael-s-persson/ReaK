/**
 * \file CRS_A465_models.hpp
 *
 * This library has constructs of the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2012
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


#ifndef RK_CRS_A465_MODELS_HPP
#define RK_CRS_A465_MODELS_HPP


#include "mbd_kte/inertia.hpp"
#include "mbd_kte/spring.hpp"
#include "mbd_kte/damper.hpp"
#include "mbd_kte/rigid_link.hpp"
#include "mbd_kte/revolute_joint.hpp"
#include "mbd_kte/prismatic_joint.hpp"
#include "mbd_kte/kte_map_chain.hpp"
#include "mbd_kte/force_actuator.hpp"
#include "mbd_kte/joint_friction.hpp"
#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "kinetostatics/motion_jacobians.hpp"
#include "mbd_kte/jacobian_joint_map.hpp"

#include "mbd_kte/manipulator_model.hpp"

#include "topologies/joint_space_topologies.hpp"
#include "topologies/se3_topologies.hpp"
#include "topologies/joint_space_limits.hpp"


namespace ReaK {


namespace robot_airship {


/**
 * This class serves to store (load / save) the data that models the CRS A465 robot. This class 
 * also provides a number of functions to create KTE chains, mass-matrix calculators, manipulator
 * models (kinematics and dynamics) as well as joint or end-effector spaces which can be used 
 * as tangent-bundle topologies in path-planning code available in ReaK.
 */
class CRS_A465_model_builder {
  public:
    
    //declare all the intermediate frames.
    shared_ptr< frame_3D<double> > robot_base;
    shared_ptr< frame_3D<double> > track_joint_end;
    shared_ptr< frame_3D<double> > arm_joint_1_base;
    shared_ptr< frame_3D<double> > arm_joint_1_end;
    shared_ptr< frame_3D<double> > arm_joint_2_base;
    shared_ptr< frame_3D<double> > arm_joint_2_end;
    shared_ptr< frame_3D<double> > arm_joint_3_base;
    shared_ptr< frame_3D<double> > arm_joint_3_end;
    shared_ptr< frame_3D<double> > arm_joint_4_base;
    shared_ptr< frame_3D<double> > arm_joint_4_end;
    shared_ptr< frame_3D<double> > arm_joint_5_base;
    shared_ptr< frame_3D<double> > arm_joint_5_end;
    shared_ptr< frame_3D<double> > arm_joint_6_base;
    shared_ptr< frame_3D<double> > arm_joint_6_end;
    shared_ptr< frame_3D<double> > arm_EE;

    //declare all the joint coordinates.
    shared_ptr< gen_coord<double> > track_joint_coord;
    shared_ptr< gen_coord<double> > arm_joint_1_coord;
    shared_ptr< gen_coord<double> > arm_joint_2_coord;
    shared_ptr< gen_coord<double> > arm_joint_3_coord;
    shared_ptr< gen_coord<double> > arm_joint_4_coord; 
    shared_ptr< gen_coord<double> > arm_joint_5_coord; 
    shared_ptr< gen_coord<double> > arm_joint_6_coord; 
  
    //declare all the joint jacobians.
    shared_ptr< jacobian_gen_3D<double> > track_joint_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_1_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_2_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_3_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_4_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_5_jacobian;
    shared_ptr< jacobian_gen_3D<double> > arm_joint_6_jacobian;
    
    //declare all the joints.
    shared_ptr< kte::prismatic_joint_3D > track_joint;
    shared_ptr< kte::revolute_joint_3D > arm_joint_1;
    shared_ptr< kte::revolute_joint_3D > arm_joint_2;
    shared_ptr< kte::revolute_joint_3D > arm_joint_3;
    shared_ptr< kte::revolute_joint_3D > arm_joint_4;
    shared_ptr< kte::revolute_joint_3D > arm_joint_5;
    shared_ptr< kte::revolute_joint_3D > arm_joint_6;
    
    //declare all motor inertias
    shared_ptr< kte::inertia_gen > track_joint_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_1_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_2_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_3_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_4_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_5_inertia;
    shared_ptr< kte::inertia_gen > arm_joint_6_inertia;
  
    //declare all force actuators
    shared_ptr< kte::driving_actuator_gen > track_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_1_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_2_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_3_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_4_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_5_actuator;
    shared_ptr< kte::driving_actuator_gen > arm_joint_6_actuator;
  
    //declare all links
    shared_ptr< kte::rigid_link_3D > link_0;
    shared_ptr< kte::rigid_link_3D > link_1;
    shared_ptr< kte::rigid_link_3D > link_2;
    shared_ptr< kte::rigid_link_3D > link_3;
    shared_ptr< kte::rigid_link_3D > link_4;
    shared_ptr< kte::rigid_link_3D > link_5;
    shared_ptr< kte::rigid_link_3D > link_6;
  
    //declare all link-dependent frames
    shared_ptr< kte::joint_dependent_frame_3D > link_0_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_1_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_2_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_3_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_4_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_5_dep_frame;
    shared_ptr< kte::joint_dependent_frame_3D > link_6_dep_frame;

    //declare all link inertias
    shared_ptr< kte::inertia_3D > link_0_inertia;
    shared_ptr< kte::inertia_3D > link_1_inertia;
    shared_ptr< kte::inertia_3D > link_2_inertia;
    shared_ptr< kte::inertia_3D > link_3_inertia;
    shared_ptr< kte::inertia_3D > link_4_inertia;
    shared_ptr< kte::inertia_3D > link_5_inertia;
    shared_ptr< kte::inertia_3D > link_6_inertia;

    vect_n<double> joint_lower_bounds;
    vect_n<double> joint_upper_bounds;
    pp::joint_limits_collection<double> joint_rate_limits;
    vect_n<double> preferred_posture;
    
    typedef pp::metric_space_array< pp::rl_joint_space_2nd_order<double>::type, 7>::type rate_limited_joint_space_type;
    typedef pp::metric_space_array< pp::joint_space_2nd_order<double>::type, 7>::type joint_space_type;
    typedef pp::se3_2nd_order_topology<double>::type end_effector_space_type;
    
    
    /**
     * Default constructor.
     */
    CRS_A465_model_builder() { };
    
    /**
     * This function will load all the KTE model data from the given file.
     * \param aFileName The absolute or relative path to the model file to load (xml format only).
     */
    void load_kte_from_file(const std::string& aFileName);
    
    /**
     * This function will load all the joint-limit data from the given file.
     * \param aFileName The absolute or relative path to the model file to load (xml format only).
     */
    void load_limits_from_file(const std::string& aFileName);
    
    /**
     * This function will create all the model data from a preset definition of the model.
     * This function may be useful to either create the preset CRS A465 model or to create a 
     * complete model that can be saved to an xml file and modified subsequently.
     */
    void create_from_preset();
    
    /**
     * This function will save all the KTE model data to the given file.
     * \param aFileName The absolute or relative path to the model file to which to save (xml format only).
     */
    void save_kte_to_file(const std::string& aFileName) const;
    
    /**
     * This function will save all the joint-limit data to the given file.
     * \param aFileName The absolute or relative path to the model file to which to save (xml format only).
     */
    void save_limits_to_file(const std::string& aFileName) const;
    
    /**
     * This function will construct a KTE chain that represents the kinematics of the CRS A465 robot.
     * \return A KTE chain that represents the kinematics of the CRS A465 robot.
     */
    shared_ptr< kte::kte_map_chain > get_kinematics_kte_chain() const;
    
    /**
     * This function will construct a KTE chain that represents the dynamics of the CRS A465 robot.
     * \return A KTE chain that represents the dynamics of the CRS A465 robot.
     */
    shared_ptr< kte::kte_map_chain > get_dynamics_kte_chain() const;
    
    /**
     * This enum type is used to select the kind of inertias that should be considered in the 
     * mass-matrix calculator (see get_mass_matrix_calculator).
     */
    enum inertia_sources {
      link_inertia = 1, ///< This means the inertias of the links should be considered.
      motor_inertia = 2 ///< This means the motor inertias of the joints should be considered.
    };
    
    /**
     * This function constructs and outputs the mass matrix calculator for the given type of 
     * inertia sources.
     * \param aInertiaSources The inertia sources to include in the mass matrix calculator (can be any bitwise-or combination of the values in the inertia_sources enum type).
     * \return The mass-matrix calculator object for the specified inertia sources.
     */
    shared_ptr< kte::mass_matrix_calc > get_mass_matrix_calculator( int aInertiaSources = link_inertia | motor_inertia ) const;
    
    /**
     * This enum type is used to specify which dependent frames (outputs of direct kinematics) should be 
     * considered when creating a manipulator kinematics model (see get_manipulator_kin_model).
     */
    enum dependent_frames {
      link_frames = 1, ///< This means that all frames attached to links should be an output of the kinematics model.
      end_effector_frame = 2 ///< This means that the end-effector frame should be an output of the kinematics model.
    };
    
    /**
     * This function constructs and outputs a manipulator kinematics model object for the specified 
     * dependent frames.
     * \param aDependentFrames The output frames for the kinematics model, can be any bitwise-or combination of values of the dependent_frames enum-type.
     * \return A manipulator kinematics model object for the specified dependent frames.
     */
    shared_ptr< kte::manipulator_kinematics_model > get_manipulator_kin_model( int aDependentFrames = end_effector_frame ) const;
    
    /**
     * This function constructs and outputs a manipulator dynamics model object for the specified 
     * dependent frames.
     * \param aDependentFrames The output frames for the dynamics model, can be any bitwise-or combination of values of the dependent_frames enum-type.
     * \return A manipulator dynamics model object for the specified dependent frames.
     */
    shared_ptr< kte::manipulator_dynamics_model > get_manipulator_dyn_model( int aDependentFrames = end_effector_frame ) const;
    
    /**
     * This function construct a rate-limited joint-space on which the joint vectors can reside. The 
     * joint space is a topology (see pp::TopologyConcept) which is also differentiable (see pp::TangentBundleConcept),
     * and can serve joint states which are stored as reach-time values (i.e. position normalized by speed, etc.).
     * \return A rate-limited joint-space corresponding to the rate-limits and joint-limits stored in this model.
     */
    rate_limited_joint_space_type get_rl_joint_space() const;
    
    /**
     * This function construct a joint-space on which the joint vectors can reside. The 
     * joint space is a topology (see pp::TopologyConcept) which is also differentiable (see pp::TangentBundleConcept).
     * \return A joint-space corresponding to the rate-limits and joint-limits stored in this model.
     */
    joint_space_type get_joint_space() const;
    
    /**
     * This function construct an end-effector space on which the end-effector state can be represented. The 
     * end-effector space is a topology (see pp::TopologyConcept) which is also differentiable (see pp::TangentBundleConcept).
     * \return An end-effector space on which the end-effector state can be represented (note that the bounds of the end-effector spaces are very approximate).
     */
    end_effector_space_type get_end_effector_space() const;
    
    
    
    
    
    
};


};
  
  
};


#endif






