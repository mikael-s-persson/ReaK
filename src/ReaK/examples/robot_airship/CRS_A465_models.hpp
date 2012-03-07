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


namespace ReaK {


namespace robot_airship {



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

    
    /**
     * Default constructor.
     */
    CRS_A465_model_builder() { };
    
    /**
     * This function will load all the model data from the given file.
     * \param aFileName The absolute or relative path to the model file to load (xml format only).
     */
    void load_from_file(const std::string& aFileName);
    
    /**
     * This function will create all the model data from a preset definition of the model.
     * This function may be useful to either create the preset CRS A465 model or to create a 
     * complete model that can be saved to an xml file and modified subsequently.
     */
    void create_from_preset();
    
    /**
     * This function will save all the model data from the given file.
     * \param aFileName The absolute or relative path to the model file to which to save (xml format only).
     */
    void save_to_file(const std::string& aFileName) const;
    
    
    shared_ptr< kte::kte_map_chain > get_kinematics_kte_chain() const;
    
    
    shared_ptr< kte::kte_map_chain > get_dynamics_kte_chain() const;
    
    enum inertia_sources {
      link_inertia = 1,
      motor_inertia = 2
    };
    
    shared_ptr< kte::mass_matrix_calc > get_mass_matrix_calculator( int aInertiaSources = link_inertia | motor_inertia ) const;
    
    
    enum dependent_frames {
      link_frames = 1,
      end_effector_frame = 2
    };
    
    shared_ptr< kte::manipulator_kinematics_model > get_manipulator_kin_model( int aDependentFrames = end_effector_frame ) const;
    
    
    shared_ptr< kte::manipulator_dynamics_model > get_manipulator_dyn_model( int aDependentFrames = end_effector_frame ) const;
    
    
    

};


};
  
  
};


#endif






