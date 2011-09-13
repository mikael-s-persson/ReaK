
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

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/state_measures.hpp"
#include "mbd_kte/free_joints.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"


int main() {
  using namespace ReaK;
  
  shared_pointer< frame_2D<double> >::type 
    global_base( new frame_2D<double>(), scoped_deleter());
  
  shared_pointer< frame_2D<double> >::type 
    airship2D_frame( new frame_2D<double>(), scoped_deleter());
  
  shared_pointer< frame_2D<double> >::type 
    airship2D_output_frame( new frame_2D<double>(), scoped_deleter());
    
  shared_pointer< jacobian_2D_2D<double> >::type
    airship2D_joint_jac( new jacobian_2D_2D<double>(), scoped_deleter());
  
  shared_pointer< kte::free_joint_2D >::type
    airship2D_joint( new kte::free_joint_2D("airship2D_joint",
                                          airship2D_frame,
					  global_base,
					  airship2D_output_frame,
					  airship2D_joint_jac), scoped_deleter());
    
  shared_pointer< kte::joint_dependent_frame_2D >::type
    airship2D_dep_frame( new kte::joint_dependent_frame_2D(airship2D_output_frame),
                         scoped_deleter());
  airship2D_dep_frame->add_joint(airship2D_frame,airship2D_joint_jac);
  
  shared_pointer< kte::inertia_2D >::type
    airship2D_inertia( new kte::inertia_2D("airship2D_inertia",
                                           airship2D_dep_frame,
					   1.0,
					   1.0), scoped_deleter());
  
  shared_pointer< kte::mass_matrix_calc >::type
    airship2D_mass_calc( new kte::mass_matrix_calc("airship2D_mass_calc"), scoped_deleter());
    
  (*airship2D_mass_calc) << airship2D_inertia;
  (*airship2D_mass_calc) << airship2D_frame;
  
  shared_pointer< kte::driving_actuator_2D >::type
    airship2D_actuator( new kte::driving_actuator_2D("airship2D_actuator",
                                                   airship2D_frame,
						   airship2D_joint), scoped_deleter());
  
  shared_pointer< kte::position_measure_2D >::type
    airship2D_position( new kte::position_measure_2D("airship2D_position",
                                                     airship2D_frame), scoped_deleter());
  
  shared_pointer< kte::rotation_measure_2D >::type
    airship2D_rotation( new kte::rotation_measure_2D("airship2D_rotation",
                                                     airship2D_frame), scoped_deleter());
  
  shared_pointer< kte::kte_map_chain >::type
    airship2D_model( new kte::kte_map_chain("airship2D_model"), scoped_deleter());
    
  (*airship2D_model) << airship2D_position 
                     << airship2D_rotation 
                     << airship2D_actuator 
                     << airship2D_joint 
                     << airship2D_inertia;
  
  {
    serialization::xml_oarchive out("airship2D_basic.xml");
    out << airship2D_frame
        << airship2D_position
        << airship2D_rotation
        << airship2D_actuator
        << airship2D_model
        << airship2D_mass_calc;
  };
  
  {
    serialization::xml_oarchive out("airship2D_init.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("initial_motion",*airship2D_frame);
  };
  
  {
    mat<double,mat_structure::diagonal> Qu(3);
    Qu(0,0) = 0.1;
    Qu(1,1) = 0.1;
    Qu(2,2) = 0.1;
    serialization::xml_oarchive out("airship2D_Qu.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("input_disturbance",Qu);
  };
  
  {
    mat<double,mat_structure::diagonal> R(4);
    R(0,0) = 0.01;
    R(1,1) = 0.01;
    R(2,2) = 0.01;
    R(3,3) = 0.01;
    serialization::xml_oarchive out("airship2D_R.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("measurement_noise",R);
  };
  
};





