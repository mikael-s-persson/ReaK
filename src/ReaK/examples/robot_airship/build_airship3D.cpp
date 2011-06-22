
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
  
  boost::shared_ptr< frame_3D<double> > 
    global_base( new frame_3D<double>(), scoped_deleter());
  
  boost::shared_ptr< frame_3D<double> > 
    airship3D_frame( new frame_3D<double>(), scoped_deleter());
  
  boost::shared_ptr< frame_3D<double> > 
    airship3D_output_frame( new frame_3D<double>(), scoped_deleter());
    
  boost::shared_ptr< jacobian_3D_3D<double> >
    airship3D_joint_jac( new jacobian_3D_3D<double>(), scoped_deleter());
  
  boost::shared_ptr< kte::free_joint_3D > 
    airship3D_joint( new kte::free_joint_3D("airship3D_joint",
                                            airship3D_frame,
					    global_base,
					    airship3D_output_frame,
					    airship3D_joint_jac), scoped_deleter());
    
  kte::jacobian_joint3D_map_3D airship3D_jacmap;
  airship3D_jacmap[airship3D_frame] = airship3D_joint_jac;
  
  boost::shared_ptr< kte::inertia_3D >
    airship3D_inertia( new kte::inertia_3D("airship3D_inertia",
                                           airship3D_output_frame,
					   1.0,
					   mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
					   kte::jacobian_joint_map_3D(),
					   kte::jacobian_joint2D_map_3D(),
					   airship3D_jacmap), scoped_deleter());
  
  boost::shared_ptr< kte::mass_matrix_calc >
    airship3D_mass_calc( new kte::mass_matrix_calc("airship3D_mass_calc"), scoped_deleter());
    
  (*airship3D_mass_calc) << airship3D_inertia;
  (*airship3D_mass_calc) << airship3D_frame;
  
  boost::shared_ptr< kte::driving_actuator_3D >
    airship3D_actuator( new kte::driving_actuator_3D("airship3D_actuator",
                                                     airship3D_frame,
						     airship3D_joint), scoped_deleter());
  
  boost::shared_ptr< kte::position_measure_3D >
    airship3D_position( new kte::position_measure_3D("airship3D_position",
                                                     airship3D_frame), scoped_deleter());
  
  boost::shared_ptr< kte::rotation_measure_3D >
    airship3D_rotation( new kte::rotation_measure_3D("airship3D_rotation",
                                                     airship3D_frame), scoped_deleter());
  
  boost::shared_ptr< kte::kte_map_chain >
    airship3D_model( new kte::kte_map_chain("airship3D_model"), scoped_deleter());
    
  (*airship3D_model) << airship3D_position 
                     << airship3D_rotation 
                     << airship3D_actuator 
                     << airship3D_joint 
                     << airship3D_inertia;
  
  {
    serialization::xml_oarchive out("airship3D_basic.xml");
    out << airship3D_frame
        << airship3D_position
        << airship3D_rotation
        << airship3D_actuator
        << airship3D_model
        << airship3D_mass_calc;
  };
  
  {
    serialization::xml_oarchive out("airship3D_init.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("initial_motion",*airship3D_frame);
  };
  
  {
    mat<double,mat_structure::diagonal> Qu(6);
    Qu(0,0) = 0.1;
    Qu(1,1) = 0.1;
    Qu(2,2) = 0.1;
    Qu(3,3) = 0.1;
    Qu(4,4) = 0.1;
    Qu(5,5) = 0.1;
    serialization::xml_oarchive out("airship3D_Qu.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("input_disturbance",Qu);
  };
  
  {
    mat<double,mat_structure::diagonal> R(7);
    R(0,0) = 0.01;
    R(1,1) = 0.01;
    R(2,2) = 0.01;
    R(3,3) = 0.01;
    R(3,3) = 0.01;
    R(3,3) = 0.01;
    R(3,3) = 0.01;
    serialization::xml_oarchive out("airship3D_R.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("measurement_noise",R);
  };
  
};






