/**
 *\file py_kte_models.cpp
 *
 * This source file defines export functions for the python bindings for classes of
 * the ReaK KTE-based Models library.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date October 2012
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



#include "base/defs.hpp"

#include "direct_kinematics_model.hpp"   // DONE.
#include "inverse_kinematics_model.hpp"  // DONE.
#include "inverse_dynamics_model.hpp"    // DONE.
#include "manip_kinematics_model.hpp"    // DONE.
#include "manip_dynamics_model.hpp"      // DONE.

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "base/py_fixes.hpp"

#include <sstream>




namespace PyReaK {
  


ReaK::vect_n<double> py_dyn_mdl_compute_output(ReaK::kte::manipulator_dynamics_model& aMdl, double aTime, const ReaK::vect_n<double>& aState) {
  ReaK::vect_n<double> aOutput;
  aMdl.computeOutput(aTime, aState, aOutput);
  return aOutput;
};

ReaK::vect_n<double> py_dyn_mdl_compute_state_rate(ReaK::kte::manipulator_dynamics_model& aMdl, double aTime, const ReaK::vect_n<double>& aState) {
  ReaK::vect_n<double> aStateRate;
  aMdl.computeStateRate(aTime, aState, aStateRate);
  return aStateRate;
};


void export_kte_models() {

  using namespace boost::python;
  
    
  class_< ReaK::kte::direct_kinematics_model,
          boost::noncopyable,
          bases< ReaK::named_object >,
          ReaK::shared_ptr< ReaK::kte::direct_kinematics_model >
        >("KTEDirectKinematics")
    .def(init<std::string>())
    .add_property("joint_positions", &ReaK::kte::direct_kinematics_model::getJointPositions, &ReaK::kte::direct_kinematics_model::setJointPositions)
    .add_property("joint_velocities", &ReaK::kte::direct_kinematics_model::getJointVelocities, &ReaK::kte::direct_kinematics_model::setJointVelocities)
    .add_property("joint_accelerations", &ReaK::kte::direct_kinematics_model::getJointAccelerations, &ReaK::kte::direct_kinematics_model::setJointAccelerations)
    .add_property("dependent_positions", &ReaK::kte::direct_kinematics_model::getDependentPositions)
    .add_property("dependent_velocities", &ReaK::kte::direct_kinematics_model::getDependentVelocities)
    .add_property("dependent_accelerations", &ReaK::kte::direct_kinematics_model::getDependentAccelerations)
    .def("coords_count", &ReaK::kte::direct_kinematics_model::getCoordsCount)
    .def("coord", &ReaK::kte::direct_kinematics_model::getCoord)
    .def("frames_2D_count", &ReaK::kte::direct_kinematics_model::getFrames2DCount)
    .def("frame_2D", &ReaK::kte::direct_kinematics_model::getFrame2D)
    .def("frames_3D_count", &ReaK::kte::direct_kinematics_model::getFrames3DCount)
    .def("frame_3D", &ReaK::kte::direct_kinematics_model::getFrame3D)
    .def("dependent_coords_count", &ReaK::kte::direct_kinematics_model::getDependentCoordsCount)
    .def("dependent_coord", &ReaK::kte::direct_kinematics_model::getDependentCoord)
    .def("dependent_frames_2D_count", &ReaK::kte::direct_kinematics_model::getDependentFrames2DCount)
    .def("dependent_frame_2D", &ReaK::kte::direct_kinematics_model::getDependentFrame2D)
    .def("dependent_frames_3D_count", &ReaK::kte::direct_kinematics_model::getDependentFrames3DCount)
    .def("dependent_frame_3D", &ReaK::kte::direct_kinematics_model::getDependentFrame3D)
    .def("do_direct_motion", &ReaK::kte::direct_kinematics_model::doDirectMotion);
    
  class_< ReaK::kte::inverse_kinematics_model,
          boost::noncopyable,
          bases< ReaK::kte::direct_kinematics_model >,
          ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model >
        >("KTEInverseKinematics")
    .def(init<std::string>())
    .def("set_dependent_positions", &ReaK::kte::inverse_kinematics_model::setDependentPositions)
    .def("set_dependent_velocities", &ReaK::kte::inverse_kinematics_model::setDependentVelocities)
    .def("set_dependent_accelerations", &ReaK::kte::inverse_kinematics_model::setDependentAccelerations)
    .def("do_inverse_motion", &ReaK::kte::inverse_kinematics_model::doInverseMotion);
    
  class_< ReaK::kte::manipulator_kinematics_model,
          boost::noncopyable,
          bases< ReaK::kte::direct_kinematics_model >,
          ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model >
        >("KTEManipKinematics")
    .def(init<std::string>())
    .def("add_dependent_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .add_property("kte_model", &ReaK::kte::manipulator_kinematics_model::getModel, &ReaK::kte::manipulator_kinematics_model::setModel);
    
  class_< ReaK::kte::inverse_dynamics_model,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model >
        >("KTEInverseDynamics")
    .def(init<std::string>())
    .add_property("joint_states", &ReaK::kte::inverse_dynamics_model::getJointStates, &ReaK::kte::inverse_dynamics_model::getJointStates)
    .def("joint_states_count", &ReaK::kte::inverse_dynamics_model::getJointStatesCount);
    
    
  class_< ReaK::kte::manipulator_dynamics_model,
          boost::noncopyable,
          bases< ReaK::kte::manipulator_kinematics_model, ReaK::kte::inverse_dynamics_model >,
          ReaK::shared_ptr< ReaK::kte::manipulator_dynamics_model >
        >("KTEManipDynamics")
    .def(init<std::string>())
    .def("add_dependent_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_frame_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_frame_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_frame_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_inertia_gen", static_cast< ReaK::kte::manipulator_dynamics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::inertia_gen >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_inertia_2D",  static_cast< ReaK::kte::manipulator_dynamics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::inertia_2D >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_inertia_3D",  static_cast< ReaK::kte::manipulator_dynamics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::inertia_3D >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_system_input", static_cast< ReaK::kte::manipulator_dynamics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::system_input >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .def("add_system_output",  static_cast< ReaK::kte::manipulator_dynamics_model& (ReaK::kte::manipulator_dynamics_model::*)(const ReaK::shared_ptr< ReaK::kte::system_output >&) >(&ReaK::kte::manipulator_dynamics_model::operator<<), return_internal_reference<>())
    .add_property("inputs", &ReaK::kte::manipulator_dynamics_model::getInput, &ReaK::kte::manipulator_dynamics_model::setInput)
    .add_property("dependent_states", &ReaK::kte::manipulator_dynamics_model::getDependentStates)
    .def("mass_calculator", &ReaK::kte::manipulator_dynamics_model::getMassCalc, return_internal_reference<>())
    .def("compute_output", &py_dyn_mdl_compute_output)
    .def("compute_state_rate", &py_dyn_mdl_compute_state_rate);
    
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::direct_kinematics_model >, 
                          ReaK::shared_ptr< ReaK::named_object > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model >, 
                          ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model >, 
                          ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::manipulator_dynamics_model >, 
                          ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::manipulator_dynamics_model >, 
                          ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model > >();
                          
#endif
  
};


};
















