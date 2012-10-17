/**
 *\file py_mbd_kte.cpp
 *
 * This source file defines export functions for the python bindings for classes of
 * the ReaK Multi-Body Dynamics with Kinetostatic Transmission Elements library.
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

#include "kte_map.hpp"                // DONE.
#include "kte_map_chain.hpp"          // DONE.
#include "kte_system_input.hpp"       // DONE.
#include "kte_system_output.hpp"      // DONE.
#include "damper.hpp"                 // DONE.
#include "driving_actuator.hpp"       // DONE.
#include "dry_revolute_joint.hpp"     // DONE.
#include "flexible_beam.hpp"          // DONE.
#include "force_actuator.hpp"         // DONE.
#include "free_joints.hpp"            // DONE.
#include "inertia.hpp"                // DONE.
#include "inertial_beam.hpp"          // DONE.
#include "jacobian_joint_map.hpp"     // DONE.
//#include "kte_ext_mappings.hpp"       // not needed.
#include "line_point_mindist.hpp"     // DONE.
#include "manipulator_model.hpp"      // DONE.
#include "mass_matrix_calculator.hpp" // DONE.
#include "plane_point_mindist.hpp"    // DONE.
#include "prismatic_joint.hpp"        // DONE.
#include "reacting_kte.hpp"           // DONE.
#include "revolute_joint.hpp"         // DONE.
#include "rigid_link.hpp"             // DONE.
#include "spring.hpp"                 // DONE.
#include "state_controls.hpp"         // DONE.
#include "state_measures.hpp"         // DONE.
#include "torsion_damper.hpp"         // DONE.
#include "torsion_spring.hpp"         // DONE.
#include "virtual_kte_interface.hpp"  // DONE.
#include "vmc_revolute_joint.hpp"     // DONE.

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "base/py_fixes.hpp"

#include <sstream>




namespace PyReaK {
  
  
template <typename T>
struct py_jac_joint_map {
    typedef typename T::key_type K;
    typedef typename T::mapped_type V;
    static V get(T const& x, K const& i) {
      typename T::const_iterator it = x.find(i);
      if( it != x.end() ) 
        return it->second;
      PyErr_SetString(PyExc_KeyError, "Key not found");
      return V();
    };
    static void set(T& x, K const& i, V const& v) {
      x[i] = v; // use map autocreation feature
    };
    static void del(T& x, K const& i) {
      if( x.find(i) != x.end() ) 
        x.erase(i);
      else 
        PyErr_SetString(PyExc_KeyError, "Key not found");
    };
};

template <typename FrameType, const std::vector< ReaK::shared_ptr< FrameType > >& (ReaK::kte::mass_matrix_calc::*MemFuncPtr)() const>
struct py_mmc_frame_vector {
  const ReaK::kte::mass_matrix_calc* mmc;
  
  py_mmc_frame_vector(const ReaK::kte::mass_matrix_calc* aMMC) : mmc(aMMC) { };
  
  static py_mmc_frame_vector create(const ReaK::kte::mass_matrix_calc& x) { return py_mmc_frame_vector(&x); };
  
  std::size_t size() const {
    return (mmc->*MemFuncPtr)().size();
  };
  
  ReaK::shared_ptr< FrameType > get(std::size_t i) const {
    return (mmc->*MemFuncPtr)()[i];
  };
};


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


template <typename FrameType, const std::vector< ReaK::shared_ptr< FrameType > >& (ReaK::kte::manipulator_kinematics_model::*MemFuncPtr)() const>
struct py_mdl_frame_vector {
  const ReaK::kte::manipulator_kinematics_model* kin_mdl;
  
  py_mdl_frame_vector(const ReaK::kte::manipulator_kinematics_model* aKinMdl) : kin_mdl(aKinMdl) { };
  
  static py_mdl_frame_vector create(const ReaK::kte::manipulator_kinematics_model& x) { return py_mdl_frame_vector(&x); };
  
  std::size_t size() const {
    return (kin_mdl->*MemFuncPtr)().size();
  };
  
  ReaK::shared_ptr< FrameType > get(std::size_t i) const {
    return (kin_mdl->*MemFuncPtr)()[i];
  };
};


void export_mbd_kte() {

  using namespace boost::python;
  
  
  
  /********************************************************************************
   *        Base classes
   * *****************************************************************************/
  
  
  class_< ReaK::kte::kte_map,
          boost::noncopyable,
          bases< ReaK::named_object >,
          ReaK::shared_ptr< ReaK::kte::kte_map >
        >("KTEMap", no_init)
    .def("do_motion",   pure_virtual(&ReaK::kte::kte_map::doMotion))
    .def("do_force",    pure_virtual(&ReaK::kte::kte_map::doForce))
    .def("clear_force", pure_virtual(&ReaK::kte::kte_map::clearForce));
  
  class_< ReaK::kte::kte_map_chain,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::kte_map_chain >
        >("KTEMapChain")
    .def(init<std::string>())
    .def("do_motion",   &ReaK::kte::kte_map_chain::doMotion)
    .def("do_force",    &ReaK::kte::kte_map_chain::doForce)
    .def("clear_force", &ReaK::kte::kte_map_chain::clearForce)
    .def("add_kte_map", &ReaK::kte::kte_map_chain::operator<<, return_internal_reference<>());
  
  class_< ReaK::kte::system_input,
          boost::noncopyable,
          bases< ReaK::named_object >,
          ReaK::shared_ptr< ReaK::kte::system_input >
        >("KTESystemInput", no_init)
    .def("get_input_count",  pure_virtual(&ReaK::kte::system_input::getInputCount))
    .def("get_input", pure_virtual(&ReaK::kte::system_input::getInput))
    .def("set_input", pure_virtual(&ReaK::kte::system_input::setInput));
    
  class_< ReaK::kte::system_output,
          boost::noncopyable,
          bases< ReaK::named_object >,
          ReaK::shared_ptr< ReaK::kte::system_output >
        >("KTESystemOutput", no_init)
    .def("get_output_count",   &ReaK::kte::system_output::getOutputCount)
    .def("get_output", pure_virtual(&ReaK::kte::system_output::getOutput));
  
  
  class_< ReaK::kte::damper_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::damper_gen >
        >("KTEDamperGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> >, 
              ReaK::shared_ptr< ReaK::gen_coord<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::damper_gen::doMotion)
    .def("do_force",    &ReaK::kte::damper_gen::doForce)
    .def("clear_force", &ReaK::kte::damper_gen::clearForce)
    .add_property("anchor1", &ReaK::kte::damper_gen::Anchor1, &ReaK::kte::damper_gen::setAnchor1)
    .add_property("anchor2", &ReaK::kte::damper_gen::Anchor2, &ReaK::kte::damper_gen::setAnchor2)
    .add_property("damping", &ReaK::kte::damper_gen::Damping, &ReaK::kte::damper_gen::setDamping);
  
  class_< ReaK::kte::damper_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::damper_2D>
        >("KTEDamper2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::damper_2D::doMotion)
    .def("do_force",    &ReaK::kte::damper_2D::doForce)
    .def("clear_force", &ReaK::kte::damper_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::damper_2D::Anchor1, &ReaK::kte::damper_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::damper_2D::Anchor2, &ReaK::kte::damper_2D::setAnchor2)
    .add_property("damping", &ReaK::kte::damper_2D::Damping, &ReaK::kte::damper_2D::setDamping);
  
  class_< ReaK::kte::damper_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::damper_3D >
        >("KTEDamper3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::damper_3D::doMotion)
    .def("do_force",    &ReaK::kte::damper_3D::doForce)
    .def("clear_force", &ReaK::kte::damper_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::damper_3D::Anchor1, &ReaK::kte::damper_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::damper_3D::Anchor2, &ReaK::kte::damper_3D::setAnchor2)
    .add_property("damping", &ReaK::kte::damper_3D::Damping, &ReaK::kte::damper_3D::setDamping);
  
  
  class_< ReaK::kte::spring_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::spring_gen >
        >("KTESpringGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> >, 
              ReaK::shared_ptr< ReaK::gen_coord<double> >, 
              double, double, double>())
    .def("do_motion",   &ReaK::kte::spring_gen::doMotion)
    .def("do_force",    &ReaK::kte::spring_gen::doForce)
    .def("clear_force", &ReaK::kte::spring_gen::clearForce)
    .add_property("anchor1", &ReaK::kte::spring_gen::Anchor1, &ReaK::kte::spring_gen::setAnchor1)
    .add_property("anchor2", &ReaK::kte::spring_gen::Anchor2, &ReaK::kte::spring_gen::setAnchor2)
    .add_property("rest_length", &ReaK::kte::spring_gen::RestLength, &ReaK::kte::spring_gen::setRestLength)
    .add_property("stiffness", &ReaK::kte::spring_gen::Stiffness, &ReaK::kte::spring_gen::setStiffness)
    .add_property("saturation", &ReaK::kte::spring_gen::Saturation, &ReaK::kte::spring_gen::setSaturation);
  
  class_< ReaK::kte::spring_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::spring_2D>
        >("KTESpring2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double, double, double>())
    .def("do_motion",   &ReaK::kte::spring_2D::doMotion)
    .def("do_force",    &ReaK::kte::spring_2D::doForce)
    .def("clear_force", &ReaK::kte::spring_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::spring_2D::Anchor1, &ReaK::kte::spring_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::spring_2D::Anchor2, &ReaK::kte::spring_2D::setAnchor2)
    .add_property("rest_length", &ReaK::kte::spring_2D::RestLength, &ReaK::kte::spring_2D::setRestLength)
    .add_property("stiffness", &ReaK::kte::spring_2D::Stiffness, &ReaK::kte::spring_2D::setStiffness)
    .add_property("saturation", &ReaK::kte::spring_2D::Saturation, &ReaK::kte::spring_2D::setSaturation);
  
  class_< ReaK::kte::spring_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::spring_3D >
        >("KTESpring3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double, double, double>())
    .def("do_motion",   &ReaK::kte::spring_3D::doMotion)
    .def("do_force",    &ReaK::kte::spring_3D::doForce)
    .def("clear_force", &ReaK::kte::spring_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::spring_3D::Anchor1, &ReaK::kte::spring_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::spring_3D::Anchor2, &ReaK::kte::spring_3D::setAnchor2)
    .add_property("rest_length", &ReaK::kte::spring_3D::RestLength, &ReaK::kte::spring_3D::setRestLength)
    .add_property("stiffness", &ReaK::kte::spring_3D::Stiffness, &ReaK::kte::spring_3D::setStiffness)
    .add_property("saturation", &ReaK::kte::spring_3D::Saturation, &ReaK::kte::spring_3D::setSaturation);
  
  
  class_< ReaK::kte::torsion_damper_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::torsion_damper_2D>
        >("KTETorsionDamper2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::torsion_damper_2D::doMotion)
    .def("do_force",    &ReaK::kte::torsion_damper_2D::doForce)
    .def("clear_force", &ReaK::kte::torsion_damper_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::torsion_damper_2D::Anchor1, &ReaK::kte::torsion_damper_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::torsion_damper_2D::Anchor2, &ReaK::kte::torsion_damper_2D::setAnchor2)
    .add_property("damping", &ReaK::kte::torsion_damper_2D::Damping, &ReaK::kte::torsion_damper_2D::setDamping);
  
  class_< ReaK::kte::torsion_damper_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::torsion_damper_3D >
        >("KTETorsionDamper3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::torsion_damper_3D::doMotion)
    .def("do_force",    &ReaK::kte::torsion_damper_3D::doForce)
    .def("clear_force", &ReaK::kte::torsion_damper_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::torsion_damper_3D::Anchor1, &ReaK::kte::torsion_damper_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::torsion_damper_3D::Anchor2, &ReaK::kte::torsion_damper_3D::setAnchor2)
    .add_property("damping", &ReaK::kte::torsion_damper_3D::Damping, &ReaK::kte::torsion_damper_3D::setDamping);
    
  
  class_< ReaK::kte::torsion_spring_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::torsion_spring_2D>
        >("KTETorsionSpring2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double, double>())
    .def("do_motion",   &ReaK::kte::torsion_spring_2D::doMotion)
    .def("do_force",    &ReaK::kte::torsion_spring_2D::doForce)
    .def("clear_force", &ReaK::kte::torsion_spring_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::torsion_spring_2D::Anchor1, &ReaK::kte::torsion_spring_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::torsion_spring_2D::Anchor2, &ReaK::kte::torsion_spring_2D::setAnchor2)
    .add_property("stiffness", &ReaK::kte::torsion_spring_2D::Stiffness, &ReaK::kte::torsion_spring_2D::setStiffness)
    .add_property("saturation", &ReaK::kte::torsion_spring_2D::Saturation, &ReaK::kte::torsion_spring_2D::setSaturation);
  
  class_< ReaK::kte::torsion_spring_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::torsion_spring_3D >
        >("KTETorsionSpring3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double, double>())
    .def("do_motion",   &ReaK::kte::torsion_spring_3D::doMotion)
    .def("do_force",    &ReaK::kte::torsion_spring_3D::doForce)
    .def("clear_force", &ReaK::kte::torsion_spring_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::torsion_spring_3D::Anchor1, &ReaK::kte::torsion_spring_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::torsion_spring_3D::Anchor2, &ReaK::kte::torsion_spring_3D::setAnchor2)
    .add_property("stiffness", &ReaK::kte::torsion_spring_3D::Stiffness, &ReaK::kte::torsion_spring_3D::setStiffness)
    .add_property("saturation", &ReaK::kte::torsion_spring_3D::Saturation, &ReaK::kte::torsion_spring_3D::setSaturation);
    
    
  class_< ReaK::kte::flexible_beam_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::flexible_beam_2D>
        >("KTEFlexibleBeam2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double, double, double>())
    .def("do_motion",   &ReaK::kte::flexible_beam_2D::doMotion)
    .def("do_force",    &ReaK::kte::flexible_beam_2D::doForce)
    .def("clear_force", &ReaK::kte::flexible_beam_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::flexible_beam_2D::Anchor1, &ReaK::kte::flexible_beam_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::flexible_beam_2D::Anchor2, &ReaK::kte::flexible_beam_2D::setAnchor2)
    .add_property("center_frame", &ReaK::kte::flexible_beam_2D::CenterFrame, &ReaK::kte::flexible_beam_2D::setCenterFrame)
    .add_property("rest_length", &ReaK::kte::flexible_beam_2D::RestLength, &ReaK::kte::flexible_beam_2D::setRestLength)
    .add_property("stiffness", &ReaK::kte::flexible_beam_2D::Stiffness, &ReaK::kte::flexible_beam_2D::setStiffness)
    .add_property("torsion_stiffness", &ReaK::kte::flexible_beam_2D::TorsionStiffness, &ReaK::kte::flexible_beam_2D::setTorsionStiffness);
  
  class_< ReaK::kte::flexible_beam_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::flexible_beam_3D >
        >("KTEFlexibleBeam3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double, double, double>())
    .def("do_motion",   &ReaK::kte::flexible_beam_3D::doMotion)
    .def("do_force",    &ReaK::kte::flexible_beam_3D::doForce)
    .def("clear_force", &ReaK::kte::flexible_beam_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::flexible_beam_3D::Anchor1, &ReaK::kte::flexible_beam_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::flexible_beam_3D::Anchor2, &ReaK::kte::flexible_beam_3D::setAnchor2)
    .add_property("center_frame", &ReaK::kte::flexible_beam_3D::CenterFrame, &ReaK::kte::flexible_beam_3D::setCenterFrame)
    .add_property("rest_length", &ReaK::kte::flexible_beam_3D::RestLength, &ReaK::kte::flexible_beam_3D::setRestLength)
    .add_property("stiffness", &ReaK::kte::flexible_beam_3D::Stiffness, &ReaK::kte::flexible_beam_3D::setStiffness)
    .add_property("torsion_stiffness", &ReaK::kte::flexible_beam_3D::TorsionStiffness, &ReaK::kte::flexible_beam_3D::setTorsionStiffness);
    
    
  class_< ReaK::kte::line_point_mindist_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::line_point_mindist_2D>
        >("KTELinePointMinDist2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::vect<double,2>, ReaK::vect<double,2> >())
    .def("do_motion",   &ReaK::kte::line_point_mindist_2D::doMotion)
    .def("do_force",    &ReaK::kte::line_point_mindist_2D::doForce)
    .def("clear_force", &ReaK::kte::line_point_mindist_2D::clearForce)
    .add_property("base_frame", &ReaK::kte::line_point_mindist_2D::BaseFrame, &ReaK::kte::line_point_mindist_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::line_point_mindist_2D::EndFrame, &ReaK::kte::line_point_mindist_2D::setEndFrame)
    .add_property("tangent", &ReaK::kte::line_point_mindist_2D::Tangent, &ReaK::kte::line_point_mindist_2D::setTangent)
    .add_property("origin_mindist", &ReaK::kte::line_point_mindist_2D::OriginMinDist, &ReaK::kte::line_point_mindist_2D::setOriginMinDist);
  
  class_< ReaK::kte::line_point_mindist_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::line_point_mindist_3D >
        >("KTELinePointMinDist3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::vect<double,3>, ReaK::vect<double,3> >())
    .def("do_motion",   &ReaK::kte::line_point_mindist_3D::doMotion)
    .def("do_force",    &ReaK::kte::line_point_mindist_3D::doForce)
    .def("clear_force", &ReaK::kte::line_point_mindist_3D::clearForce)
    .add_property("base_frame", &ReaK::kte::line_point_mindist_3D::BaseFrame, &ReaK::kte::line_point_mindist_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::line_point_mindist_3D::EndFrame, &ReaK::kte::line_point_mindist_3D::setEndFrame)
    .add_property("tangent", &ReaK::kte::line_point_mindist_3D::Tangent, &ReaK::kte::line_point_mindist_3D::setTangent)
    .add_property("origin_mindist", &ReaK::kte::line_point_mindist_3D::OriginMinDist, &ReaK::kte::line_point_mindist_3D::setOriginMinDist);
    
  class_< ReaK::kte::plane_point_mindist_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::plane_point_mindist_3D >
        >("KTEPlanePointMinDist3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::vect<double,3>, double >())
    .def("do_motion",   &ReaK::kte::plane_point_mindist_3D::doMotion)
    .def("do_force",    &ReaK::kte::plane_point_mindist_3D::doForce)
    .def("clear_force", &ReaK::kte::plane_point_mindist_3D::clearForce)
    .add_property("base_frame", &ReaK::kte::plane_point_mindist_3D::BaseFrame, &ReaK::kte::plane_point_mindist_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::plane_point_mindist_3D::EndFrame, &ReaK::kte::plane_point_mindist_3D::setEndFrame)
    .add_property("normal", &ReaK::kte::plane_point_mindist_3D::Normal, &ReaK::kte::plane_point_mindist_3D::setNormal)
    .add_property("origin_dist", &ReaK::kte::plane_point_mindist_3D::Origin, &ReaK::kte::plane_point_mindist_3D::setOrigin);
    
      
      
  class_< ReaK::kte::reacting_kte_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::reacting_kte_gen >
        >("KTEReactingGen", no_init)
    .def("apply_reaction_force",   pure_virtual(&ReaK::kte::reacting_kte_gen::applyReactionForce));
      
  class_< ReaK::kte::reacting_kte_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::reacting_kte_2D >
        >("KTEReacting2D", no_init)
    .def("apply_reaction_force",   pure_virtual(&ReaK::kte::reacting_kte_2D::applyReactionForce));
      
  class_< ReaK::kte::reacting_kte_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::reacting_kte_3D >
        >("KTEReacting3D", no_init)
    .def("apply_reaction_force",   pure_virtual(&ReaK::kte::reacting_kte_3D::applyReactionForce));
  
  
  class_< ReaK::kte::revolute_joint_2D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_gen >,
          ReaK::shared_ptr< ReaK::kte::revolute_joint_2D>
        >("KTERevoluteJoint2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_2D<double> > >())
    .def("do_motion",   &ReaK::kte::revolute_joint_2D::doMotion)
    .def("do_force",    &ReaK::kte::revolute_joint_2D::doForce)
    .def("clear_force", &ReaK::kte::revolute_joint_2D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::revolute_joint_2D::applyReactionForce)
    .add_property("angle", &ReaK::kte::revolute_joint_2D::Angle, &ReaK::kte::revolute_joint_2D::setAngle)
    .add_property("base_frame", &ReaK::kte::revolute_joint_2D::BaseFrame, &ReaK::kte::revolute_joint_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::revolute_joint_2D::EndFrame, &ReaK::kte::revolute_joint_2D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::revolute_joint_2D::Jacobian, &ReaK::kte::revolute_joint_2D::setJacobian);
  
  class_< ReaK::kte::revolute_joint_3D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_gen >,
          ReaK::shared_ptr< ReaK::kte::revolute_joint_3D >
        >("KTERevoluteJoint3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::vect<double,3>,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_3D<double> > >())
    .def("do_motion",   &ReaK::kte::revolute_joint_3D::doMotion)
    .def("do_force",    &ReaK::kte::revolute_joint_3D::doForce)
    .def("clear_force", &ReaK::kte::revolute_joint_3D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::revolute_joint_3D::applyReactionForce)
    .add_property("angle", &ReaK::kte::revolute_joint_3D::Angle, &ReaK::kte::revolute_joint_3D::setAngle)
    .add_property("axis", &ReaK::kte::revolute_joint_3D::Axis, &ReaK::kte::revolute_joint_3D::setAxis)
    .add_property("base_frame", &ReaK::kte::revolute_joint_3D::BaseFrame, &ReaK::kte::revolute_joint_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::revolute_joint_3D::EndFrame, &ReaK::kte::revolute_joint_3D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::revolute_joint_3D::Jacobian, &ReaK::kte::revolute_joint_3D::setJacobian);
    
  
  class_< ReaK::kte::dry_revolute_joint_2D,
          boost::noncopyable,
          bases< ReaK::kte::revolute_joint_2D >,
          ReaK::shared_ptr< ReaK::kte::dry_revolute_joint_2D>
        >("KTEDryRevoluteJoint2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_2D<double> >,
              double, double, double >())
    .def("do_motion",   &ReaK::kte::dry_revolute_joint_2D::doMotion)
    .def("do_force",    &ReaK::kte::dry_revolute_joint_2D::doForce)
    .def("clear_force", &ReaK::kte::dry_revolute_joint_2D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::dry_revolute_joint_2D::applyReactionForce)
    .add_property("stiction_coef", &ReaK::kte::dry_revolute_joint_2D::StictionCoefficient, &ReaK::kte::dry_revolute_joint_2D::setStictionCoefficient)
    .add_property("slip_coef", &ReaK::kte::dry_revolute_joint_2D::SlipCoefficient, &ReaK::kte::dry_revolute_joint_2D::setSlipCoefficient)
    .add_property("slip_velocity", &ReaK::kte::dry_revolute_joint_2D::SlipVelocity, &ReaK::kte::dry_revolute_joint_2D::setSlipVelocity);
  
  class_< ReaK::kte::dry_revolute_joint_3D,
          boost::noncopyable,
          bases< ReaK::kte::revolute_joint_3D >,
          ReaK::shared_ptr< ReaK::kte::dry_revolute_joint_3D >
        >("KTEDryRevoluteJoint3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::vect<double,3>,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_3D<double> >,
              double, double, double >())
    .def("do_motion",   &ReaK::kte::dry_revolute_joint_3D::doMotion)
    .def("do_force",    &ReaK::kte::dry_revolute_joint_3D::doForce)
    .def("clear_force", &ReaK::kte::dry_revolute_joint_3D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::dry_revolute_joint_3D::applyReactionForce)
    .add_property("stiction_coef", &ReaK::kte::dry_revolute_joint_3D::StictionCoefficient, &ReaK::kte::dry_revolute_joint_3D::setStictionCoefficient)
    .add_property("slip_coef", &ReaK::kte::dry_revolute_joint_3D::SlipCoefficient, &ReaK::kte::dry_revolute_joint_3D::setSlipCoefficient)
    .add_property("slip_velocity", &ReaK::kte::dry_revolute_joint_3D::SlipVelocity, &ReaK::kte::dry_revolute_joint_3D::setSlipVelocity);
    
  
  class_< ReaK::kte::prismatic_joint_2D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_gen >,
          ReaK::shared_ptr< ReaK::kte::prismatic_joint_2D>
        >("KTEPrismaticJoint2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::vect<double,2>,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_2D<double> > >())
    .def("do_motion",   &ReaK::kte::prismatic_joint_2D::doMotion)
    .def("do_force",    &ReaK::kte::prismatic_joint_2D::doForce)
    .def("clear_force", &ReaK::kte::prismatic_joint_2D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::prismatic_joint_2D::applyReactionForce)
    .add_property("coord", &ReaK::kte::prismatic_joint_2D::Coord, &ReaK::kte::prismatic_joint_2D::setCoord)
    .add_property("axis", &ReaK::kte::prismatic_joint_2D::Axis, &ReaK::kte::prismatic_joint_2D::setAxis)
    .add_property("base_frame", &ReaK::kte::prismatic_joint_2D::BaseFrame, &ReaK::kte::prismatic_joint_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::prismatic_joint_2D::EndFrame, &ReaK::kte::prismatic_joint_2D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::prismatic_joint_2D::Jacobian, &ReaK::kte::prismatic_joint_2D::setJacobian);
  
  class_< ReaK::kte::prismatic_joint_3D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_gen >,
          ReaK::shared_ptr< ReaK::kte::prismatic_joint_3D >
        >("KTEPrismaticJoint3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::vect<double,3>,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_3D<double> > >())
    .def("do_motion",   &ReaK::kte::prismatic_joint_3D::doMotion)
    .def("do_force",    &ReaK::kte::prismatic_joint_3D::doForce)
    .def("clear_force", &ReaK::kte::prismatic_joint_3D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::prismatic_joint_3D::applyReactionForce)
    .add_property("coord", &ReaK::kte::prismatic_joint_3D::Coord, &ReaK::kte::prismatic_joint_3D::setCoord)
    .add_property("axis", &ReaK::kte::prismatic_joint_3D::Axis, &ReaK::kte::prismatic_joint_3D::setAxis)
    .add_property("base_frame", &ReaK::kte::prismatic_joint_3D::BaseFrame, &ReaK::kte::prismatic_joint_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::prismatic_joint_3D::EndFrame, &ReaK::kte::prismatic_joint_3D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::prismatic_joint_3D::Jacobian, &ReaK::kte::prismatic_joint_3D::setJacobian);
    
  
  class_< ReaK::kte::free_joint_2D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_2D >,
          ReaK::shared_ptr< ReaK::kte::free_joint_2D>
        >("KTEFreeJoint2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_2D_2D<double> > >())
    .def("do_motion",   &ReaK::kte::free_joint_2D::doMotion)
    .def("do_force",    &ReaK::kte::free_joint_2D::doForce)
    .def("clear_force", &ReaK::kte::free_joint_2D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::free_joint_2D::applyReactionForce)
    .add_property("coord", &ReaK::kte::free_joint_2D::Coord, &ReaK::kte::free_joint_2D::setCoord)
    .add_property("base_frame", &ReaK::kte::free_joint_2D::BaseFrame, &ReaK::kte::free_joint_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::free_joint_2D::EndFrame, &ReaK::kte::free_joint_2D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::free_joint_2D::Jacobian, &ReaK::kte::free_joint_2D::setJacobian);
  
  class_< ReaK::kte::free_joint_3D,
          boost::noncopyable,
          bases< ReaK::kte::reacting_kte_3D >,
          ReaK::shared_ptr< ReaK::kte::free_joint_3D >
        >("KTEFreeJoint3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_3D_3D<double> > >())
    .def("do_motion",   &ReaK::kte::free_joint_3D::doMotion)
    .def("do_force",    &ReaK::kte::free_joint_3D::doForce)
    .def("clear_force", &ReaK::kte::free_joint_3D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::free_joint_3D::applyReactionForce)
    .add_property("coord", &ReaK::kte::free_joint_3D::Coord, &ReaK::kte::free_joint_3D::setCoord)
    .add_property("base_frame", &ReaK::kte::free_joint_3D::BaseFrame, &ReaK::kte::free_joint_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::free_joint_3D::EndFrame, &ReaK::kte::free_joint_3D::setEndFrame)
    .add_property("jacobian", &ReaK::kte::free_joint_3D::Jacobian, &ReaK::kte::free_joint_3D::setJacobian);
    
    
  
  class_< ReaK::kte::rigid_link_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::rigid_link_gen>
        >("KTERigidLinkGen")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              double >())
    .def("do_motion",   &ReaK::kte::rigid_link_gen::doMotion)
    .def("do_force",    &ReaK::kte::rigid_link_gen::doForce)
    .def("clear_force", &ReaK::kte::rigid_link_gen::clearForce)
    .add_property("base_frame", &ReaK::kte::rigid_link_gen::BaseFrame, &ReaK::kte::rigid_link_gen::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::rigid_link_gen::EndFrame, &ReaK::kte::rigid_link_gen::setEndFrame)
    .add_property("offset", &ReaK::kte::rigid_link_gen::Offset, &ReaK::kte::rigid_link_gen::setOffset);
  
  class_< ReaK::kte::rigid_link_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::rigid_link_2D>
        >("KTERigidLink2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::pose_2D<double> >())
    .def("do_motion",   &ReaK::kte::rigid_link_2D::doMotion)
    .def("do_force",    &ReaK::kte::rigid_link_2D::doForce)
    .def("clear_force", &ReaK::kte::rigid_link_2D::clearForce)
    .add_property("base_frame", &ReaK::kte::rigid_link_2D::BaseFrame, &ReaK::kte::rigid_link_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::rigid_link_2D::EndFrame, &ReaK::kte::rigid_link_2D::setEndFrame)
    .add_property("pose_offset", &ReaK::kte::rigid_link_2D::PoseOffset, &ReaK::kte::rigid_link_2D::setPoseOffset);
  
  class_< ReaK::kte::rigid_link_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::rigid_link_3D >
        >("KTERigidLink3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::pose_3D<double> >())
    .def("do_motion",   &ReaK::kte::rigid_link_3D::doMotion)
    .def("do_force",    &ReaK::kte::rigid_link_3D::doForce)
    .def("clear_force", &ReaK::kte::rigid_link_3D::clearForce)
    .add_property("base_frame", &ReaK::kte::rigid_link_3D::BaseFrame, &ReaK::kte::rigid_link_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::rigid_link_3D::EndFrame, &ReaK::kte::rigid_link_3D::setEndFrame)
    .add_property("pose_offset", &ReaK::kte::rigid_link_3D::PoseOffset, &ReaK::kte::rigid_link_3D::setPoseOffset);
    
    
    
    
    
  class_< ReaK::kte::position_control_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::position_control_gen >
        >("KTEPositionCtrlGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> > >())
    .def("do_motion",   &ReaK::kte::position_control_gen::doMotion)
    .def("do_force",    &ReaK::kte::position_control_gen::doForce)
    .def("clear_force", &ReaK::kte::position_control_gen::clearForce)
    .def("get_input_count",   &ReaK::kte::position_control_gen::getInputCount)
    .def("get_input", &ReaK::kte::position_control_gen::getInput)
    .def("set_input", &ReaK::kte::position_control_gen::setInput)
    .add_property("anchor", &ReaK::kte::position_control_gen::Anchor, &ReaK::kte::position_control_gen::setAnchor)
    .add_property("position_desired", &ReaK::kte::position_control_gen::PosDesired, &ReaK::kte::position_control_gen::setPosDesired);
  
  class_< ReaK::kte::position_control_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::position_control_2D >
        >("KTEPositionCtrl2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::position_control_2D::doMotion)
    .def("do_force",    &ReaK::kte::position_control_2D::doForce)
    .def("clear_force", &ReaK::kte::position_control_2D::clearForce)
    .def("get_input_count",   &ReaK::kte::position_control_2D::getInputCount)
    .def("get_input", &ReaK::kte::position_control_2D::getInput)
    .def("set_input", &ReaK::kte::position_control_2D::setInput)
    .add_property("anchor", &ReaK::kte::position_control_2D::Anchor, &ReaK::kte::position_control_2D::setAnchor)
    .add_property("position_desired", &ReaK::kte::position_control_2D::PosDesired, &ReaK::kte::position_control_2D::setPosDesired);
  
  class_< ReaK::kte::position_control_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::position_control_3D >
        >("KTEPositionCtrl3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::position_control_3D::doMotion)
    .def("do_force",    &ReaK::kte::position_control_3D::doForce)
    .def("clear_force", &ReaK::kte::position_control_3D::clearForce)
    .def("get_input_count",   &ReaK::kte::position_control_3D::getInputCount)
    .def("get_input", &ReaK::kte::position_control_3D::getInput)
    .def("set_input", &ReaK::kte::position_control_3D::setInput)
    .add_property("anchor", &ReaK::kte::position_control_3D::Anchor, &ReaK::kte::position_control_3D::setAnchor)
    .add_property("position_desired", &ReaK::kte::position_control_3D::PosDesired, &ReaK::kte::position_control_3D::setPosDesired);
  
    
  class_< ReaK::kte::rotation_control_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::rotation_control_2D >
        >("KTERotationCtrl2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::rotation_control_2D::doMotion)
    .def("do_force",    &ReaK::kte::rotation_control_2D::doForce)
    .def("clear_force", &ReaK::kte::rotation_control_2D::clearForce)
    .def("get_input_count",   &ReaK::kte::rotation_control_2D::getInputCount)
    .def("get_input", &ReaK::kte::rotation_control_2D::getInput)
    .def("set_input", &ReaK::kte::rotation_control_2D::setInput)
    .add_property("anchor", &ReaK::kte::rotation_control_2D::Anchor, &ReaK::kte::rotation_control_2D::setAnchor)
    .add_property("rotation_desired", &ReaK::kte::rotation_control_2D::AngleDesired, &ReaK::kte::rotation_control_2D::setAngleDesired);
  
  class_< ReaK::kte::rotation_control_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::rotation_control_3D >
        >("KTERotationCtrl3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::rotation_control_3D::doMotion)
    .def("do_force",    &ReaK::kte::rotation_control_3D::doForce)
    .def("clear_force", &ReaK::kte::rotation_control_3D::clearForce)
    .def("get_input_count",   &ReaK::kte::rotation_control_3D::getInputCount)
    .def("get_input", &ReaK::kte::rotation_control_3D::getInput)
    .def("set_input", &ReaK::kte::rotation_control_3D::setInput)
    .add_property("anchor", &ReaK::kte::rotation_control_3D::Anchor, &ReaK::kte::rotation_control_3D::setAnchor)
    .add_property("rotation_desired", &ReaK::kte::rotation_control_3D::QuatDesired, &ReaK::kte::rotation_control_3D::setQuatDesired);
  
    
  class_< ReaK::kte::velocity_control_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::velocity_control_gen >
        >("KTEVelocityCtrlGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_control_gen::doMotion)
    .def("do_force",    &ReaK::kte::velocity_control_gen::doForce)
    .def("clear_force", &ReaK::kte::velocity_control_gen::clearForce)
    .def("get_input_count",   &ReaK::kte::velocity_control_gen::getInputCount)
    .def("get_input", &ReaK::kte::velocity_control_gen::getInput)
    .def("set_input", &ReaK::kte::velocity_control_gen::setInput)
    .add_property("anchor", &ReaK::kte::velocity_control_gen::Anchor, &ReaK::kte::velocity_control_gen::setAnchor)
    .add_property("velocity_desired", &ReaK::kte::velocity_control_gen::VelDesired, &ReaK::kte::velocity_control_gen::setVelDesired);
  
  class_< ReaK::kte::velocity_control_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::velocity_control_2D >
        >("KTEVelocityCtrl2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_control_2D::doMotion)
    .def("do_force",    &ReaK::kte::velocity_control_2D::doForce)
    .def("clear_force", &ReaK::kte::velocity_control_2D::clearForce)
    .def("get_input_count",   &ReaK::kte::velocity_control_2D::getInputCount)
    .def("get_input", &ReaK::kte::velocity_control_2D::getInput)
    .def("set_input", &ReaK::kte::velocity_control_2D::setInput)
    .add_property("anchor", &ReaK::kte::velocity_control_2D::Anchor, &ReaK::kte::velocity_control_2D::setAnchor)
    .add_property("velocity_desired", &ReaK::kte::velocity_control_2D::VelDesired, &ReaK::kte::velocity_control_2D::setVelDesired);
  
  class_< ReaK::kte::velocity_control_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::velocity_control_3D >
        >("KTEVelocityCtrl3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_control_3D::doMotion)
    .def("do_force",    &ReaK::kte::velocity_control_3D::doForce)
    .def("clear_force", &ReaK::kte::velocity_control_3D::clearForce)
    .def("get_input_count",   &ReaK::kte::velocity_control_3D::getInputCount)
    .def("get_input", &ReaK::kte::velocity_control_3D::getInput)
    .def("set_input", &ReaK::kte::velocity_control_3D::setInput)
    .add_property("anchor", &ReaK::kte::velocity_control_3D::Anchor, &ReaK::kte::velocity_control_3D::setAnchor)
    .add_property("velocity_desired", &ReaK::kte::velocity_control_3D::VelDesired, &ReaK::kte::velocity_control_3D::setVelDesired);
  
    
  class_< ReaK::kte::ang_velocity_control_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::ang_velocity_control_2D >
        >("KTEAngVelocityCtrl2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::ang_velocity_control_2D::doMotion)
    .def("do_force",    &ReaK::kte::ang_velocity_control_2D::doForce)
    .def("clear_force", &ReaK::kte::ang_velocity_control_2D::clearForce)
    .def("get_input_count",   &ReaK::kte::ang_velocity_control_2D::getInputCount)
    .def("get_input", &ReaK::kte::ang_velocity_control_2D::getInput)
    .def("set_input", &ReaK::kte::ang_velocity_control_2D::setInput)
    .add_property("anchor", &ReaK::kte::ang_velocity_control_2D::Anchor, &ReaK::kte::ang_velocity_control_2D::setAnchor)
    .add_property("ang_velocity_desired", &ReaK::kte::ang_velocity_control_2D::AngVelDesired, &ReaK::kte::ang_velocity_control_2D::setAngVelDesired);
  
  class_< ReaK::kte::ang_velocity_control_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::ang_velocity_control_3D >
        >("KTEAngVelocityCtrl3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::ang_velocity_control_3D::doMotion)
    .def("do_force",    &ReaK::kte::ang_velocity_control_3D::doForce)
    .def("clear_force", &ReaK::kte::ang_velocity_control_3D::clearForce)
    .def("get_input_count",   &ReaK::kte::ang_velocity_control_3D::getInputCount)
    .def("get_input", &ReaK::kte::ang_velocity_control_3D::getInput)
    .def("set_input", &ReaK::kte::ang_velocity_control_3D::setInput)
    .add_property("anchor", &ReaK::kte::ang_velocity_control_3D::Anchor, &ReaK::kte::ang_velocity_control_3D::setAnchor)
    .add_property("ang_velocity_desired", &ReaK::kte::ang_velocity_control_3D::AngVelDesired, &ReaK::kte::ang_velocity_control_3D::setAngVelDesired);
  
    
    
  class_< ReaK::kte::position_measure_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::position_measure_gen >
        >("KTEPositionMeasGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> > >())
    .def("do_motion",   &ReaK::kte::position_measure_gen::doMotion)
    .def("do_force",    &ReaK::kte::position_measure_gen::doForce)
    .def("clear_force", &ReaK::kte::position_measure_gen::clearForce)
    .def("get_output_count",   &ReaK::kte::position_measure_gen::getOutputCount)
    .def("get_output", &ReaK::kte::position_measure_gen::getOutput)
    .add_property("anchor", &ReaK::kte::position_measure_gen::Anchor, &ReaK::kte::position_measure_gen::setAnchor)
    .add_property("position_measure", &ReaK::kte::position_measure_gen::PosMeasure, &ReaK::kte::position_measure_gen::setPosMeasure);
  
  class_< ReaK::kte::position_measure_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::position_measure_2D >
        >("KTEPositionMeas2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::position_measure_2D::doMotion)
    .def("do_force",    &ReaK::kte::position_measure_2D::doForce)
    .def("clear_force", &ReaK::kte::position_measure_2D::clearForce)
    .def("get_output_count",   &ReaK::kte::position_measure_2D::getOutputCount)
    .def("get_output", &ReaK::kte::position_measure_2D::getOutput)
    .add_property("anchor", &ReaK::kte::position_measure_2D::Anchor, &ReaK::kte::position_measure_2D::setAnchor)
    .add_property("position_measure", &ReaK::kte::position_measure_2D::PosMeasure, &ReaK::kte::position_measure_2D::setPosMeasure);
  
  class_< ReaK::kte::position_measure_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::position_measure_3D >
        >("KTEPositionMeas3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::position_measure_3D::doMotion)
    .def("do_force",    &ReaK::kte::position_measure_3D::doForce)
    .def("clear_force", &ReaK::kte::position_measure_3D::clearForce)
    .def("get_output_count",   &ReaK::kte::position_measure_3D::getOutputCount)
    .def("get_output", &ReaK::kte::position_measure_3D::getOutput)
    .add_property("anchor", &ReaK::kte::position_measure_3D::Anchor, &ReaK::kte::position_measure_3D::setAnchor)
    .add_property("position_measure", &ReaK::kte::position_measure_3D::PosMeasure, &ReaK::kte::position_measure_3D::setPosMeasure);
  
    
  class_< ReaK::kte::rotation_measure_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::rotation_measure_2D >
        >("KTERotationMeas2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::rotation_measure_2D::doMotion)
    .def("do_force",    &ReaK::kte::rotation_measure_2D::doForce)
    .def("clear_force", &ReaK::kte::rotation_measure_2D::clearForce)
    .def("get_output_count",   &ReaK::kte::rotation_measure_2D::getOutputCount)
    .def("get_output", &ReaK::kte::rotation_measure_2D::getOutput)
    .add_property("anchor", &ReaK::kte::rotation_measure_2D::Anchor, &ReaK::kte::rotation_measure_2D::setAnchor)
    .add_property("rotation_measure", &ReaK::kte::rotation_measure_2D::AngleMeasure, &ReaK::kte::rotation_measure_2D::setAngleMeasure);
  
  class_< ReaK::kte::rotation_measure_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::rotation_measure_3D >
        >("KTERotationMeas3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::rotation_measure_3D::doMotion)
    .def("do_force",    &ReaK::kte::rotation_measure_3D::doForce)
    .def("clear_force", &ReaK::kte::rotation_measure_3D::clearForce)
    .def("get_output_count",   &ReaK::kte::rotation_measure_3D::getOutputCount)
    .def("get_output", &ReaK::kte::rotation_measure_3D::getOutput)
    .add_property("anchor", &ReaK::kte::rotation_measure_3D::Anchor, &ReaK::kte::rotation_measure_3D::setAnchor)
    .add_property("rotation_measure", &ReaK::kte::rotation_measure_3D::QuatMeasure, &ReaK::kte::rotation_measure_3D::setQuatMeasure);
  
    
  class_< ReaK::kte::velocity_measure_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::velocity_measure_gen >
        >("KTEVelocityMeasGen")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::gen_coord<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_measure_gen::doMotion)
    .def("do_force",    &ReaK::kte::velocity_measure_gen::doForce)
    .def("clear_force", &ReaK::kte::velocity_measure_gen::clearForce)
    .def("get_output_count",   &ReaK::kte::velocity_measure_gen::getOutputCount)
    .def("get_output", &ReaK::kte::velocity_measure_gen::getOutput)
    .add_property("anchor", &ReaK::kte::velocity_measure_gen::Anchor, &ReaK::kte::velocity_measure_gen::setAnchor)
    .add_property("velocity_measure", &ReaK::kte::velocity_measure_gen::VelMeasure, &ReaK::kte::velocity_measure_gen::setVelMeasure);
  
  class_< ReaK::kte::velocity_measure_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::velocity_measure_2D >
        >("KTEVelocityMeas2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_measure_2D::doMotion)
    .def("do_force",    &ReaK::kte::velocity_measure_2D::doForce)
    .def("clear_force", &ReaK::kte::velocity_measure_2D::clearForce)
    .def("get_output_count",   &ReaK::kte::velocity_measure_2D::getOutputCount)
    .def("get_output", &ReaK::kte::velocity_measure_2D::getOutput)
    .add_property("anchor", &ReaK::kte::velocity_measure_2D::Anchor, &ReaK::kte::velocity_measure_2D::setAnchor)
    .add_property("velocity_measure", &ReaK::kte::velocity_measure_2D::VelMeasure, &ReaK::kte::velocity_measure_2D::setVelMeasure);
  
  class_< ReaK::kte::velocity_measure_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::velocity_measure_3D >
        >("KTEVelocityMeas3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::velocity_measure_3D::doMotion)
    .def("do_force",    &ReaK::kte::velocity_measure_3D::doForce)
    .def("clear_force", &ReaK::kte::velocity_measure_3D::clearForce)
    .def("get_output_count",   &ReaK::kte::velocity_measure_3D::getOutputCount)
    .def("get_output", &ReaK::kte::velocity_measure_3D::getOutput)
    .add_property("anchor", &ReaK::kte::velocity_measure_3D::Anchor, &ReaK::kte::velocity_measure_3D::setAnchor)
    .add_property("velocity_measure", &ReaK::kte::velocity_measure_3D::VelMeasure, &ReaK::kte::velocity_measure_3D::setVelMeasure);
  
    
  class_< ReaK::kte::ang_velocity_measure_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_2D >
        >("KTEAngVelocityMeas2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::ang_velocity_measure_2D::doMotion)
    .def("do_force",    &ReaK::kte::ang_velocity_measure_2D::doForce)
    .def("clear_force", &ReaK::kte::ang_velocity_measure_2D::clearForce)
    .def("get_output_count",   &ReaK::kte::ang_velocity_measure_2D::getOutputCount)
    .def("get_output", &ReaK::kte::ang_velocity_measure_2D::getOutput)
    .add_property("anchor", &ReaK::kte::ang_velocity_measure_2D::Anchor, &ReaK::kte::ang_velocity_measure_2D::setAnchor)
    .add_property("ang_velocity_measure", &ReaK::kte::ang_velocity_measure_2D::AngVelMeasure, &ReaK::kte::ang_velocity_measure_2D::setAngVelMeasure);
  
  class_< ReaK::kte::ang_velocity_measure_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map, ReaK::kte::system_output >,
          ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_3D >
        >("KTEAngVelocityMeas3D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::ang_velocity_measure_3D::doMotion)
    .def("do_force",    &ReaK::kte::ang_velocity_measure_3D::doForce)
    .def("clear_force", &ReaK::kte::ang_velocity_measure_3D::clearForce)
    .def("get_output_count",   &ReaK::kte::ang_velocity_measure_3D::getOutputCount)
    .def("get_output", &ReaK::kte::ang_velocity_measure_3D::getOutput)
    .add_property("anchor", &ReaK::kte::ang_velocity_measure_3D::Anchor, &ReaK::kte::ang_velocity_measure_3D::setAnchor)
    .add_property("ang_velocity_measure", &ReaK::kte::ang_velocity_measure_3D::AngVelMeasure, &ReaK::kte::ang_velocity_measure_3D::setAngVelMeasure);
  
  
  class_<ReaK::kte::jacobian_joint_map_gen>("JacJointMapGen")
    .def("__len__", &ReaK::kte::jacobian_joint_map_gen::size)
    .def("clear", &ReaK::kte::jacobian_joint_map_gen::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_gen>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_gen>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_gen>::del);
  
  class_<ReaK::kte::jacobian_joint_map_2D>("JacJointMap2D")
    .def("__len__", &ReaK::kte::jacobian_joint_map_2D::size)
    .def("clear", &ReaK::kte::jacobian_joint_map_2D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_2D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_2D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_2D>::del);
  
  class_<ReaK::kte::jacobian_joint_map_3D>("JacJointMap3D")
    .def("__len__", &ReaK::kte::jacobian_joint_map_3D::size)
    .def("clear", &ReaK::kte::jacobian_joint_map_3D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_3D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_3D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint_map_3D>::del);
  
  class_<ReaK::kte::jacobian_joint2D_map_gen>("JacJoint2DMapGen")
    .def("__len__", &ReaK::kte::jacobian_joint2D_map_gen::size)
    .def("clear", &ReaK::kte::jacobian_joint2D_map_gen::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_gen>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_gen>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_gen>::del);
  
  class_<ReaK::kte::jacobian_joint2D_map_2D>("JacJoint2DMap2D")
    .def("__len__", &ReaK::kte::jacobian_joint2D_map_2D::size)
    .def("clear", &ReaK::kte::jacobian_joint2D_map_2D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_2D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_2D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_2D>::del);
  
  class_<ReaK::kte::jacobian_joint2D_map_3D>("JacJoint2DMap3D")
    .def("__len__", &ReaK::kte::jacobian_joint2D_map_3D::size)
    .def("clear", &ReaK::kte::jacobian_joint2D_map_3D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_3D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_3D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint2D_map_3D>::del);
  
  class_<ReaK::kte::jacobian_joint3D_map_gen>("JacJoint3DMapGen")
    .def("__len__", &ReaK::kte::jacobian_joint3D_map_gen::size)
    .def("clear", &ReaK::kte::jacobian_joint3D_map_gen::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_gen>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_gen>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_gen>::del);
  
  class_<ReaK::kte::jacobian_joint3D_map_2D>("JacJoint3DMap2D")
    .def("__len__", &ReaK::kte::jacobian_joint3D_map_2D::size)
    .def("clear", &ReaK::kte::jacobian_joint3D_map_2D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_2D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_2D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_2D>::del);
  
  class_<ReaK::kte::jacobian_joint3D_map_3D>("JacJoint3DMap3D")
    .def("__len__", &ReaK::kte::jacobian_joint3D_map_3D::size)
    .def("clear", &ReaK::kte::jacobian_joint3D_map_3D::clear)
    .def("__getitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_3D>::get)
    .def("__setitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_3D>::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &py_jac_joint_map<ReaK::kte::jacobian_joint3D_map_3D>::del);
  
  class_< ReaK::kte::joint_dependent_gen_coord,
          bases< ReaK::shared_object >,
          ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >
        >("JointDependentGenCoord")
    .def(init<ReaK::shared_ptr< ReaK::gen_coord<double> >,
              const ReaK::kte::jacobian_joint_map_gen& >())
    .def(init<ReaK::shared_ptr< ReaK::gen_coord<double> >,
              const ReaK::kte::jacobian_joint_map_gen&,
              const ReaK::kte::jacobian_joint2D_map_gen& >())
    .def(init<ReaK::shared_ptr< ReaK::gen_coord<double> >,
              const ReaK::kte::jacobian_joint_map_gen&,
              const ReaK::kte::jacobian_joint2D_map_gen&,
              const ReaK::kte::jacobian_joint3D_map_gen& >())
    .def_readwrite("frame", &ReaK::kte::joint_dependent_gen_coord::mFrame)
    .def_readwrite("upstream_joints", &ReaK::kte::joint_dependent_gen_coord::mUpStreamJoints)
    .def_readwrite("upstream_2Djoints", &ReaK::kte::joint_dependent_gen_coord::mUpStream2DJoints)
    .def_readwrite("upstream_3Djoints", &ReaK::kte::joint_dependent_gen_coord::mUpStream3DJoints)
    .def("add_joint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&, const ReaK::shared_ptr< ReaK::jacobian_gen_gen<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::add_joint), return_internal_reference<>())
    .def("add_2Djoint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_2D_gen<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::add_joint), return_internal_reference<>())
    .def("add_3Djoint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_3D_gen<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::add_joint), return_internal_reference<>())
    .def("remove_joint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::remove_joint), return_internal_reference<>())
    .def("remove_2Djoint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::remove_joint), return_internal_reference<>())
    .def("remove_3Djoint", static_cast<ReaK::kte::joint_dependent_gen_coord& (ReaK::kte::joint_dependent_gen_coord::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&)>(&ReaK::kte::joint_dependent_gen_coord::remove_joint), return_internal_reference<>())
    ;
  
  class_< ReaK::kte::joint_dependent_frame_2D,
          bases< ReaK::shared_object >,
          ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >
        >("JointDependentFrame2D")
    .def(init<ReaK::shared_ptr< ReaK::frame_2D<double> >,
              const ReaK::kte::jacobian_joint_map_2D& >())
    .def(init<ReaK::shared_ptr< ReaK::frame_2D<double> >,
              const ReaK::kte::jacobian_joint_map_2D&,
              const ReaK::kte::jacobian_joint2D_map_2D& >())
    .def(init<ReaK::shared_ptr< ReaK::frame_2D<double> >,
              const ReaK::kte::jacobian_joint_map_2D&,
              const ReaK::kte::jacobian_joint2D_map_2D&,
              const ReaK::kte::jacobian_joint3D_map_2D& >())
    .def_readwrite("frame", &ReaK::kte::joint_dependent_frame_2D::mFrame)
    .def_readwrite("upstream_joints", &ReaK::kte::joint_dependent_frame_2D::mUpStreamJoints)
    .def_readwrite("upstream_2Djoints", &ReaK::kte::joint_dependent_frame_2D::mUpStream2DJoints)
    .def_readwrite("upstream_3Djoints", &ReaK::kte::joint_dependent_frame_2D::mUpStream3DJoints)
    .def("add_joint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&, const ReaK::shared_ptr< ReaK::jacobian_gen_2D<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::add_joint), return_internal_reference<>())
    .def("add_2Djoint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_2D_2D<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::add_joint), return_internal_reference<>())
    .def("add_3Djoint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_3D_2D<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::add_joint), return_internal_reference<>())
    .def("remove_joint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::remove_joint), return_internal_reference<>())
    .def("remove_2Djoint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::remove_joint), return_internal_reference<>())
    .def("remove_3Djoint", static_cast<ReaK::kte::joint_dependent_frame_2D& (ReaK::kte::joint_dependent_frame_2D::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&)>(&ReaK::kte::joint_dependent_frame_2D::remove_joint), return_internal_reference<>())
    ;
  
  class_< ReaK::kte::joint_dependent_frame_3D,
          bases< ReaK::shared_object >,
          ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >
        >("JointDependentFrame3D")
    .def(init<ReaK::shared_ptr< ReaK::frame_3D<double> >,
              const ReaK::kte::jacobian_joint_map_3D& >())
    .def(init<ReaK::shared_ptr< ReaK::frame_3D<double> >,
              const ReaK::kte::jacobian_joint_map_3D&,
              const ReaK::kte::jacobian_joint2D_map_3D& >())
    .def(init<ReaK::shared_ptr< ReaK::frame_3D<double> >,
              const ReaK::kte::jacobian_joint_map_3D&,
              const ReaK::kte::jacobian_joint2D_map_3D&,
              const ReaK::kte::jacobian_joint3D_map_3D& >())
    .def_readwrite("frame", &ReaK::kte::joint_dependent_frame_3D::mFrame)
    .def_readwrite("upstream_joints", &ReaK::kte::joint_dependent_frame_3D::mUpStreamJoints)
    .def_readwrite("upstream_2Djoints", &ReaK::kte::joint_dependent_frame_3D::mUpStream2DJoints)
    .def_readwrite("upstream_3Djoints", &ReaK::kte::joint_dependent_frame_3D::mUpStream3DJoints)
    .def("add_joint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&, const ReaK::shared_ptr< ReaK::jacobian_gen_3D<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::add_joint), return_internal_reference<>())
    .def("add_2Djoint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_2D_3D<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::add_joint), return_internal_reference<>())
    .def("add_3Djoint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&, const ReaK::shared_ptr< ReaK::jacobian_3D_3D<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::add_joint), return_internal_reference<>())
    .def("remove_joint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::remove_joint), return_internal_reference<>())
    .def("remove_2Djoint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::remove_joint), return_internal_reference<>())
    .def("remove_3Djoint", static_cast<ReaK::kte::joint_dependent_frame_3D& (ReaK::kte::joint_dependent_frame_3D::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&)>(&ReaK::kte::joint_dependent_frame_3D::remove_joint), return_internal_reference<>())
    ;
  
    
    
  class_< ReaK::kte::inertia_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inertia_gen >
        >("KTEInertiaGen")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >, double >())
    .def("do_motion",   &ReaK::kte::inertia_gen::doMotion)
    .def("do_force",    &ReaK::kte::inertia_gen::doForce)
    .def("clear_force", &ReaK::kte::inertia_gen::clearForce)
    .add_property("mass", &ReaK::kte::inertia_gen::Mass, &ReaK::kte::inertia_gen::setMass)
    .add_property("center_of_mass", &ReaK::kte::inertia_gen::CenterOfMass, &ReaK::kte::inertia_gen::setCenterOfMass);
  
  class_< ReaK::kte::inertia_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inertia_2D >
        >("KTEInertia2D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >, double, double >())
    .def("do_motion",   &ReaK::kte::inertia_2D::doMotion)
    .def("do_force",    &ReaK::kte::inertia_2D::doForce)
    .def("clear_force", &ReaK::kte::inertia_2D::clearForce)
    .add_property("mass", &ReaK::kte::inertia_2D::Mass, &ReaK::kte::inertia_2D::setMass)
    .add_property("moment_of_inertia", &ReaK::kte::inertia_2D::MomentOfInertia, &ReaK::kte::inertia_2D::setMomentOfInertia)
    .add_property("center_of_mass", &ReaK::kte::inertia_2D::CenterOfMass, &ReaK::kte::inertia_2D::setCenterOfMass);
  
  class_< ReaK::kte::inertia_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inertia_3D >
        >("KTEInertia3D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >, double, const ReaK::mat<double,ReaK::mat_structure::symmetric>& >())
    .def("do_motion",   &ReaK::kte::inertia_3D::doMotion)
    .def("do_force",    &ReaK::kte::inertia_3D::doForce)
    .def("clear_force", &ReaK::kte::inertia_3D::clearForce)
    .add_property("mass", &ReaK::kte::inertia_3D::Mass, &ReaK::kte::inertia_3D::setMass)
    .add_property("inertia_tensor", &ReaK::kte::inertia_3D::InertiaTensor, &ReaK::kte::inertia_3D::setInertiaTensor)
    .add_property("center_of_mass", &ReaK::kte::inertia_3D::CenterOfMass, &ReaK::kte::inertia_3D::setCenterOfMass);
  
  
  class_< ReaK::kte::force_actuator_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::force_actuator_gen >
        >("KTEForceActuatorGen")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::gen_coord<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >())
    .def("do_motion",   &ReaK::kte::force_actuator_gen::doMotion)
    .def("do_force",    &ReaK::kte::force_actuator_gen::doForce)
    .def("clear_force", &ReaK::kte::force_actuator_gen::clearForce)
    .add_property("frame", &ReaK::kte::force_actuator_gen::Frame, &ReaK::kte::force_actuator_gen::setFrame)
    .add_property("joint", &ReaK::kte::force_actuator_gen::Joint, &ReaK::kte::force_actuator_gen::setJoint);
  
  class_< ReaK::kte::force_actuator_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::force_actuator_2D >
        >("KTEForceActuator2D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::frame_2D<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_2D > >())
    .def("do_motion",   &ReaK::kte::force_actuator_2D::doMotion)
    .def("do_force",    &ReaK::kte::force_actuator_2D::doForce)
    .def("clear_force", &ReaK::kte::force_actuator_2D::clearForce)
    .add_property("frame", &ReaK::kte::force_actuator_2D::Frame, &ReaK::kte::force_actuator_2D::setFrame)
    .add_property("joint", &ReaK::kte::force_actuator_2D::Joint, &ReaK::kte::force_actuator_2D::setJoint);
  
  class_< ReaK::kte::force_actuator_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::force_actuator_3D >
        >("KTEForceActuator3D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::frame_3D<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_3D > >())
    .def("do_motion",   &ReaK::kte::force_actuator_3D::doMotion)
    .def("do_force",    &ReaK::kte::force_actuator_3D::doForce)
    .def("clear_force", &ReaK::kte::force_actuator_3D::clearForce)
    .add_property("frame", &ReaK::kte::force_actuator_3D::Frame, &ReaK::kte::force_actuator_3D::setFrame)
    .add_property("joint", &ReaK::kte::force_actuator_3D::Joint, &ReaK::kte::force_actuator_3D::setJoint);
  
  
  class_< ReaK::kte::driving_actuator_gen,
          boost::noncopyable,
          bases< ReaK::kte::force_actuator_gen, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::driving_actuator_gen >
        >("KTEDrivingActuatorGen")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::gen_coord<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >())
    .def("do_motion",   &ReaK::kte::driving_actuator_gen::doMotion)
    .def("do_force",    &ReaK::kte::driving_actuator_gen::doForce)
    .def("clear_force", &ReaK::kte::driving_actuator_gen::clearForce)
    .def("get_input_count",   &ReaK::kte::driving_actuator_gen::getInputCount)
    .def("get_input", &ReaK::kte::driving_actuator_gen::getInput)
    .def("set_input", &ReaK::kte::driving_actuator_gen::setInput)
    .add_property("drive_force", &ReaK::kte::driving_actuator_gen::DriveForce, &ReaK::kte::driving_actuator_gen::setDriveForce);
  
  class_< ReaK::kte::driving_actuator_2D,
          boost::noncopyable,
          bases< ReaK::kte::force_actuator_2D, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::driving_actuator_2D >
        >("KTEDrivingActuator2D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::frame_2D<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_2D > >())
    .def("do_motion",   &ReaK::kte::driving_actuator_2D::doMotion)
    .def("do_force",    &ReaK::kte::driving_actuator_2D::doForce)
    .def("clear_force", &ReaK::kte::driving_actuator_2D::clearForce)
    .def("get_input_count",   &ReaK::kte::driving_actuator_2D::getInputCount)
    .def("get_input", &ReaK::kte::driving_actuator_2D::getInput)
    .def("set_input", &ReaK::kte::driving_actuator_2D::setInput)
    .add_property("drive_force", &ReaK::kte::driving_actuator_2D::DriveForce, &ReaK::kte::driving_actuator_2D::setDriveForce)
    .add_property("drive_torque", &ReaK::kte::driving_actuator_2D::DriveTorque, &ReaK::kte::driving_actuator_2D::setDriveTorque);
  
  class_< ReaK::kte::driving_actuator_3D,
          boost::noncopyable,
          bases< ReaK::kte::force_actuator_3D, ReaK::kte::system_input >,
          ReaK::shared_ptr< ReaK::kte::driving_actuator_3D >
        >("KTEDrivingActuator3D")
    .def(init<std::string>())
    .def(init<std::string, ReaK::shared_ptr< ReaK::frame_3D<double> >, ReaK::shared_ptr< ReaK::kte::reacting_kte_3D > >())
    .def("do_motion",   &ReaK::kte::driving_actuator_3D::doMotion)
    .def("do_force",    &ReaK::kte::driving_actuator_3D::doForce)
    .def("clear_force", &ReaK::kte::driving_actuator_3D::clearForce)
    .def("get_input_count",   &ReaK::kte::driving_actuator_3D::getInputCount)
    .def("get_input", &ReaK::kte::driving_actuator_3D::getInput)
    .def("set_input", &ReaK::kte::driving_actuator_3D::setInput)
    .add_property("drive_force", &ReaK::kte::driving_actuator_3D::DriveForce, &ReaK::kte::driving_actuator_3D::setDriveForce)
    .add_property("drive_torque", &ReaK::kte::driving_actuator_3D::DriveTorque, &ReaK::kte::driving_actuator_3D::setDriveTorque);
  
  
    
  class_< ReaK::kte::inertial_beam_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inertial_beam_2D>
        >("KTEInertialBeam2D")
    .def(init<std::string>())
    .def(init<std::string, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              ReaK::shared_ptr< ReaK::frame_2D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::inertial_beam_2D::doMotion)
    .def("do_force",    &ReaK::kte::inertial_beam_2D::doForce)
    .def("clear_force", &ReaK::kte::inertial_beam_2D::clearForce)
    .add_property("anchor1", &ReaK::kte::inertial_beam_2D::Anchor1, &ReaK::kte::inertial_beam_2D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::inertial_beam_2D::Anchor2, &ReaK::kte::inertial_beam_2D::setAnchor2)
    .add_property("mass", &ReaK::kte::inertial_beam_2D::Mass, &ReaK::kte::inertial_beam_2D::setMass);
  
  class_< ReaK::kte::inertial_beam_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::inertial_beam_3D >
        >("KTEInertialBeam3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              ReaK::shared_ptr< ReaK::frame_3D<double> >, 
              double>())
    .def("do_motion",   &ReaK::kte::inertial_beam_3D::doMotion)
    .def("do_force",    &ReaK::kte::inertial_beam_3D::doForce)
    .def("clear_force", &ReaK::kte::inertial_beam_3D::clearForce)
    .add_property("anchor1", &ReaK::kte::inertial_beam_3D::Anchor1, &ReaK::kte::inertial_beam_3D::setAnchor1)
    .add_property("anchor2", &ReaK::kte::inertial_beam_3D::Anchor2, &ReaK::kte::inertial_beam_3D::setAnchor2)
    .add_property("mass", &ReaK::kte::inertial_beam_3D::Mass, &ReaK::kte::inertial_beam_3D::setMass);
    
    
  class_< ReaK::kte::virtual_kte_interface_gen,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_gen>
        >("KTEVirtualInterfaceGen")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::shared_ptr< ReaK::gen_coord<double> > >())
    .def("do_motion",   &ReaK::kte::virtual_kte_interface_gen::doMotion)
    .def("do_force",    &ReaK::kte::virtual_kte_interface_gen::doForce)
    .def("clear_force", &ReaK::kte::virtual_kte_interface_gen::clearForce)
    .add_property("base_frame", &ReaK::kte::virtual_kte_interface_gen::BaseFrame, &ReaK::kte::virtual_kte_interface_gen::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::virtual_kte_interface_gen::EndFrame, &ReaK::kte::virtual_kte_interface_gen::setEndFrame);
  
  class_< ReaK::kte::virtual_kte_interface_2D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_2D>
        >("KTEVirtualInterface2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> > >())
    .def("do_motion",   &ReaK::kte::virtual_kte_interface_2D::doMotion)
    .def("do_force",    &ReaK::kte::virtual_kte_interface_2D::doForce)
    .def("clear_force", &ReaK::kte::virtual_kte_interface_2D::clearForce)
    .add_property("base_frame", &ReaK::kte::virtual_kte_interface_2D::BaseFrame, &ReaK::kte::virtual_kte_interface_2D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::virtual_kte_interface_2D::EndFrame, &ReaK::kte::virtual_kte_interface_2D::setEndFrame);
  
  class_< ReaK::kte::virtual_kte_interface_3D,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_3D >
        >("KTEVirtualInterface3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> > >())
    .def("do_motion",   &ReaK::kte::virtual_kte_interface_3D::doMotion)
    .def("do_force",    &ReaK::kte::virtual_kte_interface_3D::doForce)
    .def("clear_force", &ReaK::kte::virtual_kte_interface_3D::clearForce)
    .add_property("base_frame", &ReaK::kte::virtual_kte_interface_3D::BaseFrame, &ReaK::kte::virtual_kte_interface_3D::setBaseFrame)
    .add_property("end_frame", &ReaK::kte::virtual_kte_interface_3D::EndFrame, &ReaK::kte::virtual_kte_interface_3D::setEndFrame);
    
    
  
  class_< ReaK::kte::vmc_revolute_joint_2D,
          boost::noncopyable,
          bases< ReaK::kte::revolute_joint_2D >,
          ReaK::shared_ptr< ReaK::kte::vmc_revolute_joint_2D>
        >("KTEVMCRevoluteJoint2D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::frame_2D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_2D<double> >,
              double, double, double >())
    .def("do_motion",   &ReaK::kte::vmc_revolute_joint_2D::doMotion)
    .def("do_force",    &ReaK::kte::vmc_revolute_joint_2D::doForce)
    .def("clear_force", &ReaK::kte::vmc_revolute_joint_2D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::vmc_revolute_joint_2D::applyReactionForce)
    .add_property("stiction_coef", &ReaK::kte::vmc_revolute_joint_2D::StictionCoefficient, &ReaK::kte::vmc_revolute_joint_2D::setStictionCoefficient)
    .add_property("slip_coef", &ReaK::kte::vmc_revolute_joint_2D::SlipCoefficient, &ReaK::kte::vmc_revolute_joint_2D::setSlipCoefficient)
    .add_property("slip_velocity", &ReaK::kte::vmc_revolute_joint_2D::SlipVelocity, &ReaK::kte::vmc_revolute_joint_2D::setSlipVelocity);
  
  class_< ReaK::kte::vmc_revolute_joint_3D,
          boost::noncopyable,
          bases< ReaK::kte::revolute_joint_3D >,
          ReaK::shared_ptr< ReaK::kte::vmc_revolute_joint_3D >
        >("KTEVMCRevoluteJoint3D")
    .def(init<std::string>())
    .def(init<std::string,
              ReaK::shared_ptr< ReaK::gen_coord<double> >,
              ReaK::vect<double,3>,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::frame_3D<double> >,
              ReaK::shared_ptr< ReaK::jacobian_gen_3D<double> >,
              double, double, double >())
    .def("do_motion",   &ReaK::kte::vmc_revolute_joint_3D::doMotion)
    .def("do_force",    &ReaK::kte::vmc_revolute_joint_3D::doForce)
    .def("clear_force", &ReaK::kte::vmc_revolute_joint_3D::clearForce)
    .def("apply_reaction_force", &ReaK::kte::vmc_revolute_joint_3D::applyReactionForce)
    .add_property("stiction_coef", &ReaK::kte::vmc_revolute_joint_3D::StictionCoefficient, &ReaK::kte::vmc_revolute_joint_3D::setStictionCoefficient)
    .add_property("slip_coef", &ReaK::kte::vmc_revolute_joint_3D::SlipCoefficient, &ReaK::kte::vmc_revolute_joint_3D::setSlipCoefficient)
    .add_property("slip_velocity", &ReaK::kte::vmc_revolute_joint_3D::SlipVelocity, &ReaK::kte::vmc_revolute_joint_3D::setSlipVelocity);
   
    
    
    
  class_< py_mmc_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::mass_matrix_calc::Coords > >("MMCGenCoordVector", no_init)
    .def("__len__", &py_mmc_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::mass_matrix_calc::Coords >::size)
    .def("__getitem__", &py_mmc_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::mass_matrix_calc::Coords >::get);
  
  class_< py_mmc_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::mass_matrix_calc::Frames2D > >("MMCFrame2DVector", no_init)
    .def("__len__", &py_mmc_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::mass_matrix_calc::Frames2D >::size)
    .def("__getitem__", &py_mmc_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::mass_matrix_calc::Frames2D >::get);
  
  class_< py_mmc_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::mass_matrix_calc::Frames3D > >("MMCFrame3DVector", no_init)
    .def("__len__", &py_mmc_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::mass_matrix_calc::Frames3D >::size)
    .def("__getitem__", &py_mmc_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::mass_matrix_calc::Frames3D >::get);
  
    
  class_< ReaK::kte::mass_matrix_calc,
          boost::noncopyable,
          bases< ReaK::named_object >,
          ReaK::shared_ptr< ReaK::kte::mass_matrix_calc >
        >("KTEMassMatrixCalc")
    .def(init<std::string>())
    .def("add_inertia_gen", static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::kte::inertia_gen >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("add_inertia_2D",  static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::kte::inertia_2D >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("add_inertia_3D",  static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::kte::inertia_3D >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("add_frame_gen", static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("add_frame_2D",  static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("add_frame_3D",  static_cast< ReaK::kte::mass_matrix_calc& (ReaK::kte::mass_matrix_calc::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&) >(&ReaK::kte::mass_matrix_calc::operator<<), return_internal_reference<>())
    .def("coords", &py_mmc_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::mass_matrix_calc::Coords >::create)
    .def("frames_2D", &py_mmc_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::mass_matrix_calc::Frames2D >::create)
    .def("frames_3D", &py_mmc_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::mass_matrix_calc::Frames3D >::create);
    
    
    
  class_< py_mdl_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::manipulator_kinematics_model::Coords > >("MDLGenCoordVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::manipulator_kinematics_model::Coords >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::manipulator_kinematics_model::Coords >::get);
  
  class_< py_mdl_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::manipulator_kinematics_model::Frames2D > >("MDLFrame2DVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::manipulator_kinematics_model::Frames2D >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::manipulator_kinematics_model::Frames2D >::get);
  
  class_< py_mdl_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::manipulator_kinematics_model::Frames3D > >("MDLFrame3DVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::manipulator_kinematics_model::Frames3D >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::manipulator_kinematics_model::Frames3D >::get);
  
  class_< py_mdl_frame_vector< ReaK::kte::joint_dependent_gen_coord, &ReaK::kte::manipulator_kinematics_model::DependentCoords > >("MDLDependentGenCoordVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_gen_coord, &ReaK::kte::manipulator_kinematics_model::DependentCoords >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_gen_coord, &ReaK::kte::manipulator_kinematics_model::DependentCoords >::get);
  
  class_< py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_2D, &ReaK::kte::manipulator_kinematics_model::DependentFrames2D > >("MDLDependentFrame2DVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_2D, &ReaK::kte::manipulator_kinematics_model::DependentFrames2D >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_2D, &ReaK::kte::manipulator_kinematics_model::DependentFrames2D >::get);
  
  class_< py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_3D, &ReaK::kte::manipulator_kinematics_model::DependentFrames3D > >("MDLDependentFrame3DVector", no_init)
    .def("__len__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_3D, &ReaK::kte::manipulator_kinematics_model::DependentFrames3D >::size)
    .def("__getitem__", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_3D, &ReaK::kte::manipulator_kinematics_model::DependentFrames3D >::get);
    
  class_< ReaK::kte::manipulator_kinematics_model,
          boost::noncopyable,
          bases< ReaK::kte::kte_map >,
          ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model >
        >("KTEManipKinematics")
    .def(init<std::string>())
    .def("add_dependent_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_dependent_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_gen", static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::gen_coord<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_2D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::frame_2D<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .def("add_frame_3D",  static_cast< ReaK::kte::manipulator_kinematics_model& (ReaK::kte::manipulator_kinematics_model::*)(const ReaK::shared_ptr< ReaK::frame_3D<double> >&) >(&ReaK::kte::manipulator_kinematics_model::operator<<), return_internal_reference<>())
    .add_property("kte_model", &ReaK::kte::manipulator_kinematics_model::getModel, &ReaK::kte::manipulator_kinematics_model::setModel)
    .add_property("joint_positions", &ReaK::kte::manipulator_kinematics_model::getJointPositions, &ReaK::kte::manipulator_kinematics_model::setJointPositions)
    .add_property("joint_velocities", &ReaK::kte::manipulator_kinematics_model::getJointVelocities, &ReaK::kte::manipulator_kinematics_model::setJointVelocities)
    .add_property("joint_accelerations", &ReaK::kte::manipulator_kinematics_model::getJointAccelerations, &ReaK::kte::manipulator_kinematics_model::setJointAccelerations)
    .add_property("dependent_positions", &ReaK::kte::manipulator_kinematics_model::getDependentPositions)
    .add_property("dependent_velocities", &ReaK::kte::manipulator_kinematics_model::getDependentVelocities)
    .add_property("dependent_accelerations", &ReaK::kte::manipulator_kinematics_model::getDependentAccelerations)
    .def("coords", &py_mdl_frame_vector< ReaK::gen_coord<double>, &ReaK::kte::manipulator_kinematics_model::Coords >::create)
    .def("frames_2D", &py_mdl_frame_vector< ReaK::frame_2D<double>, &ReaK::kte::manipulator_kinematics_model::Frames2D >::create)
    .def("frames_3D", &py_mdl_frame_vector< ReaK::frame_3D<double>, &ReaK::kte::manipulator_kinematics_model::Frames3D >::create)
    .def("dependent_coords", &py_mdl_frame_vector< ReaK::kte::joint_dependent_gen_coord, &ReaK::kte::manipulator_kinematics_model::DependentCoords >::create)
    .def("dependent_frames_2D", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_2D, &ReaK::kte::manipulator_kinematics_model::DependentFrames2D >::create)
    .def("dependent_frames_3D", &py_mdl_frame_vector< ReaK::kte::joint_dependent_frame_3D, &ReaK::kte::manipulator_kinematics_model::DependentFrames3D >::create);
    
    
  class_< ReaK::kte::manipulator_dynamics_model,
          boost::noncopyable,
          bases< ReaK::kte::manipulator_kinematics_model >,
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
    .add_property("joint_states", &ReaK::kte::manipulator_dynamics_model::getJointStates, &ReaK::kte::manipulator_dynamics_model::setJointStates)
    .add_property("dependent_states", &ReaK::kte::manipulator_dynamics_model::getDependentStates)
    .def("mass_calculator", &ReaK::kte::manipulator_dynamics_model::getMassCalc, return_internal_reference<>())
    .def("compute_output", &py_dyn_mdl_compute_output)
    .def("compute_state_rate", &py_dyn_mdl_compute_state_rate);
    
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::kte_map >, 
                          ReaK::shared_ptr< ReaK::named_object > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::kte_map_chain >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::system_input >, 
                          ReaK::shared_ptr< ReaK::named_object > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::system_output >, 
                          ReaK::shared_ptr< ReaK::named_object > >();
  
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::damper_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::damper_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::damper_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::spring_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::spring_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::spring_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::torsion_damper_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::torsion_damper_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::torsion_spring_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::torsion_spring_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::flexible_beam_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::flexible_beam_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::line_point_mindist_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::line_point_mindist_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::plane_point_mindist_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::reacting_kte_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::reacting_kte_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::reacting_kte_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::revolute_joint_2D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::revolute_joint_3D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::dry_revolute_joint_2D >, 
                          ReaK::shared_ptr< ReaK::kte::revolute_joint_2D > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::dry_revolute_joint_3D >, 
                          ReaK::shared_ptr< ReaK::kte::revolute_joint_3D > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::prismatic_joint_2D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::prismatic_joint_3D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_gen > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::free_joint_2D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_2D > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::free_joint_3D >, 
                          ReaK::shared_ptr< ReaK::kte::reacting_kte_3D > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rigid_link_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rigid_link_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rigid_link_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_gen >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_gen >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_control_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_control_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_gen >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::position_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::rotation_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_gen >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::velocity_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::ang_velocity_measure_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_output > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::joint_dependent_gen_coord >, 
                          ReaK::shared_ptr< ReaK::shared_object > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_2D >, 
                          ReaK::shared_ptr< ReaK::shared_object > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::joint_dependent_frame_3D >, 
                          ReaK::shared_ptr< ReaK::shared_object > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inertia_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inertia_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inertia_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::force_actuator_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::force_actuator_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::force_actuator_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_gen >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_2D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::driving_actuator_3D >, 
                          ReaK::shared_ptr< ReaK::kte::system_input > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inertial_beam_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::inertial_beam_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_gen >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_2D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::virtual_kte_interface_3D >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::vmc_revolute_joint_2D >, 
                          ReaK::shared_ptr< ReaK::kte::revolute_joint_2D > >();
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::vmc_revolute_joint_3D >, 
                          ReaK::shared_ptr< ReaK::kte::revolute_joint_3D > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::mass_matrix_calc >, 
                          ReaK::shared_ptr< ReaK::named_object > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model >, 
                          ReaK::shared_ptr< ReaK::kte::kte_map > >();
                          
  implicitly_convertible< ReaK::shared_ptr< ReaK::kte::manipulator_dynamics_model >, 
                          ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > >();
#endif
  
};


};
















