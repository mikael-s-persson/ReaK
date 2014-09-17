
/*
 *    Copyright 2014 Sven Mikael Persson
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


#include <ReaK/ctrl/mbd_kte/inertia.hpp>
#include <ReaK/ctrl/mbd_kte/driving_actuator.hpp>
#include <ReaK/ctrl/mbd_kte/state_measures.hpp>
#include <ReaK/ctrl/mbd_kte/free_joints.hpp>
#include <ReaK/ctrl/mbd_kte/kte_map_chain.hpp>

#include <ReaK/ctrl/kte_models/manip_dynamics_model.hpp>

#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/shapes/sphere.hpp>
#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/core/recorders/data_record_po.hpp>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;




int main(int argc, char ** argv) {
  
  using namespace ReaK;
  using namespace geom;
  using namespace kte;
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description mdl_options("Model I/O options");
  mdl_options.add_options()
    ("target-model", po::value< std::string >()->default_value("models/airship3D.model.rkx"), "specify the input file for the kinematic model of the target (default is 'models/airship3D.model.rkx')")
  ;
  
  po::options_description io_options = recorder::get_data_stream_options_po_desc(true, true);
  
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(mdl_options).add(io_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::string target_mdl_fname = vm["target-model"].as<std::string>();
  
  shared_ptr< frame_3D<double> > global_base;
  shared_ptr< manipulator_kinematics_model > airship3D_kin_model;
  shared_ptr< frame_3D<double> > airship3D_grasp_frame;
  
  try {
    shared_ptr< colored_model_3D > geom_mdl;
    shared_ptr< proxy_query_model_3D > proxy_mdl;
    (*serialization::open_iarchive(target_mdl_fname))
      >> global_base
      >> airship3D_kin_model
      >> airship3D_grasp_frame
      >> geom_mdl 
      >> proxy_mdl;
  } catch(std::exception& e) {
    std::cerr << "An error occurred while trying to load the target's kinematics model! Got exception: '" << e.what() << "'." << std::endl;
    return 2;
  };
  
  
  recorder::data_stream_options data_in_opt = recorder::get_data_stream_options_from_po(vm);
  
  shared_ptr< recorder::data_extractor > data_in;
  std::vector< std::string > data_in_names;
  try {
    std::tie(data_in, data_in_names) = data_in_opt.create_extractor();
  } catch(std::exception& e) {
    std::cerr << "An error occurred while trying to open the input data stream! Got exception: '" << e.what() << "'." << std::endl;
    return 3;
  };
  
  recorder::named_value_row nvr_in  = data_in->getFreshNamedValueRow();
  
  bool is_measured_pose = true;
  try {
    nvr_in["time"]; nvr_in["p_x"]; nvr_in["p_y"]; nvr_in["p_z"]; 
    nvr_in["q_0"]; nvr_in["q_1"]; nvr_in["q_2"]; nvr_in["q_3"];
  } catch(recorder::out_of_bounds& e) { RK_UNUSED(e);
    try {
      nvr_in["time"]; nvr_in["pos_x"]; nvr_in["pos_y"]; nvr_in["pos_z"]; 
      nvr_in["q0"]; nvr_in["q1"]; nvr_in["q2"]; nvr_in["q3"];
      nvr_in["vel_x"]; nvr_in["vel_y"]; nvr_in["vel_z"]; 
      nvr_in["avel_x"]; nvr_in["avel_y"]; nvr_in["avel_z"]; 
      is_measured_pose = false;
    } catch(recorder::out_of_bounds& e) { RK_UNUSED(e);
      std::cerr << "Could not recognize the input data fields!" << std::endl;
      return 4;
    };
  };
  
  recorder::data_stream_options data_out_opt = recorder::get_data_stream_options_from_po(vm, true);
  data_out_opt.names.clear();
  if(is_measured_pose) {
    data_out_opt
      .add_name("time").add_name("p_x").add_name("p_y").add_name("p_z")
      .add_name("q_0").add_name("q_1").add_name("q_2").add_name("q_3");
  } else {
    data_out_opt
      .add_name("time").add_name("pos_x").add_name("pos_y").add_name("pos_z")
      .add_name("q0").add_name("q1").add_name("q2").add_name("q3")
      .add_name("vel_x").add_name("vel_y").add_name("vel_z")
      .add_name("avel_x").add_name("avel_y").add_name("avel_z");
  };
  shared_ptr< recorder::data_recorder > data_out = data_out_opt.create_recorder();
  
  try {
    while(true) {
      (*data_in) >> nvr_in;
      
      frame_3D<double> cur_pose;
      if(is_measured_pose) {
        cur_pose.Position = vect<double,3>(nvr_in["p_x"], nvr_in["p_y"], nvr_in["p_z"]);
        cur_pose.Quat     = quaternion<double>(vect<double,4>(nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"]));
      } else {
        cur_pose.Position = vect<double,3>(nvr_in["pos_x"], nvr_in["pos_y"], nvr_in["pos_z"]);
        cur_pose.Quat     = quaternion<double>(vect<double,4>(nvr_in["q0"], nvr_in["q1"], nvr_in["q2"], nvr_in["q3"]));
        cur_pose.Velocity    = vect<double,3>(nvr_in["vel_x"], nvr_in["vel_y"], nvr_in["vel_z"]);
        cur_pose.AngVelocity = vect<double,3>(nvr_in["avel_x"], nvr_in["avel_y"], nvr_in["avel_z"]);
      };
      
      *(airship3D_kin_model->getFrame3D(0)) = cur_pose; 
      airship3D_kin_model->doDirectMotion();
      frame_3D<double> cur_grapple = airship3D_grasp_frame->getGlobalFrame();
      
      (*data_out) << nvr_in["time"] << cur_grapple.Position 
                  << cur_grapple.Quat[0] << cur_grapple.Quat[1] << cur_grapple.Quat[2] << cur_grapple.Quat[3];
      if(!is_measured_pose)
        (*data_out) << cur_grapple.Velocity << cur_grapple.AngVelocity;
      (*data_out) << recorder::data_recorder::end_value_row;
      
    };
  } catch(recorder::end_of_record& e) { RK_UNUSED(e); };
  
  (*data_out) << recorder::data_recorder::flush;
  
  return 0;
};






