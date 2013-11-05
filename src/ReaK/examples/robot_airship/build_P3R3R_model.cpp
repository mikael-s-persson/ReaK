/**
 * \file build_P3R3R_model.cpp
 *
 * This application constructs the KTE-based kinematics models for a P3R3R manipulator.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */


#include "kte_models/manip_P3R3R_arm.hpp"
#include "shapes/kte_chain_geometry.hpp"

#include "shapes/sphere.hpp"
#include "shapes/box.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "shapes/capped_cylinder.hpp"

#include "topologies/joint_space_limits.hpp"

#include "topologies/manip_P3R3R_workspaces.hpp"
#include "interpolation/linear_interp.hpp"

#include "serialization/archiver_factory.hpp"

#include <boost/tuple/tuple.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


int main(int argc, char ** argv) {
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("output-path,p", po::value< std::string >()->default_value("models"), "specify the output path (default is models)")
    ("output-name,o", po::value< std::string >()->default_value("CRS_A465"), "specify the output base-name (default is CRS_A465)")
    ("format", po::value< std::string >()->default_value("xml"), "specify the format that should be outputted (default is xml, but can also be bin or protobuf)")
  ;
  
  po::options_description output_options("Output options");
  output_options.add_options()
    ("geometry", "specify that the geometry should be outputted (default is not)")
    ("proxy", "specify that the proximity-query model (simplified bounding-geometry) should be outputted (default is not)")
    ("kte-model", "specify that the KTE kinematics model should be outputted (with limits) (default is not)")
    ("all-assembled", "specify that the KTE kinematics + geometry + proximity-query models should be outputted (with limits) together in one file (default is not)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(output_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help") || (vm.count("geometry") + vm.count("proxy") + vm.count("kte-model") + vm.count("all-assembled") < 1)) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::string output_base_name = vm["output-name"].as<std::string>();
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  std::string output_extension = ".rkx";
  if(vm["format"].as<std::string>() == "bin") {
    output_extension = ".rkb";
  } else if(vm["format"].as<std::string>() == "protobuf") {
    output_extension = ".pbuf";
  };
  
  
  
  
  using namespace ReaK;
  using namespace geom;
  using namespace kte;
  
  
  shared_ptr< kte_chain_geometry_3D > kte_geom_model;
  
  if(vm.count("geometry") + vm.count("proxy") + vm.count("all-assembled") > 0) {
    
    // CRS_A465 geometries:
    
    shared_ptr< coord_arrows_3D > robot_base_arrows( new coord_arrows_3D("robot_base_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > track_joint_arrows(new coord_arrows_3D("track_joint_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_1_arrows(new coord_arrows_3D("arm_joint_1_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_2_arrows(new coord_arrows_3D("arm_joint_2_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_3_arrows(new coord_arrows_3D("arm_joint_3_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_4_arrows(new coord_arrows_3D("arm_joint_4_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_5_arrows(new coord_arrows_3D("arm_joint_5_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    shared_ptr< coord_arrows_3D > arm_joint_6_arrows(new coord_arrows_3D("arm_joint_6_arrows", shared_ptr< frame_3D<double> >(), pose_3D<double>(),0.3));
    
    shared_ptr< capped_cylinder > link1_cyl(new capped_cylinder("link1_cyl", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.15),quaternion<double>()), 
      0.3, 0.1));
    shared_ptr< capped_cylinder > joint2_cyl(new capped_cylinder("joint2_cyl", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.3302),axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0)).getQuaternion()), 
      0.34, 0.09));
    shared_ptr< capped_cylinder > link2_cyl(new capped_cylinder("link2_cyl", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.15),quaternion<double>()), 
      0.3, 0.07));
    shared_ptr< capped_cylinder > link3_cyl(new capped_cylinder("link3_cyl", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.165),quaternion<double>()), 
      0.33, 0.07));
    shared_ptr< capped_cylinder > link5_cyl(new capped_cylinder("link5_cyl", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.0381),quaternion<double>()), 
      0.0762, 0.05));
    shared_ptr< sphere > EE_sphere(new sphere("EE_sphere", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.04,0.0,0.05),quaternion<double>()), 
      0.11));
    
    shared_ptr< box > EE_bumblebee(new box("EE_bumblebee", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.105,0.0,0.0483),quaternion<double>()), 
      vect<double,3>(0.04,0.14,0.04)));
    shared_ptr< box > EE_bumblebee_support(new box("EE_bumblebee_support", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.0525,0.0,0.0393),quaternion<double>()), 
      vect<double,3>(0.065,0.06,0.03)));
    shared_ptr< box > EE_gripper_box(new box("EE_gripper_box", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.0669),quaternion<double>()), 
      vect<double,3>(0.04,0.08,0.1338)));
    shared_ptr< box > EE_gripper_fingers(new box("EE_gripper_fingers", shared_ptr< frame_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.1588),quaternion<double>()), 
      vect<double,3>(0.02,0.01,0.05)));
    
    kte_geom_model = shared_ptr< kte_chain_geometry_3D >(new kte_chain_geometry_3D("CRS_A465_kte"));
    
    if(vm.count("geometry") || vm.count("all-assembled")) {
      (*kte_geom_model)
      .addElement("manip_P3R3R_track_base",color(0,0,0),robot_base_arrows)
      .addElement("manip_P3R3R_track_end",color(0,0,0),track_joint_arrows)
      .addElement("manip_3R3R_joint_1_end",color(0,0,0),arm_joint_1_arrows)
      .addElement("manip_3R3R_joint_2_end",color(0,0,0),arm_joint_2_arrows)
      .addElement("manip_3R3R_joint_3_end",color(0,0,0),arm_joint_3_arrows)
      .addElement("manip_3R3R_joint_4_end",color(0,0,0),arm_joint_4_arrows)
      .addElement("manip_3R3R_joint_5_end",color(0,0,0),arm_joint_5_arrows)
      .addElement("manip_3R3R_joint_6_end",color(0,0,0),arm_joint_6_arrows)
      .addElement("manip_P3R3R_track_end",color(0.2,0.2,0.2),link1_cyl)
      .addElement("manip_3R3R_joint_1_end",color(0.3,0.3,0.3),joint2_cyl)
      .addElement("manip_3R3R_joint_2_end",color(0.7,0.7,0.5),link2_cyl)
      .addElement("manip_3R3R_joint_3_end",color(0.7,0.7,0.5),link3_cyl)
      .addElement("manip_3R3R_joint_5_end",color(0.1,0.1,0.1),link5_cyl)
      //.addElement("manip_3R3R_joint_6_end",color(0.6,0.6,0.0),EE_sphere)
      .addElement("manip_3R3R_joint_6_end",color(0.6,0.6,0.0),EE_bumblebee)
      .addElement("manip_3R3R_joint_6_end",color(0.5,0.5,0.5),EE_bumblebee_support)
      .addElement("manip_3R3R_joint_6_end",color(0.1,0.1,0.1),EE_gripper_box)
      .addElement("manip_3R3R_joint_6_end",color(0.5,0.5,0.5),EE_gripper_fingers);
    };
    
    if(vm.count("proxy") || vm.count("all-assembled")) {
      (*kte_geom_model)
      .addShape("manip_3R3R_joint_1_end",joint2_cyl)
      .addShape("manip_3R3R_joint_2_end",link2_cyl)
      .addShape("manip_3R3R_joint_3_end",link3_cyl)
      .addShape("manip_3R3R_joint_5_end",link5_cyl)
      .addShape("manip_3R3R_joint_6_end",EE_sphere);
    };
  };
  
  // Save the kte_geom_model object:
  if(vm.count("geometry") || vm.count("proxy")) {
    (*serialization::open_oarchive(output_path_name + "/" + output_base_name + ".geom" + output_extension))
      << kte_geom_model;
  };
  
  
  // CRS_A465 model parameters:
  
  if(vm.count("kte-model") + vm.count("all-assembled") > 0) {
    
    shared_ptr< frame_3D<double> > base_frame(new frame_3D<double>());
    
    shared_ptr< manip_P3R3R_kinematics > kte_model(new manip_P3R3R_kinematics("CRS_A465_kte_model",
      shared_ptr< frame_3D<double> >(new frame_3D<double>(base_frame,
        vect<double,3>(0.0,-3.3,0.3),  // aGlobalToBasePlate
        axis_angle<double>(M_PI * 0.5, vect<double,3>(0.0,0.0,1.0)).getQuaternion(), // align the x-axis along the track.
        vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0)
      )),
      0.3302,  // aBaseToShoulder 
      0.3048,  // aShoulderToElbowDist
      0.1500,  // aElbowToJoint4
      0.1802,  // aJoint4ToWrist
      0.0762,  // aWristToFlange
      vect_n<double>(    // aJointLowerBounds
        0.0,
        -3.05432619099,
        -1.57079632679,
        -1.91986217719,
        -3.14159265359,
        -1.83259571459,
        -3.14159265359),
      vect_n<double>(    // aJointUpperBounds
        3.0,
        3.05432619099,
        1.57079632679,
        1.91986217719,
        3.14159265359,
        1.83259571459,
        3.14159265359)));
    
    shared_ptr< pp::joint_limits_collection<double> > joint_rate_limits(new pp::joint_limits_collection<double>("CRS_A465_joint_limits"));

    joint_rate_limits->gen_speed_limits.resize(7);
    joint_rate_limits->gen_speed_limits[0] = 0.8;
    joint_rate_limits->gen_speed_limits[1] = 3.14159265359;
    joint_rate_limits->gen_speed_limits[2] = 3.14159265359;
    joint_rate_limits->gen_speed_limits[3] = 3.14159265359;
    joint_rate_limits->gen_speed_limits[4] = 2.98451302091;
    joint_rate_limits->gen_speed_limits[5] = 3.01941960595;
    joint_rate_limits->gen_speed_limits[6] = 2.98451302091;
    
    joint_rate_limits->gen_accel_limits.resize(7);
    joint_rate_limits->gen_accel_limits[0] = 3.0;
    joint_rate_limits->gen_accel_limits[1] = 12.5663706144;
    joint_rate_limits->gen_accel_limits[2] = 12.5663706144;
    joint_rate_limits->gen_accel_limits[3] = 12.5663706144;
    joint_rate_limits->gen_accel_limits[4] = 24.9582083035;
    joint_rate_limits->gen_accel_limits[5] = 24.9582083035;
    joint_rate_limits->gen_accel_limits[6] = 24.9582083035;
    
    joint_rate_limits->gen_jerk_limits.resize(7);
    joint_rate_limits->gen_jerk_limits[0] = 12.0;
    joint_rate_limits->gen_jerk_limits[1] = 125.663706144;
    joint_rate_limits->gen_jerk_limits[2] = 125.663706144;
    joint_rate_limits->gen_jerk_limits[3] = 125.663706144;
    joint_rate_limits->gen_jerk_limits[4] = 249.582083035;
    joint_rate_limits->gen_jerk_limits[5] = 249.582083035;
    joint_rate_limits->gen_jerk_limits[6] = 249.582083035;
    
    joint_rate_limits->frame2D_speed_limits.resize(0);
    joint_rate_limits->frame2D_accel_limits.resize(0);
    joint_rate_limits->frame2D_jerk_limits.resize(0);
    joint_rate_limits->frame3D_speed_limits.resize(0);
    joint_rate_limits->frame3D_accel_limits.resize(0);
    joint_rate_limits->frame3D_jerk_limits.resize(0);
    
    
    shared_ptr< pp::manip_static_workspace< manip_P3R3R_kinematics, 0 >::rl_workspace_type >
      workspace = pp::make_manip_static_workspace< 0 >(
        pp::linear_interpolation_tag(), kte_model, joint_rate_limits, 0.1);
    
    
    
    
    // Save kte_model and joint_rate_limits together in one file (kte_model first).
    if(vm.count("kte-model")) {
      (*serialization::open_oarchive(output_path_name + "/" + output_base_name + ".kte_ik" + output_extension))
        << base_frame << kte_model << joint_rate_limits;
    };
    
    if(vm.count("all-assembled")) {
      shared_ptr< colored_model_3D > mdl_geom;
      shared_ptr< proxy_query_model_3D > mdl_prox;
      boost::tie(mdl_geom, mdl_prox) = kte_geom_model->attachToKTEChain(*(kte_model->getKTEChain()));
      
      (*serialization::open_oarchive(output_path_name + "/" + output_base_name + ".model" + output_extension))
        << base_frame << kte_model << joint_rate_limits << mdl_geom << mdl_prox;
    };
    
    
  };
  
  
  return 0;
};







