/**
 * \file build_P3R3R_model.cpp
 *
 * This application constructs the KTE-based kinematics models for a P3R3R manipulator.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */


#include <ReaK/ctrl/kte_models/manip_P3R3R_arm.hpp>

#include <ReaK/ctrl/topologies/joint_space_limits.hpp>

#include <ReaK/core/recorders/data_record_po.hpp>

#include <ReaK/core/lin_alg/mat_alg.hpp>
#include <ReaK/core/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/core/lin_alg/mat_svd_method.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;


int main(int argc, char ** argv) {
  
  using namespace ReaK;
  using namespace kte;
  using namespace recorder;
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message. This program expects, as input, a data file with robot configurations and corresponding world coordinates.")
  ;
  
  po::options_description io_options = get_data_stream_options_po_desc(true);
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  data_stream_options data_in_opt = get_data_stream_options_from_po(vm);
  
  shared_ptr< data_extractor > data_in;
  std::vector< std::string > data_in_names;
  std::tie(data_in, data_in_names) = data_in_opt.create_extractor();
  
  named_value_row nvr_in  = data_in->getFreshNamedValueRow();
  
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
  
  pose_3D<double> EE_to_marker(weak_ptr< pose_3D<double> >(), vect<double,3>(-0.033, 0.0, 0.107), quaternion<double>());
  
  std::vector< vect<double,3> > CRS_pts;
  std::vector< vect<double,3> > world_pts;
  try {
    while(true) {
      (*data_in) >> nvr_in;
      world_pts.push_back(vect<double,3>(nvr_in["wo_x"], nvr_in["wo_y"], nvr_in["wo_z"]));
      vect_n<double> CRS_jt(nvr_in["track"], nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"], nvr_in["q_4"], nvr_in["q_5"]);
      kte_model->setJointPositions(CRS_jt);
      kte_model->doDirectMotion();
      shared_ptr< frame_3D<double> > EE_fr = kte_model->getDependentFrame3D(0)->mFrame;
      EE_to_marker.Parent = EE_fr;
      CRS_pts.push_back(EE_to_marker.getGlobalPose().Position);
      std::cout << " CRS = " << CRS_pts.back() << " and World = " << world_pts.back() << std::endl;
    };
  } catch(end_of_record& e) { 
    RK_UNUSED(e);
  };
  
  
  // NOTE: This is a bit of a brute-force method, but it's not expected that there are too many points.
  
  std::vector< vect<double,3> > CRS_vects;
  std::vector< vect<double,3> > world_vects;
  for(std::size_t i = 0; i < CRS_pts.size(); ++i) {
    for(std::size_t j = i+1; j < CRS_pts.size(); ++j) {
      CRS_vects.push_back(CRS_pts[j] - CRS_pts[i]);
      world_vects.push_back(world_pts[j] - world_pts[i]);
    };
  };
  
  mat<double,mat_structure::rectangular> X(CRS_vects.size(),3);
  mat<double,mat_structure::rectangular> B(world_vects.size(),3);
  for(std::size_t i = 0; i < CRS_vects.size(); ++i) {
    X(i, 0) = CRS_vects[i][0]; X(i, 1) = CRS_vects[i][1]; X(i, 2) = CRS_vects[i][2];
    B(i, 0) = world_vects[i][0]; B(i, 1) = world_vects[i][1]; B(i, 2) = world_vects[i][2];
  };
  
  // Solve the Orthogonal Procrustes problem: (Kabsch's algorithm)
  mat<double,mat_structure::square> C( transpose_view(X) * B );
  
  mat<double,mat_structure::square> U(3);
  mat<double,mat_structure::diagonal> E(3);
  mat<double,mat_structure::square> V(3);

  decompose_SVD(C,U,E,V);
  
  mat<double,mat_structure::square> Rt(U * transpose_view(V));
  
  std::cout << "Rt = \n" << Rt << std::endl;
  
  return 0;
};







