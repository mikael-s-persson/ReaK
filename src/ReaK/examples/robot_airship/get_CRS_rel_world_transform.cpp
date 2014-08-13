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

#include <ReaK/core/kinetostatics/calibrate_frames_3D.hpp>

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
  
  std::vector< std::pair< vect<double,3>, vect<double,3> > > calib_pts;
  try {
    while(true) {
      (*data_in) >> nvr_in;
      std::pair< vect<double,3>, vect<double,3> > cur_pt;
      cur_pt.second = vect<double,3>(nvr_in["wo_x"], nvr_in["wo_y"], nvr_in["wo_z"]);
      vect_n<double> CRS_jt(nvr_in["track"], nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"], nvr_in["q_4"], nvr_in["q_5"]);
      kte_model->setJointPositions(CRS_jt);
      kte_model->doDirectMotion();
      shared_ptr< frame_3D<double> > EE_fr = kte_model->getDependentFrame3D(0)->mFrame;
      EE_to_marker.Parent = EE_fr;
      cur_pt.first = EE_to_marker.getGlobalPose().Position;
      calib_pts.push_back(cur_pt);
      std::cout << " CRS = " << cur_pt.first << " and World = " << cur_pt.second << std::endl;
    };
  } catch(end_of_record& e) { 
    RK_UNUSED(e);
  };
  
  pose_3D<double> CRS_wo = get_relative_pose_pointcloud(calib_pts);
  
  std::cout << "Rt = \n" << CRS_wo.Quat.getRotMat() << std::endl;
  std::cout << "CRS_wo = \n" << CRS_wo << std::endl;
  
  double avg_err = 0.0;
  for(std::size_t i = 0; i < calib_pts.size(); ++i) {
    vect<double,3> CRS_new_pt = CRS_wo.transformToParent(calib_pts[i].second);
    double cur_err = norm_2(calib_pts[i].first - CRS_new_pt);
    avg_err += cur_err;
    std::cout << " CRS = " << calib_pts[i].first << " New CRS = " << CRS_new_pt 
              << " ErrNorm = " << cur_err << std::endl;
  };
  avg_err /= calib_pts.size();
  std::cout << " Total Average Error Norm = " << avg_err << std::endl;
  
  return 0;
};







