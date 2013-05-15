/**
 * \file test_lab_scene.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */


#include "CRS_A465_geom_model.hpp"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/sensors/SoIdleSensor.h>

#include "shapes/oi_scene_graph.hpp"

#include "shapes/plane.hpp"
#include "shapes/box.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "shapes/capped_cylinder.hpp"

#include "proximity/proxy_query_model.hpp"

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/state_measures.hpp"
#include "mbd_kte/free_joints.hpp"

#include "mbd_kte/spring.hpp"
#include "mbd_kte/damper.hpp"
#include "mbd_kte/free_joints.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"
#include "optimization/optim_exceptions.hpp"

#include "interpolation/cubic_hermite_interp.hpp"

#include "topologies/manip_free_workspace.hpp"
#include "path_planning/rrt_path_planner.hpp"
#include "path_planning/rrtstar_path_planner.hpp"
#include "topologies/direct_kinematics_topomap.hpp"
#include "topologies/inverse_kinematics_topomap.hpp"
#include "path_planning/frame_tracer_coin3d.hpp"


struct all_robot_info {
  SoSeparator* sg_root;
  SoSeparator* pp_root;
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > kin_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > robot_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > lab_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_lab_proxy;
  ReaK::shared_ptr< ReaK::frame_3D<double> > airship_frame;
  ReaK::pose_3D<double> target_frame;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > airship_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > lab_airship_proxy;
  std::vector< SoSeparator* > solution_seps;
};

void add_new_solution_sep(void* userData, SoSensor*) {
  all_robot_info* r_info = reinterpret_cast< all_robot_info* >(userData);
  
  while(!r_info->solution_seps.empty()) {
    r_info->pp_root->addChild(r_info->solution_seps.back());
    r_info->solution_seps.back()->unref();
    r_info->solution_seps.pop_back();
  };
  
};


void keyboard_press_hdl(void* userData, SoEventCallback* eventCB) {
  all_robot_info* r_info = reinterpret_cast< all_robot_info* >(userData);
  const SoEvent* event = eventCB->getEvent();
  
  static bool IK_enabled = false;
  static bool PP_enabled = false;
  
  if( SO_KEY_PRESS_EVENT(event, Q) ) {
    r_info->builder.track_joint_coord->q += 0.01;
    if(r_info->builder.track_joint_coord->q > r_info->builder.joint_upper_bounds[0])
      r_info->builder.track_joint_coord->q = r_info->builder.joint_upper_bounds[0];
  } else if( SO_KEY_PRESS_EVENT(event, A) ) {
    r_info->builder.track_joint_coord->q -= 0.01;
    if(r_info->builder.track_joint_coord->q < r_info->builder.joint_lower_bounds[0])
      r_info->builder.track_joint_coord->q = r_info->builder.joint_lower_bounds[0];
  } else if( SO_KEY_PRESS_EVENT(event, W) ) {
    r_info->builder.arm_joint_1_coord->q += M_PI * 0.01;
    if(r_info->builder.arm_joint_1_coord->q > r_info->builder.joint_upper_bounds[1])
      r_info->builder.arm_joint_1_coord->q = r_info->builder.joint_upper_bounds[1];
  } else if( SO_KEY_PRESS_EVENT(event, S) ) {
    r_info->builder.arm_joint_1_coord->q -= M_PI * 0.01;
    if(r_info->builder.arm_joint_1_coord->q < r_info->builder.joint_lower_bounds[1])
      r_info->builder.arm_joint_1_coord->q = r_info->builder.joint_lower_bounds[1];
  } else if( SO_KEY_PRESS_EVENT(event, E) ) {
    r_info->builder.arm_joint_2_coord->q += M_PI * 0.01;
    if(r_info->builder.arm_joint_2_coord->q > r_info->builder.joint_upper_bounds[2])
      r_info->builder.arm_joint_2_coord->q = r_info->builder.joint_upper_bounds[2];
  } else if( SO_KEY_PRESS_EVENT(event, D) ) {
    r_info->builder.arm_joint_2_coord->q -= M_PI * 0.01;
    if(r_info->builder.arm_joint_2_coord->q < r_info->builder.joint_lower_bounds[2])
      r_info->builder.arm_joint_2_coord->q = r_info->builder.joint_lower_bounds[2];
  } else if( SO_KEY_PRESS_EVENT(event, R) ) {
    r_info->builder.arm_joint_3_coord->q += M_PI * 0.01;
    if(r_info->builder.arm_joint_3_coord->q > r_info->builder.joint_upper_bounds[3])
      r_info->builder.arm_joint_3_coord->q = r_info->builder.joint_upper_bounds[3];
  } else if( SO_KEY_PRESS_EVENT(event, F) ) {
    r_info->builder.arm_joint_3_coord->q -= M_PI * 0.01;
    if(r_info->builder.arm_joint_3_coord->q < r_info->builder.joint_lower_bounds[3])
      r_info->builder.arm_joint_3_coord->q = r_info->builder.joint_lower_bounds[3];
  } else if( SO_KEY_PRESS_EVENT(event, T) ) {
    r_info->builder.arm_joint_4_coord->q += M_PI * 0.01; 
    if(r_info->builder.arm_joint_4_coord->q > r_info->builder.joint_upper_bounds[4])
      r_info->builder.arm_joint_4_coord->q = r_info->builder.joint_upper_bounds[4];
  } else if( SO_KEY_PRESS_EVENT(event, G) ) {
    r_info->builder.arm_joint_4_coord->q -= M_PI * 0.01; 
    if(r_info->builder.arm_joint_4_coord->q < r_info->builder.joint_lower_bounds[4])
      r_info->builder.arm_joint_4_coord->q = r_info->builder.joint_lower_bounds[4];
  } else if( SO_KEY_PRESS_EVENT(event, Y) ) {
    r_info->builder.arm_joint_5_coord->q += M_PI * 0.01; 
    if(r_info->builder.arm_joint_5_coord->q > r_info->builder.joint_upper_bounds[5])
      r_info->builder.arm_joint_5_coord->q = r_info->builder.joint_upper_bounds[5];
  } else if( SO_KEY_PRESS_EVENT(event, H) ) {
    r_info->builder.arm_joint_5_coord->q -= M_PI * 0.01; 
    if(r_info->builder.arm_joint_5_coord->q < r_info->builder.joint_lower_bounds[5])
      r_info->builder.arm_joint_5_coord->q = r_info->builder.joint_lower_bounds[5];
  } else if( SO_KEY_PRESS_EVENT(event, U) ) {
    r_info->builder.arm_joint_6_coord->q += M_PI * 0.01;
    if(r_info->builder.arm_joint_6_coord->q > r_info->builder.joint_upper_bounds[6])
      r_info->builder.arm_joint_6_coord->q = r_info->builder.joint_upper_bounds[6];
  } else if( SO_KEY_PRESS_EVENT(event, J) ) {
    r_info->builder.arm_joint_6_coord->q -= M_PI * 0.01;
    if(r_info->builder.arm_joint_6_coord->q < r_info->builder.joint_lower_bounds[6])
      r_info->builder.arm_joint_6_coord->q = r_info->builder.joint_lower_bounds[6];
  } else if( SO_KEY_PRESS_EVENT(event, Z) ) {
    r_info->airship_frame->Position[2] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, X) ) {
    r_info->airship_frame->Position[2] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, UP_ARROW) ) {
    r_info->airship_frame->Position[1] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, DOWN_ARROW) ) {
    r_info->airship_frame->Position[1] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, LEFT_ARROW) ) {
    r_info->airship_frame->Position[0] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, RIGHT_ARROW) ) {
    r_info->airship_frame->Position[0] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, B) ) {
    r_info->airship_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(1.0,0.0,0.0));
  } else if( SO_KEY_PRESS_EVENT(event, N) ) {
    r_info->airship_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(0.0,1.0,0.0));
  } else if( SO_KEY_PRESS_EVENT(event, M) ) {
    r_info->airship_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(0.0,0.0,1.0));
  } else if( SO_KEY_PRESS_EVENT(event, P) ) {
    IK_enabled = !IK_enabled;
  } else if( SO_KEY_PRESS_EVENT(event, L) ) {
    PP_enabled = true;
  };
  
  r_info->airship_chain->doMotion();
  
  if(IK_enabled) {
    try {
      ReaK::vect_n<double> jt_sol = r_info->builder.compute_inverse_kinematics(r_info->target_frame.getGlobalPose());
      r_info->builder.track_joint_coord->q = jt_sol[0];
      r_info->builder.arm_joint_1_coord->q = jt_sol[1];
      r_info->builder.arm_joint_2_coord->q = jt_sol[2];
      r_info->builder.arm_joint_3_coord->q = jt_sol[3];
      r_info->builder.arm_joint_4_coord->q = jt_sol[4];
      r_info->builder.arm_joint_5_coord->q = jt_sol[5];
      r_info->builder.arm_joint_6_coord->q = jt_sol[6];
    } catch( ReaK::optim::infeasible_problem& e ) { RK_UNUSED(e); };
  };
  
  if(PP_enabled) {
  try {
    ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > manip_kin_mdl = r_info->builder.get_manipulator_kin_model();
    ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> > manip_jt_limits(&(r_info->builder.joint_rate_limits), ReaK::null_deleter());
    typedef ReaK::pp::manip_quasi_static_env<ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type, ReaK::pp::cubic_hermite_interpolation_tag> workspace_type;
    ReaK::shared_ptr<workspace_type> 
      workspace(new workspace_type(
        r_info->builder.get_rl_joint_space_1st(),
        manip_kin_mdl,
        manip_jt_limits,
        0.1, 1.2));
    /*typedef ReaK::pp::manip_quasi_static_env<ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_type, ReaK::pp::sap_interpolation_tag> workspace_type;
    ReaK::shared_ptr<workspace_type> 
      workspace(new workspace_type(
        r_info->builder.get_rl_joint_space(),
        manip_kin_mdl,
        manip_jt_limits,
        0.1, 1.0, 1e-6, 60));*/
    
    (*workspace) << r_info->robot_lab_proxy << r_info->robot_airship_proxy;
    
    ReaK::shared_ptr< ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type > jt_space(new ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type(r_info->builder.get_rl_joint_space_1st()));
    ReaK::shared_ptr< ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type > normal_jt_space(new ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type(r_info->builder.get_joint_space_1st()));
    
    
    ReaK::pp::topology_traits<ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type>::point_type start_point, goal_point;
    ReaK::pp::topology_traits<ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type>::point_type start_inter, goal_inter;
    
    ReaK::pp::detail::read_joint_coordinates_impl(start_inter, r_info->builder.get_joint_space(), manip_kin_mdl);
    
    get<0>(start_inter) = ReaK::vect<double,7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    get<1>(start_inter) = ReaK::vect<double,7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
//     get<2>(start_inter) = ReaK::vect<double,7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    
    start_point = manip_jt_limits->map_to_space(start_inter, r_info->builder.get_joint_space_1st(), r_info->builder.get_rl_joint_space_1st());
    
    ReaK::vect_n<double> jt_desired = r_info->builder.compute_inverse_kinematics(r_info->target_frame.getGlobalPose());
    get<0>(goal_inter) = ReaK::vect<double,7>(jt_desired);
    get<1>(goal_inter) = ReaK::vect<double,7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
//     get<2>(goal_inter) = ReaK::vect<double,7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    
    goal_point = manip_jt_limits->map_to_space(goal_inter, r_info->builder.get_joint_space_1st(), r_info->builder.get_rl_joint_space_1st());
    
    typedef ReaK::pp::frame_tracer_3D<
      ReaK::pp::manip_rl_direct_kin_map< ReaK::pp::joint_limits_collection<double>, ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type >,
      ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type, 
      ReaK::pp::identity_topo_map, 
      ReaK::pp::print_sbmp_progress<> > frame_reporter_type;
    frame_reporter_type temp_reporter(
      ReaK::pp::manip_rl_direct_kin_map< ReaK::pp::joint_limits_collection<double>, ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type >(
        manip_kin_mdl, manip_jt_limits, normal_jt_space),
      jt_space,
      0.05);
    temp_reporter.add_traced_frame(r_info->builder.arm_joint_6_end);
    
    ReaK::pp::rrtstar_path_planner<workspace_type, frame_reporter_type>
      workspace_planner(
        workspace,
        start_point,
        goal_point,
        3000,
        1000,
        ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
        ReaK::pp::UNIDIRECTIONAL_PLANNING,
        temp_reporter,
        100);
    
    workspace_planner.solve_path();
    
    r_info->solution_seps.push_back(workspace_planner.get_reporter().get_motion_graph_tracer(r_info->builder.arm_joint_6_end).get_separator());
    r_info->solution_seps.back()->ref();
    if(workspace_planner.get_reporter().get_solution_count()) {
      r_info->solution_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info->builder.arm_joint_6_end).get_separator());
      r_info->solution_seps.back()->ref();
    };
    SoIdleSensor* ask_for_new_solsep = new SoIdleSensor(add_new_solution_sep, r_info);
    ask_for_new_solsep->schedule();
    
    r_info->builder.track_joint_coord->q = get<0>(start_inter)[0];
    r_info->builder.arm_joint_1_coord->q = get<0>(start_inter)[1];
    r_info->builder.arm_joint_2_coord->q = get<0>(start_inter)[2];
    r_info->builder.arm_joint_3_coord->q = get<0>(start_inter)[3];
    r_info->builder.arm_joint_4_coord->q = get<0>(start_inter)[4];
    r_info->builder.arm_joint_5_coord->q = get<0>(start_inter)[5];
    r_info->builder.arm_joint_6_coord->q = get<0>(start_inter)[6];
    
  } catch( ReaK::optim::infeasible_problem& e ) { RK_UNUSED(e);
    
  };
  PP_enabled = false;
  };
  
  r_info->kin_chain->doMotion();
  
};


int main(int argc, char ** argv) {
  using namespace ReaK;
  
  all_robot_info r_info;
  
  //r_info.builder.create_geom_from_preset();
  r_info.builder.load_kte_and_geom("models/CRS_A465_with_geom.xml");
  r_info.builder.load_limits_from_file("models/CRS_A465_limits.xml");
  
  r_info.robot_proxy = r_info.builder.get_proximity_model();
  r_info.kin_chain = r_info.builder.get_kinematics_kte_chain();
  
  
  shared_ptr< geom::colored_model_3D > lab_geom_model;
  {
    serialization::xml_iarchive in("models/MD148_lab_model.xml");
    in >> lab_geom_model >> r_info.lab_proxy;
  };
  
  shared_ptr< geom::colored_model_3D > airship_geom_model;
  
  {
    
    shared_ptr< kte::position_measure_3D > airship3D_position;
    shared_ptr< kte::rotation_measure_3D > airship3D_rotation;
    shared_ptr< kte::driving_actuator_3D > airship3D_actuator;
    shared_ptr< kte::inertia_3D > airship3D_inertia;
    shared_ptr< kte::mass_matrix_calc > airship3D_mass_calc;
    
    serialization::xml_iarchive in("models/airship3D_with_geom.xml");
    in >> r_info.airship_frame
       >> airship3D_position
       >> airship3D_rotation
       >> airship3D_actuator
       >> airship3D_inertia
       >> r_info.airship_chain
       >> airship3D_mass_calc
       >> airship_geom_model 
       >> r_info.airship_proxy;
  };
  
  r_info.target_frame = pose_3D<double>(r_info.airship_frame,
      vect<double,3>(0.97 * std::sin(0.2 / 0.93),0.0,0.97 * std::cos(0.2 / 0.93)),
      axis_angle<double>(0.2 / 0.93 / 2.0,vect<double,3>(0.0,1.0,0.0)).getQuaternion()
      * quaternion<double>::yrot(M_PI) * quaternion<double>::zrot(0.5 * M_PI));
  r_info.target_frame.Position += r_info.target_frame.Quat * (-0.3 * vect_k);
  
  r_info.robot_lab_proxy     = ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D >(new ReaK::geom::proxy_query_pair_3D("robot_lab_proxy",r_info.robot_proxy, r_info.lab_proxy));
  r_info.robot_airship_proxy = ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D >(new ReaK::geom::proxy_query_pair_3D("robot_airship_proxy",r_info.robot_proxy, r_info.airship_proxy));
  r_info.lab_airship_proxy   = ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D >(new ReaK::geom::proxy_query_pair_3D("lab_airship_proxy",r_info.lab_proxy, r_info.airship_proxy));
  
  r_info.robot_lab_proxy->setModelPair(r_info.robot_proxy, r_info.lab_proxy);
  r_info.robot_airship_proxy->setModelPair(r_info.robot_proxy, r_info.airship_proxy);
  r_info.lab_airship_proxy->setModelPair(r_info.lab_proxy, r_info.airship_proxy);
  
  {
  QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
  
  
  {
    
    r_info.builder.track_joint_coord->q = 0.2;
    r_info.builder.arm_joint_1_coord->q = M_PI * 0.25;
    r_info.builder.arm_joint_2_coord->q = -M_PI * 0.125;
    r_info.builder.arm_joint_3_coord->q = -M_PI * 0.375;
    r_info.builder.arm_joint_4_coord->q = M_PI * 0.125; 
    r_info.builder.arm_joint_5_coord->q = M_PI * 0.25; 
    r_info.builder.arm_joint_6_coord->q = -M_PI * 0.125; 
    r_info.kin_chain->doMotion();
    
    r_info.airship_frame->Position = vect<double,3>(-0.8, -0.5, 1.4);
    r_info.airship_frame->Quat = axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0));
    r_info.airship_chain->doMotion();
    
    geom::oi_scene_graph sg;
#if 0
    geom::oi_scene_graph sg_tmp;
    sg_tmp << (*r_info.builder.get_geometric_model());
    
    sg.setCharacteristicLength(sg_tmp.computeCharacteristicLength());
    
    shared_ptr< kte::kte_map_chain > d_chain = r_info.builder.get_dynamics_kte_chain();
    
    shared_ptr< kte::spring_3D > spr1 = shared_ptr< kte::spring_3D >(new kte::spring_3D(
      "",
      r_info.builder.track_joint_end,
      r_info.builder.arm_EE,
      1.0,
      1.0));
    
    shared_ptr< kte::damper_3D > dmp1 = shared_ptr< kte::damper_3D >(new kte::damper_3D(
      "",
      r_info.builder.arm_joint_2_end,
      r_info.builder.arm_EE,
      1.0));
    
    shared_ptr< frame_3D<double> > f1 = shared_ptr< frame_3D<double> >(new frame_3D<double>());
    shared_ptr< frame_3D<double> > f2 = shared_ptr< frame_3D<double> >(new frame_3D<double>());
    shared_ptr< frame_3D<double> > f3 = shared_ptr< frame_3D<double> >(new frame_3D<double>());
    
    shared_ptr< kte::free_joint_3D > fj1 = shared_ptr< kte::free_joint_3D >(new kte::free_joint_3D(
      "",
      f1,
      f2,
      f3
    ));
    
    (*d_chain) << spr1 << dmp1 << fj1;
    
    sg << (*d_chain);
#else
    
    sg << (*r_info.builder.get_geometric_model());
    sg << (*lab_geom_model) << (*airship_geom_model);
    sg << geom::coord_arrows_3D("target_arrows",r_info.airship_frame,r_info.target_frame,0.3);
#endif
    
    r_info.sg_root = new SoSeparator;
    r_info.sg_root->ref();
    
    r_info.sg_root->addChild(sg.getSceneGraph());
    
    SoEventCallback* keypressCB = new SoEventCallback;
    keypressCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),
                                 keyboard_press_hdl, &r_info);
    r_info.sg_root->addChild(keypressCB);
    
    r_info.pp_root = new SoSeparator;
    r_info.sg_root->addChild(r_info.pp_root);
    
    sg.enableAnchorUpdates();
    
    // Use one of the convenient SoQt viewer classes.
    //SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
    SoQtExaminerViewer * eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(r_info.sg_root);
    eviewer->show();
    
    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();
    
    sg.disableAnchorUpdates();
    
    //delete cone_rot_anim;
    // Clean up resources.
    for(unsigned int i = 0; i < r_info.solution_seps.size(); ++i)
      r_info.solution_seps[i]->unref();
    delete eviewer;
    r_info.sg_root->unref();
  };
  
  
  SoQt::done();
  };
  
  return 0;
};







