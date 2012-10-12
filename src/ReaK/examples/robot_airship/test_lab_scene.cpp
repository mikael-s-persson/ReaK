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
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/events/SoKeyboardEvent.h>

#include "shapes/oi_scene_graph.hpp"

#include "shapes/plane.hpp"
#include "shapes/box.hpp"
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


struct all_robot_info {
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > kin_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > robot_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > lab_proxy;
  ReaK::geom::proxy_query_pair_3D robot_lab_proxy;
  ReaK::shared_ptr< ReaK::frame_3D<double> > airship_frame;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > airship_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > airship_proxy;
  ReaK::geom::proxy_query_pair_3D robot_airship_proxy;
  ReaK::geom::proxy_query_pair_3D lab_airship_proxy;
  SoCoordinate3* l_r_proxy_line;
  SoCoordinate3* r_a_proxy_line;
  SoCoordinate3* l_a_proxy_line;
};


void keyboard_press_hdl(void* userData, SoEventCallback* eventCB) {
  all_robot_info* r_info = reinterpret_cast< all_robot_info* >(userData);
  const SoEvent* event = eventCB->getEvent();
  
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
  };
  
  
  r_info->kin_chain->doMotion();
  r_info->airship_chain->doMotion();
  
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > lr_pline = r_info->robot_lab_proxy.findMinimumDistance();
  if(lr_pline) {
    r_info->l_r_proxy_line->point.set1Value(0, lr_pline->getLastResult().mPoint1[0], lr_pline->getLastResult().mPoint1[1], lr_pline->getLastResult().mPoint1[2]);
    r_info->l_r_proxy_line->point.set1Value(1, lr_pline->getLastResult().mPoint2[0], lr_pline->getLastResult().mPoint2[1], lr_pline->getLastResult().mPoint2[2]);
  };
  
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > ra_pline = r_info->robot_airship_proxy.findMinimumDistance();
  if(ra_pline) {
    r_info->r_a_proxy_line->point.set1Value(0, ra_pline->getLastResult().mPoint1[0], ra_pline->getLastResult().mPoint1[1], ra_pline->getLastResult().mPoint1[2]);
    r_info->r_a_proxy_line->point.set1Value(1, ra_pline->getLastResult().mPoint2[0], ra_pline->getLastResult().mPoint2[1], ra_pline->getLastResult().mPoint2[2]);
  };
  
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > la_pline = r_info->lab_airship_proxy.findMinimumDistance();
  if(la_pline) {
    r_info->l_a_proxy_line->point.set1Value(0, la_pline->getLastResult().mPoint1[0], la_pline->getLastResult().mPoint1[1], la_pline->getLastResult().mPoint1[2]);
    r_info->l_a_proxy_line->point.set1Value(1, la_pline->getLastResult().mPoint2[0], la_pline->getLastResult().mPoint2[1], la_pline->getLastResult().mPoint2[2]);
  };
};


int main(int argc, char ** argv) {
  using namespace ReaK;
  
  all_robot_info r_info;
  
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
  
  r_info.robot_lab_proxy.setModelPair(r_info.robot_proxy, r_info.lab_proxy);
  r_info.robot_airship_proxy.setModelPair(r_info.robot_proxy, r_info.airship_proxy);
  r_info.lab_airship_proxy.setModelPair(r_info.lab_proxy, r_info.airship_proxy);
  
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
    
    geom::oi_scene_graph sg_tmp;
    sg_tmp << (*r_info.builder.get_geometric_model());
    
    geom::oi_scene_graph sg;
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
    
    //sg << (*r_info.kin_chain);
    //sg << (*lab_geom_model) << (*airship_geom_model);
    
    SoSeparator* root = new SoSeparator;
    root->ref();
    
    root->addChild(sg.getSceneGraph());
    
    
    
    shared_ptr< geom::proximity_finder_3D > lr_pline = r_info.robot_lab_proxy.findMinimumDistance();
    
    SoSeparator* sep_lr_pline = new SoSeparator;
    
    SoBaseColor * col_lr_pline = new SoBaseColor;
    col_lr_pline->rgb = SbColor(1, 0, 1);
    sep_lr_pline->addChild(col_lr_pline);
    
    SoCoordinate3* coords_lr_pline = new SoCoordinate3;
    if(lr_pline) {
      coords_lr_pline->point.set1Value(0, lr_pline->getLastResult().mPoint1[0], lr_pline->getLastResult().mPoint1[1], lr_pline->getLastResult().mPoint1[2]);
      coords_lr_pline->point.set1Value(1, lr_pline->getLastResult().mPoint2[0], lr_pline->getLastResult().mPoint2[1], lr_pline->getLastResult().mPoint2[2]);
    } else {
      coords_lr_pline->point.set1Value(0, 0.0, 0.0, 0.0);
      coords_lr_pline->point.set1Value(1, 0.0, 0.0, 0.0);
    };
    sep_lr_pline->addChild(coords_lr_pline);
    r_info.l_r_proxy_line = coords_lr_pline;
    
    SoLineSet* ln_set_lr_pline = new SoLineSet;
    ln_set_lr_pline->numVertices.set1Value(0, 2);
    sep_lr_pline->addChild(ln_set_lr_pline);
    
    //root->addChild(sep_lr_pline);
    
    
    
    shared_ptr< geom::proximity_finder_3D > ra_pline = r_info.robot_airship_proxy.findMinimumDistance();
    
    SoSeparator* sep_ra_pline = new SoSeparator;
    
    SoBaseColor * col_ra_pline = new SoBaseColor;
    col_ra_pline->rgb = SbColor(1, 0, 1);
    sep_ra_pline->addChild(col_ra_pline);
    
    SoCoordinate3* coords_ra_pline = new SoCoordinate3;
    if(ra_pline) {
      coords_ra_pline->point.set1Value(0, ra_pline->getLastResult().mPoint1[0], ra_pline->getLastResult().mPoint1[1], ra_pline->getLastResult().mPoint1[2]);
      coords_ra_pline->point.set1Value(1, ra_pline->getLastResult().mPoint2[0], ra_pline->getLastResult().mPoint2[1], ra_pline->getLastResult().mPoint2[2]);
    } else {
      coords_ra_pline->point.set1Value(0, 0.0, 0.0, 0.0);
      coords_ra_pline->point.set1Value(1, 0.0, 0.0, 0.0);
    };
    sep_ra_pline->addChild(coords_ra_pline);
    r_info.r_a_proxy_line = coords_ra_pline;
    
    SoLineSet* ln_set_ra_pline = new SoLineSet;
    ln_set_ra_pline->numVertices.set1Value(0, 2);
    sep_ra_pline->addChild(ln_set_ra_pline);
    
    //root->addChild(sep_ra_pline);
    
    
    
    shared_ptr< geom::proximity_finder_3D > la_pline = r_info.lab_airship_proxy.findMinimumDistance();
    
    SoSeparator* sep_la_pline = new SoSeparator;
    
    SoBaseColor * col_la_pline = new SoBaseColor;
    col_la_pline->rgb = SbColor(1, 0, 1);
    sep_la_pline->addChild(col_la_pline);
    
    SoCoordinate3* coords_la_pline = new SoCoordinate3;
    if(la_pline) {
      coords_la_pline->point.set1Value(0, la_pline->getLastResult().mPoint1[0], la_pline->getLastResult().mPoint1[1], la_pline->getLastResult().mPoint1[2]);
      coords_la_pline->point.set1Value(1, la_pline->getLastResult().mPoint2[0], la_pline->getLastResult().mPoint2[1], la_pline->getLastResult().mPoint2[2]);
    } else {
      coords_la_pline->point.set1Value(0, 0.0, 0.0, 0.0);
      coords_la_pline->point.set1Value(1, 0.0, 0.0, 0.0);
    };
    sep_la_pline->addChild(coords_la_pline);
    r_info.l_a_proxy_line = coords_la_pline;
    
    SoLineSet* ln_set_la_pline = new SoLineSet;
    ln_set_la_pline->numVertices.set1Value(0, 2);
    sep_la_pline->addChild(ln_set_la_pline);
    
    //root->addChild(sep_la_pline);
    
    
    
    SoEventCallback* keypressCB = new SoEventCallback;
    keypressCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),
                                 keyboard_press_hdl, &r_info);
    root->addChild(keypressCB);
    
    sg.enableAnchorUpdates();
    
    // Use one of the convenient SoQt viewer classes.
    //SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
    SoQtExaminerViewer * eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(root);
    eviewer->show();
    
    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();
    
    sg.disableAnchorUpdates();
    
    //delete cone_rot_anim;
    // Clean up resources.
    delete eviewer;
    root->unref();
  };
  
  
  SoQt::done();
  };
  
  return 0;
};







