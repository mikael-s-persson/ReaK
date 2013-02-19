/**
 * \file X8_test_scene.cpp
 *
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */


#include "X8_quadrotor_model.hpp"
#include "X8_quadrotor_geom.hpp"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/events/SoKeyboardEvent.h>

#include "shapes/oi_scene_graph.hpp"
#include "proximity/proxy_query_model.hpp"

#include "shapes/coord_arrows_3D.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"
#include "optimization/optim_exceptions.hpp"

using namespace ReaK;

struct all_robot_info {
  kte::X8_quadrotor_kinematics builder;
  geom::X8_quadrotor_geom geom_builder;
  shared_ptr< kte::kte_map_chain > kin_chain;
  shared_ptr< frame_3D<double> > target_frame;
  SoSwitch* sw_proxy_show;
};


void keyboard_press_hdl(void* userData, SoEventCallback* eventCB) {
  all_robot_info* r_info = reinterpret_cast< all_robot_info* >(userData);
  const SoEvent* event = eventCB->getEvent();
  
  static bool proxy_show = false;
  
  vect_n<double> j_pos = r_info->builder.getJointPositions();
  
  if( SO_KEY_PRESS_EVENT(event, Z) ) {
    r_info->target_frame->Position[2] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, X) ) {
    r_info->target_frame->Position[2] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, UP_ARROW) ) {
    r_info->target_frame->Position[1] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, DOWN_ARROW) ) {
    r_info->target_frame->Position[1] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, LEFT_ARROW) ) {
    r_info->target_frame->Position[0] -= 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, RIGHT_ARROW) ) {
    r_info->target_frame->Position[0] += 0.01;
  } else if( SO_KEY_PRESS_EVENT(event, B) ) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(1.0,0.0,0.0));
  } else if( SO_KEY_PRESS_EVENT(event, N) ) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(0.0,1.0,0.0));
  } else if( SO_KEY_PRESS_EVENT(event, M) ) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(M_PI * 0.01,ReaK::vect<double,3>(0.0,0.0,1.0));
  } else if( SO_KEY_PRESS_EVENT(event, L) ) {
    proxy_show = !proxy_show;
    r_info->sw_proxy_show->whichChild.setValue((proxy_show ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  };
  
  r_info->kin_chain->doMotion();
};


int main(int argc, char ** argv) {
  using namespace ReaK;
  
  all_robot_info r_info;
  
  
  r_info.kin_chain = r_info.builder.getKTEChain();
  r_info.geom_builder.create_geom_from_preset(r_info.builder);
  
  r_info.target_frame = r_info.builder.getFrame3D(0);
  
  r_info.target_frame->Position = vect<double,3>(0.0, 0.0, 1.0);
  r_info.target_frame->Quat = quaternion<double>::xrot(M_PI).getQuaternion();
  
  
  shared_ptr< geom::colored_model_3D > world_geom_model;
  shared_ptr< geom::proxy_query_model_3D > world_geom_proxy;
  {
    serialization::xml_iarchive in("models/window_crossing.xml");
    in >> world_geom_model >> world_geom_proxy;
  };
  
  
  {
  QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
  
  
  {
    r_info.kin_chain->doMotion();
    
    geom::oi_scene_graph sg;

    sg.setCharacteristicLength(10.0);
    
    sg << (*r_info.geom_builder.get_geometric_model());
    sg << *world_geom_model;
    sg << geom::coord_arrows_3D("base_arrows", shared_ptr< frame_3D<double> >(new frame_3D<double>()), pose_3D<double>(), 0.3);
    
    SoSeparator* root = new SoSeparator;
    root->ref();
    
    root->addChild(sg.getSceneGraph());
    
    SoEventCallback* keypressCB = new SoEventCallback;
    keypressCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),
                                 keyboard_press_hdl, &r_info);
    root->addChild(keypressCB);
    
    
    r_info.sw_proxy_show = new SoSwitch();
    geom::oi_scene_graph sg_proxy;
    
    sg_proxy << *(r_info.geom_builder.get_proximity_model()->mShapeList[0]);
    r_info.sw_proxy_show->addChild(sg_proxy.getSceneGraph());
    r_info.sw_proxy_show->whichChild.setValue(SO_SWITCH_NONE);
    
    root->addChild(r_info.sw_proxy_show);
    
    sg.enableAnchorUpdates();
    sg_proxy.enableAnchorUpdates();
    
    // Use one of the convenient SoQt viewer classes.
    //SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
    SoQtExaminerViewer * eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(root);
    eviewer->show();
    
    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();
    
    sg_proxy.disableAnchorUpdates();
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







