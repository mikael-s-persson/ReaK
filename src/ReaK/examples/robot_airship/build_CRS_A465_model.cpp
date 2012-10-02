/**
 * \file build_CRS_A465_model.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */


#include "CRS_A465_geom_model.hpp"
#include "CRS_A465_2D_analog.hpp"


#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/nodes/SoSeparator.h>

#include "shapes/oi_scene_graph.hpp"

int main(int argc, char ** argv) {
  
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  
  builder.create_geom_from_preset();
  
  builder.save_kte_to_file("models/CRS_A465_raw_components.xml");
  builder.save_limits_to_file("models/CRS_A465_limits.xml");
  
  ReaK::robot_airship::CRS_A465_2D_model_builder builder2D;
  
  builder2D.create_from_preset();
  
  builder2D.save_kte_to_file("models/CRS_A465_2D_raw_components.xml");
  builder2D.save_limits_to_file("models/CRS_A465_2D_limits.xml");
  
  
  {
  QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
  
  
  {
    ReaK::geom::oi_scene_graph sg;
    
    builder.track_joint_coord->q = 0.2;
    builder.arm_joint_1_coord->q = M_PI * 0.25;
    builder.arm_joint_2_coord->q = -M_PI * 0.125;
    builder.arm_joint_3_coord->q = -M_PI * 0.375;
    builder.arm_joint_4_coord->q = M_PI * 0.125; 
    builder.arm_joint_5_coord->q = M_PI * 0.25; 
    builder.arm_joint_6_coord->q = -M_PI * 0.125; 
    builder.get_kinematics_kte_chain()->doMotion();
    
    sg << (*builder.get_geometric_model());
    
    // Use one of the convenient SoQt viewer classes.
    //SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
    SoQtExaminerViewer * eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(sg.getSceneGraph());
    eviewer->show();
    
    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();
    
    //delete cone_rot_anim;
    // Clean up resources.
    delete eviewer;
    
  };
  
  
  SoQt::done();
  };
  
  return 0;
};







