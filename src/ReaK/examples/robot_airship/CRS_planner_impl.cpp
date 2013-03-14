
/*
 *    Copyright 2011 Sven Mikael Persson
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

#include "CRS_planner_impl.h"


#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>



#include "CRS_A465_geom_model.hpp"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

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
#include "mbd_kte/kte_map_chain.hpp"
#include "kte_models/manip_dynamics_model.hpp"

#include "serialization/xml_archiver.hpp"

#include "CRS_workspaces.hpp"
#include "CRS_rrt_planners.hpp"
#include "CRS_rrtstar_planners.hpp"
#include "CRS_prm_planners.hpp"
#include "CRS_fadprm_planners.hpp"
#include "CRS_sbastar_planners.hpp"

#include "optimization/optim_exceptions.hpp"

#include <chrono>

struct all_robot_info {
  SoQtExaminerViewer * eviewer;
  SoSeparator* sg_root;
  ReaK::geom::oi_scene_graph* sg_robot_geom;
  SoSwitch* sw_robot_geom;
  ReaK::geom::oi_scene_graph* sg_robot_kin;
  SoSwitch* sw_robot_kin;
  ReaK::geom::oi_scene_graph* sg_airship_geom;
  SoSwitch* sw_airship_geom;
  ReaK::geom::oi_scene_graph* sg_lab_geom;
  SoSwitch* sw_lab_geom;
  ReaK::geom::oi_scene_graph* sg_proxy_geom;
  SoSwitch* sw_proxy_geom;
  SoSwitch* sw_motion_graph;
  SoSwitch* sw_solutions;
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > kin_chain;
  ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > manip_kin_mdl;
  ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> > manip_jt_limits;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > robot_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > lab_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_lab_proxy;
  ReaK::shared_ptr< ReaK::frame_3D<double> > airship_frame;
  ReaK::pose_3D<double> target_frame;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > airship_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > lab_airship_proxy;
  SoCoordinate3* l_r_proxy_line;
  SoCoordinate3* r_a_proxy_line;
  SoCoordinate3* l_a_proxy_line;
  
  std::vector< ReaK::vect<double,7> > bestsol_trajectory;
  SoTimerSensor* animation_timer;
  std::size_t animation_progress;
  std::chrono::high_resolution_clock::time_point animation_last_render;
} r_info;


void CRSPlannerGUI_animate_bestsol_trajectory(void*, SoSensor*) {
  if(r_info.animation_progress < r_info.bestsol_trajectory.size()) {
    if(std::chrono::high_resolution_clock::now() - r_info.animation_last_render >= std::chrono::milliseconds(100)) {
      const ReaK::vect<double,7>& cur_pos = r_info.bestsol_trajectory[r_info.animation_progress];
      std::cout << "animate point: " << cur_pos << std::endl;
      r_info.builder.track_joint_coord->q = cur_pos[0];
      r_info.builder.arm_joint_1_coord->q = cur_pos[1];
      r_info.builder.arm_joint_2_coord->q = cur_pos[2];
      r_info.builder.arm_joint_3_coord->q = cur_pos[3];
      r_info.builder.arm_joint_4_coord->q = cur_pos[4]; 
      r_info.builder.arm_joint_5_coord->q = cur_pos[5]; 
      r_info.builder.arm_joint_6_coord->q = cur_pos[6]; 
      r_info.kin_chain->doMotion();
      r_info.animation_progress++;
      r_info.animation_last_render = std::chrono::high_resolution_clock::now();
    };
  } else {
    r_info.animation_timer->unschedule();
    r_info.animation_progress = 0;
    r_info.animation_last_render = std::chrono::high_resolution_clock::now();
  };
};



CRSPlannerGUI::CRSPlannerGUI( QWidget * parent, Qt::WindowFlags flags ) : QMainWindow(parent,flags),
                                                                          configs() {
  setupUi(this);
  using namespace ReaK;
  
  //connect(actionRun, SIGNAL(triggered()), this, SLOT(startPlanner()));
  //connect(actionStop_4, SIGNAL(triggered()), this, SLOT(stopPlanner()));
  //connect(actionLoad_World, SIGNAL(triggered()), this, SLOT(launchPlayerStage()));
  //connect(actionLoad_Planner, SIGNAL(triggered()), this, SLOT(loadPathPlanner()));
  //connect(actionLoad_Map, SIGNAL(triggered()), this, SLOT(loadWorldMap()));
  //connect(actionClose_Map, SIGNAL(triggered()), this, SLOT(closeTestScenario()));
  //connect(actionProperties, SIGNAL(triggered()), this, SLOT(setPlannerProperties()));
  //onnect(actionSave_Results, SIGNAL(triggered()), this, SLOT(saveResults()));
  
  //actionClose_Map->setEnabled(false);
  //actionStart_Robot->setEnabled(false);
  
  configs.setupUi(this->config_dock->widget());
  connect(configs.actionStart_Robot, SIGNAL(triggered()), this, SLOT(startRobot()));
  connect(configs.actionExecutePlanner, SIGNAL(triggered()), this, SLOT(executePlanner()));
  connect(configs.actionJointChange, SIGNAL(triggered()), this, SLOT(onJointChange()));
  connect(configs.actionTargetChange, SIGNAL(triggered()), this, SLOT(onTargetChange()));
  
  connect(configs.actionRobotVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotVisible()));
  connect(configs.actionRobotKinVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotKinVisible()));
  connect(configs.actionTargetVisibleToggle, SIGNAL(triggered()), this, SLOT(onTargetVisible()));
  connect(configs.actionEnvVisibleToggle, SIGNAL(triggered()), this, SLOT(onEnvVisible()));
  connect(configs.actionProxyVisibleToggle, SIGNAL(triggered()), this, SLOT(onProxyVisible()));
  connect(configs.actionMGVisibleToggle, SIGNAL(triggered()), this, SLOT(onMGVisible()));
  connect(configs.actionSolutionsVisibleToggle, SIGNAL(triggered()), this, SLOT(onSolutionsVisible()));
  
  
  r_info.animation_progress = 0;
  r_info.animation_timer = new SoTimerSensor(CRSPlannerGUI_animate_bestsol_trajectory, NULL);
  
  //r_info.builder.create_geom_from_preset();
  r_info.builder.load_kte_and_geom("models/CRS_A465_with_geom.xml");
  r_info.builder.load_limits_from_file("models/CRS_A465_limits.xml");
  
  r_info.robot_proxy = r_info.builder.get_proximity_model();
  r_info.kin_chain = r_info.builder.get_kinematics_kte_chain();
  r_info.manip_kin_mdl = r_info.builder.get_manipulator_kin_model();
  r_info.manip_jt_limits = ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> >(&(r_info.builder.joint_rate_limits), ReaK::null_deleter());
  
  
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
  
  
  
  SoQt::init(this->centralwidget);
  
  r_info.sg_root = new SoSeparator;
  r_info.sg_root->ref();
  
  
  r_info.sw_robot_geom = new SoSwitch();
  r_info.sg_robot_geom = new geom::oi_scene_graph();
  
  (*r_info.sg_robot_geom) << (*r_info.builder.get_geometric_model());
  double charact_length = r_info.sg_robot_geom->computeCharacteristicLength();
  
  r_info.sw_robot_geom->addChild(r_info.sg_robot_geom->getSceneGraph());
  r_info.sw_robot_geom->whichChild.setValue((configs.check_show_geom->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  r_info.sg_root->addChild(r_info.sw_robot_geom);
  
  
  r_info.sw_robot_kin = new SoSwitch();
  r_info.sg_robot_kin = new geom::oi_scene_graph();
  
  r_info.sg_robot_kin->setCharacteristicLength(charact_length);
  (*r_info.sg_robot_kin) << (*r_info.kin_chain);
  
  r_info.sw_robot_kin->addChild(r_info.sg_robot_kin->getSceneGraph());
  r_info.sw_robot_kin->whichChild.setValue((configs.check_show_kinmdl->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  r_info.sg_root->addChild(r_info.sw_robot_kin);
  
  
  r_info.sw_airship_geom = new SoSwitch();
  r_info.sg_airship_geom = new geom::oi_scene_graph();
  
  (*r_info.sg_airship_geom) << (*airship_geom_model) 
                         << geom::coord_arrows_3D("target_arrows",r_info.airship_frame,r_info.target_frame,0.3);
  
  r_info.sw_airship_geom->addChild(r_info.sg_airship_geom->getSceneGraph());
  r_info.sw_airship_geom->whichChild.setValue((configs.check_show_target->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  r_info.sg_root->addChild(r_info.sw_airship_geom);
  
  
  r_info.sw_lab_geom = new SoSwitch();
  r_info.sg_lab_geom = new geom::oi_scene_graph();
  
  (*r_info.sg_lab_geom) << (*lab_geom_model);
  
  r_info.sw_lab_geom->addChild(r_info.sg_lab_geom->getSceneGraph());
  r_info.sw_lab_geom->whichChild.setValue((configs.check_show_env->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  r_info.sg_root->addChild(r_info.sw_lab_geom);
  
  
  r_info.sw_motion_graph = new SoSwitch();
  r_info.sw_motion_graph->whichChild.setValue((configs.check_show_motiongraph->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  r_info.sg_root->addChild(r_info.sw_motion_graph);
  
  
  r_info.sw_solutions = new SoSwitch();
  r_info.sw_solutions->whichChild.setValue((configs.check_show_sol->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  r_info.sg_root->addChild(r_info.sw_solutions);
  
  
  
  
  
  r_info.sw_proxy_geom = new SoSwitch();
  
  shared_ptr< geom::proximity_finder_3D > lr_pline = r_info.robot_lab_proxy->findMinimumDistance();
  
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
  
  r_info.sw_proxy_geom->addChild(sep_lr_pline);
  
  
  shared_ptr< geom::proximity_finder_3D > ra_pline = r_info.robot_airship_proxy->findMinimumDistance();
  
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
  
  r_info.sw_proxy_geom->addChild(sep_ra_pline);
  
  
  shared_ptr< geom::proximity_finder_3D > la_pline = r_info.lab_airship_proxy->findMinimumDistance();
  
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
  
  r_info.sw_proxy_geom->addChild(sep_la_pline);
  
  
  r_info.sg_root->addChild(r_info.sw_proxy_geom);
  
  
  
  onJointChange();
  onTargetChange();
  
  
  r_info.eviewer = new SoQtExaminerViewer(this->centralwidget);
  r_info.eviewer->setSceneGraph(r_info.sg_root);
  r_info.eviewer->show();
  
  r_info.sg_robot_geom->enableAnchorUpdates();
  r_info.sg_robot_kin->enableAnchorUpdates();
  r_info.sg_airship_geom->enableAnchorUpdates();
};


CRSPlannerGUI::~CRSPlannerGUI() {
  
  r_info.sg_robot_geom->disableAnchorUpdates();
  r_info.sg_robot_kin->disableAnchorUpdates();
  r_info.sg_airship_geom->disableAnchorUpdates();
  
  delete r_info.sg_robot_geom;
  delete r_info.sg_robot_kin;
  delete r_info.sg_airship_geom;
  delete r_info.sg_lab_geom;
  
  delete r_info.animation_timer;
  
  delete r_info.eviewer;
  r_info.sg_root->unref();
  SoQt::done();
};


void CRSPlannerGUI::onJointChange() {
  r_info.builder.track_joint_coord->q = double(configs.track_pos->value()) * 0.001;
  r_info.builder.arm_joint_1_coord->q = double(configs.joint1_pos->value()) * 0.001;
  r_info.builder.arm_joint_2_coord->q = double(configs.joint2_pos->value()) * 0.001;
  r_info.builder.arm_joint_3_coord->q = double(configs.joint3_pos->value()) * 0.001;
  r_info.builder.arm_joint_4_coord->q = double(configs.joint4_pos->value()) * 0.001; 
  r_info.builder.arm_joint_5_coord->q = double(configs.joint5_pos->value()) * 0.001; 
  r_info.builder.arm_joint_6_coord->q = double(configs.joint6_pos->value()) * 0.001; 
  r_info.kin_chain->doMotion();
  onProxyChange();
};

void CRSPlannerGUI::onTargetChange() {
  r_info.airship_frame->Position = ReaK::vect<double,3>(
    double(configs.target_x->value()) * 0.001, 
    double(configs.target_y->value()) * 0.001, 
    double(configs.target_z->value()) * 0.001);
  r_info.airship_frame->Quat = 
    ReaK::quaternion<double>::zrot(double(configs.target_yaw->value()) * 0.001) * 
    ReaK::quaternion<double>::yrot(double(configs.target_pitch->value()) * 0.001) * 
    ReaK::quaternion<double>::xrot(double(configs.target_roll->value()) * 0.001);
  r_info.airship_chain->doMotion();
  
  
  if(configs.check_enable_ik->isChecked()) {
    try {
      ReaK::vect_n<double> jt_sol = r_info.builder.compute_inverse_kinematics(r_info.target_frame.getGlobalPose());
      r_info.builder.track_joint_coord->q = jt_sol[0];
      r_info.builder.arm_joint_1_coord->q = jt_sol[1];
      r_info.builder.arm_joint_2_coord->q = jt_sol[2];
      r_info.builder.arm_joint_3_coord->q = jt_sol[3];
      r_info.builder.arm_joint_4_coord->q = jt_sol[4];
      r_info.builder.arm_joint_5_coord->q = jt_sol[5];
      r_info.builder.arm_joint_6_coord->q = jt_sol[6];
    } catch( ReaK::optim::infeasible_problem& e ) { RK_UNUSED(e); };
    r_info.kin_chain->doMotion();
  };
  
  onProxyChange();
};

void CRSPlannerGUI::onRobotVisible() {
  r_info.sw_robot_geom->whichChild.setValue((configs.check_show_geom->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onRobotKinVisible() {
  r_info.sw_robot_kin->whichChild.setValue((configs.check_show_kinmdl->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onTargetVisible() {
  r_info.sw_airship_geom->whichChild.setValue((configs.check_show_target->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onEnvVisible() {
  r_info.sw_lab_geom->whichChild.setValue((configs.check_show_env->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onProxyVisible() {
  r_info.sw_proxy_geom->whichChild.setValue((configs.check_show_proxy->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onMGVisible() {
  r_info.sw_motion_graph->whichChild.setValue((configs.check_show_motiongraph->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onSolutionsVisible() {
  r_info.sw_solutions->whichChild.setValue((configs.check_show_sol->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onProxyChange() {
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > lr_pline = r_info.robot_lab_proxy->findMinimumDistance();
  if(lr_pline) {
    r_info.l_r_proxy_line->point.set1Value(0, lr_pline->getLastResult().mPoint1[0], lr_pline->getLastResult().mPoint1[1], lr_pline->getLastResult().mPoint1[2]);
    r_info.l_r_proxy_line->point.set1Value(1, lr_pline->getLastResult().mPoint2[0], lr_pline->getLastResult().mPoint2[1], lr_pline->getLastResult().mPoint2[2]);
  };
  
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > ra_pline = r_info.robot_airship_proxy->findMinimumDistance();
  if(ra_pline) {
    r_info.r_a_proxy_line->point.set1Value(0, ra_pline->getLastResult().mPoint1[0], ra_pline->getLastResult().mPoint1[1], ra_pline->getLastResult().mPoint1[2]);
    r_info.r_a_proxy_line->point.set1Value(1, ra_pline->getLastResult().mPoint2[0], ra_pline->getLastResult().mPoint2[1], ra_pline->getLastResult().mPoint2[2]);
  };
  
  ReaK::shared_ptr< ReaK::geom::proximity_finder_3D > la_pline = r_info.lab_airship_proxy->findMinimumDistance();
  if(la_pline) {
    r_info.l_a_proxy_line->point.set1Value(0, la_pline->getLastResult().mPoint1[0], la_pline->getLastResult().mPoint1[1], la_pline->getLastResult().mPoint1[2]);
    r_info.l_a_proxy_line->point.set1Value(1, la_pline->getLastResult().mPoint2[0], la_pline->getLastResult().mPoint2[1], la_pline->getLastResult().mPoint2[2]);
  };
};



void CRSPlannerGUI::executePlanner() {
  
  ReaK::vect_n<double> jt_desired(7,0.0);
  if(configs.check_ik_goal->isChecked()) {
    try {
      
      jt_desired = r_info.builder.compute_inverse_kinematics(r_info.target_frame.getGlobalPose());
      
    } catch( ReaK::optim::infeasible_problem& e ) { RK_UNUSED(e);
      QMessageBox::critical(this,
                    "Inverse Kinematics Error!",
                    "The target frame cannot be reached! No inverse kinematics solution possible!",
                    QMessageBox::Ok);
      return;
    };
  } else {
    std::stringstream ss(configs.custom_goal_edit->text().toStdString());
    ss >> jt_desired;
  };
  
  ReaK::vect_n<double> jt_start;
  if(configs.check_current_start->isChecked()) {
    jt_start.resize(7);
    jt_start[0] = r_info.builder.track_joint_coord->q;
    jt_start[1] = r_info.builder.arm_joint_1_coord->q;
    jt_start[2] = r_info.builder.arm_joint_2_coord->q;
    jt_start[3] = r_info.builder.arm_joint_3_coord->q;
    jt_start[4] = r_info.builder.arm_joint_4_coord->q; 
    jt_start[5] = r_info.builder.arm_joint_5_coord->q; 
    jt_start[6] = r_info.builder.arm_joint_6_coord->q; 
  } else {
    std::stringstream ss(configs.custom_start_edit->text().toStdString());
    ss >> jt_start;
  };
  
  
  // joint-space parameters:
  std::size_t space_order = configs.order_selection->currentIndex();
  std::size_t interp_id = configs.interp_selection->currentIndex();
  
  double min_travel = configs.min_interval_spinbox->value();
  double max_travel = configs.max_interval_spinbox->value();
  
  
  // planner parameters:
  std::size_t max_vertices = configs.maxvertices_spinbox->value();
  std::size_t prog_interval = configs.progress_interval_spinbox->value();
  std::size_t max_results = configs.maxsolutions_spinbox->value();
  
  std::size_t rrt_dir = ReaK::pp::UNIDIRECTIONAL_RRT;
  if(configs.direction_selection->currentIndex())
    rrt_dir = ReaK::pp::BIDIRECTIONAL_RRT;
  
  std::size_t store_policy = ReaK::pp::ADJ_LIST_MOTION_GRAPH;
  if(configs.graph_storage_selection->currentIndex())
    store_policy = ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH;
  
  std::size_t knn_method = ReaK::pp::LINEAR_SEARCH_KNN;
  switch(configs.KNN_method_selection->currentIndex()) {
    case 1:
      if(store_policy == ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH)
        knn_method = ReaK::pp::DVP_ALT_BF2_KNN;
      else
        knn_method = ReaK::pp::DVP_BF2_TREE_KNN;
      break;
    case 2:
      if(store_policy == ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH)
        knn_method = ReaK::pp::DVP_ALT_BF4_KNN;
      else
        knn_method = ReaK::pp::DVP_BF4_TREE_KNN;
      break;
    case 3:
      if(store_policy == ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH)
        knn_method = ReaK::pp::DVP_ALT_COB2_KNN;
      else
        knn_method = ReaK::pp::DVP_COB2_TREE_KNN;
      break;
    case 4:
      if(store_policy == ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH)
        knn_method = ReaK::pp::DVP_ALT_COB4_KNN;
      else
        knn_method = ReaK::pp::DVP_COB4_TREE_KNN;
      break;
    default:
      break;
  };
  
  SoSeparator* mg_sep = NULL;
  std::vector< SoSeparator* > sol_seps;
  
  
  /*
   * Below are a few rather large MACROs that are used to generate all the code for the different planner-space
   * combinations.
   */
  
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(WORKSPACE) \
        typedef ReaK::pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_rl_o0_type > jt_space(new ReaK::robot_airship::CRS3D_jspace_rl_o0_type(r_info.builder.get_rl_joint_space_0th())); \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_o0_type > normal_jt_space(new ReaK::robot_airship::CRS3D_jspace_o0_type(r_info.builder.get_joint_space_0th())); \
         \
        ReaK::shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            r_info.manip_kin_mdl, \
            r_info.manip_jt_limits, \
            min_travel, max_travel)); \
         \
        (*workspace) << r_info.robot_lab_proxy << r_info.robot_airship_proxy; \
         \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_rl_o0_type >::point_type RLPointType; \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_o0_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = ReaK::vect<double,7>(jt_start); \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = ReaK::vect<double,7>(jt_desired); \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::robot_airship::CRS3D_rl_o0_tracer frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          ReaK::robot_airship::CRS3D_rlDK_o0_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * min_travel); \
        temp_reporter.add_traced_frame(r_info.builder.arm_joint_6_end);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(WORKSPACE) \
        typedef ReaK::pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_rl_o1_type > jt_space(new ReaK::robot_airship::CRS3D_jspace_rl_o1_type(r_info.builder.get_rl_joint_space_1st())); \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_o1_type > normal_jt_space(new ReaK::robot_airship::CRS3D_jspace_o1_type(r_info.builder.get_joint_space_1st())); \
         \
        ReaK::shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            r_info.manip_kin_mdl, \
            r_info.manip_jt_limits, \
            min_travel, max_travel)); \
         \
        (*workspace) << r_info.robot_lab_proxy << r_info.robot_airship_proxy; \
         \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_rl_o1_type >::point_type RLPointType; \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_o1_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = ReaK::vect<double,7>(jt_start); \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = ReaK::vect<double,7>(jt_desired); \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::robot_airship::CRS3D_rl_o1_tracer frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          ReaK::robot_airship::CRS3D_rlDK_o1_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * min_travel); \
        temp_reporter.add_traced_frame(r_info.builder.arm_joint_6_end);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(WORKSPACE) \
        typedef ReaK::pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_rl_o2_type > jt_space(new ReaK::robot_airship::CRS3D_jspace_rl_o2_type(r_info.builder.get_rl_joint_space())); \
        ReaK::shared_ptr< ReaK::robot_airship::CRS3D_jspace_o2_type > normal_jt_space(new ReaK::robot_airship::CRS3D_jspace_o2_type(r_info.builder.get_joint_space())); \
         \
        ReaK::shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            r_info.manip_kin_mdl, \
            r_info.manip_jt_limits, \
            min_travel, max_travel)); \
         \
        (*workspace) << r_info.robot_lab_proxy << r_info.robot_airship_proxy; \
         \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_rl_o2_type >::point_type RLPointType; \
        typedef ReaK::pp::topology_traits< ReaK::robot_airship::CRS3D_jspace_o2_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = ReaK::vect<double,7>(jt_start); \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = ReaK::vect<double,7>(jt_desired); \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::robot_airship::CRS3D_rl_o2_tracer frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          ReaK::robot_airship::CRS3D_rlDK_o2_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * min_travel); \
        temp_reporter.add_traced_frame(r_info.builder.arm_joint_6_end);
        
#define RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(WORKSPACE) \
        ReaK::pp::rrt_path_planner<WORKSPACE, frame_reporter_type> \
          workspace_planner( \
            workspace, \
            start_point, \
            goal_point, \
            max_vertices, \
            prog_interval, \
            rrt_dir, \
            store_policy, \
            knn_method, \
            temp_reporter,  \
            max_results); \
         \
        ReaK::shared_ptr< ReaK::pp::seq_path_base< SuperSpaceType > > bestsol_rlpath = workspace_planner.solve_path(); \
         \
        r_info.bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef ReaK::pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            r_info.bestsol_trajectory.push_back( get<0>(r_info.manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        r_info.builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        r_info.builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        r_info.builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        r_info.builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        r_info.builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        r_info.builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        r_info.kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(WORKSPACE) \
        ReaK::pp::rrtstar_path_planner<WORKSPACE, frame_reporter_type> \
          workspace_planner( \
            workspace, \
            start_point, \
            goal_point, \
            max_vertices, \
            prog_interval, \
            rrt_dir, \
            store_policy, \
            knn_method, \
            temp_reporter, \
            max_results); \
         \
        ReaK::shared_ptr< ReaK::pp::seq_path_base< SuperSpaceType > > bestsol_rlpath = workspace_planner.solve_path(); \
         \
        r_info.bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef ReaK::pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            r_info.bestsol_trajectory.push_back( get<0>(r_info.manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        r_info.builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        r_info.builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        r_info.builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        r_info.builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        r_info.builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        r_info.builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        r_info.kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(WORKSPACE) \
        ReaK::pp::prm_path_planner<WORKSPACE, frame_reporter_type> \
          workspace_planner( \
            workspace, \
            start_point, \
            goal_point, \
            max_vertices, \
            prog_interval, \
            store_policy, \
            knn_method, \
            temp_reporter, \
            max_results); \
         \
        ReaK::shared_ptr< ReaK::pp::seq_path_base< SuperSpaceType > > bestsol_rlpath = workspace_planner.solve_path(); \
         \
        r_info.bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef ReaK::pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            r_info.bestsol_trajectory.push_back( get<0>(r_info.manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        r_info.builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        r_info.builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        r_info.builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        r_info.builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        r_info.builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        r_info.builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        r_info.kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(WORKSPACE) \
        ReaK::pp::fadprm_path_planner<WORKSPACE, frame_reporter_type> \
          workspace_planner( \
            workspace, \
            start_point, \
            goal_point, \
            0.1, \
            max_vertices, \
            prog_interval, \
            store_policy, \
            knn_method, \
            temp_reporter, \
            max_results); \
         \
        ReaK::shared_ptr< ReaK::pp::seq_path_base< SuperSpaceType > > bestsol_rlpath = workspace_planner.solve_path(); \
         \
        r_info.bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef ReaK::pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            r_info.bestsol_trajectory.push_back( get<0>(r_info.manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        r_info.builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        r_info.builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        r_info.builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        r_info.builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        r_info.builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        r_info.builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        r_info.kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(WORKSPACE) \
        ReaK::pp::sbastar_path_planner<WORKSPACE, frame_reporter_type> \
          workspace_planner( \
            workspace, \
            start_point, \
            goal_point, \
            0.9, \
            0.9, \
            max_travel, \
            max_vertices, \
            prog_interval, \
            store_policy, \
            knn_method, \
            ReaK::pp::LAZY_COLLISION_CHECKING, \
            ReaK::pp::PLAN_WITH_VORONOI_PULL, \
            temp_reporter, \
            max_results); \
         \
        ReaK::shared_ptr< ReaK::pp::seq_path_base< SuperSpaceType > > bestsol_rlpath = workspace_planner.solve_path(); \
        std::cout << "The shortest distance is: " << workspace_planner.get_best_solution_distance() << std::endl; \
        r_info.bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef ReaK::pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            r_info.bestsol_trajectory.push_back( get<0>(r_info.manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        r_info.builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        r_info.builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        r_info.builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        r_info.builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        r_info.builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        r_info.builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        r_info.kin_chain->doMotion();
  
  
  
  
  switch(configs.planning_algo_selection->currentIndex()) {
    case 0:  // RRT
    {
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
    case 1:  // RRT*
    {
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
    case 2:  // PRM
    {
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
    case 3:  // SBA*
    { 
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
    case 4:  // FADPRM
    { 
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(ReaK::robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
    case 5:  // LSBA*
    { 
      
    };// break;
    case 6:  // ???
    { 
      
    };// break;
    default:
      QMessageBox::information(this,
                  "Planner Not Supported!",
                  "Sorry, the planning method you selected is not yet supported!",
                  QMessageBox::Ok);
  };
  
  // Check the motion-graph separator and solution separators
  //  add them to the switches.
  if(mg_sep) {
    r_info.sw_motion_graph->removeAllChildren();
    r_info.sw_motion_graph->addChild(mg_sep);
    mg_sep->unref();
  };
  
  r_info.sw_solutions->removeAllChildren();
  if(configs.check_print_allsol->isChecked()) {
    for(std::size_t i = 0; i < sol_seps.size(); ++i) {
      r_info.sw_solutions->addChild(sol_seps[i]);
      sol_seps[i]->unref();
    };
  } else if((configs.check_print_best->isChecked()) && (sol_seps.size())) {
    r_info.sw_solutions->addChild(sol_seps[0]);
    sol_seps[0]->unref();
  };
  
};

void CRSPlannerGUI::startRobot() {
  r_info.animation_progress = 0;
  r_info.animation_last_render = std::chrono::high_resolution_clock::now();
  r_info.animation_timer->schedule();
};


int main(int argc, char** argv) {
  QApplication app(argc,argv);
  CRSPlannerGUI window;
  window.show();
  // Pop up the main window.
  SoQt::show(&window);
  // Loop until exit.
  SoQt::mainLoop();
  
  return 0;
  //return app.exec();
};








