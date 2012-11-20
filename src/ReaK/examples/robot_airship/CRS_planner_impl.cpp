
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
#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"

#include "interpolation/linear_interp.hpp"
#include "interpolation/cubic_hermite_interp.hpp"
#include "interpolation/quintic_hermite_interp.hpp"
#include "interpolation/sustained_velocity_pulse.hpp"
#include "interpolation/sustained_acceleration_pulse.hpp"
#include "topologies/manip_free_workspace.hpp"
#include "path_planning/rrt_path_planner.hpp"
#include "path_planning/rrtstar_path_planner.hpp"
#include "path_planning/frame_tracer_coin3d.hpp"
#include "optimization/optim_exceptions.hpp"



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
} r_info;


typedef ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_0th_type jspace_rl_0_type;
typedef ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_1st_type jspace_rl_1_type;
typedef ReaK::robot_airship::CRS_A465_model_builder::rate_limited_joint_space_type jspace_rl_2_type;

typedef ReaK::robot_airship::CRS_A465_model_builder::joint_space_0th_type jspace_0_type;
typedef ReaK::robot_airship::CRS_A465_model_builder::joint_space_1st_type jspace_1_type;
typedef ReaK::robot_airship::CRS_A465_model_builder::joint_space_type jspace_2_type;

typedef ReaK::pp::manip_quasi_static_env<jspace_rl_0_type, ReaK::pp::linear_interpolation_tag> workspace_0_linear_type;
typedef ReaK::pp::manip_quasi_static_env<jspace_rl_1_type, ReaK::pp::linear_interpolation_tag> workspace_1_linear_type;
typedef ReaK::pp::manip_quasi_static_env<jspace_rl_2_type, ReaK::pp::linear_interpolation_tag> workspace_2_linear_type;

typedef ReaK::pp::manip_quasi_static_env<jspace_rl_1_type, ReaK::pp::cubic_hermite_interpolation_tag> workspace_1_cubic_type;
typedef ReaK::pp::manip_quasi_static_env<jspace_rl_2_type, ReaK::pp::cubic_hermite_interpolation_tag> workspace_2_cubic_type;

typedef ReaK::pp::manip_quasi_static_env<jspace_rl_2_type, ReaK::pp::quintic_hermite_interpolation_tag> workspace_2_quintic_type;

typedef ReaK::pp::manip_quasi_static_env<jspace_rl_1_type, ReaK::pp::svp_interpolation_tag> workspace_1_svp_type;
typedef ReaK::pp::manip_quasi_static_env<jspace_rl_2_type, ReaK::pp::svp_interpolation_tag> workspace_2_svp_type;

typedef ReaK::pp::manip_quasi_static_env<jspace_rl_2_type, ReaK::pp::sap_interpolation_tag> workspace_2_sap_type;


typedef ReaK::pp::manip_rl_direct_kin_map< ReaK::pp::joint_limits_collection<double>, jspace_0_type > rldk_0_type;
typedef ReaK::pp::manip_rl_direct_kin_map< ReaK::pp::joint_limits_collection<double>, jspace_1_type > rldk_1_type;
typedef ReaK::pp::manip_rl_direct_kin_map< ReaK::pp::joint_limits_collection<double>, jspace_2_type > rldk_2_type;


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
        ReaK::shared_ptr< jspace_rl_0_type > jt_space(new jspace_rl_0_type(r_info.builder.get_rl_joint_space_0th())); \
        ReaK::shared_ptr< jspace_0_type > normal_jt_space(new jspace_0_type(r_info.builder.get_joint_space_0th())); \
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
        ReaK::pp::topology_traits< jspace_rl_0_type >::point_type start_point, goal_point; \
        ReaK::pp::topology_traits< jspace_0_type >::point_type start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(get<0>(start_inter)) = jt_start[0]; \
        get<0>(get<1>(start_inter)) = jt_start[1]; \
        get<0>(get<2>(start_inter)) = jt_start[2]; \
        get<0>(get<3>(start_inter)) = jt_start[3]; \
        get<0>(get<4>(start_inter)) = jt_start[4]; \
        get<0>(get<5>(start_inter)) = jt_start[5]; \
        get<0>(get<6>(start_inter)) = jt_start[6]; \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(get<0>(goal_inter)) = jt_desired[0]; \
        get<0>(get<1>(goal_inter)) = jt_desired[1]; \
        get<0>(get<2>(goal_inter)) = jt_desired[2]; \
        get<0>(get<3>(goal_inter)) = jt_desired[3]; \
        get<0>(get<4>(goal_inter)) = jt_desired[4]; \
        get<0>(get<5>(goal_inter)) = jt_desired[5]; \
        get<0>(get<6>(goal_inter)) = jt_desired[6]; \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::pp::frame_tracer_3D< \
          rldk_0_type, jspace_rl_0_type, ReaK::pp::identity_topo_map, ReaK::pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          rldk_0_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * min_travel); \
        temp_reporter.add_traced_frame(r_info.builder.arm_joint_6_end);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(WORKSPACE) \
        ReaK::shared_ptr< jspace_rl_1_type > jt_space(new jspace_rl_1_type(r_info.builder.get_rl_joint_space_1st())); \
        ReaK::shared_ptr< jspace_1_type > normal_jt_space(new jspace_1_type(r_info.builder.get_joint_space_1st())); \
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
        ReaK::pp::topology_traits< jspace_rl_1_type >::point_type start_point, goal_point; \
        ReaK::pp::topology_traits< jspace_1_type >::point_type start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(get<0>(start_inter)) = jt_start[0]; \
        get<0>(get<1>(start_inter)) = jt_start[1]; \
        get<0>(get<2>(start_inter)) = jt_start[2]; \
        get<0>(get<3>(start_inter)) = jt_start[3]; \
        get<0>(get<4>(start_inter)) = jt_start[4]; \
        get<0>(get<5>(start_inter)) = jt_start[5]; \
        get<0>(get<6>(start_inter)) = jt_start[6]; \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(get<0>(goal_inter)) = jt_desired[0]; \
        get<0>(get<1>(goal_inter)) = jt_desired[1]; \
        get<0>(get<2>(goal_inter)) = jt_desired[2]; \
        get<0>(get<3>(goal_inter)) = jt_desired[3]; \
        get<0>(get<4>(goal_inter)) = jt_desired[4]; \
        get<0>(get<5>(goal_inter)) = jt_desired[5]; \
        get<0>(get<6>(goal_inter)) = jt_desired[6]; \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::pp::frame_tracer_3D< \
          rldk_1_type, jspace_rl_1_type, ReaK::pp::identity_topo_map, ReaK::pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          rldk_1_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * min_travel); \
        temp_reporter.add_traced_frame(r_info.builder.arm_joint_6_end);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(WORKSPACE) \
        ReaK::shared_ptr< jspace_rl_2_type > jt_space(new jspace_rl_2_type(r_info.builder.get_rl_joint_space())); \
        ReaK::shared_ptr< jspace_2_type > normal_jt_space(new jspace_2_type(r_info.builder.get_joint_space())); \
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
        ReaK::pp::topology_traits< jspace_rl_2_type >::point_type start_point, goal_point; \
        ReaK::pp::topology_traits< jspace_2_type >::point_type start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(get<0>(start_inter)) = jt_start[0]; \
        get<0>(get<1>(start_inter)) = jt_start[1]; \
        get<0>(get<2>(start_inter)) = jt_start[2]; \
        get<0>(get<3>(start_inter)) = jt_start[3]; \
        get<0>(get<4>(start_inter)) = jt_start[4]; \
        get<0>(get<5>(start_inter)) = jt_start[5]; \
        get<0>(get<6>(start_inter)) = jt_start[6]; \
        start_point = r_info.manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(get<0>(goal_inter)) = jt_desired[0]; \
        get<0>(get<1>(goal_inter)) = jt_desired[1]; \
        get<0>(get<2>(goal_inter)) = jt_desired[2]; \
        get<0>(get<3>(goal_inter)) = jt_desired[3]; \
        get<0>(get<4>(goal_inter)) = jt_desired[4]; \
        get<0>(get<5>(goal_inter)) = jt_desired[5]; \
        get<0>(get<6>(goal_inter)) = jt_desired[6]; \
        goal_point = r_info.manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef ReaK::pp::frame_tracer_3D< \
          rldk_2_type, jspace_rl_2_type, ReaK::pp::identity_topo_map, ReaK::pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          rldk_2_type(r_info.manip_kin_mdl, r_info.manip_jt_limits, normal_jt_space), \
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
        workspace_planner.solve_path(); \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(get<0>(start_inter)); \
        r_info.builder.arm_joint_1_coord->q = get<0>(get<1>(start_inter)); \
        r_info.builder.arm_joint_2_coord->q = get<0>(get<2>(start_inter)); \
        r_info.builder.arm_joint_3_coord->q = get<0>(get<3>(start_inter)); \
        r_info.builder.arm_joint_4_coord->q = get<0>(get<4>(start_inter)); \
        r_info.builder.arm_joint_5_coord->q = get<0>(get<5>(start_inter)); \
        r_info.builder.arm_joint_6_coord->q = get<0>(get<6>(start_inter)); \
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
        workspace_planner.solve_path(); \
         \
        mg_sep = workspace_planner.get_reporter().get_motion_graph_tracer(r_info.builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < workspace_planner.get_reporter().get_solution_count(); ++i) { \
          sol_seps.push_back(workspace_planner.get_reporter().get_solution_tracer(r_info.builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        r_info.builder.track_joint_coord->q = get<0>(get<0>(start_inter)); \
        r_info.builder.arm_joint_1_coord->q = get<0>(get<1>(start_inter)); \
        r_info.builder.arm_joint_2_coord->q = get<0>(get<2>(start_inter)); \
        r_info.builder.arm_joint_3_coord->q = get<0>(get<3>(start_inter)); \
        r_info.builder.arm_joint_4_coord->q = get<0>(get<4>(start_inter)); \
        r_info.builder.arm_joint_5_coord->q = get<0>(get<5>(start_inter)); \
        r_info.builder.arm_joint_6_coord->q = get<0>(get<6>(start_inter)); \
        r_info.kin_chain->doMotion();
  
  
  
  
  switch(configs.planning_algo_selection->currentIndex()) {
    case 0:  // RRT
    {
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(workspace_0_linear_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_0_linear_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_linear_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_1_linear_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_linear_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_2_linear_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_cubic_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_1_cubic_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_cubic_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_2_cubic_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_quintic_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_2_quintic_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(workspace_2_sap_type)
      };
      
    }; break;
    case 1:  // RRT*
    {
      
      if((space_order == 0) && (interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(workspace_0_linear_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_0_linear_type)
      } else if((space_order == 1) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_linear_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_1_linear_type)
      } else if((space_order == 2) && (interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_linear_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_2_linear_type)
      } else if((space_order == 1) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_cubic_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_1_cubic_type)
      } else if((space_order == 2) && (interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_cubic_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_2_cubic_type)
      } else if((space_order == 2) && (interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_quintic_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_2_quintic_type)
      } else if((space_order == 1) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(workspace_1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_1_svp_type)
      } else if((space_order == 2) && (interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_2_svp_type)
      } else if((space_order == 2) && (interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(workspace_2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(workspace_2_sap_type)
      };
      
    }; break;
    case 2:  // RRG
    { };// break;
    case 3:  // PRM
    { };// break;
    case 4:  // PRM*
    { };// break;
    case 5:  // FADPRM
    { };// break;
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
  QMessageBox::information(this,
                  "Action captured!",
                  "You have just pressed the animate best solution button!",
                  QMessageBox::Ok);
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








