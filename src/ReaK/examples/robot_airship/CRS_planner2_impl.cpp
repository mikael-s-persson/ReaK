
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

#include "CRS_planner2_impl.hpp"


#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>
#include <QMainWindow>
#include <QDir>



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

#include "serialization/archiver_factory.hpp"


#include "CRS_planners_utility.hpp"

#include "CRS_planner_data.hpp"

#include "optimization/optim_exceptions.hpp"


#include <chrono>


using namespace ReaK;




static QString last_used_path;





void CRSPlannerGUI_animate_bestsol_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);
  
  static std::size_t animation_progress = 0;
  
  if( (p->sol_anim->enabled) && ( animation_progress < p->sol_anim->bestsol_trajectory.size() ) ) {
    if(std::chrono::high_resolution_clock::now() - p->sol_anim->animation_last_render >= std::chrono::milliseconds(100)) {
      
      p->scene_data.chaser_kin_model->setJointPositions(vect_n<double>(p->sol_anim->bestsol_trajectory[animation_progress]));
      p->scene_data.chaser_kin_model->doDirectMotion();
      
      animation_progress++;
      p->sol_anim->animation_last_render = std::chrono::high_resolution_clock::now();
    };
  } else {
    p->sol_anim->animation_timer->unschedule();
    animation_progress = 0;
    p->sol_anim->animation_last_render = std::chrono::high_resolution_clock::now();
  };
};

void CRSPlannerGUI::startSolutionAnimation() {
  sol_anim->animation_last_render = std::chrono::high_resolution_clock::now();
  sol_anim->enabled = true;
  sol_anim->animation_timer->schedule();
};







void CRSPlannerGUI_animate_target_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);
  
  static shared_ptr< sat_traj_type > target_traj = p->target_anim->target_trajectory;
  static sat_traj_type::point_time_iterator cur_pit = target_traj->begin_time_travel();
  
  if(!target_traj) {
    target_traj = p->target_anim->target_trajectory;
    cur_pit = target_traj->begin_time_travel();
  };
  if( (p->target_anim->enabled) && ( cur_pit->time < target_traj->get_end_time() ) ) {
    if(std::chrono::high_resolution_clock::now() - p->target_anim->target_anim_last_render >= std::chrono::milliseconds(100)) {
      p->target_anim->target_anim_last_render = std::chrono::high_resolution_clock::now();
      cur_pit += 0.1;
      *(p->scene_data.target_kin_model->getFrame3D(0)) = get_frame_3D(cur_pit->pt); 
      p->scene_data.target_kin_model->doDirectMotion();
    };
  } else {
    p->target_anim->target_anim_timer->unschedule();
    p->target_anim->target_anim_last_render = std::chrono::high_resolution_clock::now();
    target_traj.reset();
  };
};

void CRSPlannerGUI::startTargetAnimation() {
  if( !target_anim->target_trajectory ) {
    QMessageBox::critical(this,
                  "Animation Error!",
                  "The target trajectory is missing (not loaded or erroneous)! Cannot animate target!",
                  QMessageBox::Ok);
    return;
  };
  target_anim->target_anim_last_render = std::chrono::high_resolution_clock::now();
  target_anim->enabled = true;
  target_anim->target_anim_timer->schedule();
};

void CRSPlannerGUI::loadTargetTrajectory() {
  
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Open Target Trajectory..."),
    last_used_path,
    tr("SE(3) Trajectories (*.rkx *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  QFileInfo fileInf(fileName);
  
  last_used_path = fileInf.absolutePath();
  
  *(serialization::open_iarchive(fileName.toStdString()))
    & RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory", target_anim->target_trajectory);
  
  if( target_anim->target_trajectory ) {
    configs.traj_filename_edit->setText(fileInf.baseName());
  } else {
    configs.traj_filename_edit->setText(tr("ERROR!"));
  };
};



CRSPlannerGUI::CRSPlannerGUI( QWidget * parent, Qt::WindowFlags flags ) : QMainWindow(parent,flags),
                                                                          configs() {
  setupUi(this);
  
  scene_data.load_chaser("models/CRS_A465.model.rkx");
  scene_data.load_target("models/airship3D.model.rkx");
  scene_data.load_environment("models/MD148_lab.geom.rkx");
  
  draw_data = new CRS_coin_nodes();
  sol_anim = new CRS_sol_anim_data();
  target_anim = new CRS_target_anim_data();
  
  
  
  configs.setupUi(this->config_dock->widget());
  connect(configs.actionStart_Robot, SIGNAL(triggered()), this, SLOT(startSolutionAnimation()));
  connect(configs.actionAnimateTarget, SIGNAL(triggered()), this, SLOT(startTargetAnimation()));
  connect(configs.actionExecutePlanner, SIGNAL(triggered()), this, SLOT(executePlanner()));
  connect(configs.actionJointChange, SIGNAL(triggered()), this, SLOT(onJointChange()));
  connect(configs.actionTargetChange, SIGNAL(triggered()), this, SLOT(onTargetChange()));
  
  connect(configs.actionLoadTargetTrajectory, SIGNAL(triggered()), this, SLOT(loadTargetTrajectory()));
  
  connect(configs.actionRobotVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotVisible()));
  connect(configs.actionRobotKinVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotKinVisible()));
  connect(configs.actionTargetVisibleToggle, SIGNAL(triggered()), this, SLOT(onTargetVisible()));
  connect(configs.actionEnvVisibleToggle, SIGNAL(triggered()), this, SLOT(onEnvVisible()));
  connect(configs.actionProxyVisibleToggle, SIGNAL(triggered()), this, SLOT(onProxyVisible()));
  connect(configs.actionMGVisibleToggle, SIGNAL(triggered()), this, SLOT(onMGVisible()));
  connect(configs.actionSolutionsVisibleToggle, SIGNAL(triggered()), this, SLOT(onSolutionsVisible()));
  
  connect(configs.actionUpdateAvailOptions, SIGNAL(triggered()), this, SLOT(onUpdateAvailableOptions()));
  
  
  
  connect(actionLoad_Positions, SIGNAL(triggered()), this, SLOT(loadPositions()));
  connect(actionSave_Positions, SIGNAL(triggered()), this, SLOT(savePositions()));
  connect(actionLoad_Planner, SIGNAL(triggered()), this, SLOT(loadPlannerConfig()));
  connect(actionSave_Planner, SIGNAL(triggered()), this, SLOT(savePlannerConfig()));
  
  
  
  
  sol_anim->animation_timer = new SoTimerSensor(CRSPlannerGUI_animate_bestsol_trajectory, this);
  target_anim->target_anim_timer = new SoTimerSensor(CRSPlannerGUI_animate_target_trajectory, this);
  
  
  
  SoQt::init(this->centralwidget);
  
  draw_data->sg_root = new SoSeparator;
  draw_data->sg_root->ref();
  
  
  draw_data->sw_robot_geom = new SoSwitch();
  draw_data->sg_robot_geom = new geom::oi_scene_graph();
  
  (*draw_data->sg_robot_geom) << (*scene_data.chaser_geom_model);
  double charact_length = draw_data->sg_robot_geom->computeCharacteristicLength();
  
  draw_data->sw_robot_geom->addChild(draw_data->sg_robot_geom->getSceneGraph());
  draw_data->sw_robot_geom->whichChild.setValue((configs.check_show_geom->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  draw_data->sg_root->addChild(draw_data->sw_robot_geom);
  
  
  draw_data->sw_robot_kin = new SoSwitch();
  draw_data->sg_robot_kin = new geom::oi_scene_graph();
  
  draw_data->sg_robot_kin->setCharacteristicLength(charact_length);
  (*draw_data->sg_robot_kin) << (*scene_data.chaser_kin_model->getKTEChain());
  
  draw_data->sw_robot_kin->addChild(draw_data->sg_robot_kin->getSceneGraph());
  draw_data->sw_robot_kin->whichChild.setValue((configs.check_show_kinmdl->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  draw_data->sg_root->addChild(draw_data->sw_robot_kin);
  
  
  draw_data->sw_airship_geom = new SoSwitch();
  draw_data->sg_airship_geom = new geom::oi_scene_graph();
  
  (*draw_data->sg_airship_geom) << (*scene_data.target_geom_model);
  
  draw_data->sw_airship_geom->addChild(draw_data->sg_airship_geom->getSceneGraph());
  draw_data->sw_airship_geom->whichChild.setValue((configs.check_show_target->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  draw_data->sg_root->addChild(draw_data->sw_airship_geom);
  
  
  draw_data->sw_lab_geom = new SoSwitch();
  draw_data->sg_lab_geom = new geom::oi_scene_graph();
  
  for(std::size_t i = 0; i < scene_data.env_geom_models.size(); ++i)
    (*draw_data->sg_lab_geom) << (*(scene_data.env_geom_models[i]));
  
  draw_data->sw_lab_geom->addChild(draw_data->sg_lab_geom->getSceneGraph());
  draw_data->sw_lab_geom->whichChild.setValue((configs.check_show_env->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
  draw_data->sg_root->addChild(draw_data->sw_lab_geom);
  
  
  draw_data->sw_motion_graph = new SoSwitch();
  draw_data->sw_motion_graph->whichChild.setValue((configs.check_show_motiongraph->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  draw_data->sg_root->addChild(draw_data->sw_motion_graph);
  
  draw_data->sw_solutions = new SoSwitch();
  draw_data->sw_solutions->whichChild.setValue((configs.check_show_sol->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  draw_data->sg_root->addChild(draw_data->sw_solutions);
  
  
  
  onJointChange();
  onTargetChange();
  
  
  
  plan_options.space_order = 0;
  plan_options.interp_id = 0;
  plan_options.min_travel = 0.1;
  plan_options.max_travel = 1.0;
  plan_options.planning_algo = 0;
  plan_options.max_vertices = 2000;
  plan_options.prog_interval = 500;
  plan_options.max_results = 50;
  plan_options.planning_options = 0;
  plan_options.store_policy = 0;
  plan_options.knn_method = 2;
  plan_options.init_SA_temp = 0.0;
  plan_options.init_relax = 5.0;
  updateConfigs();
  
  
  
  
  
  draw_data->eviewer = new SoQtExaminerViewer(this->centralwidget);
  draw_data->eviewer->setSceneGraph(draw_data->sg_root);
  draw_data->eviewer->show();
  
  draw_data->sg_robot_geom->enableAnchorUpdates();
  draw_data->sg_robot_kin->enableAnchorUpdates();
  draw_data->sg_airship_geom->enableAnchorUpdates();
};


CRSPlannerGUI::~CRSPlannerGUI() {
  
  draw_data->sg_robot_geom->disableAnchorUpdates();
  draw_data->sg_robot_kin->disableAnchorUpdates();
  draw_data->sg_airship_geom->disableAnchorUpdates();
  
  delete draw_data->sg_robot_geom;
  delete draw_data->sg_robot_kin;
  delete draw_data->sg_airship_geom;
  delete draw_data->sg_lab_geom;
  
  delete target_anim->target_anim_timer;
  
  delete sol_anim->animation_timer;
  
  delete draw_data->eviewer;
  draw_data->sg_root->unref();
  SoQt::done();
  
  
  
  delete draw_data;
  delete sol_anim;
  delete target_anim;
  
};


void CRSPlannerGUI::onJointChange() {
  scene_data.chaser_kin_model->setJointPositions(
    vect_n<double>(
      double(configs.track_pos->value()) * 0.001,
      double(configs.joint1_pos->value()) * 0.001,
      double(configs.joint2_pos->value()) * 0.001,
      double(configs.joint3_pos->value()) * 0.001,
      double(configs.joint4_pos->value()) * 0.001,
      double(configs.joint5_pos->value()) * 0.001,
      double(configs.joint6_pos->value()) * 0.001
    )
  );
  scene_data.chaser_kin_model->doDirectMotion();
  
  onProxyChange();
};

void CRSPlannerGUI::onTargetChange() {
  shared_ptr< frame_3D<double> > target_state = scene_data.target_kin_model->getFrame3D(0);
  target_state->Position = vect<double,3>(
    double(configs.target_x->value()) * 0.001, 
    double(configs.target_y->value()) * 0.001, 
    double(configs.target_z->value()) * 0.001);
  target_state->Quat = 
    quaternion<double>::zrot(double(configs.target_yaw->value()) * 0.001) * 
    quaternion<double>::yrot(double(configs.target_pitch->value()) * 0.001) * 
    quaternion<double>::xrot(double(configs.target_roll->value()) * 0.001);
  scene_data.target_kin_model->doDirectMotion();
  
  
  if(configs.check_enable_ik->isChecked()) {
    try {
      frame_3D<double> tf = scene_data.target_frame->getFrameRelativeTo(scene_data.chaser_kin_model->getDependentFrame3D(0)->mFrame);
      scene_data.chaser_kin_model->getDependentFrame3D(0)->mFrame->addBefore(tf);
      scene_data.chaser_kin_model->doInverseMotion();
    } catch( optim::infeasible_problem& e ) { RK_UNUSED(e); };
    scene_data.chaser_kin_model->doDirectMotion();
  };
  
  onProxyChange();
};

void CRSPlannerGUI::onRobotVisible() {
  draw_data->sw_robot_geom->whichChild.setValue((configs.check_show_geom->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onRobotKinVisible() {
  draw_data->sw_robot_kin->whichChild.setValue((configs.check_show_kinmdl->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onTargetVisible() {
  draw_data->sw_airship_geom->whichChild.setValue((configs.check_show_target->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onEnvVisible() {
  draw_data->sw_lab_geom->whichChild.setValue((configs.check_show_env->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onProxyVisible() {
  
};

void CRSPlannerGUI::onMGVisible() {
  draw_data->sw_motion_graph->whichChild.setValue((configs.check_show_motiongraph->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onSolutionsVisible() {
  draw_data->sw_solutions->whichChild.setValue((configs.check_show_sol->isChecked() ? SO_SWITCH_ALL : SO_SWITCH_NONE));
};

void CRSPlannerGUI::onProxyChange() {
  
};


void CRSPlannerGUI::onConfigsChanged() {
  
  // joint-space parameters:
  plan_options.space_order = configs.order_selection->currentIndex();
  plan_options.interp_id = configs.interp_selection->currentIndex();
  
  plan_options.min_travel = configs.min_interval_spinbox->value();
  plan_options.max_travel = configs.max_interval_spinbox->value();
  
  
  // planner parameters:
  plan_options.planning_algo = configs.planning_algo_selection->currentIndex();
  
  plan_options.max_vertices = configs.maxvertices_spinbox->value();
  plan_options.prog_interval = configs.progress_interval_spinbox->value();
  plan_options.max_results = configs.maxsolutions_spinbox->value();
  
  plan_options.planning_options = pp::UNIDIRECTIONAL_PLANNING;
  
  if(configs.check_bidir->isChecked())
    plan_options.planning_options |= pp::BIDIRECTIONAL_PLANNING;
  
  if(configs.check_lazy_collision->isChecked())
    plan_options.planning_options |= pp::LAZY_COLLISION_CHECKING;
  
  plan_options.init_SA_temp = -1.0;
  if(configs.check_voronoi_pull->isChecked()) {
    plan_options.planning_options |= pp::PLAN_WITH_VORONOI_PULL;
    plan_options.init_SA_temp = configs.init_sa_temp_spinbox->value();
    if(plan_options.init_SA_temp < 1e-6)
      plan_options.init_SA_temp = -1.0;
  };
  
  plan_options.init_relax = 0.0;
  if(configs.check_anytime_heuristic->isChecked()) {
    plan_options.planning_options |= pp::PLAN_WITH_ANYTIME_HEURISTIC;
    plan_options.init_relax = configs.init_relax_spinbox->value();
  };
  
  if(configs.check_bnb->isChecked())
    plan_options.planning_options |= pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;
  
  
  plan_options.store_policy = pp::ADJ_LIST_MOTION_GRAPH;
  if(configs.graph_storage_selection->currentIndex())
    plan_options.store_policy = pp::DVP_ADJ_LIST_MOTION_GRAPH;
  
  plan_options.knn_method = pp::LINEAR_SEARCH_KNN;
  switch(configs.KNN_method_selection->currentIndex()) {
    case 1:
      plan_options.knn_method = pp::DVP_BF2_TREE_KNN;
      break;
    case 2:
      plan_options.knn_method = pp::DVP_BF4_TREE_KNN;
      break;
    case 3:
      plan_options.knn_method = pp::DVP_COB2_TREE_KNN;
      break;
    case 4:
      plan_options.knn_method = pp::DVP_COB4_TREE_KNN;
      break;
    default:
      break;
  };
  
};


void CRSPlannerGUI::updateConfigs() {
  configs.order_selection->setCurrentIndex(plan_options.space_order);
  configs.interp_selection->setCurrentIndex(plan_options.interp_id);
  
  configs.min_interval_spinbox->setValue(plan_options.min_travel);
  configs.max_interval_spinbox->setValue(plan_options.max_travel);
  
  // planner parameters:
  configs.planning_algo_selection->setCurrentIndex(plan_options.planning_algo);
  
  configs.maxvertices_spinbox->setValue(plan_options.max_vertices);
  configs.progress_interval_spinbox->setValue(plan_options.prog_interval);
  configs.maxsolutions_spinbox->setValue(plan_options.max_results);
  
  if(plan_options.planning_options & pp::BIDIRECTIONAL_PLANNING)
    configs.check_bidir->setChecked(true);
  else
    configs.check_bidir->setChecked(false);
  
  
  if(plan_options.planning_options & pp::LAZY_COLLISION_CHECKING)
    configs.check_lazy_collision->setChecked(true);
  else
    configs.check_lazy_collision->setChecked(false);
  
  
  if(plan_options.planning_options & pp::PLAN_WITH_VORONOI_PULL) {
    configs.init_sa_temp_spinbox->setValue(plan_options.init_SA_temp);
    configs.check_voronoi_pull->setChecked(true);
  } else {
    configs.init_sa_temp_spinbox->setValue(0.0);
    configs.check_voronoi_pull->setChecked(false);
  };
  
  if(plan_options.planning_options & pp::PLAN_WITH_ANYTIME_HEURISTIC) {
    configs.init_relax_spinbox->setValue(plan_options.init_relax);
    configs.check_anytime_heuristic->setChecked(true);
  } else {
    configs.init_relax_spinbox->setValue(0.0);
    configs.check_anytime_heuristic->setChecked(false);
  };
  
  if(plan_options.planning_options & pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG)
    configs.check_bnb->setChecked(true);
  else
    configs.check_bnb->setChecked(false);
  
  
  if( plan_options.store_policy == pp::DVP_ADJ_LIST_MOTION_GRAPH ) {
    configs.graph_storage_selection->setCurrentIndex(1);
  } else {
    configs.graph_storage_selection->setCurrentIndex(0);
  };
  
  switch(plan_options.knn_method) {
    case pp::DVP_BF2_TREE_KNN:
      configs.KNN_method_selection->setCurrentIndex(1);
      break;
    case pp::DVP_BF4_TREE_KNN:
      configs.KNN_method_selection->setCurrentIndex(2);
      break;
    case pp::DVP_COB2_TREE_KNN:
      configs.KNN_method_selection->setCurrentIndex(3);
      break;
    case pp::DVP_COB4_TREE_KNN:
      configs.KNN_method_selection->setCurrentIndex(4);
      break;
    default:
      configs.KNN_method_selection->setCurrentIndex(0);
      break;
  };
  
  onUpdateAvailableOptions();
  
};



void CRSPlannerGUI::onUpdateAvailableOptions() {
  int plan_alg = configs.planning_algo_selection->currentIndex();
  
  switch(plan_alg) {
    case 1:  // RRT*
      configs.check_lazy_collision->setEnabled(false);
      configs.check_lazy_collision->setChecked(true);
      break;
    case 3:  // SBA*
      configs.check_lazy_collision->setEnabled(true);
      configs.check_lazy_collision->setChecked(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      configs.check_lazy_collision->setEnabled(false);
      configs.check_lazy_collision->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 0:  // RRT
      configs.check_bidir->setEnabled(true);
      break;
    case 1:  // RRT*
    case 2:  // PRM
    case 3:  // SBA*
    case 4:  // FADPRM
    default:
      configs.check_bidir->setEnabled(false);
      configs.check_bidir->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 3:  // SBA*
      configs.check_voronoi_pull->setEnabled(true);
      configs.check_anytime_heuristic->setEnabled(true);
      break;
    case 0:  // RRT
    case 1:  // RRT*
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      configs.check_voronoi_pull->setEnabled(false);
      configs.check_voronoi_pull->setChecked(false);
      configs.check_anytime_heuristic->setEnabled(false);
      configs.check_anytime_heuristic->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 1:  // RRT*
    case 3:  // SBA*
      configs.check_bnb->setEnabled(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      configs.check_bnb->setEnabled(false);
      configs.check_bnb->setChecked(false);
      break;
  };
  
};




void CRSPlannerGUI::savePositions() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Positions Record..."), last_used_path,
    tr("Robot-Airship Positions Record (*.rapos.rkx *.rapos.rkb *.rapos.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  vect<double,7> robot_joint_positions;
  robot_joint_positions[0] = double(configs.track_pos->value()) * 0.001;
  robot_joint_positions[1] = double(configs.joint1_pos->value()) * 0.001;
  robot_joint_positions[2] = double(configs.joint2_pos->value()) * 0.001;
  robot_joint_positions[3] = double(configs.joint3_pos->value()) * 0.001;
  robot_joint_positions[4] = double(configs.joint4_pos->value()) * 0.001;
  robot_joint_positions[5] = double(configs.joint5_pos->value()) * 0.001;
  robot_joint_positions[6] = double(configs.joint6_pos->value()) * 0.001;
  
  vect<double,3> airship_position;
  airship_position[0] = double(configs.target_x->value()) * 0.001;
  airship_position[1] = double(configs.target_y->value()) * 0.001;
  airship_position[2] = double(configs.target_z->value()) * 0.001;
  
  euler_angles_TB<double> ea;
  ea.yaw() = double(configs.target_yaw->value()) * 0.001;
  ea.pitch() = double(configs.target_pitch->value()) * 0.001;
  ea.roll() = double(configs.target_roll->value()) * 0.001;
  quaternion<double> airship_quaternion = ea.getQuaternion();
  
  
  try {
    *(serialization::open_oarchive(fileName.toStdString()))
      & RK_SERIAL_SAVE_WITH_NAME(robot_joint_positions)
      & RK_SERIAL_SAVE_WITH_NAME(airship_position)
      & RK_SERIAL_SAVE_WITH_NAME(airship_quaternion);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void CRSPlannerGUI::loadPositions() {
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Open Positions Record..."),
    last_used_path,
    tr("Robot-Airship Positions Record (*.rapos.rkx *.rapos.rkb *.rapos.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  vect<double,7> robot_joint_positions;
  vect<double,3> airship_position;
  quaternion<double> airship_quaternion;
  
  try {
    *(serialization::open_iarchive(fileName.toStdString()))
      & RK_SERIAL_LOAD_WITH_NAME(robot_joint_positions)
      & RK_SERIAL_LOAD_WITH_NAME(airship_position)
      & RK_SERIAL_LOAD_WITH_NAME(airship_quaternion);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  configs.track_pos->setValue(int(1000.0 * robot_joint_positions[0]));
  configs.joint1_pos->setValue(int(1000.0 * robot_joint_positions[1]));
  configs.joint2_pos->setValue(int(1000.0 * robot_joint_positions[2]));
  configs.joint3_pos->setValue(int(1000.0 * robot_joint_positions[3]));
  configs.joint4_pos->setValue(int(1000.0 * robot_joint_positions[4]));
  configs.joint5_pos->setValue(int(1000.0 * robot_joint_positions[5]));
  configs.joint6_pos->setValue(int(1000.0 * robot_joint_positions[6])); 
  //onJointChange();
  
  configs.target_x->setValue(int(1000.0 * airship_position[0]));
  configs.target_y->setValue(int(1000.0 * airship_position[1]));
  configs.target_z->setValue(int(1000.0 * airship_position[2]));
  euler_angles_TB<double> ea = euler_angles_TB<double>(airship_quaternion);
  configs.target_yaw->setValue(int(1000.0 * ea.yaw()));
  configs.target_pitch->setValue(int(1000.0 * ea.pitch()));
  configs.target_roll->setValue(int(1000.0 * ea.roll()));
  //onTargetChange();
  
};

void CRSPlannerGUI::savePlannerConfig() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Planner Configurations..."), last_used_path,
    tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb *.raplan.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    (*serialization::open_oarchive(fileName.toStdString())) << plan_options;
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void CRSPlannerGUI::loadPlannerConfig() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Planner Configurations..."), last_used_path,
    tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb *.raplan.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    (*serialization::open_iarchive(fileName.toStdString())) >> plan_options;
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
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








