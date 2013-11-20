
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

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

#include "shapes/oi_scene_graph.hpp"
#include "proximity/proxy_query_model.hpp"

#include "serialization/archiver_factory.hpp"

#include "CRS_planner_data.hpp"

#include "optimization/optim_exceptions.hpp"

#include "base/chrono_incl.hpp"



using namespace ReaK;




static QString last_used_path;





void CRSPlannerGUI_animate_bestsol_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);
  
  static shared_ptr< CRS_sol_anim_data::trajectory_type > manip_traj     = p->sol_anim.trajectory;
  static CRS_sol_anim_data::trajectory_type::point_time_iterator cur_pit = manip_traj->begin_time_travel();
  static ReaKaux::chrono::high_resolution_clock::time_point animation_start = ReaKaux::chrono::high_resolution_clock::now();
  
  if(!manip_traj) {
    manip_traj = p->sol_anim.trajectory;
    cur_pit = manip_traj->begin_time_travel();
    animation_start = ReaKaux::chrono::high_resolution_clock::now();
  };
  if( (p->sol_anim.enabled) && ( (*cur_pit).time < manip_traj->get_end_time() ) ) {
    if( (*cur_pit).time <= 0.001 * (ReaKaux::chrono::duration_cast<ReaKaux::chrono::milliseconds>(ReaKaux::chrono::high_resolution_clock::now() - animation_start)).count() ) {
      cur_pit += 0.1;
      p->ct_config.sceneData.chaser_kin_model->setJointPositions(vect_n<double>((*cur_pit).pt));
      p->ct_config.sceneData.chaser_kin_model->doDirectMotion();
    };
  } else {
    p->sol_anim.animation_timer->unschedule();
    animation_start = ReaKaux::chrono::high_resolution_clock::now();
    manip_traj.reset();
  };
};

void CRSPlannerGUI::startSolutionAnimation() {
  
  if( space_config.is_temporal ) {
    startCompleteAnimation();
    return;
  };
  
  if( !sol_anim.trajectory || !ct_config.sceneData.chaser_kin_model ) {
    QMessageBox::critical(this,
                  "Animation Error!",
                  "The best-solution trajectory is missing (not loaded or erroneous)! Cannot animate chaser!",
                  QMessageBox::Ok);
    return;
  };
  sol_anim.enabled = true;
  sol_anim.animation_timer->schedule();
};


void CRSPlannerGUI::stopSolutionAnimation() {
  
  if( space_config.is_temporal ) {
    stopCompleteAnimation();
    return;
  };
  
  sol_anim.enabled = false;
};





void CRSPlannerGUI_animate_target_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);
  
  static shared_ptr< CRS_target_anim_data::trajectory_type > target_traj    = p->target_anim.trajectory;
  static CRS_target_anim_data::trajectory_type::point_time_iterator cur_pit = target_traj->begin_time_travel();
  static ReaKaux::chrono::high_resolution_clock::time_point animation_start = ReaKaux::chrono::high_resolution_clock::now();
  
  if(!target_traj) {
    target_traj = p->target_anim.trajectory;
    cur_pit = target_traj->begin_time_travel();
    animation_start = ReaKaux::chrono::high_resolution_clock::now();
  };
  if( (p->target_anim.enabled) && ( (*cur_pit).time < target_traj->get_end_time() ) ) {
    if( (*cur_pit).time <= 0.001 * (ReaKaux::chrono::duration_cast<ReaKaux::chrono::milliseconds>(ReaKaux::chrono::high_resolution_clock::now() - animation_start)).count() ) {
      cur_pit += 0.1;
      *(p->ct_config.sceneData.target_kin_model->getFrame3D(0)) = get_frame_3D((*cur_pit).pt); 
      p->ct_config.sceneData.target_kin_model->doDirectMotion();
      
      if( p->ct_config.sceneData.chaser_kin_model && p->ct_interact.isIKEnabled() && !(p->sol_anim.enabled)) {
        try {
          frame_3D<double> tf = p->ct_config.sceneData.target_frame->getFrameRelativeTo(p->ct_config.sceneData.chaser_kin_model->getDependentFrame3D(0)->mFrame);
          p->ct_config.sceneData.chaser_kin_model->getDependentFrame3D(0)->mFrame->addBefore(tf);
          p->ct_config.sceneData.chaser_kin_model->doInverseMotion();
        } catch( optim::infeasible_problem& e ) { RK_UNUSED(e); };
        p->ct_config.sceneData.chaser_kin_model->doDirectMotion();
      };
    
    };
  } else {
    p->target_anim.animation_timer->unschedule();
    animation_start = ReaKaux::chrono::high_resolution_clock::now();
    target_traj.reset();
  };
};

void CRSPlannerGUI::startTargetAnimation() {
  if( !target_anim.trajectory || !ct_config.sceneData.target_kin_model ) {
    QMessageBox::critical(this,
                  "Animation Error!",
                  "The target trajectory is missing (not loaded or erroneous)! Cannot animate target!",
                  QMessageBox::Ok);
    return;
  };
  target_anim.enabled = true;
  target_anim.animation_timer->schedule();
};

void CRSPlannerGUI::stopTargetAnimation() {
  target_anim.enabled = false;
};

void CRSPlannerGUI::loadTargetTrajectory(QString fileName) {
  *(serialization::open_iarchive(fileName.toStdString()))
    & RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory", target_anim.trajectory);
};




void CRSPlannerGUI::startCompleteAnimation() {
  
  if( !sol_anim.trajectory || !target_anim.trajectory || !ct_config.sceneData.chaser_kin_model || !ct_config.sceneData.target_kin_model ) {
    QMessageBox::critical(this,
                  "Animation Error!",
                  "One of the trajectories is missing (not loaded or erroneous)! Cannot animate chaser and target!",
                  QMessageBox::Ok);
    return;
  };
  sol_anim.enabled = true;
  target_anim.enabled = true;
  sol_anim.animation_timer->schedule();
  target_anim.animation_timer->schedule();
};

void CRSPlannerGUI::stopCompleteAnimation() {
  sol_anim.enabled = false;
  target_anim.enabled = false;
};



void CRSPlannerGUI::runPlanner() {
  if( space_config.is_temporal ) {
    executeDynamicPlanner();
  } else {
    executePlanner();
  };
};



CRSPlannerGUI::CRSPlannerGUI( QWidget * parent, Qt::WindowFlags flags ) : 
    QMainWindow(parent,flags),
    view3d_menu(this),
    ct_config(&view3d_menu, this),
    ct_interact(&ct_config.sceneData, this),
    space_config(this), 
    plan_alg_config(this) {
  setupUi(this);
  
  
  addDockWidget(Qt::LeftDockWidgetArea, &ct_config);
  addDockWidget(Qt::LeftDockWidgetArea, &ct_interact);
  addDockWidget(Qt::LeftDockWidgetArea, &space_config);
  addDockWidget(Qt::LeftDockWidgetArea, &plan_alg_config);
  
  
  tabifyDockWidget(&space_config, &ct_config);
  tabifyDockWidget(&ct_config, &plan_alg_config);
  tabifyDockWidget(&plan_alg_config, &ct_interact);
  
  
  connect(actionStart_Sol_Anim, SIGNAL(triggered()), this, SLOT(startSolutionAnimation()));
  connect(actionStop_Sol_Anim, SIGNAL(triggered()), this, SLOT(stopSolutionAnimation()));
  connect(&ct_interact, SIGNAL(onStartTargetAnimation()), this, SLOT(startTargetAnimation()));
  connect(&ct_interact, SIGNAL(onStopTargetAnimation()), this, SLOT(stopTargetAnimation()));
  connect(&ct_interact, SIGNAL(onLoadTargetTrajectory(QString)), this, SLOT(loadTargetTrajectory(QString)));
  
  connect(actionRun_Planner, SIGNAL(triggered()), this, SLOT(runPlanner()));
//   connect(configs.actionExecutePlanner, SIGNAL(triggered()), this, SLOT(executePlanner()));
//   connect(configs.actionExecuteDynamicPlanner, SIGNAL(triggered()), this, SLOT(executeDynamicPlanner()));
  
  SoQt::init(this->centralwidget);
  
  menubar->addMenu(&view3d_menu);
  view3d_menu.setViewer(new SoQtExaminerViewer(this->centralwidget));
  
  view3d_menu.getGeometryGroup("Chaser Geometry",true);
  view3d_menu.getGeometryGroup("Chaser KTE Chain",false);
  view3d_menu.getGeometryGroup("Target Geometry",true);
  view3d_menu.getGeometryGroup("Environment",true);
  view3d_menu.getDisplayGroup("Motion-Graph",true);
  view3d_menu.getDisplayGroup("Solution(s)",true);
  
  sol_anim.animation_timer    = new SoTimerSensor(CRSPlannerGUI_animate_bestsol_trajectory, this);
  target_anim.animation_timer = new SoTimerSensor(CRSPlannerGUI_animate_target_trajectory, this);
  
};


CRSPlannerGUI::~CRSPlannerGUI() {
  
  delete target_anim.animation_timer;
  delete sol_anim.animation_timer;
  
  
  view3d_menu.removeGeometryGroup("Chaser Geometry");
  view3d_menu.removeGeometryGroup("Chaser KTE Chain");
  view3d_menu.removeGeometryGroup("Target Geometry");
  view3d_menu.removeGeometryGroup("Environment");
  view3d_menu.removeDisplayGroup("Motion-Graph");
  view3d_menu.removeDisplayGroup("Solution(s)");
  
  view3d_menu.setViewer(NULL);
  
  SoQt::done();
  
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








