
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

#include "CRS_planner_impl.hpp"

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

#include "interpolation/discrete_point_trajectory.hpp"
#include "path_planning/trajectory_base.hpp"

#include "CRS_planner_data.hpp"

#include "optimization/optim_exceptions.hpp"

#include "base/chrono_incl.hpp"


#include <boost/asio.hpp>
#include <boost/bind.hpp>

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
  
  static ReaKaux::chrono::high_resolution_clock::time_point animation_start = ReaKaux::chrono::high_resolution_clock::now();
  
  static shared_ptr< CRS_target_anim_data::trajectory_type > target_traj;
  
  static shared_ptr< CRS_target_anim_data::trajectory_type::point_time_iterator > cur_pit_ptr;
  
  if(!target_traj) {
    animation_start = ReaKaux::chrono::high_resolution_clock::now();
    target_traj = p->target_anim.trajectory;
    if( target_traj ) 
      cur_pit_ptr.reset(new CRS_target_anim_data::trajectory_type::point_time_iterator(target_traj->begin_time_travel()));
    else
      return;
  };
  
  if( target_traj && p->target_anim.enabled && ( (*(*cur_pit_ptr)).time < target_traj->get_end_time() ) ) {
    CRS_target_anim_data::trajectory_type::point_time_iterator& cur_pit = *cur_pit_ptr;
    if( (*cur_pit).time <= 0.001 * (ReaKaux::chrono::duration_cast<ReaKaux::chrono::milliseconds>(ReaKaux::chrono::high_resolution_clock::now() - animation_start)).count() ) {
      cur_pit += 0.1;
      *(p->ct_config.sceneData.target_kin_model->getFrame3D(0)) = get_frame_3D((*cur_pit).pt); 
      p->ct_config.sceneData.target_kin_model->doDirectMotion();
      
      if( p->ct_config.sceneData.chaser_kin_model && p->ct_interact.isIKEnabled() && !(p->sol_anim.enabled)) {
        try {
          *(p->ct_config.sceneData.chaser_kin_model->getDependentFrame3D(0)->mFrame) = *(p->ct_config.sceneData.target_frame);
          p->ct_config.sceneData.chaser_kin_model->doInverseMotion();
        } catch( optim::infeasible_problem& e ) { RK_UNUSED(e); };
        p->ct_config.sceneData.chaser_kin_model->doDirectMotion();
      };
    };
  } else {
    p->target_anim.animation_timer->unschedule();
    cur_pit_ptr.reset();
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
  typedef pp::se3_1st_order_topology<double>::type sat_state_space_type;
  typedef pp::temporal_space<sat_state_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_space_type;
  typedef pp::discrete_point_trajectory< temporal_space_type > trajectory_type;
  
  typedef pp::trajectory_wrapper<trajectory_type> wrapped_traj_type;
  
  shared_ptr< wrapped_traj_type > tmp_traj(new wrapped_traj_type());
  
  *(serialization::open_iarchive(fileName.toStdString()))
    & RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory", tmp_traj->get_underlying_trajectory());
  target_anim.trajectory = tmp_traj;
  
//   *(serialization::open_iarchive(fileName.toStdString()))
//     & RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory", target_anim.trajectory);
};




void CRSPlannerGUI::startCompleteAnimation() {
  
  if( !sol_anim.trajectory || !target_anim.trajectory 
      || !ct_config.sceneData.chaser_kin_model || !ct_config.sceneData.target_kin_model ) {
    QMessageBox::critical(this,
                  "Animation Error!",
                  "One of the trajectories is missing (not loaded or erroneous)! Cannot animate chaser and target!",
                  QMessageBox::Ok);
    return;
  };
  sol_anim.enabled = true;
  sol_anim.animation_timer->schedule();
  target_anim.enabled = true;
  target_anim.animation_timer->schedule();
};

void CRSPlannerGUI::stopCompleteAnimation() {
  sol_anim.enabled = false;
  target_anim.enabled = false;
};



void CRSPlannerGUI::runPlanner() {
  try {
    if( space_config.is_temporal ) {
      executeDynamicPlanner();
    } else {
      executePlanner();
    };
  } catch(std::exception& e) {
    std::stringstream ss;
    ss << "An exception was raised during the planning:\nwhat(): " << e.what();
    QMessageBox::critical(this,
                  "Motion-Planning Error!",
                  QString::fromStdString(ss.str()),
                  QMessageBox::Ok);
  };
};


void CRSPlannerGUI::onShowRunDialog() {
  run_dialog.setModal(false);
  run_dialog.show();
};

void CRSPlannerGUI::threadedPlanningFunction(int mode) {
  
  if((mode == 2) || (mode == 3)) {
    run_dialog.publishConsoleMessage("Starting state estimation...\n");
    target_anim.trajectory.reset();
    target_pred_config.startStatePrediction();
    if(!target_anim.trajectory) {
      QMessageBox::critical(this,
                            "State-Prediction Error!",
                            "The live state-estimation of the target failed to produce a viable predicted trajectory!",
                            QMessageBox::Ok);
      run_dialog.onReset();
      return;
    };
    run_dialog.publishConsoleMessage("State estimation done!\nSwitched to predicted target trajectory!\n");
  };
  
  run_dialog.onInitializationDone();
  
  try {
    if( mode > 0 ) {
      executeDynamicPlanner();
    } else {
      executePlanner();
    };
  } catch(std::exception& e) {
    std::stringstream ss;
    ss << "An exception was raised during the planning:\nwhat(): " << e.what();
    QMessageBox::critical(this,
                  "Motion-Planning Error!",
                  QString::fromStdString(ss.str()),
                  QMessageBox::Ok);
    run_dialog.onReset();
  };
  
  run_dialog.onPlanningDone();
  
  if(sol_anim.trajectory)
    run_dialog.onLaunchOpportunity();
  
//   onInitializationDone();
//   onPlanningDone();
//   onLaunchOpportunity();
//   onLaunchStarted();
//   onCaptureReached();
//   onReset();
  
};

void CRSPlannerGUI::onStartPlanning(int mode) {
  if( (mode > 0) && !space_config.is_temporal ) {
    QMessageBox::critical(this,
                          "Motion-Planning Error!",
                          "Trying to run a dynamic planner with a non-temporal planning space! Please configure a temporal planning space!",
                          QMessageBox::Ok);
    run_dialog.onReset();
    return;
  };
  
//   mode == 0;  // static-sim (executePlanner)
//   mode == 1;  // dynamic-sim (loaded target-trajectory + executeDynamicPlanner + startCompleteAnimation)
//   mode == 2;  // live-sim (predicted target-trajectory + executeDynamicPlanner + startCompleteAnimation)
//   mode == 3;  // live-run (predicted target-trajectory + executeDynamicPlanner + startCompleteAnimation + executeSolTrajectory)
  
  if(planning_thr) {
    // TODO Signal the stop.
    planning_thr->join();
    delete planning_thr;
    planning_thr = NULL;
  };
  
  planning_thr = new ReaKaux::thread(boost::bind(&CRSPlannerGUI::threadedPlanningFunction, this, mode));
  
};

void CRSPlannerGUI::onStopPlanning() {
  
  target_pred_config.stopStatePrediction();
  
  if(planning_thr) {
    // TODO Signal the stop.
    planning_thr->join();
    delete planning_thr;
    planning_thr = NULL;
  };
  
};

void CRSPlannerGUI::executeSolutionTrajectory() {
  
  // Setup the UDP streaming to the robot (FIXME remove the hard-coded address / port)
  boost::asio::io_service io_service;
  boost::asio::ip::udp::endpoint endpoint(boost::asio::ip::address_v4::from_string("192.168.0.3"), 17050);
  boost::asio::ip::udp::socket socket((io_service));
  boost::asio::basic_streambuf<> udp_buf;
  socket.open(boost::asio::ip::udp::v4());
  socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
  std::ostream udp_buf_stream(&(udp_buf));
  
  shared_ptr< CRS_sol_anim_data::trajectory_type > manip_traj = sol_anim.trajectory;
  CRS_sol_anim_data::trajectory_type::point_time_iterator cur_pit = manip_traj->begin_time_travel();
  
  // Synchronize the trajectory to the current time (from target trajectory).
  cur_pit += current_target_anim_time - (*cur_pit).time;
  
  ReaKaux::chrono::high_resolution_clock::time_point exec_start = ReaKaux::chrono::high_resolution_clock::now();
  
  // UDP output for the robot:
  vect<double,7> cur_pt = (*cur_pit).pt;
  for(std::size_t i = 0; i < 7; ++i)
    udp_buf_stream.write(reinterpret_cast<char*>(&cur_pt[i]),sizeof(double));
  std::size_t len = socket.send_to(udp_buf.data(), endpoint);
  udp_buf.consume(len);
  
  
  while( exec_robot_thr && ( (*cur_pit).time < manip_traj->get_end_time() ) ) {
    cur_pit += 0.001;
    ReaKaux::this_thread::sleep_until(exec_start + ReaKaux::chrono::microseconds(1000));
    
    //UDP output for the robot:
    cur_pt = (*cur_pit).pt;
    for(std::size_t i = 0; i < 7; ++i)
      udp_buf_stream.write(reinterpret_cast<char*>(&cur_pt[i]),sizeof(double));
    len = socket.send_to(udp_buf.data(), endpoint);
    udp_buf.consume(len);
  };
  //UDP output for the robot:
  double terminating_value = -100.0;
  for(std::size_t i = 0; i < 7; ++i)
    udp_buf_stream.write(reinterpret_cast<char*>(&terminating_value),sizeof(double));
  len = socket.send_to(udp_buf.data(), endpoint);
  udp_buf.consume(len);
  
  if( (*cur_pit).time >= manip_traj->get_end_time() ) {
    run_dialog.onCaptureReached();
  };
  
  exec_robot_enabled = false;
  
};



void CRSPlannerGUI::onLaunch(int mode) {
  
  if(this->exec_robot_thr) {
    this->exec_robot_enabled = false;
    this->exec_robot_thr->join();
    delete this->exec_robot_thr;
    this->exec_robot_thr = NULL;
  };
  
  if( mode == 3 ) {
    this->exec_robot_enabled = true;
    this->exec_robot_thr = new ReaKaux::thread(boost::bind(&CRSPlannerGUI::executeSolutionTrajectory, this));
  };
  
  if(mode == 0)
    startSolutionAnimation();
  else
    startCompleteAnimation();
  
};

void CRSPlannerGUI::onAbort() {
  
  // first, stop the robot
  if(this->exec_robot_thr) {
    this->exec_robot_enabled = false;
    this->exec_robot_thr->join();
    delete this->exec_robot_thr;
    this->exec_robot_thr = NULL;
  };
  
  // then, stop the state prediction and/or planning
  this->onStopPlanning();
  
};




CRSPlannerGUI::CRSPlannerGUI( QWidget * parent, Qt::WindowFlags flags ) : 
    QMainWindow(parent,flags),
    view3d_menu(this),
    ct_config(&view3d_menu, this),
    ct_interact(&ct_config.sceneData, this),
    space_config(this), 
    plan_alg_config(this),
    target_pred_config(&target_anim, &current_target_anim_time, this),
    run_dialog(this),
    planning_thr(NULL),
    exec_robot_thr(NULL) {
  setupUi(this);
  
  
  addDockWidget(Qt::LeftDockWidgetArea, &ct_config);
  addDockWidget(Qt::LeftDockWidgetArea, &ct_interact);
  addDockWidget(Qt::LeftDockWidgetArea, &space_config);
  addDockWidget(Qt::LeftDockWidgetArea, &plan_alg_config);
  addDockWidget(Qt::LeftDockWidgetArea, &target_pred_config);
  
  
  tabifyDockWidget(&space_config, &ct_config);
  tabifyDockWidget(&ct_config, &plan_alg_config);
  tabifyDockWidget(&plan_alg_config, &ct_interact);
  tabifyDockWidget(&ct_interact, &target_pred_config);
  
  
  connect(actionStart_Sol_Anim, SIGNAL(triggered()), this, SLOT(startSolutionAnimation()));
  connect(actionStop_Sol_Anim, SIGNAL(triggered()), this, SLOT(stopSolutionAnimation()));
  connect(actionStart_Target_Anim, SIGNAL(triggered()), this, SLOT(startTargetAnimation()));
  connect(actionStop_Target_Anim, SIGNAL(triggered()), this, SLOT(stopTargetAnimation()));
  connect(&ct_interact, SIGNAL(onLoadTargetTrajectory(QString)), this, SLOT(loadTargetTrajectory(QString)));
  
  connect(&ct_config, SIGNAL(onChaserLoaded()), &ct_interact, SLOT(loadJointPosFromModel()));
  connect(&ct_config, SIGNAL(onTargetLoaded()), &ct_interact, SLOT(loadTargetPosFromModel()));
  
  connect(actionRun_Planner, SIGNAL(triggered()), this, SLOT(runPlanner()));
  connect(actionRun_Dialog, SIGNAL(triggered()), this, SLOT(onShowRunDialog()));
  
  connect(&run_dialog, SIGNAL(triggeredStartPlanning(int)), this, SLOT(onStartPlanning(int)));
  connect(&run_dialog, SIGNAL(triggeredStopPlanning()), this, SLOT(onStopPlanning()));
  connect(&run_dialog, SIGNAL(triggeredLaunch(int)), this, SLOT(onLaunch(int)));
  connect(&run_dialog, SIGNAL(triggeredAbort()), this, SLOT(onAbort()));
  
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








