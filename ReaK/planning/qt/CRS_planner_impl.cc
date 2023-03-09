
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
#include <QDir>
#include <QFileDialog>
#include <QMainWindow>
#include <QMessageBox>

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/SoDB.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/mbd/coin3D/oi_scene_graph.h"

#include "ReaK/core/serialization/archiver_factory.h"

#include "ReaK/topologies/interpolation/discrete_point_trajectory.h"
#include "ReaK/topologies/interpolation/trajectory_base.h"

#include "CRS_planner_data.hpp"

#include "ReaK/math/optimization/optim_exceptions.h"

#include <chrono>

#include <functional>
#include "boost/asio.hpp"

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

// Configuration options
ABSL_FLAG(std::string, scene_data, "", "File-name of a complete model file.");
ABSL_FLAG(std::string, chaser_target_pos, "",
          "File-name of a chaser-target positions file.");
ABSL_FLAG(std::string, planning_cfg, "",
          "File-name of a planning algorithm configuration file.");
ABSL_FLAG(std::string, space_cfg, "",
          "File-name of a robot state-space configuration file.");
ABSL_FLAG(std::string, predictor_cfg, "",
          "File-name of a target predictor configuration file.");
ABSL_FLAG(std::string, jtctrl_log_cfg, "",
          "File-name of a configuration file for the chaser joint logging.");
ABSL_FLAG(
    std::string, jtctrl_network_cfg, "",
    "File-name of a configuration file for the chaser joint control network.");
ABSL_FLAG(
    std::string, target_meas_log_cfg, "",
    "File-name of a configuration file for the target measurement logging.");
ABSL_FLAG(
    std::string, target_est_log_cfg, "",
    "File-name of a configuration file for the target estimates logging.");
ABSL_FLAG(
    std::string, target_pred_log_cfg, "",
    "File-name of a configuration file for the target predictions logging.");

using namespace ReaK;

static QString last_used_path;

// #define RK_CRSPLANNER_USE_RAW_ASIO_SOCKET

namespace {

union double_to_ulong {
  double d;
  uint64_t ui64;
  uint32_t ui32[2];
};

template <typename UnionT>
void ntoh_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = ntohl(value.ui32[0]);
  value.ui32[0] = ntohl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP
  // platforms.
};

template <typename UnionT>
void hton_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = htonl(value.ui32[0]);
  value.ui32[0] = htonl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP
  // platforms.
};

#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
void write_double_to_net_stream(std::ostream& out, double value) {
  double_to_ulong tmp;
  tmp.d = value;
  hton_2ui32(tmp);
  out.write(reinterpret_cast<char*>(&tmp), sizeof(double));
  //   out.write(reinterpret_cast<char*>(&value),sizeof(double));
};
#endif
};  // namespace

void CRSPlannerGUI_animate_bestsol_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);

  static std::shared_ptr<CRS_sol_anim_data::trajectory_type> manip_traj;
  static CRS_sol_anim_data::trajectory_type::point_time_iterator cur_pit;
  static std::chrono::high_resolution_clock::time_point animation_start;

  if (!manip_traj) {
    manip_traj = p->sol_anim.trajectory;
    cur_pit = manip_traj->begin_time_travel();

    // Synchronize the trajectory to the current time (from target trajectory).
    //  at the start of planning, starting point is at time:  (delay + target-start-time)
    //  therefore, start of solution, at time start-time, should satisfy:  target-start-time  ==  start-time - delay
    //  and then, because we don't start animation at the start-time, we can adjust to the current time with:
    //    t-offset = target-curr-time - (start-time - delay)
    //     double t_offset = p->current_target_anim_time - ((*cur_pit).time -
    //     p->plan_alg_config.planOptions.start_delay);
    if (p->current_target_anim_time > (*cur_pit).time)
      cur_pit += p->current_target_anim_time - (*cur_pit).time;

    animation_start = std::chrono::high_resolution_clock::now() -
                      std::chrono::duration_cast<
                          std::chrono::high_resolution_clock::duration>(
                          std::chrono::duration<double>(
                              p->current_target_anim_time -
                              p->plan_alg_config.planOptions.start_delay));
  };
  if ((p->sol_anim.enabled) && ((*cur_pit).time < manip_traj->get_end_time())) {
    if ((*cur_pit).time <=
        0.001 *
            (std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::high_resolution_clock::now() - animation_start))
                .count()) {
      //     if( ((*cur_pit).time - manip_traj->get_start_time()) <= 0.001 *
      //     (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()
      //     - animation_start)).count() ) {
      cur_pit += 0.1;
      p->ct_config.sceneData.chaser_kin_model->setJointPositions(
          vect_n<double>((*cur_pit).pt));
      p->ct_config.sceneData.chaser_kin_model->doDirectMotion();
    };
  } else {
    p->sol_anim.animation_timer->unschedule();
    animation_start = std::chrono::high_resolution_clock::now();
    manip_traj.reset();
    emit p->notifyCaptureReached();
  };
};

void CRSPlannerGUI::startSolutionAnimation() {

  if (space_config.is_temporal) {
    startCompleteAnimation();
    return;
  };

  if (!sol_anim.trajectory || !ct_config.sceneData.chaser_kin_model) {
    QMessageBox::critical(this, "Animation Error!",
                          "The best-solution trajectory is missing (not loaded "
                          "or erroneous)! Cannot animate chaser!",
                          QMessageBox::Ok);
    return;
  };
  sol_anim.enabled = true;
  sol_anim.animation_timer->schedule();
};

void CRSPlannerGUI::stopSolutionAnimation() {

  if (space_config.is_temporal) {
    stopCompleteAnimation();
    return;
  };

  sol_anim.enabled = false;
};

void CRSPlannerGUI_animate_target_trajectory(void* pv, SoSensor*) {
  CRSPlannerGUI* p = static_cast<CRSPlannerGUI*>(pv);

  static std::shared_ptr<CRS_target_anim_data::trajectory_type> target_traj;
  static CRS_target_anim_data::trajectory_type::point_time_iterator cur_pit;
  static std::chrono::high_resolution_clock::time_point animation_start;

  if (!target_traj) {
    target_traj = p->target_anim.trajectory;
    cur_pit = target_traj->begin_time_travel();

    // Synchronize the trajectory to the current time (from target trajectory).
    double t_offset = p->current_target_anim_time - (*cur_pit).time;
    if (t_offset > 0.0)
      cur_pit += t_offset;

    animation_start =
        std::chrono::high_resolution_clock::now() -
        std::chrono::duration_cast<
            std::chrono::high_resolution_clock::duration>(
            std::chrono::duration<double>(p->current_target_anim_time));
  };

  if (target_traj && p->target_anim.enabled &&
      ((*cur_pit).time < target_traj->get_end_time())) {
    if ((*cur_pit).time <=
        0.001 *
            (std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::high_resolution_clock::now() - animation_start))
                .count()) {
      //     if( ((*cur_pit).time - target_traj->get_start_time()) <= 0.001 *
      //     (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()
      //     - animation_start)).count() ) {
      cur_pit += 0.1;
      *(p->ct_config.sceneData.target_kin_model->getFrame3D(0)) =
          get_frame_3D((*cur_pit).pt);
      //       std::cout << " Position = " << p->ct_config.sceneData.target_kin_model->getFrame3D(0)->Position
      //                 << " Rotation = " << p->ct_config.sceneData.target_kin_model->getFrame3D(0)->Quat << std::endl;
      p->ct_config.sceneData.target_kin_model->doDirectMotion();

      if (p->ct_config.sceneData.chaser_kin_model &&
          p->ct_interact.isIKEnabled() && !(p->sol_anim.enabled)) {
        try {
          *(p->ct_config.sceneData.chaser_kin_model->getDependentFrame3D(0)
                ->mFrame) = *(p->ct_config.sceneData.target_frame);
          p->ct_config.sceneData.chaser_kin_model->doInverseMotion();
        } catch (optim::infeasible_problem& e) {
          RK_UNUSED(e);
        };
        p->ct_config.sceneData.chaser_kin_model->doDirectMotion();
      };
    };
  } else {
    p->target_anim.animation_timer->unschedule();
    target_traj.reset();
  };
};

void CRSPlannerGUI::startTargetAnimation() {
  if (!target_anim.trajectory || !ct_config.sceneData.target_kin_model) {
    QMessageBox::critical(this, "Animation Error!",
                          "The target trajectory is missing (not loaded or "
                          "erroneous)! Cannot animate target!",
                          QMessageBox::Ok);
    return;
  };
  target_anim.enabled = true;
  target_anim.animation_timer->schedule();
  current_target_anim_time = target_anim.trajectory->get_start_time();
};

void CRSPlannerGUI::stopTargetAnimation() {
  target_anim.enabled = false;
};

void CRSPlannerGUI::loadTargetTrajectory(QString fileName) {
  typedef pp::se3_1st_order_topology<double>::type sat_state_space_type;
  typedef pp::temporal_space<sat_state_space_type, pp::time_poisson_topology,
                             pp::time_distance_only>
      temporal_space_type;
  typedef pp::discrete_point_trajectory<temporal_space_type> trajectory_type;

  typedef pp::trajectory_wrapper<trajectory_type> wrapped_traj_type;

  std::shared_ptr<wrapped_traj_type> tmp_traj(new wrapped_traj_type());

  *(serialization::open_iarchive(fileName.toStdString())) &
      RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory",
                                tmp_traj->get_underlying_trajectory());
  target_anim.trajectory = tmp_traj;

  //   *(serialization::open_iarchive(fileName.toStdString()))
  //     & RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory", target_anim.trajectory);
};

void CRSPlannerGUI::startCompleteAnimation() {

  if (!sol_anim.trajectory || !target_anim.trajectory ||
      !ct_config.sceneData.chaser_kin_model ||
      !ct_config.sceneData.target_kin_model) {
    QMessageBox::critical(this, "Animation Error!",
                          "One of the trajectories is missing (not loaded or "
                          "erroneous)! Cannot animate chaser and target!",
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
    if (space_config.is_temporal) {
      executeDynamicPlanner();
    } else {
      executePlanner();
    };
  } catch (std::exception& e) {
    std::stringstream ss;
    ss << "An exception was raised during the planning:\nwhat(): " << e.what();
    QMessageBox::critical(
        this, "Motion-Planning Error!",
        "An exception was raised during the planning:\nwhat(): " +
            QString(e.what()),
        QMessageBox::Ok);
  };
};

void CRSPlannerGUI::onShowRunDialog() {
  run_dialog.setModal(false);
  run_dialog.show();
};

void CRSPlannerGUI::threadedPlanningFunction(int mode) {

  if ((mode == 2) || (mode == 3)) {

    emit notifyConsoleMessage("Starting state estimation...\n");

    target_anim.trajectory.reset();
    target_pred_config.startStatePrediction();
    if (!target_anim.trajectory) {
      emit notifyConsoleMessage(
          "State-Prediction Error!\nThe live state-estimation of the target "
          "failed to produce a "
          "viable predicted trajectory!");
      emit notifyReset();
      return;
    };

    emit notifyConsoleMessage(
        "State estimation done!\nSwitched to predicted target trajectory!\n");
  };

  emit notifyInitializationDone();

  try {
    if (mode > 0) {
      executeDynamicPlanner();
    } else {
      executePlanner();
    };
  } catch (std::exception& e) {
    emit notifyConsoleMessage("Motion-Planning Error!\nwhat(): " +
                              QString(e.what()));
    emit notifyReset();
  };

  emit notifyPlanningDone();

  if (sol_anim.trajectory) {
    emit notifyLaunchOpportunity();
  };
};

void CRSPlannerGUI::onStartPlanning(int mode) {
  if ((mode > 0) && !space_config.is_temporal) {
    QMessageBox::critical(
        this, "Motion-Planning Error!",
        "Trying to run a dynamic planner with a non-temporal planning space! "
        "Please configure a temporal planning space!",
        QMessageBox::Ok);
    emit notifyReset();
    return;
  };

  //   mode == 0;  // static-sim (executePlanner)
  //   mode == 1;  // dynamic-sim (loaded target-trajectory + executeDynamicPlanner + startCompleteAnimation)
  //   mode == 2;  // live-sim (predicted target-trajectory + executeDynamicPlanner + startCompleteAnimation)
  //   mode == 3;  // live-run (predicted target-trajectory + executeDynamicPlanner + startCompleteAnimation +
  //   executeSolTrajectory)

  if (stop_planner)
    stop_planner();
  if (planning_thr.joinable())
    planning_thr.join();

  planning_thr = std::thread(
      std::bind(&CRSPlannerGUI::threadedPlanningFunction, this, mode));
};

void CRSPlannerGUI::onStopPlanning() {

  target_pred_config.stopStatePrediction();

  if (stop_planner)
    stop_planner();
  if (planning_thr.joinable())
    planning_thr.join();
};

void CRSPlannerGUI::executeSolutionTrajectory() {

  if (!sol_anim.trajectory) {
    emit notifyConsoleMessage(
        "Execution Error! Cannot execute planned path, because no such plan "
        "exists!");
    return;
  };

#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
  // Setup the UDP streaming to the robot (FIXME remove the hard-coded address / port)
  boost::asio::io_service io_service;
  //   boost::asio::ip::udp::endpoint endpoint(boost::asio::ip::address_v4::from_string("127.0.0.1"), 17050);
  boost::asio::ip::udp::endpoint endpoint(
      boost::asio::ip::address_v4::from_string("192.168.0.3"), 17050);
  boost::asio::ip::udp::socket socket((io_service));
  boost::asio::basic_streambuf<> udp_buf;
  std::ostream udp_buf_stream(&(udp_buf));
  socket.open(boost::asio::ip::udp::v4());
  socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
#else
  std::shared_ptr<recorder::data_recorder> jtctrl_net =
      jtctrl_network_opt.create_recorder();
#endif

  // Setup the recording of the control signals to a space-separated file:
  std::shared_ptr<recorder::data_recorder> jtctrl_out =
      jtctrl_log_opt.create_recorder();

  std::shared_ptr<CRS_sol_anim_data::trajectory_type> manip_traj =
      sol_anim.trajectory;
  CRS_sol_anim_data::trajectory_type::point_time_iterator cur_pit =
      manip_traj->begin_time_travel();

  std::chrono::high_resolution_clock::time_point exec_start =
      std::chrono::high_resolution_clock::now();

  // UDP output for the robot:
  vect<double, 7> cur_pt = (*cur_pit).pt;
#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
  for (std::size_t i = 0; i < 7; ++i)
    write_double_to_net_stream(udp_buf_stream, cur_pt[i]);
  std::size_t len = socket.send_to(udp_buf.data(), endpoint);
  udp_buf.consume(len);
#else
  (*jtctrl_net) << cur_pt << recorder::data_recorder::end_value_row;
#endif

  double cur_time = (*cur_pit).time;

  // Synchronize the trajectory to the current time (from target trajectory).
  if (current_target_anim_time > cur_time) {
    cur_pit += current_target_anim_time - cur_time;
  } else {

    std::chrono::high_resolution_clock::time_point delayed_start =
        exec_start +
        std::chrono::duration_cast<
            std::chrono::high_resolution_clock::duration>(
            std::chrono::duration<double>(cur_time - current_target_anim_time));

    while (exec_robot_enabled && (exec_start < delayed_start) &&
           (cur_time > current_target_anim_time)) {

      std::this_thread::sleep_until(exec_start +
                                    std::chrono::microseconds(1000));

      exec_start = std::chrono::high_resolution_clock::now();
      if (jtctrl_out)
        (*jtctrl_out) << 0.0;
#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
      for (std::size_t i = 0; i < 7; ++i) {
        write_double_to_net_stream(udp_buf_stream, cur_pt[i]);
        if (jtctrl_out)
          (*jtctrl_out) << cur_pt[i];
      };
      std::size_t len = socket.send_to(udp_buf.data(), endpoint);
      udp_buf.consume(len);
      if (jtctrl_out)
        (*jtctrl_out) << recorder::data_recorder::end_value_row;
#else
      (*jtctrl_net) << cur_pt << recorder::data_recorder::end_value_row;
      if (jtctrl_out)
        (*jtctrl_out) << cur_pt << recorder::data_recorder::end_value_row;
#endif
    };
  };

  exec_start = std::chrono::high_resolution_clock::now();

  while (exec_robot_enabled && ((*cur_pit).time < manip_traj->get_end_time())) {
    cur_pit += 0.001;
    std::this_thread::sleep_until(exec_start + std::chrono::microseconds(1000));
    //     cur_pit += 0.01;
    //     std::this_thread::sleep_until(exec_start + std::chrono::microseconds(10000));

    exec_start = std::chrono::high_resolution_clock::now();

    // UDP output for the robot:
    if (jtctrl_out)
      (*jtctrl_out) << (*cur_pit).time;
    cur_pt = (*cur_pit).pt;

#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
    for (std::size_t i = 0; i < 7; ++i) {
      write_double_to_net_stream(udp_buf_stream, cur_pt[i]);
      if (jtctrl_out)
        (*jtctrl_out) << cur_pt[i];
    };
    len = socket.send_to(udp_buf.data(), endpoint);
    udp_buf.consume(len);
    if (jtctrl_out)
      (*jtctrl_out) << recorder::data_recorder::end_value_row;
#else
    (*jtctrl_net) << cur_pt << recorder::data_recorder::end_value_row;
    if (jtctrl_out)
      (*jtctrl_out) << cur_pt << recorder::data_recorder::end_value_row;
#endif
  };

  // UDP output for the robot:
  double terminating_value = -100.0;
#ifdef RK_CRSPLANNER_USE_RAW_ASIO_SOCKET
  for (std::size_t i = 0; i < 7; ++i)
    write_double_to_net_stream(udp_buf_stream, terminating_value);
  len = socket.send_to(udp_buf.data(), endpoint);
  udp_buf.consume(len);
#else
  (*jtctrl_net) << terminating_value << terminating_value << terminating_value
                << terminating_value << terminating_value << terminating_value
                << terminating_value << recorder::data_recorder::end_value_row;
#endif

  if (jtctrl_out)
    (*jtctrl_out) << recorder::data_recorder::flush;

  if ((*cur_pit).time >= manip_traj->get_end_time()) {
    emit notifyCaptureReached();
  };
  exec_robot_enabled = false;
};

void CRSPlannerGUI::onLaunch(int mode) {

  // first, stop the robot
  if (this->exec_robot_thr.joinable()) {
    this->exec_robot_enabled = false;
    this->exec_robot_thr.join();
  };

  if (!sol_anim.trajectory || !ct_config.sceneData.chaser_kin_model) {
    QMessageBox::critical(this, "Launch Error!",
                          "The best-solution trajectory is missing (not loaded "
                          "or erroneous)! Cannot launch chaser!",
                          QMessageBox::Ok);
    return;
  };

  if (mode == 0) {
    current_target_anim_time = 0.0;
    startSolutionAnimation();
    return;
  } else if (mode == 1) {
    current_target_anim_time = target_anim.trajectory->get_start_time();
  } else if (mode == 3) {
    exec_robot_enabled = true;
    exec_robot_thr =
        std::thread(std::bind(&CRSPlannerGUI::executeSolutionTrajectory, this));
  };
  startCompleteAnimation();
};

void CRSPlannerGUI::onAbort() {

  // first, stop the robot
  if (this->exec_robot_thr.joinable()) {
    this->exec_robot_enabled = false;
    this->exec_robot_thr.join();
  };

  // then, stop the state prediction and/or planning
  this->onStopPlanning();

  // finally, stop any animation that might have been started
  this->stopCompleteAnimation();

  // finally, restore robot to set position.
  this->ct_interact.onJointChange();
};

CRSPlannerGUI::CRSPlannerGUI(QWidget* parent, Qt::WindowFlags flags)
    : QMainWindow(parent, flags),
      view3d_menu(this),
      ct_config(&view3d_menu, this),
      ct_interact(&ct_config.sceneData, this),
      space_config(this),
      plan_alg_config(this),
      target_pred_config(&target_anim.trajectory, &current_target_anim_time,
                         this),
      run_dialog(this),
      planning_thr(),
      exec_robot_thr() {
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

  connect(actionStart_Sol_Anim, SIGNAL(triggered()), this,
          SLOT(startSolutionAnimation()));
  connect(actionStop_Sol_Anim, SIGNAL(triggered()), this,
          SLOT(stopSolutionAnimation()));
  connect(actionStart_Target_Anim, SIGNAL(triggered()), this,
          SLOT(startTargetAnimation()));
  connect(actionStop_Target_Anim, SIGNAL(triggered()), this,
          SLOT(stopTargetAnimation()));
  connect(&ct_interact, SIGNAL(onLoadTargetTrajectory(QString)), this,
          SLOT(loadTargetTrajectory(QString)));

  connect(&ct_config, SIGNAL(onChaserLoaded()), &ct_interact,
          SLOT(loadJointPosFromModel()));
  connect(&ct_config, SIGNAL(onTargetLoaded()), &ct_interact,
          SLOT(loadTargetPosFromModel()));

  connect(actionRun_Planner, SIGNAL(triggered()), this, SLOT(runPlanner()));
  connect(actionRun_Dialog, SIGNAL(triggered()), this, SLOT(onShowRunDialog()));

  connect(&run_dialog, SIGNAL(triggeredStartPlanning(int)), this,
          SLOT(onStartPlanning(int)), Qt::QueuedConnection);
  connect(&run_dialog, SIGNAL(triggeredStopPlanning()), this,
          SLOT(onStopPlanning()), Qt::QueuedConnection);
  connect(&run_dialog, SIGNAL(triggeredLaunch(int)), this, SLOT(onLaunch(int)),
          Qt::QueuedConnection);
  connect(&run_dialog, SIGNAL(triggeredAbort()), this, SLOT(onAbort()),
          Qt::QueuedConnection);

  connect(this, SIGNAL(notifyCaptureReached()), &run_dialog,
          SLOT(onCaptureReached()), Qt::QueuedConnection);
  connect(this, SIGNAL(notifyConsoleMessage(QString)), &run_dialog,
          SLOT(publishConsoleMessage(QString)), Qt::QueuedConnection);
  connect(this, SIGNAL(notifyReset()), &run_dialog, SLOT(onReset()),
          Qt::QueuedConnection);
  connect(this, SIGNAL(notifyInitializationDone()), &run_dialog,
          SLOT(onInitializationDone()), Qt::QueuedConnection);
  connect(this, SIGNAL(notifyPlanningDone()), &run_dialog,
          SLOT(onPlanningDone()), Qt::QueuedConnection);
  connect(this, SIGNAL(notifyLaunchOpportunity()), &run_dialog,
          SLOT(onLaunchOpportunity()), Qt::QueuedConnection);

  SoQt::init(this->centralwidget);

  menubar->addMenu(&view3d_menu);
  view3d_menu.setViewer(new SoQtExaminerViewer(this->centralwidget));

  view3d_menu.getGeometryGroup("Chaser Geometry", true);
  view3d_menu.getGeometryGroup("Chaser KTE Chain", false);
  view3d_menu.getGeometryGroup("Target Geometry", true);
  view3d_menu.getGeometryGroup("Environment", true);
  view3d_menu.getDisplayGroup("Motion-Graph", true);
  view3d_menu.getDisplayGroup("Solution(s)", true);

  sol_anim.animation_timer =
      new SoTimerSensor(CRSPlannerGUI_animate_bestsol_trajectory, this);
  target_anim.animation_timer =
      new SoTimerSensor(CRSPlannerGUI_animate_target_trajectory, this);

  jtctrl_log_opt.kind = recorder::data_stream_options::space_separated;
  jtctrl_log_opt.file_name =
      std::string("exp_results/robot_airship/jtctrl_$d_$t.ssv");

  jtctrl_log_opt.add_name("time")
      .add_name("q_0")
      .add_name("q_1")
      .add_name("q_2")
      .add_name("q_3")
      .add_name("q_4")
      .add_name("q_5")
      .add_name("q_6");

  jtctrl_network_opt.kind = recorder::data_stream_options::raw_udp_stream;
  jtctrl_network_opt.file_name = std::string("192.168.0.3:17050");
  //   jtctrl_network_opt.file_name = std::string("127.0.0.1:17050");
  jtctrl_network_opt.set_unbuffered();
  jtctrl_network_opt.time_sync_name = "time";

  jtctrl_network_opt.add_name("q_0")
      .add_name("q_1")
      .add_name("q_2")
      .add_name("q_3")
      .add_name("q_4")
      .add_name("q_5")
      .add_name("q_6");
};

CRSPlannerGUI::~CRSPlannerGUI() {

  // first, stop the robot
  if (this->exec_robot_thr.joinable()) {
    this->exec_robot_enabled = false;
    this->exec_robot_thr.join();
  };

  // then, stop the state prediction and/or planning
  this->onStopPlanning();

  // finally, stop any animation that might have been started
  this->stopCompleteAnimation();

  delete target_anim.animation_timer;
  delete sol_anim.animation_timer;

  view3d_menu.removeGeometryGroup("Chaser Geometry");
  view3d_menu.removeGeometryGroup("Chaser KTE Chain");
  view3d_menu.removeGeometryGroup("Target Geometry");
  view3d_menu.removeGeometryGroup("Environment");
  view3d_menu.removeDisplayGroup("Motion-Graph");
  view3d_menu.removeDisplayGroup("Solution(s)");

  view3d_menu.setViewer(nullptr);

  SoQt::done();
};

int main(int argc, char** argv) {

  absl::ParseCommandLine(argc, argv);

  QApplication app(argc, argv);

  CRSPlannerGUI window;

  window.show();

  // Pop up the main window.
  SoQt::show(&window);

  if (!absl::GetFlag(FLAGS_scene_data).empty()) {
    window.ct_config.loadCompleteModel(absl::GetFlag(FLAGS_scene_data));
  }

  if (!absl::GetFlag(FLAGS_chaser_target_pos).empty()) {
    window.ct_interact.loadChaserTargetPositions(
        absl::GetFlag(FLAGS_chaser_target_pos));
  }

  if (!absl::GetFlag(FLAGS_planning_cfg).empty()) {
    window.plan_alg_config.loadPlannerConfiguration(
        absl::GetFlag(FLAGS_planning_cfg));
  }

  if (!absl::GetFlag(FLAGS_space_cfg).empty()) {
    window.space_config.loadSpaceConfiguration(absl::GetFlag(FLAGS_space_cfg));
  }

  if (!absl::GetFlag(FLAGS_predictor_cfg).empty()) {
    window.target_pred_config.loadPredictorConfigurations(
        absl::GetFlag(FLAGS_predictor_cfg));
  }

  if (!absl::GetFlag(FLAGS_jtctrl_log_cfg).empty()) {
    window.jtctrl_log_opt.load_all_configs(absl::GetFlag(FLAGS_jtctrl_log_cfg));
  }

  if (!absl::GetFlag(FLAGS_jtctrl_network_cfg).empty()) {
    window.jtctrl_network_opt.load_all_configs(
        absl::GetFlag(FLAGS_jtctrl_network_cfg));
  }

  if (!absl::GetFlag(FLAGS_target_meas_log_cfg).empty()) {
    window.target_pred_config.meas_out_opt.load_all_configs(
        absl::GetFlag(FLAGS_target_meas_log_cfg));
  }

  if (!absl::GetFlag(FLAGS_target_est_log_cfg).empty()) {
    window.target_pred_config.est_out_opt.load_all_configs(
        absl::GetFlag(FLAGS_target_est_log_cfg));
  }

  if (!absl::GetFlag(FLAGS_target_pred_log_cfg).empty()) {
    window.target_pred_config.pred_out_opt.load_all_configs(
        absl::GetFlag(FLAGS_target_pred_log_cfg));
  }

  // Loop until exit.
  SoQt::mainLoop();

  return 0;
};
