
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

#ifndef REAK_CRS_PLANNER_IMPL_H
#define REAK_CRS_PLANNER_IMPL_H


#include "ui_CRS_planner_window.h"
#include "ui_CRS_planner_config.h"

#include "CRS_planner_data.hpp"
#include <ReaK/planning/path_planning/path_planner_options.hpp>
#include <ReaK/mbd/models/chaser_target_model_data.hpp>

#include <ReaK/mbd/qt/rk_view3d_menu.hpp>
#include <ReaK/core/base/thread_incl.hpp>
#include <ReaK/core/base/atomic_incl.hpp>
#include <ReaK/core/recorders/data_record_options.hpp>

#include <ReaK/mbd/qt/chaser_target_config_widget.hpp>
#include <ReaK/mbd/qt/chaser_target_interact_widget.hpp>
#include <ReaK/planning/qt/planner_alg_config_widget.hpp>
#include <ReaK/planning/qt/manip_space_config_widget.hpp>
#include <ReaK/control/qt/target_pred_config_widget.hpp>

#include "CRS_run_dialog.hpp"


class SoSensor;


class CRSPlannerGUI : public QMainWindow, private Ui::CRSPlannerWindow {
  Q_OBJECT

public:
  CRSPlannerGUI( QWidget* parent = 0, Qt::WindowFlags flags = 0 );
  ~CRSPlannerGUI();

private slots:

  void runPlanner();

  void executePlanner();
  void executeDynamicPlanner();

  void startSolutionAnimation();
  void startTargetAnimation();

  void stopSolutionAnimation();
  void stopTargetAnimation();

  void loadTargetTrajectory( QString );


  void onShowRunDialog();

  void onStartPlanning( int mode );
  void onStopPlanning();

  void onLaunch( int mode );

  void onAbort();


  friend void CRSPlannerGUI_animate_bestsol_trajectory( void*, SoSensor* );
  friend void CRSPlannerGUI_animate_target_trajectory( void*, SoSensor* );

signals:

  void notifyCaptureReached();

  void notifyConsoleMessage( QString );

  void notifyReset();

  void notifyInitializationDone();

  void notifyPlanningDone();

  void notifyLaunchOpportunity();

private:
  void startCompleteAnimation();
  void stopCompleteAnimation();

  // non-copyable:
  CRSPlannerGUI( const CRSPlannerGUI& );
  CRSPlannerGUI& operator=( const CRSPlannerGUI& );

public:
  CRS_sol_anim_data sol_anim;
  CRS_target_anim_data target_anim;
  ReaKaux::atomic< double > current_target_anim_time;

  ReaK::qt::View3DMenu view3d_menu;

  ReaK::qt::ChaserTargetConfigWidget ct_config;
  ReaK::qt::ChaserTargetInteractWidget ct_interact;
  ReaK::qt::ManipSpaceConfigWidget space_config;
  ReaK::qt::PlannerAlgConfigWidget plan_alg_config;
  ReaK::qt::TargetPredConfigWidget target_pred_config;

  ReaK::qt::CRSRunDialogWidget run_dialog;

  ReaKaux::function< void() > stop_planner;
  void threadedPlanningFunction( int mode );

  ReaKaux::thread planning_thr;

  void executeSolutionTrajectory();

  ReaKaux::atomic< bool > exec_robot_enabled;
  ReaKaux::thread exec_robot_thr;
  ReaK::recorder::data_stream_options jtctrl_log_opt;
  ReaK::recorder::data_stream_options jtctrl_network_opt;
};


#endif
