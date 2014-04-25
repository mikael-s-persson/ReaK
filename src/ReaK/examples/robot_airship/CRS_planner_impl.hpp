
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
#include "path_planning/path_planner_options.hpp"
#include "kte_models/chaser_target_model_data.hpp"

#include "rk_view3d_menu.hpp"

#include "chaser_target_config_widget.hpp"
#include "chaser_target_interact_widget.hpp"
#include "planner_alg_config_widget.hpp"
#include "manip_space_config_widget.hpp"
#include "target_pred_config_widget.hpp"

#include "CRS_run_dialog.hpp"


class SoSensor;



class CRSPlannerGUI : public QMainWindow, private Ui::CRSPlannerWindow {
    Q_OBJECT
  
  public:
    CRSPlannerGUI( QWidget * parent = 0, Qt::WindowFlags flags = 0 );
    ~CRSPlannerGUI();
    
  private slots:
    
    void runPlanner();
    
    void executePlanner();
    void executeDynamicPlanner();
    
    void startSolutionAnimation();
    void startTargetAnimation();
    
    void stopSolutionAnimation();
    void stopTargetAnimation();
    
    void loadTargetTrajectory(QString);
    
    
    void onShowRunDialog();
    
    void onStartPlanning(int mode);
    void onStopPlanning();
    
    void onLaunch(int mode);
    
    void onAbort();
    
    
    friend void CRSPlannerGUI_animate_bestsol_trajectory(void*, SoSensor*);
    friend void CRSPlannerGUI_animate_target_trajectory(void*, SoSensor*);
    
    
  private:
    
    void startCompleteAnimation();
    void stopCompleteAnimation();
    
    // non-copyable:
    CRSPlannerGUI( const CRSPlannerGUI& );
    CRSPlannerGUI& operator=( const CRSPlannerGUI& );
    
    CRS_sol_anim_data    sol_anim;
    CRS_target_anim_data target_anim;
    double current_target_anim_time;
    
    ReaK::rkqt::View3DMenu view3d_menu;
    
    ReaK::rkqt::ChaserTargetConfigWidget   ct_config;
    ReaK::rkqt::ChaserTargetInteractWidget ct_interact;
    ReaK::rkqt::ManipSpaceConfigWidget     space_config;
    ReaK::rkqt::PlannerAlgConfigWidget     plan_alg_config;
    ReaK::rkqt::TargetPredConfigWidget     target_pred_config;
    
    ReaK::rkqt::CRSRunDialogWidget run_dialog;
    
    void threadedPlanningFunction(int mode);
    
    ReaKaux::thread* planning_thr;
    
    void executeSolutionTrajectory();
    
    volatile bool exec_robot_enabled;
    ReaKaux::thread* exec_robot_thr;
    
    
};



#endif














