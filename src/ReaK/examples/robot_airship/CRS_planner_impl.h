
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




struct CRS_coin_nodes;
struct CRS_model_data;
struct CRS_sol_anim_data;
struct CRS_target_anim_data;
struct CRS_planning_options;


class SoSensor;



class CRSPlannerGUI : public QMainWindow, private Ui::CRSPlannerWindow {
    Q_OBJECT
  
  public:
    CRSPlannerGUI( QWidget * parent = 0, Qt::WindowFlags flags = 0 );
    ~CRSPlannerGUI();
    
  private slots:
    
    void executePlanner();
    void startSolutionAnimation();
    void startTargetAnimation();
    void onJointChange();
    void onTargetChange();
    
    void onRobotVisible();
    void onRobotKinVisible();
    void onTargetVisible();
    void onEnvVisible();
    void onProxyVisible();
    void onMGVisible();
    void onSolutionsVisible();
    
    void onUpdateAvailableOptions();
    
    void onConfigsChanged();
    void updateConfigs();
    
    void savePositions();
    void loadPositions();
    
    void savePlannerConfig();
    void loadPlannerConfig();
    
    void loadTargetTrajectory();
    
    
    friend void CRSPlannerGUI_animate_bestsol_trajectory(void*, SoSensor*);
    friend void CRSPlannerGUI_animate_target_trajectory(void*, SoSensor*);
    
  private:
    
    CRSPlannerGUI( const CRSPlannerGUI& );
    CRSPlannerGUI& operator=( const CRSPlannerGUI& );
    
    void onProxyChange();
    
    Ui::CRSPlannerConfig configs;
    
    CRS_coin_nodes* draw_data;
    CRS_model_data* mdl_data;
    CRS_sol_anim_data* sol_anim;
    CRS_target_anim_data* target_anim;
    CRS_planning_options* plan_options;
    
};



#endif














