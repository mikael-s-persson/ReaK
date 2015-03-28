
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


#include "ui_planner_3Dview.h"
#include "ui_planner_space_config.h"
#include "ui_planner_alg_config.h"

#include <QDockWidget>

class SoSeparator;
class SoQtExaminerViewer;


class Planner3DWindow : public QMainWindow, private Ui::Planner3DView {
  Q_OBJECT

public:
  Planner3DWindow( QWidget* parent = 0, Qt::WindowFlags flags = 0 );
  ~Planner3DWindow();

private slots:

  void executePlanner();
#if 0
    void startRobot();
    void onJointChange();
    void onTargetChange();
    
    void onRobotVisible();
    void onRobotKinVisible();
    void onTargetVisible();
    void onEnvVisible();
    void onProxyVisible();
    void onMGVisible();
    void onSolutionsVisible();
#endif

  void onLoadTopology();
  void onLoadRobotModel();

private:
  void refreshTopoData();

  void onProxyChange();

  Ui::PlannerSpaceConfig space_configs;
  QDockWidget* space_configs_dock;
  QWidget* space_configs_widget;
  Ui::PlannerAlgConfig alg_configs;
  QDockWidget* alg_configs_dock;
  QWidget* alg_configs_widget;


  SoQtExaminerViewer* eviewer;
  SoSeparator* sg_root;
};


#endif
