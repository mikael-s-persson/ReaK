/**
 *
 *
 *
 */

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

#ifndef REAK_PLANNER_ALG_CONFIG_WIDGET_HPP
#define REAK_PLANNER_ALG_CONFIG_WIDGET_HPP

#include <ReaK/planning/path_planning/path_planner_options.hpp>

#include <QDockWidget>

namespace Ui {
class PlannerAlgConfig;
};

namespace ReaK {

namespace qt {

class PlannerAlgConfigWidget : public QDockWidget {
  Q_OBJECT
private:
  Ui::PlannerAlgConfig* ui;

public:
  PlannerAlgConfigWidget( QWidget* parent = NULL, Qt::WindowFlags flags = 0 );
  virtual ~PlannerAlgConfigWidget();

private slots:

  void onUpdateAvailableOptions( int );

  void onConfigsChanged();
  void updateConfigs();

  void savePlannerConfig();
  void loadPlannerConfig();

public:
  pp::planning_option_collection planOptions;

  bool outputTiming() const;
  bool outputNodeCounter() const;
  bool outputMotionGraph() const;
  bool outputSolution() const;

  void savePlannerConfiguration( const std::string& aFilename );
  void loadPlannerConfiguration( const std::string& aFilename );
};
};
};

#endif
