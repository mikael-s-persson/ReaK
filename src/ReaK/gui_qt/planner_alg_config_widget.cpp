
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

#include "planner_alg_config_widget.hpp"

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QScrollArea>

#include "serialization/archiver_factory.hpp"

#include "shapes/oi_scene_graph.hpp"
#include "proximity/proxy_query_model.hpp"


namespace ReaK {
  
namespace rkqt {


static QString last_used_path;



PlannerAlgConfigWidget::PlannerAlgConfigWidget(QWidget * parent, Qt::WindowFlags flags) :
                                               QDockWidget(tr("Planner"), parent, flags),
                                               Ui::PlannerAlgConfig(),
                                               planOptions()
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  setupUi(dock_wid);
  
  connect(this->planning_algo_selection, SIGNAL(currentIndexChanged(int)), this, SLOT(onUpdateAvailableOptions(int)));
  connect(this->actionValuesChanged, SIGNAL(triggered()), this, SLOT(onConfigsChanged()));
  connect(this->load_button, SIGNAL(clicked()), this, SLOT(loadPlannerConfig()));
  connect(this->save_button, SIGNAL(clicked()), this, SLOT(savePlannerConfig()));
  
  
  planOptions.planning_algo = 0;
  planOptions.max_vertices = 2000;
  planOptions.prog_interval = 500;
  planOptions.max_results = 50;
  planOptions.planning_options = 0;
  planOptions.store_policy = 0;
  planOptions.knn_method = 2;
  planOptions.init_SA_temp = 0.0;
  planOptions.init_relax = 5.0;
  planOptions.max_random_walk = 1.0;
  planOptions.start_delay = 20.0;
  updateConfigs();
  
};

PlannerAlgConfigWidget::~PlannerAlgConfigWidget() {
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
};




void PlannerAlgConfigWidget::onConfigsChanged() {
  
  
  planOptions.max_random_walk = this->max_expansion_spinbox->value();
  
  // planner parameters:
  planOptions.planning_algo = this->planning_algo_selection->currentIndex();
  
  planOptions.max_vertices = this->maxvertices_spinbox->value();
  planOptions.prog_interval = this->progress_interval_spinbox->value();
  planOptions.max_results = this->maxsolutions_spinbox->value();
  
  planOptions.planning_options = pp::UNIDIRECTIONAL_PLANNING;
  
  if(this->check_bidir->isChecked())
    planOptions.planning_options |= pp::BIDIRECTIONAL_PLANNING;
  
  if(this->check_lazy_collision->isChecked())
    planOptions.planning_options |= pp::LAZY_COLLISION_CHECKING;
  
  planOptions.init_SA_temp = -1.0;
  if(this->check_voronoi_pull->isChecked()) {
    planOptions.planning_options |= pp::PLAN_WITH_VORONOI_PULL;
    planOptions.init_SA_temp = this->init_sa_temp_spinbox->value();
    if(planOptions.init_SA_temp < 1e-6)
      planOptions.init_SA_temp = -1.0;
  };
  
  planOptions.init_relax = 0.0;
  if(this->check_anytime_heuristic->isChecked()) {
    planOptions.planning_options |= pp::PLAN_WITH_ANYTIME_HEURISTIC;
    planOptions.init_relax = this->init_relax_spinbox->value();
  };
  
  if(this->check_bnb->isChecked())
    planOptions.planning_options |= pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;
  
  planOptions.start_delay = this->start_delay_spinbox->value();
  
  planOptions.store_policy = pp::ADJ_LIST_MOTION_GRAPH;
  if(this->graph_storage_selection->currentIndex())
    planOptions.store_policy = pp::DVP_ADJ_LIST_MOTION_GRAPH;
  
  planOptions.knn_method = pp::LINEAR_SEARCH_KNN;
  switch(this->KNN_method_selection->currentIndex()) {
    case 1:
      planOptions.knn_method = pp::DVP_BF2_TREE_KNN;
      break;
    case 2:
      planOptions.knn_method = pp::DVP_BF4_TREE_KNN;
      break;
    case 3:
      planOptions.knn_method = pp::DVP_COB2_TREE_KNN;
      break;
    case 4:
      planOptions.knn_method = pp::DVP_COB4_TREE_KNN;
      break;
    default:
      break;
  };
  
};


void PlannerAlgConfigWidget::updateConfigs() {
  
  
  this->max_expansion_spinbox->setValue(planOptions.max_random_walk);
  
  // planner parameters:
  this->planning_algo_selection->setCurrentIndex(planOptions.planning_algo);
  
  this->maxvertices_spinbox->setValue(planOptions.max_vertices);
  this->progress_interval_spinbox->setValue(planOptions.prog_interval);
  this->maxsolutions_spinbox->setValue(planOptions.max_results);
  
  if(planOptions.planning_options & pp::BIDIRECTIONAL_PLANNING)
    this->check_bidir->setChecked(true);
  else
    this->check_bidir->setChecked(false);
  
  
  if(planOptions.planning_options & pp::LAZY_COLLISION_CHECKING)
    this->check_lazy_collision->setChecked(true);
  else
    this->check_lazy_collision->setChecked(false);
  
  
  if(planOptions.planning_options & pp::PLAN_WITH_VORONOI_PULL) {
    this->init_sa_temp_spinbox->setValue(planOptions.init_SA_temp);
    this->check_voronoi_pull->setChecked(true);
  } else {
    this->init_sa_temp_spinbox->setValue(0.0);
    this->check_voronoi_pull->setChecked(false);
  };
  
  if(planOptions.planning_options & pp::PLAN_WITH_ANYTIME_HEURISTIC) {
    this->init_relax_spinbox->setValue(planOptions.init_relax);
    this->check_anytime_heuristic->setChecked(true);
  } else {
    this->init_relax_spinbox->setValue(0.0);
    this->check_anytime_heuristic->setChecked(false);
  };
  
  if(planOptions.planning_options & pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG)
    this->check_bnb->setChecked(true);
  else
    this->check_bnb->setChecked(false);
  
  this->start_delay_spinbox->setValue(planOptions.start_delay);
  
  
  if( planOptions.store_policy == pp::DVP_ADJ_LIST_MOTION_GRAPH ) {
    this->graph_storage_selection->setCurrentIndex(1);
  } else {
    this->graph_storage_selection->setCurrentIndex(0);
  };
  
  switch(planOptions.knn_method) {
    case pp::DVP_BF2_TREE_KNN:
      this->KNN_method_selection->setCurrentIndex(1);
      break;
    case pp::DVP_BF4_TREE_KNN:
      this->KNN_method_selection->setCurrentIndex(2);
      break;
    case pp::DVP_COB2_TREE_KNN:
      this->KNN_method_selection->setCurrentIndex(3);
      break;
    case pp::DVP_COB4_TREE_KNN:
      this->KNN_method_selection->setCurrentIndex(4);
      break;
    default:
      this->KNN_method_selection->setCurrentIndex(0);
      break;
  };
  
  onUpdateAvailableOptions(this->planning_algo_selection->currentIndex());
  
};



void PlannerAlgConfigWidget::onUpdateAvailableOptions(int plan_alg) {
  
  switch(plan_alg) {
    case 0:  // RRT
    case 1:  // RRT*
      this->max_expansion_spinbox->setEnabled(false);
      break;
    case 2:  // PRM
    case 3:  // SBA*
    case 4:  // FADPRM
    default:
      this->max_expansion_spinbox->setEnabled(true);
      break;
  };
  
  switch(plan_alg) {
    case 1:  // RRT*
      this->check_lazy_collision->setEnabled(false);
      this->check_lazy_collision->setChecked(true);
      break;
    case 3:  // SBA*
      this->check_lazy_collision->setEnabled(true);
      this->check_lazy_collision->setChecked(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      this->check_lazy_collision->setEnabled(false);
      this->check_lazy_collision->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 0:  // RRT
      this->check_bidir->setEnabled(true);
      break;
    case 1:  // RRT*
    case 2:  // PRM
    case 3:  // SBA*
    case 4:  // FADPRM
    default:
      this->check_bidir->setEnabled(false);
      this->check_bidir->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 3:  // SBA*
      this->check_voronoi_pull->setEnabled(true);
      this->check_anytime_heuristic->setEnabled(true);
      break;
    case 0:  // RRT
    case 1:  // RRT*
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      this->check_voronoi_pull->setEnabled(false);
      this->check_voronoi_pull->setChecked(false);
      this->check_anytime_heuristic->setEnabled(false);
      this->check_anytime_heuristic->setChecked(false);
      break;
  };
  
  switch(plan_alg) {
    case 1:  // RRT*
    case 3:  // SBA*
      this->check_bnb->setEnabled(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      this->check_bnb->setEnabled(false);
      this->check_bnb->setChecked(false);
      break;
  };
  
};


void PlannerAlgConfigWidget::savePlannerConfig() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Planner Configurations..."), last_used_path,
    tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb *.raplan.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    (*serialization::open_oarchive(fileName.toStdString())) << planOptions;
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void PlannerAlgConfigWidget::loadPlannerConfig() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Planner Configurations..."), last_used_path,
    tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb *.raplan.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    (*serialization::open_iarchive(fileName.toStdString())) >> planOptions;
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
};



bool PlannerAlgConfigWidget::outputTiming() const {
  return check_print_timing->isChecked();
};

bool PlannerAlgConfigWidget::outputNodeCounter() const {
  return check_print_counter->isChecked();
};

bool PlannerAlgConfigWidget::outputMotionGraph() const {
  return check_print_graph->isChecked();
};

bool PlannerAlgConfigWidget::outputSolution() const {
  return check_print_best->isChecked();
};




};

};














