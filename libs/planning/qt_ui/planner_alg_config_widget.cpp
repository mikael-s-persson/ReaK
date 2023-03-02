
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

#include <ReaK/planning/qt/planner_alg_config_widget.hpp>

#include "ui_planner_alg_config.h"

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QScrollArea>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/geometry/proximity/proxy_query_model.hpp>
#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>

namespace ReaK {

namespace qt {

static QString last_used_path;

PlannerAlgConfigWidget::PlannerAlgConfigWidget(QWidget* parent,
                                               Qt::WindowFlags flags)
    : QDockWidget(tr("Planner"), parent, flags),
      ui(new Ui::PlannerAlgConfig()),
      planOptions() {
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  ui->setupUi(dock_wid);

  connect(ui->planning_algo_selection, SIGNAL(currentIndexChanged(int)), this,
          SLOT(onUpdateAvailableOptions(int)));
  connect(ui->actionValuesChanged, SIGNAL(triggered()), this,
          SLOT(onConfigsChanged()));
  connect(ui->load_button, SIGNAL(clicked()), this, SLOT(loadPlannerConfig()));
  connect(ui->save_button, SIGNAL(clicked()), this, SLOT(savePlannerConfig()));

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
  delete ui;
};

void PlannerAlgConfigWidget::onConfigsChanged() {

  planOptions.max_random_walk = ui->max_expansion_spinbox->value();

  // planner parameters:
  planOptions.planning_algo = ui->planning_algo_selection->currentIndex();

  planOptions.max_vertices = ui->maxvertices_spinbox->value();
  planOptions.prog_interval = ui->progress_interval_spinbox->value();
  planOptions.max_results = ui->maxsolutions_spinbox->value();

  planOptions.planning_options = pp::UNIDIRECTIONAL_PLANNING;

  if (ui->check_bidir->isChecked())
    planOptions.planning_options |= pp::BIDIRECTIONAL_PLANNING;

  if (ui->check_lazy_collision->isChecked())
    planOptions.planning_options |= pp::LAZY_COLLISION_CHECKING;

  planOptions.init_SA_temp = -1.0;
  if (ui->check_voronoi_pull->isChecked()) {
    planOptions.planning_options |= pp::PLAN_WITH_VORONOI_PULL;
    planOptions.init_SA_temp = ui->init_sa_temp_spinbox->value();
    if (planOptions.init_SA_temp < 1e-6)
      planOptions.init_SA_temp = -1.0;
  };

  planOptions.init_relax = 0.0;
  if (ui->check_anytime_heuristic->isChecked()) {
    planOptions.planning_options |= pp::PLAN_WITH_ANYTIME_HEURISTIC;
    planOptions.init_relax = ui->init_relax_spinbox->value();
  };

  if (ui->check_bnb->isChecked())
    planOptions.planning_options |= pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;

  planOptions.start_delay = ui->start_delay_spinbox->value();

  planOptions.store_policy = pp::ADJ_LIST_MOTION_GRAPH;
  if (ui->graph_storage_selection->currentIndex())
    planOptions.store_policy = pp::DVP_ADJ_LIST_MOTION_GRAPH;

  planOptions.knn_method = pp::LINEAR_SEARCH_KNN;
  switch (ui->KNN_method_selection->currentIndex()) {
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

  ui->actionValuesChanged->disconnect(this, SLOT(onConfigsChanged()));

  ui->max_expansion_spinbox->setValue(planOptions.max_random_walk);

  // planner parameters:
  ui->planning_algo_selection->setCurrentIndex(planOptions.planning_algo);

  ui->maxvertices_spinbox->setValue(planOptions.max_vertices);
  ui->progress_interval_spinbox->setValue(planOptions.prog_interval);
  ui->maxsolutions_spinbox->setValue(planOptions.max_results);

  if (planOptions.planning_options & pp::BIDIRECTIONAL_PLANNING)
    ui->check_bidir->setChecked(true);
  else
    ui->check_bidir->setChecked(false);

  if (planOptions.planning_options & pp::LAZY_COLLISION_CHECKING)
    ui->check_lazy_collision->setChecked(true);
  else
    ui->check_lazy_collision->setChecked(false);

  if (planOptions.planning_options & pp::PLAN_WITH_VORONOI_PULL) {
    ui->init_sa_temp_spinbox->setValue(planOptions.init_SA_temp);
    ui->check_voronoi_pull->setChecked(true);
  } else {
    ui->init_sa_temp_spinbox->setValue(0.0);
    ui->check_voronoi_pull->setChecked(false);
  };

  if (planOptions.planning_options & pp::PLAN_WITH_ANYTIME_HEURISTIC) {
    ui->init_relax_spinbox->setValue(planOptions.init_relax);
    ui->check_anytime_heuristic->setChecked(true);
  } else {
    ui->init_relax_spinbox->setValue(0.0);
    ui->check_anytime_heuristic->setChecked(false);
  };

  if (planOptions.planning_options & pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG)
    ui->check_bnb->setChecked(true);
  else
    ui->check_bnb->setChecked(false);

  ui->start_delay_spinbox->setValue(planOptions.start_delay);

  if (planOptions.store_policy == pp::DVP_ADJ_LIST_MOTION_GRAPH) {
    ui->graph_storage_selection->setCurrentIndex(1);
  } else {
    ui->graph_storage_selection->setCurrentIndex(0);
  };

  switch (planOptions.knn_method) {
    case pp::DVP_BF2_TREE_KNN:
      ui->KNN_method_selection->setCurrentIndex(1);
      break;
    case pp::DVP_BF4_TREE_KNN:
      ui->KNN_method_selection->setCurrentIndex(2);
      break;
    case pp::DVP_COB2_TREE_KNN:
      ui->KNN_method_selection->setCurrentIndex(3);
      break;
    case pp::DVP_COB4_TREE_KNN:
      ui->KNN_method_selection->setCurrentIndex(4);
      break;
    default:
      ui->KNN_method_selection->setCurrentIndex(0);
      break;
  };

  connect(ui->actionValuesChanged, SIGNAL(triggered()), this,
          SLOT(onConfigsChanged()));

  onUpdateAvailableOptions(ui->planning_algo_selection->currentIndex());
};

void PlannerAlgConfigWidget::onUpdateAvailableOptions(int plan_alg) {

  switch (plan_alg) {
    case 0:  // RRT
    case 1:  // RRT*
      ui->max_expansion_spinbox->setEnabled(false);
      break;
    case 2:  // PRM
    case 3:  // SBA*
    case 4:  // FADPRM
    default:
      ui->max_expansion_spinbox->setEnabled(true);
      break;
  };

  switch (plan_alg) {
    case 1:  // RRT*
      ui->check_lazy_collision->setEnabled(false);
      ui->check_lazy_collision->setChecked(true);
      break;
    case 3:  // SBA*
      ui->check_lazy_collision->setEnabled(true);
      ui->check_lazy_collision->setChecked(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      ui->check_lazy_collision->setEnabled(false);
      ui->check_lazy_collision->setChecked(false);
      break;
  };

  switch (plan_alg) {
    case 0:  // RRT
      ui->check_bidir->setEnabled(true);
      break;
    case 1:  // RRT*
    case 2:  // PRM
    case 3:  // SBA*
    case 4:  // FADPRM
    default:
      ui->check_bidir->setEnabled(false);
      ui->check_bidir->setChecked(false);
      break;
  };

  switch (plan_alg) {
    case 3:  // SBA*
      ui->check_voronoi_pull->setEnabled(true);
      ui->check_anytime_heuristic->setEnabled(true);
      break;
    case 0:  // RRT
    case 1:  // RRT*
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      ui->check_voronoi_pull->setEnabled(false);
      ui->check_voronoi_pull->setChecked(false);
      ui->check_anytime_heuristic->setEnabled(false);
      ui->check_anytime_heuristic->setChecked(false);
      break;
  };

  switch (plan_alg) {
    case 1:  // RRT*
    case 3:  // SBA*
      ui->check_bnb->setEnabled(true);
      break;
    case 0:  // RRT
    case 2:  // PRM
    case 4:  // FADPRM
    default:
      ui->check_bnb->setEnabled(false);
      ui->check_bnb->setChecked(false);
      break;
  };
};

void PlannerAlgConfigWidget::savePlannerConfig() {
  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save Planner Configurations..."), last_used_path,
      tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb "
         "*.raplan.pbuf)"));

  if (fileName == tr(""))
    return;

  last_used_path = QFileInfo(fileName).absolutePath();

  savePlannerConfiguration(fileName.toStdString());
};

void PlannerAlgConfigWidget::loadPlannerConfig() {
  QString fileName = QFileDialog::getOpenFileName(
      this, tr("Open Planner Configurations..."), last_used_path,
      tr("Robot-Airship Planner Configurations (*.raplan.rkx *.raplan.rkb "
         "*.raplan.pbuf)"));

  if (fileName == tr(""))
    return;

  last_used_path = QFileInfo(fileName).absolutePath();

  loadPlannerConfiguration(fileName.toStdString());
};

void PlannerAlgConfigWidget::savePlannerConfiguration(
    const std::string& aFilename) {

  onConfigsChanged();

  try {
    (*serialization::open_oarchive(aFilename)) << planOptions;
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  };
};

void PlannerAlgConfigWidget::loadPlannerConfiguration(
    const std::string& aFilename) {

  try {
    (*serialization::open_iarchive(aFilename)) >> planOptions;
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  };

  updateConfigs();
};

bool PlannerAlgConfigWidget::outputTiming() const {
  return ui->check_print_timing->isChecked();
};

bool PlannerAlgConfigWidget::outputNodeCounter() const {
  return ui->check_print_counter->isChecked();
};

bool PlannerAlgConfigWidget::outputMotionGraph() const {
  return ui->check_print_graph->isChecked();
};

bool PlannerAlgConfigWidget::outputSolution() const {
  return ui->check_print_best->isChecked();
};
};  // namespace qt
};  // namespace ReaK
