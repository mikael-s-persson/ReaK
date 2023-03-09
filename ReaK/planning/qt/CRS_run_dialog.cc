
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include "CRS_run_dialog.hpp"
#include "ReaK/core/base/defs.h"

#include "ui_CRS_run_dialog.h"

#include <QDialog>
#include <QFileInfo>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QScrollBar>

namespace ReaK {

namespace qt {

CRSRunDialogWidget::CRSRunDialogWidget(QWidget* parent, Qt::WindowFlags flags)
    : QDialog(parent, flags),
      ui(new Ui::CRSRunDialog()),
      flashing_button_timer(this) {
  ui->setupUi(this);

  connect(ui->start_button, SIGNAL(clicked()), this, SLOT(onStartPressed()));
  connect(ui->stop_button, SIGNAL(clicked()), this, SLOT(onStopPressed()));
  connect(ui->proceed_button, SIGNAL(clicked()), this, SLOT(onLaunchPressed()));
  connect(ui->abort_button, SIGNAL(clicked()), this, SLOT(onAbortPressed()));

  flashing_button_timer.setInterval(500);
  flashing_button_timer.setSingleShot(false);
  connect(&(flashing_button_timer), SIGNAL(timeout()), this,
          SLOT(flipLaunchButtonColor()));
};

CRSRunDialogWidget::~CRSRunDialogWidget() {
  delete ui;
};

void CRSRunDialogWidget::onStartPressed() {

  ui->init_label->setText("Initializing...");
  ui->init_label->setStyleSheet("color: red;");

  ui->static_sim_check->setEnabled(false);
  ui->dynamic_sim_check->setEnabled(false);
  ui->live_sim_check->setEnabled(false);
  ui->live_run_check->setEnabled(false);

  QString mes = "Starting ";

  int mode = 0;
  if (ui->static_sim_check->isChecked()) {
    mode = 0;
    mes += "static simulation...\n";
  };
  if (ui->dynamic_sim_check->isChecked()) {
    mode = 1;
    mes += "dynamic simulation...\n";
  };
  if (ui->live_sim_check->isChecked()) {
    mode = 2;
    mes += "live simulation...\n";
  };
  if (ui->live_run_check->isChecked()) {
    mode = 3;
    mes += "live run...\n";
  };

  ui->stop_button->setEnabled(true);
  ui->start_button->setEnabled(false);
  ui->proceed_button->setEnabled(false);

  publishConsoleMessage(mes);

  emit triggeredStartPlanning(mode);
};

void CRSRunDialogWidget::onStopPressed() {

  publishConsoleMessage("Stopping planner...\n");

  emit triggeredStopPlanning();

  ui->static_sim_check->setEnabled(true);
  ui->dynamic_sim_check->setEnabled(true);
  ui->live_sim_check->setEnabled(true);
  ui->live_run_check->setEnabled(true);

  ui->stop_button->setEnabled(false);
  ui->start_button->setEnabled(true);
  ui->proceed_button->setEnabled(false);
};

void CRSRunDialogWidget::onLaunchPressed() {

  QString mes = "Starting ";

  int mode = 0;
  if (ui->static_sim_check->isChecked()) {
    mode = 0;
    mes += "solution animation...\n";
  };
  if (ui->dynamic_sim_check->isChecked()) {
    mode = 1;
    mes += "solution animation...\n";
  };
  if (ui->live_sim_check->isChecked()) {
    mode = 2;
    mes += "solution animation...\n";
  };
  if (ui->live_run_check->isChecked()) {
    mode = 3;
    mes +=
        "solution live execution...\nKEEP HAND ON EMERGENCY STOP BUTTON!!!\n";
  };

  ui->stop_button->setEnabled(false);
  ui->start_button->setEnabled(false);

  publishConsoleMessage(mes);

  emit triggeredLaunch(mode);
};

void CRSRunDialogWidget::onAbortPressed() {

  publishConsoleMessage("ABORT! ABORT! ABORT!\n");

  emit triggeredAbort();

  onReset();
};

void CRSRunDialogWidget::flipLaunchButtonColor() {

  QString cur_ss = ui->proceed_button->styleSheet();
  if (cur_ss == "color: navy; background-color: lime;") {
    ui->proceed_button->setStyleSheet("color: navy; background-color: green;");
  } else {
    ui->proceed_button->setStyleSheet("color: navy; background-color: lime;");
  };
};

void CRSRunDialogWidget::onInitializationDone() {

  publishConsoleMessage("Initializing Done!\n");

  ui->init_label->setText("Initializing... Done!");
  ui->init_label->setStyleSheet("color: green;");
  ui->planning_label->setText("Planning...");
  ui->planning_label->setStyleSheet("color: red;");
  ui->launch_label->setText("Launching...");
  ui->launch_label->setStyleSheet("color: red;");
};

void CRSRunDialogWidget::onPlanningDone() {

  publishConsoleMessage("Planning Done!\n");

  ui->planning_label->setText("Planning... Done!");
  ui->planning_label->setStyleSheet("color: green;");

  ui->start_button->setEnabled(true);
  ui->stop_button->setEnabled(false);
};

void CRSRunDialogWidget::onLaunchOpportunity() {

  publishConsoleMessage("Launch opportunity found!\n");

  ui->launch_label->setStyleSheet("color: green;");

  ui->proceed_button->setEnabled(true);

  flashing_button_timer.start();
};

void CRSRunDialogWidget::onLaunchStarted() {

  publishConsoleMessage("Executing solution...\n");

  ui->launch_label->setText("Launching... Done");
  ui->launch_label->setStyleSheet("color: green;");

  ui->stop_button->setEnabled(false);
  ui->start_button->setEnabled(false);
  ui->proceed_button->setEnabled(false);

  flashing_button_timer.stop();
  ui->proceed_button->setStyleSheet("color: navy; background-color: lime;");
};

void CRSRunDialogWidget::onCaptureReached() {

  publishConsoleMessage("Capture solution reached!\n");

  ui->launch_label->setText("Captured!");
  ui->launch_label->setStyleSheet("color: green;");

  ui->stop_button->setEnabled(false);
  ui->start_button->setEnabled(true);
  ui->proceed_button->setEnabled(false);
  ui->proceed_button->setStyleSheet("color: navy; background-color: lime;");

  flashing_button_timer.stop();
};

void CRSRunDialogWidget::onReset() {

  ui->static_sim_check->setEnabled(true);
  ui->dynamic_sim_check->setEnabled(true);
  ui->live_sim_check->setEnabled(true);
  ui->live_run_check->setEnabled(true);

  ui->start_button->setEnabled(true);
  ui->stop_button->setEnabled(false);
  ui->proceed_button->setEnabled(false);
  ui->proceed_button->setStyleSheet("color: navy; background-color: lime;");

  ui->init_label->setText("");
  ui->init_label->setStyleSheet("color: red;");

  ui->planning_label->setText("");
  ui->planning_label->setStyleSheet("color: red;");

  ui->launch_label->setText("");
  ui->launch_label->setStyleSheet("color: red;");

  ui->capture_label->setText("");
  ui->capture_label->setStyleSheet("color: red;");

  //   ui->status_text->clear();

  publishConsoleMessage("\nReset complete! Ready to work...\n");

  flashing_button_timer.stop();
};

void CRSRunDialogWidget::publishConsoleMessage(QString aMessage) {
  ui->status_text->appendPlainText(aMessage);
  ui->status_text->verticalScrollBar()->setSliderPosition(
      ui->status_text->verticalScrollBar()->maximum());
};
};  // namespace qt
};  // namespace ReaK
