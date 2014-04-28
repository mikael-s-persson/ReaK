
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

#include <QDialog>
#include <QLineEdit>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>


namespace ReaK {
  
namespace rkqt {


CRSRunDialogWidget::CRSRunDialogWidget(QWidget * parent, Qt::WindowFlags flags) :
  QDialog(parent, flags),
  Ui::CRSRunDialog(),
  flashing_button_timer(this) {
  setupUi(this);
  
  connect(this->start_button, SIGNAL(clicked()), this, SLOT(onStartPressed()));
  connect(this->stop_button, SIGNAL(clicked()), this, SLOT(onStopPressed()));
  connect(this->proceed_button, SIGNAL(clicked()), this, SLOT(onLaunchPressed()));
  connect(this->abort_button, SIGNAL(clicked()), this, SLOT(onAbortPressed()));
  
  this->flashing_button_timer.setInterval(500);
  this->flashing_button_timer.setSingleShot(false);
  connect(&(this->flashing_button_timer), SIGNAL(timeout()), this, SLOT(flipLaunchButtonColor()));
  
};

CRSRunDialogWidget::~CRSRunDialogWidget() { };


void CRSRunDialogWidget::onStartPressed() {
  
  this->init_label->setText("Initializing...");
  this->init_label->setStyleSheet("color: red;");
  
  this->static_sim_check->setEnabled(false);
  this->dynamic_sim_check->setEnabled(false);
  this->live_sim_check->setEnabled(false);
  this->live_run_check->setEnabled(false);
  
  std::string mes = "Starting ";
  
  int mode = 0;
  if ( this->static_sim_check->isChecked() ) {
    mode = 0;
    mes += "static simulation...\n";
  };
  if ( this->dynamic_sim_check->isChecked() ) {
    mode = 1;
    mes += "dynamic simulation...\n";
  };
  if ( this->live_sim_check->isChecked() ) {
    mode = 2;
    mes += "live simulation...\n";
  };
  if ( this->live_run_check->isChecked() ) {
    mode = 3;
    mes += "live run...\n";
  };
  
  this->stop_button->setEnabled(true);
  this->start_button->setEnabled(false);
  
  this->publishConsoleMessage(mes);
  
  emit triggeredStartPlanning(mode);
  
};

void CRSRunDialogWidget::onStopPressed() {
  
  this->publishConsoleMessage("Stopping planner...\n");
  
  emit triggeredStopPlanning();
  
  this->static_sim_check->setEnabled(true);
  this->dynamic_sim_check->setEnabled(true);
  this->live_sim_check->setEnabled(true);
  this->live_run_check->setEnabled(true);
  
  this->stop_button->setEnabled(false);
  this->start_button->setEnabled(true);
  this->proceed_button->setEnabled(false);
  
};

void CRSRunDialogWidget::onLaunchPressed() {
  
  std::string mes = "Starting ";
  
  int mode = 0;
  if ( this->static_sim_check->isChecked() ) {
    mode = 0;
    mes += "solution animation...\n";
  };
  if ( this->dynamic_sim_check->isChecked() ) {
    mode = 1;
    mes += "solution animation...\n";
  };
  if ( this->live_sim_check->isChecked() ) {
    mode = 2;
    mes += "solution animation...\n";
  };
  if ( this->live_run_check->isChecked() ) {
    mode = 3;
    mes += "solution live execution...\nKEEP HAND ON EMERGENCY STOP BUTTON!!!\n";
  };
  
  this->publishConsoleMessage(mes);
  
  emit triggeredLaunch(mode);
  
};

void CRSRunDialogWidget::onAbortPressed() {
  
  this->publishConsoleMessage("ABORT! ABORT! ABORT!\n");
  
  emit triggeredAbort();
  
  onReset();
  
};


void CRSRunDialogWidget::flipLaunchButtonColor() {
  
  QString cur_ss = this->proceed_button->styleSheet();
  if( cur_ss == "color: navy; background-color: lime;" ) {
    this->proceed_button->setStyleSheet("color: navy; background-color: green;");
  } else {
    this->proceed_button->setStyleSheet("color: navy; background-color: lime;");
  };
  
};


void CRSRunDialogWidget::onInitializationDone() {
  
  this->publishConsoleMessage("Initializing Done!\n");
  
  this->init_label->setText("Initializing... Done!");
  this->init_label->setStyleSheet("color: green;");
  this->planning_label->setText("Planning...");
  this->planning_label->setStyleSheet("color: red;");
  this->launch_label->setText("Launching...");
  this->launch_label->setStyleSheet("color: red;");
  
};

void CRSRunDialogWidget::onPlanningDone() {
  
  this->publishConsoleMessage("Planning Done!\n");
  
  this->planning_label->setText("Planning... Done!");
  this->planning_label->setStyleSheet("color: green;");
  
  this->stop_button->setEnabled(false);
  this->proceed_button->setEnabled(true);
  
};

void CRSRunDialogWidget::onLaunchOpportunity() {
  
  this->publishConsoleMessage("Launch opportunity found!\n");
  
  this->launch_label->setStyleSheet("color: green;");
  
  this->proceed_button->setEnabled(true);
  
  this->flashing_button_timer.start();
  
};

void CRSRunDialogWidget::onLaunchStarted() {
  
  this->publishConsoleMessage("Executing solution...\n");
  
  this->launch_label->setText("Launching... Done");
  this->launch_label->setStyleSheet("color: green;");
  
  this->proceed_button->setEnabled(false);
  
  this->flashing_button_timer.stop();
  
};

void CRSRunDialogWidget::onCaptureReached() {
  
  this->publishConsoleMessage("Capture solution reached!\n");
  
  this->launch_label->setText("Captured!");
  this->launch_label->setStyleSheet("color: green;");
  
  this->flashing_button_timer.stop();
  
};

void CRSRunDialogWidget::onReset() {
  
  this->static_sim_check->setEnabled(true);
  this->dynamic_sim_check->setEnabled(true);
  this->live_sim_check->setEnabled(true);
  this->live_run_check->setEnabled(true);
  
  this->start_button->setEnabled(true);
  this->stop_button->setEnabled(false);
  this->proceed_button->setEnabled(false);
  
  this->init_label->setText("");
  this->init_label->setStyleSheet("color: red;");
  
  this->planning_label->setText("");
  this->planning_label->setStyleSheet("color: red;");
  
  this->launch_label->setText("");
  this->launch_label->setStyleSheet("color: red;");
  
  this->capture_label->setText("");
  this->capture_label->setStyleSheet("color: red;");
  
  this->status_text->clear();
  
  this->flashing_button_timer.stop();
  
};

void CRSRunDialogWidget::publishConsoleMessage(const std::string& aMessage) {
  this->status_text->appendPlainText(QString::fromStdString(aMessage));
};


};

};


