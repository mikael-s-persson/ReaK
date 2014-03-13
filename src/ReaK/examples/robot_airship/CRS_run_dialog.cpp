
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
  
  int mode = 0;
  if ( this->static_sim_check->isChecked() )
    mode = 0;
  if ( this->dynamic_sim_check->isChecked() )
    mode = 1;
  if ( this->live_sim_check->isChecked() )
    mode = 2;
  if ( this->live_run_check->isChecked() )
    mode = 3;
  
  this->stop_button->setEnabled(true);
  this->start_button->setEnabled(false);
  
  emit triggeredStartPlanning(mode);
  
};

void CRSRunDialogWidget::onStopPressed() {
  
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
  
  int mode = 0;
  if ( this->static_sim_check->isChecked() )
    mode = 0;
  if ( this->dynamic_sim_check->isChecked() )
    mode = 1;
  if ( this->live_sim_check->isChecked() )
    mode = 2;
  if ( this->live_run_check->isChecked() )
    mode = 3;
  
  emit triggeredLaunch(mode);
  
};

void CRSRunDialogWidget::onAbortPressed() {
  
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
  
  this->init_label->setText("Initializing... Done!");
  this->init_label->setStyleSheet("color: green;");
  this->planning_label->setText("Planning...");
  this->planning_label->setStyleSheet("color: red;");
  this->launch_label->setText("Launching...");
  this->launch_label->setStyleSheet("color: red;");
  
};

void CRSRunDialogWidget::onPlanningDone() {
  
  this->planning_label->setText("Planning... Done!");
  this->planning_label->setStyleSheet("color: green;");
  
  this->stop_button->setEnabled(false);
  this->proceed_button->setEnabled(true);
  
};

void CRSRunDialogWidget::onLaunchOpportunity() {
  
  this->launch_label->setStyleSheet("color: green;");
  
  this->proceed_button->setEnabled(true);
  
  this->flashing_button_timer.start();
  
};

void CRSRunDialogWidget::onLaunchStarted() {
  
  this->launch_label->setText("Launching... Done");
  this->launch_label->setStyleSheet("color: green;");
  
  this->proceed_button->setEnabled(false);
  
  this->flashing_button_timer.stop();
  
};

void CRSRunDialogWidget::onCaptureReached() {
  
  this->launch_label->setText("Captured!");
  this->launch_label->setStyleSheet("color: green;");
  
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
  
  
};



};

};


