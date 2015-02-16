
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

#include <ReaK/planning/qt/manip_space_config_widget.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include "ui_manip_space_config.h"

#include <QDockWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QScrollArea>

namespace ReaK {
  
namespace qt {

static QString last_used_path;


ManipSpaceConfigWidget::ManipSpaceConfigWidget(QWidget * parent, Qt::WindowFlags flags) :
                                               QDockWidget(tr("Space"), parent, flags),
                                               ui(new Ui::ManipSpaceConfig())
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  ui->setupUi(dock_wid);
  
  connect(ui->actionValuesChanged, SIGNAL(triggered()), this, SLOT(updateInternalValues()));
  connect(ui->load_button, SIGNAL(clicked()), this, SLOT(loadSpaceConfig()));
  connect(ui->save_button, SIGNAL(clicked()), this, SLOT(saveSpaceConfig()));
  
  updateInternalValues();
  
};

ManipSpaceConfigWidget::~ManipSpaceConfigWidget() { 
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
  delete ui;
};


void ManipSpaceConfigWidget::updateInternalValues() {
  
  space_order = ui->order_selection->currentIndex();
  interp_id = ui->interp_selection->currentIndex();
  min_travel = ui->min_interval_spinbox->value();
  max_travel = ui->max_interval_spinbox->value();
  
  is_temporal = ui->temporal_space_check->isChecked();
  is_rate_limited = ui->rate_limited_check->isChecked();
  
  output_space_order = ui->output_space_selection->currentIndex();
  
};

void ManipSpaceConfigWidget::updateExternalValues() {
  
  ui->actionValuesChanged->disconnect(this, SLOT(updateInternalValues()));
  
  ui->order_selection->setCurrentIndex(space_order);
  ui->interp_selection->setCurrentIndex(interp_id);
  ui->min_interval_spinbox->setValue(min_travel);
  ui->max_interval_spinbox->setValue(max_travel);
  
  ui->temporal_space_check->setChecked(is_temporal);
  ui->rate_limited_check->setChecked(is_rate_limited);
  
  ui->output_space_selection->setCurrentIndex(output_space_order);
  
  connect(ui->actionValuesChanged, SIGNAL(triggered()), this, SLOT(updateInternalValues()));
  
};

void ManipSpaceConfigWidget::saveSpaceConfig() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Space Configurations..."), last_used_path,
    tr("Robot Space Configurations (*.rspace.rkx *.rspace.rkb *.rspace.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  saveSpaceConfiguration(fileName.toStdString());
};

void ManipSpaceConfigWidget::loadSpaceConfig() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Space Configurations..."), last_used_path,
    tr("Robot Space Configurations (*.rspace.rkx *.rspace.rkb *.rspace.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  loadSpaceConfiguration(fileName.toStdString());
};

void ManipSpaceConfigWidget::saveSpaceConfiguration(const std::string& aFilename) {
  
  updateInternalValues();
  
  try {
    shared_ptr< serialization::oarchive > p_ao = serialization::open_oarchive(aFilename); 
    (*p_ao) & RK_SERIAL_SAVE_WITH_NAME(space_order)
            & RK_SERIAL_SAVE_WITH_NAME(interp_id)
            & RK_SERIAL_SAVE_WITH_NAME(min_travel)
            & RK_SERIAL_SAVE_WITH_NAME(max_travel)
            & RK_SERIAL_SAVE_WITH_NAME(is_temporal)
            & RK_SERIAL_SAVE_WITH_NAME(is_rate_limited)
            & RK_SERIAL_SAVE_WITH_NAME(output_space_order);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void ManipSpaceConfigWidget::loadSpaceConfiguration(const std::string& aFilename) {
  
  try {
    shared_ptr< serialization::iarchive > p_ai = serialization::open_iarchive(aFilename); 
    (*p_ai) & RK_SERIAL_LOAD_WITH_NAME(space_order)
            & RK_SERIAL_LOAD_WITH_NAME(interp_id)
            & RK_SERIAL_LOAD_WITH_NAME(min_travel)
            & RK_SERIAL_LOAD_WITH_NAME(max_travel)
            & RK_SERIAL_LOAD_WITH_NAME(is_temporal)
            & RK_SERIAL_LOAD_WITH_NAME(is_rate_limited)
            & RK_SERIAL_LOAD_WITH_NAME(output_space_order);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateExternalValues();
  
};


};

};














