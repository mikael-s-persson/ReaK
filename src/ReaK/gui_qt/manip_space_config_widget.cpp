
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

#include <ReaK/gui_qt/manip_space_config_widget.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <QDockWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QScrollArea>

namespace ReaK {
  
namespace rkqt {

static QString last_used_path;


ManipSpaceConfigWidget::ManipSpaceConfigWidget(QWidget * parent, Qt::WindowFlags flags) :
                                               QDockWidget(tr("Space"), parent, flags),
                                               Ui::ManipSpaceConfig()
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  setupUi(dock_wid);
  
  connect(this->actionValuesChanged, SIGNAL(triggered()), this, SLOT(updateInternalValues()));
  connect(this->load_button, SIGNAL(clicked()), this, SLOT(loadSpaceConfig()));
  connect(this->save_button, SIGNAL(clicked()), this, SLOT(saveSpaceConfig()));
  
  updateInternalValues();
  
};

ManipSpaceConfigWidget::~ManipSpaceConfigWidget() { 
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
};


void ManipSpaceConfigWidget::updateInternalValues() {
  
  space_order = this->order_selection->currentIndex();
  interp_id = this->interp_selection->currentIndex();
  min_travel = this->min_interval_spinbox->value();
  max_travel = this->max_interval_spinbox->value();
  
  is_temporal = this->temporal_space_check->isChecked();
  is_rate_limited = this->rate_limited_check->isChecked();
  
  output_space_order = this->output_space_selection->currentIndex();
  
};

void ManipSpaceConfigWidget::updateExternalValues() {
  
  this->order_selection->setCurrentIndex(space_order);
  this->interp_selection->setCurrentIndex(interp_id);
  this->min_interval_spinbox->setValue(min_travel);
  this->max_interval_spinbox->setValue(max_travel);
  
  this->temporal_space_check->setChecked(is_temporal);
  this->rate_limited_check->setChecked(is_rate_limited);
  
  this->output_space_selection->setCurrentIndex(output_space_order);
  
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














