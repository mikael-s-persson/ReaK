
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

#include <ReaK/gui_qt/data_stream_opt_widget.hpp>

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QMainWindow>
#include <QProcess>
#include <QScrollArea>

#if 0
#include <ReaK/core/serialization/archiver_factory.hpp>
#endif

#include <ReaK/core/recorders/data_record_options.hpp>


namespace ReaK {
  
namespace rkqt {



DataStreamOptWidget::DataStreamOptWidget(QWidget * parent, Qt::WindowFlags flags) :
  QDockWidget(tr("Predictor"), parent, flags),
  Ui::DataStreamOpt()
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  setupUi(dock_wid);
  
  connect(this->actionUpdateURI, SIGNAL(triggered()), this, SLOT(onUpdateURIAndDataOpt()));
  connect(this->Columns_edit, SIGNAL(textChanged()), this, SLOT(onUpdateFieldsAndDataOpt()));
  connect(this->URI_edit, SIGNAL(textChanged(QString)), this, SLOT(onUpdateFieldsAndDataOpt()));
  
  this->URI_edit->setPlainText(tr("file:./data.ssv"));
  
  onUpdateFieldsAndDataOpt();
  
};

DataStreamOptWidget::~DataStreamOptWidget() {
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
};



void DataStreamOptWidget::onConfigsChanged() {
  
  typedef ctrl::satellite_predictor_options sat_opt;
  
  sat_options.system_kind = 0;
  switch(this->kf_model_selection->currentIndex()) {
    case 0:  // IEKF
      sat_options.system_kind |= sat_opt::invariant;
      break;
    case 1:  // IMKF
      sat_options.system_kind |= sat_opt::invar_mom;
      break;
    case 2:  // IMKFv2
      sat_options.system_kind |= sat_opt::invar_mom2;
      break;
    case 3:  // IMKF_em
      sat_options.system_kind |= sat_opt::invar_mom_em;
      break;
    case 4:  // IMKF_emd
      sat_options.system_kind |= sat_opt::invar_mom_emd;
      break;
    case 5:  // IMKF_emdJ
      sat_options.system_kind |= sat_opt::invar_mom_emdJ;
      break;
    case 6:  // TSOSAIKF_em
      sat_options.system_kind |= sat_opt::invar_mom_em | sat_opt::TSOSAKF;
      break;
    case 7:  // TSOSAIKF_emd
      sat_options.system_kind |= sat_opt::invar_mom_emd | sat_opt::TSOSAKF;
      break;
    case 8:  // TSOSAIKF_emdJ
      sat_options.system_kind |= sat_opt::invar_mom_emdJ | sat_opt::TSOSAKF;
      break;
  };
  
  if(this->gyro_check->isChecked())
    sat_options.system_kind |= sat_opt::gyro_measures;
  if(this->IMU_check->isChecked())
    sat_options.system_kind |= sat_opt::IMU_measures;
  
  switch(this->predict_assumption_selection->currentIndex()) {
    case 0:  // No future measurements
      sat_options.predict_assumption = sat_opt::no_measurements;
      break;
    case 1:  // Maximum-likelihood Measurements
      sat_options.predict_assumption = sat_opt::most_likely_measurements;
      break;
    case 2:  // Full certainty
      sat_options.predict_assumption = sat_opt::full_certainty;
      break;
  };
  
  
  sat_options.time_step = getTimeStep();
  
  sat_options.mass = getMass();
  sat_options.inertia_tensor = getInertiaTensor();
  
  sat_options.input_disturbance = getInputDisturbance();
  sat_options.measurement_noise = getMeasurementNoise();
  
  sat_options.IMU_orientation = getIMUOrientation();
  sat_options.IMU_location = getIMULocation();
  sat_options.earth_orientation = getEarthOrientation();
  sat_options.mag_field_direction = getMagFieldDirection();
  
  sat_options.initial_motion = frame_3D<double>();
  
  sat_options.predict_time_horizon = getTimeHorizon();
  sat_options.predict_Pnorm_threshold = getPThreshold();
  
};


void DataStreamOptWidget::updateConfigs() {
  
  typedef ctrl::satellite_predictor_options sat_opt;
  
  switch(sat_options.system_kind & 15) {
    case sat_opt::invariant:  // IEKF
      this->kf_model_selection->setCurrentIndex(0);
      break;
    case sat_opt::invar_mom:  // IMKF
      this->kf_model_selection->setCurrentIndex(1);
      break;
    case sat_opt::invar_mom2:  // IMKFv2
      this->kf_model_selection->setCurrentIndex(2);
      break;
    case sat_opt::invar_mom_em:  // IMKF_em
      if(sat_options.system_kind & sat_opt::TSOSAKF)
        this->kf_model_selection->setCurrentIndex(6);
      else
        this->kf_model_selection->setCurrentIndex(3);
      break;
    case sat_opt::invar_mom_emd:  // IMKF_emd
      if(sat_options.system_kind & sat_opt::TSOSAKF)
        this->kf_model_selection->setCurrentIndex(7);
      else
        this->kf_model_selection->setCurrentIndex(4);
      break;
    case sat_opt::invar_mom_emdJ:  // IMKF_emdJ
      if(sat_options.system_kind & sat_opt::TSOSAKF)
        this->kf_model_selection->setCurrentIndex(8);
      else
        this->kf_model_selection->setCurrentIndex(5);
      break;
  };
  
  if(sat_options.system_kind & sat_opt::gyro_measures)
    this->gyro_check->setChecked(true);
  if(sat_options.system_kind & sat_opt::IMU_measures)
    this->IMU_check->setChecked(true);
  
  switch(sat_options.predict_assumption) {
    case sat_opt::no_measurements:
      this->predict_assumption_selection->setCurrentIndex(0);
      break;
    case sat_opt::most_likely_measurements:
      this->predict_assumption_selection->setCurrentIndex(1);
      break;
    case sat_opt::full_certainty:
      this->predict_assumption_selection->setCurrentIndex(2);
      break;
  };
  
  this->time_step_spin->setValue(sat_options.time_step);
  
  this->mass_spin->setValue(sat_options.mass);
  inertia_storage->inertia_tensor = sat_options.inertia_tensor;
  
  this->Qf_spin->setValue((sat_options.input_disturbance(0,0) + sat_options.input_disturbance(1,1) + sat_options.input_disturbance(2,2)) / 3.0);
  this->Qt_spin->setValue((sat_options.input_disturbance(3,3) + sat_options.input_disturbance(4,4) + sat_options.input_disturbance(5,5)) / 3.0);
  
  this->Rpos_spin->setValue((sat_options.measurement_noise(0,0) + sat_options.measurement_noise(1,1) + sat_options.measurement_noise(2,2)) / 3.0);
  this->Rang_spin->setValue((sat_options.measurement_noise(3,3) + sat_options.measurement_noise(4,4) + sat_options.measurement_noise(5,5)) / 3.0);
  if(sat_options.measurement_noise.get_col_count() > 6) {
    this->Rgyro_spin->setValue((sat_options.measurement_noise(6,6) + sat_options.measurement_noise(7,7) + sat_options.measurement_noise(8,8)) / 3.0);
    if(sat_options.measurement_noise.get_col_count() > 9) {
      this->Racc_spin->setValue((sat_options.measurement_noise(9,9) + sat_options.measurement_noise(10,10) + sat_options.measurement_noise(11,11)) / 3.0);
      this->Rmag_spin->setValue((sat_options.measurement_noise(12,12) + sat_options.measurement_noise(13,13) + sat_options.measurement_noise(14,14)) / 3.0);
    };
  };
  
  IMU_storage->IMU_orientation = sat_options.IMU_orientation;
  IMU_storage->IMU_location = sat_options.IMU_location;
  IMU_storage->earth_orientation = sat_options.earth_orientation;
  IMU_storage->mag_field_direction = sat_options.mag_field_direction;
  
  this->horizon_spin->setValue(sat_options.predict_time_horizon);
  this->Pthreshold_spin->setValue(sat_options.predict_Pnorm_threshold);
  
};



void DataStreamOptWidget::onUpdateURIAndDataOpt() {
  
};

void DataStreamOptWidget::onUpdateFieldsAndDataOpt() {
  
};


void DataStreamOptWidget::onUpdateAvailableOptions() {
  
  switch(filter_method) {
    case 0:  // IEKF
      this->IMU_check->setChecked(false);
      this->IMU_check->setEnabled(false);
      break;
    case 1:  // IMKF
    case 2:  // IMKFv2
    default:
      this->IMU_check->setEnabled(false); // always false (for now!)
      break;
  };
  
};


#if 0

static QString last_used_path;

void TargetPredConfigWidget::saveInertiaTensor() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Inertia Information..."), last_used_path,
    tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    sat_options.save_mass_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void TargetPredConfigWidget::editInertiaTensor() {
  ot_inertia_win.show();
};

void TargetPredConfigWidget::loadInertiaTensor() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Inertia Information..."), last_used_path,
    tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    sat_options.load_mass_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
};

#endif


};

};





