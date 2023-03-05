
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

#include <ReaK/control/qt/target_pred_config_widget.hpp>

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMainWindow>
#include <QMessageBox>
#include <QProcess>
#include <QScrollArea>

#include <ReaK/core/serialization/archiver_factory.hpp>
#include <atomic>
#include <functional>

#include <ReaK/control/estimators/belief_state_predictor.hpp>
#include <ReaK/control/estimators/covar_topology.hpp>
#include <ReaK/control/estimators/gaussian_belief_space.hpp>
#include <ReaK/control/estimators/maximum_likelihood_mapping.hpp>
#include <ReaK/control/estimators/tsos_aug_inv_kalman_filter.hpp>
#include <ReaK/control/systems/augmented_to_state_mapping.hpp>
#include <ReaK/topologies/interpolation/constant_trajectory.hpp>
#include <ReaK/topologies/interpolation/transformed_trajectory.hpp>

#include <ReaK/core/recorders/data_record_options.hpp>

#include "ui_target_predictor_config.h"

#include <ctime>
#include <tuple>
#include <type_traits>

namespace ReaK::qt {

namespace detail {

struct inertia_tensor_storage_impl : public shared_object {

  mat<double, mat_structure::symmetric> inertia_tensor;

  inertia_tensor_storage_impl()
      : inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0) {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(inertia_tensor);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(inertia_tensor);
  }
  RK_RTTI_MAKE_CONCRETE_1BASE(inertia_tensor_storage_impl, 0xBEEF0001, 1,
                              "inertia_tensor_storage_impl", shared_object)
};

struct IMU_config_storage_impl : public shared_object {

  unit_quat<double> IMU_orientation;
  vect<double, 3> IMU_location;
  unit_quat<double> earth_orientation;
  vect<double, 3> mag_field_direction;

  IMU_config_storage_impl()
      : IMU_orientation(),
        IMU_location(),
        earth_orientation(),
        mag_field_direction(1.0, 0.0, 0.0) {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(IMU_orientation) &
        RK_SERIAL_SAVE_WITH_NAME(IMU_location) &
        RK_SERIAL_SAVE_WITH_NAME(earth_orientation) &
        RK_SERIAL_SAVE_WITH_NAME(mag_field_direction);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(IMU_orientation) &
        RK_SERIAL_LOAD_WITH_NAME(IMU_location) &
        RK_SERIAL_LOAD_WITH_NAME(earth_orientation) &
        RK_SERIAL_LOAD_WITH_NAME(mag_field_direction);
  }
  RK_RTTI_MAKE_CONCRETE_1BASE(IMU_config_storage_impl, 0xBEEF0002, 1,
                              "IMU_config_storage_impl", shared_object)
};
}  // namespace detail

static QString last_used_path;

TargetPredConfigWidget::TargetPredConfigWidget(
    std::shared_ptr<trajectory_type>* aTargetAnimTraj,
    std::atomic<double>* aCurrentTargetAnimTime, QWidget* parent,
    Qt::WindowFlags flags)
    : QDockWidget(tr("Predictor"), parent, flags),
      ui(new Ui::TargetPredConfig()),
      inertia_storage(new detail::inertia_tensor_storage_impl()),
      IMU_storage(new detail::IMU_config_storage_impl()),
      objtree_sch_bld(),
      ot_inertia_graph(std::shared_ptr<serialization::object_graph>(
          new serialization::object_graph())),
      ot_inertia_root(add_vertex(*ot_inertia_graph)),
      ot_inertia_win(this, Qt::WindowFlags(Qt::Popup | Qt::Dialog)),
      ot_IMU_graph(std::shared_ptr<serialization::object_graph>(
          new serialization::object_graph())),
      ot_IMU_root(add_vertex(*ot_IMU_graph)),
      ot_IMU_win(this, Qt::WindowFlags(Qt::Popup | Qt::Dialog)),
      target_anim_traj(aTargetAnimTraj),
      current_target_anim_time(aCurrentTargetAnimTime) {
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  ui->setupUi(dock_wid);

  connect(ui->kf_model_selection, SIGNAL(currentIndexChanged(int)), this,
          SLOT(onUpdateAvailableOptions(int)));
  connect(ui->actionValuesChanged, SIGNAL(triggered()), this,
          SLOT(onConfigsChanged()));
  connect(ui->load_button, SIGNAL(clicked()), this,
          SLOT(loadPredictorConfig()));
  connect(ui->save_button, SIGNAL(clicked()), this,
          SLOT(savePredictorConfig()));

  connect(ui->I_load_button, SIGNAL(clicked()), this,
          SLOT(loadInertiaTensor()));
  connect(ui->I_edit_button, SIGNAL(clicked()), this,
          SLOT(editInertiaTensor()));
  connect(ui->I_save_button, SIGNAL(clicked()), this,
          SLOT(saveInertiaTensor()));

  connect(ui->IMU_load_button, SIGNAL(clicked()), this, SLOT(loadIMUConfig()));
  connect(ui->IMU_edit_button, SIGNAL(clicked()), this, SLOT(editIMUConfig()));
  connect(ui->IMU_save_button, SIGNAL(clicked()), this, SLOT(saveIMUConfig()));

  objtree_sch_bld << inertia_storage << IMU_storage;

  {
    ot_inertia_widget = new ObjectTreeWidget(ot_inertia_graph, ot_inertia_root);
    ot_inertia_propedit = new PropEditorWidget(&(ot_inertia_widget->mdl));
    ot_inertia_edit = &(ot_inertia_propedit->mdl.get_object_editor());
    ot_inertia_edit->add_new_object(inertia_storage);
    ot_inertia_widget->mdl.refreshObjTree();

    ot_inertia_win.resize(400, 500);
    ot_inertia_win.move(100, 100);
    ot_inertia_win.setWindowTitle("Edit Inertia Tensor");
    ot_inertia_win.addDockWidget(Qt::RightDockWidgetArea, ot_inertia_widget);
    ot_inertia_win.addDockWidget(Qt::RightDockWidgetArea, ot_inertia_propedit);
  }

  {
    ot_IMU_widget = new ObjectTreeWidget(ot_IMU_graph, ot_IMU_root);
    ot_IMU_propedit = new PropEditorWidget(&(ot_IMU_widget->mdl));
    ot_IMU_edit = &(ot_IMU_propedit->mdl.get_object_editor());
    ot_IMU_edit->add_new_object(IMU_storage);
    ot_IMU_widget->mdl.refreshObjTree();

    ot_IMU_win.resize(400, 500);
    ot_IMU_win.move(100, 100);
    ot_IMU_win.setWindowTitle("Edit IMU Configurations");
    ot_IMU_win.addDockWidget(Qt::RightDockWidgetArea, ot_IMU_widget);
    ot_IMU_win.addDockWidget(Qt::RightDockWidgetArea, ot_IMU_propedit);
  }

  meas_out_opt.kind = recorder::data_stream_options::space_separated;
  meas_out_opt.file_name = "exp_results/robot_airship/";
  est_out_opt = meas_out_opt;
  pred_out_opt = meas_out_opt;
  meas_out_opt.file_name += std::string("meas_$d_$t.ssv");
  est_out_opt.file_name += std::string("est_$d_$t.ssv");
  pred_out_opt.file_name += std::string("pred_$d_$t.ssv");

  updateConfigs();
}

TargetPredConfigWidget::~TargetPredConfigWidget() {
  stopStatePrediction();

  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
  delete ui;
}

double TargetPredConfigWidget::getTimeStep() const {
  return ui->time_step_spin->value();
}

double TargetPredConfigWidget::getMass() const {
  return ui->mass_spin->value();
}

const mat<double, mat_structure::symmetric>&
TargetPredConfigWidget::getInertiaTensor() const {
  return inertia_storage->inertia_tensor;
}

mat<double, mat_structure::diagonal>
TargetPredConfigWidget::getInputDisturbance() const {
  mat<double, mat_structure::diagonal> result(6, true);
  result(2, 2) = result(1, 1) = result(0, 0) = ui->Qf_spin->value();
  result(5, 5) = result(4, 4) = result(3, 3) = ui->Qt_spin->value();
  return result;
}

mat<double, mat_structure::diagonal>
TargetPredConfigWidget::getMeasurementNoise() const {

  std::size_t m_noise_size = 6;
  if (ui->gyro_check->isChecked()) {
    m_noise_size = 9;
  }
  if (ui->IMU_check->isChecked()) {
    m_noise_size = 15;
  }

  mat<double, mat_structure::diagonal> result(m_noise_size, true);
  result(2, 2) = result(1, 1) = result(0, 0) = ui->Rpos_spin->value();
  result(5, 5) = result(4, 4) = result(3, 3) = ui->Rang_spin->value();
  if (ui->gyro_check->isChecked()) {
    result(8, 8) = result(7, 7) = result(6, 6) = ui->Rgyro_spin->value();
  }
  if (ui->IMU_check->isChecked()) {
    result(11, 11) = result(10, 10) = result(9, 9) = ui->Racc_spin->value();
    result(14, 14) = result(13, 13) = result(12, 12) = ui->Rmag_spin->value();
  }

  return result;
}

const unit_quat<double>& TargetPredConfigWidget::getIMUOrientation() const {
  return IMU_storage->IMU_orientation;
}

const vect<double, 3>& TargetPredConfigWidget::getIMULocation() const {
  return IMU_storage->IMU_location;
}

const unit_quat<double>& TargetPredConfigWidget::getEarthOrientation() const {
  return IMU_storage->earth_orientation;
}

const vect<double, 3>& TargetPredConfigWidget::getMagFieldDirection() const {
  return IMU_storage->mag_field_direction;
}

double TargetPredConfigWidget::getTimeHorizon() const {
  return ui->horizon_spin->value();
}

double TargetPredConfigWidget::getPThreshold() const {
  return ui->Pthreshold_spin->value();
}

std::string TargetPredConfigWidget::getServerAddress() const {
  return ui->ip_addr_edit->text().toStdString();
}

int TargetPredConfigWidget::getPortNumber() const {
  return ui->port_spin->value();
}

bool TargetPredConfigWidget::useRawUDP() const {
  return ui->raw_udp_radio->isChecked();
}

bool TargetPredConfigWidget::useUDP() const {
  return ui->udp_radio->isChecked();
}

bool TargetPredConfigWidget::useTCP() const {
  return ui->tcp_radio->isChecked();
}

std::string TargetPredConfigWidget::getStartScript() const {
  return ui->start_script_edit->text().toStdString();
}

void TargetPredConfigWidget::onConfigsChanged() {

  using sat_opt = ctrl::satellite_predictor_options;

  sat_options.system_kind = 0;
  switch (ui->kf_model_selection->currentIndex()) {
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
  }

  if (ui->gyro_check->isChecked()) {
    sat_options.system_kind |= sat_opt::gyro_measures;
  }
  if (ui->IMU_check->isChecked()) {
    sat_options.system_kind |= sat_opt::IMU_measures;
  }

  switch (ui->predict_assumption_selection->currentIndex()) {
    case 0:  // No future measurements
      sat_options.predict_assumption = sat_opt::no_measurements;
      break;
    case 1:  // Maximum-likelihood Measurements
      sat_options.predict_assumption = sat_opt::most_likely_measurements;
      break;
    case 2:  // Full certainty
      sat_options.predict_assumption = sat_opt::full_certainty;
      break;
  }

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
}

void TargetPredConfigWidget::updateConfigs() {

  ui->kf_model_selection->disconnect(this, SLOT(onUpdateAvailableOptions(int)));
  ui->actionValuesChanged->disconnect(this, SLOT(onConfigsChanged()));

  using sat_opt = ctrl::satellite_predictor_options;

  switch (sat_options.system_kind & 15) {
    case sat_opt::invariant:  // IEKF
      ui->kf_model_selection->setCurrentIndex(0);
      break;
    case sat_opt::invar_mom:  // IMKF
      ui->kf_model_selection->setCurrentIndex(1);
      break;
    case sat_opt::invar_mom2:  // IMKFv2
      ui->kf_model_selection->setCurrentIndex(2);
      break;
    case sat_opt::invar_mom_em:  // IMKF_em
      if (sat_options.system_kind & sat_opt::TSOSAKF)
        ui->kf_model_selection->setCurrentIndex(6);
      else
        ui->kf_model_selection->setCurrentIndex(3);
      break;
    case sat_opt::invar_mom_emd:  // IMKF_emd
      if (sat_options.system_kind & sat_opt::TSOSAKF)
        ui->kf_model_selection->setCurrentIndex(7);
      else
        ui->kf_model_selection->setCurrentIndex(4);
      break;
    case sat_opt::invar_mom_emdJ:  // IMKF_emdJ
      if (sat_options.system_kind & sat_opt::TSOSAKF)
        ui->kf_model_selection->setCurrentIndex(8);
      else
        ui->kf_model_selection->setCurrentIndex(5);
      break;
  }

  ui->gyro_check->setChecked(sat_options.system_kind & sat_opt::gyro_measures);
  ui->IMU_check->setChecked((sat_options.system_kind & sat_opt::IMU_measures) ==
                            sat_opt::IMU_measures);

  switch (sat_options.predict_assumption) {
    case sat_opt::no_measurements:
      ui->predict_assumption_selection->setCurrentIndex(0);
      break;
    case sat_opt::most_likely_measurements:
      ui->predict_assumption_selection->setCurrentIndex(1);
      break;
    case sat_opt::full_certainty:
      ui->predict_assumption_selection->setCurrentIndex(2);
      break;
  }

  ui->time_step_spin->setValue(sat_options.time_step);

  ui->mass_spin->setValue(sat_options.mass);
  inertia_storage->inertia_tensor = sat_options.inertia_tensor;
  ot_inertia_widget->mdl.refreshObjTree();

  ui->Qf_spin->setValue((sat_options.input_disturbance(0, 0) +
                         sat_options.input_disturbance(1, 1) +
                         sat_options.input_disturbance(2, 2)) /
                        3.0);
  ui->Qt_spin->setValue((sat_options.input_disturbance(3, 3) +
                         sat_options.input_disturbance(4, 4) +
                         sat_options.input_disturbance(5, 5)) /
                        3.0);

  ui->Rpos_spin->setValue((sat_options.measurement_noise(0, 0) +
                           sat_options.measurement_noise(1, 1) +
                           sat_options.measurement_noise(2, 2)) /
                          3.0);
  ui->Rang_spin->setValue((sat_options.measurement_noise(3, 3) +
                           sat_options.measurement_noise(4, 4) +
                           sat_options.measurement_noise(5, 5)) /
                          3.0);
  if (sat_options.measurement_noise.get_col_count() > 6) {
    ui->Rgyro_spin->setValue((sat_options.measurement_noise(6, 6) +
                              sat_options.measurement_noise(7, 7) +
                              sat_options.measurement_noise(8, 8)) /
                             3.0);
    if (sat_options.measurement_noise.get_col_count() > 9) {
      ui->Racc_spin->setValue((sat_options.measurement_noise(9, 9) +
                               sat_options.measurement_noise(10, 10) +
                               sat_options.measurement_noise(11, 11)) /
                              3.0);
      ui->Rmag_spin->setValue((sat_options.measurement_noise(12, 12) +
                               sat_options.measurement_noise(13, 13) +
                               sat_options.measurement_noise(14, 14)) /
                              3.0);
    }
  }

  IMU_storage->IMU_orientation = sat_options.IMU_orientation;
  IMU_storage->IMU_location = sat_options.IMU_location;
  IMU_storage->earth_orientation = sat_options.earth_orientation;
  IMU_storage->mag_field_direction = sat_options.mag_field_direction;
  ot_IMU_widget->mdl.refreshObjTree();

  ui->horizon_spin->setValue(sat_options.predict_time_horizon);
  ui->Pthreshold_spin->setValue(sat_options.predict_Pnorm_threshold);

  connect(ui->kf_model_selection, SIGNAL(currentIndexChanged(int)), this,
          SLOT(onUpdateAvailableOptions(int)));
  connect(ui->actionValuesChanged, SIGNAL(triggered()), this,
          SLOT(onConfigsChanged()));
}

void TargetPredConfigWidget::onUpdateAvailableOptions(int filter_method) {
  switch (filter_method) {
    case 0:  // IEKF
      ui->IMU_check->setChecked(false);
      ui->IMU_check->setEnabled(false);
      break;
    case 1:  // IMKF
    case 2:  // IMKFv2
    default:
      ui->IMU_check->setEnabled(false);  // always false (for now!)
      break;
  }
}

void TargetPredConfigWidget::savePredictorConfig() {
  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save Predictor Configurations..."), last_used_path,
      tr("Target Predictor Configurations (*.tpred.rkx *.tpred.rkb "
         "*.tpred.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  savePredictorConfigurations(fileName.toStdString());
}

void TargetPredConfigWidget::loadPredictorConfig() {
  QString fileName = QFileDialog::getOpenFileName(
      this, tr("Open Predictor Configurations..."), last_used_path,
      tr("Target Predictor Configurations (*.tpred.rkx *.tpred.rkb "
         "*.tpred.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  loadPredictorConfigurations(fileName.toStdString());
}

void TargetPredConfigWidget::savePredictorConfigurations(
    const std::string& aFilename) {

  onConfigsChanged();

  try {
    sat_options.save_all_configs(aFilename);
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }
}

void TargetPredConfigWidget::loadPredictorConfigurations(
    const std::string& aFilename) {

  try {
    sat_options.load_all_configs(aFilename);
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }

  updateConfigs();
}

ctrl::satellite_predictor_options
TargetPredConfigWidget::getSatPredictorOptions() const {
  return sat_options;
}

void TargetPredConfigWidget::saveInertiaTensor() {
  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save Inertia Information..."), last_used_path,
      tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  onConfigsChanged();

  try {
    sat_options.save_mass_configs(fileName.toStdString());
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }
}

void TargetPredConfigWidget::editInertiaTensor() {
  ot_inertia_win.show();
}

void TargetPredConfigWidget::loadInertiaTensor() {
  QString fileName = QFileDialog::getOpenFileName(
      this, tr("Open Inertia Information..."), last_used_path,
      tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  try {
    sat_options.load_mass_configs(fileName.toStdString());
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }

  updateConfigs();
}

void TargetPredConfigWidget::saveIMUConfig() {
  QString fileName = QFileDialog::getSaveFileName(
      this, tr("Save IMU Configurations..."), last_used_path,
      tr("Target IMU Configurations (*.rkx *.xml *.rkb *.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  onConfigsChanged();

  try {
    sat_options.save_IMU_configs(fileName.toStdString());
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }
}

void TargetPredConfigWidget::editIMUConfig() {
  ot_IMU_win.show();
}

void TargetPredConfigWidget::loadIMUConfig() {
  QString fileName = QFileDialog::getOpenFileName(
      this, tr("Open IMU Configurations..."), last_used_path,
      tr("Target IMU Configurations (*.rkx *.xml *.rkb *.pbuf)"));

  if (fileName == tr("")) {
    return;
  }

  last_used_path = QFileInfo(fileName).absolutePath();

  try {
    sat_options.load_IMU_configs(fileName.toStdString());
  } catch (...) {
    QMessageBox::information(this, "File Type Not Supported!",
                             "Sorry, this file-type is not supported!",
                             QMessageBox::Ok);
    return;
  }

  updateConfigs();
}

namespace {

using sat3D_state_type = ReaK::ctrl::satellite_model_options::state_type;
using sat3D_cov_type = ReaK::ctrl::satellite_model_options::covar_type;
using sat3D_cov_matrix_type = sat3D_cov_type::matrix_type;

const sat3D_state_type& get_sat3D_state(const sat3D_state_type& x) {
  return x;
}

template <typename StateTuple>
const sat3D_state_type& get_sat3D_state(const StateTuple& x) {
  using ReaK::get;
  return get<0>(x);
}

struct sat3D_meas_est_pred_to_recorders {
  std::shared_ptr<ReaK::recorder::data_recorder> meas_rec;
  std::shared_ptr<ReaK::recorder::data_recorder> est_rec;
  std::shared_ptr<ReaK::recorder::data_recorder> pred_rec;

  sat3D_meas_est_pred_to_recorders(
      const std::shared_ptr<ReaK::recorder::data_recorder>& aMeasRec,
      const std::shared_ptr<ReaK::recorder::data_recorder>& aEstRec,
      const std::shared_ptr<ReaK::recorder::data_recorder>& aPredRec)
      : meas_rec(aMeasRec), est_rec(aEstRec), pred_rec(aPredRec) {}
  ~sat3D_meas_est_pred_to_recorders() {
    (*meas_rec) << ReaK::recorder::data_recorder::flush;
    (*est_rec) << ReaK::recorder::data_recorder::flush;
    (*pred_rec) << ReaK::recorder::data_recorder::flush;
  }

  template <typename BeliefStateType, typename InputBeliefType,
            typename OutputBeliefType>
  void add_record(const BeliefStateType& b, const BeliefStateType& b_pred,
                  const InputBeliefType& b_u, const OutputBeliefType& b_z,
                  double time) {
    using namespace ReaK;
    using ReaK::to_vect;

    const sat3D_state_type& x_mean = get_sat3D_state(b.get_mean_state());
    (*est_rec) << time << get_position(x_mean) << get_quaternion(x_mean)
               << get_velocity(x_mean) << get_ang_velocity(x_mean);

    vect_n<double> all_x = to_vect<double>(b.get_mean_state());
    for (std::size_t l = 13; l < all_x.size(); ++l) {
      (*est_rec) << all_x[l];
    }

    const sat3D_state_type& x_pred_mean =
        get_sat3D_state(b_pred.get_mean_state());
    (*pred_rec) << time << get_position(x_pred_mean)
                << get_quaternion(x_pred_mean) << get_velocity(x_pred_mean)
                << get_ang_velocity(x_pred_mean);

    vect_n<double> all_x_pred = to_vect<double>(b_pred.get_mean_state());
    for (std::size_t l = 13; l < all_x_pred.size(); ++l) {
      (*pred_rec) << all_x_pred[l];
    }

    const vect_n<double>& z = b_z.get_mean_state();
    (*meas_rec) << time << z << b_u.get_mean_state();

    axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) *
                               quaternion<double>(z[range(3, 7)]));
    (*est_rec) << (get_position(x_mean) - z[range(0, 3)])
               << (aa_diff.angle() * aa_diff.axis())
               << vect<double, 3>(0.0, 0.0, 0.0);
    axis_angle<double> aa_pred_diff(
        invert(get_quaternion(x_pred_mean).as_rotation()) *
        quaternion<double>(z[range(3, 7)]));
    (*pred_rec) << (get_position(x_pred_mean) - z[range(0, 3)])
                << (aa_pred_diff.angle() * aa_pred_diff.axis())
                << vect<double, 3>(0.0, 0.0, 0.0);
    if (z.size() >= 10) {
      (*est_rec) << (get_ang_velocity(x_mean) - z[range(7, 10)]);
      (*pred_rec) << (get_ang_velocity(x_pred_mean) - z[range(7, 10)]);
    } else {
      (*est_rec) << vect<double, 3>(0.0, 0.0, 0.0);
      (*pred_rec) << vect<double, 3>(0.0, 0.0, 0.0);
    }

    const sat3D_cov_matrix_type& P_xx = b.get_covariance().get_matrix();
    for (std::size_t l = 0; l < P_xx.get_row_count(); ++l) {
      (*est_rec) << P_xx(l, l);
    }

    const sat3D_cov_matrix_type& P_pred_xx =
        b_pred.get_covariance().get_matrix();
    for (std::size_t l = 0; l < P_pred_xx.get_row_count(); ++l) {
      (*pred_rec) << P_pred_xx(l, l);
    }

    (*meas_rec) << recorder::data_recorder::end_value_row;
    (*est_rec) << recorder::data_recorder::end_value_row;
    (*pred_rec) << recorder::data_recorder::end_value_row;
  }
};

static std::atomic<bool> prediction_should_stop(false);
static std::thread prediction_executer = std::thread();

template <typename Sat3DSystemType>
struct prediction_updater {

  using TempSpaceType = typename Sat3DSystemType::temporal_state_space_type;
  using BeliefSpaceType = typename Sat3DSystemType::belief_space_type;
  using TempBeliefSpaceType =
      typename Sat3DSystemType::temporal_belief_space_type;
  using CovarType = typename Sat3DSystemType::covar_type;
  using CovarMatType = typename CovarType::matrix_type;

  using StateBeliefType = typename Sat3DSystemType::state_belief_type;
  using InputBeliefType = typename Sat3DSystemType::input_belief_type;
  using OutputBeliefType = typename Sat3DSystemType::output_belief_type;

  using TempBeliefPointType =
      typename pp::topology_traits<TempBeliefSpaceType>::point_type;

  using InputTrajType =
      pp::constant_trajectory<pp::vector_topology<vect_n<double>>>;

  using PredFactoryType = typename ctrl::try_TSOSAIKF_belief_transfer_factory<
      Sat3DSystemType>::type;
  using BeliefPredTrajType =
      ctrl::belief_predicted_trajectory<BeliefSpaceType, PredFactoryType,
                                        InputTrajType>;

  std::promise<std::shared_ptr<BeliefPredTrajType>>
      predictor_promise; /* shared, sync'd by the promise-future mechanism */
  std::shared_ptr<Sat3DSystemType> satellite3D_system; /* not shared */
  ctrl::satellite_predictor_options
      sat_options; /* not shared, don't care that it's expensive to copy, I want to avoid shared states */

  std::shared_ptr<recorder::data_extractor> data_in; /* not shared */
  sat3D_meas_est_pred_to_recorders data_logger;      /* not shared */

  std::atomic<double>* current_target_anim_time; /* shared, atomic operations */

  prediction_updater(
      std::promise<std::shared_ptr<BeliefPredTrajType>>& aPredictorPromise,
      std::shared_ptr<Sat3DSystemType> aSatSys,
      const ctrl::satellite_predictor_options& aSatOptions,
      std::shared_ptr<recorder::data_extractor> aDataIn,
      sat3D_meas_est_pred_to_recorders aDataLogger,
      std::atomic<double>* aCurrentTargetAnimTime)
      : predictor_promise(),
        satellite3D_system(aSatSys),
        sat_options(aSatOptions),
        data_in(aDataIn),
        data_logger(aDataLogger),
        current_target_anim_time(aCurrentTargetAnimTime) {
    predictor_promise.swap(aPredictorPromise);
  }

  int operator()() {
    using namespace ctrl;
    using namespace pp;

    bool promise_fulfilled = false;

    try {

      std::shared_ptr<TempSpaceType> sat_temp_space =
          satellite3D_system->get_temporal_state_space(
              0.0, sat_options.predict_time_horizon);
      std::shared_ptr<TempBeliefSpaceType> sat_temp_belief_space =
          satellite3D_system->get_temporal_belief_space(
              0.0, sat_options.predict_time_horizon);

      InputBeliefType b_u;
      b_u.set_mean_state(vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
      b_u.set_covariance(
          CovarType(CovarMatType(sat_options.input_disturbance)));

      OutputBeliefType b_z;
      b_z.set_mean_state(vect_n<double>(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0));
      b_z.set_covariance(
          CovarType(CovarMatType(sat_options.measurement_noise)));

      StateBeliefType b = satellite3D_system->get_zero_state_belief(10.0);
      if ((b.get_covariance().get_matrix().get_row_count() > 12) &&
          (sat_options.system_kind &
           ctrl::satellite_predictor_options::TSOSAKF) &&
          (sat_options.steady_param_covariance.get_row_count() + 12 ==
           b.get_covariance().get_matrix().get_row_count())) {
        mat<double, mat_structure::square> P(b.get_covariance().get_matrix());
        set_block(P, sat_options.steady_param_covariance, 12, 12);
        b.set_covariance(CovarType(CovarMatType(P)));
      }

      double last_time = 0.0;
      (*current_target_anim_time) = last_time;
      double current_Pnorm =
          norm_2(b.get_covariance().get_matrix()(range(0, 12), range(0, 12)));
      recorder::named_value_row nvr_in = data_in->getFreshNamedValueRow();
      double init_time = -1000.0;

      while ((current_Pnorm > sat_options.predict_Pnorm_threshold) &&
             !prediction_should_stop) {

        (*data_in) >> nvr_in;

        vect_n<double> z(nvr_in["p_x"], nvr_in["p_y"], nvr_in["p_z"],
                         nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"],
                         nvr_in["q_3"]);

        try {

          vect<double, 3> w(nvr_in["w_x"], nvr_in["w_y"], nvr_in["w_z"]);
          z.resize(10);
          z[7] = w[0];
          z[8] = w[1];
          z[9] = w[2];

          vect<double, 3> a(nvr_in["acc_x"], nvr_in["acc_y"], nvr_in["acc_z"]);
          vect<double, 3> m(nvr_in["mag_x"], nvr_in["mag_y"], nvr_in["mag_z"]);
          z.resize(16);
          z[10] = a[0];
          a[11] = a[1];
          z[12] = a[2];
          z[13] = m[0];
          a[14] = m[1];
          z[15] = m[2];

        } catch (recorder::out_of_bounds& e) {
          RK_UNUSED(e);
        }

        b_z.set_mean_state(z);
        b_u.set_mean_state(vect_n<double>(nvr_in["f_x"], nvr_in["f_y"],
                                          nvr_in["f_z"], nvr_in["t_x"],
                                          nvr_in["t_y"], nvr_in["t_z"]));

        tsos_aug_inv_kalman_filter_step(*satellite3D_system,
                                        sat_temp_space->get_space_topology(), b,
                                        b_u, b_z, last_time);

        current_Pnorm =
            norm_2(b.get_covariance().get_matrix()(range(0, 12), range(0, 12)));

        last_time = nvr_in["time"];
        (*current_target_anim_time) = last_time;
        if (init_time < -900.0) {
          init_time = last_time;
        }

        data_logger.add_record(b, b, b_u, b_z, last_time);
      }

      typename BeliefPredTrajType::assumption pred_assumpt =
          BeliefPredTrajType::no_measurements;
      switch (sat_options.predict_assumption) {
        case satellite_predictor_options::
            no_measurements:  // No future measurements
          pred_assumpt = BeliefPredTrajType::no_measurements;
          break;
        case satellite_predictor_options::
            most_likely_measurements:  // Maximum-likelihood Measurements
          pred_assumpt = BeliefPredTrajType::most_likely_measurements;
          break;
        case satellite_predictor_options::full_certainty:  // Full certainty
          // FIXME: make this into what it really should be (what is that? ... no sure)
          pred_assumpt = BeliefPredTrajType::most_likely_measurements;
          break;
      }

      std::shared_ptr<BeliefPredTrajType> predictor(new BeliefPredTrajType(
          sat_temp_belief_space, TempBeliefPointType(last_time, b),
          InputTrajType(vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
          PredFactoryType(satellite3D_system,
                          CovarMatType(sat_options.input_disturbance),
                          CovarMatType(sat_options.measurement_noise)),
          pred_assumpt));

      predictor_promise.set_value(predictor);
      promise_fulfilled = true;

      TempBeliefPointType b_pred;

      while (!prediction_should_stop) {

        (*data_in) >> nvr_in;

        vect_n<double> z(nvr_in["p_x"], nvr_in["p_y"], nvr_in["p_z"],
                         nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"],
                         nvr_in["q_3"]);

        try {

          vect<double, 3> w(nvr_in["w_x"], nvr_in["w_y"], nvr_in["w_z"]);
          z.resize(10);
          z[7] = w[0];
          z[8] = w[1];
          z[9] = w[2];

          vect<double, 3> a(nvr_in["acc_x"], nvr_in["acc_y"], nvr_in["acc_z"]);
          vect<double, 3> m(nvr_in["mag_x"], nvr_in["mag_y"], nvr_in["mag_z"]);
          z.resize(16);
          z[10] = a[0];
          a[11] = a[1];
          z[12] = a[2];
          z[13] = m[0];
          a[14] = m[1];
          z[15] = m[2];

        } catch (ReaK::recorder::out_of_bounds& e) {
          RK_UNUSED(e);
        }

        b_z.set_mean_state(z);
        b_u.set_mean_state(vect_n<double>(nvr_in["f_x"], nvr_in["f_y"],
                                          nvr_in["f_z"], nvr_in["t_x"],
                                          nvr_in["t_y"], nvr_in["t_z"]));

        tsos_aug_inv_kalman_filter_step(*satellite3D_system,
                                        sat_temp_space->get_space_topology(), b,
                                        b_u, b_z, last_time);

        last_time = nvr_in["time"];
        (*current_target_anim_time) = last_time;

        b_pred = predictor->get_point_at_time(last_time);

        data_logger.add_record(b, b_pred.pt, b_u, b_z, last_time);
      }
    } catch (std::exception& e) {
      if (!promise_fulfilled) {
        predictor_promise.set_exception(ReaKaux::make_exception_ptr(
            std::runtime_error("Could not fulfill the promise of a predictor "
                               "due to some error during estimation!")));
      }
      return 0;
    }
    return 0;
  }
};
}  // namespace

template <typename Sat3DSystemType, typename TempSpaceType,
          typename BeliefPredTrajType>
std::enable_if_t<ctrl::is_augmented_ss_system_v<Sat3DSystemType>,
                 std::shared_ptr<TargetPredConfigWidget::trajectory_type>>
construct_wrapped_trajectory(
    std::shared_ptr<TempSpaceType> sat_temp_space,
    const std::shared_ptr<BeliefPredTrajType>& predictor,
    const ctrl::satellite_predictor_options& sat_options,
    double extended_time_horizon) {
  using namespace ctrl;
  using namespace pp;

  using BaseSpaceType = se3_1st_order_topology_t<double>;
  using TemporalBaseSpaceType =
      temporal_space<BaseSpaceType, time_poisson_topology, time_distance_only>;

#define RK_D_INF std::numeric_limits<double>::infinity()
  auto sat_base_temp_space = std::make_shared<TemporalBaseSpaceType>(
      "satellite3D_temporal_space",
      make_se3_space("satellite3D_state_space",
                     vect<double, 3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
                     vect<double, 3>(RK_D_INF, RK_D_INF, RK_D_INF), RK_D_INF,
                     RK_D_INF),
      time_poisson_topology("satellite3D_time_space", sat_options.time_step,
                            sat_options.predict_time_horizon * 0.5));
#undef RK_D_INF

  using MLTrajType = transformed_trajectory<TempSpaceType, BeliefPredTrajType,
                                            maximum_likelihood_map>;
  using BaseTrajType = transformed_trajectory<TemporalBaseSpaceType, MLTrajType,
                                              augmented_to_state_map>;
  using StateTrajType = trajectory_base<TemporalBaseSpaceType>;
  using WrappedStateTrajType = trajectory_wrapper<BaseTrajType>;

  predictor->set_minimal_horizon(extended_time_horizon);

  auto ML_traj = std::make_shared<MLTrajType>(sat_temp_space, predictor);

  return std::make_shared<WrappedStateTrajType>(
      "sat3D_predicted_traj", BaseTrajType(sat_base_temp_space, ML_traj));
}

template <typename Sat3DSystemType, typename TempSpaceType,
          typename BeliefPredTrajType>
std::enable_if_t<!ctrl::is_augmented_ss_system_v<Sat3DSystemType>,
                 std::shared_ptr<TargetPredConfigWidget::trajectory_type>>
construct_wrapped_trajectory(
    std::shared_ptr<TempSpaceType> sat_temp_space,
    const std::shared_ptr<BeliefPredTrajType>& predictor,
    const ctrl::satellite_predictor_options& sat_options,
    double extended_time_horizon) {
  using namespace pp;
  using namespace ctrl;

  using MLTrajType = transformed_trajectory<TempSpaceType, BeliefPredTrajType,
                                            maximum_likelihood_map>;
  using WrappedStateTrajType = trajectory_wrapper<MLTrajType>;

  predictor->set_minimal_horizon(extended_time_horizon);
  auto ML_traj = std::make_shared<MLTrajType>(sat_temp_space, predictor);

  return std::make_shared<WrappedStateTrajType>("sat3D_predicted_traj",
                                                *ML_traj);
}

template <typename Sat3DSystemType>
std::shared_ptr<TargetPredConfigWidget::trajectory_type>
start_state_predictions(std::shared_ptr<Sat3DSystemType> satellite3D_system,
                        const ctrl::satellite_predictor_options& sat_options,
                        std::shared_ptr<recorder::data_extractor> data_in,
                        sat3D_meas_est_pred_to_recorders data_logger,
                        std::atomic<double>* current_target_anim_time) {
  using namespace ctrl;
  using namespace pp;

  using BeliefSpaceType = typename Sat3DSystemType::belief_space_type;
  using InputTrajType = constant_trajectory<vector_topology<vect_n<double>>>;
  using PredFactoryType =
      typename try_TSOSAIKF_belief_transfer_factory<Sat3DSystemType>::type;
  using BeliefPredTrajType =
      belief_predicted_trajectory<BeliefSpaceType, PredFactoryType,
                                  InputTrajType>;

  prediction_should_stop = true;
  if (prediction_executer.joinable()) {
    prediction_executer.join();
  }

  std::promise<std::shared_ptr<BeliefPredTrajType>> predictor_promise;
  auto predictor_future = predictor_promise.get_future();

  prediction_should_stop = false;
  prediction_executer = std::thread(prediction_updater<Sat3DSystemType>(
      predictor_promise, satellite3D_system, sat_options, data_in, data_logger,
      current_target_anim_time));

  try {
    return construct_wrapped_trajectory<Sat3DSystemType>(
        satellite3D_system->get_temporal_state_space(
            0.0, sat_options.predict_time_horizon),
        predictor_future.get(), sat_options,
        sat_options.predict_time_horizon + (*current_target_anim_time));
  } catch (std::runtime_error& e) {
    RK_UNUSED(e);
    return std::shared_ptr<TargetPredConfigWidget::trajectory_type>();
  }
}

void TargetPredConfigWidget::startStatePrediction() {
  using namespace ctrl;

  // ---------- Create input data stream ----------

  recorder::data_stream_options data_in_opt;

  if (useTCP()) {
    data_in_opt.kind = recorder::data_stream_options::tcp_stream;
  } else if (useUDP()) {
    data_in_opt.kind = recorder::data_stream_options::udp_stream;
  } else {
    data_in_opt.kind = recorder::data_stream_options::raw_udp_stream;
  }

  {
    std::stringstream ss;
    ss << getServerAddress() << ":" << getPortNumber();
    data_in_opt.file_name = ss.str();
  }

  data_in_opt.time_sync_name = "time";

  sat_options.imbue_names_for_received_meas(data_in_opt);

  std::vector<std::string> names_in;
  std::shared_ptr<recorder::data_extractor> data_in;
  std::tie(data_in, names_in) = data_in_opt.create_extractor();

  sat_options.imbue_names_for_received_meas(meas_out_opt);
  sat_options.imbue_names_for_state_estimates(est_out_opt);
  sat_options.imbue_names_for_state_estimates(pred_out_opt);

  sat3D_meas_est_pred_to_recorders data_logger(meas_out_opt.create_recorder(),
                                               est_out_opt.create_recorder(),
                                               pred_out_opt.create_recorder());

  // ---------- Call impl function with the appropriate system ----------

  using sat_opt = satellite_predictor_options;

  if (sat_options.system_kind & sat_opt::gyro_measures) {
    switch (sat_options.system_kind & 15) {
      case sat_opt::invariant:   // IEKF
      case sat_opt::invar_mom:   // IMKF
      case sat_opt::invar_mom2:  // IMKFv2
        *target_anim_traj = start_state_predictions(
            sat_options.get_gyro_sat_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_em:  // IMKF_em
        *target_anim_traj = start_state_predictions(
            sat_options.get_gyro_em_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_emd:  // IMKF_emd
        *target_anim_traj = start_state_predictions(
            sat_options.get_gyro_emd_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_emdJ:  // IMKF_emdJ
        *target_anim_traj = start_state_predictions(
            sat_options.get_gyro_emdJ_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
    }
  } else {
    switch (sat_options.system_kind & 15) {
      case sat_opt::invariant:   // IEKF
      case sat_opt::invar_mom:   // IMKF
      case sat_opt::invar_mom2:  // IMKFv2
        *target_anim_traj = start_state_predictions(
            sat_options.get_base_sat_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_em:  // IMKF_em
        *target_anim_traj = start_state_predictions(
            sat_options.get_em_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_emd:  // IMKF_emd
        *target_anim_traj = start_state_predictions(
            sat_options.get_emd_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
      case sat_opt::invar_mom_emdJ:  // IMKF_emdJ
        *target_anim_traj = start_state_predictions(
            sat_options.get_emdJ_airship_system(), sat_options, data_in,
            data_logger, current_target_anim_time);
        break;
    }
  }
}

void TargetPredConfigWidget::stopStatePrediction() {
  prediction_should_stop = true;
  if (prediction_executer.joinable())
    prediction_executer.join();
}

}  // namespace ReaK::qt
