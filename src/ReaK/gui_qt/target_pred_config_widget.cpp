
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

#include "target_pred_config_widget.hpp"

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QMainWindow>
#include <QProcess>

#include "serialization/archiver_factory.hpp"

#include "recorders/tcp_recorder.hpp"
#include "recorders/udp_recorder.hpp"

namespace ReaK {
  
namespace rkqt {


namespace detail {
  
  struct inertia_tensor_storage_impl : public shared_object {
    
    mat<double,mat_structure::symmetric> inertia_tensor;
    
    inertia_tensor_storage_impl() : inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0) { };
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(inertia_tensor);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(inertia_tensor);
    };
    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_tensor_storage_impl,0xBEEF0001,1,"inertia_tensor_storage_impl",shared_object)
  };
    
  struct IMU_config_storage_impl : public shared_object {
    
    unit_quat<double> IMU_orientation;
    vect<double,3> IMU_location;
    unit_quat<double> earth_orientation;
    vect<double,3> mag_field_direction;
    
    IMU_config_storage_impl() : IMU_orientation(), IMU_location(),
                                earth_orientation(), mag_field_direction(1.0,0.0,0.0) { };
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
        & RK_SERIAL_SAVE_WITH_NAME(IMU_location)
        & RK_SERIAL_SAVE_WITH_NAME(earth_orientation)
        & RK_SERIAL_SAVE_WITH_NAME(mag_field_direction);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
        & RK_SERIAL_LOAD_WITH_NAME(IMU_location)
        & RK_SERIAL_LOAD_WITH_NAME(earth_orientation)
        & RK_SERIAL_LOAD_WITH_NAME(mag_field_direction);
    };
    RK_RTTI_MAKE_CONCRETE_1BASE(IMU_config_storage_impl,0xBEEF0002,1,"IMU_config_storage_impl",shared_object)
  };
  
};



static QString last_used_path;





TargetPredConfigWidget::TargetPredConfigWidget(QWidget * parent, Qt::WindowFlags flags) :
  QDockWidget(tr("Predictor"), parent, flags),
  Ui::TargetPredConfig(),
  inertia_storage(new detail::inertia_tensor_storage_impl()),
  IMU_storage(new detail::IMU_config_storage_impl()),
  objtree_sch_bld(),
  ot_inertia_graph(shared_ptr< serialization::object_graph >(new serialization::object_graph())),
  ot_inertia_root(add_vertex(*ot_inertia_graph)),
  ot_inertia_win(this, Qt::WindowFlags(Qt::Popup | Qt::Dialog)),
  ot_IMU_graph(shared_ptr< serialization::object_graph >(new serialization::object_graph())),
  ot_IMU_root(add_vertex(*ot_IMU_graph)),
  ot_IMU_win(this, Qt::WindowFlags(Qt::Popup | Qt::Dialog))
{
  this->QDockWidget::setWidget(new QWidget(this));
  setupUi(this->QDockWidget::widget());
  
  connect(this->kf_model_selection, SIGNAL(currentIndexChanged(int)), this, SLOT(onUpdateAvailableOptions(int)));
  connect(this->actionValuesChanged, SIGNAL(triggered()), this, SLOT(onConfigsChanged()));
  connect(this->load_button, SIGNAL(clicked()), this, SLOT(loadPredictorConfig()));
  connect(this->save_button, SIGNAL(clicked()), this, SLOT(savePredictorConfig()));
  
  connect(this->I_load_button, SIGNAL(clicked()), this, SLOT(loadInertiaTensor()));
  connect(this->I_edit_button, SIGNAL(clicked()), this, SLOT(editInertiaTensor()));
  connect(this->I_save_button, SIGNAL(clicked()), this, SLOT(saveInertiaTensor()));
  
  connect(this->IMU_load_button, SIGNAL(clicked()), this, SLOT(loadIMUConfig()));
  connect(this->IMU_edit_button, SIGNAL(clicked()), this, SLOT(editIMUConfig()));
  connect(this->IMU_save_button, SIGNAL(clicked()), this, SLOT(saveIMUConfig()));
  
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
  };
  
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
  };
  
  
  updateConfigs();
  
};

TargetPredConfigWidget::~TargetPredConfigWidget() {
  delete this->QDockWidget::widget();
};


double TargetPredConfigWidget::getTimeStep() const {
  return this->time_step_spin->value();
};


double TargetPredConfigWidget::getMass() const {
  return this->mass_spin->value();
};


const mat<double,mat_structure::symmetric>& TargetPredConfigWidget::getInertiaTensor() const { return inertia_storage->inertia_tensor; };
    

mat<double,mat_structure::diagonal> TargetPredConfigWidget::getInputDisturbance() const {
  mat<double,mat_structure::diagonal> result(6,true);
  result(2,2) = result(1,1) = result(0,0) = this->Qf_spin->value();
  result(5,5) = result(4,4) = result(3,3) = this->Qt_spin->value();
  return result;
};

mat<double,mat_structure::diagonal> TargetPredConfigWidget::getMeasurementNoise() const {
  
  std::size_t m_noise_size = 6;
  if( this->gyro_check->isChecked() )
    m_noise_size = 9;
  if( this->IMU_check->isChecked() )
    m_noise_size = 15;
  
  mat<double,mat_structure::diagonal> result(m_noise_size,true);
  result(2,2) = result(1,1) = result(0,0) = this->Rpos_spin->value();
  result(5,5) = result(4,4) = result(3,3) = this->Rang_spin->value();
  if( this->gyro_check->isChecked() ) {
    result(8,8) = result(7,7) = result(6,6) = this->Rgyro_spin->value();
  };
  if( this->IMU_check->isChecked() ) {
    result(11,11) = result(10,10) = result(9,9) = this->Racc_spin->value();
    result(14,14) = result(13,13) = result(12,12) = this->Rmag_spin->value();
  };
  
  return result;
};


const unit_quat<double>& TargetPredConfigWidget::getIMUOrientation() const { 
  return IMU_storage->IMU_orientation;
};

const vect<double,3>& TargetPredConfigWidget::getIMULocation() const { 
  return IMU_storage->IMU_location;
};

const unit_quat<double>& TargetPredConfigWidget::getEarthOrientation() const { 
  return IMU_storage->earth_orientation;
};

const vect<double,3>& TargetPredConfigWidget::getMagFieldDirection() const { 
  return IMU_storage->mag_field_direction;
};


double TargetPredConfigWidget::getTimeHorizon() const {
  return this->horizon_spin->value();
};

double TargetPredConfigWidget::getPThreshold() const {
  return this->Pthreshold_spin->value();
};


std::string TargetPredConfigWidget::getServerAddress() const {
  return this->ip_addr_edit->text().toStdString();
};

int TargetPredConfigWidget::getPortNumber() const {
  return this->port_spin->value();
};

bool TargetPredConfigWidget::useUDP() const {
  return this->udp_radio->isChecked();
};

bool TargetPredConfigWidget::useTCP() const {
  return this->tcp_radio->isChecked();
};

std::string TargetPredConfigWidget::getStartScript() const {
  return this->start_script_edit->text().toStdString();
};


void TargetPredConfigWidget::onConfigsChanged() {
  
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


void TargetPredConfigWidget::updateConfigs() {
  
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



void TargetPredConfigWidget::onUpdateAvailableOptions(int filter_method) {
  
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


void TargetPredConfigWidget::savePredictorConfig() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Predictor Configurations..."), last_used_path,
    tr("Target Predictor Configurations (*.tpred.rkx *.tpred.rkb *.tpred.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    sat_options.save_all_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void TargetPredConfigWidget::loadPredictorConfig() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Predictor Configurations..."), last_used_path,
    tr("Target Predictor Configurations (*.tpred.rkx *.tpred.rkb *.tpred.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    sat_options.load_all_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
};


ctrl::satellite_predictor_options TargetPredConfigWidget::getSatPredictorOptions() const { return sat_options; };


    
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


void TargetPredConfigWidget::saveIMUConfig() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save IMU Configurations..."), last_used_path,
    tr("Target IMU Configurations (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    sat_options.save_IMU_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void TargetPredConfigWidget::editIMUConfig() {
  ot_IMU_win.show();
};

void TargetPredConfigWidget::loadIMUConfig() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open IMU Configurations..."), last_used_path,
    tr("Target IMU Configurations (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    sat_options.load_IMU_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
};


namespace {

struct prediction_updater {
  
  typedef satellite_predict_data::system_type SysType;
  typedef satellite_predict_data::temp_state_space_type TempStateSpaceType;
  typedef satellite_predict_data::belief_space_type BeliefSpaceType;
  typedef pp::topology_traits< BeliefSpaceType >::point_type BeliefPointType;
  typedef satellite_predict_data::belief_pred_traj_type BeliefPredTrajType;
  typedef BeliefPredTrajType::topology TempBeliefSpaceType;
  typedef pp::topology_traits< TempBeliefSpaceType >::point_type TempBeliefPointType;
  
  typedef satellite_predict_data::input_type InputType;
  typedef ctrl::discrete_sss_traits< SysType >::output_type OutputType;
  typedef ctrl::covariance_matrix< vect_n<double> > IOCovType;
  typedef ctrl::gaussian_belief_state< InputType,  IOCovType > InputBeliefType;
  typedef ctrl::gaussian_belief_state< OutputType, IOCovType > OutputBeliefType;
  
  shared_ptr< BeliefPredTrajType > predictor;
  
  shared_ptr< SysType > satellite3D_system;
  shared_ptr< TempStateSpaceType > sat_temp_space;
  
  recorder::named_value_row nvr_in;
  shared_ptr< recorder::data_extractor > data_in;
  
  BeliefPointType b;
  InputBeliefType b_u;
  OutputBeliefType b_z;
  
  double last_time;
  double diff_tolerance;
  
  static volatile bool should_stop;
  
  prediction_updater(
    shared_ptr< BeliefPredTrajType > aPredictor,
    shared_ptr< SysType > aSatSys,
    shared_ptr< TempStateSpaceType > aSatTempSpace,
    recorder::named_value_row aNVRIn,
    shared_ptr< recorder::data_extractor > aDataIn,
    BeliefPointType aB,
    InputBeliefType aBU,
    OutputBeliefType aBZ,
    double aLastTime, 
    double aDiffTolerance
  ) : 
    predictor(aPredictor),
    satellite3D_system(aSatSys),
    sat_temp_space(aSatTempSpace), 
    nvr_in(aNVRIn), data_in(aDataIn),
    b(aB), b_u(aBU), b_z(aBZ),
    last_time(aLastTime), diff_tolerance(aDiffTolerance) { };
  
  int operator()() {
    
    TempBeliefPointType b_pred;
    
    try {
      while ( !should_stop ) {
        
        (*data_in) >> nvr_in;
        
        b_z.set_mean_state(vect_n<double>(
          nvr_in["p_x"], nvr_in["p_y"], nvr_in["p_z"],
          nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"]));
        
        invariant_kalman_filter_step(*satellite3D_system, 
                                     sat_temp_space->get_space_topology(), 
                                     b, b_u, b_z, last_time);
        
        last_time = nvr_in["time"];
        
        b_pred = predictor->get_point_at_time(last_time);
        
        if( predictor->get_temporal_space().get_space_topology().distance(b, b_pred.pt) > diff_tolerance )
          predictor->set_initial_point(TempBeliefPointType(last_time, b));
        
      };
    } catch (std::exception& e) {
      return 0;
    };
    return 0;
  };
  
  static shared_ptr< ReaKaux::thread > executer;
  
  static QProcess side_script;
  
};

volatile bool prediction_updater::should_stop = false;
shared_ptr< ReaKaux::thread > prediction_updater::executer = shared_ptr< ReaKaux::thread >();
QProcess prediction_updater::side_script;

};


void TargetPredConfigWidget::startStatePrediction() {
  
  std::string start_script = getStartScript();
  if( start_script != "" ) {
    prediction_updater::side_script.kill();
    prediction_updater::side_script.start(QString::fromStdString(start_script));
  };
  
  using namespace ctrl;
  
  typedef satellite_predict_data::system_type SysType;
  
  switch(this->kf_model_selection->currentIndex()) {
    case 0:  // IEKF
      satellite3D_system = shared_ptr< SysType >(
        new satellite3D_inv_dt_system("satellite3D_inv", getMass(), getInertiaTensor(), getTimeStep())
      );
      break;
    case 1:  // IMKF
      satellite3D_system = shared_ptr< SysType >(
        new satellite3D_imdt_sys("satellite3D_invmom", getMass(), getInertiaTensor(), getTimeStep())
      );
      break;
    case 2:  // IMKFv2
      satellite3D_system = shared_ptr< SysType >(
        new satellite3D_imdt_sys("satellite3D_invmom", getMass(), getInertiaTensor(), getTimeStep(), 2)
      );
      break;
  };
  
  
  typedef satellite_predict_data::state_space_type StateSpaceType;
  typedef pp::topology_traits< StateSpaceType >::point_type StatePointType;
  
  typedef satellite_predict_data::temp_state_space_type TempStateSpaceType;
  
  shared_ptr< TempStateSpaceType > sat_temp_space(new TempStateSpaceType(
    "satellite3D_temporal_space",
    StateSpaceType(),
    pp::time_poisson_topology("satellite3D_time_space", getTimeStep(), getTimeHorizon())));
  
  
  typedef satellite_predict_data::input_type InputType;
  typedef satellite_predict_data::input_traj_type InputTrajType;
  
  InputTrajType in_cst_traj( InputType(0.0,0.0,0.0,0.0,0.0,0.0) );
  
  
  typedef satellite_predict_data::pred_factory_type PredFactoryType;
  typedef PredFactoryType::matrix_type NoiseMatType;
  
  
  typedef satellite_predict_data::covar_type StateCovarType;
  typedef satellite_predict_data::covar_space_type CovarSpaceType;
  typedef satellite_predict_data::belief_space_type BeliefSpaceType;
  typedef pp::topology_traits< BeliefSpaceType >::point_type BeliefPointType;
  
  
  StatePointType x_init;
  set_position(x_init, vect<double,3>(0.0, 0.0, 0.0));
  set_velocity(x_init, vect<double,3>(0.0, 0.0, 0.0));
  set_quaternion(x_init, unit_quat<double>());
  set_ang_velocity(x_init, vect<double,3>(0.0, 0.0, 0.0));
  
  BeliefPointType b_init(x_init, StateCovarType(StateCovarType::matrix_type(mat<double,mat_structure::diagonal>(12,100.0))));
  
  
  typedef satellite_predict_data::belief_pred_traj_type BeliefPredTrajType;
  
  BeliefPredTrajType::assumption pred_assumpt = BeliefPredTrajType::no_measurements;
  switch(this->predict_assumption_selection->currentIndex()) {
    case 0:  // No future measurements
      pred_assumpt = BeliefPredTrajType::no_measurements;
      break;
    case 1:  // Maximum-likelihood Measurements
      pred_assumpt = BeliefPredTrajType::most_likely_measurements;
      break;
    case 2:  // Full certainty
      // FIXME: make this into what it really should be (what is that? ... no sure)
      pred_assumpt = BeliefPredTrajType::most_likely_measurements;
      break;
  };
  
//   typedef pp::temporal_space<BeliefTopology, pp::time_poisson_topology, pp::time_distance_only> topology;
  typedef BeliefPredTrajType::topology TempBeliefSpaceType;
  typedef pp::topology_traits< TempBeliefSpaceType >::point_type TempBeliefPointType;
  
  shared_ptr< TempBeliefSpaceType > sat_temp_belief_space(new TempBeliefSpaceType(
    "satellite3D_temporal_belief_space",
    BeliefSpaceType(shared_ptr< StateSpaceType > (new StateSpaceType()), 
                    shared_ptr< CovarSpaceType >(new CovarSpaceType(12)), "satellite3D_belief_space"),
    pp::time_poisson_topology("satellite3D_time_space", getTimeStep(), getTimeHorizon())));
  
  
  pred_anim_data.predictor = shared_ptr< BeliefPredTrajType >(new BeliefPredTrajType(
    sat_temp_belief_space, 
    TempBeliefPointType(0.0, b_init), 
    in_cst_traj,
    PredFactoryType(satellite3D_system, NoiseMatType(getInputDisturbance()), NoiseMatType(getMeasurementNoise())),
    pred_assumpt
  ));
  
  
  typedef satellite_predict_data::ML_traj_type MLTrajType;
  
  pred_anim_data.trajectory = shared_ptr< MLTrajType >(new MLTrajType(sat_temp_space, pred_anim_data.predictor));
  
  typedef discrete_sss_traits< SysType >::output_type OutputType;
  typedef covariance_matrix< vect_n<double> > IOCovType;
  typedef IOCovType::matrix_type IOCovMatrixType;
  typedef gaussian_belief_state< InputType,  IOCovType > InputBeliefType;
  typedef gaussian_belief_state< OutputType, IOCovType > OutputBeliefType;
  
  // NOTE: at this point, I have constructed the trajectories needed for planning.
  //  Now, I need to start streaming measurements into it, and do so until I have good enough
  //  estimates to be able to construct a real predicted trajectory.
  
  
  shared_ptr< recorder::data_extractor > data_in;
  {
    std::stringstream ss;
    ss << getServerAddress() << ":" << getPortNumber();
    if( useTCP() )
      data_in = shared_ptr< recorder::data_extractor >(new recorder::tcp_extractor(ss.str()));
    else 
      data_in = shared_ptr< recorder::data_extractor >(new recorder::udp_extractor(ss.str()));
  };
  
  std::vector<std::string> names_in(data_in->getColCount(), "");
  for(std::size_t i = 0; i < names_in.size(); ++i)
    (*data_in) >> names_in[i];
  
  if( ( names_in.size() < 8 ) ||
      ( names_in[0] != "time" ) ||
      ( names_in[1] != "p_x" ) ||
      ( names_in[2] != "p_y" ) ||
      ( names_in[3] != "p_z" ) ||
      ( names_in[4] != "q_0" ) ||
      ( names_in[5] != "q_1" ) ||
      ( names_in[6] != "q_2" ) ||
      ( names_in[7] != "q_3" ) ) {
    QMessageBox::information(this,
                "Data stream seems to be corrupt!",
                "The meta-data obtained from the data-stream is corrupt or of the wrong format!",
                QMessageBox::Ok);
    return;
  };
  
  recorder::named_value_row nvr_in  = data_in->getFreshNamedValueRow();
  
  InputBeliefType b_u;
  OutputBeliefType b_z;
  BeliefPointType b = b_init;
  
  b_u.set_mean_state(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0));
  b_u.set_covariance(IOCovType(IOCovMatrixType(getInputDisturbance())));
  b_z.set_mean_state(vect_n<double>(0.0,0.0,0.0,1.0,0.0,0.0,0.0));
  b_z.set_covariance(IOCovType(IOCovMatrixType(getMeasurementNoise())));
  
  double last_time = 0.0;
  double current_Pnorm = norm_2(b.get_covariance().get_matrix());
  
  try {
    while ( current_Pnorm > getPThreshold() ) {
      
      (*data_in) >> nvr_in;
      
      b_z.set_mean_state(vect_n<double>(
        nvr_in["p_x"], nvr_in["p_y"], nvr_in["p_z"],
        nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"]));
      
      invariant_kalman_filter_step(*satellite3D_system, 
                                  sat_temp_space->get_space_topology(), 
                                  b, b_u, b_z, last_time);
      
//       const StatePointType& x_mean = b.get_mean_state();
//       result_vect[range(0,2)]   = get_position(x_mean);
//       result_vect[range(3,6)]   = get_quaternion(x_mean);
//       result_vect[range(7,9)]   = get_velocity(x_mean);
//       result_vect[range(10,12)] = get_ang_velocity(x_mean);
      
      current_Pnorm = norm_2(b.get_covariance().get_matrix());
      last_time = nvr_in["time"];
      
    };
  } catch (std::exception& e) {
    QMessageBox::information(this,
                "Error!",
                "An error occurred during the initial phase of measurement streaming!",
                QMessageBox::Ok);
    return;
  };
  
  pred_anim_data.predictor->set_initial_point(TempBeliefPointType(last_time, b));
  
  // Start a thread that will update the predictor as new data comes in.
  
  if( prediction_updater::executer ) {
    prediction_updater::should_stop = true;
    prediction_updater::executer->join();
    prediction_updater::executer.reset();
  };
  
  prediction_updater::should_stop = false;
  prediction_updater::executer = shared_ptr< ReaKaux::thread >(new ReaKaux::thread(
    prediction_updater(
    pred_anim_data.predictor,
    satellite3D_system,
    sat_temp_space, nvr_in, data_in,
    b, b_u, b_z, last_time,  getPThreshold())));
  
};


void TargetPredConfigWidget::stopStatePrediction() {
  
  if( prediction_updater::executer ) {
    prediction_updater::should_stop = true;
    prediction_updater::executer->join();
    prediction_updater::executer.reset();
  };
  
  prediction_updater::side_script.kill();
  
};


};

};














