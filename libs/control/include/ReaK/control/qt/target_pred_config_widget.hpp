/**
 * 
 * 
 * 
 */

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

#ifndef REAK_TARGET_PRED_CONFIG_WIDGET_HPP
#define REAK_TARGET_PRED_CONFIG_WIDGET_HPP

#include <ReaK/math/lin_alg/mat_alg_symmetric.hpp>
#include <ReaK/math/lin_alg/mat_alg_diagonal.hpp>

#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/core/base/atomic_incl.hpp>
#include <ReaK/core/recorders/data_record_options.hpp>

#include <ReaK/core/qt/rk_object_tree_widget.hpp>
#include <ReaK/core/qt/rk_prop_editor_widget.hpp>
#include <ReaK/core/serialization/scheme_builder.hpp>

#include <ReaK/control/systems/satellite_modeling_options.hpp>
#include <ReaK/topologies/interpolation/trajectory_base.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>
#include <ReaK/topologies/spaces/temporal_space.hpp>

#include <QDockWidget>
#include <QMainWindow>

namespace Ui {
  class TargetPredConfig;
};

namespace ReaK {
  
namespace qt {
  
/* forward-decl */
namespace detail {
  struct inertia_tensor_storage_impl;
  struct IMU_config_storage_impl;
};

class TargetPredConfigWidget : public QDockWidget {
    Q_OBJECT
  private:
    Ui::TargetPredConfig* ui;
  public:
    
    typedef pp::se3_1st_order_topology<double>::type sat_state_space_type;
    typedef pp::temporal_space<sat_state_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_space_type;
    typedef pp::trajectory_base< temporal_space_type > trajectory_type;
    
    TargetPredConfigWidget(shared_ptr< trajectory_type >* aTargetAnimTraj, 
                           ReaKaux::atomic<double>* aCurrentTargetAnimTime, 
                           QWidget * parent = NULL, Qt::WindowFlags flags = 0);
    virtual ~TargetPredConfigWidget();
    
  private slots:
    
    void onUpdateAvailableOptions(int);
    
    void onConfigsChanged();
    void updateConfigs();
    
    void savePredictorConfig();
    void loadPredictorConfig();
    
    void saveInertiaTensor();
    void editInertiaTensor();
    void loadInertiaTensor();
    
    void saveIMUConfig();
    void editIMUConfig();
    void loadIMUConfig();
    
  private:
    
    shared_ptr<detail::inertia_tensor_storage_impl> inertia_storage;
    shared_ptr<detail::IMU_config_storage_impl> IMU_storage;
    
    serialization::scheme_builder objtree_sch_bld;
    
    shared_ptr< serialization::object_graph > ot_inertia_graph;
    serialization::object_node_desc ot_inertia_root;
    ObjectTreeWidget* ot_inertia_widget;
    PropEditorWidget* ot_inertia_propedit;
    serialization::objtree_editor* ot_inertia_edit;
    QMainWindow ot_inertia_win;
    
    shared_ptr< serialization::object_graph > ot_IMU_graph;
    serialization::object_node_desc ot_IMU_root;
    ObjectTreeWidget* ot_IMU_widget;
    PropEditorWidget* ot_IMU_propedit;
    serialization::objtree_editor* ot_IMU_edit;
    QMainWindow ot_IMU_win;
    
    shared_ptr< ctrl::satellite3D_inv_dt_system > satellite3D_system;
    
    ctrl::satellite_predictor_options sat_options;
    
    shared_ptr< trajectory_type >* target_anim_traj;
    ReaKaux::atomic<double>* current_target_anim_time;
    
  public:
    
    recorder::data_stream_options meas_out_opt;
    recorder::data_stream_options est_out_opt;
    recorder::data_stream_options pred_out_opt;
    
    double getTimeStep() const;
    double getMass() const;
    const mat<double,mat_structure::symmetric>& getInertiaTensor() const;
    
    mat<double,mat_structure::diagonal> getInputDisturbance() const;
    mat<double,mat_structure::diagonal> getMeasurementNoise() const;
    
    const unit_quat<double>& getIMUOrientation() const;
    const vect<double,3>&    getIMULocation() const;
    const unit_quat<double>& getEarthOrientation() const;
    const vect<double,3>&    getMagFieldDirection() const;
    
    double getTimeHorizon() const;
    double getPThreshold() const;
    
    ctrl::satellite_predictor_options getSatPredictorOptions() const;
    
    std::string getServerAddress() const;
    int getPortNumber() const;
    bool useRawUDP() const;
    bool useUDP() const;
    bool useTCP() const;
    std::string getStartScript() const;
    
    void startStatePrediction();
    void stopStatePrediction();
    
    void savePredictorConfigurations(const std::string& aFilename);
    void loadPredictorConfigurations(const std::string& aFilename);
    
    
};

};

};

#endif














