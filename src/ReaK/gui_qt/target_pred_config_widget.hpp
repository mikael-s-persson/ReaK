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

#include "path_planning/path_planner_options.hpp"

#include "lin_alg/mat_alg_symmetric.hpp"
#include "lin_alg/mat_alg_diagonal.hpp"

#include "lin_alg/vect_alg.hpp"
#include "kinetostatics/quat_alg.hpp"

#include "ui_target_predictor_config.h"

#include "rk_object_tree_widget.hpp"
#include "rk_prop_editor_widget.hpp"
#include "serialization/scheme_builder.hpp"

#include <QDockWidget>

#include <QMainWindow>

namespace ReaK {
  
namespace rkqt {

class TargetPredConfigWidget : public QDockWidget, private Ui::TargetPredConfig {
    Q_OBJECT
  
  public:
    TargetPredConfigWidget(QWidget * parent = NULL, Qt::WindowFlags flags = 0);
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
    
    mat<double,mat_structure::symmetric> inertia_tensor;
    
    unit_quat<double> IMU_orientation;
    vect<double,3> IMU_location;
    unit_quat<double> earth_orientation;
    vect<double,3> mag_field_direction;
    
    
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
    
  public:
    
    double getTimeStep() const;
    double getMass() const;
    const mat<double,mat_structure::symmetric>& getInertiaTensor() const { return inertia_tensor; };
    
    mat<double,mat_structure::diagonal> getInputDisturbance() const;
    mat<double,mat_structure::diagonal> getMeasurementNoise() const;
    
    const unit_quat<double>& getIMUOrientation() const { return IMU_orientation; };
    const vect<double,3>&    getIMULocation() const { return IMU_location; };
    const unit_quat<double>& getEarthOrientation() const { return earth_orientation; };
    const vect<double,3>&    getMagFieldDirection() const { return mag_field_direction; };
    
    double getTimeHorizon() const;
    double getPThreshold() const;
    
    std::string getServerAddress() const;
    int getPortNumber() const;
    bool useUDP() const;
    bool useTCP() const;
    std::string getStartScript() const;
    
};

};

};

#endif














