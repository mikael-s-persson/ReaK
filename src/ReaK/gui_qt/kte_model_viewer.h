
/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_KTE_MODEL_VIEWER_H
#define REAK_KTE_MODEL_VIEWER_H


#include "ui_kte_model_viewer.h"
#include "rk_object_tree_widget.hpp"
#include "rk_prop_editor_widget.hpp"
#include "rk_view3d_menu.hpp"
#include "serialization/scheme_builder.hpp"

#include <map>
#include <string>

#include <QMainWindow>
#include <QDockWidget>

class SoSeparator;
class SoQtExaminerViewer;

namespace ReaK {
  namespace geom { 
    class kte_chain_geometry_3D;
    class oi_scene_graph;
    class colored_model_3D;
    class proxy_query_model_3D;
  };
  
  template <typename T>
  class pose_3D;
  
  namespace pp {
    template <typename T>
    struct joint_limits_collection;
  };
  
  namespace kte {
    class direct_kinematics_model;
    class inverse_kinematics_model;
    class inverse_dynamics_model;
  };
};


class KTEModelViewerEditor : public QMainWindow, private Ui::KTEModelView {
    Q_OBJECT
  
  public:
    KTEModelViewerEditor( QWidget * parent = 0, Qt::WindowFlags flags = 0 );
    ~KTEModelViewerEditor();
    
  private slots:
    
    void onLoad();
    void onSave();
    void onCloseAll();
    
    void onRefreshView();
    
  private:
    
    std::string loadKinModelFromArchive(ReaK::serialization::iarchive& in);
    std::string addGeometricModel(const ReaK::shared_ptr< ReaK::geom::colored_model_3D >& mdl_geom);
    std::string addProximityModel(const ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D >& mdl_prox);
    void loadFromArchive(ReaK::serialization::iarchive& in, QString fileContentExt);
    
    void saveKTEChainGeometry(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::geom::kte_chain_geometry_3D >& kte_geom);
    void saveGeometricModel(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::geom::colored_model_3D >& geom_mdl);
    void saveProximityModel(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D >& prox_mdl);
    void saveDirectKinModel(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::kte::direct_kinematics_model >& dk_mdl);
    void saveInverseKinModel(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model >& ik_mdl);
    void saveDynamicsModel(ReaK::serialization::oarchive& out, const std::string& mdl_name, const ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model >& dyn_mdl);
    
    ReaK::shared_ptr< ReaK::serialization::object_graph > objtree_graph;
    ReaK::serialization::object_node_desc objtree_root;
    
    ReaK::serialization::scheme_builder objtree_sch_bld;
    
    ReaK::rkqt::ObjectTreeWidget objtree;
    ReaK::rkqt::PropEditorWidget propedit;
    
    ReaK::serialization::objtree_editor& objtree_edit;
    
    ReaK::rkqt::View3DMenu view3d_menu;
    
    std::map< std::string, ReaK::shared_ptr< ReaK::geom::kte_chain_geometry_3D > > kte_geometries;
    std::map< std::string, ReaK::shared_ptr< ReaK::geom::colored_model_3D > > geom_models;
    std::map< std::string, ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > > proxy_models;
    
    std::map< std::string, ReaK::shared_ptr< ReaK::pose_3D<double> > > mdl_base_frames;
    std::map< std::string, ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> > > mdl_jt_limits;
    std::map< std::string, ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > > dk_models;
    std::map< std::string, ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model > > ik_models;
    std::map< std::string, ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model > > dyn_models;
    
    std::map< std::string, std::string > mdl_to_base;
    std::map< std::string, std::string > mdl_to_jt_lim;
    
};



#endif














