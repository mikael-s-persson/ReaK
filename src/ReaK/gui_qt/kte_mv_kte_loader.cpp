
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

#include "kte_model_viewer.h"

#include <QMessageBox>


#include "shapes/oi_scene_graph.hpp"
#include "shapes/kte_chain_geometry.hpp"
#include "proximity/proxy_query_model.hpp"

#include "mbd_kte/kte_map_chain.hpp"
#include "kte_models/direct_kinematics_model.hpp"
#include "kte_models/inverse_kinematics_model.hpp"
#include "kte_models/inverse_dynamics_model.hpp"
#include "kte_models/manip_3R3R_arm.hpp"
#include "kte_models/manip_3R_arm.hpp"
#include "kte_models/manip_ERA_arm.hpp"
#include "kte_models/manip_P3R3R_arm.hpp"
#include "kte_models/manip_SCARA_arm.hpp"
#include "kte_models/manip_SSRMS_arm.hpp"

#include "topologies/joint_space_limits.hpp"

#include "serialization/archiver.hpp"

void KTEModelViewerEditor::loadFromArchive(ReaK::serialization::iarchive& in, QString fileContentExt) {
  
  if( (fileContentExt == tr("kte_ik")) || (fileContentExt == tr("kte_dk")) ) {
    ReaK::shared_ptr< ReaK::pose_3D<double> > base_frame;
    ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > kin_mdl;
    ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> > mdl_jt_lim;
    
    try {
      in >> base_frame;
      in >> kin_mdl;
      in >> mdl_jt_lim;
    } catch(std::exception& e) {
      QMessageBox::information(this,
                tr("Failed to load file!"),
                tr("Sorry, this file could not be loaded! Got error: '") + tr(e.what()) + tr("'"),
                QMessageBox::Ok);
      return;
    };
    
    if( !base_frame || !kin_mdl || !mdl_jt_lim ) {
      QString culprit;
      if( !base_frame )
        culprit += tr(" base-frame");
      if( !kin_mdl )
        culprit += tr(" kinematics-model");
      if( !mdl_jt_lim )
        culprit += tr(" joint-limits");
      QMessageBox::information(this,
                tr("Failed to load file!"),
                tr("Sorry, some required objects could not be loaded from this file! Culprits are ") + culprit,
                QMessageBox::Ok);
      return;
    };
    
    std::string kin_mdl_name = kin_mdl->getName();
    objtree_sch_bld << base_frame << kin_mdl << mdl_jt_lim;
    ReaK::serialization::object_node_desc base_frame_id = objtree_edit.add_new_object(base_frame);
    ReaK::serialization::object_node_desc kin_mdl_id    = objtree_edit.add_new_object(kin_mdl);
    ReaK::serialization::object_node_desc mdl_jt_lim_id = objtree_edit.add_new_object(mdl_jt_lim);
    objtree.mdl.refreshObjTree();
    
    mdl_base_frames[objtree_edit.get_object_name(base_frame_id)] = base_frame;
    dk_models[objtree_edit.get_object_name(kin_mdl_id)] = kin_mdl;
    mdl_jt_limits[objtree_edit.get_object_name(mdl_jt_lim_id)] = mdl_jt_lim;
    
    ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model > kin_ik_mdl = ReaK::rtti::rk_dynamic_ptr_cast< ReaK::kte::inverse_kinematics_model >(kin_mdl);
    if(kin_ik_mdl) 
      ik_models[objtree_edit.get_object_name(kin_mdl_id)] = kin_ik_mdl;
    
    
    
  } else if( fileContentExt == tr("kte_dyn") ) {
    
    
  } else if( fileContentExt == tr("geom") ) {
    
    
  } else if( fileContentExt == tr("model") ) {
    
    
  };
  
};





