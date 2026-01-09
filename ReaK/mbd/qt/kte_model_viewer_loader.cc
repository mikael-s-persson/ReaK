
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

#include "kte_model_viewer.hpp"

#include <QMessageBox>

#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/mbd/coin3D/oi_scene_graph.h"
#include "ReaK/mbd/kte/kte_chain_geometry.h"

#include "ReaK/mbd/kte/kte_map_chain.h"
#include "ReaK/mbd/models/direct_kinematics_model.h"
#include "ReaK/mbd/models/inverse_dynamics_model.h"
#include "ReaK/mbd/models/inverse_kinematics_model.h"
#include "ReaK/mbd/models/joint_space_limits.h"
#include "ReaK/mbd/models/manip_3R3R_arm.h"
#include "ReaK/mbd/models/manip_3R_arm.h"
#include "ReaK/mbd/models/manip_ERA_arm.h"
#include "ReaK/mbd/models/manip_P3R3R_arm.h"
#include "ReaK/mbd/models/manip_SCARA_arm.h"
#include "ReaK/mbd/models/manip_SSRMS_arm.h"

#include "ReaK/core/serialization/archiver.h"

namespace ReaK {

namespace qt {

std::string KTEModelViewerEditor::loadKinModelFromArchive(
    serialization::iarchive& in) {
  std::shared_ptr<pose_3D<double>> base_frame;
  std::shared_ptr<kte::direct_kinematics_model> kin_mdl;
  std::shared_ptr<kte::joint_limits_collection<double>> mdl_jt_lim;

  try {
    in >> base_frame;
    in >> kin_mdl;
    in >> mdl_jt_lim;
  } catch (std::exception& e) {
    QMessageBox::information(
        this, tr("Failed to load file!"),
        tr("Sorry, this file could not be loaded! Got error: '") +
            tr(e.what()) + tr("'"),
        QMessageBox::Ok);
    return std::string();
  };

  if (!base_frame || !kin_mdl || !mdl_jt_lim) {
    QString culprit;
    if (!base_frame)
      culprit += tr(" base-frame");
    if (!kin_mdl)
      culprit += tr(" kinematics-model");
    if (!mdl_jt_lim)
      culprit += tr(" joint-limits");
    QMessageBox::information(this, tr("Failed to load file!"),
                             tr("Sorry, some required objects could not be "
                                "loaded from this file! Culprits are ") +
                                 culprit,
                             QMessageBox::Ok);
    return std::string();
  };

  objtree_sch_bld << base_frame << kin_mdl << mdl_jt_lim;
  serialization::object_node_desc base_frame_id =
      objtree_edit.add_new_object(base_frame);
  serialization::object_node_desc kin_mdl_id =
      objtree_edit.add_new_object(kin_mdl);
  serialization::object_node_desc mdl_jt_lim_id =
      objtree_edit.add_new_object(mdl_jt_lim);
  objtree.mdl.refreshObjTree();

  dk_models[objtree_edit.get_object_name(kin_mdl_id)] = kin_mdl;
  mdl_to_base[objtree_edit.get_object_name(kin_mdl_id)] =
      objtree_edit.get_object_name(base_frame_id);
  mdl_base_frames[objtree_edit.get_object_name(base_frame_id)] = base_frame;
  mdl_to_jt_lim[objtree_edit.get_object_name(kin_mdl_id)] =
      objtree_edit.get_object_name(mdl_jt_lim_id);
  mdl_jt_limits[objtree_edit.get_object_name(mdl_jt_lim_id)] = mdl_jt_lim;

  std::shared_ptr<kte::inverse_kinematics_model> kin_ik_mdl =
      rtti::rk_dynamic_ptr_cast<kte::inverse_kinematics_model>(kin_mdl);
  if (kin_ik_mdl)
    ik_models[objtree_edit.get_object_name(kin_mdl_id)] = kin_ik_mdl;

  std::shared_ptr<const kte::kte_map_chain> kin_mdl_chain =
      kin_mdl->getKTEChain();
  if (kin_mdl_chain) {
    kin_mdl->doDirectMotion();
    std::shared_ptr<geom::oi_scene_graph> psg =
        view3d_menu.getGeometryGroup("Kinematics Models");
    (*psg) << (*kin_mdl_chain);
  };

  return objtree_edit.get_object_name(kin_mdl_id);
};

std::string KTEModelViewerEditor::addGeometricModel(
    const std::shared_ptr<geom::colored_model_3D>& mdl_geom) {
  if (!mdl_geom)
    return std::string();

  objtree_sch_bld << mdl_geom;
  serialization::object_node_desc mdl_geom_id =
      objtree_edit.add_new_object(mdl_geom);

  std::string sg_name = objtree_edit.get_object_name(mdl_geom_id);
  geom_models[sg_name] = mdl_geom;

  std::shared_ptr<geom::oi_scene_graph> psg =
      view3d_menu.getGeometryGroup("Geometric Models");
  (*psg) << (*mdl_geom);

  return sg_name;
};

std::string KTEModelViewerEditor::addProximityModel(
    const std::shared_ptr<geom::proxy_query_model_3D>& mdl_prox) {
  if (!mdl_prox)
    return std::string();

  objtree_sch_bld << mdl_prox;
  serialization::object_node_desc mdl_prox_id =
      objtree_edit.add_new_object(mdl_prox);

  std::string sg_name = objtree_edit.get_object_name(mdl_prox_id);
  proxy_models[sg_name] = mdl_prox;

  std::shared_ptr<geom::oi_scene_graph> psg =
      view3d_menu.getGeometryGroup("Proximity Models");
  (*psg) << (*mdl_prox);

  return sg_name;
};

void KTEModelViewerEditor::loadFromArchive(serialization::iarchive& in,
                                           QString fileContentExt) {

  if ((fileContentExt == tr("kte_ik")) || (fileContentExt == tr("kte_dk"))) {

    loadKinModelFromArchive(in);

  } else if (fileContentExt == tr("kte_dyn")) {

  } else if (fileContentExt == tr("geom")) {
    std::shared_ptr<kte::kte_chain_geometry_3D> kte_geom;

    try {
      in >> kte_geom;
    } catch (std::exception& e) {
      QMessageBox::information(
          this, tr("Failed to load file!"),
          tr("Sorry, this file could not be loaded! Got error: '") +
              tr(e.what()) + tr("'"),
          QMessageBox::Ok);
      return;
    };

    if (!kte_geom) {
      QMessageBox::information(
          this, tr("Failed to load file!"),
          tr("Sorry, KTE geometry could not be loaded from this file!"),
          QMessageBox::Ok);
      return;
    };

    std::string kin_mdl_name = kte_geom->get_name();
    objtree_sch_bld << kte_geom;
    serialization::object_node_desc kte_geom_id =
        objtree_edit.add_new_object(kte_geom);

    kte_geometries[objtree_edit.get_object_name(kte_geom_id)] = kte_geom;
    for (std::map<std::string,
                  std::shared_ptr<kte::direct_kinematics_model>>::iterator it =
             dk_models.begin();
         it != dk_models.end(); ++it) {

      std::shared_ptr<geom::colored_model_3D> mdl_geom;
      std::shared_ptr<geom::proxy_query_model_3D> mdl_prox;
      std::tie(mdl_geom, mdl_prox) =
          kte_geom->attachToKTEChain(*(it->second->getKTEChain()));

      addGeometricModel(mdl_geom);
      addProximityModel(mdl_prox);
    };
    objtree.mdl.refreshObjTree();

  } else if (fileContentExt == tr("model")) {

    std::string kte_mdl_str = loadKinModelFromArchive(in);

    std::shared_ptr<geom::colored_model_3D> mdl_geom;
    std::shared_ptr<geom::proxy_query_model_3D> mdl_prox;

    try {
      in >> mdl_geom;
      in >> mdl_prox;
    } catch (std::exception& e) {
      QMessageBox::information(
          this, tr("Failed to load file!"),
          tr("Sorry, this file could not be loaded! Got error: '") +
              tr(e.what()) + tr("'"),
          QMessageBox::Ok);
      return;
    };

    if (!mdl_geom || !mdl_prox) {
      QString culprit;
      if (!mdl_geom)
        culprit += tr(" geometry");
      if (!mdl_prox)
        culprit += tr(" proxy-model");
      QMessageBox::information(this, tr("Failed to load file!"),
                               tr("Sorry, some required objects could not be "
                                  "loaded from this file! Culprits are ") +
                                   culprit,
                               QMessageBox::Ok);
      return;
    };

    addGeometricModel(mdl_geom);
    addProximityModel(mdl_prox);
  };
};

void KTEModelViewerEditor::saveKTEChainGeometry(
    serialization::oarchive& out, const std::string&,
    const std::shared_ptr<kte::kte_chain_geometry_3D>& kte_geom) {
  out << kte_geom;
};

void KTEModelViewerEditor::saveGeometricModel(
    serialization::oarchive& out, const std::string&,
    const std::shared_ptr<geom::colored_model_3D>& geom_mdl) {
  out << geom_mdl;
};

void KTEModelViewerEditor::saveProximityModel(
    serialization::oarchive& out, const std::string&,
    const std::shared_ptr<geom::proxy_query_model_3D>& prox_mdl) {
  out << prox_mdl;
};

void KTEModelViewerEditor::saveDirectKinModel(
    serialization::oarchive& out, const std::string& mdl_name,
    const std::shared_ptr<kte::direct_kinematics_model>& dk_mdl) {
  std::shared_ptr<pose_3D<double>> base_frame =
      mdl_base_frames[mdl_to_base[mdl_name]];
  std::shared_ptr<kte::joint_limits_collection<double>> jt_limits =
      mdl_jt_limits[mdl_to_jt_lim[mdl_name]];

  out << base_frame << dk_mdl << jt_limits;
};

void KTEModelViewerEditor::saveInverseKinModel(
    serialization::oarchive& out, const std::string& mdl_name,
    const std::shared_ptr<kte::inverse_kinematics_model>& ik_mdl) {
  std::shared_ptr<pose_3D<double>> base_frame =
      mdl_base_frames[mdl_to_base[mdl_name]];
  std::shared_ptr<kte::joint_limits_collection<double>> jt_limits =
      mdl_jt_limits[mdl_to_jt_lim[mdl_name]];

  out << base_frame << ik_mdl << jt_limits;
};

void KTEModelViewerEditor::saveDynamicsModel(
    serialization::oarchive& out, const std::string& mdl_name,
    const std::shared_ptr<kte::inverse_dynamics_model>& dyn_mdl) {
  std::shared_ptr<pose_3D<double>> base_frame =
      mdl_base_frames[mdl_to_base[mdl_name]];
  std::shared_ptr<kte::joint_limits_collection<double>> jt_limits =
      mdl_jt_limits[mdl_to_jt_lim[mdl_name]];

  out << base_frame << dyn_mdl << jt_limits;
};
};  // namespace qt
};  // namespace ReaK
