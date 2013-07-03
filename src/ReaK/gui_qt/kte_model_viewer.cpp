
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

#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>
#include <QMainWindow>
#include <QDir>



#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

#include "shapes/oi_scene_graph.hpp"
#include "shapes/kte_chain_geometry.hpp"
#include "proximity/proxy_query_model.hpp"

#include "mbd_kte/kte_map_chain.hpp"
#include "kte_models/direct_kinematics_model.hpp"
#include "kte_models/inverse_kinematics_model.hpp"
#include "kte_models/inverse_dynamics_model.hpp"
#include "kte_models/manip_P3R3R_arm.hpp"

#include "topologies/joint_space_limits.hpp"

#include "serialization/xml_archiver.hpp"
#include "serialization/bin_archiver.hpp"
#include "serialization/protobuf_archiver.hpp"

#include "optimization/optim_exceptions.hpp"



static QString last_used_path;



KTEModelViewerEditor::KTEModelViewerEditor(QWidget * parent, Qt::WindowFlags flags ) : 
                                           QMainWindow(parent,flags),
                                           objtree_graph(ReaK::shared_ptr< ReaK::serialization::object_graph >(new ReaK::serialization::object_graph())),
                                           objtree_root(add_vertex(*objtree_graph)),
                                           objtree_sch_bld(),
                                           objtree(objtree_graph, objtree_root),
                                           propedit(&(objtree.mdl)),
                                           objtree_edit(propedit.mdl.get_object_editor()),
                                           view3d_menu(this) {
  setupUi(this);
  using namespace ReaK;
  
  addDockWidget(Qt::LeftDockWidgetArea, &objtree);
  
  addDockWidget(Qt::LeftDockWidgetArea, &propedit);
  
  tabifyDockWidget(&objtree, &propedit);
  
  
  menubar->addMenu(&view3d_menu);
  
  
  /*
  connect(configs.actionStart_Robot, SIGNAL(triggered()), this, SLOT(startRobot()));
  connect(configs.actionExecutePlanner, SIGNAL(triggered()), this, SLOT(executePlanner()));
  connect(configs.actionJointChange, SIGNAL(triggered()), this, SLOT(onJointChange()));
  connect(configs.actionTargetChange, SIGNAL(triggered()), this, SLOT(onTargetChange()));
  
  connect(configs.actionRobotVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotVisible()));
  connect(configs.actionRobotKinVisibleToggle, SIGNAL(triggered()), this, SLOT(onRobotKinVisible()));
  connect(configs.actionTargetVisibleToggle, SIGNAL(triggered()), this, SLOT(onTargetVisible()));
  connect(configs.actionEnvVisibleToggle, SIGNAL(triggered()), this, SLOT(onEnvVisible()));
  connect(configs.actionProxyVisibleToggle, SIGNAL(triggered()), this, SLOT(onProxyVisible()));
  connect(configs.actionMGVisibleToggle, SIGNAL(triggered()), this, SLOT(onMGVisible()));
  connect(configs.actionSolutionsVisibleToggle, SIGNAL(triggered()), this, SLOT(onSolutionsVisible()));
  */
  
  connect(actionLoad, SIGNAL(triggered()), this, SLOT(onLoad()));
  connect(actionSave, SIGNAL(triggered()), this, SLOT(onSave()));
  connect(actionClose_all, SIGNAL(triggered()), this, SLOT(onCloseAll()));
  
  actionLoad->setIcon( style()->standardIcon( QStyle::SP_DialogOpenButton, NULL, actionLoad->menu() ) );
  actionSave->setIcon( style()->standardIcon( QStyle::SP_DialogSaveButton, NULL, actionSave->menu() ) );
  actionClose_all->setIcon( style()->standardIcon( QStyle::SP_DialogCloseButton, NULL, actionClose_all->menu() ) );
  actionQuit->setIcon( style()->standardIcon( QStyle::SP_DialogCloseButton, NULL, actionQuit->menu() ) );
  
  
  connect(&propedit.mdl, SIGNAL(sourceDataChanged()), this, SLOT(onRefreshView()));
  
  
  SoQt::init(this->centralwidget);
  
  sg_root = new SoSeparator;
  sg_root->ref();
  
  view3d_menu.setRoot(sg_root);
  
  // insert some code here.  Add geometries.
  
  eviewer = new SoQtExaminerViewer(this->centralwidget);
  eviewer->setSceneGraph(sg_root);
  eviewer->show();
  
#if 0
  
  r_info.sg_robot_geom->enableAnchorUpdates();
  r_info.sg_robot_kin->enableAnchorUpdates();
  r_info.sg_airship_geom->enableAnchorUpdates();
  
#endif
};


KTEModelViewerEditor::~KTEModelViewerEditor() {
  
  // insert some code here.  Delete geometries.
  
  delete eviewer;
  sg_root->unref();
  view3d_menu.setRoot(NULL);
  scene_graphs.clear();
  
  SoQt::done();
  
};



void KTEModelViewerEditor::onLoad() {
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Load a Model..."),
    last_used_path,
    tr("\
Kinematics Model (*.kte_dk.rkx *.kte_dk.rkb *.kte_dk.pbuf *.kte_ik.rkx *.kte_ik.rkb *.kte_ik.pbuf);;\
Dynamics Model (*.kte_dyn.rkx *.kte_dyn.rkb *.kte_dyn.pbuf);;\
KTE-chain Geometry Specification (*.geom.rkx *.geom.rkb *.geom.pbuf);;\
Geometric Model (*.geom_mdl.rkx *.geom_mdl.rkb *.geom_mdl.pbuf);;\
Proximity Model (*.prox_mdl.rkx *.prox_mdl.rkb *.prox_mdl.pbuf);;\
Complete Model (*.model.rkx *.model.rkb *.model.pbuf)"));
  
  if(fileName.isEmpty())
    return;
  
  QFileInfo fileInf(fileName);
  
  last_used_path = fileInf.absolutePath();
  
  QString fileExt = fileInf.completeSuffix();
  QString fileType = fileInf.suffix();
  QString fileContentExt = fileExt;
  fileContentExt.chop(fileType.length() + 1);
  
  
  if( fileType == tr("rkx") ) {
    ReaK::serialization::xml_iarchive in(fileName.toStdString());
    loadFromArchive(in, fileContentExt);
  } else if( fileType == tr("rkb") ) {
    ReaK::serialization::bin_iarchive in(fileName.toStdString());
    loadFromArchive(in, fileContentExt);
  } else if( fileType == tr("pbuf") ) {
    ReaK::serialization::protobuf_iarchive in(fileName.toStdString());
    loadFromArchive(in, fileContentExt);
  } else {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
  };
  
};


void KTEModelViewerEditor::onSave() {
  
  ReaK::serialization::object_node_desc current_obj_id = propedit.mdl.get_current_node();
  std::string current_obj_name = objtree_edit.get_object_name(current_obj_id);
  
  QString diag_title;
  QString diag_filter;
  
  std::map< std::string, ReaK::shared_ptr< ReaK::geom::kte_chain_geometry_3D > >::iterator kte_geom_it = kte_geometries.find(current_obj_name);
  if(kte_geom_it != kte_geometries.end()) {
    diag_title = tr("Save KTE-chain Geometry...");
    diag_filter = tr("\
KTE-chain Geometry Specification - XML (*.geom.rkx);;\
KTE-chain Geometry Specification - Binary (*.geom.rkb);;\
KTE-chain Geometry Specification - Protobuf (*.geom.pbuf)");
  };
  
  std::map< std::string, ReaK::shared_ptr< ReaK::geom::colored_model_3D > >::iterator geom_it = geom_models.find(current_obj_name);
  if(geom_it != geom_models.end()) {
    diag_title = tr("Save Geometric Model...");
    diag_filter = tr("\
Geometric Model - XML (*.geom_mdl.rkx);;\
Geometric Model - Binary (*.geom_mdl.rkb);;\
Geometric Model - Protobuf (*.geom_mdl.pbuf)");
  };
  
  std::map< std::string, ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > >::iterator prox_it = proxy_models.find(current_obj_name);
  if(prox_it != proxy_models.end()) {
    diag_title = tr("Save Proximity Model...");
    diag_filter = tr("\
Proximity Model - XML (*.prox_mdl.rkx);;\
Proximity Model - Binary (*.prox_mdl.rkb);;\
Proximity Model - Protobuf (*.prox_mdl.pbuf)");
  };
  
  std::map< std::string, ReaK::shared_ptr< ReaK::kte::inverse_dynamics_model > >::iterator dyn_it = dyn_models.find(current_obj_name);
  std::map< std::string, ReaK::shared_ptr< ReaK::kte::inverse_kinematics_model > >::iterator ik_it = ik_models.find(current_obj_name);
  std::map< std::string, ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > >::iterator dk_it = dk_models.find(current_obj_name);
  
  if(dyn_it != dyn_models.end()) {
    diag_title = tr("Save Dynamics Model...");
    diag_filter = tr("\
Dynamics Model - XML (*.kte_dyn.rkx);;\
Dynamics Model - Binary (*.kte_dyn.rkb);;\
Dynamics Model - Protobuf (*.kte_dyn.pbuf)");
  } else if(ik_it != ik_models.end()) {
    diag_title = tr("Save Kinematics Model...");
    diag_filter = tr("\
Inverse Kinematics Model - XML (*.kte_ik.rkx);;\
Inverse Kinematics Model - Binary (*.kte_ik.rkb);;\
Inverse Kinematics Model - Protobuf (*.kte_ik.pbuf);;\
Direct Kinematics Model - XML (*.kte_dk.rkx);;\
Direct Kinematics Model - Binary (*.kte_dk.rkb);;\
Direct Kinematics Model - Protobuf (*.kte_dk.pbuf)");
  } else if(dk_it != dk_models.end()) {
    diag_title = tr("Save Kinematics Model...");
    diag_filter = tr("\
Direct Kinematics Model - XML (*.kte_dk.rkx);;\
Direct Kinematics Model - Binary (*.kte_dk.rkb);;\
Direct Kinematics Model - Protobuf (*.kte_dk.pbuf)");
  };
  
  if( diag_title.isEmpty() ) {
    diag_title = tr("Save in ReaK Archive...");
    diag_filter = tr("\
ReaK Archive - XML (*.rkx);;\
ReaK Archive - Binary (*.rkb);;\
ReaK Archive - Protobuf (*.pbuf)");
  };
  
  QString fileName = QFileDialog::getSaveFileName(
    this, 
    diag_title,
    last_used_path,
    diag_filter);
  
  if( fileName == tr("") )
    return;
  
  QFileInfo fileInf(fileName);
  last_used_path = fileInf.absolutePath();
  QString fileExt = fileInf.suffix();
  
  ReaK::shared_ptr< ReaK::serialization::oarchive > out_ar;
  if( fileExt == tr("rkb") ) {
    out_ar = ReaK::shared_ptr< ReaK::serialization::oarchive >(new ReaK::serialization::bin_oarchive(fileName.toStdString()));
  } else if( fileExt == tr("pbuf") ) {
    out_ar = ReaK::shared_ptr< ReaK::serialization::oarchive >(new ReaK::serialization::protobuf_oarchive(fileName.toStdString()));
  } else {
    out_ar = ReaK::shared_ptr< ReaK::serialization::oarchive >(new ReaK::serialization::xml_oarchive(fileName.toStdString()));
  };
  
  if(kte_geom_it != kte_geometries.end()) {
    saveKTEChainGeometry(*out_ar, current_obj_name, kte_geom_it->second);
  } else if(geom_it != geom_models.end()) {
    saveGeometricModel(*out_ar, current_obj_name, geom_it->second);
  } else if(prox_it != proxy_models.end()) {
    saveProximityModel(*out_ar, current_obj_name, prox_it->second);
  } else if(dyn_it != dyn_models.end()) {
    saveDynamicsModel(*out_ar, current_obj_name, dyn_it->second);
  } else if(ik_it != ik_models.end()) {
    saveInverseKinModel(*out_ar, current_obj_name, ik_it->second);
  } else if(dk_it != dk_models.end()) {
    saveDirectKinModel(*out_ar, current_obj_name, dk_it->second);
  } else {
    ReaK::shared_ptr< ReaK::serialization::serializable > current_obj_ptr = (*objtree_graph)[current_obj_id].p_obj;
    (*out_ar) << current_obj_ptr;
  };
  
};

void KTEModelViewerEditor::onCloseAll() {
  
  objtree_edit.remove_object(objtree_root);
  objtree.mdl.refreshObjTree();
  propedit.mdl.selectObject(objtree_root);
  
  
  view3d_menu.removeDisplayGroup("Kinematics Models");
  view3d_menu.removeDisplayGroup("Geometric Models");
  view3d_menu.removeDisplayGroup("Proximity Models");
  
  sg_root->removeAllChildren();
  
  scene_graphs.clear();
  
  kte_geometries.clear();
  geom_models.clear();
  proxy_models.clear();
  
  mdl_base_frames.clear();
  mdl_jt_limits.clear();
  dk_models.clear();
  ik_models.clear();
  dyn_models.clear();
  
  mdl_to_base.clear();
  mdl_to_jt_lim.clear();
  
};

void KTEModelViewerEditor::onRefreshView() {
  
  for(std::map< std::string, ReaK::shared_ptr< ReaK::kte::direct_kinematics_model > >::iterator it = dk_models.begin(); it != dk_models.end(); ++it) {
    it->second->doDirectMotion();
  };
  
};




int main(int argc, char** argv) {
  QApplication app(argc,argv);
  
  last_used_path = QDir::currentPath();
  
  KTEModelViewerEditor window;
  window.show();
  // Pop up the main window.
  SoQt::show(&window);
  // Loop until exit.
  SoQt::mainLoop();
  
  return 0;
  //return app.exec();
};








