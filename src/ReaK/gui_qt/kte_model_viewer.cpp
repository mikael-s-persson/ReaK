
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

#include "topologies/joint_space_limits.hpp"

#include "serialization/xml_archiver.hpp"
#include "serialization/bin_archiver.hpp"
#include "serialization/protobuf_archiver.hpp"

#include "optimization/optim_exceptions.hpp"


struct env_element {
  ReaK::shared_ptr< ReaK::frame_3D<double> > mdl_base;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > mdl_kin_chain;
  ReaK::geom::oi_scene_graph* mdl_kin_chain_sg;
  SoSwitch* mdl_kin_chain_sw;
  ReaK::shared_ptr< ReaK::geom::colored_model_3D > mdl_render;
  ReaK::geom::oi_scene_graph* mdl_render_sg;
  SoSwitch* mdl_render_switch;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > mdl_proxy;
};


static QString last_used_path;



KTEModelViewerEditor::KTEModelViewerEditor(QWidget * parent, Qt::WindowFlags flags ) : 
                                           QMainWindow(parent,flags),
                                           objtree_graph(ReaK::shared_ptr< ReaK::serialization::object_graph >(new ReaK::serialization::object_graph())),
                                           objtree_root(add_vertex(*objtree_graph)),
                                           objtree_out_arc(objtree_graph, objtree_root),
                                           objtree_sch_bld(),
                                           objtree(objtree_graph, objtree_root),
                                           propedit(&(objtree.mdl)) {
  setupUi(this);
  using namespace ReaK;
  
  addDockWidget(Qt::LeftDockWidgetArea, &objtree);
  
  addDockWidget(Qt::LeftDockWidgetArea, &propedit);
  
  tabifyDockWidget(&objtree, &propedit);
  
  
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
  connect(actionSave_as, SIGNAL(triggered()), this, SLOT(onSaveAs()));
  connect(actionClose_all, SIGNAL(triggered()), this, SLOT(onCloseAll()));
  
  actionLoad->setIcon( style()->standardIcon( QStyle::SP_DialogOpenButton, NULL, actionLoad->menu() ) );
  actionSave->setIcon( style()->standardIcon( QStyle::SP_DialogSaveButton, NULL, actionSave->menu() ) );
  actionSave_as->setIcon( style()->standardIcon( QStyle::SP_DialogSaveButton, NULL, actionSave_as->menu() ) );
  actionClose_all->setIcon( style()->standardIcon( QStyle::SP_DialogCloseButton, NULL, actionClose_all->menu() ) );
  actionQuit->setIcon( style()->standardIcon( QStyle::SP_DialogCloseButton, NULL, actionQuit->menu() ) );
  
  
  
  SoQt::init(this->centralwidget);
  
  sg_root = new SoSeparator;
  sg_root->ref();
  
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
  
  SoQt::done();
  
};


void KTEModelViewerEditor::loadFromArchive(ReaK::serialization::iarchive& in, QString fileContentExt) {
  
  if( (fileContentExt == tr("kte_ik")) || (fileContentExt == tr("kte_dk")) ) {
    ReaK::shared_ptr< ReaK::pose_3D<double> > base_frame;
    
    
    
    
  } else if( fileContentExt == tr("kte_dyn") ) {
    
    
  } else if( fileContentExt == tr("geom") ) {
    
    
  } else if( fileContentExt == tr("model") ) {
    
    
  };
  
};


void KTEModelViewerEditor::onLoad() {
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Load a Model..."),
    last_used_path,
    tr("\
Kinematics Model (*.kte_dk.rkx *.kte_dk.rkb *.kte_dk.pbuf *.kte_ik.rkx *.kte_ik.rkb *.kte_ik.pbuf);;\
Dynamics Model (*.kte_dyn.rkx *.kte_dyn.rkb *.kte_dyn.pbuf);;\
Geometric Model (*.geom.rkx *.geom.rkb *.geom.pbuf *.proxy.rkx *.proxy.rkb *.proxy.pbuf);;\
Complete Model (*.model.rkx *.model.rkb *.model.pbuf)"));
  
  if(fileName.isEmpty())
    return;
  
  QFileInfo fileInf(fileName);
  
  last_used_path = fileInf.absolutePath();
  
  QString fileExt = fileInf.completeSuffix();
  QString fileType = fileInf.suffix();
  QString fileContentExt = fileExt.left(fileContentExt.lastIndexOf('.'));
  
  
  if( fileExt == tr("rkx") ) {
    ReaK::serialization::xml_iarchive in(fileName.toStdString());
    loadFromArchive(in, fileContentExt);
  } else if( fileExt == tr("rkb") ) {
    ReaK::serialization::bin_iarchive in(fileName.toStdString());
    loadFromArchive(in, fileContentExt);
  } else if( fileExt == tr("pbuf") ) {
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
  
  
};

void KTEModelViewerEditor::onSaveAs() {
  
  
};

void KTEModelViewerEditor::onCloseAll() {
  
  
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








