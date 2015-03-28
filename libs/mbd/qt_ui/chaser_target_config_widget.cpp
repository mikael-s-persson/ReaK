
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

#include <ReaK/mbd/qt/chaser_target_config_widget.hpp>

#include <QDockWidget>
#include <QListView>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QScrollArea>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include "ui_chaser_target_mdl_config.h"

namespace ReaK {

namespace qt {


static QString last_used_path;


ChaserTargetConfigWidget::ChaserTargetConfigWidget( View3DMenu* aView3dMenu, QWidget* parent, Qt::WindowFlags flags )
    : QDockWidget( tr( "Models" ), parent, flags ), ui( new Ui::ChaserTargetMdlConfig() ), view3d_menu( aView3dMenu ),
      sceneData() {
  QScrollArea* dock_scroll = new QScrollArea( this );
  dock_scroll->setWidgetResizable( true );
  QWidget* dock_wid = new QWidget( this );
  dock_scroll->setWidget( dock_wid );
  this->QDockWidget::setWidget( dock_scroll );
  ui->setupUi( dock_wid );

  connect( ui->actionLoadChaserMdl, SIGNAL( triggered() ), this, SLOT( loadChaserMdl() ) );
  connect( ui->actionEditChaserMdl, SIGNAL( triggered() ), this, SLOT( editChaserMdl() ) );
  connect( ui->actionSaveChaserMdl, SIGNAL( triggered() ), this, SLOT( saveChaserMdl() ) );

  connect( ui->actionLoadTargetMdl, SIGNAL( triggered() ), this, SLOT( loadTargetMdl() ) );
  connect( ui->actionEditTargetMdl, SIGNAL( triggered() ), this, SLOT( editTargetMdl() ) );
  connect( ui->actionSaveTargetMdl, SIGNAL( triggered() ), this, SLOT( saveTargetMdl() ) );

  connect( ui->actionEnvGeomAdd, SIGNAL( triggered() ), this, SLOT( addEnvMdl() ) );
  connect( ui->actionEnvGeomEdit, SIGNAL( triggered() ), this, SLOT( editEnvMdl() ) );
  connect( ui->actionEnvGeomClear, SIGNAL( triggered() ), this, SLOT( clearEnvMdls() ) );
  connect( ui->actionEnvGeomSave, SIGNAL( triggered() ), this, SLOT( saveEnvMdl() ) );

  connect( ui->actionLoadCompleteMdl, SIGNAL( triggered() ), this, SLOT( loadCompleteMdl() ) );
  connect( ui->actionEditCompleteMdl, SIGNAL( triggered() ), this, SLOT( editCompleteMdl() ) );
  connect( ui->actionSaveCompleteMdl, SIGNAL( triggered() ), this, SLOT( saveCompleteMdl() ) );
};

ChaserTargetConfigWidget::~ChaserTargetConfigWidget() {
  delete static_cast< QScrollArea* >( this->QDockWidget::widget() )->widget();
  delete this->QDockWidget::widget();
  delete ui;
};


void ChaserTargetConfigWidget::loadChaserMdl() {

  QString fileName
    = QFileDialog::getOpenFileName( this, tr( "Open Chaser Kinematic Model..." ), last_used_path,
                                    tr( "Chaser Kinematic Model (*.model.rkx *.model.rkb *.model.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.load_chaser( fileName.toStdString() ); // "models/CRS_A465.model.rkx"

    if( !sceneData.chaser_kin_model ) {
      QMessageBox::information( this, "Error!", "An error occurred when loading the file! No chaser model was found!",
                                QMessageBox::Ok );
      return;
    };

    ui->chaser_filename_edit->setText( QString::fromStdString( sceneData.chaser_kin_model->getName() ) );

    if( view3d_menu ) {
      shared_ptr< geom::oi_scene_graph > psg = view3d_menu->getGeometryGroup( "Chaser Geometry" );
      psg->clearAll();
      ( *psg ) << ( *sceneData.chaser_geom_model );

      shared_ptr< geom::oi_scene_graph > psg_kte = view3d_menu->getGeometryGroup( "Chaser KTE Chain" );
      psg_kte->clearAll();
      psg_kte->setCharacteristicLength( psg->computeCharacteristicLength() );
      ( *psg_kte ) << ( *sceneData.chaser_kin_model->getKTEChain() );
    };

    emit onChaserLoaded();

  } catch( std::exception& e ) {
    std::cout << "Got exception during loading chaser: " << e.what() << std::endl;
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };
};

void ChaserTargetConfigWidget::editChaserMdl(){

};

void ChaserTargetConfigWidget::saveChaserMdl() {

  QString fileName
    = QFileDialog::getSaveFileName( this, tr( "Save Chaser Kinematic Model..." ), last_used_path,
                                    tr( "Chaser Kinematic Model (*.model.rkx *.model.rkb *.model.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.save_chaser( fileName.toStdString() );
  } catch( std::exception& e ) {
    std::cout << "Got exception during saving chaser: " << e.what() << std::endl;
    QMessageBox::information( this, "Error!", "An error occurred while saving the chaser model to file!",
                              QMessageBox::Ok );
    return;
  };
};


void ChaserTargetConfigWidget::loadTargetMdl() {

  QString fileName = QFileDialog::getOpenFileName( this, tr( "Open Target Model..." ), last_used_path,
                                                   tr( "Target Model (*.model.rkx *.model.rkb *.model.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.load_target( fileName.toStdString() ); // "models/airship3D.model.rkx"

    if( !sceneData.target_kin_model ) {
      QMessageBox::information( this, "Error!", "An error occurred when loading the file! No target model was found!",
                                QMessageBox::Ok );
      return;
    };

    ui->target_filename_edit->setText( QString::fromStdString( sceneData.target_kin_model->getName() ) );

    if( view3d_menu ) {
      shared_ptr< geom::oi_scene_graph > psg = view3d_menu->getGeometryGroup( "Target Geometry" );
      psg->clearAll();
      ( *psg ) << ( *sceneData.target_geom_model );
    };

    emit onTargetLoaded();

  } catch( std::exception& e ) {
    std::cout << "Got exception during loading target: " << e.what() << std::endl;
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };
};

void ChaserTargetConfigWidget::editTargetMdl(){

};

void ChaserTargetConfigWidget::saveTargetMdl() {

  QString fileName = QFileDialog::getSaveFileName( this, tr( "Save Target Model..." ), last_used_path,
                                                   tr( "Target Model (*.model.rkx *.model.rkb *.model.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.save_target( fileName.toStdString() );
  } catch( std::exception& e ) {
    std::cout << "Got exception during saving target: " << e.what() << std::endl;
    QMessageBox::information( this, "Error!", "An error occurred while saving the target model to file!",
                              QMessageBox::Ok );
    return;
  };
};


void ChaserTargetConfigWidget::addEnvMdl() {

  QString fileName = QFileDialog::getOpenFileName( this, tr( "Open Environment Geometry..." ), last_used_path,
                                                   tr( "Environment Geometry (*.geom.rkx *.geom.rkb *.geom.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.load_environment( fileName.toStdString() ); // "models/MD148_lab.geom.rkx"

    for( std::size_t i = ui->env_geoms_list->count(); i < sceneData.env_geom_models.size(); ++i )
      ui->env_geoms_list->addItem( QString::fromStdString( sceneData.env_geom_models[i]->getName() ) );

    if( view3d_menu ) {
      shared_ptr< geom::oi_scene_graph > psg = view3d_menu->getGeometryGroup( "Environment" );
      psg->clearAll();
      for( std::size_t i = 0; i < sceneData.env_geom_models.size(); ++i )
        ( *psg ) << ( *( sceneData.env_geom_models[i] ) );
    };

  } catch( std::exception& e ) {
    std::cout << "Got exception during loading environment element: " << e.what() << std::endl;
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };
};

void ChaserTargetConfigWidget::editEnvMdl(){

};

void ChaserTargetConfigWidget::clearEnvMdls() {
  if( view3d_menu )
    view3d_menu->getGeometryGroup( "Environment" )->clearAll();
  ui->env_geoms_list->clear();
  sceneData.clear_environment();
};

void ChaserTargetConfigWidget::saveEnvMdl() {

  if( ui->env_geoms_list->count() == 0 ) {
    QMessageBox::information( this, "Error!", "There are no environment geometries!", QMessageBox::Ok );
    return;
  };

  QString fileName = QFileDialog::getSaveFileName( this, tr( "Save Environment Geometry..." ), last_used_path,
                                                   tr( "Environment Geometry (*.geom.rkx *.geom.rkb *.geom.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  try {
    sceneData.save_environment( ui->env_geoms_list->currentRow(), fileName.toStdString() );
  } catch( std::exception& e ) {
    std::cout << "Got exception during saving environment element: " << e.what() << std::endl;
    QMessageBox::information(
      this, "Error!", "An error occurred while saving the environment geometry element to file!", QMessageBox::Ok );
    return;
  };
};


void ChaserTargetConfigWidget::loadCompleteMdl() {

  QString fileName = QFileDialog::getOpenFileName( this, tr( "Open Chaser-Target Scenario..." ), last_used_path,
                                                   tr( "Chaser-Target Scenario (*.rkx *.rkb *.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  loadCompleteModel( fileName.toStdString() );
};

void ChaserTargetConfigWidget::editCompleteMdl(){

};

void ChaserTargetConfigWidget::saveCompleteMdl() {

  QString fileName = QFileDialog::getSaveFileName( this, tr( "Save Chaser-Target Scenario..." ), last_used_path,
                                                   tr( "Chaser-Target Scenario (*.rkx *.rkb *.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  saveCompleteModel( fileName.toStdString() );
};

void ChaserTargetConfigWidget::loadCompleteModel( const std::string& aFilename ) {

  try {
    ( *serialization::open_iarchive( aFilename ) ) >> sceneData;

    ui->complete_filename_edit->setText( QString::fromStdString( aFilename ) );

    if( !sceneData.chaser_kin_model ) {
      QMessageBox::information( this, "Error!", "An error occurred when loading the file! No chaser model was found!",
                                QMessageBox::Ok );
    } else {
      ui->chaser_filename_edit->setText( QString::fromStdString( sceneData.chaser_kin_model->getName() ) );

      if( view3d_menu ) {
        shared_ptr< geom::oi_scene_graph > psg_chase = view3d_menu->getGeometryGroup( "Chaser Geometry" );
        psg_chase->clearAll();
        ( *psg_chase ) << ( *sceneData.chaser_geom_model );

        shared_ptr< geom::oi_scene_graph > psg_kte = view3d_menu->getGeometryGroup( "Chaser KTE Chain" );
        psg_kte->clearAll();
        psg_kte->setCharacteristicLength( psg_chase->computeCharacteristicLength() );
        ( *psg_kte ) << ( *sceneData.chaser_kin_model->getKTEChain() );
      };

      emit onChaserLoaded();
    };

    if( !sceneData.target_kin_model ) {
      QMessageBox::information( this, "Error!", "An error occurred when loading the file! No target model was found!",
                                QMessageBox::Ok );
    } else {
      ui->target_filename_edit->setText( QString::fromStdString( sceneData.target_kin_model->getName() ) );

      if( view3d_menu ) {
        shared_ptr< geom::oi_scene_graph > psg_target = view3d_menu->getGeometryGroup( "Target Geometry" );
        psg_target->clearAll();
        ( *psg_target ) << ( *sceneData.target_geom_model );
      };

      emit onTargetLoaded();
    };

    ui->env_geoms_list->clear();
    for( std::size_t i = 0; i < sceneData.env_geom_models.size(); ++i )
      ui->env_geoms_list->addItem( QString::fromStdString( sceneData.env_geom_models[i]->getName() ) );

    if( view3d_menu ) {

      shared_ptr< geom::oi_scene_graph > psg_env = view3d_menu->getGeometryGroup( "Environment" );
      psg_env->clearAll();
      for( std::size_t i = 0; i < sceneData.env_geom_models.size(); ++i )
        ( *psg_env ) << ( *( sceneData.env_geom_models[i] ) );
    };

  } catch( std::exception& e ) {
    std::cout << "Got exception during loading complete scenario: " << e.what() << std::endl;
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };
};

void ChaserTargetConfigWidget::saveCompleteModel( const std::string& aFilename ) {

  try {
    ( *serialization::open_oarchive( aFilename ) ) << sceneData;
    ui->complete_filename_edit->setText( QString::fromStdString( aFilename ) );
  } catch( std::exception& e ) {
    std::cout << "Got exception during saving complete scenario: " << e.what() << std::endl;
    QMessageBox::information( this, "Error!", "An error occurred while saving the chaser model to file!",
                              QMessageBox::Ok );
    return;
  };
};
};
};
