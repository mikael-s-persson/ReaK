
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

#include <ReaK/mbd/qt/chaser_target_interact_widget.hpp>

#include <QDockWidget>
#include <QMessageBox>
#include <QFileInfo>
#include <QFileDialog>
#include <QString>
#include <QScrollArea>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/math/optimization/optim_exceptions.hpp>

#include "ui_chaser_target_interact.h"

namespace ReaK {

namespace qt {

static QString last_used_path;


ChaserTargetInteractWidget::ChaserTargetInteractWidget( kte::chaser_target_data* aPSceneData, QWidget* parent,
                                                        Qt::WindowFlags flags )
    : QDockWidget( tr( "Interact" ), parent, flags ), ui( new Ui::ChaserTargetInteract() ), pSceneData( aPSceneData ) {
  QScrollArea* dock_scroll = new QScrollArea( this );
  dock_scroll->setWidgetResizable( true );
  QWidget* dock_wid = new QWidget( this );
  dock_scroll->setWidget( dock_wid );
  this->QDockWidget::setWidget( dock_scroll );
  ui->setupUi( dock_wid );

  connect( ui->actionJointChange, SIGNAL( triggered() ), this, SLOT( onJointChange() ) );
  connect( ui->actionTargetChange, SIGNAL( triggered() ), this, SLOT( onTargetChange() ) );

  connect( ui->load_traj_button, SIGNAL( clicked() ), this, SLOT( loadTargetTrajectory() ) );

  connect( ui->load_button, SIGNAL( clicked() ), this, SLOT( loadPositions() ) );
  connect( ui->save_button, SIGNAL( clicked() ), this, SLOT( savePositions() ) );
};

ChaserTargetInteractWidget::~ChaserTargetInteractWidget() {
  delete static_cast< QScrollArea* >( this->QDockWidget::widget() )->widget();
  delete this->QDockWidget::widget();
  delete ui;
};


void ChaserTargetInteractWidget::onJointChange() {
  if( !pSceneData->chaser_kin_model )
    return;

  pSceneData->chaser_kin_model->setJointPositions(
    vect_n< double >( double( ui->track_pos->value() ) * 0.001, double( ui->joint1_pos->value() ) * 0.001,
                      double( ui->joint2_pos->value() ) * 0.001, double( ui->joint3_pos->value() ) * 0.001,
                      double( ui->joint4_pos->value() ) * 0.001, double( ui->joint5_pos->value() ) * 0.001,
                      double( ui->joint6_pos->value() ) * 0.001 ) );
  pSceneData->chaser_kin_model->doDirectMotion();
};


void ChaserTargetInteractWidget::onTargetChange() {
  if( !pSceneData->target_kin_model )
    return;

  shared_ptr< frame_3D< double > > target_state = pSceneData->target_kin_model->getFrame3D( 0 );
  target_state->Position
    = vect< double, 3 >( double( ui->target_x->value() ) * 0.001, double( ui->target_y->value() ) * 0.001,
                         double( ui->target_z->value() ) * 0.001 );
  target_state->Quat = quaternion< double >::zrot( double( ui->target_yaw->value() ) * 0.001 )
                       * quaternion< double >::yrot( double( ui->target_pitch->value() ) * 0.001 )
                       * quaternion< double >::xrot( double( ui->target_roll->value() ) * 0.001 );
  pSceneData->target_kin_model->doDirectMotion();


  if( pSceneData->chaser_kin_model && ui->enable_ik_check->isChecked() ) {
    try {
      *( pSceneData->chaser_kin_model->getDependentFrame3D( 0 )->mFrame ) = *( pSceneData->target_frame );
      pSceneData->chaser_kin_model->doInverseMotion();
    } catch( optim::infeasible_problem& e ) {
      RK_UNUSED( e );
    };
    pSceneData->chaser_kin_model->doDirectMotion();
  };
};


void ChaserTargetInteractWidget::loadTargetTrajectory() {

  QString fileName = QFileDialog::getOpenFileName( this, tr( "Open Target Trajectory..." ), last_used_path,
                                                   tr( "SE(3) Trajectories (*.rkx *.rkb *.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  QFileInfo fileInf( fileName );

  last_used_path = fileInf.absolutePath();

  emit onLoadTargetTrajectory( fileName );

  ui->traj_filename_edit->setText( fileInf.baseName() );
};


void ChaserTargetInteractWidget::savePositions() {
  QString fileName
    = QFileDialog::getSaveFileName( this, tr( "Save Positions Record..." ), last_used_path,
                                    tr( "Chaser-Target Positions Record (*.ctpos.rkx *.ctpos.rkb *.ctpos.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  saveChaserTargetPositions( fileName.toStdString() );
};

void ChaserTargetInteractWidget::loadPositions() {
  QString fileName
    = QFileDialog::getOpenFileName( this, tr( "Open Positions Record..." ), last_used_path,
                                    tr( "Chaser-Target Positions Record (*.ctpos.rkx *.ctpos.rkb *.ctpos.pbuf)" ) );

  if( fileName == tr( "" ) )
    return;

  last_used_path = QFileInfo( fileName ).absolutePath();

  loadChaserTargetPositions( fileName.toStdString() );
};

void ChaserTargetInteractWidget::saveChaserTargetPositions( const std::string& aFilename ) {

  vect< double, 7 > chaser_positions;
  chaser_positions[0] = double( ui->track_pos->value() ) * 0.001;
  chaser_positions[1] = double( ui->joint1_pos->value() ) * 0.001;
  chaser_positions[2] = double( ui->joint2_pos->value() ) * 0.001;
  chaser_positions[3] = double( ui->joint3_pos->value() ) * 0.001;
  chaser_positions[4] = double( ui->joint4_pos->value() ) * 0.001;
  chaser_positions[5] = double( ui->joint5_pos->value() ) * 0.001;
  chaser_positions[6] = double( ui->joint6_pos->value() ) * 0.001;

  vect< double, 3 > target_position;
  target_position[0] = double( ui->target_x->value() ) * 0.001;
  target_position[1] = double( ui->target_y->value() ) * 0.001;
  target_position[2] = double( ui->target_z->value() ) * 0.001;

  euler_angles_TB< double > ea;
  ea.yaw() = double( ui->target_yaw->value() ) * 0.001;
  ea.pitch() = double( ui->target_pitch->value() ) * 0.001;
  ea.roll() = double( ui->target_roll->value() ) * 0.001;
  quaternion< double > target_quaternion = ea.getQuaternion();

  try {
    *( serialization::open_oarchive( aFilename ) ) & RK_SERIAL_SAVE_WITH_NAME( chaser_positions )
      & RK_SERIAL_SAVE_WITH_NAME( target_position ) & RK_SERIAL_SAVE_WITH_NAME( target_quaternion );
  } catch( ... ) {
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };
};

void ChaserTargetInteractWidget::loadChaserTargetPositions( const std::string& aFilename ) {

  vect< double, 7 > chaser_positions;
  vect< double, 3 > target_position;
  quaternion< double > target_quaternion;

  try {
    *( serialization::open_iarchive( aFilename ) ) & RK_SERIAL_LOAD_WITH_NAME( chaser_positions )
      & RK_SERIAL_LOAD_WITH_NAME( target_position ) & RK_SERIAL_LOAD_WITH_NAME( target_quaternion );
  } catch( ... ) {
    QMessageBox::information( this, "File Type Not Supported!", "Sorry, this file-type is not supported!",
                              QMessageBox::Ok );
    return;
  };

  ui->track_pos->setValue( int( 1000.0 * chaser_positions[0] ) );
  ui->joint1_pos->setValue( int( 1000.0 * chaser_positions[1] ) );
  ui->joint2_pos->setValue( int( 1000.0 * chaser_positions[2] ) );
  ui->joint3_pos->setValue( int( 1000.0 * chaser_positions[3] ) );
  ui->joint4_pos->setValue( int( 1000.0 * chaser_positions[4] ) );
  ui->joint5_pos->setValue( int( 1000.0 * chaser_positions[5] ) );
  ui->joint6_pos->setValue( int( 1000.0 * chaser_positions[6] ) );
  onJointChange();

  ui->target_x->setValue( int( 1000.0 * target_position[0] ) );
  ui->target_y->setValue( int( 1000.0 * target_position[1] ) );
  ui->target_z->setValue( int( 1000.0 * target_position[2] ) );
  euler_angles_TB< double > ea = euler_angles_TB< double >( target_quaternion );
  ui->target_yaw->setValue( int( 1000.0 * ea.yaw() ) );
  ui->target_pitch->setValue( int( 1000.0 * ea.pitch() ) );
  ui->target_roll->setValue( int( 1000.0 * ea.roll() ) );
  onTargetChange();
};


void ChaserTargetInteractWidget::loadJointPosFromModel() {
  if( !pSceneData->chaser_kin_model )
    return;

  vect_n< double > v = pSceneData->chaser_kin_model->getJointPositions();

  ui->track_pos->setValue( 1000.0 * v[0] );
  ui->joint1_pos->setValue( 1000.0 * v[1] );
  ui->joint2_pos->setValue( 1000.0 * v[2] );
  ui->joint3_pos->setValue( 1000.0 * v[3] );
  ui->joint4_pos->setValue( 1000.0 * v[4] );
  ui->joint5_pos->setValue( 1000.0 * v[5] );
  ui->joint6_pos->setValue( 1000.0 * v[6] );

  onJointChange();
};


void ChaserTargetInteractWidget::loadTargetPosFromModel() {
  if( !pSceneData->target_kin_model )
    return;

  shared_ptr< frame_3D< double > > target_state = pSceneData->target_kin_model->getFrame3D( 0 );
  ui->target_x->setValue( int( 1000.0 * target_state->Position[0] ) );
  ui->target_y->setValue( int( 1000.0 * target_state->Position[1] ) );
  ui->target_z->setValue( int( 1000.0 * target_state->Position[2] ) );
  euler_angles_TB< double > ea = euler_angles_TB< double >( target_state->Quat );
  ui->target_yaw->setValue( int( 1000.0 * ea.yaw() ) );
  ui->target_pitch->setValue( int( 1000.0 * ea.pitch() ) );
  ui->target_roll->setValue( int( 1000.0 * ea.roll() ) );

  onTargetChange();
};

bool ChaserTargetInteractWidget::isIKEnabled() const { return ui->enable_ik_check->isChecked(); };
};
};
