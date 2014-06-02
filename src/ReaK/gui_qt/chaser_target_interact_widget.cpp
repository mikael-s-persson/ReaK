
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

#include <ReaK/gui_qt/chaser_target_interact_widget.hpp>

#include <QDockWidget>
#include <QMessageBox>
#include <QFileInfo>
#include <QFileDialog>
#include <QString>
#include <QScrollArea>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/core/optimization/optim_exceptions.hpp>


namespace ReaK {
  
namespace rkqt {

static QString last_used_path;


ChaserTargetInteractWidget::ChaserTargetInteractWidget(
  kte::chaser_target_data* aPSceneData, 
  QWidget* parent, 
  Qt::WindowFlags flags) :
  QDockWidget(tr("Interact"), parent, flags),
  Ui::ChaserTargetInteract(),
  pSceneData(aPSceneData)
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  setupUi(dock_wid);
  
  connect(actionJointChange, SIGNAL(triggered()), this, SLOT(onJointChange()));
  connect(actionTargetChange, SIGNAL(triggered()), this, SLOT(onTargetChange()));
  
  connect(load_traj_button, SIGNAL(clicked()), this, SLOT(loadTargetTrajectory()));
  
  connect(load_button, SIGNAL(clicked()), this, SLOT(loadPositions()));
  connect(save_button, SIGNAL(clicked()), this, SLOT(savePositions()));
  
};

ChaserTargetInteractWidget::~ChaserTargetInteractWidget() { 
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
};


void ChaserTargetInteractWidget::onJointChange() {
  if( !pSceneData->chaser_kin_model )
    return;
  
  pSceneData->chaser_kin_model->setJointPositions(
    vect_n<double>(
      double(this->track_pos->value()) * 0.001,
      double(this->joint1_pos->value()) * 0.001,
      double(this->joint2_pos->value()) * 0.001,
      double(this->joint3_pos->value()) * 0.001,
      double(this->joint4_pos->value()) * 0.001,
      double(this->joint5_pos->value()) * 0.001,
      double(this->joint6_pos->value()) * 0.001
    )
  );
  pSceneData->chaser_kin_model->doDirectMotion();
};


void ChaserTargetInteractWidget::onTargetChange() {
  if( !pSceneData->target_kin_model )
    return;
  
  shared_ptr< frame_3D<double> > target_state = pSceneData->target_kin_model->getFrame3D(0);
  target_state->Position = vect<double,3>(
    double(this->target_x->value()) * 0.001, 
    double(this->target_y->value()) * 0.001, 
    double(this->target_z->value()) * 0.001);
  target_state->Quat = 
    quaternion<double>::zrot(double(this->target_yaw->value()) * 0.001) * 
    quaternion<double>::yrot(double(this->target_pitch->value()) * 0.001) * 
    quaternion<double>::xrot(double(this->target_roll->value()) * 0.001);
  pSceneData->target_kin_model->doDirectMotion();
  
  
  if( pSceneData->chaser_kin_model && this->enable_ik_check->isChecked()) {
    try {
      *(pSceneData->chaser_kin_model->getDependentFrame3D(0)->mFrame) = *(pSceneData->target_frame);
      pSceneData->chaser_kin_model->doInverseMotion();
    } catch( optim::infeasible_problem& e ) { RK_UNUSED(e); };
    pSceneData->chaser_kin_model->doDirectMotion();
  };
  
};


void ChaserTargetInteractWidget::loadTargetTrajectory() {
  
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Open Target Trajectory..."),
    last_used_path,
    tr("SE(3) Trajectories (*.rkx *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  QFileInfo fileInf(fileName);
  
  last_used_path = fileInf.absolutePath();
  
  emit onLoadTargetTrajectory(fileName);
  
  this->traj_filename_edit->setText(fileInf.baseName());
  
};



void ChaserTargetInteractWidget::savePositions() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Positions Record..."), last_used_path,
    tr("Chaser-Target Positions Record (*.ctpos.rkx *.ctpos.rkb *.ctpos.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  vect<double,7> chaser_positions;
  chaser_positions[0] = double(this->track_pos->value()) * 0.001;
  chaser_positions[1] = double(this->joint1_pos->value()) * 0.001;
  chaser_positions[2] = double(this->joint2_pos->value()) * 0.001;
  chaser_positions[3] = double(this->joint3_pos->value()) * 0.001;
  chaser_positions[4] = double(this->joint4_pos->value()) * 0.001;
  chaser_positions[5] = double(this->joint5_pos->value()) * 0.001;
  chaser_positions[6] = double(this->joint6_pos->value()) * 0.001;
  
  vect<double,3> target_position;
  target_position[0] = double(this->target_x->value()) * 0.001;
  target_position[1] = double(this->target_y->value()) * 0.001;
  target_position[2] = double(this->target_z->value()) * 0.001;
  
  euler_angles_TB<double> ea;
  ea.yaw() = double(this->target_yaw->value()) * 0.001;
  ea.pitch() = double(this->target_pitch->value()) * 0.001;
  ea.roll() = double(this->target_roll->value()) * 0.001;
  quaternion<double> target_quaternion = ea.getQuaternion();
  
  
  try {
    *(serialization::open_oarchive(fileName.toStdString()))
      & RK_SERIAL_SAVE_WITH_NAME(chaser_positions)
      & RK_SERIAL_SAVE_WITH_NAME(target_position)
      & RK_SERIAL_SAVE_WITH_NAME(target_quaternion);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void ChaserTargetInteractWidget::loadPositions() {
  QString fileName = QFileDialog::getOpenFileName(
    this, 
    tr("Open Positions Record..."),
    last_used_path,
    tr("Chaser-Target Positions Record (*.ctpos.rkx *.ctpos.rkb *.ctpos.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  vect<double,7> chaser_positions;
  vect<double,3> target_position;
  quaternion<double> target_quaternion;
  
  try {
    *(serialization::open_iarchive(fileName.toStdString()))
      & RK_SERIAL_LOAD_WITH_NAME(chaser_positions)
      & RK_SERIAL_LOAD_WITH_NAME(target_position)
      & RK_SERIAL_LOAD_WITH_NAME(target_quaternion);
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  this->track_pos->setValue(int(1000.0 * chaser_positions[0]));
  this->joint1_pos->setValue(int(1000.0 * chaser_positions[1]));
  this->joint2_pos->setValue(int(1000.0 * chaser_positions[2]));
  this->joint3_pos->setValue(int(1000.0 * chaser_positions[3]));
  this->joint4_pos->setValue(int(1000.0 * chaser_positions[4]));
  this->joint5_pos->setValue(int(1000.0 * chaser_positions[5]));
  this->joint6_pos->setValue(int(1000.0 * chaser_positions[6])); 
  onJointChange();
  
  this->target_x->setValue(int(1000.0 * target_position[0]));
  this->target_y->setValue(int(1000.0 * target_position[1]));
  this->target_z->setValue(int(1000.0 * target_position[2]));
  euler_angles_TB<double> ea = euler_angles_TB<double>(target_quaternion);
  this->target_yaw->setValue(int(1000.0 * ea.yaw()));
  this->target_pitch->setValue(int(1000.0 * ea.pitch()));
  this->target_roll->setValue(int(1000.0 * ea.roll()));
  onTargetChange();
  
};



void ChaserTargetInteractWidget::loadJointPosFromModel() {
  if( !pSceneData->chaser_kin_model )
    return;
  
  vect_n<double> v = pSceneData->chaser_kin_model->getJointPositions();
  
  this->track_pos->setValue(1000.0 * v[0]);
  this->joint1_pos->setValue(1000.0 * v[1]);
  this->joint2_pos->setValue(1000.0 * v[2]);
  this->joint3_pos->setValue(1000.0 * v[3]);
  this->joint4_pos->setValue(1000.0 * v[4]);
  this->joint5_pos->setValue(1000.0 * v[5]);
  this->joint6_pos->setValue(1000.0 * v[6]);
  
  onJointChange();
};


void ChaserTargetInteractWidget::loadTargetPosFromModel() {
  if( !pSceneData->target_kin_model )
    return;
  
  shared_ptr< frame_3D<double> > target_state = pSceneData->target_kin_model->getFrame3D(0);
  this->target_x->setValue(int(1000.0 * target_state->Position[0]));
  this->target_y->setValue(int(1000.0 * target_state->Position[1]));
  this->target_z->setValue(int(1000.0 * target_state->Position[2]));
  euler_angles_TB<double> ea = euler_angles_TB<double>(target_state->Quat);
  this->target_yaw->setValue(int(1000.0 * ea.yaw()));
  this->target_pitch->setValue(int(1000.0 * ea.pitch()));
  this->target_roll->setValue(int(1000.0 * ea.roll()));
  
  
  onTargetChange();
};




bool ChaserTargetInteractWidget::isIKEnabled() const {
  return this->enable_ik_check->isChecked();
};




};

};














