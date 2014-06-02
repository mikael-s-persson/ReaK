/**
 * 
 * 
 * 
 */

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

#ifndef REAK_CHASER_TARGET_INTERACT_WIDGET_HPP
#define REAK_CHASER_TARGET_INTERACT_WIDGET_HPP

#include <ReaK/ctrl/kte_models/chaser_target_model_data.hpp>

#include "ui_chaser_target_interact.h"
#include <QDockWidget>

namespace ReaK {
  
namespace rkqt {

class ChaserTargetInteractWidget : public QDockWidget, private Ui::ChaserTargetInteract {
    Q_OBJECT
  
  public:
    ChaserTargetInteractWidget(kte::chaser_target_data* aPSceneData, QWidget * parent = NULL, Qt::WindowFlags flags = 0);
    virtual ~ChaserTargetInteractWidget();
    
  private slots:
    
    void savePositions();
    void loadPositions();
    
    void onJointChange();
    void onTargetChange();
    
    void loadTargetTrajectory();
    
  public slots:
    
    void loadJointPosFromModel();
    void loadTargetPosFromModel();
    
  signals:
    
    void onLoadTargetTrajectory(QString);
    
  public:
    
    kte::chaser_target_data* pSceneData;
    
    bool isIKEnabled() const;
    
};

};

};

#endif














