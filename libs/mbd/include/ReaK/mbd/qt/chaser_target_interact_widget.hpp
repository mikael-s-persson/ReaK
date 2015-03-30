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

#include <ReaK/mbd/models/chaser_target_model_data.hpp>

#include <QDockWidget>

namespace Ui {
class ChaserTargetInteract;
};

namespace ReaK {

namespace qt {

class ChaserTargetInteractWidget : public QDockWidget {
  Q_OBJECT

private:
  Ui::ChaserTargetInteract* ui;

public:
  ChaserTargetInteractWidget( kte::chaser_target_data* aPSceneData, QWidget* parent = nullptr,
                              Qt::WindowFlags flags = 0 );
  virtual ~ChaserTargetInteractWidget();

public slots:

  void savePositions();
  void loadPositions();

  void onJointChange();
  void onTargetChange();

  void loadTargetTrajectory();

public slots:

  void loadJointPosFromModel();
  void loadTargetPosFromModel();

signals:

  void onLoadTargetTrajectory( QString );

public:
  kte::chaser_target_data* pSceneData;

  bool isIKEnabled() const;

  void saveChaserTargetPositions( const std::string& aFilename );
  void loadChaserTargetPositions( const std::string& aFilename );
};
};
};

#endif
