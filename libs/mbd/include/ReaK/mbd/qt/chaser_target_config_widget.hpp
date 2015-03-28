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

#ifndef REAK_CHASER_TARGET_CONFIG_WIDGET_HPP
#define REAK_CHASER_TARGET_CONFIG_WIDGET_HPP

#include <ReaK/mbd/models/chaser_target_model_data.hpp>

#include "rk_view3d_menu.hpp"

#include <QDockWidget>

namespace Ui {
class ChaserTargetMdlConfig;
};

namespace ReaK {

namespace qt {

class ChaserTargetConfigWidget : public QDockWidget {
  Q_OBJECT

private:
  Ui::ChaserTargetMdlConfig* ui;

public:
  ChaserTargetConfigWidget( View3DMenu* aView3dMenu = NULL, QWidget* parent = NULL, Qt::WindowFlags flags = 0 );
  virtual ~ChaserTargetConfigWidget();

private slots:

  void loadChaserMdl();
  void editChaserMdl();
  void saveChaserMdl();

  void loadTargetMdl();
  void editTargetMdl();
  void saveTargetMdl();

  void addEnvMdl();
  void editEnvMdl();
  void clearEnvMdls();
  void saveEnvMdl();

  void loadCompleteMdl();
  void editCompleteMdl();
  void saveCompleteMdl();

public:
signals:

  void onTargetLoaded();
  void onChaserLoaded();

public:
  View3DMenu* view3d_menu;

  kte::chaser_target_data sceneData;

  void loadCompleteModel( const std::string& aFilename );
  void saveCompleteModel( const std::string& aFilename );
};
};
};

#endif
