/**
 *
 *
 *
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_DATA_STREAM_OPT_WIDGET_HPP
#define REAK_DATA_STREAM_OPT_WIDGET_HPP

#include <ReaK/core/recorders/data_record_options.hpp>

#include <QDockWidget>

#include <QMainWindow>

namespace Ui {
class DataStreamOpt;
};

namespace ReaK::qt {

class DataStreamOptWidget : public QDockWidget {
  Q_OBJECT

 public:
  DataStreamOptWidget(
      recorder::data_stream_options aDataOpt = recorder::data_stream_options(),
      QWidget* parent = nullptr, Qt::WindowFlags flags = 0);
  virtual ~DataStreamOptWidget();

 private slots:

  void onUpdateURIAndDataOpt();
  void onUpdateFieldsAndDataOpt();
  void onUpdateDataOptNames();

 private:
  Ui::DataStreamOpt* ui;

 public:
  recorder::data_stream_options data_opt;

  void setEnabledNameEdits(bool aEnabled);
  void setEnabledFreqSetting(bool aEnabled);
  void setEnabledBufferSetting(bool aEnabled);
  void setEnabledTimeSyncSetting(bool aEnabled);

  enum stream_types { allow_all, file_only, network_only } mode;

  void setMode(stream_types aMode) {
    mode = aMode;
    onUpdateFieldsAndDataOpt();
    onUpdateDataOptNames();
  }
};

}  // namespace ReaK::qt

#endif
