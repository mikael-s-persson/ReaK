
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

#include <ReaK/gui_qt/data_stream_opt_widget.hpp>

#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QMainWindow>
#include <QProcess>
#include <QScrollArea>

#if 0
#include <ReaK/core/serialization/archiver_factory.hpp>
#endif

#include <ReaK/core/recorders/data_record_options.hpp>


namespace ReaK {
  
namespace rkqt {



DataStreamOptWidget::DataStreamOptWidget(recorder::data_stream_options aDataOpt, QWidget * parent, Qt::WindowFlags flags) :
  QDockWidget(tr("Predictor"), parent, flags),
  Ui::DataStreamOpt(),
  data_opt(aDataOpt), mode(allow_all)
{
  QScrollArea* dock_scroll = new QScrollArea(this);
  dock_scroll->setWidgetResizable(true);
  QWidget* dock_wid = new QWidget(this);
  dock_scroll->setWidget(dock_wid);
  this->QDockWidget::setWidget(dock_scroll);
  setupUi(dock_wid);
  
  connect(this->actionUpdateURI, SIGNAL(triggered()), this, SLOT(onUpdateURIAndDataOpt()));
  connect(this->Columns_edit, SIGNAL(textChanged()), this, SLOT(onUpdateDataOptNames()));
  connect(this->URI_edit, SIGNAL(textChanged(QString)), this, SLOT(onUpdateFieldsAndDataOpt()));
  
  this->URI_edit->setText(tr("file:./data.ssv"));
  
  onUpdateFieldsAndDataOpt();
  onUpdateDataOptNames();
  
};

DataStreamOptWidget::~DataStreamOptWidget() {
  delete static_cast<QScrollArea*>(this->QDockWidget::widget())->widget();
  delete this->QDockWidget::widget();
};


void DataStreamOptWidget::onUpdateURIAndDataOpt() {
  // fill data_opt with the values of the fields.
  switch(mode) {
    case network_only:
      if(this->Format_combo->currentIndex() < 5) {
        QMessageBox::information(this,
                    "Format not allowed!",
                    "Sorry, only network streams are allowed! Please choose a network data-stream type.",
                    QMessageBox::Ok);
        onUpdateFieldsAndDataOpt();
        return;
      };
      break;
    case file_only:
      if(this->Format_combo->currentIndex() >= 5) {
        QMessageBox::information(this,
                    "Format not allowed!",
                    "Sorry, only file streams are allowed! Please choose a file format.",
                    QMessageBox::Ok);
        onUpdateFieldsAndDataOpt();
        return;
      };
    case allow_all:
    default:
      break;
  };
  switch(this->Format_combo->currentIndex()) {
    case 1:
      data_opt.kind = recorder::data_stream_options::binary;
      break;
    case 2:
      data_opt.kind = recorder::data_stream_options::space_separated;
      break;
    case 3:
      data_opt.kind = recorder::data_stream_options::comma_separated;
      break;
    case 4:
      data_opt.kind = recorder::data_stream_options::tab_separated;
      break;
    case 5:
      data_opt.kind = recorder::data_stream_options::tcp_stream;
      break;
    case 6:
      data_opt.kind = recorder::data_stream_options::udp_stream;
      break;
    case 7:
      data_opt.kind = recorder::data_stream_options::raw_udp_stream;
      break;
    default:
      data_opt.kind = recorder::data_stream_options::vector_stream;
      break;
  };
  switch(data_opt.kind) {
    case recorder::data_stream_options::tcp_stream:
    case recorder::data_stream_options::udp_stream:
    case recorder::data_stream_options::raw_udp_stream:
      this->Filename_label->setText(tr("IP Addr. / Host-name:"));
      this->Filename_button->setEnabled(false);
      this->Port_spin->setEnabled(true);
      break;
    case recorder::data_stream_options::binary:
    case recorder::data_stream_options::space_separated:
    case recorder::data_stream_options::comma_separated:
    case recorder::data_stream_options::tab_separated:
    default:
      this->Filename_label->setText(tr("File Name:"));
      this->Filename_button->setEnabled(true);
      this->Port_spin->setEnabled(false);
      break;
  };
  std::string fname = this->Filename_edit->text().toStdString();
  std::size_t portnum = this->Port_spin->value();
  if(this->Port_spin->isEnabled()) {
    std::stringstream ss(fname);
    ss << ":" << portnum;
    fname = ss.str();
  };
  data_opt.file_name = fname;
  
  data_opt.flush_rate = this->Freq_spin->value();
  data_opt.buffer_size = this->Buffer_spin->value();
  data_opt.time_sync_name = this->TimeSync_edit->text().toStdString();
  
  // create URI string to set the URI edit with.
  this->URI_edit->setText(QString::fromStdString(data_opt.get_URI()));
};

void DataStreamOptWidget::onUpdateFieldsAndDataOpt() {
  // fill data_opt with the values of the URI.
  data_opt.set_from_URI(this->URI_edit->text().toStdString());
  // fill the fields with the values from data_opt.
  switch(data_opt.kind) {
    case recorder::data_stream_options::binary:
      this->Format_combo->setCurrentIndex(1);
      break;
    case recorder::data_stream_options::space_separated:
      this->Format_combo->setCurrentIndex(2);
      break;
    case recorder::data_stream_options::comma_separated:
      this->Format_combo->setCurrentIndex(3);
      break;
    case recorder::data_stream_options::tab_separated:
      this->Format_combo->setCurrentIndex(4);
      break;
    case recorder::data_stream_options::tcp_stream:
      this->Format_combo->setCurrentIndex(5);
      break;
    case recorder::data_stream_options::udp_stream:
      this->Format_combo->setCurrentIndex(6);
      break;
    case recorder::data_stream_options::raw_udp_stream:
      this->Format_combo->setCurrentIndex(7);
      break;
    default:
      this->Format_combo->setCurrentIndex(0);
      break;
  };
  switch(data_opt.kind) {
    case recorder::data_stream_options::tcp_stream:
    case recorder::data_stream_options::udp_stream:
    case recorder::data_stream_options::raw_udp_stream:
      this->Filename_label->setText(tr("IP Addr. / Host-name:"));
      this->Filename_button->setEnabled(false);
      this->Port_spin->setEnabled(true);
      break;
    case recorder::data_stream_options::binary:
    case recorder::data_stream_options::space_separated:
    case recorder::data_stream_options::comma_separated:
    case recorder::data_stream_options::tab_separated:
    default:
      this->Filename_label->setText(tr("File Name:"));
      this->Filename_button->setEnabled(true);
      this->Port_spin->setEnabled(false);
      break;
  };
  
  std::string fname = "localhost";
  std::size_t portnum = 17000;
  std::stringstream ss(data_opt.file_name);
  std::getline(ss, fname, ':');
  ss >> portnum;
  this->Filename_edit->setText(QString::fromStdString(fname));
  this->Port_spin->setValue(portnum);
  
  this->Freq_spin->setValue(data_opt.flush_rate);
  this->Buffer_spin->setValue(data_opt.buffer_size);
  this->TimeSync_edit->setText(QString::fromStdString(data_opt.time_sync_name));
  
  if(!this->Columns_edit->isEnabled()) {
    // fill the Columns_edit with the names from data_opt.
    std::stringstream ss;
    for(std::size_t i = 0; i < data_opt.names.size(); ++i) {
      ss << data_opt.names[i] << '\n';
    };
    this->Columns_edit->setPlainText(QString::fromStdString(ss.str()));
  };
  
};

void DataStreamOptWidget::onUpdateDataOptNames() {
  // fill the data_opt names with the names from Columns_edit.
  std::stringstream ss(this->Columns_edit->toPlainText().toStdString());
  data_opt.names.clear();
  for(std::size_t i = 0; i < data_opt.names.size(); ++i) {
    std::string tmp;
    std::getline(ss, tmp, '\n');
    data_opt.add_name(tmp);
  };
};




#if 0

static QString last_used_path;

void TargetPredConfigWidget::saveInertiaTensor() {
  QString fileName = QFileDialog::getSaveFileName(
    this, tr("Save Inertia Information..."), last_used_path,
    tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  onConfigsChanged();
  
  try {
    sat_options.save_mass_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
};

void TargetPredConfigWidget::editInertiaTensor() {
  ot_inertia_win.show();
};

void TargetPredConfigWidget::loadInertiaTensor() {
  QString fileName = QFileDialog::getOpenFileName(
    this, tr("Open Inertia Information..."), last_used_path,
    tr("Target Inertia Information (*.rkx *.xml *.rkb *.pbuf)"));
  
  if( fileName == tr("") )
    return;
  
  last_used_path = QFileInfo(fileName).absolutePath();
  
  try {
    sat_options.load_mass_configs(fileName.toStdString());
  } catch(...) {
    QMessageBox::information(this,
                "File Type Not Supported!",
                "Sorry, this file-type is not supported!",
                QMessageBox::Ok);
    return;
  };
  
  updateConfigs();
};

#endif


};

};





