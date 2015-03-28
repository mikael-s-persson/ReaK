
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

#include <ReaK/core/qt/data_stream_opt_widget.hpp>

#include "ui_data_stream_opt.h"

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

namespace qt {


DataStreamOptWidget::DataStreamOptWidget( recorder::data_stream_options aDataOpt, QWidget* parent,
                                          Qt::WindowFlags flags )
    : QDockWidget( tr( "Data Stream" ), parent, flags ), ui( new Ui::DataStreamOpt() ), data_opt( aDataOpt ),
      mode( allow_all ) {
  QScrollArea* dock_scroll = new QScrollArea( this );
  dock_scroll->setWidgetResizable( true );
  QWidget* dock_wid = new QWidget( this );
  dock_scroll->setWidget( dock_wid );
  this->QDockWidget::setWidget( dock_scroll );
  ui->setupUi( dock_wid );

  connect( ui->actionUpdateURI, SIGNAL( triggered() ), this, SLOT( onUpdateURIAndDataOpt() ) );
  connect( ui->Columns_edit, SIGNAL( textChanged() ), this, SLOT( onUpdateDataOptNames() ) );
  connect( ui->URI_edit, SIGNAL( textChanged( QString ) ), this, SLOT( onUpdateFieldsAndDataOpt() ) );

  ui->URI_edit->setText( tr( "file:./data.ssv" ) );

  onUpdateFieldsAndDataOpt();
  onUpdateDataOptNames();
};

DataStreamOptWidget::~DataStreamOptWidget() {
  delete static_cast< QScrollArea* >( this->QDockWidget::widget() )->widget();
  delete this->QDockWidget::widget();
  delete ui;
};


void DataStreamOptWidget::onUpdateURIAndDataOpt() {
  // fill data_opt with the values of the fields.
  switch( mode ) {
    case network_only:
      if( ui->Format_combo->currentIndex() < 5 ) {
        QMessageBox::information( this, "Format not allowed!",
                                  "Sorry, only network streams are allowed! Please choose a network data-stream type.",
                                  QMessageBox::Ok );
        onUpdateFieldsAndDataOpt();
        return;
      };
      break;
    case file_only:
      if( ui->Format_combo->currentIndex() >= 5 ) {
        QMessageBox::information( this, "Format not allowed!",
                                  "Sorry, only file streams are allowed! Please choose a file format.",
                                  QMessageBox::Ok );
        onUpdateFieldsAndDataOpt();
        return;
      };
    case allow_all:
    default:
      break;
  };
  switch( ui->Format_combo->currentIndex() ) {
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
  switch( data_opt.kind ) {
    case recorder::data_stream_options::tcp_stream:
    case recorder::data_stream_options::udp_stream:
    case recorder::data_stream_options::raw_udp_stream:
      ui->Filename_label->setText( tr( "IP Addr. / Host-name:" ) );
      ui->Filename_button->setEnabled( false );
      ui->Port_spin->setEnabled( true );
      break;
    case recorder::data_stream_options::binary:
    case recorder::data_stream_options::space_separated:
    case recorder::data_stream_options::comma_separated:
    case recorder::data_stream_options::tab_separated:
    default:
      ui->Filename_label->setText( tr( "File Name:" ) );
      ui->Filename_button->setEnabled( true );
      ui->Port_spin->setEnabled( false );
      break;
  };
  std::string fname = ui->Filename_edit->text().toStdString();
  std::size_t portnum = ui->Port_spin->value();
  if( ui->Port_spin->isEnabled() ) {
    std::stringstream ss( fname );
    ss << ":" << portnum;
    fname = ss.str();
  };
  data_opt.file_name = fname;

  data_opt.flush_rate = ui->Freq_spin->value();
  data_opt.buffer_size = ui->Buffer_spin->value();
  data_opt.time_sync_name = ui->TimeSync_edit->text().toStdString();

  // create URI string to set the URI edit with.
  ui->URI_edit->setText( QString::fromStdString( data_opt.get_URI() ) );
};

void DataStreamOptWidget::onUpdateFieldsAndDataOpt() {
  // fill data_opt with the values of the URI.
  data_opt.set_from_URI( ui->URI_edit->text().toStdString() );
  // fill the fields with the values from data_opt.
  switch( data_opt.kind ) {
    case recorder::data_stream_options::binary:
      ui->Format_combo->setCurrentIndex( 1 );
      break;
    case recorder::data_stream_options::space_separated:
      ui->Format_combo->setCurrentIndex( 2 );
      break;
    case recorder::data_stream_options::comma_separated:
      ui->Format_combo->setCurrentIndex( 3 );
      break;
    case recorder::data_stream_options::tab_separated:
      ui->Format_combo->setCurrentIndex( 4 );
      break;
    case recorder::data_stream_options::tcp_stream:
      ui->Format_combo->setCurrentIndex( 5 );
      break;
    case recorder::data_stream_options::udp_stream:
      ui->Format_combo->setCurrentIndex( 6 );
      break;
    case recorder::data_stream_options::raw_udp_stream:
      ui->Format_combo->setCurrentIndex( 7 );
      break;
    default:
      ui->Format_combo->setCurrentIndex( 0 );
      break;
  };
  switch( data_opt.kind ) {
    case recorder::data_stream_options::tcp_stream:
    case recorder::data_stream_options::udp_stream:
    case recorder::data_stream_options::raw_udp_stream:
      ui->Filename_label->setText( tr( "IP Addr. / Host-name:" ) );
      ui->Filename_button->setEnabled( false );
      ui->Port_spin->setEnabled( true );
      break;
    case recorder::data_stream_options::binary:
    case recorder::data_stream_options::space_separated:
    case recorder::data_stream_options::comma_separated:
    case recorder::data_stream_options::tab_separated:
    default:
      ui->Filename_label->setText( tr( "File Name:" ) );
      ui->Filename_button->setEnabled( true );
      ui->Port_spin->setEnabled( false );
      break;
  };

  std::string fname = "localhost";
  std::size_t portnum = 17000;
  std::stringstream ss( data_opt.file_name );
  std::getline( ss, fname, ':' );
  ss >> portnum;
  ui->Filename_edit->setText( QString::fromStdString( fname ) );
  ui->Port_spin->setValue( portnum );

  ui->Freq_spin->setValue( data_opt.flush_rate );
  ui->Buffer_spin->setValue( data_opt.buffer_size );
  ui->TimeSync_edit->setText( QString::fromStdString( data_opt.time_sync_name ) );

  if( !ui->Columns_edit->isEnabled() ) {
    // fill the Columns_edit with the names from data_opt.
    std::stringstream ss;
    for( std::size_t i = 0; i < data_opt.names.size(); ++i ) {
      ss << data_opt.names[i] << '\n';
    };
    ui->Columns_edit->setPlainText( QString::fromStdString( ss.str() ) );
  };
};

void DataStreamOptWidget::onUpdateDataOptNames() {
  // fill the data_opt names with the names from Columns_edit.
  std::stringstream ss( ui->Columns_edit->toPlainText().toStdString() );
  data_opt.names.clear();
  for( std::size_t i = 0; i < data_opt.names.size(); ++i ) {
    std::string tmp;
    std::getline( ss, tmp, '\n' );
    data_opt.add_name( tmp );
  };
};

void DataStreamOptWidget::setEnabledNameEdits( bool aEnabled ) {
  ui->Columns_edit->setEnabled( aEnabled );
  onUpdateFieldsAndDataOpt();
  onUpdateDataOptNames();
};

void DataStreamOptWidget::setEnabledFreqSetting( bool aEnabled ) { ui->Freq_spin->setEnabled( aEnabled ); };

void DataStreamOptWidget::setEnabledBufferSetting( bool aEnabled ) { ui->Buffer_spin->setEnabled( aEnabled ); };

void DataStreamOptWidget::setEnabledTimeSyncSetting( bool aEnabled ) { ui->TimeSync_edit->setEnabled( aEnabled ); };


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
