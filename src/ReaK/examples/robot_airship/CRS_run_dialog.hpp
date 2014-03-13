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

#ifndef REAK_CRS_RUN_DIALOG_HPP
#define REAK_CRS_RUN_DIALOG_HPP

#include "ui_CRS_run_dialog.h"

#include <QDialog>
#include <QTimer>

namespace ReaK {
  
namespace rkqt {


class CRSRunDialogWidget : public QDialog, private Ui::CRSRunDialog {
    Q_OBJECT
  
  public:
    CRSRunDialogWidget(QWidget * parent = NULL, Qt::WindowFlags flags = 0);
    virtual ~CRSRunDialogWidget();
    
  private slots:
    
    void onStartPressed();
    void onStopPressed();
    void onLaunchPressed();
    void onAbortPressed();
    
    void flipLaunchButtonColor();
    
  public slots:
    
    void onInitializationDone();
    void onPlanningDone();
    
    void onLaunchOpportunity();
    
    void onLaunchStarted();
    void onCaptureReached();
    
    void onReset();
    
  signals:
    
    void triggeredStartPlanning(int mode);
    void triggeredStopPlanning();
    
    void triggeredLaunch(int mode);
    
    void triggeredAbort();
    
  private:
    
    QTimer flashing_button_timer;
    
  public:
    
    
};

};

};

#endif




