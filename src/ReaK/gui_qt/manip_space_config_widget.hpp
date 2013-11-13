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

#ifndef REAK_MANIP_SPACE_CONFIG_WIDGET_HPP
#define REAK_MANIP_SPACE_CONFIG_WIDGET_HPP

#include "ui_manip_space_config.h"
#include <QDockWidget>

namespace ReaK {
  
namespace rkqt {

class ManipSpaceConfigWidget : public QDockWidget, private Ui::ManipSpaceConfig {
    Q_OBJECT
  
  public:
    ManipSpaceConfigWidget(QWidget * parent = NULL, Qt::WindowFlags flags = 0);
    virtual ~ManipSpaceConfigWidget();
    
  private slots:
    
    void onUpdateAvailableOptions();
    
  public:
    
    int space_order;
    int interp_id;
    double min_travel;
    double max_travel;
    bool is_temporal;
    bool is_rate_limited;
    
    int output_space_order;
    
};

};

};

#endif














