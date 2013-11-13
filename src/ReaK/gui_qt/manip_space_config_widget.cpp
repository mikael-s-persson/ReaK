
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

#include "manip_space_config_widget.hpp"

#include <QDockWidget>


namespace ReaK {
  
namespace rkqt {


ManipSpaceConfigWidget::ManipSpaceConfigWidget(QWidget * parent, Qt::WindowFlags flags) :
                                               QDockWidget(tr("Space"), parent, flags),
                                               Ui::ManipSpaceConfig()
{
  this->QDockWidget::setWidget(new QWidget(this));
  setupUi(this->QDockWidget::widget());
  
  connect(this->actionValuesChanged, SIGNAL(triggered()), this, SLOT(onUpdateAvailableOptions()));
  
  onUpdateAvailableOptions();
  
};

ManipSpaceConfigWidget::~ManipSpaceConfigWidget() { 
  delete this->QDockWidget::widget();
};


void ManipSpaceConfigWidget::onUpdateAvailableOptions() {
  
  space_order = this->order_selection->currentIndex();
  interp_id = this->interp_selection->currentIndex();
  min_travel = this->min_interval_spinbox->value();
  max_travel = this->max_interval_spinbox->value();
  
  is_temporal = this->temporal_space_check->isChecked();
  is_rate_limited = this->rate_limited_check->isChecked();
  
  output_space_order = this->output_space_selection->currentIndex();
  
};




};

};














