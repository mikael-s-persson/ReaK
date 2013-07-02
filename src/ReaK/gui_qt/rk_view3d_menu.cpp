
/*
 *    Copyright 2013 Sven Mikael Persson
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

#include "rk_view3d_menu.hpp"

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>


namespace ReaK {
  
namespace rkqt {

View3DMenu::View3DMenu( QWidget * parent, SoSeparator* aRoot ) : QMenu(tr("View"), parent), root_sep(aRoot) {
  // addAction("Some General Option");
  addSeparator();
  
  if(root_sep)
    root_sep->ref();
};

View3DMenu::~View3DMenu() {
  if(root_sep)
    root_sep->unref();
};


void View3DMenu::toggleDisplayGroup(bool isChecked) {
  QAction* snder = static_cast<QAction*>(sender());
  
  std::string snder_name = snder->text().toStdString();
  std::map< std::string, display_group >::iterator it = display_items.find(snder_name);
  if(it == display_items.end())
    return;
  
  it->second.display_switch->whichChild.setValue((isChecked ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  
};

void View3DMenu::setRoot(SoSeparator* aRoot) {
  if((root_sep) && (aRoot)) {
    for(std::map< std::string, display_group >::iterator it = display_items.begin(); it != display_items.end(); ++it) {
      aRoot->addChild(it->second.display_switch);
      root_sep->removeChild(it->second.display_switch);
    };
  };
  if(root_sep)
    root_sep->unref();
  root_sep = aRoot;
  if(root_sep)
    root_sep->ref();
};

SoSwitch* View3DMenu::getDisplayGroup(const std::string& aGroupName, bool initChecked) {
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it != display_items.end())
    return it->second.display_switch;
  if(!root_sep)
    return NULL;
  
  display_group& dg = display_items[aGroupName];
  dg.display_switch = new SoSwitch;
  dg.display_switch->whichChild.setValue((initChecked ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  root_sep->addChild(dg.display_switch);
  dg.selector = addAction(QString::fromStdString(aGroupName));
  dg.selector->setCheckable(true);
  dg.selector->setChecked(initChecked);
  
  connect(dg.selector, SIGNAL(toggled(bool)), this, SLOT(toggleDisplayGroup(bool)));
  
  return dg.display_switch;
};

void View3DMenu::removeDisplayGroup(const std::string& aGroupName) {
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it == display_items.end())
    return;
  
  root_sep->removeChild(it->second.display_switch);
  removeAction(it->second.selector);
  display_items.erase(it);
  
};




};

};



