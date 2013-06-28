
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

#include <Inventor/nodes/SoSwitch.h>


namespace ReaK {
  
namespace rkqt {

View3DMenu::View3DMenu( QWidget * parent ) : QMenu(tr("View"), parent) {
  // addAction("Some General Option");
  addSeparator();
};

View3DMenu::~View3DMenu() {
  for(std::map< std::string, display_group >::iterator it = display_items.begin(); it != display_items.end(); ++it)
    for(std::vector< SoSwitch* >::iterator sit = it->second.switches.begin(); sit != it->second.switches.end(); ++sit)
      (*sit)->unref();
};


void View3DMenu::addDisplayGroup(const std::string& aGroupName, bool initChecked) {
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it != display_items.end())
    return;
  
  display_group& dg = display_items[aGroupName];
  dg = display_group(addAction(QString::fromStdString(aGroupName)));
  dg.selector->setCheckable(true);
  dg.selector->setChecked(initChecked);
  
};

void View3DMenu::addSwitchToGroup(const std::string& aGroupName, SoSwitch* aSwitch) {
  if( ! aSwitch )
    return;
  
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it == display_items.end())
    addDisplayGroup(aGroupName, true);
  
  aSwitch->ref();
  display_items[aGroupName].switches.push_back(aSwitch);
  
};

void View3DMenu::removeDisplayGroup(const std::string& aGroupName) {
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it == display_items.end())
    return;
  
  for(std::vector< SoSwitch* >::iterator sit = it->second.switches.begin(); sit != it->second.switches.end(); ++sit)
    (*sit)->unref();
  
  removeAction(it->second.selector);
  display_items.erase(it);
  
};

void View3DMenu::removeSwitchFromGroup(const std::string& aGroupName, SoSwitch* aSwitch) {
  std::map< std::string, display_group >::iterator it = display_items.find(aGroupName);
  if(it == display_items.end())
    return;
  
  std::vector< SoSwitch* >::iterator sit = std::find(it->second.switches.begin(), it->second.switches.end(), aSwitch);
  if( sit == it->second.switches.end() )
    return;
  
  (*sit)->unref();
  it->second.switches.erase(sit);
  
};




};

};



