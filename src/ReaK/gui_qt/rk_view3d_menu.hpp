
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

#ifndef REAK_VIEW3D_MENU_HPP
#define REAK_VIEW3D_MENU_HPP

#include "base/defs.hpp"

#include <QMenu>

#include <string>

class SoQtExaminerViewer;
class SoSwitch;
class SoSeparator;

namespace ReaK { 

namespace geom {
  class oi_scene_graph;
}; 

};


namespace ReaK {
  
  
namespace rkqt {

class View3DMenu : public QMenu {
    Q_OBJECT
  
  public:
    View3DMenu( QWidget * parent = 0, SoQtExaminerViewer* aViewer = NULL );
    ~View3DMenu();
    
  private slots:
    
    void toggleDisplayGroup(bool isChecked);
    
  public:
    
    void setViewer(SoQtExaminerViewer* aViewer);
    SoQtExaminerViewer* getViewer() const { return qtviewer; };
    
    SoSeparator* getRoot() const { return root_sep; };
    
    SoSwitch* getDisplayGroup(const std::string& aGroupName, bool initChecked = true);
    void removeDisplayGroup(const std::string& aGroupName);
    
    shared_ptr<geom::oi_scene_graph> getGeometryGroup(const std::string& aGroupName, bool initChecked = true);
    void removeGeometryGroup(const std::string& aGroupName);
    
  private:
    
    struct display_group {
      QAction* selector;
      SoSwitch* display_switch;
      
      display_group(QAction* aSelector = NULL, SoSwitch* aDisplaySwitch = NULL) : selector(aSelector), display_switch(aDisplaySwitch) { };
    };
    
    struct geometry_group {
      QAction* selector;
      shared_ptr<geom::oi_scene_graph> geom_scene;
      
      geometry_group(QAction* aSelector = NULL, const shared_ptr<geom::oi_scene_graph>& aGeomScene = shared_ptr<geom::oi_scene_graph>()) : 
                     selector(aSelector), geom_scene(aGeomScene) { };
    };
    
    SoQtExaminerViewer* qtviewer;
    SoSeparator* root_sep;
    std::map< std::string, display_group >  display_items;
    std::map< std::string, geometry_group > geometry_items;
    
};


};


};


#endif














