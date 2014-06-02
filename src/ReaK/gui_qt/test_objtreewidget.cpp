
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

#include <QApplication>
#include <QMainWindow>

#include "rk_object_tree_widget.hpp"
#include "rk_prop_editor_widget.hpp"

#include <ReaK/ctrl/kte_models/manip_3R_arm.hpp>
#include <ReaK/core/serialization/scheme_builder.hpp>


int main(int argc, char *argv[]) {
  QApplication app(argc, argv); 
  
  QMainWindow window;
  window.resize(400, 500);
  window.move(100, 100);  
  window.setWindowTitle("Test");
  
  ReaK::shared_ptr< ReaK::serialization::object_graph > obj_g = ReaK::shared_ptr< ReaK::serialization::object_graph >(new ReaK::serialization::object_graph());
  ReaK::serialization::object_node_desc root = add_vertex(*obj_g);
  
  ReaK::shared_ptr< ReaK::kte::manip_3R_2D_kinematics > mdl = ReaK::shared_ptr< ReaK::kte::manip_3R_2D_kinematics >(new ReaK::kte::manip_3R_2D_kinematics("3R_model"));
  
  ReaK::serialization::objtree_oarchive out_arc(obj_g, root);
  out_arc << mdl;
  ReaK::serialization::scheme_builder sch_bld;
  sch_bld << mdl;
  
  ReaK::rkqt::ObjectTreeWidget objtree(obj_g, root);
  window.addDockWidget(Qt::RightDockWidgetArea, &objtree);
  
  ReaK::rkqt::PropEditorWidget propedit(&(objtree.mdl));
  window.addDockWidget(Qt::RightDockWidgetArea, &propedit);
  
  window.show();
  
  int result = app.exec();
  
  ReaK::serialization::get_global_schemes().clear();
  return result;
};













