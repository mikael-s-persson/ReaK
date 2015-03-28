
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

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/core/qt/rk_object_tree_widget.hpp>
#include <ReaK/core/qt/rk_prop_editor_widget.hpp>

#include <ReaK/core/serialization/scheme_builder.hpp>

class test_class : public ReaK::named_object {
public:
  double mass;
  int count;
  std::string last_name;

  test_class( const std::string& aName = "" ) : ReaK::named_object() { this->setName( aName ); };

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    ReaK::named_object::save( A, ReaK::named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( mass ) & RK_SERIAL_SAVE_WITH_NAME( count ) & RK_SERIAL_SAVE_WITH_NAME( last_name );
  };
  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    ReaK::named_object::load( A, ReaK::named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( mass ) & RK_SERIAL_LOAD_WITH_NAME( count ) & RK_SERIAL_LOAD_WITH_NAME( last_name );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( test_class, 0xDEADBEEF, 1, "test_class", ReaK::named_object )
};


int main( int argc, char* argv[] ) {
  QApplication app( argc, argv );

  QMainWindow window;
  window.resize( 400, 500 );
  window.move( 100, 100 );
  window.setWindowTitle( "Test" );

  ReaK::shared_ptr< ReaK::serialization::object_graph > obj_g
    = ReaK::shared_ptr< ReaK::serialization::object_graph >( new ReaK::serialization::object_graph() );
  ReaK::serialization::object_node_desc root = add_vertex( *obj_g );

  ReaK::shared_ptr< test_class > mdl = ReaK::shared_ptr< test_class >( new test_class( "model" ) );

  ReaK::serialization::objtree_oarchive out_arc( obj_g, root );
  out_arc << mdl;
  ReaK::serialization::scheme_builder sch_bld;
  sch_bld << mdl;

  ReaK::qt::ObjectTreeWidget objtree( obj_g, root );
  window.addDockWidget( Qt::RightDockWidgetArea, &objtree );

  ReaK::qt::PropEditorWidget propedit( &( objtree.mdl ) );
  window.addDockWidget( Qt::RightDockWidgetArea, &propedit );

  window.show();

  int result = app.exec();

  ReaK::serialization::get_global_schemes().clear();
  return result;
};
