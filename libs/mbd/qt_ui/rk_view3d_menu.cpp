
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

#include <ReaK/mbd/qt/rk_view3d_menu.hpp>

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>


namespace ReaK {

namespace qt {

View3DMenu::View3DMenu( QWidget* parent, SoQtExaminerViewer* aViewer )
    : QMenu( tr( "View" ), parent ), qtviewer( aViewer ), root_sep( nullptr ) {
  // addAction("Some General Option");
  addSeparator();

  if( qtviewer ) {
    root_sep = new SoSeparator;
    root_sep->ref();
    qtviewer->setSceneGraph( root_sep );
    qtviewer->show();
  };
};

View3DMenu::~View3DMenu() {
  if( root_sep )
    root_sep->unref();
  delete qtviewer;
};


void View3DMenu::toggleDisplayGroup( bool isChecked ) {
  QAction* snder = static_cast< QAction* >( sender() );

  std::string snder_name = snder->text().toStdString();
  std::map< std::string, display_group >::iterator itd = display_items.find( snder_name );
  if( itd != display_items.end() ) {
    itd->second.display_switch->whichChild.setValue( ( isChecked ? SO_SWITCH_ALL : SO_SWITCH_NONE ) );
    return;
  };

  std::map< std::string, geometry_group >::iterator itg = geometry_items.find( snder_name );
  if( itg != geometry_items.end() ) {
    itg->second.geom_scene->setVisibility( isChecked );
    return;
  };
};


void View3DMenu::setViewer( SoQtExaminerViewer* aViewer ) {
  SoSeparator* newRoot = nullptr;
  if( aViewer ) {
    newRoot = new SoSeparator;
    newRoot->ref();
    aViewer->setSceneGraph( newRoot );
    aViewer->show();
  };
  if( ( root_sep ) && ( newRoot ) ) {
    for( std::map< std::string, display_group >::iterator it = display_items.begin(); it != display_items.end();
         ++it ) {
      newRoot->addChild( it->second.display_switch );
      root_sep->removeChild( it->second.display_switch );
    };
    for( std::map< std::string, geometry_group >::iterator it = geometry_items.begin(); it != geometry_items.end();
         ++it ) {
      newRoot->addChild( it->second.geom_scene->getSceneGraph() );
      root_sep->removeChild( it->second.geom_scene->getSceneGraph() );
    };
  };
  if( root_sep )
    root_sep->unref();
  delete qtviewer;
  qtviewer = aViewer;
  root_sep = newRoot; // already been ref'ed (above).
};


SoSwitch* View3DMenu::getDisplayGroup( const std::string& aGroupName, bool initChecked ) {
  std::map< std::string, display_group >::iterator it = display_items.find( aGroupName );
  if( it != display_items.end() )
    return it->second.display_switch;
  if( !root_sep )
    return nullptr;

  display_group& dg = display_items[aGroupName];
  dg.display_switch = new SoSwitch;
  dg.display_switch->whichChild.setValue( ( initChecked ? SO_SWITCH_ALL : SO_SWITCH_NONE ) );
  root_sep->addChild( dg.display_switch );
  dg.selector = addAction( QString::fromStdString( aGroupName ) );
  dg.selector->setCheckable( true );
  dg.selector->setChecked( initChecked );

  connect( dg.selector, SIGNAL( toggled(bool)), this, SLOT( toggleDisplayGroup(bool)) );

  return dg.display_switch;
};

void View3DMenu::removeDisplayGroup( const std::string& aGroupName ) {
  std::map< std::string, display_group >::iterator it = display_items.find( aGroupName );
  if( it == display_items.end() )
    return;

  if( root_sep )
    root_sep->removeChild( it->second.display_switch );
  removeAction( it->second.selector );
  display_items.erase( it );
};


shared_ptr< geom::oi_scene_graph > View3DMenu::getGeometryGroup( const std::string& aGroupName, bool initChecked ) {
  std::map< std::string, geometry_group >::iterator it = geometry_items.find( aGroupName );
  if( it != geometry_items.end() )
    return it->second.geom_scene;
  if( !root_sep )
    return nullptr;

  geometry_group& gg = geometry_items[aGroupName];
  gg.geom_scene = shared_ptr< geom::oi_scene_graph >( new geom::oi_scene_graph );
  gg.geom_scene->enableAnchorUpdates();
  gg.geom_scene->setVisibility( initChecked );
  root_sep->addChild( gg.geom_scene->getSceneGraph() );
  gg.selector = addAction( QString::fromStdString( aGroupName ) );
  gg.selector->setCheckable( true );
  gg.selector->setChecked( initChecked );

  connect( gg.selector, SIGNAL( toggled(bool)), this, SLOT( toggleDisplayGroup(bool)) );

  return gg.geom_scene;
};

void View3DMenu::removeGeometryGroup( const std::string& aGroupName ) {
  std::map< std::string, geometry_group >::iterator it = geometry_items.find( aGroupName );
  if( it == geometry_items.end() )
    return;

  if( root_sep )
    root_sep->removeChild( it->second.geom_scene->getSceneGraph() );
  removeAction( it->second.selector );
  geometry_items.erase( it );
};
};
};
