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

#ifndef REAK_PROP_EDITOR_WIDGET_HPP
#define REAK_PROP_EDITOR_WIDGET_HPP

#include "objtree_qtmodel.hpp"
#include "obj_properties_qtmodel.hpp"

#include <QDockWidget>

namespace Ui {
class RKPropEditorWidget;
};

namespace ReaK {

namespace qt {

class PropEditorWidget : public QDockWidget {
  Q_OBJECT

public:
  PropEditorWidget( ObjTreeQtModel* aObjTreeMdl, QWidget* parent = NULL, Qt::WindowFlags flags = 0 );
  virtual ~PropEditorWidget();

private slots:

  void xmlSrcChanged();
  void onTextChanged();
  void applyButtonClick();
  void cancelButtonClick();
  void objNameChanged( const std::string& newName );

signals:

  void editedXMLSrc( const std::string& newSrc );

private:
  Ui::RKPropEditorWidget* ui;

public:
  ObjPropertiesQtModel mdl;
  ObjPropertiesQtDelegate delegate;
};
};
};

#endif
