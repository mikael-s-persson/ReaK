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

#ifndef REAK_CORE_QT_RK_OBJECT_TREE_WIDGET_H_
#define REAK_CORE_QT_RK_OBJECT_TREE_WIDGET_H_

#include "ReaK/core/qt/objtree_qtmodel.h"

#include <QDockWidget>

namespace Ui {
class RKObjectTreeWidget;
};

namespace ReaK::qt {

class ObjectTreeWidget : public QDockWidget {
  Q_OBJECT

 public:
  ObjectTreeWidget(
      const std::shared_ptr<serialization::object_graph>& aObjGraph,
      serialization::object_node_desc aRoot = 0, QWidget* parent = nullptr,
      Qt::WindowFlags flags = 0);
  virtual ~ObjectTreeWidget();

 private slots:

  void recordPreviousSelection();
  void restorePreviousSelection();

 private:
  Ui::RKObjectTreeWidget* ui;

 public:
  ObjTreeQtModel mdl;
  QString prev_selection;
};

}  // namespace ReaK::qt

#endif  // REAK_CORE_QT_RK_OBJECT_TREE_WIDGET_H_
