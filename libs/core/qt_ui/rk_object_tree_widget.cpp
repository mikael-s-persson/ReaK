
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

#include "ReaK/core/qt/rk_object_tree_widget.hpp"

#include "ui_rk_object_tree.h"

#include <QDockWidget>
#include <QTreeView>

namespace ReaK {

namespace qt {

ObjectTreeWidget::ObjectTreeWidget(
    const std::shared_ptr<serialization::object_graph>& aObjGraph,
    serialization::object_node_desc aRoot, QWidget* parent,
    Qt::WindowFlags flags)
    : QDockWidget(parent, flags),
      ui(new Ui::RKObjectTreeWidget()),
      mdl(aObjGraph, aRoot) {
  ui->setupUi(this);
  ui->treeView->setModel(&mdl);
  ui->treeView->setRootIndex(QModelIndex());

  connect(ui->treeView, SIGNAL(clicked(QModelIndex)), &mdl,
          SLOT(itemSelected(QModelIndex)));
  connect(&mdl, SIGNAL(aboutToResetModel()), this,
          SLOT(recordPreviousSelection()));
  connect(&mdl, SIGNAL(justAfterModelReset()), this,
          SLOT(restorePreviousSelection()));
};

ObjectTreeWidget::~ObjectTreeWidget() {
  delete ui;
};

void ObjectTreeWidget::recordPreviousSelection() {
  prev_selection =
      mdl.data(ui->treeView->currentIndex(), Qt::DisplayRole).toString();
};

void ObjectTreeWidget::restorePreviousSelection() {
  ui->treeView->keyboardSearch(prev_selection);
};
};  // namespace qt
};  // namespace ReaK
