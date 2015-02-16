
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

#include <ReaK/core/qt/rk_prop_editor_widget.hpp>

#include <QDockWidget>
#include <QTreeView>

#include "ui_rk_prop_editor.h"

namespace ReaK {
  
namespace qt {


PropEditorWidget::PropEditorWidget(ObjTreeQtModel* aObjTreeMdl, QWidget * parent, Qt::WindowFlags flags) :
                                   QDockWidget(parent, flags),
                                   ui(new Ui::RKPropEditorWidget()),
                                   mdl(aObjTreeMdl->get_object_graph(),aObjTreeMdl->get_root_node()),
                                   delegate(&mdl)
{
  ui->setupUi(this);
  ui->tableView->setModel(&mdl);
  ui->tableView->setRootIndex(QModelIndex());
  ui->tableView->setItemDelegate(&delegate);
  ui->tableTab->setAttribute(Qt::WA_AlwaysShowToolTips, true);
  ui->sourceTab->setAttribute(Qt::WA_AlwaysShowToolTips, true);
  
  connect(aObjTreeMdl, SIGNAL(objectNodeSelected(serialization::object_node_desc)), &mdl, SLOT(selectObject(serialization::object_node_desc)));
  connect(&mdl, SIGNAL(objectTreeChanged()), aObjTreeMdl, SLOT(treeChanged()));

  connect(&mdl, SIGNAL(objectNameChanged(std::string)), this, SLOT(objNameChanged(std::string)));
  connect(&mdl, SIGNAL(sourceDataChanged()), this, SLOT(xmlSrcChanged()));
  connect(this, SIGNAL(editedXMLSrc(std::string)), &mdl, SLOT(sourceDataEdited(std::string)));
  
  connect(ui->applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClick()));
  connect(ui->cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClick()));
  
  connect(ui->textEdit, SIGNAL(textChanged()), this, SLOT(onTextChanged()));
};

PropEditorWidget::~PropEditorWidget() { 
  delete ui;
};

void PropEditorWidget::xmlSrcChanged() {
  ui->textEdit->setText(QString::fromStdString(mdl.get_current_fields().get_complete_src()));
  ui->applyButton->setEnabled(false);
  ui->cancelButton->setEnabled(false);
};

void PropEditorWidget::onTextChanged() {
  ui->applyButton->setEnabled(true);
  ui->cancelButton->setEnabled(true);
};

void PropEditorWidget::applyButtonClick() {
  emit editedXMLSrc(ui->textEdit->toPlainText().toStdString());
  ui->applyButton->setEnabled(false);
  ui->cancelButton->setEnabled(false);
};

void PropEditorWidget::cancelButtonClick() {
  ui->textEdit->setText(QString::fromStdString(mdl.get_current_fields().get_complete_src()));
  ui->applyButton->setEnabled(false);
  ui->cancelButton->setEnabled(false);
};

void PropEditorWidget::objNameChanged(const std::string& newName) {
  if(newName != "")
    this->setWindowTitle(QString("Property Editor - ") + QString::fromStdString(newName));
  else
    this->setWindowTitle(QString("Property Editor"));
};


};

};


