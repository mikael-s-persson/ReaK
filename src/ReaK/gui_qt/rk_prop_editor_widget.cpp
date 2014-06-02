
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

#include <ReaK/gui_qt/rk_prop_editor_widget.hpp>

#include <QDockWidget>
#include <QTreeView>

namespace ReaK {
  
namespace rkqt {


PropEditorWidget::PropEditorWidget(ObjTreeQtModel* aObjTreeMdl, QWidget * parent, Qt::WindowFlags flags) :
                                   QDockWidget(parent, flags),
                                   Ui::RKPropEditorWidget(),
                                   mdl(aObjTreeMdl->get_object_graph(),aObjTreeMdl->get_root_node()),
                                   delegate(&mdl)
{
  setupUi(this);
  this->tableView->setModel(&mdl);
  this->tableView->setRootIndex(QModelIndex());
  this->tableView->setItemDelegate(&delegate);
  this->tableTab->setAttribute(Qt::WA_AlwaysShowToolTips, true);
  this->sourceTab->setAttribute(Qt::WA_AlwaysShowToolTips, true);
  
  connect(aObjTreeMdl, SIGNAL(objectNodeSelected(serialization::object_node_desc)), &mdl, SLOT(selectObject(serialization::object_node_desc)));
  connect(&mdl, SIGNAL(objectTreeChanged()), aObjTreeMdl, SLOT(treeChanged()));

  connect(&mdl, SIGNAL(objectNameChanged(std::string)), this, SLOT(objNameChanged(std::string)));
  connect(&mdl, SIGNAL(sourceDataChanged()), this, SLOT(xmlSrcChanged()));
  connect(this, SIGNAL(editedXMLSrc(std::string)), &mdl, SLOT(sourceDataEdited(std::string)));
  
  connect(this->applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClick()));
  connect(this->cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClick()));
  
  connect(this->textEdit, SIGNAL(textChanged()), this, SLOT(onTextChanged()));
};

void PropEditorWidget::xmlSrcChanged() {
  this->textEdit->setText(QString::fromStdString(mdl.get_current_fields().get_complete_src()));
  this->applyButton->setEnabled(false);
  this->cancelButton->setEnabled(false);
};

void PropEditorWidget::onTextChanged() {
  this->applyButton->setEnabled(true);
  this->cancelButton->setEnabled(true);
};

void PropEditorWidget::applyButtonClick() {
  emit editedXMLSrc(this->textEdit->toPlainText().toStdString());
  this->applyButton->setEnabled(false);
  this->cancelButton->setEnabled(false);
};

void PropEditorWidget::cancelButtonClick() {
  this->textEdit->setText(QString::fromStdString(mdl.get_current_fields().get_complete_src()));
  this->applyButton->setEnabled(false);
  this->cancelButton->setEnabled(false);
};

void PropEditorWidget::objNameChanged(const std::string& newName) {
  if(newName != "")
    this->setWindowTitle(QString("Property Editor - ") + QString::fromStdString(newName));
  else
    this->setWindowTitle(QString("Property Editor"));
};


};

};














