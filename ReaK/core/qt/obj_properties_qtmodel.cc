
/*
 *    Copyright 2012 Sven Mikael Persson
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

#include "ReaK/core/qt/obj_properties_qtmodel.h"

#include "ReaK/core/base/named_object.h"

#include <QMessageBox>
#include "ReaK/core/serialization/scheme_builder.h"
#include "ReaK/core/serialization/type_schemes.h"

#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSpinBox>

#include "ui_rk_class_select.h"

namespace ReaK {

namespace qt {

//     std::shared_ptr< serialization::object_graph > obj_graph;
//     serialization::object_node_desc current_node;

ObjPropertiesQtModel::ObjPropertiesQtModel(
    const std::shared_ptr<serialization::object_graph>& aObjGraph,
    serialization::object_node_desc aCurrentNode)
    : obj_graph_editor(aObjGraph),
      current_node(aCurrentNode),
      current_fields(obj_graph_editor.create_field_editor(aCurrentNode)){};

ObjPropertiesQtModel::~ObjPropertiesQtModel(){};

int ObjPropertiesQtModel::rowCount(const QModelIndex&) const {
  // compute the total row-count:
  return current_fields.get_total_field_count();
};

int ObjPropertiesQtModel::columnCount(const QModelIndex&) const {
  return 2;
};

QVariant ObjPropertiesQtModel::data(const QModelIndex& index, int role) const {
  if (!index.isValid())
    return QVariant();

  std::size_t cur_count = current_fields.get_total_field_count();

  if (index.row() >= int(cur_count) || index.row() < 0)
    return QVariant();

  if (role == Qt::DisplayRole) {
    if (index.column() == 0)
      return QString(current_fields.get_field(index.row()).first.c_str());
    else if (index.column() == 1)
      return QString(current_fields.get_field_value(index.row()).c_str());
  };
  return QVariant();
};

QVariant ObjPropertiesQtModel::headerData(int section,
                                          Qt::Orientation orientation,
                                          int role) const {
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
    if (section == 0)
      return QVariant(QString("Property"));
    else if (section == 1)
      return QVariant(QString("Value"));
  };
  return QVariant();
};

Qt::ItemFlags ObjPropertiesQtModel::flags(const QModelIndex& index) const {
  if (!index.isValid())
    return 0;

  if (index.column() == 0)
    return Qt::ItemIsEnabled;
  else if (index.column() == 1)
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
  else
    return 0;
};

bool ObjPropertiesQtModel::setData(const QModelIndex& index,
                                   const QVariant& value, int role) {
  if ((!index.isValid()) || (role != Qt::EditRole) || (index.column() != 1))
    return false;

  int row = index.row();

  std::pair<std::string, std::shared_ptr<serialization::type_scheme>> fld =
      current_fields.get_field(row);
  rtti::so_type* p_type = fld.second->getObjectType();
  QVariant tmp = value;

  bool convertWorked = false;
  if (p_type == serialization::primitive_scheme<int>::getStaticObjectType()) {
    convertWorked = tmp.convert(QVariant::Int);
  } else if ((p_type == serialization::primitive_scheme<
                            unsigned char>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            unsigned int>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            long unsigned int>::getStaticObjectType())) {
    convertWorked = tmp.convert(QVariant::UInt);
  } else if ((p_type ==
              serialization::primitive_scheme<float>::getStaticObjectType()) ||
             (p_type ==
              serialization::primitive_scheme<double>::getStaticObjectType())) {
    convertWorked = tmp.convert(QVariant::Double);
  } else if (p_type ==
             serialization::primitive_scheme<char>::getStaticObjectType()) {
    convertWorked = tmp.convert(QVariant::Char);
  } else if (p_type ==
             serialization::primitive_scheme<bool>::getStaticObjectType()) {
    convertWorked = tmp.convert(QVariant::Bool);
  } else if (p_type == serialization::primitive_scheme<
                           std::string>::getStaticObjectType()) {
    convertWorked = tmp.convert(QVariant::String);
  } else if (p_type ==
             serialization::serializable_ptr_scheme::getStaticObjectType()) {
    convertWorked = tmp.convert(QVariant::String);
  };
  if (!convertWorked)
    return false;

  current_fields.set_field_value(row, tmp.toString().toStdString());

  emit(dataChanged(index, index));
  emit sourceDataChanged();
  if (p_type ==
      serialization::primitive_scheme<std::string>::getStaticObjectType())
    emit objectNameChanged(current_fields.get_object_name());
  if (p_type == serialization::serializable_ptr_scheme::getStaticObjectType())
    emit objectTreeChanged();

  return true;
};

void ObjPropertiesQtModel::setDataNewObject(
    const QModelIndex& index, const std::shared_ptr<shared_object>& aNewPtr) {
  if ((!index.isValid()) || (index.column() != 1))
    return;

  current_fields.set_field_newptr(index.row(), aNewPtr);

  emit dataChanged(index, index);
  emit sourceDataChanged();
  emit objectTreeChanged();
};

void ObjPropertiesQtModel::selectObject(serialization::object_node_desc u) {
  beginResetModel();
  current_node = u;
  current_fields = obj_graph_editor.create_field_editor(u);
  endResetModel();

  emit objectNameChanged(current_fields.get_object_name());
  emit sourceDataChanged();
};

void ObjPropertiesQtModel::sourceDataEdited(const std::string& newSrc) {
  current_fields.set_complete_src(newSrc);
  emit objectNameChanged(current_fields.get_object_name());
  emit sourceDataChanged();
  emit objectTreeChanged();
};

ObjPropertiesQtDelegate::ObjPropertiesQtDelegate(
    ObjPropertiesQtModel* aParentModel)
    : QStyledItemDelegate(), parentModel(aParentModel){};

ObjPropertiesQtDelegate::~ObjPropertiesQtDelegate(){};

QWidget* ObjPropertiesQtDelegate::createEditor(
    QWidget* parent, const QStyleOptionViewItem& option,
    const QModelIndex& index) const {
  if ((!index.isValid()) || (index.column() != 1))
    return QStyledItemDelegate::createEditor(parent, option, index);

  std::pair<std::string, std::shared_ptr<serialization::type_scheme>> fld =
      parentModel->current_fields.get_field(index.row());
  rtti::so_type* p_type = fld.second->getObjectType();

  if (p_type == serialization::primitive_scheme<bool>::getStaticObjectType()) {
    QComboBox* box = new QComboBox(parent);
    box->addItem(QVariant(QBool(false)).toString(), QVariant(QBool(false)));
    box->addItem(QVariant(QBool(true)).toString(), QVariant(QBool(true)));
    return box;
  } else if (p_type ==
             serialization::primitive_scheme<int>::getStaticObjectType()) {
    QSpinBox* box = new QSpinBox(parent);
    box->setMinimum(std::numeric_limits<int>::min());
    box->setMaximum(std::numeric_limits<int>::max());
    return box;
  } else if ((p_type == serialization::primitive_scheme<
                            unsigned char>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            unsigned int>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            long unsigned int>::getStaticObjectType())) {
    QSpinBox* box = new QSpinBox(parent);
    box->setMinimum(0);
    box->setMaximum(std::numeric_limits<int>::max());
    return box;
  } else if ((p_type ==
              serialization::primitive_scheme<float>::getStaticObjectType()) ||
             (p_type ==
              serialization::primitive_scheme<double>::getStaticObjectType())) {
    QDoubleSpinBox* box = new QDoubleSpinBox(parent);
    box->setMinimum(-std::numeric_limits<double>::infinity());
    box->setMaximum(std::numeric_limits<double>::infinity());
    box->setSingleStep(0.01);
    return box;
  } else if (p_type ==
             serialization::primitive_scheme<char>::getStaticObjectType()) {
    QSpinBox* box = new QSpinBox(parent);
    box->setMinimum(std::numeric_limits<char>::min());
    box->setMaximum(std::numeric_limits<char>::max());
    return box;
  } else if (p_type == serialization::primitive_scheme<
                           std::string>::getStaticObjectType()) {
    return QStyledItemDelegate::createEditor(parent, option, index);
  } else if (p_type ==
             serialization::serializable_ptr_scheme::getStaticObjectType()) {
    QComboBox* box = new QComboBox(parent);
    rtti::so_type* type_ptr =
        rtti::so_type_repo::getInstance().findType(fld.second->get_type_ID());
    if (!type_ptr)
      return box;
    std::vector<std::string> v_tmp =
        parentModel->obj_graph_editor.get_objects_derived_from(type_ptr);
    for (std::size_t i = 0; i < v_tmp.size(); ++i)
      box->addItem(QString::fromStdString(v_tmp[i]));
    box->addItem("New object...");
    return box;
  };

  return QStyledItemDelegate::createEditor(parent, option, index);
};

void ObjPropertiesQtDelegate::setEditorData(QWidget* editor,
                                            const QModelIndex& index) const {
  if ((!index.isValid()) || (index.column() != 1))
    return QStyledItemDelegate::setEditorData(editor, index);

  std::pair<std::string, std::shared_ptr<serialization::type_scheme>> fld =
      parentModel->current_fields.get_field(index.row());
  rtti::so_type* p_type = fld.second->getObjectType();

  if (p_type == serialization::primitive_scheme<bool>::getStaticObjectType()) {
    QComboBox* box = static_cast<QComboBox*>(editor);
    bool value = parentModel->data(index, Qt::DisplayRole).toBool();
    box->setCurrentIndex(box->findData(QVariant(value)));
    return;
  } else if (p_type ==
             serialization::primitive_scheme<int>::getStaticObjectType()) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->setValue(parentModel->data(index, Qt::DisplayRole).toInt());
    return;
  } else if ((p_type == serialization::primitive_scheme<
                            unsigned char>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            unsigned int>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            long unsigned int>::getStaticObjectType())) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->setValue(parentModel->data(index, Qt::DisplayRole).toInt());
    return;
  } else if ((p_type ==
              serialization::primitive_scheme<float>::getStaticObjectType()) ||
             (p_type ==
              serialization::primitive_scheme<double>::getStaticObjectType())) {
    QDoubleSpinBox* box = static_cast<QDoubleSpinBox*>(editor);
    box->setValue(parentModel->data(index, Qt::DisplayRole).toDouble());
    return;
  } else if (p_type ==
             serialization::primitive_scheme<char>::getStaticObjectType()) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->setValue(parentModel->data(index, Qt::DisplayRole).toInt());
    return;
  } else if (p_type == serialization::primitive_scheme<
                           std::string>::getStaticObjectType()) {
    return QStyledItemDelegate::setEditorData(editor, index);
  } else if (p_type ==
             serialization::serializable_ptr_scheme::getStaticObjectType()) {
    QComboBox* box = static_cast<QComboBox*>(editor);
    box->setCurrentIndex(
        box->findText(parentModel->data(index, Qt::DisplayRole).toString()));
    return;
  };

  QStyledItemDelegate::setEditorData(editor, index);
  return;
};

void ObjPropertiesQtDelegate::setModelData(QWidget* editor,
                                           QAbstractItemModel* model,
                                           const QModelIndex& index) const {
  if ((!index.isValid()) || (index.column() != 1))
    return QStyledItemDelegate::setModelData(editor, model, index);

  std::pair<std::string, std::shared_ptr<serialization::type_scheme>> fld =
      parentModel->current_fields.get_field(index.row());
  rtti::so_type* p_type = fld.second->getObjectType();

  if (p_type == serialization::primitive_scheme<bool>::getStaticObjectType()) {
    QComboBox* box = static_cast<QComboBox*>(editor);
    parentModel->setData(index, box->itemData(box->currentIndex()),
                         Qt::EditRole);
    return;
  } else if (p_type ==
             serialization::primitive_scheme<int>::getStaticObjectType()) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->interpretText();
    parentModel->setData(index, QVariant(box->value()), Qt::EditRole);
    return;
  } else if ((p_type == serialization::primitive_scheme<
                            unsigned char>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            unsigned int>::getStaticObjectType()) ||
             (p_type == serialization::primitive_scheme<
                            long unsigned int>::getStaticObjectType())) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->interpretText();
    parentModel->setData(index, QVariant(box->value()), Qt::EditRole);
    return;
  } else if ((p_type ==
              serialization::primitive_scheme<float>::getStaticObjectType()) ||
             (p_type ==
              serialization::primitive_scheme<double>::getStaticObjectType())) {
    QDoubleSpinBox* box = static_cast<QDoubleSpinBox*>(editor);
    box->interpretText();
    parentModel->setData(index, QVariant(box->value()), Qt::EditRole);
    return;
  } else if (p_type ==
             serialization::primitive_scheme<char>::getStaticObjectType()) {
    QSpinBox* box = static_cast<QSpinBox*>(editor);
    box->interpretText();
    parentModel->setData(index, QVariant(box->value()), Qt::EditRole);
    return;
  } else if (p_type == serialization::primitive_scheme<
                           std::string>::getStaticObjectType()) {
    return QStyledItemDelegate::setModelData(editor, model, index);
  } else if (p_type ==
             serialization::serializable_ptr_scheme::getStaticObjectType()) {
    QComboBox* box = static_cast<QComboBox*>(editor);
    QString value = box->itemText(box->currentIndex());
    if (value != "New object...")
      parentModel->setData(index, value, Qt::EditRole);
    else {
      // This requires the creation of a new object.
      // First, get the list of derived types.
      Ui::RKClassSelectWidget diag_w;
      QDialog diag;
      diag_w.setupUi(&diag);

      rtti::so_type* basetype_ptr =
          rtti::getRKSharedObjTypeRepo().findType(fld.second->get_type_ID());
      if (!basetype_ptr)
        return;

      std::stack<rtti::so_type*> btype_stack;
      btype_stack.push(basetype_ptr);
      while (!btype_stack.empty()) {
        rtti::so_type* tmp_bt = btype_stack.top();
        btype_stack.pop();
        if (tmp_bt->isConcrete()) {
          QListWidgetItem* itm_ptr =
              new QListWidgetItem(QString::fromStdString(tmp_bt->TypeName()));
          itm_ptr->setData(
              Qt::UserRole,
              QVariant(reinterpret_cast<quint64>(tmp_bt->TypeID_begin())));
          diag_w.listWidget->addItem(itm_ptr);
        };
        for (std::size_t i = 0; i < tmp_bt->getDirectDescendantCount(); ++i) {
          rtti::so_type* tmp = tmp_bt->getDirectDescendant(i);
          if (tmp)
            btype_stack.push(tmp);
        };
      };

      // Execute the "derived class select" dialog.
      int diag_result = diag.exec();
      if ((diag_result == QDialog::Rejected) ||
          (diag_w.listWidget->currentIndex().row() < 0))
        return;
      const unsigned int* newobj_typeID = reinterpret_cast<const unsigned int*>(
          diag_w.listWidget->currentItem()->data(Qt::UserRole).toULongLong());
      rtti::so_type* newobj_typeptr =
          basetype_ptr->findDescendant(newobj_typeID);
      if (!newobj_typeptr)
        return;

      // Execute the "set initial values" dialog.
      std::shared_ptr<shared_object> newobj_ptr =
          newobj_typeptr->CreateObject();
      parentModel->setDataNewObject(index, newobj_ptr);
    };
    return;
  };

  QStyledItemDelegate::setEditorData(editor, index);
  return;
};

void ObjPropertiesQtDelegate::updateEditorGeometry(
    QWidget* editor, const QStyleOptionViewItem& option,
    const QModelIndex& index) const {
  return QStyledItemDelegate::updateEditorGeometry(editor, option, index);
};
};  // namespace qt

};  // namespace ReaK
