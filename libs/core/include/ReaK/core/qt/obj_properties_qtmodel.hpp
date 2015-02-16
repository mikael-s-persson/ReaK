/**
 * \file obj_properties_qtmodel.hpp
 *
 * This library declares the class 
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date December 2012
 */

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

#ifndef REAK_OBJ_PROPERTIES_QTMODEL_HPP
#define REAK_OBJ_PROPERTIES_QTMODEL_HPP

#include <ReaK/core/serialization/objtree_archiver.hpp>

#include <QAbstractTableModel>
#include <QStyledItemDelegate>
#include <QVariant>
#include <QString>
#include <QModelIndex>

namespace ReaK {

namespace qt {
  
class ObjPropertiesQtDelegate; // forward-declare.


class ObjPropertiesQtModel : public QAbstractTableModel {
  Q_OBJECT
  
  public:
    ObjPropertiesQtModel(const shared_ptr< serialization::object_graph >& aObjGraph, serialization::object_node_desc aCurrentNode);
    ~ObjPropertiesQtModel();
    
    // All the virtual functions required by QAbstractItemModel:
    int rowCount(const QModelIndex& parent = QModelIndex()) const;
    int columnCount(const QModelIndex& parent = QModelIndex()) const;
    QVariant data(const QModelIndex& index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    Qt::ItemFlags flags(const QModelIndex& index) const;
    bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole);
    
    void setDataNewObject(const QModelIndex& index, const shared_ptr< shared_object >& aNewPtr);
    
  private:
    serialization::objtree_editor obj_graph_editor;
    serialization::object_node_desc current_node;
    serialization::xml_field_editor current_fields;
    
  public slots:
    
    void selectObject(serialization::object_node_desc u);
    void sourceDataEdited(const std::string& newSrc);
    
  signals:
    
    void objectNameChanged(const std::string& newName);
    void sourceDataChanged();
    void objectTreeChanged();
    
  public:
    friend class ObjPropertiesQtDelegate;
    
    /**
     * This function returns the shared-pointer to the object-graph wrapped by this model.
     * \return The shared-pointer to the object-graph wrapped by this model.
     */ 
    const serialization::objtree_editor& get_object_editor() const { return obj_graph_editor; };
    serialization::objtree_editor& get_object_editor() { return obj_graph_editor; };
    /**
     * This function returns the root node of the object-graph wrapped by this model.
     * \return The root node of the object-graph wrapped by this model.
     */
    serialization::object_node_desc get_current_node() const { return current_node; };
    
    /**
     * This function return the current field-editor object (by const-reference) which 
     * describes the property fields of the current object.
     * \return The current field-editor object (by const-reference).
     */
    const serialization::xml_field_editor& get_current_fields() const { return current_fields; };
    
};

class ObjPropertiesQtDelegate : public QStyledItemDelegate {
    Q_OBJECT
  
  public:
    ObjPropertiesQtDelegate(ObjPropertiesQtModel* aParent);
    virtual ~ObjPropertiesQtDelegate();
    
    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
    void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;
    
  private:
    ObjPropertiesQtModel* parentModel;
    
  public:
    
};


};

}; //ReaK

#endif














