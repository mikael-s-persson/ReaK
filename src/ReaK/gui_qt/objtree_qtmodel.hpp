/**
 * \file objtree_qtmodel.hpp
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

#ifndef REAK_OBJTREE_QTMODEL_HPP
#define REAK_OBJTREE_QTMODEL_HPP

#include <ReaK/core/serialization/objtree_archiver.hpp>

#include <QAbstractItemModel>
#include <QVariant>
#include <QString>
#include <QModelIndex>

namespace ReaK {

namespace rkqt {

class ObjTreeQtModel : public QAbstractItemModel {
  Q_OBJECT
  
  public:
    ObjTreeQtModel(const shared_ptr< serialization::object_graph >& aObjGraph, serialization::object_node_desc aRoot = 0);
    ~ObjTreeQtModel();
    
    // All the virtual functions required by QAbstractItemModel:
    QVariant data(const QModelIndex &index, int role) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;
    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    
    /**
     * This function refreshes the object-tree representation to synchronize it to the object-graph.
     */
    void refreshObjTree();
    
  private:
    shared_ptr< serialization::object_graph > obj_graph;
    serialization::object_node_desc root_node;  
    
    struct objtree_indirect_node {
      serialization::object_node_desc g_node;
      objtree_indirect_node(serialization::object_node_desc aGNode = serialization::object_node_desc()) : g_node(aGNode) { };
    };
    typedef boost::adjacency_list_BC< boost::vecBC, boost::vecBC, boost::bidirectionalS, objtree_indirect_node > obj_indirect_tree_type;
    typedef boost::graph_traits< obj_indirect_tree_type >::vertex_descriptor obj_indirect_node_desc;
    obj_indirect_tree_type obj_tree;
    
    void updateObjTreeNode(obj_indirect_node_desc u);
    
  public slots:
    
    void itemSelected(const QModelIndex &index);
    void treeChanged();
    
  signals:
    
    void aboutToResetModel();
    void justAfterModelReset();
    void objectNodeSelected(serialization::object_node_desc u);
    
  public:
    /**
     * This function returns the shared-pointer to the object-graph wrapped by this model.
     * \return The shared-pointer to the object-graph wrapped by this model.
     */
    shared_ptr< serialization::object_graph > get_object_graph() const { return obj_graph; };
    /**
     * This function returns the root node of the object-graph wrapped by this model.
     * \return The root node of the object-graph wrapped by this model.
     */
    serialization::object_node_desc get_root_node() const { return root_node; };
    
};


};

}; //ReaK

#endif






