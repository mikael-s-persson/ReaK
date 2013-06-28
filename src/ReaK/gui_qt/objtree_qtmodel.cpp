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

#include "objtree_qtmodel.hpp"

#include <base/named_object.hpp>
#include "graph_alg/bgl_tree_adaptor.hpp"

namespace ReaK {

namespace rkqt {
  
  
void ObjTreeQtModel::updateObjTreeNode(obj_indirect_node_desc u) {
  typedef boost::graph_traits< serialization::object_graph >::out_edge_iterator OutEdgeIter;
  OutEdgeIter eo, eo_end;
  boost::tie(eo, eo_end) = out_edges(obj_tree[u].g_node, *obj_graph);
  for(; eo != eo_end; ++eo)
    updateObjTreeNode(add_child_vertex(u, objtree_indirect_node(target(*eo,*obj_graph)), obj_tree).first);
};

void ObjTreeQtModel::refreshObjTree() {
  emit aboutToResetModel();
  beginResetModel();
  obj_tree.clear();
  obj_indirect_node_desc t_root = create_root(objtree_indirect_node(root_node), obj_tree);
  updateObjTreeNode(t_root);
  endResetModel();
  emit justAfterModelReset();
};

 
//     shared_ptr< serialization::object_graph > obj_graph;
//     serialization::object_node_desc root_node;  

ObjTreeQtModel::ObjTreeQtModel(const shared_ptr< serialization::object_graph >& aObjGraph, 
                               serialization::object_node_desc aRoot) :
                               obj_graph(aObjGraph), root_node(aRoot) {
  obj_indirect_node_desc t_root = create_root(objtree_indirect_node(root_node), obj_tree);
  updateObjTreeNode(t_root);
};

ObjTreeQtModel::~ObjTreeQtModel() { };

QVariant ObjTreeQtModel::data(const QModelIndex &index, int role) const {
  if(!index.isValid())
    return QVariant();
  
  if(role != Qt::DisplayRole)
    return QVariant();
  
  serialization::object_node_desc item(obj_tree[obj_indirect_node_desc(index.internalId())].g_node);
  
  if(index.column() == 0) {
    return QVariant(QString::fromStdString(serialization::get_objtree_name(*obj_graph, item)));
  } else {  // column 1
    shared_ptr< serialization::serializable > item_ptr = (*obj_graph)[item].p_obj;
    std::string s_tmp = item_ptr->getObjectType()->TypeName();
    if(s_tmp.length() > 40) {
      // Try to eliminate the template arguments.
      std::string::iterator it = std::find(s_tmp.begin(), s_tmp.end(), '<');
      if(it == s_tmp.end())
        return QVariant(QString(s_tmp.c_str()));
      std::string::reverse_iterator rit = std::find(s_tmp.rbegin(), s_tmp.rend(), '>');
      return QVariant(QString((std::string(s_tmp.begin(), it) + "<..>" + std::string(rit.base(), s_tmp.end())).c_str()));
    } else {
      return QVariant(QString(s_tmp.c_str()));
    };
  };
  return QVariant();
};

Qt::ItemFlags ObjTreeQtModel::flags(const QModelIndex &index) const {
  if (!index.isValid())
    return 0;
  
  if(index.column() < 1)
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
  else
    return Qt::ItemIsEnabled;
};

QVariant ObjTreeQtModel::headerData(int section, Qt::Orientation orientation, int role) const {
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
    if(section == 0)
      return QVariant(QString("Object Name"));
    else if(section == 1)
      return QVariant(QString("Type"));
  };
  
  return QVariant();  
};

QModelIndex ObjTreeQtModel::index(int row, int column, const QModelIndex &parent) const {
  if (!hasIndex(row, column, parent))
    return QModelIndex();
  
  obj_indirect_node_desc parentItem = 0;
  
  if(parent.isValid())
    parentItem = obj_indirect_node_desc(parent.internalId());
  
  if(row >= int(out_degree(parentItem, obj_tree)))
    return QModelIndex();
  
  typedef boost::graph_traits< obj_indirect_tree_type >::out_edge_iterator OutEdgeIter;
  OutEdgeIter eo, eo_end;
  boost::tie(eo, eo_end) = out_edges(parentItem, obj_tree);
  std::advance(eo, row);
  obj_indirect_node_desc child = target(*eo, obj_tree);
  return createIndex(row, column, quint32(child));
};

QModelIndex ObjTreeQtModel::parent(const QModelIndex &index) const {
  if(!index.isValid())
    return QModelIndex();
  
  obj_indirect_node_desc childItem(index.internalId());
  typedef boost::graph_traits< obj_indirect_tree_type >::in_edge_iterator InEdgeIter;
  InEdgeIter ei, ei_end;
  boost::tie(ei,ei_end) = in_edges(childItem, obj_tree);
  if(ei == ei_end)
    return QModelIndex();
  obj_indirect_node_desc parentItem = source(*ei, obj_tree);
  
  if(parentItem == root_node)
    return QModelIndex();
  
  // find the row in the parent of the parent node:
  typedef boost::graph_traits< obj_indirect_tree_type >::out_edge_iterator OutEdgeIter;
  OutEdgeIter eo, eo_end;
  obj_indirect_node_desc grandparent = source(*(in_edges(parentItem, obj_tree).first), obj_tree);
  boost::tie(eo, eo_end) = out_edges(grandparent, obj_tree);
  int i = 0;
  for(; eo != eo_end; ++eo, ++i) 
    if(target(*eo, obj_tree) == parentItem)
      break;
  
  return createIndex(i, 0, quint32(parentItem));
};

int ObjTreeQtModel::rowCount(const QModelIndex &parent) const {
  if(!parent.isValid()) 
    return out_degree(obj_indirect_node_desc(0), obj_tree);
  else
    return out_degree(obj_indirect_node_desc(parent.internalId()), obj_tree);
};

int ObjTreeQtModel::columnCount(const QModelIndex &parent) const {
  return 2;
};

void ObjTreeQtModel::treeChanged() {
  refreshObjTree();
};

void ObjTreeQtModel::itemSelected(const QModelIndex &index) {
  if(index.isValid())
    emit objectNodeSelected(obj_tree[obj_indirect_node_desc(index.internalId())].g_node);
};

};

}; //ReaK




