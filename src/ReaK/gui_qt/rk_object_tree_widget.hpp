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

#ifndef REAK_OBJECT_TREE_WIDGET_HPP
#define REAK_OBJECT_TREE_WIDGET_HPP

#include "objtree_qtmodel.hpp"

#include "ui_rk_object_tree.h"
#include <QDockWidget>

namespace ReaK {
  
namespace rkqt {

class ObjectTreeWidget : public QDockWidget, private Ui::RKObjectTreeWidget {
    Q_OBJECT
  
  public:
    ObjectTreeWidget(const shared_ptr< serialization::object_graph >& aObjGraph, 
                     serialization::object_node_desc aRoot = 0,
                     QWidget * parent = NULL, Qt::WindowFlags flags = 0);
    virtual ~ObjectTreeWidget() { };
    
  private slots:
    
    void recordPreviousSelection();
    void restorePreviousSelection();
    
  public:
    
    ObjTreeQtModel mdl;
    QString prev_selection;
    
};

};

};

#endif














