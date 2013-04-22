/**
 * \file bst_inorder_iterator.hpp
 * 
 * This library implements a simple binary search tree in-order iterator for a BGL tree.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_BST_INORDER_ITERATOR_HPP
#define REAK_BST_INORDER_ITERATOR_HPP

#include "base/defs.hpp"

#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include "graph_alg/tree_concepts.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_raw_property_graph.hpp"

namespace ReaK {

namespace graph {



template <typename CompleteBinaryTree>
class vertex_iterator {
  public:
    typedef vertex_iterator<CompleteBinaryTree> self;
    typedef CompleteBinaryTree tree_type;
    
    typedef std::ptrdiff_t difference_type;
    typedef typename tree_type::vertex_bundled value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::bidirectional_iterator_tag iterator_category;
    
  private:
    typedef typename boost::graph_traits<tree_type>::vertex_descriptor vertex_type;
    typedef typename tree_traits<tree_type>::child_vertex_iterator child_vertex_iter;
    typedef typename boost::graph_traits<tree_type>::in_edge_iterator in_edge_iter;
    
    enum traversal_status {
      OnLeftBranch,
      OnMiddleBranch,
      OnRightBranch
    };
    
    tree_type* m_tree;
    vertex_type m_u;
    traversal_status m_status;
    
    vertex_iterator(tree_type* aTree, vertex_type aU, traversal_status aStatus) : m_tree(aTree), m_u(aU), m_status(aStatus) { };
    
    static vertex_type go_down_left(tree_type* aTree, vertex_type aStart) {
      std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(aStart, *aTree);
      // look down-left until we reach the leaf:
      while(cur_children.first != cur_children.second) {
        aStart = *(cur_children.first);
        cur_children = child_vertices(aStart, *aTree);
      };
      return aStart;
    };
    
    static vertex_type go_down_right(tree_type* aTree, vertex_type aStart) {
      std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(aStart, *aTree);
      // look down-right until we reach the leaf:
      while((cur_children.first != cur_children.second) && (++(cur_children.first) != cur_children.second)) {
        aStart = *(cur_children.first);
        cur_children = child_vertices(aStart, *aTree);
      };
      return aStart;
    };
    
    void move_up_to_next() {
      in_edge_iter ei, ei_end;
      child_vertex_iter vil, vi_end;
      m_status = OnRightBranch;
      while(true) {
        boost::tie(ei, ei_end) = in_edges(m_u, *m_tree);
        if(ei == ei_end) // at the root, go to the end (root, rightbranch)
          return;
        vertex_type v = source(*ei, *m_tree);
        boost::tie(vil, vi_end) = child_vertices(v, *m_tree);
        if(*vil == m_u) { // u is the left child of v.
          m_u = v;
          m_status = OnMiddleBranch;
          return;
        }
        // u must be the right child of v. keep going.
        m_u = v;
      };
    };
    
  public:
    
    static self begin(tree_type* aTree) {
      if(aTree)
        return self(aTree, go_down_left(aTree, get_root_vertex(*aTree)), OnLeftBranch);
      else
        return self(NULL, boost::graph_traits<tree_type>::null_vertex(), OnRightBranch);
    };
    
    static self end(tree_type* aTree) {
      if(aTree)
        return self(aTree, get_root_vertex(*aTree), OnRightBranch);
      else
        return self(NULL, boost::graph_traits<tree_type>::null_vertex(), OnRightBranch);
    };
    
    friend bool operator==( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree == rhs.m_tree) && (lhs.m_u == rhs.m_u));
    };
    friend bool operator!=( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree != rhs.m_tree) || (lhs.m_u != rhs.m_u));
    };
    
    self& operator++() { 
      if(!m_tree)
        return *this;
      switch(m_status) {
        case OnLeftBranch:
          // on a left-leaf, must move up to the parent as a middle-point
          std::pair<in_edge_iter, in_edge_iter> eis = in_edges(m_u, *m_tree);
          if(eis.first == eis.second) { // at the root, go to the end (root, rightbranch).
            m_status = OnRightBranch;
            return *this;
          };
          m_u = source(*(eis.first), *m_tree);
          m_status = OnMiddleBranch;
          break;
        case OnMiddleBranch:
          // on a middle-point, must move down to the right once, and then left to the bottom.
          std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
          ++(cur_children.first);
          
          break;
        case OnRightBranch:
          
          break;
      };
      return *this;
    };
    self operator++(int) { self result(*this); return ++result; };
    self& operator--() { --base; return *this; };
    self operator--(int) { self result(*this); return --result; };
    
    value_type operator[](difference_type i) const { return base + i; };
    reference operator*() { return base; };
    pointer operator->() { return &base; };
};








};

};


#endif









