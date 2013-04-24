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
  
  
namespace detail {
  
  
  enum bst_traversal_status {
    OnLeftBranch,
    OnMiddleBranch,
    OnRightBranch
  };
  
  
  template <typename TreeType, typename VertexType>
  VertexType bst_go_down_left(const TreeType& aTree, VertexType aStart) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(aStart, aTree);
    // look down-left until we reach the leaf:
    while(cur_children.first != cur_children.second) {
      aStart = *(cur_children.first);
      cur_children = child_vertices(aStart, aTree);
    };
    return aStart;
  };
  
  template <typename TreeType, typename VertexType>
  VertexType bst_go_down_right(const TreeType& aTree, VertexType aStart) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(aStart, *aTree);
    // look down-right until we reach the leaf:
    while((cur_children.first != cur_children.second) && (++(cur_children.first) != cur_children.second)) {
      aStart = *(cur_children.first);
      cur_children = child_vertices(aStart, *aTree);
    };
    return aStart;
  };
  
  template <typename TreeType, typename VertexType>
  void bst_move_up_to_next(const TreeType& aTree, VertexType& aU, bst_traversal_status& aStatus) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    typedef typename boost::graph_traits<TreeType>::in_edge_iterator in_edge_iter;
    aStatus = OnRightBranch;
    while(true) {
      in_edge_iter ei, ei_end;
      boost::tie(ei, ei_end) = in_edges(aU, aTree);
      if(ei == ei_end) // at the root, go to the end (root, rightbranch)
        return;
      vertex_type v = source(*ei, aTree);
      child_vertex_iter vil, vi_end;
      boost::tie(vil, vi_end) = child_vertices(v, aTree);
      if(*vil == aU) { // u is the left child of v.
        aU = v;
        aStatus = OnMiddleBranch;
        return;
      }
      // u must be the right child of v. keep going.
      aU = v;
    };
  };
  
  template <typename TreeType, typename VertexType>
  void bst_move_up_to_prev(const TreeType& aTree, VertexType& aU, bst_traversal_status& aStatus) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    typedef typename boost::graph_traits<TreeType>::in_edge_iterator in_edge_iter;
    child_vertex_iter vil, vi_end;
    aStatus = OnLeftBranch;
    vertex_type u = aU;
    while(true) {
      in_edge_iter ei, ei_end;
      boost::tie(ei, ei_end) = in_edges(u, aTree);
      if(ei == ei_end) // at the root, so, aU must be the beginning node.
        return;
      vertex_type v = source(*ei, aTree);
      child_vertex_iter vil, vi_end;
      boost::tie(vil, vi_end) = child_vertices(v, aTree);
      if(*vil == u) { // u is the left child of v.
        u = v;
        continue;
      };
      // u must be the right child of v. keep going.
      aU = v;
      aStatus = OnMiddleBranch;
      return;
    };
  };
  
};



template <typename CompleteBinaryTree, typename ValueType>
class bst_inorder_iterator {
  public:
    typedef bst_inorder_iterator<CompleteBinaryTree> self;
    typedef CompleteBinaryTree tree_type;
    
    typedef std::ptrdiff_t difference_type;
    typedef ValueType value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::bidirectional_iterator_tag iterator_category;
    
  private:
    typedef typename boost::graph_traits<tree_type>::vertex_descriptor vertex_type;
    typedef typename tree_traits<tree_type>::child_vertex_iterator child_vertex_iter;
    typedef typename boost::graph_traits<tree_type>::in_edge_iterator in_edge_iter;
    
    
    tree_type* m_tree;
    vertex_type m_u;
    detail::bst_traversal_status m_status;
    
    bst_inorder_iterator(tree_type* aTree, vertex_type aU, detail::bst_traversal_status aStatus) : m_tree(aTree), m_u(aU), m_status(aStatus) { };
    
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
    
    void move_up_to_prev() {
      in_edge_iter ei, ei_end;
      child_vertex_iter vil, vi_end;
      m_status = OnLeftBranch;
      vertex_type u = m_u;
      while(true) {
        boost::tie(ei, ei_end) = in_edges(u, *m_tree);
        if(ei == ei_end) // at the root, so, m_u must be the beginning node.
          return;
        vertex_type v = source(*ei, *m_tree);
        boost::tie(vil, vi_end) = child_vertices(v, *m_tree);
        if(*vil == u) { // u is the left child of v.
          u = v;
          continue;
        };
        // u must be the right child of v. keep going.
        m_u = v;
        m_status = OnMiddleBranch;
        return;
      };
    };
    
  public:
    
    vertex_type base() const { return m_u; };
    
    bst_inorder_iterator(tree_type* aTree, vertex_type aU) : m_tree(aTree), m_u(aU), m_status(detail::OnRightBranch) {
      // must figure out what the case is.
      if((!m_tree) || (m_u == boost::graph_traits<tree_type>::null_vertex()))
        return;
      if(m_u == get_root_vertex(*m_tree)) {
        m_status = detail::OnMiddleBranch;
        return;
      };
      // first check if there are any children:
      std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
      if(cur_children.first != cur_children.second) {
        m_status = detail::OnMiddleBranch; // not on leaf.
        return;
      };
      // then, check if m_u is a left or right child of its parent:
      in_edge_iter ei, ei_end;
      boost::tie(ei, ei_end) = in_edges(m_u, *m_tree);
      vertex_type v = source(*ei, *m_tree);
      child_vertex_iter vil, vi_end;
      boost::tie(vil, vi_end) = child_vertices(v, *m_tree);
      if(*vil == m_u) // u is the left child of v.
        m_status = detail::OnLeftBranch;
      else // u must be the right child of v.
        m_status = detail::OnRightBranch;
    };
    
    static self begin(tree_type* aTree) {
      if(aTree)
        return self(aTree, go_down_left(aTree, get_root_vertex(*aTree)), detail::OnLeftBranch);
      else
        return self(NULL, boost::graph_traits<tree_type>::null_vertex(), detail::OnRightBranch);
    };
    
    static self end(tree_type* aTree) {
      if(aTree)
        return self(aTree, get_root_vertex(*aTree), detail::OnRightBranch);
      else
        return self(NULL, boost::graph_traits<tree_type>::null_vertex(), detail::OnRightBranch);
    };
    
    friend bool operator==( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree == rhs.m_tree) && (lhs.m_u == rhs.m_u) && (lhs.m_status == rhs.m_status));
    };
    friend bool operator!=( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree != rhs.m_tree) || (lhs.m_u != rhs.m_u) || (lhs.m_status != rhs.m_status));
    };
    
    self& operator++() { 
      if(!m_tree)
        return *this;
      switch(m_status) {
        case detail::OnLeftBranch:
          // on a left-leaf, must move up to the parent as a middle-point
          std::pair<in_edge_iter, in_edge_iter> eis = in_edges(m_u, *m_tree);
          if(eis.first == eis.second) { // at the root, go to the end (root, rightbranch).
            m_status = detail::OnRightBranch;
            return *this;
          };
          m_u = source(*(eis.first), *m_tree);
          m_status = detail::OnMiddleBranch;
          break;
        case detail::OnMiddleBranch:
          // on a middle-point, must move down to the right once, and then left to the bottom.
          std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
          ++(cur_children.first);
          if(cur_children.first != cur_children.second) {
            // go to the right child.
            m_u = go_down_left(m_tree, *(cur_children.first));
            if(m_u == *(cur_children.first))
              m_status = detail::OnRightBranch;
            else
              m_status = detail::OnLeftBranch;
            break;
          };
          // this means that we must move up to the next value (no right child here).
        case detail::OnRightBranch:
          detail::bst_move_up_to_next(*m_tree, m_u, m_status);
          break;
      };
      return *this;
    };
    
    self& operator--() { 
      if(!m_tree)
        return *this;
      switch(m_status) {
        case detail::OnRightBranch:
          if(m_u == get_root_vertex(*m_tree)) {
            // go to the left child, and down the right:
            m_u = go_down_right(m_tree, m_u);
            if(m_u == get_root_vertex(*m_tree))
              m_status = detail::OnMiddleBranch;
            else {
              std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
              if(cur_children.first == cur_children.second) 
                m_status = detail::OnRightBranch;  // on the leaf.
              else
                m_status = detail::OnMiddleBranch; // not on leaf.
            };
            break;
          };
          // this means that we are either on a right-leaf or on a mis-labeled middle-node, 
          // in either case, try the middle-branch case:
        case detail::OnMiddleBranch:
          // on a middle-point or right-point, must move down to the left once (if exists), and then right to the bottom.
          std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
          if(cur_children.first != cur_children.second) {
            // go to the left child, and down the right:
            m_u = go_down_right(m_tree, *(cur_children.first));
            if(m_u == *(cur_children.first))
              m_status = detail::OnMiddleBranch;
            else {
              cur_children = child_vertices(m_u, *m_tree);
              if(cur_children.first == cur_children.second) 
                m_status = detail::OnRightBranch;  // on the leaf.
              else
                m_status = detail::OnMiddleBranch; // not on leaf.
            };
            break;
          };
          // this means that we must move up to the previous value (no left child here).
        case detail::OnLeftBranch:
          detail::bst_move_up_to_prev(*m_tree, m_u, m_status);
          break;
      };
      return *this;
    };
    
    self operator++(int) { self result(*this); return ++result; };
    self operator--(int) { self result(*this); return --result; };
    
    reference operator*() { return (*m_tree)[m_u]; };
    pointer operator->() { return &(*m_tree)[m_u]; };
};








};

};


#endif









