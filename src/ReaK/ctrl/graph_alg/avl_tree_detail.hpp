/**
 * \file avl_tree_detail.hpp
 * 
 * This library implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN) 
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which 
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * 
 * \todo Currently, there is no support for set/map (i.e., there is no enforcement of unique-keys).
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

#ifndef REAK_AVL_TREE_DETAIL_HPP
#define REAK_AVL_TREE_DETAIL_HPP

#include <ReaK/core/base/defs.hpp>

#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

// BGL-Extra includes:
#include <boost/graph/tree_traits.hpp>
#include <boost/graph/tree_adaptor.hpp>
#include <boost/graph/bst_inorder_iterator.hpp>


#include <vector>
#include <queue>
#include <iterator>


namespace ReaK {

namespace graph {
  
namespace detail {
  
  enum avl_container_style {
    avl_set_style,
    avl_multiset_style,
    avl_map_style,
    avl_multimap_style
  };
  
  
  template <typename TreeType,
            typename Compare,
            avl_container_style ContainerStyle>
  struct avl_tree_helper {
    typedef typename TreeType::vertex_bundled value_type;
    typedef value_type key_type;
    typedef value_type mapped_type;
    
    typedef Compare value_compare;
    typedef Compare key_compare;
    
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    
    static const bool allow_duplicates = false;
    
    static const key_type& value_to_key(const value_type& rhs) { return rhs; };
    
    static const mapped_type& value_to_mapped(const value_type& rhs) { return rhs; };
    static mapped_type& value_to_mapped(value_type& rhs) { return rhs; };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    template <typename U, typename V>
    static value_type keymap_to_value(U&& k, V&&) { 
      return value_type(std::forward<U>(k));
    };
#endif
    static value_type keymap_to_value(const key_type& k, const mapped_type&) { 
      return value_type(k);
    };
    
  };
  
  
  template <typename TreeType,
            typename Compare>
  struct avl_tree_helper<TreeType, Compare, avl_multiset_style> {
    typedef typename TreeType::vertex_bundled value_type;
    typedef value_type key_type;
    typedef value_type mapped_type;
    
    typedef Compare value_compare;
    typedef Compare key_compare;
    
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    
    static const bool allow_duplicates = true;
    
    static const key_type& value_to_key(const value_type& rhs) { return rhs; };
    
    static const mapped_type& value_to_mapped(const value_type& rhs) { return rhs; };
    static mapped_type& value_to_mapped(value_type& rhs) { return rhs; };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    template <typename U, typename V>
    static value_type keymap_to_value(U&& k, V&&) { 
      return value_type(std::forward<U>(k));
    };
#endif
    static value_type keymap_to_value(const key_type& k, const mapped_type&) { 
      return value_type(k);
    };
    
  };
  
  
  template <typename TreeType,
            typename Compare>
  struct avl_tree_helper<TreeType, Compare, avl_map_style> {
    typedef typename TreeType::vertex_bundled value_type;
    typedef typename value_type::first_type key_type;
    typedef typename value_type::second_type mapped_type;
    
    struct value_compare {
      Compare m_compare;
      
      value_compare(Compare aComp) : m_compare(aComp) { };
      
      bool operator()(const value_type& lhs, const value_type& rhs) const {
        return m_compare(lhs.first, rhs.first);
      };
      bool operator()(const value_type& lhs, const key_type& rhs) const {
        return m_compare(lhs.first, rhs);
      };
      bool operator()(const key_type& lhs, const value_type& rhs) const {
        return m_compare(lhs, rhs.first);
      };
      bool operator()(const key_type& lhs, const key_type& rhs) const {
        return m_compare(lhs, rhs);
      };
      
    };
    typedef value_compare key_compare;
    
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    
    static const bool allow_duplicates = false;
    
    static const key_type& value_to_key(const value_type& rhs) { return rhs.first; };
    
    static const mapped_type& value_to_mapped(const value_type& rhs) { return rhs.second; };
    static mapped_type& value_to_mapped(value_type& rhs) { return rhs.second; };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    template <typename U, typename V>
    static value_type keymap_to_value(U&& k, V&& m) { 
      return value_type(std::forward<U>(k), std::forward<V>(m));
    };
#endif
    static value_type keymap_to_value(const key_type& k, const mapped_type& m) { 
      return value_type(k, m);
    };
    
  };
  
  
  template <typename TreeType,
            typename Compare>
  struct avl_tree_helper<TreeType, Compare, avl_multimap_style> {
    typedef typename TreeType::vertex_bundled value_type;
    typedef typename value_type::first_type key_type;
    typedef typename value_type::second_type mapped_type;
    
    struct value_compare {
      Compare m_compare;
      
      value_compare(Compare aComp) : m_compare(aComp) { };
      
      bool operator()(const value_type& lhs, const value_type& rhs) const {
        return m_compare(lhs.first, rhs.first);
      };
      bool operator()(const value_type& lhs, const key_type& rhs) const {
        return m_compare(lhs.first, rhs);
      };
      bool operator()(const key_type& lhs, const value_type& rhs) const {
        return m_compare(lhs, rhs.first);
      };
      bool operator()(const key_type& lhs, const key_type& rhs) const {
        return m_compare(lhs, rhs);
      };
      
    };
    typedef value_compare key_compare;
    
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    
    static const bool allow_duplicates = true;
    
    static const key_type& value_to_key(const value_type& rhs) { return rhs.first; };
    
    static const mapped_type& value_to_mapped(const value_type& rhs) { return rhs.second; };
    static mapped_type& value_to_mapped(value_type& rhs) { return rhs.second; };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    template <typename U, typename V>
    static value_type keymap_to_value(U&& k, V&& m) { 
      return value_type(std::forward<U>(k), std::forward<V>(m));
    };
#endif
    static value_type keymap_to_value(const key_type& k, const mapped_type& m) { 
      return value_type(k, m);
    };
    
  };
  
  
};


/**
 * This class implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN) 
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which 
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * \tparam TreeType The tree type to be used to store the entries of this AVL tree.
 * \tparam Compare The comparison functor type that can compare two elements. 
 */
template <typename TreeType,
          typename Compare,
          detail::avl_container_style ContainerStyle>
class avl_tree_impl
{
  public:
    
    typedef avl_tree_impl<TreeType, Compare, ContainerStyle> self;
    
    typedef detail::avl_tree_helper<TreeType, Compare, ContainerStyle> helper_type;
    
    typedef typename helper_type::value_type value_type;
    typedef typename helper_type::key_type key_type;
    typedef typename helper_type::mapped_type mapped_type;
    typedef typename helper_type::value_compare value_compare;
    typedef typename helper_type::key_compare key_compare;
    
    typedef typename helper_type::reference reference;
    typedef typename helper_type::const_reference const_reference;
    typedef typename helper_type::pointer pointer;
    typedef typename helper_type::const_pointer const_pointer;
    
    typedef boost::bst_inorder_iterator<const TreeType, const value_type> iterator;
    typedef iterator const_iterator;
    
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::allocator<value_type> allocator_type;
    
    
  private:
    
    typedef TreeType tree_indexer;
    
    typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<tree_indexer>::edge_descriptor edge_type;
    typedef typename boost::graph_traits<tree_indexer>::out_edge_iterator out_edge_iter;
    typedef typename boost::graph_traits<tree_indexer>::in_edge_iterator in_edge_iter;
    typedef typename boost::tree_traits<tree_indexer>::child_vertex_iterator child_vertex_iter;
    
    typedef typename tree_indexer::vertex_property_type vertex_property;
    typedef typename tree_indexer::edge_property_type edge_property;
    
    //typedef detail::compare_pair_first< distance_type, vertex_type, std::less< distance_type > > priority_compare_type;
    //typedef std::vector< std::pair< distance_type, vertex_type > > priority_queue_type;
    
    tree_indexer m_tree;   ///< Tree storage.
    vertex_type m_root;    ///< Root node of the tree.
    value_compare m_compare;  ///< The comparison functor.
    
    struct construction_task {
      typedef typename std::vector<vertex_property>::iterator PropIter;
      
      vertex_type node;
      typename std::vector<vertex_property>::iterator first;
      typename std::vector<vertex_property>::iterator last;
      
      construction_task(vertex_type aNode, PropIter aBegin, PropIter aEnd) : 
        node(aNode), first(aBegin), last(aEnd) { };
    };
    
    
#if 0
    /* for debugging purposes */
    void print_complete_tree() const {
      using std::swap;
      
      if(m_root == tree_indexer::null_vertex()) {
        std::cout << " <empty>" << std::endl;
        return;
      };
      
      std::queue<vertex_type> parents;
      parents.push(m_root);
      std::cout << std::setw(8) << helper_type::value_to_key(m_tree[m_root]) << std::endl;
      bool finished = false;
      while(!finished) {
        // print parents and collect children:
        std::queue<vertex_type> children;
        finished = true;
        while(!parents.empty()) {
          vertex_type cur_parent = parents.front(); parents.pop();
          
          if((cur_parent == tree_indexer::null_vertex()) || (out_degree(cur_parent, m_tree) == 0)) {
            std::cout << " <empty> <empty>";
            children.push(tree_indexer::null_vertex());
            children.push(tree_indexer::null_vertex());
            continue; // this is a leaf node.
          };
          finished = false;
          child_vertex_iter vil, vi_end;
          boost::tie(vil, vi_end) = child_vertices(cur_parent, m_tree);
          children.push(*vil);
          std::cout << std::setw(8) << helper_type::value_to_key(m_tree[*vil]);
          if(out_degree(cur_parent, m_tree) == 1) {
            std::cout << " <empty>";
            children.push(tree_indexer::null_vertex());
            continue;
          };
          ++vil;
          children.push(*vil);
          std::cout << std::setw(8) << helper_type::value_to_key(m_tree[*vil]);
        };
        std::cout << std::endl;
        // swap parents for children:
        swap(parents,children);
      };
    };
#endif
    
    
    
    /* NOTE Invalidates vertices */
    /* NOTE This is a non-recursive version of the construct-node algorithm */
    /* Does not require persistent vertices */
    /* This is the main tree construction function. It takes the vertices in the iterator range and organizes them 
     * as a sub-tree below the aParentNode node. */
    void construct_node(vertex_type aCurrentNode,
                        typename std::vector<vertex_property>::iterator aBegin, 
                        typename std::vector<vertex_property>::iterator aEnd) {
      typedef typename std::vector<vertex_property>::iterator PropIter;
      using std::swap;
      
      std::queue<construction_task> tasks;
      tasks.push(construction_task(aCurrentNode, aBegin, aEnd));
      
      // breadth-first construction of the tree from aCurrentNode and everything below it.
      while(!tasks.empty()) {
        construction_task cur_task = tasks.front(); tasks.pop();
        
        // partition the elements around the median element:
        PropIter cur_task_median = cur_task.first + std::distance(cur_task.first, cur_task.last) / 2;
        std::nth_element(cur_task.first, cur_task_median, cur_task.last, m_compare);
        
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        m_tree[cur_task.node] = std::move(*cur_task_median);
#else
        m_tree[cur_task.node] = *cur_task_median;
#endif
        
        // first take care of the [first, median) range.
        if(cur_task.first == cur_task_median) 
          continue; // this is a leaf node.
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(cur_task.node, m_tree);
        tasks.push(construction_task(new_child.first, cur_task.first, cur_task_median));
        
        // then, take care of the (median, last) range.
        ++cur_task_median; // exclude the median.
        if(cur_task_median == cur_task.last)
          continue; // there is no right-child.
        new_child = add_child_vertex(cur_task.node, m_tree);
        tasks.push(construction_task(new_child.first, cur_task_median, cur_task.last));
        
      };
    };
    
    bool has_left_child(vertex_type u) const { return out_degree(u,m_tree) >= 1; };
    bool has_right_child(vertex_type u) const { return out_degree(u,m_tree) >= 2; };
    
    vertex_type get_left_child(vertex_type u) const {
      return *(child_vertices(u, m_tree).first);
    };
    vertex_type get_right_child(vertex_type u) const {
      return *(++(child_vertices(u, m_tree).first));
    };
    
    vertex_type find_lower_bound(const key_type& k, vertex_type cur_node) const {
      vertex_type cur_lb = tree_indexer::null_vertex();
      while(true) {
        if(!m_compare(m_tree[cur_node], k)) {
          cur_lb = cur_node;
          // descend on the left side.
          if( !has_left_child(cur_node) )
            break;
          cur_node = get_left_child(cur_node);
        } else {
          // descend on the right side.
          if( !has_right_child(cur_node) )
            break;
          cur_node = get_right_child(cur_node);
        };
      };
      return cur_lb;
    };
    
    vertex_type find_upper_bound(const key_type& k, vertex_type cur_node) const {
      vertex_type cur_ub = tree_indexer::null_vertex();
      while(true) {
        if(m_compare(k, m_tree[cur_node])) {
          cur_ub = cur_node;
          // descend on the left side.
          if( !has_left_child(cur_node) )
            break;
          cur_node = get_left_child(cur_node);
        } else {
          // descend on the right side.
          if( !has_right_child(cur_node) )
            break;
          cur_node = get_right_child(cur_node);
        };
      };
      return cur_ub;
    };
    
    std::pair< std::size_t, std::size_t> get_minmax_depth(vertex_type aNode) const {
      std::pair< std::size_t, std::size_t> result(std::numeric_limits<std::size_t>::max(),0);
      
      // assume that a breadth-first traversal is more efficient:
      std::queue< std::pair< vertex_type, std::size_t > > tasks;
      tasks.push(std::pair< vertex_type, std::size_t >(aNode, 0));
      while(!tasks.empty()) {
        std::pair< vertex_type, std::size_t > cur_task = tasks.front(); tasks.pop();
        if( !has_left_child(cur_task.first) ) {
          if(cur_task.second < result.first)
            result.first = cur_task.second;
          if(cur_task.second > result.second)
            result.second = cur_task.second;
          continue;
        };
        ++(cur_task.second);
        tasks.push(std::pair< vertex_type, std::size_t >(get_left_child(cur_task.first), cur_task.second));
        if( has_right_child(cur_task.first) )
          tasks.push(std::pair< vertex_type, std::size_t >(get_right_child(cur_task.first), cur_task.second));
      };
      
      return result;
    };
    
    // collect all the keys below and including the given root-node. Appends all the keys to the list.
    void collect_nodes(vertex_type aRootNode, std::vector< vertex_property >& aList) {
      // assume that a breadth-first traversal is more efficient:
      std::queue< vertex_type > tasks;
      tasks.push(aRootNode);
      while(!tasks.empty()) {
        vertex_type cur_task = tasks.front(); tasks.pop();
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        aList.push_back(std::move(m_tree[cur_task]));
#else
        aList.push_back(m_tree[cur_task]);
#endif
        if( !has_left_child(cur_task) )
          continue;
        tasks.push(get_left_child(cur_task));
        if( has_right_child(cur_task) )
          tasks.push(get_right_child(cur_task));
      };
      if( has_right_child(aRootNode) )
        remove_branch(get_right_child(aRootNode), m_tree);
      if( has_left_child(aRootNode) )
        remove_branch(get_left_child(aRootNode), m_tree);
    };
    
    // re-balance the sub-tree rooted at the given node. Uses a "collect and re-construct" approach.
    void rebalance_subtree(vertex_type aRootNode) {
      // collect and remove.
      std::vector< vertex_property > collected_nodes;
      collect_nodes(aRootNode, collected_nodes);
      // re-construct:
      construct_node(aRootNode, collected_nodes.begin(), collected_nodes.end());
    };
    
    // re-balance the sub-tree rooted at the given node. Uses a "collect and re-construct" approach.
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    void rebalance_and_insert_subtree(vertex_type aRootNode, vertex_property&& aAdded) {
#else
    void rebalance_and_insert_subtree(vertex_type aRootNode, const vertex_property& aAdded) {
#endif
      // collect and remove.
      std::vector< vertex_property > collected_nodes;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      collected_nodes.push_back(std::move(aAdded));
#else
      collected_nodes.push_back(aAdded);
#endif
      collect_nodes(aRootNode, collected_nodes);
      // re-construct:
      construct_node(aRootNode, collected_nodes.begin(), collected_nodes.end());
    };
    
    
    std::pair< vertex_type, bool > find_imbalance_impl(vertex_type u, int depth_change) const {
      // if we were to change the depth below u by the given depth-change, then would it 
      // create an imbalance somewhere? if yes, what vertex stems the sub-tree that must be re-balanced?
      std::pair< std::size_t, std::size_t> c_depth = get_minmax_depth(u);
      // assume c_depth.first == c_depth.second, i.e., the sub-tree is full and balanced (otherwise, what's the point of calling this function):
      c_depth.first  += depth_change;
      c_depth.second += depth_change;
      vertex_type orig_u = u;
      
      in_edge_iter ei, ei_end;
      boost::tie(ei, ei_end) = in_edges(u, m_tree);
      while(ei != ei_end) { // until you hit the root node.
        vertex_type p = source(*ei, m_tree);
        child_vertex_iter vil, vi_end;
        boost::tie(vil, vi_end) = child_vertices(p, m_tree);
        std::pair< std::size_t, std::size_t> o_depth;
        if(*vil != u)
          o_depth = get_minmax_depth(*vil);
        else if(out_degree(p, m_tree) > 1)
          o_depth = get_minmax_depth(*(++vil));
        else
          o_depth = std::pair< std::size_t, std::size_t>(0,0);
        
        // check if the other branch already has compatible depth-bounds:
        if((c_depth.first >= o_depth.first) && (c_depth.second <= o_depth.second)) 
          return std::pair< vertex_type, bool >(orig_u, false);
        
        // if not, then check if the other branch already has incompatible depth-bounds:
        if((o_depth.first < c_depth.second - 1) || 
           (o_depth.second > c_depth.first + 1)) // then, the sub-tree spanned from "p" will be unbalanced:
          return std::pair< vertex_type, bool >(p, true);
        
        // otherwise, keep on looking:
        if(o_depth.first < c_depth.first)
          c_depth.first = o_depth.first;
        if(o_depth.second > c_depth.second)
          c_depth.second = o_depth.second;
        c_depth.first  += 1;  // move up one level.
        c_depth.second += 1;
        u = p;
        boost::tie(ei, ei_end) = in_edges(u, m_tree);
      };
      return std::pair< vertex_type, bool >(orig_u, false); // could not find an imbalance anywhere.
    };
    
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< iterator, bool > insert_after_terminal_impl(vertex_type aBefore, vertex_property&& aValue) {
#else
    std::pair< iterator, bool > insert_after_terminal_impl(vertex_type aBefore, const vertex_property& aValue) {
#endif
      using std::swap;
      if(out_degree(aBefore, m_tree)) { // then the before node is terminal but not a leaf, safe to insert:
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(aBefore, std::move(aValue), m_tree);
#else
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(aBefore, aValue, m_tree);
#endif
        return std::pair< iterator, bool >(iterator(&m_tree, new_child.first), true);
      };
      // otherwise, a new child will be needed below the 'aBefore' node.
      // find an unbalanced sub-tree.
      std::pair< vertex_type, bool > imbal_u = find_imbalance_impl(aBefore, 1);
      if(imbal_u.second) {
        // must re-balance and insert at sub-tree from imbal_u.first:
        key_type tmp_k = helper_type::value_to_key(aValue);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        rebalance_and_insert_subtree(imbal_u.first, std::move(aValue));
#else
        rebalance_and_insert_subtree(imbal_u.first, aValue);
#endif
        return std::pair< iterator, bool >(iterator(&m_tree, find_lower_bound(tmp_k, imbal_u.first)), true);
      } else {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(aBefore, std::move(aValue), m_tree);
#else
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(aBefore, aValue, m_tree);
#endif
        swap(m_tree[new_child.first], m_tree[aBefore]);
        return std::pair< iterator, bool >(iterator(&m_tree, aBefore), true);
      };
    };
    
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_type insert_in_range_impl(iterator tvi, iterator tvi_end, vertex_property&& aValue) {
#else
    vertex_type insert_in_range_impl(iterator tvi, iterator tvi_end, const vertex_property& aValue) {
#endif
      // lets try to find a sweet spot (balance-wise) to place the new node:
      vertex_type tmp_u = tvi.base(); 
      for(; tvi != tvi_end; ++tvi) {
        if(out_degree(tvi.base(), m_tree) == 1) { // then, this is definitely a good spot:
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          std::pair<vertex_type, edge_type> new_child = add_child_vertex(tvi.base(), std::move(aValue), m_tree);
#else
          std::pair<vertex_type, edge_type> new_child = add_child_vertex(tvi.base(), aValue, m_tree);
#endif
          return new_child.first;
        };
      };
      // if we reach this point, then, try to go up to some re-balanceable node.
      boost::graph::detail::bst_traversal_status dummy_status = boost::graph::detail::OnRightBranch;
      boost::graph::detail::bst_move_up_to_next(m_tree, tmp_u, dummy_status);
      if(tmp_u != m_root) {
        std::pair< vertex_type, bool > imbal_u = find_imbalance_impl(tmp_u, 1);
        if(imbal_u.second) {
          // must re-balance and insert at sub-tree from imbal_u.first:
          key_type tmp_k = helper_type::value_to_key(aValue);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          rebalance_and_insert_subtree(imbal_u.first, std::move(aValue));
#else
          rebalance_and_insert_subtree(imbal_u.first, aValue);
#endif
          return find_lower_bound(tmp_k, imbal_u.first);
        } else {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          rebalance_and_insert_subtree(tmp_u, std::move(aValue));
#else
          rebalance_and_insert_subtree(tmp_u, aValue);
#endif
          return tmp_u;
        };
      } else
        return tree_indexer::null_vertex();
    };
    
    
    
    /* 
     * The vertex aBefore is the "lower-bound" (either the first equal element or the first greater element if value is not 
     * already in the tree).
     * The vertex aAfter is the "upper-bound" (the first greater element).
     * The value is the property of the vertex to be inserted.
     */
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< iterator, bool > insert_impl(vertex_type aBefore, vertex_type aAfter, vertex_property aValue) { // by-value.
#else
    std::pair< iterator, bool > insert_impl(vertex_type aBefore, vertex_type aAfter, const vertex_property& aValue) {
#endif
      using std::swap;
      
      if((aAfter == tree_indexer::null_vertex()) && (aBefore != tree_indexer::null_vertex())) {
        // this means that all the vertices remaining after aBefore are all equal to the given value.
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        vertex_type tmp_u = insert_in_range_impl(iterator(&m_tree, aBefore), iterator::end(&m_tree), std::move(aValue));
#else
        vertex_type tmp_u = insert_in_range_impl(iterator(&m_tree, aBefore), iterator::end(&m_tree), aValue);
#endif
        if(tmp_u == tree_indexer::null_vertex())
          aBefore = tmp_u; // this will cause it to be inserted at the end (see below).
        else
          return std::pair< iterator, bool >(iterator(&m_tree, tmp_u), true);
      };
      
      if((aBefore == tree_indexer::null_vertex()) && (aAfter == tree_indexer::null_vertex())) {
        // this means that the vertex must be appended to the end. (value is greater than all elements in the tree)
        aBefore = boost::graph::detail::bst_go_down_right(m_tree, m_root);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        return insert_after_terminal_impl(aBefore, std::move(aValue));
#else
        return insert_after_terminal_impl(aBefore, aValue);
#endif
      };
      
      
      // then, the special case where the equal range is empty.
      if(aBefore == aAfter) {
        if(out_degree(aAfter, m_tree) == 0) {
          // this means we could add a left-child to this node, see if it causes an imbalance:
          std::pair< vertex_type, bool > imbal_u = find_imbalance_impl(aAfter, 1);
          if(imbal_u.second) {
            // must re-balance and insert at sub-tree from imbal_u.first:
            key_type tmp_k = helper_type::value_to_key(aValue);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
            rebalance_and_insert_subtree(imbal_u.first, std::move(aValue));
#else
            rebalance_and_insert_subtree(imbal_u.first, aValue);
#endif
            return std::pair< iterator, bool >(iterator(&m_tree, find_lower_bound(tmp_k, imbal_u.first)), true);
          } else {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
            std::pair<vertex_type, edge_type> new_child = add_child_vertex(aAfter, std::move(aValue), m_tree);
#else
            std::pair<vertex_type, edge_type> new_child = add_child_vertex(aAfter, aValue, m_tree);
#endif
            return std::pair< iterator, bool >(iterator(&m_tree, new_child.first), true);
          };
        };
        
        if(out_degree(aAfter, m_tree) == 1) {
          // this means that the left-child must be a leaf (otherwise it would be imbalanced), so, we can rotate the tree.
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          std::pair<vertex_type, edge_type> new_child = add_child_vertex(aAfter, std::move(aValue), m_tree);
#else
          std::pair<vertex_type, edge_type> new_child = add_child_vertex(aAfter, aValue, m_tree);
#endif
//           swap(m_tree[new_child.first], m_tree[*(child_vertices(aAfter,m_tree).first)]);
          swap(m_tree[aAfter], m_tree[new_child.first]);
          return std::pair< iterator, bool >(iterator(&m_tree, *(child_vertices(aAfter,m_tree).first)), true);
        };
        
        // else, aAfter is full, lets try after the previous node (which must be a leaf or non-full node):
        aBefore = boost::graph::detail::bst_go_down_right(m_tree, *(child_vertices(aAfter,m_tree).first));
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        return insert_after_terminal_impl(aBefore, std::move(aValue));
#else
        return insert_after_terminal_impl(aBefore, aValue);
#endif
      };
      
      // then, the general case: the equal range is somewhere in the middle of the set.
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      vertex_type tmp_u = insert_in_range_impl(iterator(&m_tree, aBefore), iterator(&m_tree, aAfter), std::move(aValue));
#else
      vertex_type tmp_u = insert_in_range_impl(iterator(&m_tree, aBefore), iterator(&m_tree, aAfter), aValue);
#endif
      if(tmp_u == tree_indexer::null_vertex())
        return std::pair< iterator, bool >(iterator::end(&m_tree), false); // this case is impossible.
      else
        return std::pair< iterator, bool >(iterator(&m_tree, tmp_u), true);
      // at this point, the vertex must have been inserted one way or another (worst-case: the whole tree got re-constructed).
    };
    
    
    /*
     * Erase all vertices in the range [before, after)  (excl. after).
     * The vertex aBefore is the "lower-bound" (either the first equal element or the first greater element if value is not 
     * already in the tree).
     * The vertex aAfter is the "upper-bound" (the first greater element).
     */
    iterator erase_impl(iterator aBefore, iterator aAfter) {
      iterator befAfter = aAfter; --befAfter;
      if((aBefore == aAfter) || (m_root == boost::graph_traits<tree_indexer>::null_vertex()))
        return aBefore;
      
      key_type tmp_after_key;
      bool after_at_end = false;
      if(aAfter != iterator::end(&m_tree))
        tmp_after_key = helper_type::value_to_key(m_tree[aAfter.base()]);
      else
        after_at_end = true;
      
      // find the biggest possible sub-tree that contains all elements from the range.
      vertex_type mid_root = m_root;
      while(true) {
        if( m_compare(m_tree[mid_root], m_tree[aBefore.base()]) ) { // if root is less than lower-bound than go right.
          mid_root = get_right_child(mid_root);
          continue;
        };
        if( m_compare(m_tree[befAfter.base()], m_tree[mid_root]) ) { // if root is greater than upper-bound than go left.
          mid_root = get_left_child(mid_root);
          continue;
        };
        break; // if the root is the highest root that is neither less than 'before' nor greater than 'after'
      };
      
      // collect and remove.
      std::vector< vertex_property > collected_nodes;
      iterator it_near(&m_tree, boost::graph::detail::bst_go_down_left(m_tree, mid_root));
      iterator it_far(&m_tree, boost::graph::detail::bst_go_down_right(m_tree, mid_root));
      
      while(true) {
        if(it_near == aBefore) {
          if(it_far.base() == befAfter.base())
            break;
          it_near = aAfter;
          continue;
        };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        collected_nodes.push_back(std::move(m_tree[it_near.base()]));
#else
        collected_nodes.push_back(m_tree[it_near.base()]);
#endif
        if(it_near == it_far)
          break;
        ++it_near;
      };
      
      if( has_right_child(mid_root) )
        remove_branch(get_right_child(mid_root), m_tree);
      remove_branch(get_left_child(mid_root), m_tree);
      
      if((m_root == mid_root) && (collected_nodes.empty())) {
        remove_branch(m_root, m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
        return iterator::end(&m_tree);
      };
      
      if(collected_nodes.empty()) {
        in_edge_iter ei, ei_end;
        boost::tie(ei, ei_end) = in_edges(mid_root, m_tree);
        vertex_type tmp_parent = source(*ei, m_tree);
        remove_branch(mid_root, m_tree);
        mid_root = tmp_parent;
      } else {
        
        // re-construct:
        construct_node(mid_root, collected_nodes.begin(), collected_nodes.end());
      };
      
      while(true) {
        std::pair< vertex_type, bool > imbal_u = find_imbalance_impl(mid_root, 0);
        if(!imbal_u.second)
          break;
        rebalance_subtree(imbal_u.first);
        mid_root = imbal_u.first;
      };
      
      if(after_at_end)
        return iterator::end(&m_tree);
      else
        return iterator(&m_tree, find_lower_bound(tmp_after_key, mid_root));
    };
    
    
    
  public:
    
    /**
     * Creates a AVL-tree with no elements.
     */
    explicit avl_tree_impl(const allocator_type& = allocator_type()) : 
                           m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(Compare()) { };
    
    /**
     * Creates a AVL-tree with no elements.
     */
    explicit avl_tree_impl(const Compare& comp, const allocator_type& = allocator_type()) : 
                           m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(comp) { };
    
    /**
     * Builds a AVL-tree from an iterator range.
     * \param aFirst An input iterator (start of range).
     * \param aLast An input iterator (end of range).
     * \param comp A comparison functor.
     */
    template <typename InputIterator>
    avl_tree_impl(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(comp) {
      m_root = create_root(m_tree);
      std::vector<value_type> v(aFirst, aLast);
      construct_node(m_root, v.begin(), v.end());
    };
          
    /**
     * Builds a AVL-tree from an initializer_list.
     * \param aList An std::initializer_list.
     * \param comp A comparison functor.
     */  
    avl_tree_impl(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(comp) {
      m_root = create_root(m_tree);
      std::vector<value_type> v(aList);
      construct_node(m_root, v.begin(), v.end());
    };
    
    
    /**
     * Checks if the AVL-tree is empty.
     * \return True if the AVL-tree is empty.
     */
    bool empty() const
#ifndef BOOST_NO_CXX11_NOEXCEPT
    noexcept
#endif
    { return (num_vertices(m_tree) == 0); };
    /**
     * Returns the size of the AVL-tree (the number of vertices it contains).
     * \return The size of the AVL-tree (the number of vertices it contains).
     */
    size_type size() const
#ifndef BOOST_NO_CXX11_NOEXCEPT
    noexcept
#endif
    { return num_vertices(m_tree); };
    /**
     * Returns the maximum size of the AVL-tree.
     * \return The maximum size of the AVL-tree.
     */
    size_type max_size() const
#ifndef BOOST_NO_CXX11_NOEXCEPT
    noexcept
#endif
    { return std::numeric_limits<size_type>::max(); };
    
    /**
     * Standard swap function.
     */
    void swap(self& rhs) {
      using std::swap;
      swap(this->m_tree, rhs.m_tree);
      swap(this->m_root, rhs.m_root);
      swap(this->m_compare, rhs.m_compare);
    };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) {
      lhs.swap(rhs);
    };
    
    /**
     * Returns the depth of the tree.
     * \return The depth of the tree.
     */
    size_type depth() const { 
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return 0;
      else
        return get_minmax_depth(m_root).second; 
    };
    
    /** Returns the comparison object with which the tree was constructed.  */
    key_compare key_comp() const { return m_compare; };
    /** Returns the comparison object with which the tree was constructed.  */
    value_compare value_comp() const { return m_compare; };
    
    /** Returns the allocator object with which the tree was constructed.  */
    allocator_type get_allocator() const { return allocator_type(); };
    
    
    
    /**
     * Returns a read-only (constant) iterator that points to the first element in the set. 
     * Iteration is done in ascending order according to the keys.
     */
    const_iterator begin() const {
      return const_iterator::begin(&m_tree);
    };
    
    /**
     * Returns a read-only (constant) iterator that points one past the last element in the set. 
     * Iteration is done in ascending order according to the keys.
     */
    const_iterator end() const {
      return const_iterator::end(&m_tree);
    };
    
    /**
     * Returns a read-only (constant) iterator that points to the first element in the set. 
     * Iteration is done in ascending order according to the keys.
     */
    const_iterator cbegin() const { return begin(); };
    /**
     * Returns a read-only (constant) iterator that points one past the last element in the set. 
     * Iteration is done in ascending order according to the keys.
     */
    const_iterator cend() const { return end(); };
    
    
    /**
     * Returns a read-only (constant) reverse iterator that points to the last element in the set. 
     * Iteration is done in descending order according to the keys.
     */
    const_reverse_iterator rbegin() const {
      return const_reverse_iterator(const_iterator::end(&m_tree));
    };
    
    /**
     * Returns a read-only (constant) reverse iterator that points to the last element in the set. 
     * Iteration is done in descending order according to the keys.
     */
    const_reverse_iterator rend() const {
      return const_reverse_iterator(const_iterator::begin(&m_tree));
    };
    
    /**
     * Returns a read-only (constant) reverse iterator that points to the last element in the set. 
     * Iteration is done in descending order according to the keys.
     */
    const_reverse_iterator crbegin() const { return rbegin(); };
    /**
     * Returns a read-only (constant) reverse iterator that points to the last element in the set. 
     * Iteration is done in descending order according to the keys.
     */
    const_reverse_iterator crend() const { return rend(); };
    
    
    /**
     * Tries to locate a key in a set.
     * \param aKey Key to be located.
     * \return Iterator pointing to sought-after element, or end() if not found.
     */
    const_iterator find(const key_type& aKey) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return const_iterator::end(&m_tree);
      vertex_type u = find_lower_bound(aKey, m_root);
      if(!m_compare(aKey, m_tree[u]))
        return const_iterator(&m_tree, u);
      else
        return const_iterator::end(&m_tree);
    };
    
    /**
     * Finds the beginning of a subsequence matching given key.
     * \param aKey Key to be located.
     * \return Iterator pointing to first element equal to or greater than key, or end().
     */
    const_iterator lower_bound(const key_type& aKey) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return const_iterator::end(&m_tree);
      vertex_type u = find_lower_bound(aKey, m_root);
      if(u == tree_indexer::null_vertex())
        return end();
      return const_iterator(&m_tree, u);
    };
    
    /**
     * Finds the end of a subsequence matching given key.
     * \param aKey Key to be located.
     * \return Iterator pointing to the first element greater than key, or end().
     */
    const_iterator upper_bound(const key_type& aKey) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return const_iterator::end(&m_tree);
      vertex_type u = find_upper_bound(aKey, m_root);
      if(u == tree_indexer::null_vertex())
        return end();
      return const_iterator(&m_tree, u);
    };
    
    /**
     * Finds a subsequence matching given key.
     * \param aKey Key to be located.
     * \return Pair of iterators that possibly points to the subsequence matching given key.
     */
    std::pair<const_iterator,const_iterator> equal_range(const key_type& aKey) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return std::pair<const_iterator,const_iterator>(const_iterator::end(&m_tree),const_iterator::end(&m_tree));
      return std::pair<const_iterator,const_iterator>(lower_bound(aKey),upper_bound(aKey));
    };
    
    /**
     * Finds the number of keys.
     * \param aKey Key to be located.
     * \return Number of elements with specified key.
     */
    size_type count(const key_type& aKey) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        return 0;
      std::pair<const_iterator,const_iterator> er = equal_range(aKey);
      return std::distance(er.first, er.second);
    };
    
    
    
    
    /**
     * Attempts to insert an element into the set.
     * \param aValue Element to be inserted.
     * \return A pair, of which the first element is an iterator that points to the possibly inserted element, 
     *         and the second is a bool that is true if the element was actually inserted.
     */
    std::pair< iterator, bool > insert(const value_type& aValue) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex()) {
        m_root = create_root(aValue, m_tree);
        return std::pair<iterator, bool>(begin(), true);
      };
      key_type k = helper_type::value_to_key(aValue);
      return insert_impl(find_lower_bound(k, m_root), find_upper_bound(k, m_root), aValue);
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Attempts to insert an element into the set.
     * \param aValue Element to be inserted.
     * \return A pair, of which the first element is an iterator that points to the possibly inserted element, 
     *         and the second is a bool that is true if the element was actually inserted.
     */
    template <typename P>
    std::pair< iterator, bool > insert(P&& aValue) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex()) {
        m_root = create_root(value_type(std::forward<P>(aValue)), m_tree);
        return std::pair<iterator, bool>(begin(), true);
      };
      key_type k = helper_type::value_to_key(aValue);
      return insert_impl(find_lower_bound(k, m_root), find_upper_bound(k, m_root), value_type(std::forward<P>(aValue)));
    };
#endif
    
    /**
     * Attempts to insert an element into the set.
     * \param aPosition An iterator that serves as a hint as to where the element should be inserted.
     * \param aValue Element to be inserted.
     * \return An iterator that points to the element with key of x (may or may not be the element passed in).
     */
    iterator insert(const_iterator aPosition, const value_type& aValue) {
      if(aPosition == end())
        return insert(aValue).first;
      if(m_compare(*aPosition, aValue)) { // if position is less than value.
        const_iterator p2 = aPosition; ++p2;
        if(!m_compare(*p2, aValue)) { // if p2 if not less than value.
          // find first iterator that is greater than value:
          const_iterator p3 = p2;
          while((p3 != end()) && (!m_compare(aValue, *p3)))
            ++p3;
          return insert_impl(p2.base(), p3.base(), aValue).first;
        };
        // else, it means that the hint is wrong:
        return insert(aValue).first;
      };
      // else, position is not less than value.
      if(aPosition == begin())
        return insert_impl(aPosition.base(), aPosition.base(), aValue).first;
      const_iterator p4 = aPosition; --p4;
      if(!m_compare(*p4, aValue)) {
        // hint is wrong:
        return insert(aValue).first;
      };
      // else, p4 is less than value. Find first iterator that is greater than value:
      const_iterator p5 = aPosition;
      while((p5 != end()) && (!m_compare(aValue, *p5)))
        ++p5;
      return insert_impl(aPosition.base(), p5.base(), aValue).first;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Attempts to insert an element into the set.
     * \param aPosition An iterator that serves as a hint as to where the element should be inserted.
     * \param aValue Element to be inserted.
     * \return An iterator that points to the element with key of x (may or may not be the element passed in).
     */
    template <typename P>
    iterator insert(const_iterator aPosition, P&& aValue) {
      if(aPosition == end())
        return insert(std::forward<P>(aValue)).first;
      if(m_compare(*aPosition, aValue)) { // if position is less than value.
        const_iterator p2 = aPosition; ++p2;
        if(!m_compare(*p2, aValue)) { // if p2 if not less than value.
          // find first iterator that is greater than value:
          const_iterator p3 = p2;
          while((p3 != end()) && (!m_compare(aValue, *p3)))
            ++p3;
          return insert_impl(p2.base(), p3.base(), value_type(std::forward<P>(aValue))).first;
        };
        // else, it means that the hint is wrong:
        return insert(std::forward<P>(aValue)).first;
      };
      // else, position is not less than value.
      if(aPosition == begin())
        return insert_impl(aPosition.base(), aPosition.base(), value_type(std::forward<P>(aValue))).first;
      const_iterator p4 = aPosition; --p4;
      if(!m_compare(*p4, aValue)) {
        // hint is wrong:
        return insert(std::forward<P>(aValue)).first;
      };
      // else, p4 is less than value. Find first iterator that is greater than value:
      const_iterator p5 = aPosition;
      while((p5 != end()) && (!m_compare(aValue, *p5)))
        ++p5;
      return insert_impl(aPosition.base(), p5.base(), value_type(std::forward<P>(aValue))).first;
    };
#endif
    
    /**
     * A template function that attempts to insert a range of elements.
     * \tparam InputIterator An input-iterator type that can be dereferenced to a value-type rvalue (or rvalue-ref (C++11)).
     * \param aFirst Iterator pointing to the start of the range to be inserted.
     * \param aLast Iterator pointing to the end of the range.
     */
    template <typename InputIterator>
    void insert(InputIterator aFirst, InputIterator aLast) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex()) {
        m_root = create_root(m_tree);
        std::vector<value_type> v(aFirst, aLast);
        construct_node(m_root, v.begin(), v.end());
        return;
      };
      for(; aFirst != aLast; ++aFirst)
        insert(*aFirst);
    };
    
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
    /**
     * Attempts to insert a list of elements into the set.
     * \param aList A std::initializer_list of elements to be inserted.
     */
    void insert(std::initializer_list<value_type> aList) {
      insert(aList.begin(), aList.end());
    };
#endif
    
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
    /**
     * Attempts to emplace-create an element into the set.
     * \param args The constructor arguments required to create the element.
     * \return A pair, of which the first element is an iterator that points to the possibly inserted element, 
     *         and the second is a bool that is true if the element was actually inserted.
     * \note This function does not really create the element emplace since it must first form the element to find 
     *       its appropriate position in the set.
     */
    template <typename... Args>
    std::pair< iterator, bool > emplace(Args&&... args) {
      value_type tmp_val(std::forward<Args>(args)...);
      return insert(std::move(tmp_val));
    };
    
    /**
     * Attempts to emplace-create an element into the set.
     * \param aPosition An iterator that serves as a hint as to where the element should be inserted.
     * \param args The constructor arguments required to create the element.
     * \return A pair, of which the first element is an iterator that points to the possibly inserted element, 
     *         and the second is a bool that is true if the element was actually inserted.
     * \note This function does not really create the element emplace since it must first form the element to find 
     *       its appropriate position in the set.
     */
    template <typename... Args>
    iterator emplace_hint(const_iterator aPosition, Args&&... args) {
      value_type tmp_val(std::forward<Args>(args)...);
      return insert(aPosition, std::move(tmp_val));
    };
#endif
    
    /**
     * Erases a [first,last) range of elements from a set.
     * \param aFirst Iterator pointing to the start of the range to be erased.
     * \param aLast Iterator pointing to the end of the range to be erased.
     * \return New iterator to aLast.
     */
    iterator erase(const_iterator aFirst, const_iterator aLast) {
      if(aFirst == aLast)
        return aLast;
      return erase_impl(aFirst, aLast);
    };
    
    /**
     * Erases an element from a set.
     * \param aPosition An iterator pointing to the element to be erased.
     * \return An iterator pointing to the element immediately following position 
     *         prior to the element being erased. If no such element exists, end() is returned.
     */
    iterator erase(const_iterator aPosition) {
      const_iterator p2 = aPosition; ++p2;
      return erase(aPosition,p2);
    };
    
    /**
     * Erases elements according to the provided key.
     * \param aKey Key of element to be erased.
     * \return The number of elements erased.
     */
    size_type erase(const key_type& aKey) {
      std::pair<const_iterator,const_iterator> er = equal_range(aKey);
      size_type result = std::distance(er.first, er.second);
      if(result == 0)
        return 0;
      erase(er.first, er.second);
      return result;
    };
    
    /**
     * Erases all elements in a set.
     */
#ifndef BOOST_NO_CXX11_NOEXCEPT
    void clear() noexcept {
#else
    void clear() {
#endif
      if( num_vertices(m_tree) != 0 ) {
        remove_branch(m_root,m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
      };
    };
    
    
    /**
     * Subscript ( [] ) access to map data. Allows for easy lookup with the subscript ( [] ) operator. 
     * Returns data associated with the key specified in subscript. If the key does not exist, a pair 
     * with that key is created using default values, which is then returned.
     * \param k The key for which data should be retrieved.
     * \return A reference to the data of the (key,data) pair.
     * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
     */
    mapped_type& operator[](const key_type& k) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex()) {
        iterator it = insert(helper_type::keymap_to_value(k, mapped_type())).first;
        return helper_type::value_to_mapped(m_tree[it.base()]);
      };
      vertex_type u = find_lower_bound(k, m_root);
      if(!m_compare(k, m_tree[u]))
        return helper_type::value_to_mapped(m_tree[u]);
      else { // insert the key:
        iterator it = insert(const_iterator(&m_tree, u), helper_type::keymap_to_value(k, mapped_type()));
        return helper_type::value_to_mapped(m_tree[it.base()]);
      };
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Subscript ( [] ) access to map data. Allows for easy lookup with the subscript ( [] ) operator. 
     * Returns data associated with the key specified in subscript. If the key does not exist, a pair 
     * with that key is created using default values, which is then returned.
     * \param k The key for which data should be retrieved.
     * \return A reference to the data of the (key,data) pair.
     * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
     */
    mapped_type& operator[](key_type&& k) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex()) {
        iterator it = insert(helper_type::keymap_to_value(std::move(k), mapped_type())).first;
        return helper_type::value_to_mapped(m_tree[it.base()]);
      };
      vertex_type u = find_lower_bound(k, m_root);
      if(!m_compare(k, m_tree[u]))
        return helper_type::value_to_mapped(m_tree[u]);
      else { // insert the key:
        iterator it = insert(const_iterator(&m_tree, u), helper_type::keymap_to_value(std::move(k), mapped_type()));
        return helper_type::value_to_mapped(m_tree[it.base()]);
      };
    };
#endif
    
    /**
     * Access to map data.
     * \param k The key for which data should be retrieved.
     * \return A reference to the data whose key is equivalent to k, if such a data is present in the map.
     * \throw std::out_of_range If no such data is present.
     * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
     */
    mapped_type& at(const key_type& k) {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        throw std::out_of_range("The key does not match an element of the map!");
      vertex_type u = find_lower_bound(k, m_root);
      if(!m_compare(k, m_tree[u]))
        return helper_type::value_to_mapped(m_tree[u]);
      else 
        throw std::out_of_range("The key does not match an element of the map!");
    };
    
    /**
     * Access to map const data.
     * \param k The key for which data should be retrieved.
     * \return A const reference to the data whose key is equivalent to k, if such a data is present in the map.
     * \throw std::out_of_range If no such data is present.
     * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
     */
    const mapped_type& at(const key_type& k) const {
      if(m_root == boost::graph_traits<tree_indexer>::null_vertex())
        throw std::out_of_range("The key does not match an element of the map!");
      vertex_type u = find_lower_bound(k, m_root);
      if(!m_compare(k, m_tree[u]))
        return helper_type::value_to_mapped(m_tree[u]);
      else 
        throw std::out_of_range("The key does not match an element of the map!");
    };
    
    
    
    /* for private use only (useful to keep a map / set synchronized to another dynamic structure. */
    struct mutation_visitor {
      self* m_parent;
      
      mutation_visitor(self* aParent) : m_parent(aParent) { };
      
      void remove_vertex(vertex_type v, tree_indexer&) const {
        m_parent->erase(v);
      };
      
      void add_vertex(const vertex_property& vp, tree_indexer&) const {
        m_parent->insert(vp);
      };
      
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      void add_vertex(vertex_property&& vp, tree_indexer&) const {
        m_parent->insert(std::move(vp));
      };
#endif
    };
    
};





};

};


#endif


















