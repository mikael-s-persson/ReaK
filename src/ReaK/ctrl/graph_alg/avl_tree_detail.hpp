/**
 * \file avl_tree_detail.hpp
 * 
 * This library implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN) 
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which 
 * only requires a strict weak ordering (comparison function, analogous to less-than).
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

#include "base/defs.hpp"

#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include <unordered_map>
#include <vector>
#include <stack>

#include "graph_alg/tree_concepts.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_raw_property_graph.hpp"

#include "graph_alg/bst_inorder_iterator.hpp"

#include <iterator>


namespace ReaK {

namespace graph {


/**
 * This class implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN) 
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which 
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * \tparam TreeType The tree type to be used to store the entries of this AVL tree.
 * \tparam Compare The comparison functor type that can compare two elements. 
 */
template <typename TreeType,
          typename Compare>
class avl_tree_impl
{
  public:
    
    typedef avl_tree_impl<TreeType, Compare> self;
    
    typedef typename TreeType::vertex_bundled value_type;
    typedef Compare value_compare;
    typedef value_type key_type;
    typedef Compare key_compare;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    
    typedef bst_inorder_iterator<const TreeType, const value_type> iterator;
    typedef iterator const_iterator;
    
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::allocator<typename TreeType::vertex_bundled> allocator_type;
    
    
  private:
    
    typedef TreeType tree_indexer;
    
    typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<tree_indexer>::edge_descriptor edge_type;
    typedef typename boost::graph_traits<tree_indexer>::out_edge_iterator out_edge_iter;
    typedef typename boost::graph_traits<tree_indexer>::in_edge_iterator in_edge_iter;
    
    typedef typename tree_indexer::vertex_property_type vertex_property;
    typedef typename tree_indexer::edge_property_type edge_property;
    
    //typedef detail::compare_pair_first< distance_type, vertex_type, std::less< distance_type > > priority_compare_type;
    //typedef std::vector< std::pair< distance_type, vertex_type > > priority_queue_type;
    
    tree_indexer m_tree;   ///< Tree storage.
    vertex_type m_root;    ///< Root node of the tree.
    value_compare m_compare;  ///< The comparison functor.
    
    //non-copyable.
    avl_tree_impl(const self&);
    self& operator=(const self&); 
    
    
    struct construction_task {
      typedef typename std::vector<vertex_property>::iterator PropIter;
      
      vertex_type node;
      typename std::vector<vertex_property>::iterator first;
      typename std::vector<vertex_property>::iterator last;
      
      construction_task(vertex_type aNode, PropIter aBegin, PropIter aEnd) : 
        node(aNode), first(aBegin), last(aEnd) { };
    };
    
    
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
        
#ifdef RK_ENABLE_CXX0X_FEATURES
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
        std::pair<vertex_type, edge_type> new_child = add_child_vertex(cur_task.node, m_tree);
        tasks.push(construction_task(new_child.first, cur_task_median, cur_task.last));
        
      };
    };
    
    bool has_left_child(vertex_type u) const { return in_degree(u,m_tree) >= 1; };
    bool has_right_child(vertex_type u) const { return in_degree(u,m_tree) >= 2; };
    
    vertex_type get_left_child(vertex_type u) const {
      return *(child_vertices(u, m_tree).first);
    };
    vertex_type get_right_child(vertex_type u) const {
      return *(++(child_vertices(u, m_tree).first));
    };
    
    template <typename Key>
    vertex_type find_lower_bound(const Key& k, vertex_type cur_node) const {
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
    
    template <typename Key>
    vertex_type find_upper_bound(const Key& k, vertex_type cur_node) const {
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
#ifdef RK_ENABLE_CXX0X_FEATURES
        aList.emplace_back(std::move(m_tree[cur_task]));
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
      remove_branch(get_left_child(aRootNode), m_tree);
    };
    
    // collect all the keys below and including the given root-node. Appends all the keys to the list.
    void collect_nodes(vertex_type aRootNode, std::vector< vertex_property >& aList, vertex_type aExcluded) {
      // assume that a breadth-first traversal is more efficient:
      std::queue< vertex_type > tasks;
      tasks.push(aRootNode);
      while(!tasks.empty()) {
        vertex_type cur_task = tasks.front(); tasks.pop();
        if(cur_task != aExcluded) {
#ifdef RK_ENABLE_CXX0X_FEATURES
          aList.emplace_back(std::move(m_tree[cur_task]));
#else
          aList.push_back(m_tree[cur_task]);
#endif
        };
        if( !has_left_child(cur_task) )
          continue;
        tasks.push(get_left_child(cur_task));
        if( has_right_child(cur_task))
          tasks.push(get_right_child(cur_task));
      };
      if( has_right_child(aRootNode) )
        remove_branch(get_right_child(aRootNode), m_tree);
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
    void rebalance_and_remove_subtree(vertex_type aRootNode, vertex_type aExcluded) {
      // collect and remove.
      std::vector< vertex_property > collected_nodes;
      collect_nodes(aRootNode, collected_nodes, aExcluded);
      // re-construct:
      construct_node(aRootNode, collected_nodes.begin(), collected_nodes.end());
    };
    
    // re-balance the sub-tree rooted at the given node. Uses a "collect and re-construct" approach.
#ifdef RK_ENABLE_CXX0X_FEATURES
    void rebalance_and_insert_subtree(vertex_type aRootNode, vertex_property&& aAdded) {
#else
    void rebalance_and_insert_subtree(vertex_type aRootNode, const vertex_property& aAdded) {
#endif
      // collect and remove.
      std::vector< vertex_property > collected_nodes;
#ifdef RK_ENABLE_CXX0X_FEATURES
      collected_nodes.emplace_back(std::move(aAdded));
#else
      collected_nodes.push_back(aAdded);
#endif
      collect_nodes(aRootNode, collected_nodes);
      // re-construct:
      construct_node(aRootNode, collected_nodes.begin(), collected_nodes.end());
    };
    
    
    
  public:
    
    /**
     * Creates a AVL-tree with no elements.
     */
    avl_tree_impl(const value_compare& comp = value_compare(), const allocator_type& = allocator_type()) : 
                  m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(comp) { };
    
    /**
     * Builds a AVL-tree from an iterator range.
     * \param aFirst An input iterator (start of range).
     * \param aLast An input iterator (end of range).
     * \param comp A comparison functor.
     */
    template <typename InputIterator>
    avl_tree_impl(InputIterator aFirst, InputIterator aLast, const value_compare& comp = value_compare(), const allocator_type& = allocator_type()) : 
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
    avl_tree_impl(std::initializer_list< value_type > aList, const value_compare& comp = value_compare(), const allocator_type& = allocator_type()) : 
                  m_tree(), m_root(boost::graph_traits<tree_indexer>::null_vertex()), m_compare(comp) {
      m_root = create_root(m_tree);
      std::vector<value_type> v(aList);
      construct_node(m_root, v.begin(), v.end());
    };
    
    
    /**
     * Checks if the AVL-tree is empty.
     * \return True if the AVL-tree is empty.
     */
    bool empty() const { return (num_vertices(m_tree) == 0); };
    /**
     * Returns the size of the AVL-tree (the number of vertices it contains).
     * \return The size of the AVL-tree (the number of vertices it contains).
     */
    size_type size() const { return num_vertices(m_tree); };
    /**
     * Returns the maximum size of the AVL-tree.
     * \return The maximum size of the AVL-tree.
     */
    size_type max_size() const { return std::numeric_limits<size_type>::max(); };
    
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
    size_type depth() const { return get_minmax_depth(m_root).second; };
    
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
      vertex_type u = find_lower_bound(aKey, m_root);
      return const_iterator(&m_tree, u);
    };
    
    /**
     * Finds the end of a subsequence matching given key.
     * \param aKey Key to be located.
     * \return Iterator pointing to the first element greater than key, or end().
     */
    const_iterator upper_bound(const key_type& aKey) const {
      vertex_type u = find_upper_bound(aKey, m_root);
      return const_iterator(&m_tree, u);
    };
    
    /**
     * Finds a subsequence matching given key.
     * \param aKey Key to be located.
     * \return Pair of iterators that possibly points to the subsequence matching given key.
     */
    std::pair<const_iterator,const_iterator> equal_range(const key_type& aKey) const {
      return std::pair<const_iterator,const_iterator>(lower_bound(aKey),upper_bound(aKey));
    };
    
    /**
     * Finds the number of keys.
     * \param aKey Key to be located.
     * \return Number of elements with specified key.
     */
    size_type count(const key_type& aKey) const {
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
      
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Attempts to insert an element into the set.
     * \param aValue Element to be inserted.
     * \return A pair, of which the first element is an iterator that points to the possibly inserted element, 
     *         and the second is a bool that is true if the element was actually inserted.
     */
    std::pair< iterator, bool > insert(value_type&& aValue) {
      
    };
#endif
    
    /**
     * Attempts to insert an element into the set.
     * \param aPosition An iterator that serves as a hint as to where the element should be inserted.
     * \param aValue Element to be inserted.
     * \return An iterator that points to the element with key of x (may or may not be the element passed in).
     */
    iterator insert(const_iterator aPosition, const value_type& aValue) {
      
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Attempts to insert an element into the set.
     * \param aPosition An iterator that serves as a hint as to where the element should be inserted.
     * \param aValue Element to be inserted.
     * \return An iterator that points to the element with key of x (may or may not be the element passed in).
     */
    iterator insert(const_iterator aPosition, value_type&& aValue) {
      
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
      
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Attempts to insert a list of elements into the set.
     * \param aList A std::initializer_list of elements to be inserted.
     */
    void insert(std::initializer_list<value_type> aList) {
      insert(aList.begin(), aList.end());
    };
    
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
     * Erases an element from a set.
     * \param aPosition An iterator pointing to the element to be erased.
     * \return An iterator pointing to the element immediately following position 
     *         prior to the element being erased. If no such element exists, end() is returned.
     */
    iterator erase(const_iterator aPosition) {
      
    };
    
    /**
     * Erases a [first,last) range of elements from a set.
     * \param aFirst Iterator pointing to the start of the range to be erased.
     * \param aLast Iterator pointing to the end of the range to be erased.
     * \return New iterator to aLast.
     */
    iterator erase(const_iterator aFirst, const_iterator aLast) {
      
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
#ifdef RK_ENABLE_CXX0X_FEATURES
    void clear() noexcept {
#else
    void clear() {
#endif
      if( num_vertices(m_tree) == 0 ) {
        remove_branch(m_root,m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
      };
    };
    
    
    
    
    
#if 0
    // keeping the code below just for quick reference.
    
    /**
     * Inserts a vertex into the tree.
     * \param up The vertex-property to be added to the DVP-tree.
     */
#ifdef RK_ENABLE_CXX0X_FEATURES
    void insert(vertex_property up) {
#else
    void insert(const vertex_property& up) {
#endif
      if(num_vertices(m_tree) == 0) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        m_root = create_root(std::move(up), m_tree); 
#else
        m_root = create_root(up, m_tree); 
#endif
        return;
      };
      point_type u_pt = get(m_position, up); 
      vertex_type u_realleaf = get_leaf(u_pt,m_root);
      if(u_realleaf == m_root) { //if the root is the leaf, it requires special attention since no parent exists.
        std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
        prop_list.push_back(std::move(up));
#else
        prop_list.push_back(up);
#endif
        remove_branch(u_realleaf, back_inserter(prop_list), m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
        u_realleaf = m_root;
        construct_node(u_realleaf, 0.0, prop_list.begin(), prop_list.end()); 
        return;
      };
      vertex_type u_leaf = source(*(in_edges(u_realleaf,m_tree).first),m_tree);
      if((out_degree(u_leaf,m_tree) < Arity) || (!is_leaf_node(u_leaf))) {
        // leaf node is not full of children, an additional child can be added 
        //  (must be reconstructed to keep ordering, but this is a trivial operation O(Arity)).
        //OR 
        // if leaf is not really a leaf, then it means that this sub-tree is definitely not balanced and not full either,
        //  then all the Keys ought to be collected and u_leaf ought to be reconstructed.
        update_mu_upwards(u_pt,u_leaf);
        distance_type e_dist = 0.0;
        vertex_type u_leaf_parent;
        if(u_leaf != m_root) {
          e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(u_leaf,m_tree).first)));
          u_leaf_parent = source(*(in_edges(u_leaf,m_tree).first),m_tree);
        } else 
          u_leaf_parent = boost::graph_traits<tree_indexer>::null_vertex();
        std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
        prop_list.push_back(std::move(up));
#else
        prop_list.push_back(up);
#endif
        remove_branch(u_leaf, back_inserter(prop_list), m_tree);
        construct_node(u_leaf_parent, e_dist, prop_list.begin(), prop_list.end()); 
      } else {
        //if it is a full-leaf, then this is a leaf node, and it is balanced but full, 
        // we should then find a non-full parent.
        vertex_type p = u_leaf;   
        int actual_depth_limit = 1;
        int last_depth_limit = actual_depth_limit;
        while((p != m_root) && (is_node_full(p,last_depth_limit))) {
          p = source(*(in_edges(p,m_tree).first),m_tree);
          last_depth_limit = ++actual_depth_limit;
        };
        bool is_p_full = false; 
        if(p == m_root)
          is_p_full = is_node_full(p,last_depth_limit);
        if((!is_p_full) && (last_depth_limit >= 0)) {
          //this means that we can add our key to the sub-tree of p and reconstruct from there.
          update_mu_upwards(u_pt,p);
          distance_type e_dist = 0.0;
          vertex_type p_parent;
          if(p != m_root) {
            e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(p,m_tree).first)));
            p_parent = source(*(in_edges(p,m_tree).first),m_tree);
          } else {
            p_parent = boost::graph_traits<tree_indexer>::null_vertex();
          };
          std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
          prop_list.push_back(std::move(up));
#else
          prop_list.push_back(up);
#endif
          remove_branch(p, back_inserter(prop_list), m_tree);
          construct_node(p_parent, e_dist, prop_list.begin(), prop_list.end());
        } else {
          //this means that either the root node is full or there are branches of the tree that are deeper than u_realleaf, 
          // and thus, in either case, u_realleaf should be expanded.
          edge_type l_p;
          edge_property ep;
          put(m_mu, ep, m_distance(u_pt, get(m_position, get(boost::vertex_raw_property,m_tree,u_realleaf)), *m_space));
#ifdef RK_ENABLE_CXX0X_FEATURES
          boost::tie(p, l_p) = add_child_vertex(u_realleaf, std::move(up), std::move(ep), m_tree);
#else
          boost::tie(p, l_p) = add_child_vertex(u_realleaf, up, ep, m_tree);
#endif
          update_mu_upwards(u_pt,u_realleaf);
        };
      };
    };
    /**
     * Inserts a range of vertices.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     */
    template <typename ForwardIterator>
    void insert(ForwardIterator aBegin, ForwardIterator aEnd) { 
      std::for_each(aBegin,aEnd,boost::bind(&self::insert,this,_1));
      //TODO: There's got to be a better way to insert many elements (most likely a similar strategy to the erase multiple function).
    };
    
    
    /**
     * Erases the given vertex from the DVP-tree.
     * \param u_node The vertex to be removed from the DVP-tree.
     */
    void erase(vertex_type u_node) { 
      if(num_vertices(m_tree) == 0) 
        return;
      if( (u_node == m_root) && (num_vertices(m_tree) == 1) ) {
        std::vector<vertex_property> prop_list;
        remove_branch(m_root, back_inserter(prop_list), m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
        return;
      };
      distance_type e_dist = 0.0;
      vertex_type u_parent = boost::graph_traits<tree_indexer>::null_vertex();
      if(u_node != m_root) {
        e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(u_node,m_tree).first)));
        u_parent = source(*(in_edges(u_node,m_tree).first), m_tree);
      };
      
      out_edge_iter ei, ei_end;
      std::vector<vertex_property> prop_list;
      if( (out_degree(u_node, m_tree) > 0) ||
          (u_parent == boost::graph_traits<tree_indexer>::null_vertex()) ) {
        remove_branch(u_node, back_inserter(prop_list), m_tree);
      } else {
        remove_branch(u_node, back_inserter(prop_list), m_tree);
        u_node = u_parent;
        if(u_parent == m_root)
          u_parent = boost::graph_traits<tree_indexer>::null_vertex();
        else
          u_parent = source(*(in_edges(u_node,m_tree).first), m_tree);
        remove_branch(u_node, back_inserter(prop_list), m_tree);
      };
      construct_node(u_parent, e_dist, prop_list.begin() + 1 /* skip first node (u_node) */, prop_list.end());
    };
    
    /**
     * Erases the given key-value and position from the DVP-tree.
     * Note that this function is not recommended if the vertex descriptor is known, as it will require a key-lookup.
     * \param u_key The key-value to be removed from the DVP-tree.
     * \param u_pt The position-value corresponding to the key-value.
     */
    void erase(key_type u_key, const point_type& u_pt) {
      vertex_type u_node;
      try {
        u_node = get_vertex(u_key, u_pt, m_root);
      } catch (int err) {
        return;
      };
      erase(u_node);
    };
    
    /**
     * Erases the given vertex-range from the DVP-tree.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices (by tree vertex descriptors).
     * \param aBegin The start of the range from which to take the vertices to be erased.
     * \param aEnd The end of the range from which to take the vertices to be erased (one-past-last).
     */
    template <typename ForwardIterator>
    void erase(ForwardIterator aBegin, ForwardIterator aEnd) { 
      if(num_vertices(m_tree) == 0) return;
      
      typedef std::list< std::pair< vertex_type,   // the node at the trunk of the sub-tree to be re-balanced.
                                    std::vector<vertex_type>  // the list of nodes below the re-balanced trunk.
                                  > > vertex_listing;
      vertex_listing v_lists; //will hold a list of unique nodes and all their non-erased 
      
      // First, generate the vertex-listings in preparation for the deletion.
      for(ForwardIterator first = aBegin; first != aEnd; ++first) {
        vertex_type removal_trunk = *first;
        if(out_degree(*first, m_tree) == 0)
          removal_trunk = source(*(in_edges(*first,m_tree).first), m_tree);
        put(m_key, get(boost::vertex_raw_property,m_tree,*first), reinterpret_cast<key_type>(-1)); // mark as invalid, for deletion.
        
        bool already_collected = false;
        for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ++it) {
          // is removal_trunk contained in the *it listing?
          typename std::vector<vertex_type>::iterator it_key = std::binary_search(it->second.begin(), it->second.end(), removal_trunk);
          if( it_key != it->second.end() ) {
            already_collected = true;
            break;
          };
        };
        
        if(!already_collected) {
          std::vector<vertex_type> removal_list;
          removal_list.push_back(removal_trunk);
          collect_vertices(removal_list, removal_trunk);
          std::sort(removal_list.begin(), removal_list.end());
          
          for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ) {
            // is the *it trunk contained in the removal_trunk?
            typename std::vector<vertex_type>::iterator it_key = std::binary_search(removal_list.begin(), removal_list.end(), it->first);
            if( it_key != removal_list.end() )
              it = v_lists.erase(it);
            else
              ++it;
          };
          
          v_lists.push_back( std::make_pair(removal_trunk, std::vector<vertex_type>()) );
          v_lists.back().second.swap(removal_list);
        };
        
      };
      
      for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ++it) {
        
        distance_type e_dist = 0.0;
        vertex_type u_parent = boost::graph_traits<tree_indexer>::null_vertex();
        if(it->first != m_root) {
          e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(it->first,m_tree).first)));
          u_parent = source(*(in_edges(it->first, m_tree).first), m_tree);
        };
      
        out_edge_iter ei, ei_end;
        std::vector<vertex_property> prop_list;
        remove_branch(it->first, back_inserter(prop_list), m_tree);
        prop_list.erase( remove_if(prop_list.begin(), prop_list.end(), boost::bind(is_vertex_prop_valid, m_key, _1)), prop_list.end());
        construct_node(u_parent, e_dist, prop_list.begin(), prop_list.end());
      };      
    };
    
    
#endif
    
    
    struct mutation_visitor {
      self* m_parent;
      
      mutation_visitor(self* aParent) : m_parent(aParent) { };
      
      void remove_vertex(vertex_type v, tree_indexer&) const {
        m_parent->erase(v);
      };
      
      void add_vertex(const vertex_property& vp, tree_indexer&) const {
        m_parent->insert(vp);
      };
      
#ifdef RK_ENABLE_CXX0X_FEATURES
      void add_vertex(vertex_property&& vp, tree_indexer&) const {
        m_parent->insert(std::move(vp));
      };
#endif
    };
    
};





};

};


#endif


















