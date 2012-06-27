/**
 * \file linked_tree.hpp
 * 
 * This library provides a class that implements a linked tree structure. This is a 
 * classic tree implementation in which each node contain a list of edges to its children.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2012
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

#ifndef REAK_LINKED_TREE_HPP
#define REAK_LINKED_TREE_HPP


#include <vector>
#include <stdexcept>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>

#include "bgl_tree_adaptor.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

namespace ReaK {

namespace graph {
  

/**
 * This traits class template is used to obtain the types (and meta-values) that describe 
 * the basic types used in a linked-tree with the given out-edge-list and vertex-list storage 
 * policies. This traits class is useful to obtain type-erased (or type-agnostic) vertex and edge
 * descriptors. Note, this traits class is essentially the linked-tree equivalent of the BGL adjacency_list_traits class.
 */
template <typename OutEdgeListS = boost::vecS, 
          typename VertexListS = boost::vecS>
struct linked_tree_traits {
  /** This meta-value tells if the vertex storage is random-access, or not. */
  typedef typename boost::detail::is_random_access<VertexListS>::type is_rand_access;
  /** This meta-value tells if the edges are bidirectional, or not. */
  typedef boost::mpl::true_ is_bidir;
  /** This meta-value tells if the edges are directional, or not. */
  typedef boost::mpl::true_ is_directed;
  
  /** This tag gives the edges' directional categorization. */
  typedef boost::bidirectional_tag directed_category;
  
  /** This meta-value tells if the parallel edges are allowed, or not. */
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  
  /** This type is used to describe the number of vertices. */
  typedef std::size_t vertices_size_type;
  typedef void* vertex_ptr;
  /** This type is used to describe a vertex in the tree. */
  typedef typename boost::mpl::if_< is_rand_access,
    vertices_size_type, vertex_ptr>::type vertex_descriptor;
  
  /** This type is used to describe an edge in the tree. */
  struct edge_descriptor {
    vertex_descriptor source;
    void* edge_id;
    
    edge_descriptor(vertex_descriptor aSrc, void* aEdgeId) : source(aSrc), edge_id(aEdgeId) { };
  };
  /** This type is used to describe the number of edges. */
  typedef std::size_t edges_size_type;
  
};



/**
 * This class implements a D-Ary Breadth-first tree that is tailored 
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
 * \tparam OutEdgeList A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexList A type tag to choose the storage policy for the vertices.
 * \tparam VertexProperties A POD type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A POD type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS,
          typename VertexListS,
          typename VertexProperties,
          typename EdgeProperties = boost::no_property >
class linked_tree
{
  public:
    typedef linked_tree<OutEdgeListS, VertexListS, VertexProperties, EdgeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::directed_category directed_category;
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::edge_parallel_category edge_parallel_category;
    struct traversal_category : 
      boost::bidirectional_graph_tag, 
      boost::vertex_list_graph_tag { };
    
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::vertex_descriptor vertex_descriptor;
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::edge_descriptor edge_descriptor;
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::vertices_size_type vertices_size_type;
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::edges_size_type edges_size_type;
    typedef edges_size_type degree_size_type;
    
    struct edge_value_type {
      vertex_descriptor target;
      edge_property_type data;
      
      edge_value_type(vertex_descriptor aTarget, const edge_property_type& aData) : target(aTarget), data(aData) { };
#ifdef RK_ENABLE_CXX0X_FEATURES
      edge_value_type(vertex_descriptor aTarget, edge_property_type&& aData) : target(aTarget), data(std::move(aData)) { };
#endif
    };
    typedef typename boost::container_gen<OutEdgeListS, edge_value_type>::type out_edge_container;
    
    struct vertex_value_type {
      vertex_property_type data;
      out_edge_container out_edges;
      edge_descriptor in_edge;
      
      vertex_value_type(const vertex_property_type& aData) : data(aData), out_edges(), in_edge(NULL) { };
#ifdef RK_ENABLE_CXX0X_FEATURES
      vertex_value_type(vertex_property_type&& aData) : data(std::move(aData)), out_edges(), in_edge(NULL) { };
#endif
    };
    typedef typename boost::container_gen<VertexListS, vertex_value_type>::type vertex_container;
    
    typedef typename container_type::difference_type vertex_descriptor;
    
    
    struct out_edge_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef edge_descriptor value_type;
      typedef edge_descriptor* pointer;
      typedef edge_descriptor reference;
      typedef typename out_edge_container::iterator edge_iter_impl;
      typedef typename std::iterator_traits<edge_iter_impl>::iterator_category iterator_category;
      
      vertex_descriptor source;
      edge_iter_impl edge_it;
      out_edge_iterator(vertex_descriptor aSrc, edge_iter_impl aEdgeIt) : source(aSrc), edge_it(aEdgeIt) { };
      
      friend bool operator==( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return (lhs.source == rhs.source) && (lhs.edge_it == rhs.edge_it); };
      friend bool operator!=( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return (lhs.source != rhs.source) || (lhs.edge_it != rhs.edge_it); };
      friend bool operator >( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { 
	return (lhs.source > rhs.source) ||
              ((lhs.source == rhs.source) && (lhs.edge_it > rhs.edge_it)); 
      };
      friend bool operator >=(const out_edge_iterator& lhs, const out_edge_iterator& rhs) { 
	return (lhs.source > rhs.source) ||
              ((lhs.source == rhs.source) && (lhs.edge_it >= rhs.edge_it)); 
      };
      friend bool operator <( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { 
	return (lhs.source < rhs.source) ||
              ((lhs.source == rhs.source) && (lhs.edge_it < rhs.edge_it)); 
      };
      friend bool operator <=(const out_edge_iterator& lhs, const out_edge_iterator& rhs) { 
	return (lhs.source < rhs.source) ||
              ((lhs.source == rhs.source) && (lhs.edge_it <= rhs.edge_it)); 
      };
      
      out_edge_iterator& operator++() { ++edge_it; return *this; };
      out_edge_iterator operator++(int) { out_edge_iterator result(*this); ++edge_it; return result; };
      out_edge_iterator& operator--() { --edge_it; return *this; };
      out_edge_iterator operator--(int) { out_edge_iterator result(*this); --edge_it; return result; };
      
      friend out_edge_iterator operator+(const out_edge_iterator& lhs, difference_type i) {
	return out_edge_iterator(lhs.source, std::advance(lhs.edge_it,i));
      };
      friend out_edge_iterator operator+(difference_type i, const out_edge_iterator& rhs) {
	return out_edge_iterator(rhs.source, std::advance(rhs.edge_it,i));
      };
      friend out_edge_iterator operator-(const out_edge_iterator& lhs, difference_type i) {
	return out_edge_iterator(lhs.source, std::advance(lhs.edge_it,-i));
      };
      friend difference_type operator-(const out_edge_iterator& lhs, const out_edge_iterator& rhs) {
        if(lhs.source == rhs.source)
	  return std::distance(rhs.edge_it, lhs.edge_it);
	throw invalid_argument("Attempted to compute the difference between out-edge iterators from different vertices!");
      };
      
      out_edge_iterator& operator +=(difference_type i) { edge_it = std::advance(edge_it,i); return *this; };
      out_edge_iterator& operator -=(difference_type i) { edge_it = std::advance(edge_it,-i); return *this; };
      
      value_type operator[](difference_type i) const { return edge_descriptor(source, &(*(*this + i))); };
      reference operator*() { return edge_descriptor(source, &(*edge_it)); };
    };
    
    typedef edge_descriptor* in_edge_iterator;
    typedef edge_descriptor* edge_iterator;
    
    struct vertex_iterator_with_iter {
      typedef std::ptrdiff_t difference_type;
      typedef void* value_type;
      typedef void** pointer;
      typedef void* reference;
      typedef typename vertex_container::iterator vertex_iter_impl;
      typedef typename std::iterator_traits<out_edge_iterator>::iterator_category iterator_category;
      
      vertex_iter_impl v_it;
      
      explicit vertex_iterator_with_iter(vertex_iter_impl aVIt) : v_it(aVIt) { };
      
      friend bool operator==( const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { return (lhs.v_it == rhs.v_it); };
      friend bool operator!=( const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { return (lhs.v_it != rhs.v_it); };
      friend bool operator >( const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { 
	return (lhs.v_it > rhs.v_it);
      };
      friend bool operator >=(const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { 
	return (lhs.v_it >= rhs.v_it);
      };
      friend bool operator <( const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { 
	return (lhs.v_it < rhs.v_it);
      };
      friend bool operator <=(const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) { 
	return (lhs.v_it <= rhs.v_it);
      };
      
      vertex_iterator_with_iter& operator++() { ++v_it; return *this; };
      vertex_iterator_with_iter operator++(int) { vertex_iterator_with_iter result(*this); ++v_it; return result; };
      vertex_iterator_with_iter& operator--() { --v_it; return *this; };
      vertex_iterator_with_iter operator--(int) { vertex_iterator_with_iter result(*this); --v_it; return result; };
      
      friend vertex_iterator_with_iter operator+(const vertex_iterator_with_iter& lhs, difference_type i) {
	return vertex_iterator_with_iter(*v_list, std::advance(lhs.v_it,i));
      };
      friend vertex_iterator_with_iter operator+(difference_type i, const vertex_iterator_with_iter& rhs) {
	return vertex_iterator_with_iter(*v_list, std::advance(rhs.v_it,i));
      };
      friend vertex_iterator_with_iter operator-(const vertex_iterator_with_iter& lhs, difference_type i) {
	return vertex_iterator_with_iter(*v_list, std::advance(lhs.v_it,-i));
      };
      friend difference_type operator-(const vertex_iterator_with_iter& lhs, const vertex_iterator_with_iter& rhs) {
        return std::distance(rhs.v_it, lhs.v_it);
      };
      
      vertex_iterator_with_iter& operator +=(difference_type i) { v_it = std::advance(v_it,i); return *this; };
      vertex_iterator_with_iter& operator -=(difference_type i) { v_it = std::advance(v_it,-i); return *this; };
      
      value_type operator[](difference_type i) const { 
	return &(*std::advance(v_it,i));
      };
      reference operator*() { 
	return &(*v_it);
      };
      
      static vertex_iterator_with_iter begin(vertex_container& c) { return vertex_iterator_with_iter(c.begin()); };
      static vertex_iterator_with_iter end(vertex_container& c) { return vertex_iterator_with_iter(c.end()); };
    };
    
    struct vertex_iterator_with_index {
      typedef std::ptrdiff_t difference_type;
      typedef std::size_t value_type;
      typedef std::size_t* pointer;
      typedef std::size_t reference;
      typedef typename std::iterator_traits<out_edge_iterator>::iterator_category iterator_category;
      
      std::size_t v_it;
      
      explicit vertex_iterator_with_index(std::size_t aVIt) : v_it(aVIt) { };
      
      friend bool operator==( const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { return (lhs.v_it == rhs.v_it); };
      friend bool operator!=( const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { return (lhs.v_it != rhs.v_it); };
      friend bool operator >( const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { 
	return (lhs.v_it > rhs.v_it);
      };
      friend bool operator >=(const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { 
	return (lhs.v_it >= rhs.v_it);
      };
      friend bool operator <( const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { 
	return (lhs.v_it < rhs.v_it);
      };
      friend bool operator <=(const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) { 
	return (lhs.v_it <= rhs.v_it);
      };
      
      vertex_iterator_with_index& operator++() { ++v_it; return *this; };
      vertex_iterator_with_index operator++(int) { vertex_iterator_with_index result(*this); ++v_it; return result; };
      vertex_iterator_with_index& operator--() { --v_it; return *this; };
      vertex_iterator_with_index operator--(int) { vertex_iterator_with_index result(*this); --v_it; return result; };
      
      friend vertex_iterator_with_index operator+(const vertex_iterator_with_index& lhs, difference_type i) {
	return vertex_iterator_with_index(lhs.v_it + i);
      };
      friend vertex_iterator_with_index operator+(difference_type i, const vertex_iterator_with_index& rhs) {
	return vertex_iterator_with_index(rhs.v_it + i);
      };
      friend vertex_iterator_with_index operator-(const vertex_iterator_with_index& lhs, difference_type i) {
	return vertex_iterator_with_index(lhs.v_it - i);
      };
      friend difference_type operator-(const vertex_iterator_with_index& lhs, const vertex_iterator_with_index& rhs) {
        return lhs.v_it - rhs.v_it;
      };
      
      vertex_iterator_with_index& operator +=(difference_type i) { v_it += i; return *this; };
      vertex_iterator_with_index& operator -=(difference_type i) { v_it -= i; return *this; };
      
      value_type operator[](difference_type i) const { 
	return v_it + i;
      };
      reference operator*() { 
	return v_it;
      };
      
      static vertex_iterator_with_index begin(vertex_container&) { return vertex_iterator_with_index(0); };
      static vertex_iterator_with_index end(vertex_container& c) { return vertex_iterator_with_index(c.size()); };
    };
    
    typedef typename linked_tree_traits<OutEdgeListS, VertexListS>::is_rand_access vertex_rand_access;
    
    typedef typename boost::mpl::if_<
      typename linked_tree_traits<OutEdgeListS, VertexListS>::is_rand_access,
      vertex_iterator_with_index,
      vertex_iterator_with_iter >::type vertex_iterator;
    
    struct child_vertex_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef vertex_descriptor value_type;
      typedef vertex_descriptor* pointer;
      typedef vertex_descriptor& reference;
      typedef typename out_edge_container::iterator edge_iter_impl;
      typedef typename std::iterator_traits<out_edge_iterator>::iterator_category iterator_category;
      
      edge_iter_impl base;
      
      explicit child_vertex_iterator(edge_iter_impl aEdgeIt) : base(aEdgeIt) { };
      
      friend bool operator==( const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { return (lhs.base == rhs.base); };
      friend bool operator!=( const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { return (lhs.base != rhs.base); };
      friend bool operator >( const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { 
	return (lhs.base > rhs.base);
      };
      friend bool operator >=(const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { 
	return (lhs.base >= rhs.base);
      };
      friend bool operator <( const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { 
	return (lhs.base < rhs.base);
      };
      friend bool operator <=(const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) { 
	return (lhs.base <= rhs.base);
      };
      
      child_vertex_iterator& operator++() { ++base; return *this; };
      child_vertex_iterator operator++(int) { child_vertex_iterator result(*this); ++base; return result; };
      child_vertex_iterator& operator--() { --base; return *this; };
      child_vertex_iterator operator--(int) { child_vertex_iterator result(*this); --base; return result; };
      
      friend child_vertex_iterator operator+(const child_vertex_iterator& lhs, difference_type i) {
	return child_vertex_iterator(std::advance(lhs.base,i));
      };
      friend child_vertex_iterator operator+(difference_type i, const child_vertex_iterator& rhs) {
	return child_vertex_iterator(std::advance(rhs.base,i));
      };
      friend child_vertex_iterator operator-(const child_vertex_iterator& lhs, difference_type i) {
	return child_vertex_iterator(std::advance(lhs.base,-i));
      };
      friend difference_type operator-(const child_vertex_iterator& lhs, const child_vertex_iterator& rhs) {
        return std::distance(rhs.base, lhs.base);
      };
      
      child_vertex_iterator& operator +=(difference_type i) { base = std::advance(base,i); return *this; };
      child_vertex_iterator& operator -=(difference_type i) { base = std::advance(base,-i); return *this; };
      
      value_type operator[](difference_type i) const { return (*std::advance(base,i))->target; };
      reference operator*() { return (*base)->target; };
    };
    
  private:
    
    container_type m_vertices;
    vertex_descriptor m_root;
    
    typedef typename boost::mpl::if_< vertex_rand_access,
      std::vector< vertex_descriptor >,
      int >::type avail_node_queue;
    
    avail_node_queue m_avail_nodes;
    
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    vertex_value_type& >::type get_vertex_value(vertex_descriptor v) {
      return m_vertices[v];
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    vertex_value_type& >::type get_vertex_value(vertex_descriptor v) {
      return *static_cast<vertex_value_type*>(v);
    };
    
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    const vertex_value_type& >::type get_vertex_value(vertex_descriptor v) const {
      return m_vertices[v];
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    const vertex_value_type& >::type get_vertex_value(vertex_descriptor v) const {
      return *static_cast<const vertex_value_type*>(v);
    };
    
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    vertex_descriptor >::type add_new_vertex(const vertex_property_type& vp) {
      if(m_avail_nodes.empty()) {
        m_vertices.push_back(vertex_value_type(vp));
        return m_vertices.size() - 1;
      } else {
	vertex_descriptor result = m_avail_nodes.front();
	m_vertices[result] = vertex_value_type(vp);
	std::pop_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
	m_avail_nodes.pop_back();
	return result;
      };
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    vertex_descriptor >::type add_new_vertex(const vertex_property_type& vp) {
      m_vertices.push_back(vertex_value_type(vp));
      return &m_vertices.back();
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    vertex_descriptor >::type add_new_vertex(vertex_property_type&& vp) {
      if(m_avail_nodes.empty()) {
        m_vertices.push_back(vertex_value_type(std::move(vp)));
        return m_vertices.size() - 1;
      } else {
	vertex_descriptor result = m_avail_nodes.front();
	m_vertices[result] = vertex_value_type(std::move(vp));
	std::pop_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
	m_avail_nodes.pop_back();
	return result;
      };
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    vertex_descriptor >::type add_new_vertex(vertex_property_type&& vp) {
      m_vertices.push_back(vertex_value_type(std::move(vp)));
      return &m_vertices.back();
    };
#endif
    
    
    
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    void >::type add_root_vertex(const vertex_property_type& vp) {
      if(m_avail_nodes.empty()) {
        m_vertices.push_back(vertex_value_type(vp));
	m_root = m_vertices.size() - 1;
      } else {
	m_root = m_avail_nodes.front();
	m_vertices[result] = vertex_value_type(vp);
	std::pop_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
	m_avail_nodes.pop_back();
      };
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    void >::type add_root_vertex(const vertex_property_type& vp) {
      m_vertices.push_back(vertex_value_type(vp));
      m_root = &m_vertices.back();
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    void >::type add_root_vertex(vertex_property_type&& vp) {
      if(m_avail_nodes.empty()) {
        m_vertices.push_back(vertex_value_type(std::move(vp)));
        m_root = m_vertices.size() - 1;
      } else {
	m_root = m_avail_nodes.front();
	m_vertices[result] = vertex_value_type(std::move(vp));
	std::pop_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
	m_avail_nodes.pop_back();
      };
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    void >::type add_root_vertex(vertex_property_type&& vp) {
      m_vertices.push_back(vertex_value_type(std::move(vp)));
      m_root = &m_vertices.back();
    };
#endif
    
    
    // Random-access version of the remove function.
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    void >::type remove_branch_recursion(vertex_descriptor v) {
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(typename out_edge_container::iterator ei = m_vertices[v].out_edges.begin(); 
	  ei != m_vertices[v].out_edges.end(); ++ei)
	remove_branch_recursion<IsRandAccess>(ei->target);
      m_vertices[v].in_edge = null_edge();
      m_vertices[v].out_edges.clear();
      m_avail_nodes.push_back(v);
      std::push_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
    };
    
    template <typename IsRandAccess>
    typename boost::enable_if< IsRandAccess,
    void >::type remove_branch_impl(vertex_descriptor v) {
      if(v == m_root) {
	//remove_branch_recursion<IsRandAccess>(v);
        m_vertices.clear();
        m_root = null_vertex();
        m_avail_nodes.clear();
	return;
      };
      vertex_descriptor parent = m_vertices[v].in_edge.source;
      remove_branch_recursion<IsRandAccess>(v);
      // find the edge to be removed:
      typename out_edge_container::iterator ei = m_vertices[parent].out_edges.begin();
      while( (ei != m_vertices[parent].out_edges.end()) && (ei->target != v) )
	++ei;
      m_vertices[parent].out_edges.erase(ei); // remove the edge.
      // update the edge-descriptors in the children nodes.
      for(ei = m_vertices[parent].out_edges.begin(); ei != m_vertices[parent].out_edges.end(); ++ei)
	m_vertices[ei->target].in_edge.edge_id = &(*ei);
    };
    
    // Non-random-access version of the remove function.
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    void >::type remove_branch_recursion(vertex_descriptor v) {
      vertex_value_type* v_ptr = static_cast<vertex_value_type*>(v);
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(typename out_edge_container::iterator ei = v_ptr->out_edges.begin(); 
	  ei != v_ptr->out_edges.end(); ++ei)
	remove_branch_recursion<IsRandAccess>(ei->target);
      v_ptr->in_edge = null_edge();
      v_ptr->out_edges.clear();
    };
    
    template <typename IsRandAccess>
    typename boost::disable_if< IsRandAccess,
    void >::type remove_branch_impl(vertex_descriptor v) {
      if(v == m_root) {
	//remove_branch_recursion<IsRandAccess>(v);
        m_vertices.clear();
        m_root = null_vertex();
        m_avail_nodes.clear();
	return;
      };
      vertex_value_type* parent = static_cast<vertex_value_type*>((static_cast<vertex_value_type*>(v))->in_edge.source);
      remove_branch_recursion<IsRandAccess>(v);
      // find the edge to be removed:
      typename out_edge_container::iterator ei = parent->out_edges.begin();
      while( (ei != parent->out_edges.end()) && (ei->target != v) )
	++ei;
      parent->out_edges.erase(ei); // remove the edge.
      // update the edge-descriptors in the children nodes.
      for(ei = parent->out_edges.begin(); ei != parent->out_edges.end(); ++ei)
	(static_cast<vertex_value_type*>(ei->target))->in_edge.edge_id = &(*ei);
    };
    
    
    
    
    // Random-access version of the remove function.
    template <typename IsRandAccess, typename OutputIter>
    typename boost::enable_if< IsRandAccess,
    void >::type remove_branch_recursion(vertex_descriptor v, OutputIter it_out) {
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
#ifdef RK_ENABLE_CXX0X_FEATURES
      *(it_out++) = std::move(m_vertices[v].data);
#else
      *(it_out++) = m_vertices[v].data;
#endif
      for(typename out_edge_container::iterator ei = m_vertices[v].out_edges.begin(); 
	  ei != m_vertices[v].out_edges.end(); ++ei)
	it_out = remove_branch_recursion<IsRandAccess>(ei->target,it_out);
      m_vertices[v].in_edge = null_edge();
      m_vertices[v].out_edges.clear();
      m_avail_nodes.push_back(v);
      std::push_heap(m_avail_nodes.begin(),m_avail_nodes.end(),std::greater<vertex_descriptor>());
    };
    
    template <typename IsRandAccess, typename OutputIter>
    typename boost::enable_if< IsRandAccess,
    void >::type remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      if(v == m_root) {
	it_out = remove_branch_recursion<IsRandAccess>(v,it_out);
        m_vertices.clear();
        m_root = null_vertex();
        m_avail_nodes.clear();
	return;
      };
      vertex_descriptor parent = m_vertices[v].in_edge.source;
      it_out = remove_branch_recursion<IsRandAccess>(v,it_out);
      // find the edge to be removed:
      typename out_edge_container::iterator ei = m_vertices[parent].out_edges.begin();
      while( (ei != m_vertices[parent].out_edges.end()) && (ei->target != v) )
	++ei;
      m_vertices[parent].out_edges.erase(ei); // remove the edge.
      // update the edge-descriptors in the children nodes.
      for(ei = m_vertices[parent].out_edges.begin(); ei != m_vertices[parent].out_edges.end(); ++ei)
	m_vertices[ei->target].in_edge.edge_id = &(*ei);
    };
    
    // Non-random-access version of the remove function.
    template <typename IsRandAccess, typename OutputIter>
    typename boost::disable_if< IsRandAccess,
    void >::type remove_branch_recursion(vertex_descriptor v, OutputIter it_out) {
      vertex_value_type* v_ptr = static_cast<vertex_value_type*>(v);
#ifdef RK_ENABLE_CXX0X_FEATURES
      *(it_out++) = std::move(v_ptr->data);
#else
      *(it_out++) = v_ptr->data;
#endif
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(typename out_edge_container::iterator ei = v_ptr->out_edges.begin(); 
	  ei != v_ptr->out_edges.end(); ++ei)
	it_out = remove_branch_recursion<IsRandAccess>(ei->target,it_out);
      v_ptr->in_edge = null_edge();
      v_ptr->out_edges.clear();
    };
    
    template <typename IsRandAccess, typename OutputIter>
    typename boost::disable_if< IsRandAccess,
    OutputIter >::type remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      if(v == m_root) {
	it_out = remove_branch_recursion<IsRandAccess>(v,it_out); 
        m_vertices.clear();
        m_root = null_vertex();
        m_avail_nodes.clear();
	return it_out;
      };
      vertex_value_type* parent = static_cast<vertex_value_type*>((static_cast<vertex_value_type*>(v))->in_edge.source);
      it_out = remove_branch_recursion<IsRandAccess>(v,it_out);
      // find the edge to be removed:
      typename out_edge_container::iterator ei = parent->out_edges.begin();
      while( (ei != parent->out_edges.end()) && (ei->target != v) )
	++ei;
      parent->out_edges.erase(ei); // remove the edge.
      // update the edge-descriptors in the children nodes.
      for(ei = parent->out_edges.begin(); ei != parent->out_edges.end(); ++ei)
	(static_cast<vertex_value_type*>(ei->target))->in_edge.edge_id = &(*ei);
    };
    
    
    
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    std::size_t get_depth(const vertex_value_type& aNode) const {
      std::size_t max_depth = 0;
      typedef typename out_edge_container::const_iterator EdgeIter;
      for(EdgeIter ei = aNode.out_edges.begin(); ei != aNode.out_edges.end(); ++ei) {
	std::size_t temp = get_depth(get_vertex_value<vertex_rand_access>(ei->target));
	if(temp > max_depth)
	  max_depth = temp;
      };
      return max_depth + 1;
    };
    
  public:
    
    /**
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() {
      if(vertex_rand_access::type::value)
        return reinterpret_cast<vertex_descriptor>(-1);
      else
	return reinterpret_cast<vertex_descriptor>(NULL);
    };
    
    /**
     * This static member function outputs the null-edge (invalid edge descriptor).
     * \return A null-edge descriptor (invalid edge descriptor).
     */
    static edge_descriptor null_edge() { 
      return edge_descriptor(reinterpret_cast<vertex_descriptor>(-1), NULL);
    };
    
    /**
     * Constructs an empty linked-tree.
     */
    linked_tree() : m_vertices(), m_root(null_vertex()) { };
    
    /**
     * Checks if the tree is empty.
     * \return True if the tree is empty.
     */
    bool empty() const { return m_vertices.empty(); };
    
    /**
     * Returns the size of the tree (the number of vertices it contains).
     * \return The size of the tree (the number of vertices it contains).
     */
    std::size_t size() const { return m_vertices.size(); };
    
    /**
     * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
     * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
     */
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    /**
     * Returns the depth of the tree.
     * \note This operation must recurse through all the branches of the tree (depth-first), and is 
     * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack) 
     * w.r.t. the depth of tree).
     * \return The depth of the tree.
     */
    std::size_t depth() const { 
      if(m_root != null_vertex())
        return get_depth(get_vertex_value<vertex_rand_access>(m_root))
      return 0;
    };
    
    /**
     * Standard swap function.
     */
    friend
    void swap(self& lhs, self& rhs) {
      using std::swap;
      lhs.m_vertices.swap(rhs.m_vertices);
      swap(lhs.m_root, rhs.m_root);
      lhs.m_avail_nodes.swap(rhs.m_avail_nodes);
    };
    
    /**
     * Clears the tree of all vertices and edges.
     */
    void clear() { 
      m_vertices.clear();
      m_root = null_vertex();
      m_avail_nodes.clear();
    };
    
    /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[]( const vertex_descriptor& v_i) {
      return get_vertex_value<vertex_rand_access>(v_i).data;
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[]( const vertex_descriptor& v_i) const {
      return get_vertex_value<vertex_rand_access>(v_i).data;
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[]( const edge_descriptor& e_i) {
      return static_cast<edge_value_type*>(e_i.edge_id)->data;
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[]( const edge_descriptor& e_i) const {
      return static_cast<edge_value_type*>(e_i.edge_id)->data;
    };
    
    /**
     * Indexing function. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \param g The tree from which to draw the vertex.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    friend
    vertex_property_type& get_property(const vertex_descriptor& v_i, self& g) {
      return g[v_i];
    };
    /**
     * Indexing function. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \param g The tree from which to draw the vertex.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    friend
    const vertex_property_type& get_property( const vertex_descriptor& v_i, const self& g) {
      return g[v_i];
    };
    /**
     * Indexing function. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \param g The tree from which to draw the edge.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    friend
    edge_property_type& get_property(const edge_descriptor& e_i, self& g) {
      return g[e_i];
    };
    /**
     * Indexing function. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \param g The tree from which to draw the edge.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    friend
    const edge_property_type& get_property( const edge_descriptor& e_i, const self& g) {
      return g[e_i];
    };
    
    /**
     * Get function. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param g The tree from which to draw the vertex.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    friend const vertex_property_type& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };
    
    /**
     * Put function. Sets the vertex-property associated to the given vertex descriptor to the given value.
     * \param g The tree from which to draw the vertex.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \param value The vertex-property value to assign to the given vertex.
     */
    friend void put( self& g, const vertex_descriptor& v_i, const vertex_property_type& value) {
      g[v_i] = value;
    };
    
    /**
     * Get function. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param g The tree from which to draw the edge.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    friend const edge_property_type& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };
    
    /**
     * Put function. Sets the edge-property associated to the given edge descriptor to the given value.
     * \param g The tree from which to draw the edge.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \param value The edge-property value to assign to the given edge.
     */
    friend void put( self& g, const edge_descriptor& e_i, const edge_property_type& value) {
      g[e_i] = value;
    };
    
    template <typename T, typename Bundle>
    friend 
    typename boost::property_map< self, T Bundle::* >::type
    get( T Bundle::* p, self& g) {
      return typename boost::property_map< self, T Bundle::* >::type(&g, p);
    };
    
    template <typename T, typename Bundle>
    friend 
    typename boost::property_map< self, T Bundle::* >::const_type
    get( T Bundle::* p, const self& g) {
      return typename boost::property_map< self, T Bundle::* >::const_type(&g, p);
    };
    
    /**
     * Checks if a given vertex descriptor leads to a valid vertex of the tree.
     * \param v_i The vertex descriptor to test for validity.
     * \return True if the given vertex is valid.
     */
    bool is_valid(const vertex_descriptor& v_i) const {
      return (v_i != null_vertex()) && !((v_i != m_root) && (get_vertex_value<vertex_rand_access>(v_i).in_edge == null_edge())) ;
    };
    
    /**
     * Returns the out-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The out-degree of the given vertex descriptor.
     */
    edges_size_type get_raw_out_degree( const vertex_descriptor& v_i) const {
      if(is_valid(v_i))
        return get_vertex_value<vertex_rand_access>(v_i).out_edges.size();
      return 0;
    };
    
    /**
     * Returns the out-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The out-degree of the given vertex descriptor.
     */
    edges_size_type get_out_degree( const vertex_descriptor& v_i) const {
      return get_raw_out_degree(v_i);
    };
    
    /**
     * Returns the in-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
     */
    edges_size_type get_in_degree( const vertex_descriptor& v_i) const {
      if(is_valid(v_i) && (v_i != m_root))
        return 1;
      return 0;
    };
    
    /**
     * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The edge iterator range for the out-edges of a given vertex descriptor.
     */
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor v) {
      return std::make_pair(out_edge_iterator(v,get_vertex_value<vertex_rand_access>(v).out_edges.begin()),
			    out_edge_iterator(v,get_vertex_value<vertex_rand_access>(v).out_edges.end()));
    };
    
    /**
     * Returns the edge iterator range for the in-edges of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The edge iterator range for the in-edges of a given vertex descriptor.
     */
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) {
      return std::make_pair(&(get_vertex_value<vertex_rand_access>(v).in_edge),
			    &(get_vertex_value<vertex_rand_access>(v).in_edge) + 1);
    };
    
    /**
     * Returns the vertex iterator range for all the vertices of the tree.
     * \return The vertex iterator range for all the vertices of the tree.
     */
    std::pair< vertex_iterator, vertex_iterator > vertices() {
      return std::make_pair(vertex_iterator::begin(m_vertices),
			    vertex_iterator::end(m_vertices));
    };
    
    /**
     * Returns the edge descriptor for the edge between two given vertex descriptors.
     * \param u The vertex descriptor of the source vertex.
     * \param v The vertex descriptor of the target vertex.
     * \return The edge descriptor for the given vertex descriptor pair.
     */
    edge_descriptor get_edge( vertex_descriptor u, vertex_descriptor v) const {
      typename out_edge_container::const_iterator ei = get_vertex_value<vertex_rand_access>(u).out_edges.begin();
      while( (ei != get_vertex_value<vertex_rand_access>(u).out_edges.end()) && (ei->target != v) )
	++ei;
      return edge_descriptor(u, &(*ei));
    };
    
    /**
     * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
     * \param v The vertex descriptor whose children are sought.
     * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
     */
    std::pair< child_vertex_iterator, child_vertex_iterator > child_vertices(vertex_descriptor v) {
      return std::make_pair(child_vertex_iterator(get_vertex_value<vertex_rand_access>(v).out_edges.begin()),
			    child_vertex_iterator(get_vertex_value<vertex_rand_access>(v).out_edges.end()));
    };
    
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values.
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value for the newly created vertex.
     * \param ep The property value for the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     const vertex_property_type& vp = vertex_property_type(), 
							     const edge_property_type& ep = edge_property_type()) {
      if(is_valid(v))
	throw std::range_error("Cannot add child-node to an empty node!");
      // create a new node.
      vertex_descriptor new_node = add_new_vertex<vertex_rand_access>(vp);
      // create a new edge.
      (*this)[v].out_edges.push_back(edge_value_type(new_node,ep));
      (*this)[new_node].in_edge = edge_descriptor(v,&((*this)[v].out_edges.back()));
      return std::make_pair(new_node, (*this)[new_node].in_edge);
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values, by move-semantics (C++11).
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value to be moved into the newly created vertex.
     * \param ep The property value to be moved into the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     vertex_property_type&& vp, 
							     edge_property_type&& ep = edge_property_type()) {
      if(is_valid(v))
	throw std::range_error("Cannot add child-node to an empty node!");
      // create a new node.
      vertex_descriptor new_node = add_new_vertex<vertex_rand_access>(std::move(vp));
      // create a new edge.
      (*this)[v].out_edges.push_back(edge_value_type(new_node,std::move(ep)));
      (*this)[new_node].in_edge = edge_descriptor(v,&((*this)[v].out_edges.back()));
      return std::make_pair(new_node, (*this)[new_node].in_edge);
    };
#endif
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
    void remove_branch(vertex_descriptor v) {
      if(is_valid(v))
	return;  // vertex is already deleted.
      remove_branch_impl<vertex_rand_access>(v);
    };
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex, while 
     * recording the vertex-properties of all the removed vertices into an output-iterator.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \return The output-iterator after the collection of all the removed vertices.
     * \note The first vertex-property to figure in the output range is that of the vertex v.
     */
    template <typename OutputIter>
    OutputIter remove_branch(vertex_descriptor v, OutputIter it_out) {
      if(is_valid(v))
	return it_out;  // vertex is already deleted.
      return remove_branch_impl<vertex_rand_access>(v,i_out);
    };
    
    /**
     * Returns the vertex-descriptor of the root of the tree.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor get_root_vertex() const {
      return m_root; 
    };
    
    /**
     * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
     * \param vp The vertex-property to assign to the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(m_vertex_count)
	clear();
      add_root_vertex<vertex_rand_access>(vp);
      return m_root;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
     * \param vp The vertex-property to move into the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(vertex_property_type&& vp) {
      if(m_vertex_count)
	clear();
      add_root_vertex<vertex_rand_access>(std::move(vp));
      return m_root;
    };
#endif
    
    
    
};


/**
 * This is the tree-storage specifier for a linked-tree of the given edge and vertex storage policies.
 */
template <typename OutEdgeListS, typename VertexListS>
struct linked_tree_storage { };


template <typename VertexDescriptor, typename EdgeDescriptor, typename OutEdgeListS, typename VertexListS>
struct tree_storage<VertexDescriptor, EdgeDescriptor, linked_tree_storage<OutEdgeListS, VertexListS> > {
  typedef linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor> type;
};


template <typename OutEdgeListS, typename VertexListS>
struct tree_storage_traits< linked_tree_storage<OutEdgeListS, VertexListS> > :
  linked_tree_traits<OutEdgeListS, VertexListS> { };





/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  source( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor& e,
	  const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>&) {
  return e.source;
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  target( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor& e,
	  const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>&) {
  typedef typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_value_type ValueType;
  return (static_cast<ValueType*>(e.edge_id))->target;
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair<
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::out_edge_iterator,
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::out_edge_iterator >
  out_edges( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
	     linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.out_edges(v);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::size_t
  out_degree( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
	      const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair<
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::in_edge_iterator,
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::in_edge_iterator >
  in_edges( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
	    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.in_edges(v);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::size_t
  in_degree( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
             const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.get_in_degree(v);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::size_t
  degree( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
          const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.get_in_degree(v) + g.get_out_degree(v);
};



/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair<
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_iterator,
 typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_iterator >
  vertices( linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.vertices();
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertices_size_type
  num_vertices( const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.size();
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor
  edge( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& u,
	const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
        const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.get_edge(u,v);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  get_root_vertex( const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.get_root_vertex();
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_iterator,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_iterator >
  child_vertices( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
                  linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.child_vertices(v);
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  create_root( linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.create_root_vertex();
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor >
  add_child_vertex( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.add_child(v);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
void remove_branch( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.remove_branch(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  create_root( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type& vp, 
	       linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.create_root_vertex(vp);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor
  create_root( typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type&& vp, 
	       linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.create_root_vertex(std::move(vp));
};
#endif

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor >
  add_child_vertex( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
		    const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type& vp,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.add_child(v,vp);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor >
  add_child_vertex( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
		    const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type& vp,
		    const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_property_type& ep,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.add_child(v,vp,ep);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor >
  add_child_vertex( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
		    typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type&& vp,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.add_child(v,std::move(vp));
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
std::pair< 
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor,
typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor >
  add_child_vertex( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
		    typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_property_type&& vp,
		    typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_property_type&& ep,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.add_child(v,std::move(vp),std::move(ep));
};
#endif

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS,
	  typename OutputIter>
inline
OutputIter remove_branch( const typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor& v,
		    OutputIter it_out,
                    linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.remove_branch(v,it_out);
};



/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/


template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
bool is_vertex_valid( typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::vertex_descriptor u,
                      const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.is_valid(u);
};

template <typename VertexDescriptor, typename EdgeDescriptor, 
          typename OutEdgeListS, typename VertexListS>
inline
bool is_edge_valid( typename linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>::edge_descriptor e,
                    const linked_tree<OutEdgeListS, VertexListS, VertexDescriptor, EdgeDescriptor>& g) {
  return g.is_valid(e.source) && g.is_valid(target(e,g));
};




};

};


#endif


















