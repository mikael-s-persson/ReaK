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

#ifndef REAK_D_ARY_BF_TREE_HPP
#define REAK_D_ARY_BF_TREE_HPP

#include <boost/graph/properties.hpp>

#include <vector>
#include <stdexcept>
#include <map>
#include <iterator>

#include "bgl_tree_adaptor.hpp"

#include <boost/graph/adjacency_list.hpp>

namespace ReaK {

namespace graph {
  


template <typename OutEdgeListS = boost::vecS, 
          typename VertexListS = boost::vecS>
struct linked_tree_traits {
  typedef typename boost::detail::is_random_access<VertexListS>::type is_rand_access;
  typedef boost::mpl::true_ is_bidir;
  typedef boost::mpl::true_ is_directed;
  
  typedef boost::bidirectional_tag directed_category;
  
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  
  typedef std::size_t vertices_size_type;
  typedef void* vertex_ptr;
  typedef typename boost::mpl::if_< is_rand_access,
    vertices_size_type, vertex_ptr>::type vertex_descriptor;
  
  struct edge_descriptor {
    vertex_descriptor source;
    void* edge_id;
    
    edge_descriptor(vertex_descriptor aSrc, void* aEdgeId) : source(aSrc), edge_id(aEdgeId) { };
  };
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
    };
    typedef typename boost::container_gen<OutEdgeListS, edge_value_type>::type out_edge_container;
    
    struct vertex_value_type {
      vertex_property_type data;
      out_edge_container out_edges;
      edge_descriptor in_edge;
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
    vertices_size_type m_vertex_count;
    vertex_descriptor m_root;
    
  public:
    
    
    static vertex_descriptor null_vertex() { 
      return reinterpret_cast<vertex_descriptor>(-1);
    };
    
    static edge_descriptor null_edge() { 
      return edge_descriptor(reinterpret_cast<vertex_descriptor>(-1), NULL);
    };
    
    /**
     * Construct the D-ary BF-tree with a given reserved depth.
     * \param aDepth The depth of the graph to reserve space for.
     */
    linked_tree() : m_vertex_count(0), m_root(null_vertex()) { };
    
    /**
     * Checks if the DVP-tree is empty.
     * \return True if the DVP-tree is empty.
     */
    bool empty() const { return m_vertex_count == 0; };
    
    /**
     * Returns the size of the DVP-tree (the number of vertices it contains.
     * \return The size of the DVP-tree (the number of vertices it contains.
     */
    std::size_t size() const { return m_vertex_count; };
    
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    std::size_t depth() const { // TODO
      vertices_size_type vert_count = 1;
      vertices_size_type accum = 1;
      vertices_size_type depth_count = 0;
      for(; vert_count < m_vertices.size(); ++depth_count) {
	accum *= Arity;
	vert_count += accum;
      };
      return depth_count; 
    };
    
    void swap(self& rhs) {
      using std::swap;
      m_vertices.swap(rhs.m_vertices);
      swap(m_vertex_count, rhs.m_vertex_count);
      swap(m_root, rhs.m_root);
    };
    
    void clear() { 
      m_vertices.clear();
      m_vertices.resize(1);
      m_vertex_count = 0;
      m_root = null_vertex();
    };
    
    vertex_property_type& operator[]( const vertex_descriptor& v_i) {
      return m_vertices[v_i].v;
    };
    const vertex_property_type& operator[]( const vertex_descriptor& v_i) const {
      return m_vertices[v_i].v;
    };
    edge_property_type& operator[]( const edge_descriptor& e_i) {
      return m_vertices[e_i.source_vertex].e[e_i.edge_index];
    };
    const edge_property_type& operator[]( const edge_descriptor& e_i) const {
      return m_vertices[e_i.source_vertex].e[e_i.edge_index];
    };
    
    friend
    vertex_property_type& get_property(const vertex_descriptor& v_i, self& g) {
      return g[v_i];
    };
    friend
    const vertex_property_type& get_property( const vertex_descriptor& v_i, const self& g) {
      return g[v_i];
    };
    friend
    edge_property_type& get_property(const edge_descriptor& e_i, self& g) {
      return g[e_i];
    };
    friend
    const edge_property_type& get_property( const edge_descriptor& e_i, const self& g) {
      return g[e_i];
    };
    
    friend const vertex_property_type& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };
    
    friend void put( self& g, const vertex_descriptor& v_i, const vertex_property_type& value) {
      g[v_i] = value;
    };
    
    friend const edge_property_type& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };
    
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
    
    
    bool is_valid(const vertex_descriptor& v_i) const {
      return (v_i < m_vertices.size()) && ( m_vertices[v_i].out_degree >= 0 );
    };
    
    edges_size_type get_raw_out_degree( const vertex_descriptor& v_i) const {
      if( m_vertices[v_i].out_degree < 0 )
	return 0;
      else 
	return m_vertices[v_i].out_degree;
    };
    
    edges_size_type get_out_degree( const vertex_descriptor& v_i) const {
      if( m_vertices[v_i].out_degree < 0 )
	return 0;
      else {
	edges_size_type result = 0;
	for( edges_size_type i = 0; i < edges_size_type(m_vertices[v_i].out_degree); ++i) {
	  if( ( v_i * Arity + i + 1 < m_vertices.size() ) &&
	      ( m_vertices[v_i * Arity + i + 1].out_degree >= 0 ) )
	    ++result;
	};
	return result;
      };
    };
    
    edges_size_type get_in_degree( const vertex_descriptor& v_i) const {
      if(( v_i == 0 ) || ( m_vertices[v_i].out_degree < 0 ))
	return 0;
      else
	return 1;
    };
    
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     const vertex_property_type& vp = vertex_property_type(), 
							     const edge_property_type& ep = edge_property_type()) {
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) ) 
	throw std::range_error("Cannot add child-node to an empty node!");
      int new_edge = 0;
      for(; new_edge < m_vertices[v].out_degree; ++new_edge)
	if( m_vertices[Arity * v + 1 + new_edge].out_degree < 0 )
	  break;
      if( new_edge == Arity ) 
	throw std::range_error("Cannot add child-node to a full node!");
      vertex_descriptor result = Arity * v + 1 + new_edge;
      if( result >= vertex_descriptor(m_vertices.size()) )
	m_vertices.resize(result + 1);
      m_vertices[result].out_degree = 0;
      m_vertices[result].v = vp;
      m_vertices[v].e[new_edge] = ep;
      if( new_edge == m_vertices[v].out_degree )
        ++(m_vertices[v].out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     vertex_property_type&& vp, 
							     edge_property_type&& ep = edge_property_type()) {
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) ) 
	throw std::range_error("Cannot add child-node to an empty node!");
      int new_edge = 0;
      for(; new_edge < m_vertices[v].out_degree; ++new_edge)
	if( m_vertices[Arity * v + 1 + new_edge].out_degree < 0 )
	  break;
      if( new_edge == Arity ) 
	throw std::range_error("Cannot add child-node to a full node!");
      vertex_descriptor result = Arity * v + 1 + new_edge;
      if( result >= vertex_descriptor(m_vertices.size()) )
	m_vertices.resize(result + 1);
      m_vertices[result].out_degree = 0;
      m_vertices[result].v = std::move(vp);
      m_vertices[v].e[new_edge] = std::move(ep);
      if( new_edge == m_vertices[v].out_degree )
        ++(m_vertices[v].out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
    };
#endif
    
    void update_out_degree(vertex_descriptor v) {
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) )
	return;  // vertex is already updated.
      for( int i = Arity; i > 0; --i) {
	vertex_descriptor u = Arity * v + i;
	if( ( u < vertex_descriptor(m_vertices.size()) ) &&
	    ( m_vertices[u].out_degree >= 0 ) ) {
	  m_vertices[v].out_degree = i;
	  return;
	};
      };
      m_vertices[v].out_degree = 0;
    };
    
    void remove_branch(vertex_descriptor v) {
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) )
	return;  // vertex is already deleted.
      --m_vertex_count;
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      int max_child = m_vertices[v].out_degree;
      for( int i = 0; i < max_child; ++i)
	remove_branch(Arity * v + 1 + i);
      m_vertices[v].out_degree = -1;
      if( v != 0 )  // if the node is not the root one, then update the out-degree of the parent node:
	update_out_degree( (v - 1) / Arity );
      // remove empty vertices from the end of the container:
      if( v == vertex_descriptor(m_vertices.size() - 1) ) {
	while( (v > 0) && ( m_vertices[v].out_degree < 0 ) )
	  --v;
	++v;
	m_vertices.erase(m_vertices.begin() + v, m_vertices.end());
      };
    };
    
    template <typename OutputIter>
    void remove_branch(vertex_descriptor v, OutputIter it_out) {
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) )
	return;  // vertex is already deleted.
      --m_vertex_count;
#ifdef RK_ENABLE_CXX0X_FEATURES
      *(it_out++) = std::move(m_vertices[v].v);
#else
      *(it_out++) = m_vertices[v].v;
#endif
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      int max_child = m_vertices[v].out_degree;
      for( int i = 0; i < max_child; ++i)
	remove_branch(Arity * v + 1 + i, it_out);
      m_vertices[v].out_degree = -1;
      if( v != 0 )  // if the node is not the root one, then update the out-degree of the parent node:
	update_out_degree( (v - 1) / Arity );
      // remove empty vertices from the end of the container:
      if( v == vertex_descriptor(m_vertices.size() - 1) ) {
	while( (v > 0) && ( m_vertices[v].out_degree < 0 ) )
	  --v;
	++v;
	m_vertices.erase(m_vertices.begin() + v, m_vertices.end());
      };
    };
    
    vertex_descriptor get_root_vertex() const {
      return 0; 
    };
    
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(m_vertices[0].out_degree >= 0)
	remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].v = vp;
      ++m_vertex_count;
      return 0;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    vertex_descriptor create_root_vertex(vertex_property_type&& vp) {
      if(m_vertices[0].out_degree >= 0)
	remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].v = std::move(vp);
      ++m_vertex_count;
      return 0;
    };
#endif
    
    
    
};


template <std::size_t Arity = 2>
struct d_ary_bf_tree_storage { };


template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity>
struct tree_storage<VertexDescriptor, EdgeDescriptor, d_ary_bf_tree_storage<Arity> > {
  typedef d_ary_bf_tree<VertexDescriptor, Arity, EdgeDescriptor> type;
};


template <std::size_t Arity>
struct tree_storage_traits< d_ary_bf_tree_storage<Arity> > {
  typedef boost::mpl::true_ is_rand_access;
  typedef boost::mpl::true_ is_bidir;
  typedef boost::mpl::true_ is_directed;
  
  typedef typename boost::mpl::if_< is_bidir,
    boost::bidirectional_tag,
    typename boost::mpl::if_< is_directed,
      directed_tag, undirected_tag
    >::type
  >::type directed_category;
  
  typedef disallow_parallel_edge_tag edge_parallel_category;
  
  typedef std::size_t vertices_size_type;
  typedef void* vertex_ptr;
  typedef typename d_ary_bf_tree<int,Arity>::vertex_descriptor vertex_descriptor;  // the value-type doesn't affect the vertex_descriptor type (int is a dummy type here).
  typedef typename d_ary_bf_tree<int,Arity>::edge_descriptor edge_descriptor;  // the value-type doesn't affect the edge_descriptor type (int is a dummy type here).
  typedef std::size_t edges_size_type;
  
};

};




/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  source( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor& e,
	  const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  return e.source_vertex;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  target( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor& e,
	  const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  return Arity * e.source_vertex + 1 + e.edge_index;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator >
  out_edges( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	     const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator OutIter;
  return std::make_pair(OutIter(v,0),OutIter(v,g.get_raw_out_degree(v)));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::size_t
  out_degree( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	      const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::in_edge_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::in_edge_iterator >
  in_edges( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	    const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::in_edge_iterator InIter;
  return std::make_pair(InIter(v),InIter());
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::size_t
  in_degree( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
             const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.get_in_degree(v);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::size_t
  degree( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
          const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.get_in_degree(v) + g.get_out_degree(v);
};


/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::adjacency_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::adjacency_iterator >
  adjacent_vertices( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	             const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::adjacency_iterator AdjIter;
  return std::make_pair(AdjIter(((v - 1) / Arity) * Arity + 1),
			AdjIter(((v - 1) / Arity + 1) * Arity));
};


/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator >
  vertices( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator VIter;
  return std::make_pair(VIter(0),
			VIter(g.capacity()));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertices_size_type
  num_vertices( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.size();
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_iterator >
  edges( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_iterator EIter;
  return std::make_pair(EIter(0,0),
			EIter(g.capacity(),Arity));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertices_size_type
  num_edges( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.size() - 1;
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor
  edge( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor&,
	const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
        const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor Edge;
  return Edge((v - 1) / Arity, (v - 1) % Arity);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  get_root_vertex( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.get_root_vertex();
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator >
  child_vertices( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
                  const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_iterator VIter;
  return std::make_pair(VIter(Arity * v + 1),VIter(Arity * (v + 1)));
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  create_root( d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.create_root_vertex();
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor >
  add_child_vertex( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.add_child(v);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
void remove_branch( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.remove_branch(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  create_root( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type& vp, 
	       d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.create_root_vertex(vp);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  create_root( typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type&& vp, 
	       d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.create_root_vertex(std::move(vp));
};
#endif

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor >
  add_child_vertex( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
		    const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type& vp,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.add_child(v,vp);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor >
  add_child_vertex( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
		    const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type& vp,
		    const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_property_type& ep,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.add_child(v,vp,ep);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor >
  add_child_vertex( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
		    typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type&& vp,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.add_child(v,std::move(vp));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor,
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor >
  add_child_vertex( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
		    typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_property_type&& vp,
		    typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_property_type&& ep,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.add_child(v,std::move(vp),std::move(ep));
};
#endif

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  typename OutputIter>
inline
void remove_branch( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
		    OutputIter it_out,
                    d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.remove_branch(v,it_out);
};



/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
bool is_vertex_valid( typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor u,
                      const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.is_valid(u);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
bool is_edge_valid( typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor e,
                    const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.is_valid(target(e,g));
};




};

};


#endif


















