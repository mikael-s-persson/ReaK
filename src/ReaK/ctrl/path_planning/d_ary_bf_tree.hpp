/**
 * \file d_ary_bf_tree.hpp
 * 
 * This library provides a class that implements a D-Ary Breadth-first tree that is tailored 
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
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

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include <map>
#include <unordered_map>
#include <vector>
#include "metric_space_concept.hpp"
#include "global_rng.hpp"
#include <lin_alg/vect_alg.hpp>


namespace ReaK {

namespace pp {


/**
 * This class implements a D-Ary Breadth-first tree that is tailored 
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
 * \tparam VertexProperties A POD type to be attached to each vertex in the tree.
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam EdgeProperties A POD type to be attached to each edge in the tree.
 */
template <typename VertexProperties,
          std::size_t Arity = 2,
          typename EdgeProperties = boost::no_property >
class d_ary_bf_tree
{
  public:
    typedef d_ary_bf_tree<VertexProperties, Arity, EdgeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    struct value_type {
      int out_degree;
      vertex_property_type v;
      edge_property_type e[Arity];
      
      value_type() : out_degree(-1) { };
    };
    
    typedef std::vector< value_type > container_type;
    
    typedef typename container_type::difference_type vertex_descriptor;
    
    typedef typename container_type::size_type vertices_size_type;
    typedef vertices_size_type edge_size_type;
    typedef vertices_size_type degree_size_type;
    
    
    struct edge_descriptor {
      vertex_descriptor source_vertex;
      std::ptrdiff_t edge_index;
      edge_descriptor(vertex_descriptor aSrc = 0, std::size_t aEdgeId = 0) : source_vertex(aSrc), edge_index(aEdgeId) { };
      
      friend bool operator==( const edge_descriptor& lhs, const edge_descriptor& rhs) { 
	return ((lhs.source_vertex == rhs.source_vertex) && (lhs.edge_index == rhs.edge_index)); 
      };
      friend bool operator!=( const edge_descriptor& lhs, const edge_descriptor& rhs) { 
	return ((lhs.source_vertex != rhs.source_vertex) || (lhs.edge_index != rhs.edge_index)); 
      };
      friend bool operator <( const edge_descriptor& lhs, const edge_descriptor& rhs) {
	if( lhs.source_vertex == rhs.source_vertex )
	  return ( lhs.edge_index < rhs.edge_index );
	else 
	  return ( lhs.source_vertex < rhs.source_vertex );
      };
      friend bool operator<=( const edge_descriptor& lhs, const edge_descriptor& rhs) {
	if( lhs.source_vertex == rhs.source_vertex )
	  return ( lhs.edge_index <= rhs.edge_index );
	else 
	  return ( lhs.source_vertex < rhs.source_vertex );
      };
      friend bool operator >( const edge_descriptor& lhs, const edge_descriptor& rhs) {
	if( lhs.source_vertex == rhs.source_vertex )
	  return ( lhs.edge_index > rhs.edge_index );
	else 
	  return ( lhs.source_vertex > rhs.source_vertex );
      };
      friend bool operator>=( const edge_descriptor& lhs, const edge_descriptor& rhs) {
	if( lhs.source_vertex == rhs.source_vertex )
	  return ( lhs.edge_index >= rhs.edge_index );
	else 
	  return ( lhs.source_vertex > rhs.source_vertex );
      };
    };
    
    struct out_edge_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef edge_descriptor value_type;
      typedef edge_descriptor* pointer;
      typedef edge_descriptor& reference;
      typedef std::random_access_iterator_tag iterator_category;
      
      edge_descriptor base;
      out_edge_iterator(const edge_descriptor& aBase = edge_descriptor()) : base(aBase) { };
      out_edge_iterator(vertex_descriptor aSrc, std::size_t aEdgeId) : base(aSrc, aEdgeId) { };
      
      friend bool operator==( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base == rhs.base; };
      friend bool operator!=( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base != rhs.base; };
      friend bool operator >( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base > rhs.base; };
      friend bool operator >=(const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base >= rhs.base; };
      friend bool operator <( const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base < rhs.base; };
      friend bool operator <=(const out_edge_iterator& lhs, const out_edge_iterator& rhs) { return lhs.base <= rhs.base; };
      
      out_edge_iterator& operator++() { ++base.edge_index; return *this; };
      out_edge_iterator operator++(int) { out_edge_iterator result(*this); ++base.edge_index; return result; };
      out_edge_iterator& operator--() { --base.edge_index; return *this; };
      out_edge_iterator operator--(int) { out_edge_iterator result(*this); --base.edge_index; return result; };
      
      friend out_edge_iterator operator+(const out_edge_iterator& lhs, difference_type i) {
	return out_edge_iterator(edge_descriptor(lhs.base.source_vertex, lhs.base.edge_index + i));
      };
      friend out_edge_iterator operator+(difference_type i, const out_edge_iterator& rhs) {
	return out_edge_iterator(edge_descriptor(rhs.base.source_vertex, rhs.base.edge_index + i));
      };
      friend out_edge_iterator operator-(const out_edge_iterator& lhs, difference_type i) {
	return out_edge_iterator(edge_descriptor(lhs.base.source_vertex, lhs.base.edge_index - i));
      };
      friend difference_type operator-(const out_edge_iterator& lhs, const out_edge_iterator& rhs) {
        if(lhs.base.source_vertex == rhs.base.source_vertex)
	  return difference_type(lhs.base.edge_index - rhs.base.edge_index);
	else
	  return Arity;
      };
      
      out_edge_iterator& operator +=(difference_type i) { base.edge_index += i; return *this; };
      out_edge_iterator& operator -=(difference_type i) { base.edge_index -= i; return *this; };
      
      value_type operator[](difference_type i) const { return edge_descriptor(base.source_vertex, base.edge_index + i); };
      reference operator*() { return base; };
      pointer operator->() { return &base; };
    };
    
    struct in_edge_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef edge_descriptor value_type;
      typedef edge_descriptor* pointer;
      typedef edge_descriptor& reference;
      typedef std::random_access_iterator_tag iterator_category;
      
      edge_descriptor base;
      in_edge_iterator(const vertex_descriptor& aBase = vertex_descriptor(-1)) : base((aBase - 1) / Arity, (aBase - 1) % Arity ) { 
	if((base.edge_index < 0) || (base.source_vertex < 0)) {
	  base.source_vertex = -1;
	  base.edge_index = -1;
	};
      };
      
      friend bool operator==( const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return lhs.base == rhs.base; };
      friend bool operator!=( const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return lhs.base != rhs.base; };
      friend bool operator >( const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return (lhs.base.edge_index < 0) && (rhs.base.edge_index >= 0); };
      friend bool operator >=(const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return (lhs.base.edge_index < 0); };
      friend bool operator <( const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return (rhs.base.edge_index < 0) && (lhs.base.edge_index >= 0); };
      friend bool operator <=(const in_edge_iterator& lhs, const in_edge_iterator& rhs) { return (rhs.base.edge_index < 0); };
      
      in_edge_iterator& operator++() { base.edge_index = -1; return *this; };
      in_edge_iterator operator++(int) { in_edge_iterator result(*this); base.edge_index = -1; return result; };
      in_edge_iterator& operator--() { base.edge_index = -1; return *this; };
      in_edge_iterator operator--(int) { in_edge_iterator result(*this); base.edge_index = -1; return result; };
      
      friend in_edge_iterator operator+(const in_edge_iterator& lhs, difference_type i) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return *this;
      };
      friend in_edge_iterator operator+(difference_type i, const in_edge_iterator& rhs) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return *this;
      };
      friend in_edge_iterator operator-(const in_edge_iterator& lhs, difference_type i) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return *this;
      };
      friend difference_type operator-(const in_edge_iterator& lhs, const in_edge_iterator& rhs) {
        if(lhs.base.edge_index == rhs.base.edge_index)
	  return 0;
	else if(rhs.base.edge_index == -1)
	  return 1;
	else
	  return -1;
      };
      
      in_edge_iterator& operator +=(difference_type i) { 
	if( i != 0 )
	  base.edge_index = -1; 
	return *this; 
      };
      in_edge_iterator& operator -=(difference_type i) { 
	if( i != 0 )
	  base.edge_index = -1;
	return *this; 
      };
      
      value_type operator[](difference_type i) const { 
	if( i == 0 )
	  return base;
	else
	  throw std::range_error("Tried to index past the only incident edge in a BF Tree structure!");
      };
      reference operator*() { return base; };
      pointer operator->() { return &base; };
    };
    
    struct edge_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef edge_descriptor value_type;
      typedef edge_descriptor* pointer;
      typedef edge_descriptor& reference;
      typedef std::random_access_iterator_tag iterator_category;
      
      edge_descriptor base;
      edge_iterator(const edge_descriptor& aBase = edge_descriptor()) : base(aBase) { };
      out_edge_iterator(vertex_descriptor aSrc, std::size_t aEdgeId) : base(aSrc, aEdgeId) { };
      
      friend bool operator==( const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base == rhs.base; };
      friend bool operator!=( const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base != rhs.base; };
      friend bool operator >( const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base > rhs.base; };
      friend bool operator >=(const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base >= rhs.base; };
      friend bool operator <( const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base < rhs.base; };
      friend bool operator <=(const edge_iterator& lhs, const edge_iterator& rhs) { return lhs.base <= rhs.base; };
      
      edge_iterator& operator++() { 
	++base.edge_index; 
	if(base.edge_index == Arity) {
	  ++base.source_vertex;
	  base.edge_index = 0;
	};
	return *this;
      };
      edge_iterator operator++(int) { 
	out_edge_iterator result(*this); 
	++base.edge_index; 
	if(base.edge_index == Arity) {
	  ++base.source_vertex;
	  base.edge_index = 0;
	};
	return result; 
      };
      edge_iterator& operator--() { 
	if(base.edge_index == 0) {
	  --base.source_vertex;
	  base.edge_index = Arity;
	};
	--base.edge_index; 
	return *this;
      };
      edge_iterator operator--(int) { 
	out_edge_iterator result(*this); 
	if(base.edge_index == 0) {
	  --base.source_vertex;
	  base.edge_index = Arity;
	};
	--base.edge_index; 
	return result; 
      };
      
      friend edge_iterator operator+(const edge_iterator& lhs, difference_type i) {
	edge_descriptor new_base(lhs.base.source_vertex, lhs.base.edge_index + i);
	new_base.source_vertex += new_base.edge_index / Arity;
	new_base.edge_index %= Arity;
	return edge_iterator(new_base);
      };
      friend edge_iterator operator+(difference_type i, const edge_iterator& rhs) {
	edge_descriptor new_base(lhs.base.source_vertex, lhs.base.edge_index + i);
	new_base.source_vertex += new_base.edge_index / Arity;
	new_base.edge_index %= Arity;
	return edge_iterator(new_base);
      };
      friend edge_iterator operator-(const edge_iterator& lhs, difference_type i) {
	edge_descriptor new_base(lhs.base.source_vertex, lhs.base.edge_index - i);
	new_base.source_vertex += new_base.edge_index / Arity;
	new_base.edge_index %= Arity;
	return edge_iterator(new_base);
      };
      friend difference_type operator-(const edge_iterator& lhs, const edge_iterator& rhs) {
	return difference_type( (lhs.base.source_vertex - rhs.source_vertex) * Arity + lhs.base.edge_index - rhs.base.edge_index);
      };
      
      edge_iterator& operator +=(difference_type i) { 
	base.edge_index += i; 
	base.source_vertex += base.edge_index / Arity;
	base.edge_index %= Arity;
	return *this;
      };
      edge_iterator& operator -=(difference_type i) { 
	base.edge_index -= i; 
	base.source_vertex += base.edge_index / Arity;
	base.edge_index %= Arity;
	return *this;
      };
      
      value_type operator[](difference_type i) const { return (*this + i).base; };
      reference operator*() { return base; };
      pointer operator->() { return &base; };
    };
    
    struct vertex_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef vertex_descriptor value_type;
      typedef vertex_descriptor* pointer;
      typedef vertex_descriptor& reference;
      typedef std::random_access_iterator_tag iterator_category;
      
      vertex_descriptor base;
      vertex_iterator(const vertex_descriptor& aBase = vertex_descriptor()) : base(aBase) { };
      
      friend bool operator==( const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base == rhs.base; };
      friend bool operator!=( const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base != rhs.base; };
      friend bool operator >( const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base > rhs.base; };
      friend bool operator >=(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base >= rhs.base; };
      friend bool operator <( const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base < rhs.base; };
      friend bool operator <=(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base <= rhs.base; };
      
      vertex_iterator& operator++() { ++base; return *this; };
      vertex_iterator operator++(int) { out_edge_iterator result(*this); ++base; return result; };
      vertex_iterator& operator--() { --base; return *this; };
      vertex_iterator operator--(int) { out_edge_iterator result(*this); --base; return result; };
      
      friend vertex_iterator operator+(const vertex_iterator& lhs, difference_type i) { return vertex_iterator(lhs.base + i); };
      friend vertex_iterator operator+(difference_type i, const vertex_iterator& rhs) { return vertex_iterator(rhs.base + i); };
      friend vertex_iterator operator-(const vertex_iterator& lhs, difference_type i) { return vertex_iterator(lhs.base - i); };
      friend difference_type operator-(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base - rhs.base; };
      
      vertex_iterator& operator +=(difference_type i) { base += i; return *this; };
      vertex_iterator& operator -=(difference_type i) { base -= i; return *this; };
      
      value_type operator[](difference_type i) const { return base + i; };
      reference operator*() { return base; };
      pointer operator->() { return &base; };
    };
    
    typedef vertex_iterator adjacency_iterator;
    
    typedef boost::directed_tag directed_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;
    
    struct traversal_category : 
      public boost::incidence_graph_tag,
      public boost::adjacency_graph_tag,
      public boost::bidirectional_graph_tag,
      public boost::vertex_list_graph_tag,
      public boost::edge_list_graph_tag { };
    
  private:
    
    container_type m_vertices;
    vertices_size_type m_vertex_count;
    
  public:
    
    /**
     * Construct the D-ary BF-tree with a given reserved depth.
     * \param aDepth The depth of the graph to reserve space for.
     */
    d_ary_bf_tree(vertices_size_type aDepth = 0) : m_vertex_count(0) {
      vertices_size_type vert_count = 1;
      vertices_size_type accum = 1;
      for(vertices_size_type i = 0; i < aDepth; ++i) {
	accum *= Arity;
	vert_count += accum;
      };
      m_vertices.resize(vert_count);
    };
    
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
    
    std::size_t capacity() const { return m_vertices.size(); };
    
    std::size_t depth() const { 
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
    };
    
    void clear() { 
      m_vertices.clear();
      m_vertices.resize(1);
      m_vertex_count = 0;
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
    
    bool is_valid(const vertex_descriptor& v_i) const {
      return ( m_vertices[v_i].out_degree >= 0 );
    };
    
    edge_size_type get_out_degree( const vertex_descriptor& v_i) const {
      if( m_vertices[v_i].out_degree < 0 )
	return 0;
      else
	return m_vertices[v_i].out_degree;
    };
    
    edge_size_type get_in_degree( const vertex_descriptor& v_i) const {
      if(( v_i == 0 ) || ( m_vertices[v_i].out_degree < 0 ))
	return 0;
      else
	return 1;
    };
    
    vertex_descriptor add_child(const vertex_descriptor& v) {
      if( m_vertices[v].out_degree < 0 ) 
	throw std::range_error("Cannot add child-node to an empty node!");
      if( m_vertices[v].out_degree == Arity ) 
	throw std::range_error("Cannot add child-node to a full node!");
      vertex_descriptor result = Arity * v + 1 + m_vertices[v].out_degree;
      if( result >= m_vertices.size() )
	m_vertices.resize(result + 1);
      m_vertices[result].out_degree = 0;
      ++(m_vertices[v].out_degree);
      
    };
    
};


/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  source( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor& e,
	  const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  return e.source_vertex;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor
  target( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor& e,
	  const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  return Arity * e.source_vertex + 1 + e.edge_index;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator >
  out_edges( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	     const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::out_edge_iterator OutIter;
  return std::make_pair(OutIter(v,0),OutIter(v,g.get_out_degree(v)));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
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
std::size_t
  in_degree( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
	     const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.get_in_degree(v);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
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
std::pair<
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::adjacency_iterator,
 typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::adjacency_iterator >
  in_edges( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
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
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertices_size_type
  num_vertices( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.capacity();
};


/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
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
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertices_size_type
  num_edges( const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>& g) {
  return g.capacity() * Arity;
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor
  edge( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor&,
	const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
        const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor Edge;
  return Edge((v - 1) / Arity, (v - 1) % Arity);
};




};

};


#endif


















