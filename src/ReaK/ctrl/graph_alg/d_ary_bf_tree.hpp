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

#include <boost/graph/properties.hpp>

#include <vector>
#include <stdexcept>
#include <map>
#include <iterator>

#include "bgl_tree_adaptor.hpp"

namespace ReaK {

namespace graph {


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
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    struct value_type {
      int out_degree;
      vertex_property_type v;
      edge_property_type e[Arity];
      
      value_type() : out_degree(-1) { };
    };
    
    typedef std::vector< value_type > container_type;
    
    typedef typename container_type::difference_type vertex_descriptor;
    
    typedef typename container_type::size_type vertices_size_type;
    typedef vertices_size_type edges_size_type;
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
      
      value_type operator[](difference_type i) const { return edge_descriptor(base.source_vertex, base.edge_index + i); };
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
      
      in_edge_iterator& operator++() { base.source_vertex = -1; base.edge_index = -1; return *this; };
      in_edge_iterator operator++(int) { in_edge_iterator result(*this); base.source_vertex = -1; base.edge_index = -1; return result; };
      in_edge_iterator& operator--() { base.source_vertex = -1; base.edge_index = -1; return *this; };
      in_edge_iterator operator--(int) { in_edge_iterator result(*this); base.source_vertex = -1; base.edge_index = -1; return result; };
      
      friend in_edge_iterator operator+(const in_edge_iterator& lhs, difference_type i) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return lhs;
      };
      friend in_edge_iterator operator+(difference_type i, const in_edge_iterator& rhs) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return rhs;
      };
      friend in_edge_iterator operator-(const in_edge_iterator& lhs, difference_type i) {
	if( i != 0 )
	  return in_edge_iterator();
	else
	  return lhs;
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
	if( i != 0 ) {
	  base.source_vertex = -1;
	  base.edge_index = -1;
	};
	return *this; 
      };
      in_edge_iterator& operator -=(difference_type i) { 
	if( i != 0 ) {
	  base.source_vertex = -1;
	  base.edge_index = -1;
	};
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
      edge_iterator(vertex_descriptor aSrc, std::size_t aEdgeId) : base(aSrc, aEdgeId) { };
      
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
	edge_iterator result(*this); 
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
	edge_iterator result(*this); 
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
	edge_descriptor new_base(rhs.base.source_vertex, rhs.base.edge_index + i);
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
      
      value_type operator[](difference_type i) const { return (*this + i).base; };
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
      vertex_iterator operator++(int) { vertex_iterator result(*this); ++base; return result; };
      vertex_iterator& operator--() { --base; return *this; };
      vertex_iterator operator--(int) { vertex_iterator result(*this); --base; return result; };
      
      friend vertex_iterator operator+(const vertex_iterator& lhs, difference_type i) { return vertex_iterator(lhs.base + i); };
      friend vertex_iterator operator+(difference_type i, const vertex_iterator& rhs) { return vertex_iterator(rhs.base + i); };
      friend vertex_iterator operator-(const vertex_iterator& lhs, difference_type i) { return vertex_iterator(lhs.base - i); };
      friend difference_type operator-(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base - rhs.base; };
      
      vertex_iterator& operator +=(difference_type i) { base += i; return *this; };
      vertex_iterator& operator -=(difference_type i) { base -= i; return *this; };
      
      value_type operator[](difference_type i) const { return base + i; };
      reference operator*() { return base; };
      pointer operator->() { return &base; };
    };
    
    typedef vertex_iterator adjacency_iterator;
    typedef vertex_iterator child_vertex_iterator;
    
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
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() { 
      return reinterpret_cast<vertex_descriptor>(-1);
    };
    
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
     * Checks if the tree is empty.
     * \return True if the tree is empty.
     */
    bool empty() const { return m_vertex_count == 0; };
    
    /**
     * Returns the size of the tree (the number of vertices it contains).
     * \return The size of the tree (the number of vertices it contains).
     */
    std::size_t size() const { return m_vertex_count; };
    
    /**
     * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
     * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
     */
    std::size_t capacity() const { return m_vertices.size(); };
    
    /**
     * Returns the depth of the tree.
     * \return The depth of the tree.
     */
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
    
    /**
     * Standard swap function.
     */
    friend
    void swap(self& lhs, self& rhs) {
      using std::swap;
      lhs.m_vertices.swap(rhs.m_vertices);
      swap(lhs.m_vertex_count, rhs.m_vertex_count);
    };
    
    /**
     * Clears the tree of all vertices and edges.
     */
    void clear() { 
      m_vertices.clear();
      m_vertices.resize(1);
      m_vertex_count = 0;
    };
    
    /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[]( const vertex_descriptor& v_i) {
      return m_vertices[v_i].v;
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[]( const vertex_descriptor& v_i) const {
      return m_vertices[v_i].v;
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[]( const edge_descriptor& e_i) {
      return m_vertices[e_i.source_vertex].e[e_i.edge_index];
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[]( const edge_descriptor& e_i) const {
      return m_vertices[e_i.source_vertex].e[e_i.edge_index];
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
      return (std::size_t(v_i) < m_vertices.size()) && ( m_vertices[v_i].out_degree >= 0 );
    };
    
    /**
     * Returns the raw out-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The out-degree of the given vertex descriptor.
     */
    edges_size_type get_raw_out_degree( const vertex_descriptor& v_i) const {
      if( m_vertices[v_i].out_degree < 0 )
	return 0;
      else 
	return m_vertices[v_i].out_degree;
    };
    
    /**
     * Returns the out-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The out-degree of the given vertex descriptor.
     */
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
    
    /**
     * Returns the in-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
     */
    edges_size_type get_in_degree( const vertex_descriptor& v_i) const {
      if(( v_i == 0 ) || ( m_vertices[v_i].out_degree < 0 ))
	return 0;
      else
	return 1;
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
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
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
      if( (v >= vertex_descriptor(m_vertices.size())) || (m_vertices[v].out_degree < 0) )
	return it_out;  // vertex is already deleted.
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
	it_out = remove_branch(Arity * v + 1 + i, it_out);
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
      return it_out;
    };
    
    /**
     * Returns the vertex-descriptor of the root of the tree.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor get_root_vertex() const {
      return 0; 
    };
    
    /**
     * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
     * \param vp The vertex-property to assign to the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(m_vertices[0].out_degree >= 0)
	remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].v = vp;
      ++m_vertex_count;
      return 0;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
     * \param vp The vertex-property to move into the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
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


/**
 * This is the tree-storage specifier for a D-ary BF-tree of a given Arity.
 */
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
      boost::directed_tag, boost::undirected_tag
    >::type
  >::type directed_category;
  
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  
  typedef std::size_t vertices_size_type;
  typedef void* vertex_ptr;
  typedef typename d_ary_bf_tree<int,Arity>::vertex_descriptor vertex_descriptor;  // the value-type doesn't affect the vertex_descriptor type (int is a dummy type here).
  typedef typename d_ary_bf_tree<int,Arity>::edge_descriptor edge_descriptor;  // the value-type doesn't affect the edge_descriptor type (int is a dummy type here).
  typedef std::size_t edges_size_type;
  
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
std::pair<
  typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor,
  bool >
  edge( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor&,
	const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
        const d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>&) {
  typedef typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::edge_descriptor Edge;
  return std::make_pair(Edge((v - 1) / Arity, (v - 1) % Arity),true);
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
  return std::make_pair(VIter(Arity * v + 1),VIter(Arity * (v + 1) + 1));
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
OutputIter remove_branch( const typename d_ary_bf_tree<VertexProperties,Arity,EdgeProperties>::vertex_descriptor& v,
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


















