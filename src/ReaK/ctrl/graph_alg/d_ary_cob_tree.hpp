/**
 * \file d_ary_cob_tree.hpp
 * 
 * This library provides a class that implements a D-Ary Cache-Oblivious B-tree (COB-tree) that is tailored 
 * to store elements of a tree as if they were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues, much than a classic breadth-first layout because a recursive 
 * break-down into sub-trees allow for better cache performance for any memory hierarchy (i.e. 
 * cache-oblivious). Ideally, for the least amount of wasted memory, the tree 
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

#ifndef REAK_D_ARY_COB_TREE_HPP
#define REAK_D_ARY_COB_TREE_HPP

#include <boost/graph/properties.hpp>
#include <boost/config.hpp>

#include <vector>
#include <stdexcept>
#include <map>
#include <unordered_map>
#include <iterator>

namespace ReaK {

namespace graph {
  
  
namespace detail {
  
  
  template <std::size_t Arity, std::size_t CuttingDepth>
  struct get_cob_block_property {
    BOOST_STATIC_CONSTANT(std::size_t, fanout = (Arity * get_cob_block_property<Arity, CuttingDepth - 1>::fanout));
    BOOST_STATIC_CONSTANT(std::size_t, vertex_count = (get_cob_block_property<Arity, CuttingDepth - 1>::fanout + get_cob_block_property<Arity, CuttingDepth - 1>::vertex_count));
  };
  
  template <std::size_t Arity>
  struct get_cob_block_property<Arity, 1> {
    BOOST_STATIC_CONSTANT(std::size_t, fanout = Arity);
    BOOST_STATIC_CONSTANT(std::size_t, vertex_count = 1);
  };
  
  template <std::size_t Arity>
  struct get_cob_block_property<Arity, 0> {
    char cannot_instantiate_this_specialization[0];
  };
  
  
};


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
          typename EdgeProperties = boost::no_property,
	  std::size_t CuttingDepth = 8 >
class d_ary_cob_tree
{
  public:
    typedef d_ary_cob_tree<VertexProperties, Arity, EdgeProperties, CuttingDepth> self;
    
    BOOST_STATIC_CONSTANT(std::size_t, BlockFanout = (detail::get_cob_block_property< Arity, CuttingDepth >::fanout));
    BOOST_STATIC_CONSTANT(std::size_t, BlockVertexCount = (detail::get_cob_block_property< Arity, CuttingDepth >::vertex_count));
    
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
    
    typedef std::unordered_map< std::size_t, container_type > block_map_type;
    
    typedef typename container_type::size_type vertices_size_type;
    typedef typename container_type::difference_type vertex_index_type;
    typedef vertices_size_type edges_size_type;
    typedef vertices_size_type degree_size_type;
    
    struct vertex_descriptor {
      
      std::size_t block_id;
      vertex_index_type vertex_id;
      
      vertex_descriptor(std::size_t aBlockID = 0, 
			vertex_index_type aVertexID = 0) : 
			block_id(aBlockID), vertex_id(aVertexID) { };
      
      friend bool operator==( const vertex_descriptor& lhs, const vertex_descriptor& rhs) { 
	return ((lhs.block_id == rhs.block_id) && (lhs.vertex_id == rhs.vertex_id)); 
      };
      friend bool operator!=( const vertex_descriptor& lhs, const vertex_descriptor& rhs) { 
	return ((lhs.block_id != rhs.block_id) || (lhs.vertex_id != rhs.vertex_id)); 
      };
      friend bool operator <( const vertex_descriptor& lhs, const vertex_descriptor& rhs) {
	if( lhs.block_id == rhs.block_id )
	  return ( lhs.vertex_id < rhs.vertex_id );
	else 
	  return ( lhs.block_id < rhs.block_id );
      };
      friend bool operator<=( const vertex_descriptor& lhs, const vertex_descriptor& rhs) {
	if( lhs.block_id == rhs.block_id )
	  return ( lhs.vertex_id <= rhs.vertex_id );
	else 
	  return ( lhs.block_id < rhs.block_id );
      };
      friend bool operator >( const vertex_descriptor& lhs, const vertex_descriptor& rhs) {
	if( lhs.block_id == rhs.block_id )
	  return ( lhs.vertex_id > rhs.vertex_id );
	else 
	  return ( lhs.block_id > rhs.block_id );
      };
      friend bool operator>=( const vertex_descriptor& lhs, const vertex_descriptor& rhs) {
	if( lhs.block_id == rhs.block_id )
	  return ( lhs.vertex_id >= rhs.vertex_id );
	else 
	  return ( lhs.block_id > rhs.block_id );
      };
      
      bool is_root() const { 
	return ((block_id == 0) && (vertex_id == 0));
      };
      
      vertex_descriptor get_child(std::size_t i) const {
	vertex_index_type result = Arity * vertex_id + 1 + i;
	if( result >= vertex_index_type(BlockVertexCount) ) {
	  result -= BlockVertexCount;
	  return vertex_descriptor(BlockFanout * block_id + 1 + result, 0);
	};
	return vertex_descriptor(block_id, result);
      };
      
      vertex_descriptor get_parent() const {
	if(vertex_id == 0)
	  return vertex_descriptor((block_id - 1) / BlockFanout, ( ( (block_id - 1) % BlockFanout ) / Arity + BlockVertexCount ) - BlockFanout / Arity);
	return vertex_descriptor(block_id, (vertex_id - 1) / Arity);
      };
      
      std::ptrdiff_t get_in_edge() const {
	if(vertex_id == 0)
	  return ( (block_id - 1) % BlockFanout ) % Arity;
	return (vertex_id - 1) % Arity;
      };
      
    };
    
    struct edge_descriptor {
      vertex_descriptor source_vertex;
      std::ptrdiff_t edge_index;
      edge_descriptor(vertex_descriptor aSrc = vertex_descriptor(), std::size_t aEdgeId = 0) : source_vertex(aSrc), edge_index(aEdgeId) { };
      
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
      in_edge_iterator() : base( vertex_descriptor(), -1) { };
      in_edge_iterator(const vertex_descriptor& aBase) : 
          base(aBase.get_parent(), aBase.get_in_edge()) { };
      
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
	  throw std::range_error("Tried to index past the only incident edge in a COB-Tree structure!");
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
      
      edge_iterator& operator +=(difference_type i) { 
	base.edge_index += i; 
	base.source_vertex.vertex_id += base.edge_index / Arity;
	base.source_vertex.block_id += base.source_vertex.vertex_id / BlockVertexCount;
	base.source_vertex.vertex_id %= BlockVertexCount;
	base.edge_index %= Arity;
	return *this;
      };
      edge_iterator& operator -=(difference_type i) { 
	base.edge_index -= i; 
	base.source_vertex.vertex_id += base.edge_index / Arity;
	base.source_vertex.block_id += base.source_vertex.vertex_id / BlockVertexCount;
	base.source_vertex.vertex_id %= BlockVertexCount;
	base.edge_index %= Arity;
	return *this;
      };
      
      edge_iterator& operator++() { 
	return *this += 1;
      };
      edge_iterator operator++(int) { 
	edge_iterator result(*this); 
	*this += 1;
	return result; 
      };
      edge_iterator& operator--() {
	return *this -= 1;
      };
      edge_iterator operator--(int) { 
	edge_iterator result(*this); 
	*this -= 1;
	return result; 
      };
      
      friend edge_iterator operator+(const edge_iterator& lhs, difference_type i) {
	edge_iterator result(lhs);
	result += i;
	return result;
      };
      friend edge_iterator operator+(difference_type i, const edge_iterator& rhs) {
	edge_iterator result(rhs);
	result += i;
	return result;
      };
      friend edge_iterator operator-(const edge_iterator& lhs, difference_type i) {
	edge_iterator result(lhs);
	result -= i;
	return result;
      };
      friend difference_type operator-(const edge_iterator& lhs, const edge_iterator& rhs) {
	return (difference_type(lhs.base.source_vertex.block_id) - difference_type(rhs.base.source_vertex.block_id)) * BlockVertexCount
	     + (difference_type(lhs.base.source_vertex.vertex_id) - difference_type(rhs.source_vertex.vertex_id)) * Arity
	      + difference_type(lhs.base.edge_index) - difference_type(rhs.base.edge_index);
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
      
      vertex_iterator& operator +=(difference_type i) { 
	base.vertex_id += i; 
	base.block_id += base.vertex_id / BlockVertexCount;
	base.vertex_id %= BlockVertexCount;
	return *this; 
      };
      vertex_iterator& operator -=(difference_type i) { 
	base.vertex_id -= i; 
	base.block_id += base.vertex_id / BlockVertexCount;
	base.vertex_id %= BlockVertexCount;
	return *this; 
      };
      
      vertex_iterator& operator++() { 
	return *this += 1;
      };
      vertex_iterator operator++(int) { 
	vertex_iterator result(*this); 
        *this += 1;
	return result; 
      };
      vertex_iterator& operator--() { 
	return *this -= 1;
      };
      vertex_iterator operator--(int) { 
	vertex_iterator result(*this); 
	*this -= 1;
	return result; 
      };
      
      friend vertex_iterator operator+(const vertex_iterator& lhs, difference_type i) { 
	vertex_iterator result(lhs);
	result += i;
	return result; 
      };
      friend vertex_iterator operator+(difference_type i, const vertex_iterator& rhs) { 
	vertex_iterator result(rhs);
	result += i;
	return result; 
      };
      friend vertex_iterator operator-(const vertex_iterator& lhs, difference_type i) { 
	vertex_iterator result(lhs);
	result -= i;
	return result; 
      };
      friend difference_type operator-(const vertex_iterator& lhs, const vertex_iterator& rhs) { 
	return (difference_type(lhs.base.block_id) - difference_type(rhs.base.block_id)) * BlockVertexCount
	       + difference_type(lhs.base.vertex_id) - difference_type(rhs.base.vertex_id); 
      };
      
      reference operator*() { return base; };
      value_type operator[](difference_type i) const { return *(*this + i); };
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
    
    mutable block_map_type m_vertices;  //NOTE This member must be mutable to account for the fact that indexing in the map might add an entry to it.
    vertices_size_type m_vertex_count;
    
  public:
    
    static vertex_descriptor null_vertex() { 
      return vertex_descriptor(reinterpret_cast<std::size_t>(-1));
    };
    
    /**
     * Construct the D-ary BF-tree with a given reserved depth.
     * \param aDepth The depth of the graph to reserve space for.
     */
    d_ary_cob_tree(vertices_size_type aDepth = 0) : m_vertex_count(0) {
      std::size_t block_id = 0;
      std::size_t block_count = 1;
      while(true) {
        vertices_size_type vert_count = 1;
        vertices_size_type accum = 1;
        for(vertices_size_type i = 0; i < aDepth; ++i) {
	  accum *= Arity;
	  vert_count += accum;
        };
	for(std::size_t i = 0; i < block_count; ++i, ++block_id)
          m_vertices[block_id].resize(vert_count);
	block_count *= BlockFanout;
	if(aDepth <= CuttingDepth)
	  break;
	aDepth -= CuttingDepth;
      };
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
    
    std::size_t capacity() const { return m_vertices.size() * BlockVertexCount; };
    
    std::size_t depth() const { 
      vertices_size_type block_level = 0;
      vertices_size_type accum = 1;
      vertices_size_type depth_count = 0;
      while(true) {
	vertices_size_type max_vertex_count = 0;
	for(vertices_size_type i = 0; i < accum; ++i) {
	  if( m_vertices[i].size() > max_vertex_count ) 
	    max_vertex_count = m_vertices[i].size();
	};
	if(max_vertex_count + BlockFanout / Arity > BlockVertexCount)
	  depth_count += CuttingDepth;
	else {
	  vertices_size_type accum2 = 1;
	  for(vertices_size_type vert_count = 1; vert_count < max_vertex_count; ++depth_count) {
	    accum2 *= Arity;
	    vert_count += accum2;
          };
	  return depth_count;
	};
	accum *= BlockFanout;
	++block_level;
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
      m_vertices[0].resize(1);
      m_vertex_count = 0;
    };
    
    vertex_property_type& operator[]( const vertex_descriptor& v_i) {
      return m_vertices[v_i.block_id][v_i.vertex_id].v;
    };
    const vertex_property_type& operator[]( const vertex_descriptor& v_i) const {
      return m_vertices[v_i.block_id][v_i.vertex_id].v;
    };
    edge_property_type& operator[]( const edge_descriptor& e_i) {
      return m_vertices[e_i.source_vertex.block_id][e_i.source_vertex.vertex_id].e[e_i.edge_index];
    };
    const edge_property_type& operator[]( const edge_descriptor& e_i) const {
      return m_vertices[e_i.source_vertex.block_id][e_i.source_vertex.vertex_id].e[e_i.edge_index];
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
      return (v_i.vertex_id < vertex_index_type( m_vertices[v_i.block_id].size() )) && 
             (m_vertices[v_i.block_id][v_i.vertex_id].out_degree >= 0);
    };
    
    edges_size_type get_out_degree( const vertex_descriptor& v_i) const {
      if( (v_i.vertex_id >= vertex_index_type(m_vertices[v_i.block_id].size())) || ((m_vertices[v_i.block_id][v_i.vertex_id].out_degree) < 0) )
	return 0;
      else {
	edges_size_type result = 0;
	for(edges_size_type i = 0; i < edges_size_type(m_vertices[v_i.block_id][v_i.vertex_id].out_degree); ++i) {
	  vertex_descriptor v_c = v_i.get_child(i);
	  if( is_valid(v_c) )
	    ++result;
	};
	return result;
      };
    };
    
    edges_size_type get_raw_out_degree( const vertex_descriptor& v_i) const {
      if( (v_i.vertex_id >= vertex_index_type(m_vertices[v_i.block_id].size())) || ((m_vertices[v_i.block_id][v_i.vertex_id].out_degree) < 0) )
	return 0;
      else
	return m_vertices[v_i.block_id][v_i.vertex_id].out_degree;
    };
    
    edges_size_type get_in_degree( const vertex_descriptor& v_i) const {
      if(( ( v_i.block_id == 0 ) && ( v_i.vertex_id == 0 ) ) || 
	 ( !is_valid(v_i) ))
	return 0;
      else
	return 1;
    };
    
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     const vertex_property_type& vp = vertex_property_type(), 
							     const edge_property_type& ep = edge_property_type()) {
      if( (v.vertex_id >= vertex_index_type(m_vertices[v.block_id].size())) || ((m_vertices[v.block_id][v.vertex_id].out_degree) < 0) )
	throw std::range_error("Cannot add child-node to an empty node!");
      int new_edge = 0;
      for(; new_edge < m_vertices[v.block_id][v.vertex_id].out_degree; ++new_edge) {
	vertex_descriptor v_c = v.get_child(new_edge);
	if( (m_vertices[v_c.block_id][v_c.vertex_id].out_degree) < 0 )
	  break;
      };
      if( new_edge == Arity ) 
	throw std::range_error("Cannot add child-node to a full node!");
      vertex_descriptor result = v.get_child(new_edge);
      if( result.vertex_id >= vertex_index_type(m_vertices[result.block_id].size()) )
	m_vertices[result.block_id].resize(result.vertex_id + 1);
      m_vertices[result.block_id][result.vertex_id].out_degree = 0;
      m_vertices[result.block_id][result.vertex_id].v = vp;
      m_vertices[v.block_id][v.vertex_id].e[new_edge] = ep;
      if( new_edge == m_vertices[v.block_id][v.vertex_id].out_degree )
        ++(m_vertices[v.block_id][v.vertex_id].out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
							     vertex_property_type&& vp, 
							     edge_property_type&& ep = edge_property_type()) {
      if( (v.vertex_id >= vertex_index_type(m_vertices[v.block_id].size())) || ((m_vertices[v.block_id][v.vertex_id].out_degree) < 0) )
	throw std::range_error("Cannot add child-node to an empty node!");
      int new_edge = 0;
      for(; new_edge < m_vertices[v.block_id][v.vertex_id].out_degree; ++new_edge) {
	vertex_descriptor v_c = v.get_child(new_edge);
	if( (m_vertices[v_c.block_id][v_c.vertex_id].out_degree) < 0 )
	  break;
      };
      if( new_edge == Arity ) 
	throw std::range_error("Cannot add child-node to a full node!");
      vertex_descriptor result = v.get_child(new_edge);
      if( result.vertex_id >= vertex_index_type(m_vertices[result.block_id].size()) )
	m_vertices[result.block_id].resize(result.vertex_id + 1);
      m_vertices[result.block_id][result.vertex_id].out_degree = 0;
      m_vertices[result.block_id][result.vertex_id].v = std::move(vp);
      m_vertices[v.block_id][v.vertex_id].e[new_edge] = std::move(ep);
      if( new_edge == m_vertices[v.block_id][v.vertex_id].out_degree )
        ++(m_vertices[v.block_id][v.vertex_id].out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
    };
#endif
    
    void update_out_degree(const vertex_descriptor& v) {
      if( !is_valid(v) )
	return;  // vertex is already updated.
      for( int i = Arity; i > 0; --i) {
	vertex_descriptor u = v.get_child( i - 1 );
	if( is_valid(u) ) {
	  m_vertices[v.block_id][v.vertex_id].out_degree = i;
	  return;
	};
      };
      m_vertices[v.block_id][v.vertex_id].out_degree = 0;
    };
    
    void remove_branch(const vertex_descriptor& v) {
      if( !is_valid(v) )
	return;  // vertex is already deleted.
      --m_vertex_count;
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      int max_child = m_vertices[v.block_id][v.vertex_id].out_degree;
      for( int i = 0; i < max_child; ++i)
	remove_branch(v.get_child(i));
      m_vertices[v.block_id][v.vertex_id].out_degree = -1;
      if( v != vertex_descriptor(0,0) )  // if the node is not the root one, then update the out-degree of the parent node:
	update_out_degree( v.get_parent() );
      // remove empty vertices from the end of the container:
      if( v.vertex_id == vertex_index_type(m_vertices[v.block_id].size()) - 1 ) {
	vertex_iterator v_it(v);
        while( (v_it->vertex_id > 0) && ( (m_vertices[v_it->block_id][v_it->vertex_id].out_degree) < 0 ) )
	  --v_it;
	//if( (m_vertices[v_it->block_id][v_it->vertex_id].out_degree) >= 0 )
	  ++v_it;
        m_vertices[v_it->block_id].erase(m_vertices[v_it->block_id].begin() + v_it->vertex_id, m_vertices[v_it->block_id].end());
      };
    };
    
    
    template <typename OutputIter>
    void remove_branch(vertex_descriptor v, OutputIter it_out) {
      if( !is_valid(v) )
	return;  // vertex is already deleted.
      --m_vertex_count; 
#ifdef RK_ENABLE_CXX0X_FEATURES
      *(it_out++) = std::move(m_vertices[v.block_id][v.vertex_id].v);
#else
      *(it_out++) = m_vertices[v.block_id][v.vertex_id].v;
#endif
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      int max_child = m_vertices[v.block_id][v.vertex_id].out_degree;
      for( int i = 0; i < max_child; ++i)
	remove_branch(v.get_child(i),it_out);
      m_vertices[v.block_id][v.vertex_id].out_degree = -1;
      if( v != vertex_descriptor(0,0) )  // if the node is not the root one, then update the out-degree of the parent node:
	update_out_degree( v.get_parent() );
      // remove empty vertices from the end of the container:
      if( v.vertex_id == vertex_index_type(m_vertices[v.block_id].size()) - 1 ) {
	vertex_iterator v_it(v);
        while( (v_it->vertex_id > 0) && ( (m_vertices[v_it->block_id][v_it->vertex_id].out_degree) < 0 ) )
	  --v_it;
	//if( (m_vertices[v_it->block_id][v_it->vertex_id].out_degree) >= 0 )
	  ++v_it;
        m_vertices[v_it->block_id].erase(m_vertices[v_it->block_id].begin() + v_it->vertex_id, m_vertices[v_it->block_id].end());
      };
    };
    
    vertex_descriptor get_root_vertex() const {
      return vertex_descriptor(0,0); 
    };
    
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(m_vertices[0].size() == 0)
	m_vertices[0].resize(1);
      if(m_vertices[0][0].out_degree >= 0)
	remove_branch(vertex_descriptor(0,0));
      m_vertices[0][0].out_degree = 0;
      m_vertices[0][0].v = vp;
      ++m_vertex_count;
      return vertex_descriptor(0,0);
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    vertex_descriptor create_root_vertex(vertex_property_type&& vp) {
      if(m_vertices[0].size() == 0)
	m_vertices[0].resize(1);
      if(m_vertices[0][0].out_degree >= 0)
	remove_branch(vertex_descriptor(0,0));
      m_vertices[0][0].out_degree = 0;
      m_vertices[0][0].v = std::move(vp);
      ++m_vertex_count;
      return vertex_descriptor(0,0);
    };
#endif
    
    
};



template <std::size_t Arity = 2, std::size_t CuttingDepth = 8>
struct d_ary_cob_tree_storage { };


template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity, std::size_t CuttingDepth>
struct tree_storage<VertexDescriptor, EdgeDescriptor, d_ary_cob_tree_storage<Arity, CuttingDepth> > {
  typedef d_ary_cob_tree<VertexDescriptor, Arity, EdgeDescriptor, CuttingDepth> type;
};




template <std::size_t Arity, std::size_t CuttingDepth>
struct tree_storage_traits< d_ary_cob_tree_storage<Arity, CuttingDepth> > {
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
  typedef typename d_ary_cob_tree_storage<int,Arity,boost::no_property,CuttingDepth>::vertex_descriptor vertex_descriptor;  // the value-type doesn't affect the vertex_descriptor type (int is a dummy type here).
  typedef typename d_ary_cob_tree_storage<int,Arity,boost::no_property,CuttingDepth>::edge_descriptor edge_descriptor;  // the value-type doesn't affect the edge_descriptor type (int is a dummy type here).
  typedef std::size_t edges_size_type;
  
};



/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  source( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor& e,
	  const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  return e.source_vertex;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  target( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor& e,
	  const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  return e.source_vertex.get_child(e.edge_index);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair<
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::out_edge_iterator,
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::out_edge_iterator >
  out_edges( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
	     const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::out_edge_iterator OutIter;
  return std::make_pair(OutIter(v,0),OutIter(v,g.get_raw_out_degree(v)));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::size_t
  out_degree( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
	      const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair<
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::in_edge_iterator,
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::in_edge_iterator >
  in_edges( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
	    const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::in_edge_iterator InIter;
  return std::make_pair(InIter(v),InIter());
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::size_t
  in_degree( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
             const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.get_in_degree(v);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::size_t
  degree( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
          const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.get_in_degree(v) + g.get_out_degree(v);
};


/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair<
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::adjacency_iterator,
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::adjacency_iterator >
  adjacent_vertices( typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor v,
	             const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::adjacency_iterator AdjIter;
  v = v.get_parent();
  return std::make_pair(AdjIter(v.get_child(0)),
			AdjIter(v.get_child(Arity)));
};


/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair<
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator,
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator >
  vertices( const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator VIter;
  return std::make_pair(VIter(0),
			VIter(g.capacity()));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertices_size_type
  num_vertices( const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.size();
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair<
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_iterator,
 typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_iterator >
  edges( const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_iterator EIter;
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor Vertex;
  return std::make_pair(EIter(Vertex(0,0), 0),
			EIter(Vertex(0,0), 0) + g.capacity());
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertices_size_type
  num_edges( const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.size() - 1;
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor
  edge( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor&,
	const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
        const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor Edge;
  return Edge(v.get_parent(), v.get_in_edge());
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  get_root_vertex( const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.get_root_vertex();
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator >
  child_vertices( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
                  const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>&) {
  typedef typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_iterator VIter;
  return std::make_pair(VIter(Arity * v + 1),VIter(Arity * (v + 1)));
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  create_root( d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.create_root_vertex();
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor >
  add_child_vertex( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.add_child(v);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
void remove_branch( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.remove_branch(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  create_root( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type& vp, 
	       d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.create_root_vertex(vp);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor
  create_root( typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type&& vp, 
	       d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.create_root_vertex(std::move(vp));
};
#endif

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor >
  add_child_vertex( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
		    const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type& vp,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.add_child(v,vp);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor >
  add_child_vertex( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
		    const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type& vp,
		    const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_property_type& ep,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.add_child(v,vp,ep);
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor >
  add_child_vertex( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
		    typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type&& vp,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.add_child(v,std::move(vp));
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
std::pair< 
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor,
typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor >
  add_child_vertex( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
		    typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_property_type&& vp,
		    typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_property_type&& ep,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.add_child(v,std::move(vp),std::move(ep));
};
#endif

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth,
	  typename OutputIter>
inline
void remove_branch( const typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor& v,
		    OutputIter it_out,
                    d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.remove_branch(v,it_out);
};


/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth >
inline
bool is_vertex_valid( typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::vertex_descriptor u,
                      const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.is_valid(u);
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties,
	  std::size_t CuttingDepth>
inline
bool is_edge_valid( typename d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>::edge_descriptor e,
                    const d_ary_cob_tree<VertexProperties,Arity,EdgeProperties,CuttingDepth>& g) {
  return g.is_valid(target(e,g));
};



};

};


#endif


















