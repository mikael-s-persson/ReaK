/**
 * \file adj_list_tree_overlay.hpp
 * 
 * This library provides a class 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_ADJ_LIST_OVERLAY_HPP
#define REAK_ADJ_LIST_OVERLAY_HPP

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>
#include "bgl_more_property_tags.hpp"

#include <map>
#include <unordered_map>
#include <vector>
#include <queue>

#include "d_ary_bf_tree.hpp"

#include <boost/graph/adjacency_list.hpp>



namespace ReaK {

namespace graph {


template <typename AdjListOnTreeType>
class alt_tree_view;  // forward-declaration.


template <typename AdjListOnTreeType>
class alt_graph_view;  // forward-declaration.

  
namespace detail {
  
/*
 * This class template implements an adjacency-list (in BGL-style) that over-shadows (i.e., a multi-graph)
 * a tree to store the underlying vertices according to the layout that the tree uses (e.g., B-tree (BFL), 
 * or COB-tree). An additional benefit of this class template is that the B-tree can act without being an 
 * indirect to a graph, and the two graphs are always kept consistent, which would be difficult to do externally.
 */
template <typename AdjListType,
	  typename VertexProperty,
	  typename TreeVertexProperty, 
	  typename TreeEdgeProperty,
	  typename TreeStorageTag>
class adj_list_on_tree {
  public:
    typedef adj_list_on_tree<AdjListType, VertexProperty, TreeVertexProperty, TreeEdgeProperty, TreeStorageTag> self;
    
    typedef AdjListType adj_list_type;
    
    typedef typename boost::graph_traits< adj_list_type >::vertex_descriptor adj_vertex_descriptor;
    
    struct tree_vertex_properties {
      adj_vertex_descriptor partner_node;
      TreeVertexProperty tree_data;
      VertexProperty user_data;
      
      friend bool operator <(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node < rhs.partner_node;
      };
      friend bool operator <=(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node <= rhs.partner_node;
      };
      friend bool operator >(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node > rhs.partner_node;
      };
      friend bool operator >=(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node >= rhs.partner_node;
      };
      friend bool operator ==(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node == rhs.partner_node;
      };
      friend bool operator !=(const tree_vertex_properties& lhs, const tree_vertex_properties& rhs) {
	return lhs.partner_node != rhs.partner_node;
      };
    };
    typedef TreeEdgeProperty tree_edge_properties;
    
    typedef typename tree_storage< tree_vertex_properties, tree_edge_properties, TreeStorageTag>::type tree_type;
    
    typedef VertexProperty adj_vertex_bundled;
    typedef typename adj_list_type::edge_bundled adj_edge_bundled;
    typedef TreeVertexProperty tree_vertex_bundled;
    typedef TreeEdgeProperty tree_edge_bundled;
    
  private:
    adj_list_type m_adj_list;
    tree_type m_tree;
    
    typedef std::priority_queue< 
      adj_vertex_descriptor,
      std::vector< adj_vertex_descriptor >,
      std::greater< adj_vertex_descriptor > > pnode_avail_queue;
    
    pnode_avail_queue m_available_pnodes;
    
  public:
    friend class alt_tree_view< self >;
    friend class alt_graph_view< self >;
    
};

};






template <typename OutEdgeListS = boost::vecS,
	  typename DirectedS = boost::directedS,
	  typename VertexProperty = boost::no_property,
	  typename AdjEdgeProperty = boost::no_property,
	  typename AdjEdgeListS = boost::vecS,
	  typename TreeStorageTag = ReaK::graph::d_ary_bf_tree_storage<2> >
struct adj_list_on_tree_tag {
  
  typedef typename tree_storage_traits< TreeStorageTag >::vertex_descriptor tree_vertex_descriptor;
  typedef typename tree_storage_traits< TreeStorageTag >::edge_descriptor tree_edge_descriptor;
    
  struct adj_vertex_properties {
    tree_vertex_descriptor tree_vertex;
  };
    
  typedef AdjEdgeProperty adj_edge_properties;
    
  typedef boost::adjacency_list< OutEdgeListS, boost::vecS, DirectedS, 
                                 adj_vertex_properties, adj_edge_properties,
				 boost::no_property, AdjEdgeListS > adj_list_type;
  
  template <typename TreeVertexProperty, typename TreeEdgeProperty>
  struct alt {
    typedef detail::adj_list_on_tree<adj_list_type, VertexProperty, 
                                     TreeVertexProperty, TreeEdgeProperty, TreeStorageTag> type;
  };
  
};



template <typename TreeVertexProperty, 
          typename TreeEdgeProperty, 
	  typename OutEdgeListS,
	  typename DirectedS,
	  typename VertexProperty,
	  typename AdjEdgeProperty,
	  typename AdjEdgeListS,
	  typename TreeStorageTag>
struct tree_storage< TreeVertexProperty, TreeEdgeProperty, 
                     adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> > {
  typedef adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> alt_tag_type;
  typedef typename alt_tag_type::template alt<TreeVertexProperty,TreeEdgeProperty>::type type;
};


template <typename OutEdgeListS,
	  typename DirectedS,
	  typename VertexProperty,
	  typename AdjEdgeProperty,
	  typename AdjEdgeListS,
	  typename TreeStorageTag>
struct tree_storage_traits< adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> > :
  tree_storage_traits< TreeStorageTag > { };





template <typename AdjListOnTreeType>
class alt_tree_view {
  private:
    ReaK::shared_ptr< AdjListOnTreeType > m_alt;
    
  public:
    typedef alt_tree_view< AdjListOnTreeType > self;
    
    typedef typename AdjListOnTreeType::tree_type tree_type;
    
    // Graph traits:
    typedef typename boost::graph_traits< tree_type >::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits< tree_type >::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits< tree_type >::directed_category directed_category;
    typedef typename boost::graph_traits< tree_type >::edge_parallel_category edge_parallel_category;
    typedef typename boost::graph_traits< tree_type >::traversal_category traversal_category;
    
    static vertex_descriptor null_vertex() { return boost::graph_traits< tree_type >::null_vertex(); };
    
    // IncidenceGraph traits:
    typedef typename boost::graph_traits< tree_type >::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits< tree_type >::degree_size_type degree_size_type;
    
    // BidirectionalGraph traits:
    typedef typename boost::graph_traits< tree_type >::in_edge_iterator in_edge_iterator;
    
    // VertexListGraph traits:
    typedef typename boost::graph_traits< tree_type >::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits< tree_type >::vertices_size_type vertices_size_type;
    
    // EdgeListGraph traits:
    typedef typename boost::graph_traits< tree_type >::edge_iterator edge_iterator;
    typedef typename boost::graph_traits< tree_type >::edges_size_type edges_size_type;
    
    // AdjacencyGraph traits:
    typedef typename boost::graph_traits< tree_type >::adjacency_iterator adjacency_iterator;
    
    // PropertyGraph traits:
    typedef typename tree_type::edge_property_type edge_property_type;
    typedef typename tree_type::vertex_property_type vertex_property_type;
    
    typedef typename AdjListOnTreeType::tree_vertex_bundled vertex_bundled;
    typedef typename AdjListOnTreeType::tree_edge_bundled edge_bundled;
    
    friend class alt_graph_view<AdjListOnTreeType>;
    
    alt_tree_view() : m_alt(new AdjListOnTreeType()) { };
    
    explicit alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g);
    
    
    
    // Bundled Property-map functions (used by the boost::property_map< self, T Bundle::* > classes).
    
    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return m_alt->m_tree[v_i].tree_data;
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return m_alt->m_tree[v_i].tree_data;
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_alt->m_tree[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_alt->m_tree[e_i];
    };
    
    friend const vertex_bundled& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };
    
    friend void put( self& g, const vertex_descriptor& v_i, const vertex_bundled& value) {
      g[v_i] = value;
    };
    
    friend const edge_bundled& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };
    
    friend void put( self& g, const edge_descriptor& e_i, const edge_bundled& value) {
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
    
    
    const tree_type& get_tree() const { return m_alt->m_tree; };
    
    
    // MutableTreeConcept
    
    vertex_descriptor create_root() const {
      char this_function_is_disabled_for_this_type[0];
      //return create_root(m_alt->m_tree);
    };
    
    std::pair< vertex_descriptor, edge_descriptor> 
      add_child_vertex(vertex_descriptor v) const {
      char this_function_is_disabled_for_this_type[0];
      //return add_child_vertex(v,m_alt->m_tree);
    };
    
    void remove_branch(vertex_descriptor v) const {
      char this_function_is_disabled_for_this_type[0];
      //return remove_branch(v,m_alt->m_tree);
    };
    
    
    // MutablePropertyTreeConcept
    
    vertex_descriptor create_root(const vertex_property_type& vp) {
      return create_root(vp, m_alt->m_tree);
    };
    
    std::pair< vertex_descriptor, edge_descriptor> 
      add_child_vertex(vertex_descriptor v, 
		       const vertex_property_type& vp, 
		       const edge_property_type& ep = edge_property_type()) {
      return add_child_vertex(v, vp, ep, m_alt->m_tree);
    };
    
    template <typename OutputIter>
    void remove_branch(vertex_descriptor v, OutputIter it_out) {
      return remove_branch(v, it_out, m_alt->m_tree);
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    vertex_descriptor create_root(vertex_property_type&& vp) {
      return create_root(std::move(vp), m_alt->m_tree);
    };
    
    std::pair< vertex_descriptor, edge_descriptor> 
      add_child_vertex(vertex_descriptor v, vertex_property_type&& vp) {
      return add_child_vertex(v, std::move(vp), m_alt->m_tree);
    };
    
    std::pair< vertex_descriptor, edge_descriptor> 
      add_child_vertex(vertex_descriptor v, vertex_property_type&& vp, edge_property_type&& ep) {
      return add_child_vertex(v, std::move(vp), std::move(ep), m_alt->m_tree);
    };
#endif
    
    
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
class alt_graph_view {
  private:
    ReaK::shared_ptr< AdjListOnTreeType > m_alt;
    GraphMutationVisitor m_vis;
    
  public:
    typedef alt_graph_view< AdjListOnTreeType, GraphMutationVisitor > self;
    
    typedef typename AdjListOnTreeType::adj_list_type graph_type;
    
    // Graph traits:
    typedef typename boost::graph_traits< graph_type >::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits< graph_type >::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits< graph_type >::directed_category directed_category;
    typedef typename boost::graph_traits< graph_type >::edge_parallel_category edge_parallel_category;
    typedef typename boost::graph_traits< graph_type >::traversal_category traversal_category;
    
    static vertex_descriptor null_vertex() { return boost::graph_traits< graph_type >::null_vertex(); };
    
    // IncidenceGraph traits:
    typedef typename boost::graph_traits< graph_type >::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits< graph_type >::degree_size_type degree_size_type;
    
    // BidirectionalGraph traits:
    typedef typename boost::graph_traits< graph_type >::in_edge_iterator in_edge_iterator;
    
    // VertexListGraph traits:
    typedef typename boost::graph_traits< graph_type >::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits< graph_type >::vertices_size_type vertices_size_type;
    
    // EdgeListGraph traits:
    typedef typename boost::graph_traits< graph_type >::edge_iterator edge_iterator;
    typedef typename boost::graph_traits< graph_type >::edges_size_type edges_size_type;
    
    // AdjacencyGraph traits:
    typedef typename boost::graph_traits< graph_type >::adjacency_iterator adjacency_iterator;
    
    // PropertyGraph traits:
    typedef typename graph_type::edge_property_type edge_property_type;
    typedef typename graph_type::vertex_property_type vertex_property_type;
    
    typedef typename AdjListOnTreeType::adj_vertex_bundled vertex_bundled;
    typedef typename AdjListOnTreeType::adj_edge_bundled edge_bundled;
    
    friend class alt_tree_view<AdjListOnTreeType>;
    
    alt_graph_view() : m_alt(new AdjListOnTreeType()), m_vis() { };
    
    explicit alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t, GraphMutationVisitor aVis);
    
    
    
    
    // Bundled Property-map functions (used by the boost::property_map< self, T Bundle::* > classes).
    
    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return m_alt->m_tree[ m_alt->m_adj_list[v_i].tree_vertex ].user_data;
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return m_alt->m_tree[ m_alt->m_adj_list[v_i].tree_vertex ].user_data;
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_alt->m_adj_list[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_alt->m_adj_list[e_i];
    };
    
    friend const vertex_bundled& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };
    
    friend void put( self& g, const vertex_descriptor& v_i, const vertex_bundled& value) {
      g[v_i] = value;
    };
    
    friend const edge_bundled& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };
    
    friend void put( self& g, const edge_descriptor& e_i, const edge_bundled& value) {
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
    
    
    
    
    // IncidenceGraph concept
    
    std::pair<out_edge_iterator, out_edge_iterator> out_edges(vertex_descriptor v) const {
      return out_edges(v, m_alt->m_adj_list);
    };
    vertex_descriptor source(edge_descriptor e) const {
      return source(e, m_alt->m_adj_list);
    };
    vertex_descriptor target(edge_descriptor e) const {
      return target(e, m_alt->m_adj_list);
    };
    degree_size_type out_degree(vertex_descriptor v) const {
      return out_degree(v, m_alt->m_adj_list);
    };
    
    // BidirectionalGraph concept
    
    std::pair<in_edge_iterator, in_edge_iterator> in_edges(vertex_descriptor v) const {
      return in_edges(v, m_alt->m_adj_list);
    };
    degree_size_type in_degree(vertex_descriptor v) const {
      return in_degree(v, m_alt->m_adj_list);
    };
    degree_size_type degree(edge_descriptor e) const {
      return degree(e, m_alt->m_adj_list);
    };
    
    // AdjacencyGraph concept
    
    std::pair<adjacency_iterator, adjacency_iterator> adjacent_vertices(vertex_descriptor u) const {
      return adjacent_vertices(u, m_alt->m_adj_list);
    };
    
    // VertexListGraph concept
    
    std::pair< vertex_iterator, vertex_iterator > vertices() const {
      return vertices(m_alt->m_adj_list);
    };
    vertices_size_type num_vertices() const {
      return num_vertices(m_alt->m_adj_list);
    };
    
    // EdgeListGraph concept
    
    std::pair< edge_iterator, edge_iterator > edges() const {
      return edges(m_alt->m_adj_list);
    };
    edges_size_type num_edges() const {
      return num_edges(m_alt->m_adj_list);
    };
    
    // AdjacencyMatrix concept
    
    std::pair<edge_descriptor,bool> edge(vertex_descriptor u, vertex_descriptor v) const {
      return edge(u,v,m_alt->m_adj_list);
    };
    
    // MutableGraph concept
    
    vertex_descriptor add_vertex() {
      
    };
    
    void clear_vertex(vertex_descriptor v) {
      
    };
    
    void remove_vertex(vertex_descriptor v) {
      
    };
    
    std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v) {
      
    };
    
    void remove_edge(vertex_descriptor u, vertex_descriptor v) {
      
    };
    
    void remove_edge(edge_descriptor e) {
      
    };
    
    void remove_edge(edge_iterator e_iter) {
      
    };

    // MutablePropertyGraph concept
    
    vertex_descriptor add_vertex(const vertex_property_type& vp) {
      
    };
    
    void remove_vertex(vertex_descriptor v, vertex_property_type& vp) {
      
    };
    
    std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v, const edge_property_type& ep) {
      
    };
    
    void remove_edge(vertex_descriptor u, vertex_descriptor v, edge_property_type& ep) {
      
    };
    
    void remove_edge(edge_descriptor e, edge_property_type& ep) {
      
    };
    
    void remove_edge(edge_iterator e_iter, edge_property_type& ep) {
      
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    vertex_descriptor add_vertex(vertex_property_type&& vp) {
      
    };
    
    std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v, edge_property_type&& ep) {
      
    };
#endif
    
};

template <typename AdjListOnTreeType>
alt_tree_view<AdjListOnTreeType>::alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g) : m_alt(g.m_alt) { };

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t, GraphMutationVisitor aVis ) : m_alt(t.m_alt), m_vis(aVis) { };



template <typename AdjListOnTreeType, typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> make_graph_view(const alt_tree_view< AdjListOnTreeType >& t, GraphMutationVisitor aVis) {
  return alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>(t,aVis);
};








/******************************************************************************************
 * *************************************************************************************
 *                             Tree View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::out_edge_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::out_edge_iterator > 
  out_edges(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
	    const alt_tree_view<AdjListOnTreeType>& g) {
  return out_edges(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor 
  source(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e, 
	 const alt_tree_view<AdjListOnTreeType>& g) {
  return source(e, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor 
  target(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e, 
	 const alt_tree_view<AdjListOnTreeType>& g) {
  return target(e, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type 
  out_degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
	     const alt_tree_view<AdjListOnTreeType>& g) {
  return out_degree(v, g.get_tree());
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::in_edge_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::in_edge_iterator > 
  in_edges(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
	   const alt_tree_view<AdjListOnTreeType>& g) {
  return in_edges(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type 
  in_degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
	    const alt_tree_view<AdjListOnTreeType>& g) {
  return in_degree(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type 
  degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e, 
	 const alt_tree_view<AdjListOnTreeType>& g) {
  return degree(e, g.get_tree());
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::adjacency_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::adjacency_iterator > 
  adjacent_vertices(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor u, 
		    const alt_tree_view<AdjListOnTreeType>& g) {
  return adjacent_vertices(u, g.get_tree());
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator > 
  vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return vertices(g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertices_size_type 
  num_vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_vertices(g.get_tree());
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_iterator > 
  edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return edges(g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edges_size_type 
  num_edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_edges(g.get_tree());
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor, bool> 
  edge(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor u,
       typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
       const alt_tree_view<AdjListOnTreeType>& g) {
  return edge(u,v,g.get_tree());
};

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  get_root_vertex(const alt_tree_view<AdjListOnTreeType>& g) {
  return get_root_vertex(g.get_tree());
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator > 
  child_vertices(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
		 const alt_tree_view<AdjListOnTreeType>& g) {
  return child_vertices(v,g.get_tree());
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  create_root(const alt_tree_view<AdjListOnTreeType>& g) {
  return g.create_root();
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor> 
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
		   const alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex(v);
};

template <typename AdjListOnTreeType>
void remove_branch(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
		   const alt_tree_view<AdjListOnTreeType>& g) {
  return g.remove_branch(v);
};


/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/

    
template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor 
  create_root(const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp,
	      alt_tree_view<AdjListOnTreeType>& g) const {
  return g.create_root(vp);
};
    
template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor> 
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
		   const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp, 
		   alt_tree_view<AdjListOnTreeType>& g) const {
  return g.add_child_vertex(v, vp);
};
    
template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor> 
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
		   const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp, 
		   const typename alt_tree_view<AdjListOnTreeType>::edge_property_type& ep,
		   alt_tree_view<AdjListOnTreeType>& g) const {
  return g.add_child_vertex(v, vp, ep);
};

template <typename AdjListOnTreeType, typename OutputIter>
void remove_branch(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
		   OutputIter it_out,
		   alt_tree_view<AdjListOnTreeType>& g) const {
  g.remove_branch(v,it_out);
};
    
#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor 
  create_root(typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp,
	      alt_tree_view<AdjListOnTreeType>& g) const {
  return g.create_root(std::move(vp));
};
    
template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor> 
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
		   typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp,
		   alt_tree_view<AdjListOnTreeType>& g) const {
  return g.add_child_vertex(v, std::move(vp));
};
    
template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor, 
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor> 
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v, 
		   typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp, 
		   typename alt_tree_view<AdjListOnTreeType>::edge_property_type&& ep,
		   alt_tree_view<AdjListOnTreeType>& g) const {
  return g.add_child_vertex(v, std::move(vp), std::move(ep));
};
#endif

/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
bool is_vertex_valid(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor u,
		     const alt_tree_view<AdjListOnTreeType>& g) {
  return is_vertex_valid(u,g.get_tree());
};

template <typename AdjListOnTreeType>
bool is_edge_valid(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e,
		   const alt_tree_view<AdjListOnTreeType>& g) {
  return is_edge_valid(e,g.get_tree());
};







/******************************************************************************************
 * *************************************************************************************
 *                             Graph View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::out_edge_iterator, 
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::out_edge_iterator > 
  out_edges(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	    const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.out_edges(v);
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor 
  source(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor e, 
	 const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.source(e);
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor 
  target(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor e, 
	 const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.target(e);
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::degree_size_type 
  out_degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	     const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.out_degree(v);
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::in_edge_iterator, 
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::in_edge_iterator > 
  in_edges(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	   const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.in_edges(v);
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::degree_size_type 
  in_degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	    const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.in_degree(v);
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::degree_size_type 
  degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor e, 
	 const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.degree(e);
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::adjacency_iterator, 
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::adjacency_iterator > 
  adjacent_vertices(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
		    const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.adjacent_vertices(u);
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_iterator, 
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_iterator > 
  vertices(const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.vertices();
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertices_size_type 
  num_vertices(const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.num_vertices();
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_iterator, 
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_iterator > 
  edges(const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.edges();
};

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edges_size_type 
  num_edges(const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.num_edges();
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor, bool> 
  edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u,
       typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v,
       const alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.edge(u,v);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor 
  add_vertex(alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_vertex();
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void clear_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v,
                  alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.clear_vertex(v);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v,
                   alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_vertex(v);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor, bool> 
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
	   typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v,
	   alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_edge(u,v);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
		 typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v,
		 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(u,v);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor e,
                 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(e);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_iterator e_iter,
                 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(e_iter);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor 
  add_vertex(const typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::vertex_property_type& vp,
             alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_vertex(vp);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
		   typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::vertex_property_type& vp,
		   alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_vertex(v,vp);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor, bool> 
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
	   typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	   const typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::edge_property_type& ep,
	   alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_edge(u,v,ep);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
		 typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
		 typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::edge_property_type& ep,
		 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(u,v,ep);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor e, 
		 typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::edge_property_type& ep,
		 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(e,ep);
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_iterator e_iter, 
		 typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::edge_property_type& ep,
		 alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  g.remove_edge(e_iter,ep);
};
    
#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor 
  add_vertex(typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::vertex_property_type&& vp,
             alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_vertex(std::move(vp));
};
    
template <typename AdjListOnTreeType, typename GraphMutationVisitor>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::edge_descriptor, bool> 
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor u, 
	   typename boost::graph_traits< alt_graph_view<AdjListOnTreeType,GraphMutationVisitor> >::vertex_descriptor v, 
	   typename alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>::edge_property_type&& ep,
	   alt_graph_view<AdjListOnTreeType,GraphMutationVisitor>& g) {
  return g.add_edge(u,v,std::move(ep));
};
#endif




};

};



#endif


















