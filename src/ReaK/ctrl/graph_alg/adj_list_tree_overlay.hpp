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
    friend class alt_graph_view<AdjListOnTreeType>;
    
    alt_tree_view() : m_alt(new AdjListOnTreeType()) { };
    
    explicit alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g);
    
    
};

template <typename AdjListOnTreeType>
class alt_graph_view {
  private:
    ReaK::shared_ptr< AdjListOnTreeType > m_alt;
    
  public:
    friend class alt_tree_view<AdjListOnTreeType>;
    
    alt_graph_view() : m_alt(new AdjListOnTreeType()) { };
    
    explicit alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t);
    
    
    
};

template <typename AdjListOnTreeType>
alt_tree_view<AdjListOnTreeType>::alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g) : m_alt(g.m_alt) { };

template <typename AdjListOnTreeType>
alt_graph_view<AdjListOnTreeType>::alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t) : m_alt(t.m_alt) { };












/******************************************************************************************
 * *************************************************************************************
 *                             Tree View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

std::pair<o_e_iter, o_e_iter> out_edges(v, g);
v_desc source(e, g);
v_desc target(e, g);
d_size out_degree(v, g);

/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

std::pair<i_e_iter, i_e_iter> in_edges(v, g);
d_size in_degree(v, g);
d_size degree(e, g);

/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

std::pair<adj_iter, adj_iter> adjacent_vertices(u, g);

/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

std::pair< v_iter, v_iter > vertices(g);
v_size num_vertices(g);

/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

std::pair< e_iter, e_iter > edges(g);
e_size num_edges(g);

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

e_desc edge(u,v,g);

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

v_desc get_root_vertex(g);
v_desc create_root(g);
std::pair<v_desc, e_desc> add_child_vertex(v,g);
void remove_branch(v,g);
std::pair<v_iter, v_iter> child_vertices(v,g);

/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/

bool is_vertex_valid(u,g);
bool is_edge_valid(e,g);


/*******************************************************************************************
 *                  PropertyGraph concept
 ******************************************************************************************/

p_map get(property,g);
value get(property,g,u);
value get(property,g,e);
void put(property,g,u,value);
void put(property,g,e,value);









/******************************************************************************************
 * *************************************************************************************
 *                             Graph View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

std::pair<o_e_iter, o_e_iter> out_edges(v, g);
v_desc source(e, g);
v_desc target(e, g);
d_size out_degree(v, g);

/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

std::pair<i_e_iter, i_e_iter> in_edges(v, g);
d_size in_degree(v, g);
d_size degree(e, g);

/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

std::pair<adj_iter, adj_iter> adjacent_vertices(u, g);

/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

std::pair< v_iter, v_iter > vertices(g);
v_size num_vertices(g);

/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

std::pair< e_iter, e_iter > edges(g);
e_size num_edges(g);

/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

std::pair<e_desc,bool> edge(u,v,g);

/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

v_desc add_vertex(g);
void clear_vertex(v, g);
void remove_vertex(v, g);
std::pair<e_desc, bool> add_edge(u, v, g);
void remove_edge(u, v, g);
void remove_edge(e, g);
void remove_edge(e_iter, g);

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

v_desc add_vertex(vp, g)
std::pair<e_desc, bool> add_edge(u, v, ep, g)

/*******************************************************************************************
 *                  PropertyGraph concept
 ******************************************************************************************/

p_map get(property,g);
value get(property,g,u);
value get(property,g,e);
void put(property,g,u,value);
void put(property,g,e,value);






};

};



namespace boost {
  

template <typename AdjListOnTreeType>
struct graph_traits< alt_graph_view<AdjListOnTreeType> > {
  // Graph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::edge_descriptor edge_descriptor;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::directed_category directed_category;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::edge_parallel_category edge_parallel_category;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::traversal_category traversal_category;
  
  static vertex_descriptor null_vertex() { return graph_traits< typename AdjListOnTreeType::adj_list_type >::null_vertex(); };
  
  // IncidenceGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::degree_size_type degree_size_type;
  
  // BidirectionalGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::in_edge_iterator in_edge_iterator;
  
  // VertexListGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::vertex_iterator vertex_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::vertices_size_type vertices_size_type;
  
  // EdgeListGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::edge_iterator edge_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::edges_size_type edges_size_type;
  
  // AdjacencyGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::adj_list_type >::adjacency_iterator adjacency_iterator;
  
};


// TODO NOTE HERE TODO NOTE
template <typename AdjListOnTreeType, typename TODO_CREATE_THIS>
struct property_map< alt_graph_view<AdjListOnTreeType>, TODO_CREATE_THIS > {
  typedef ... type;
  typedef ... const_type;
};



template <typename AdjListOnTreeType>
struct graph_traits< alt_tree_view<AdjListOnTreeType> > {
  // Graph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::edge_descriptor edge_descriptor;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::directed_category directed_category;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::edge_parallel_category edge_parallel_category;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::traversal_category traversal_category;
  
  static vertex_descriptor null_vertex() { return graph_traits< typename AdjListOnTreeType::tree_type >::null_vertex(); };
  
  // IncidenceGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::degree_size_type degree_size_type;
  
  // BidirectionalGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::in_edge_iterator in_edge_iterator;
  
  // VertexListGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::vertex_iterator vertex_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::vertices_size_type vertices_size_type;
  
  // EdgeListGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::edge_iterator edge_iterator;
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::edges_size_type edges_size_type;
  
  // AdjacencyGraph traits:
  typedef typename graph_traits< typename AdjListOnTreeType::tree_type >::adjacency_iterator adjacency_iterator;
  
};

  

// TODO NOTE HERE TODO NOTE
template <typename AdjListOnTreeType, typename TODO_CREATE_THIS>
struct property_map< alt_graph_view<AdjListOnTreeType>, TODO_CREATE_THIS > {
  typedef ... type;
  typedef ... const_type;
};
  
  
};




#endif


















