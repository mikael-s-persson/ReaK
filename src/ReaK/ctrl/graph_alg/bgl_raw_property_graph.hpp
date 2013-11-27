/**
 * \file bgl_raw_property_graph.hpp
 * 
 * This library provides traits and concepts definitions for a RawPropertyGraphConcept.
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

#ifndef REAK_BGL_RAW_PROPERTY_GRAPH_HPP
#define REAK_BGL_RAW_PROPERTY_GRAPH_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list_BC.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>
    
namespace boost {


enum vertex_raw_prop_to_bundle_t { vertex_raw_prop_to_bundle };
enum vertex_raw_property_t { vertex_raw_property };

BOOST_INSTALL_PROPERTY(vertex, raw_prop_to_bundle);
BOOST_INSTALL_PROPERTY(vertex, raw_property);

enum edge_raw_prop_to_bundle_t { edge_raw_prop_to_bundle };
enum edge_raw_property_t { edge_raw_property };

BOOST_INSTALL_PROPERTY(edge, raw_prop_to_bundle);
BOOST_INSTALL_PROPERTY(edge, raw_property);




template <typename Graph>
struct property_map<Graph, vertex_raw_property_t> { 
  typedef whole_bundle_property_map<Graph, vertex_bundle_t> type;
  typedef whole_bundle_property_map<const Graph, vertex_bundle_t> const_type;
};

template <typename Graph>
struct property_map<Graph, edge_raw_property_t> { 
  typedef whole_bundle_property_map<Graph, edge_bundle_t> type;
  typedef whole_bundle_property_map<const Graph, edge_bundle_t> const_type;
};

template <typename Graph>
struct property_map<Graph, vertex_raw_prop_to_bundle_t> {
  typedef self_property_map< typename Graph::vertex_bundled > type;
  typedef self_property_map< const typename Graph::vertex_bundled > const_type;
};

template <typename Graph>
struct property_map<Graph, edge_raw_prop_to_bundle_t> {  
  typedef self_property_map< typename Graph::edge_bundled > type;
  typedef self_property_map< const typename Graph::edge_bundled > const_type;
};




template <typename Graph>
struct RawPropertyGraphConcept {
  typename graph_traits<Graph>::vertex_descriptor v;
  typename graph_traits<Graph>::edge_descriptor e;
  Graph g;
  
  typedef typename property_map<Graph, vertex_raw_prop_to_bundle_t>::type vertex_raw_prop_to_bundle_map;
  typedef typename property_map<Graph, vertex_raw_property_t>::type vertex_raw_prop_map;
  typedef typename property_map<Graph, edge_raw_prop_to_bundle_t>::type edge_raw_prop_to_bundle_map;
  typedef typename property_map<Graph, edge_raw_property_t>::type edge_raw_prop_map;
  
  vertex_raw_prop_to_bundle_map vrp_to_bundle;
  vertex_raw_prop_map vrp_map;
  edge_raw_prop_to_bundle_map erp_to_bundle;
  edge_raw_prop_map erp_map;
  
  typedef typename property_traits< vertex_raw_prop_to_bundle_map >::key_type v_prop_type;
  typedef typename property_traits< vertex_raw_prop_to_bundle_map >::value_type v_bundle_type;
  typedef typename property_traits< vertex_raw_prop_to_bundle_map >::key_type e_prop_type;
  typedef typename property_traits< vertex_raw_prop_to_bundle_map >::value_type e_bundle_type;
  
  v_prop_type vrp;
  v_bundle_type vb;
  e_prop_type erp;
  e_bundle_type eb;
  
  BOOST_CONCEPT_ASSERT((boost::GraphConcept<Graph>));
  
  BOOST_CONCEPT_USAGE(RawPropertyGraphConcept) 
  {
    vrp_to_bundle = get(vertex_raw_prop_to_bundle, g);
    vrp_map = get(vertex_raw_property, g);
    erp_to_bundle = get(edge_raw_prop_to_bundle, g);
    erp_map = get(edge_raw_property, g);
    
    vb = get(vertex_raw_prop_to_bundle, g, vrp);
    vrp = get(vertex_raw_property, g, v);
    eb = get(edge_raw_prop_to_bundle, g, erp);
    erp = get(edge_raw_property, g, e);
    
    put(vertex_raw_prop_to_bundle, g, vrp, vb);
    put(vertex_raw_property, g, v, vrp);
    put(edge_raw_prop_to_bundle, g, erp, eb);
    put(edge_raw_property, g, e, erp);
    
    vb = get(vrp_to_bundle, vrp);
    vrp = get(vrp_map, v);
    eb = get(erp_to_bundle, erp);
    erp = get(erp_map, e);
    
    put(vrp_to_bundle, vrp, vb);
    put(vrp_map, v, vrp);
    put(erp_to_bundle, erp, eb);
    put(erp_map, e, erp);
  };
  
  
};


template <typename Graph>
struct is_raw_property_graph : boost::mpl::false_ { };

template <typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperty, typename EdgeProperty>
struct is_raw_property_graph< adjacency_list_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty> > : boost::mpl::true_ { };



template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, vertex_raw_property_t> >::type::type get(vertex_raw_property_t, Graph& g) {
  typedef typename property_map<Graph, vertex_raw_property_t>::type result_type;
  return result_type(&g);
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, vertex_raw_property_t> >::type::const_type get(vertex_raw_property_t, const Graph& g) {
  typedef typename property_map<Graph, vertex_raw_property_t>::const_type result_type;
  return result_type(&g);
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, edge_raw_property_t> >::type::type get(edge_raw_property_t, Graph& g) {
  typedef typename property_map<Graph, edge_raw_property_t>::type result_type;
  return result_type(&g);
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, edge_raw_property_t> >::type::const_type get(edge_raw_property_t, const Graph& g) {
  typedef typename property_map<Graph, edge_raw_property_t>::const_type result_type;
  return result_type(&g);
};



template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, vertex_raw_prop_to_bundle_t> >::type::type get(vertex_raw_prop_to_bundle_t, Graph& g) {
  return typename property_map<Graph, vertex_raw_prop_to_bundle_t>::type();
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, vertex_raw_prop_to_bundle_t> >::type::const_type get(vertex_raw_prop_to_bundle_t, const Graph& g) {
  return typename property_map<Graph, vertex_raw_prop_to_bundle_t>::const_type();
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, edge_raw_prop_to_bundle_t> >::type::type get(edge_raw_prop_to_bundle_t, Graph& g) {
  return typename property_map<Graph, edge_raw_prop_to_bundle_t>::type();
};

template <typename Graph>
typename enable_if< is_raw_property_graph< Graph >,
property_map<Graph, edge_raw_prop_to_bundle_t> >::type::const_type get(edge_raw_prop_to_bundle_t, const Graph& g) {
  return typename property_map<Graph, edge_raw_prop_to_bundle_t>::const_type();
};




template <typename Graph, typename Vertex>
typename enable_if< is_raw_property_graph< Graph >,
Graph >::type::vertex_bundled& get(vertex_raw_property_t, Graph& g, Vertex v) {
  return g[v];
};

template <typename Graph, typename Vertex>
const typename enable_if< is_raw_property_graph< Graph >,
Graph >::type::vertex_bundled& get(vertex_raw_property_t, const Graph& g, Vertex v) {
  return g[v];
};

template <typename Graph, typename Edge>
typename enable_if< is_raw_property_graph< Graph >,
Graph >::type::edge_bundled& get(edge_raw_property_t, Graph& g, Edge e) {
  return g[e];
};

template <typename Graph, typename Edge>
const typename enable_if< is_raw_property_graph<Graph>,
typename Graph::edge_bundled >::type& get(edge_raw_property_t, const Graph& g, Edge e) {
  return g[e];
};


template <typename Graph, typename Bundle>
Bundle& get(vertex_raw_prop_to_bundle_t, const Graph&, Bundle& b) {
  return b;
};

template <typename Graph, typename Bundle>
Bundle& get(edge_raw_prop_to_bundle_t, const Graph&, Bundle& b) {
  return b;
};



#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES

template <typename Graph, typename Vertex, typename Bundle>
void put(vertex_raw_property_t, Graph& g, Vertex v, const Bundle& value) {
  g[v] = value;
};

template <typename Graph, typename Edge, typename Bundle>
void put(edge_raw_property_t, Graph& g, Edge e, const Bundle& value) {
  g[e] = value;
};

#else

template <typename Graph, typename Vertex, typename Bundle>
void put(vertex_raw_property_t, Graph& g, Vertex v, Bundle&& value) {
  g[v] = std::forward<Bundle>(value);
};

template <typename Graph, typename Edge, typename Bundle>
void put(edge_raw_property_t, Graph& g, Edge e, Bundle&& value) {
  g[e] = std::forward<Bundle>(value);
};

#endif


#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES

template <typename Graph, typename Bundle1, typename Bundle2>
void put(vertex_raw_prop_to_bundle_t, const Graph&, Bundle1& b1, const Bundle2& b2) {
  b1 = b2;
};

template <typename Graph, typename Bundle1, typename Bundle2>
void put(edge_raw_prop_to_bundle_t, const Graph&, Bundle1& b1, const Bundle2& b2) {
  b1 = b2;
};

#else

template <typename Graph, typename Bundle1, typename Bundle2>
void put(vertex_raw_prop_to_bundle_t, const Graph&, Bundle1& b1, Bundle2&& b2) {
  b1 = std::forward<Bundle2>(b2);
};

template <typename Graph, typename Bundle1, typename Bundle2>
void put(edge_raw_prop_to_bundle_t, const Graph&, Bundle1& b1, Bundle2&& b2) {
  b1 = std::forward<Bundle2>(b2);
};

#endif


};


#endif


















