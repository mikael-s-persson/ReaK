/**
 * \file pooled_adjacency_list.hpp
 *
 * This library provides a class that implements a pooled-memory adjacency-list, wrapping the BGL's 
 * adjacency-list implementation. The point of a pooled-memory adjacency-list is simply to be able 
 * to retain most of the benefits of the vector storage of the vertex-list without invalidating 
 * vertex-descriptors and vertex-iterators when removing nodes from the graph. In other words, 
 * instead of removing nodes from the graph, it simply puts them in a "dead" state and revives 
 * them later when nodes are added to the graph. So, as long as the rate of insertion is greater 
 * than or equal to the rate of removal of nodes, there isn't much waste of memory and vertex 
 * descriptors or iterators are never invalidated.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2013
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

#ifndef REAK_POOLED_ADJACENCY_LIST_HPP
#define REAK_POOLED_ADJACENCY_LIST_HPP

#include "base/defs.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/variant.hpp>

#include "bgl_more_property_tags.hpp"
#include "bgl_more_property_maps.hpp"
#include "bgl_raw_property_graph.hpp"

#include <vector>
#include <queue>
#include <algorithm>
#include <functional>


namespace boost {



template <typename DirectedS = directedS,
          typename VertexProperty = no_property,
          typename EdgeProperty = no_property,
          typename GraphProperty = no_property,
          typename EdgeListS = listS>
class pooled_adjacency_list {
  public:
    typedef pooled_adjacency_list< DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS > self;
    
    struct hole_descriptor {
      std::size_t value;
      explicit hole_descriptor(std::size_t aValue = 0) : value(aValue) { };
    };
    
    // PropertyGraph traits:
    typedef EdgeProperty edge_property_type;
    typedef VertexProperty vertex_property_type;
    typedef variant<VertexProperty, hole_descriptor> raw_vertex_property_type;
    
    typedef edge_property_type edge_bundled;
    typedef VertexProperty vertex_bundled;
    
    typedef adjacency_list< vecS, vecS, DirectedS, raw_vertex_property_type, edge_property_type, GraphProperty, EdgeListS > graph_type;
    
    // Graph traits:
    typedef typename graph_traits< graph_type >::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits< graph_type >::edge_descriptor edge_descriptor;
    typedef typename graph_traits< graph_type >::directed_category directed_category;
    typedef typename graph_traits< graph_type >::edge_parallel_category edge_parallel_category;
    typedef typename graph_traits< graph_type >::traversal_category traversal_category;
    
    static vertex_descriptor null_vertex() { return graph_traits< graph_type >::null_vertex(); };
    
    // IncidenceGraph traits:
    typedef typename graph_traits< graph_type >::out_edge_iterator out_edge_iterator;
    typedef typename graph_traits< graph_type >::degree_size_type degree_size_type;
    
    // BidirectionalGraph traits:
    typedef typename graph_traits< graph_type >::in_edge_iterator in_edge_iterator;
    
    // VertexListGraph traits:
    struct vertex_iterator {
      typedef std::ptrdiff_t difference_type;
      typedef vertex_descriptor value_type;
      typedef vertex_descriptor* pointer;
      typedef vertex_descriptor& reference;
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef typename graph_traits< graph_type >::vertex_iterator base_iterator;
      
      base_iterator base;
      const graph_type* p_graph;
      vertex_iterator(const graph_type* aPGraph = NULL,
                      const base_iterator& aBase = base_iterator()) : 
                      base(aBase), p_graph(aPGraph) { };
      
      friend bool operator==(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base == rhs.base; };
      friend bool operator!=(const vertex_iterator& lhs, const vertex_iterator& rhs) { return lhs.base != rhs.base; };
      
      vertex_iterator& operator++() { 
        ++base;
        while( (base != vertices(*p_graph).second) && ((*p_graph)[*base].which() == 1) )
          ++base;
        return *this;
      };
      vertex_iterator operator++(int) { vertex_iterator result(*this); ++(*this); return result; };
      vertex_iterator& operator--() { 
        --base;
        while( (base != vertices(*p_graph).first) && ((*p_graph)[*base].which() == 1) )
          --base;
        return *this;
      };
      vertex_iterator operator--(int) { vertex_iterator result(*this); --(*this); return result; };
      
      value_type operator*() const { return *base; };
//       pointer operator->() { return &(*base); };
    };
    
    
    typedef typename graph_traits< graph_type >::vertices_size_type vertices_size_type;
    
    // EdgeListGraph traits:
    typedef typename graph_traits< graph_type >::edge_iterator edge_iterator;
    typedef typename graph_traits< graph_type >::edges_size_type edges_size_type;
    
    // AdjacencyGraph traits:
    typedef typename graph_traits< graph_type >::adjacency_iterator adjacency_iterator;
    
    typedef typename graph_type::graph_bundled graph_bundled;
    
  private:
    
    graph_type m_graph;
    
    hole_descriptor m_first_hole;
    vertices_size_type m_num_vertices;
    
  public:
    
    pooled_adjacency_list() : m_graph(), m_first_hole(null_vertex()), m_num_vertices(0) { };
    
    
    // Bundled Property-map functions (used by the property_map< self, T Bundle::* > classes).
    
    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return get<vertex_bundled>(m_graph[v_i]);
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return get<vertex_bundled>(m_graph[v_i]);
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_graph[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_graph[e_i];
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
    typename property_map< self, T Bundle::* >::type
    get( T Bundle::* p, self& g) {
      return typename property_map< self, T Bundle::* >::type(&g, p);
    };

    template <typename T, typename Bundle>
    friend
    typename property_map< self, T Bundle::* >::const_type
    get( T Bundle::* p, const self& g) {
      return typename property_map< self, T Bundle::* >::const_type(&g, p);
    };




    // IncidenceGraph concept

    std::pair<out_edge_iterator, out_edge_iterator> out_edges_impl(vertex_descriptor v) const {
      return out_edges(v, m_graph);
    };
    vertex_descriptor source_impl(edge_descriptor e) const {
      return source(e, m_graph);
    };
    vertex_descriptor target_impl(edge_descriptor e) const {
      return target(e, m_graph);
    };
    degree_size_type out_degree_impl(vertex_descriptor v) const {
      return out_degree(v, m_graph);
    };

    // BidirectionalGraph concept

    std::pair<in_edge_iterator, in_edge_iterator> in_edges_impl(vertex_descriptor v) const {
      return in_edges(v, m_graph);
    };
    degree_size_type in_degree_impl(vertex_descriptor v) const {
      return in_degree(v, m_graph);
    };
    degree_size_type degree_impl(edge_descriptor e) const {
      return degree(e, m_graph);
    };

    // AdjacencyGraph concept

    std::pair<adjacency_iterator, adjacency_iterator> adjacent_vertices_impl(vertex_descriptor u) const {
      return adjacent_vertices(u, m_graph);
    };

    // VertexListGraph concept

    std::pair< vertex_iterator, vertex_iterator > vertices_impl() const {
      typedef typename graph_traits< graph_type >::vertex_iterator VIter;
      std::pair< VIter, VIter > tmp = vertices(m_graph);
      return std::pair< vertex_iterator, vertex_iterator >(vertex_iterator(&m_graph,tmp.first), vertex_iterator(&m_graph,tmp.second));
    };
    vertices_size_type num_vertices_impl() const {
      return m_num_vertices;
    };

    // EdgeListGraph concept

    std::pair< edge_iterator, edge_iterator > edges_impl() const {
      return edges(m_graph);
    };
    edges_size_type num_edges_impl() const {
      return num_edges(m_graph);
    };

    // AdjacencyMatrix concept

    std::pair<edge_descriptor,bool> edge_impl(vertex_descriptor u, vertex_descriptor v) const {
      return edge(u,v,m_graph);
    };

    // MutableGraph concept
    
    vertex_descriptor add_vertex_impl() {
      vertex_descriptor v;
      if( m_first_hole.value == null_vertex() ) {
        // add a new node in the graph (this is safe, it will not invalidate vertex descriptors).
        v = add_vertex(m_graph);
      } else {
        // resurrect a node from the graveyard.
        v = m_first_hole.value;
        m_first_hole = get< hole_descriptor >(m_graph[v]);
      };
      m_graph[v] = vertex_property_type();
      ++m_num_vertices;
      return v;
    };

    void clear_vertex_impl(vertex_descriptor v) {
      clear_vertex(v, m_graph);
    };

    void remove_vertex_impl(vertex_descriptor v) {
      clear_vertex(v, m_graph);
      m_graph[v] = m_first_hole;
      m_first_hole.value = v;
      --m_num_vertices;
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v) {
      return add_edge(u, v, m_graph);
    };

    void remove_edge_impl(vertex_descriptor u, vertex_descriptor v) {
      remove_edge(u, v, m_graph);
    };

    void remove_edge_impl(edge_descriptor e) {
      remove_edge(e, m_graph);
    };

    void remove_edge_impl(edge_iterator e_iter) {
      remove_edge(e_iter, m_graph);
    };

    // MutablePropertyGraph concept

    vertex_descriptor add_vertex_impl(const vertex_property_type& vp) {
      vertex_descriptor v;
      if( m_first_hole.value == null_vertex() ) {
        // add a new node in the graph (this is safe, it will not invalidate vertex descriptors).
        v = add_vertex(m_graph);
      } else {
        // resurrect a node from the graveyard.
        v = m_first_hole.value;
        m_first_hole = get< hole_descriptor >(m_graph[v]);
      };
      m_graph[v] = vp;
      ++m_num_vertices;
      return v;
    };

    void remove_vertex_impl(vertex_descriptor v, vertex_property_type& vp) {
      vp = m_graph[v];
      remove_vertex_impl(v);
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v, const edge_property_type& ep) {
      return add_edge(u, v, ep, m_graph);
    };

    void remove_edge_impl(vertex_descriptor u, vertex_descriptor v, edge_property_type& ep) {
      remove_edge(u, v, ep, m_graph);
    };

    void remove_edge_impl(edge_descriptor e, edge_property_type& ep) {
      remove_edge(e, ep, m_graph);
    };

    void remove_edge_impl(edge_iterator e_iter, edge_property_type& ep) {
      remove_edge(e_iter, ep, m_graph);
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_vertex_impl(vertex_property_type&& vp) {
      vertex_descriptor v;
      if( m_first_hole.value == null_vertex() ) {
        // add a new node in the graph (this is safe, it will not invalidate vertex descriptors).
        v = add_vertex(m_graph);
      } else {
        // resurrect a node from the graveyard.
        v = m_first_hole.value;
        m_first_hole = get< hole_descriptor >(m_graph[v]);
      };
      m_graph[v] = std::move(vp);
      ++m_num_vertices;
      return v;
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v, edge_property_type&& ep) {
      return add_edge(u, v, std::move(ep), m_graph);
    };
#endif

    // NonCompactGraphConcept

    bool is_vertex_valid_impl(vertex_descriptor u) const {
      if( ( u != null_vertex() ) && ( m_graph[u].which() == 0 ) )
        return true;
      else
        return false;
    };

    bool is_edge_valid_impl(edge_descriptor e) const {
      return true;
    };


#if 0
    // TODO This will require a bit more work since we cannot establish that a given vertex is invalid or not.
    
    // this could go into NonCompactGraphConcept
    void repack_graph() {
      if(m_available_nodes.empty())
        return;
      // here, the idea is to essentially empty out the graveyard,
      if( is_same< vertices_size_type, vertex_descriptor >::type::value ) {
        // if vertex descriptors are into a random-access container, then a simple remove loop won't cut it.
        // the idea here is to achieve a behavior similar to the std::remove function, but by taking
        // elements from the back of the vertex list and swapping them with the ones in the graveyard.
        using std::swap;
        vertex_descriptor v_end = reinterpret_cast<vertex_descriptor>(num_vertices(m_graph));
        while(!m_available_nodes.empty()) {
          vertex_descriptor v = m_available_nodes.top();
          while( ( v_end > v ) && ( !is_vertex_valid(--v_end) ) )
            /* nothing */;
          if(v > v_end) {
            ++v_end;
            break;
          };
          swap(m_graph[v], m_graph[v_end]);
          out_edge_iterator ei, ei_end;
          for(boost::tie(ei,ei_end) = out_edges(v_end, m_graph); ei != ei_end; ++ei) {
            std::pair<edge_descriptor, bool> e = add_edge(v, target(*ei, m_graph), m_graph);
            swap(m_graph[e],m_graph[*ei]);
          };
          if( boost::is_same< directed_category, boost::directed_tag >::type::value ) {
            in_edge_iterator in_ei, in_ei_end;
            for(boost::tie(in_ei, in_ei_end) = in_edges(v_end, m_graph); in_ei != in_ei_end; ++in_ei) {
              std::pair<edge_descriptor, bool> e = add_edge(source(*in_ei, m_graph), v, m_graph);
              swap(m_graph[e], m_graph[*in_ei]);
            };
          };
          clear_vertex(v_end, m_graph);
          m_available_nodes.pop();
        };
        while(!m_available_nodes.empty())
          m_available_nodes.pop();
        // at this point v_end is the first invalid vertex remaining.
        vertex_descriptor v_rm = reinterpret_cast<vertex_descriptor>(num_vertices(m_graph));
        while(v_rm > v_end)
          remove_vertex(--v_rm, m_graph);

      } else {
        // if vertex descriptors are not random-access, they should not get invalidated by removal of other vertices.
        while(!m_available_nodes.empty()) {
          clear_vertex(m_available_nodes.top(), m_graph);
          remove_vertex(m_available_nodes.top(), m_graph);
          m_available_nodes.pop();
        };
      };
    };
#endif

};




/******************************************************************************************
 * *************************************************************************************
 *                             Graph functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::out_edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::out_edge_iterator >
  out_edges(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
            const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.out_edges_impl(v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  source(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.source_impl(e);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  target(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.target_impl(e);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  out_degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
             const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.out_degree_impl(v);
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::in_edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::in_edge_iterator >
  in_edges(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.in_edges_impl(v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  in_degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
            const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.in_degree_impl(v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.degree_impl(e);
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::adjacency_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::adjacency_iterator >
  adjacent_vertices(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                    const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.adjacent_vertices_impl(u);
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_iterator >
  vertices(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.vertices_impl();
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertices_size_type
  num_vertices(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.num_vertices_impl();
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator >
  edges(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.edges_impl();
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edges_size_type
  num_edges(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.num_edges_impl();
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
       typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
       const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.edge_impl(u,v);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_vertex_impl();
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void clear_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                  pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.clear_vertex_impl(v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                   pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_vertex_impl(v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_edge_impl(u,v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                 typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(u,v);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(e);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator e_iter,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(e_iter);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type& vp,
             pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_vertex_impl(vp);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                   typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type& vp,
                   pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_vertex_impl(v,vp);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_edge_impl(u,v,ep);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                 typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(u,v,ep);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(e,ep);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator e_iter,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  g.remove_edge_impl(e_iter,ep);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type&& vp,
             pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_vertex_impl(std::move(vp));
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type&& ep,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.add_edge_impl(u,v,std::move(ep));
};
#endif

/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
bool is_vertex_valid(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                     const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.is_vertex_valid_impl(u);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
bool is_edge_valid(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                   const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.is_edge_valid_impl(e);
};





template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
struct is_raw_property_graph< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> > : boost::mpl::true_ { };


template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, typename T, typename Bundle>
struct property_map< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


};


#endif


















