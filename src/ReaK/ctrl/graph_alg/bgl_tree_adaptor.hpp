/**
 * \file bgl_tree_adaptor.hpp
 * 
 * This library provides function templates to adapt a Mutable-Graph from the Boost Graph Library 
 * such that it has the mutable interface of a tree structure.
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

#ifndef REAK_BGL_TREE_ADAPTOR_HPP
#define REAK_BGL_TREE_ADAPTOR_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <vector>


namespace ReaK {

namespace graph {

struct bgl_tree_storage { };

template <typename VertexProperty, 
          typename EdgeProperty = boost::no_property, 
	  typename TreeStorage = bgl_tree_storage>
struct tree_storage {
  typedef boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS,
                                 VertexProperty,
		                 EdgeProperty,
		                 boost::vecS> type;
};

template <typename TreeStorage = bgl_tree_storage>
struct tree_storage_traits :
  boost::adjacency_list_traits< boost::vecS, boost::listS, boost::bidirectionalS, boost::vecS > { };

template <typename VertexProperty, typename EdgeProperty>
struct tree_traits< boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS,
                                           VertexProperty, EdgeProperty, boost::vecS> > {
  class child_vertex_iterator {
    public:
      typedef boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS,
                                   VertexProperty, EdgeProperty, boost::vecS> graph_type;
      typedef typename boost::graph_traits< graph_type >::out_edge_iterator out_edge_iter;
    private:
      out_edge_iter ei;
      graph_type const * g;
    public:
      
      child_vertex_iterator(const out_edge_iter& aEi, graph_type const * aG) : ei(aEi), g(aG) { };
      
      typedef child_vertex_iterator self;
      
      typedef std::ptrdiff_t difference_type;
      typedef typename boost::graph_traits< graph_type >::vertex_descriptor value_type;
      typedef value_type* pointer;
      typedef value_type& reference;
      typedef typename std::iterator_traits<out_edge_iter>::iterator_category iterator_category;
      
      friend bool operator==( const self& lhs, const self& rhs) { return lhs.ei == rhs.ei; };
      friend bool operator!=( const self& lhs, const self& rhs) { return lhs.ei != rhs.ei; };
      friend bool operator >( const self& lhs, const self& rhs) { return lhs.ei > rhs.ei; };
      friend bool operator >=(const self& lhs, const self& rhs) { return lhs.ei >= rhs.ei; };
      friend bool operator <( const self& lhs, const self& rhs) { return lhs.ei < rhs.ei; };
      friend bool operator <=(const self& lhs, const self& rhs) { return lhs.ei <= rhs.ei; };
      
      self& operator++() { ++ei; return *this; };
      self operator++(int) { self result(*this); ++ei; return result; };
      self& operator--() { --ei; return *this; };
      self operator--(int) { self result(*this); --ei; return result; };
      
      friend self operator+(const self& lhs, difference_type i) {
	return self(lhs.ei + i, g);
      };
      friend self operator+(difference_type i, const self& rhs) {
	return self(rhs.ei + i, g);
      };
      friend self operator-(const self& lhs, difference_type i) {
	return self(lhs.ei - i, g);
      };
      friend difference_type operator-(const self& lhs, const self& rhs) {
        return difference_type(lhs.ei - rhs.ei);
      };
      
      self& operator +=(difference_type i) { base.ei += i; return *this; };
      self& operator -=(difference_type i) { base.ei -= i; return *this; };
      
      value_type operator[](difference_type i) const { return target(*(ei + i), *g); };
      reference operator*() { return target(*ei, *g); };
      pointer operator->() { return target(*ei, *g); };
      
  };
};



};

};


namespace boost {


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

template <typename Graph >
inline
typename graph_traits<Graph>::vertex_descriptor
  get_root_vertex( const Graph& g) {
  return *(vertices(g).first);
};

template <typename Graph >
inline
std::pair<typename ReaK::graph::tree_traits<Graph>::child_vertex_iterator, 
          typename ReaK::graph::tree_traits<Graph>::child_vertex_iterator>
  child_vertices(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g) {
  typedef typename ReaK::graph::tree_traits<Graph>::child_vertex_iterator CVIter;
  typename graph_traits<Graph>::out_edge_iterator ei,ei_end;
  tie(ei,ei_end) = out_edges(v,g);
  return std::make_pair(CVIter(ei,&g),CVIter(ei_end,&g));
};
  
  
/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

template <typename Graph >
inline
typename graph_traits<Graph>::vertex_descriptor
  create_root( Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(g);
};

template <typename Graph >
inline
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    Graph& g) {
  BOOST_CONCEPT_ASSERT((MutableGraphConcept<Graph>));
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v,e);
};

template <typename Graph >
inline
void remove_branch( const typename graph_traits<Graph>::vertex_descriptor& u,
                    Graph& g) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::vector<Vertex> v_list; v_list.reserve(out_degree(u, g));
  EdgeIter ei, ei_end;
  for( tie(ei, ei_end) = out_edges(u,g); ei != ei_end; ++ei)
    v_list.push_back(target(*ei,g));
  for( typename std::vector<Vertex>::iterator it = v_list.begin(); it != v_list.end(); ++it)
    remove_branch(*it, g);
  clear_vertex(u, g);
  remove_vertex(u, g);
};




/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename graph_traits<Graph>::vertex_descriptor create_root(const typename Graph::vertex_property_type& vp, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = create_root(g);
  g[v] = vp;
  return v;
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
typename graph_traits<Graph>::vertex_descriptor
  create_root( typename Graph::vertex_property_type&& vp, 
	       Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = create_root(g);
  g[v] = std::move(vp);
  return v;
};
#endif

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& v,
		    const typename Graph::vertex_property_type& vp, Graph& g) {
  std::pair< typename graph_traits<Graph>::vertex_descriptor, 
             typename graph_traits<Graph>::edge_descriptor > result = add_child_vertex(v, g);
  g[result.first] = vp;
  return result;
};

template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& v,
		    const typename Graph::vertex_property_type& vp,
		    const typename Graph::edge_property_type& ep, Graph& g) {
  std::pair< typename graph_traits<Graph>::vertex_descriptor, 
             typename graph_traits<Graph>::edge_descriptor > result = add_child_vertex(v, g);
  g[result.first] = vp;
  g[result.second] = ep;
  return result;
};

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename VertexProperties,
          std::size_t Arity,
          typename EdgeProperties >
inline
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& v,
		    typename Graph::vertex_property_type&& vp, Graph& g) {
  std::pair< typename graph_traits<Graph>::vertex_descriptor, 
             typename graph_traits<Graph>::edge_descriptor > result = add_child_vertex(v, g);
  g[result.first] = std::move(vp);
  return result;
};

template <typename Graph>
inline
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& v,
		    typename Graph::vertex_property_type&& vp,
		    typename Graph::edge_property_type&& ep, Graph& g) {
  std::pair< typename graph_traits<Graph>::vertex_descriptor, 
             typename graph_traits<Graph>::edge_descriptor > result = add_child_vertex(v, g);
  g[result.first] = std::move(vp);
  g[result.second] = std::move(ep);
  return result;
};
#endif

template <typename Graph,
          typename OutputIter>
inline
void remove_branch( const typename graph_traits<Graph>::vertex_descriptor& u,
		    OutputIter it_out, Graph& g) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::vector<Vertex> v_list; v_list.reserve(out_degree(u, g));
  EdgeIter ei, ei_end;
  for( tie(ei, ei_end) = out_edges(u,g); ei != ei_end; ++ei)
    v_list.push_back(target(*ei,g));
  *(it_out++) = g[u];
  for( typename std::vector<Vertex>::iterator it = v_list.begin(); it != v_list.end(); ++it)
    remove_branch(*it, g);
  clear_vertex(u, g);
  remove_vertex(u, g);
};




/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/


template <typename Graph>
inline
bool is_vertex_valid( typename graph_traits<Graph>::vertex_descriptor, const Graph&) {
  return true;
};

template <typename Graph>
inline
bool is_edge_valid( typename graph_traits<Graph>::edge_descriptor, const Graph&) {
  return true;
};



};


#endif


















