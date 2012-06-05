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

};

};


namespace boost {

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template <typename Graph >
inline
typename graph_traits<Graph>::vertex_descriptor
  get_root_vertex( const Graph& g) {
  return *(vertices(g).first);
};

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


















