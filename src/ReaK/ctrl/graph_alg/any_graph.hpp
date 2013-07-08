/**
 * \file any_graph.hpp
 *
 * This library provides a type-erasure base-class to represent any graph.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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
 *    If not, sMOTION_ee <http://www.gnu.org/licenses/>.
 */

#ifndef RK_ANY_GRAPH_HPP
#define RK_ANY_GRAPH_HPP

#include <string>

#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <boost/any.hpp>
#include <boost/range/any_range.hpp>

#include <boost/type_traits.hpp>


#include "bgl_more_property_tags.hpp"


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {


class any_graph {
  public:
    typedef boost::any any_descriptor;
    
    struct vertex_descriptor {
      any_descriptor base;
      explicit vertex_descriptor(const any_descriptor& aBase) : base(aBase) { };
    };
    
    struct edge_descriptor {
      any_descriptor base;
      explicit edge_descriptor(const any_descriptor& aBase) : base(aBase) { };
    };
    
    typedef boost::any any_bundled;
    typedef any_bundled vertex_bundled;
    typedef any_bundled edge_bundled;
    
    typedef boost::any_range<edge_descriptor, boost::forward_pass_traversal_tag, edge_descriptor, std::ptrdiff_t > edge_range;
    typedef edge_range out_edge_range;
    typedef edge_range in_edge_range;
    
    typedef edge_range::iterator edge_iterator;
    typedef edge_iterator out_edge_iterator;
    typedef edge_iterator in_edge_iterator;
    
    typedef boost::any_range<vertex_descriptor, boost::forward_pass_traversal_tag, vertex_descriptor, std::ptrdiff_t > vertex_range;
    
    typedef vertex_range::iterator vertex_iterator;
    
    typedef std::size_t degree_size_type;
    typedef std::size_t vertices_size_type;
    typedef std::size_t edges_size_type;
    
  protected:
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const any_descriptor& aElement) const = 0;
    virtual boost::any get_property_by_any(const std::string& aProperty, const any_descriptor& aElement) const = 0;
    
    virtual vertex_bundled get_vertex_bundled(const vertex_descriptor& aVertex) const = 0;
    virtual edge_bundled get_edge_bundled(const edge_descriptor& aEdge) const = 0;
    
  public:
    
    vertex_bundled operator[](const vertex_descriptor& aVertex) const { return get_vertex_bundled(aVertex); };
    edge_bundled operator[](const edge_descriptor& aEdge) const { return get_edge_bundled(aEdge); };
    
    virtual edge_range edges() const = 0;
    virtual edges_size_type num_edges() const = 0;
    
    virtual vertex_descriptor source(const edge_descriptor& aEdge) const = 0;
    virtual vertex_descriptor target(const edge_descriptor& aEdge) const = 0;
    
    virtual std::pair<edge_descriptor, bool> edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const = 0;
    
    virtual out_edge_range out_edges(const vertex_descriptor& aVertex) const = 0;
    virtual degree_size_type out_degree(const vertex_descriptor& aVertex) const = 0;
    
    virtual in_edge_range in_edges(const vertex_descriptor& aVertex) const = 0;
    virtual degree_size_type in_degree(const vertex_descriptor& aVertex) const = 0;
    virtual degree_size_type degree(const edge_descriptor& aEdge) const = 0;
    
    virtual vertex_range vertices() const = 0;
    virtual vertices_size_type num_vertices() const = 0;
    
    virtual vertex_descriptor add_vertex() const = 0;
    virtual void clear_vertex(const vertex_descriptor& aVertex) const = 0;
    virtual void remove_vertex(const vertex_descriptor& aVertex) const = 0;
    
    virtual std::pair<edge_descriptor, bool> add_edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const = 0;
    virtual void remove_edge(const edge_descriptor& aEdge) const = 0;
    
  public:
    
    
    template <typename ValueType>
    struct property_map_by_ptr : public boost::put_get_helper<ValueType&, property_map_by_ptr<ValueType> > {
      typedef property_map_by_ptr<ValueType> self;
      
      typedef ValueType value_type;
      typedef ValueType& reference;
      typedef boost::any key_type;
      typedef typename boost::mpl::if_< boost::is_const<ValueType>,
        boost::readable_property_map_tag,
        boost::lvalue_property_map_tag >::type category;
      
      const any_graph* p_parent;
      std::string prop_name;
      
      property_map_by_ptr(const any_graph* aParentPtr, const std::string& aProperty) : p_parent(aParentPtr), prop_name(aProperty) { };
      reference operator[](const key_type& k) const { 
        return *(static_cast< ValueType* >(p_parent->get_property_by_ptr(prop_name, k)));
      };
      
    };
    
    template <typename ValueType>
    struct property_map_by_any : public boost::put_get_helper<ValueType, property_map_by_any<ValueType> > {
      typedef property_map_by_any<ValueType> self;
      
      typedef ValueType value_type;
      typedef ValueType reference;
      typedef boost::any key_type;
      typedef boost::readable_property_map_tag category;
      
      const any_graph* p_parent;
      std::string prop_name;
      
      property_map_by_any(const any_graph* aParentPtr, const std::string& aProperty) : p_parent(aParentPtr), prop_name(aProperty) { };
      reference operator[](const key_type& k) const { 
        return boost::any_cast< ValueType >(p_parent->get_property_by_any(prop_name, aElement));
      };
      
    };
    
    
};


template <typename ValueType>
boost::enable_if< boost::is_reference< ValueType >,
any_graph::property_map_by_ptr< typename boost::remove_reference< ValueType >::type > >::type 
  get(const std::string& aProperty, const any_graph& aGraph) {
  return any_graph::property_map_by_ptr< typename boost::remove_reference< ValueType >::type >(&aGraph, aProperty);
};

template <typename ValueType>
boost::disable_if< boost::is_reference< ValueType >,
any_graph::property_map_by_any< ValueType > >::type get(const std::string& aProperty, const any_graph& aGraph) {
  return any_graph::property_map_by_any< ValueType >(&aGraph, aProperty);
};



/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

std::pair< any_graph::out_edge_iterator, any_graph::out_edge_iterator > out_edges(any_graph::vertex_descriptor v, const any_graph& g) {
  any_graph::out_edge_range er = g.out_edges(v);
  return std::pair< any_graph::out_edge_iterator, any_graph::out_edge_iterator >(begin(er), end(er));
};

any_graph::vertex_descriptor source(any_graph::edge_descriptor e, const any_graph& g) {
  return g.source(e);
};

any_graph::vertex_descriptor target(any_graph::edge_descriptor e, const any_graph& g) {
  return g.target(e);
};

any_graph::degree_size_type out_degree(any_graph::vertex_descriptor v, const any_graph& g) {
  return g.out_degree(v);
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

std::pair< any_graph::in_edge_iterator, any_graph::in_edge_iterator > in_edges(any_graph::vertex_descriptor v, const any_graph& g) {
  any_graph::in_edge_range er = g.in_edges(v);
  return std::pair< any_graph::in_edge_iterator, any_graph::in_edge_iterator >(begin(er), end(er));
};

any_graph::degree_size_type in_degree(any_graph::vertex_descriptor v, const any_graph& g) {
  return g.in_degree(v);
};

any_graph::degree_size_type degree(any_graph::edge_descriptor e, const any_graph& g) {
  return g.degree(e);
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

std::pair< any_graph::vertex_iterator, any_graph::vertex_iterator > vertices(const any_graph& g) {
  any_graph::vertex_range vr = g.vertices();
  return std::pair< any_graph::vertex_iterator, any_graph::vertex_iterator >(begin(vr), end(vr));
};

any_graph::vertices_size_type num_vertices(const any_graph& g) {
  return g.num_vertices();
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

std::pair< any_graph::edge_iterator, any_graph::edge_iterator > edges(const any_graph& g) {
  any_graph::edge_range er = g.edges();
  return std::pair< any_graph::edge_iterator, any_graph::edge_iterator >(begin(er), end(er));
};

any_graph::edges_size_type num_edges(const any_graph& g) {
  return g.num_edges();
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

std::pair<any_graph::edge_descriptor, bool> edge(any_graph::vertex_descriptor u, any_graph::vertex_descriptor v, const any_graph& g) {
  return g.edge(u, v);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

any_graph::vertex_descriptor add_vertex(any_graph& g) {
  return g.add_vertex();
};

void clear_vertex(any_graph::vertex_descriptor v, any_graph& g) {
  g.clear_vertex(v);
};

void remove_vertex(any_graph::vertex_descriptor v, any_graph& g) {
  g.remove_vertex(v);
};

std::pair<any_graph::edge_descriptor, bool> add_edge(any_graph::vertex_descriptor u, any_graph::vertex_descriptor v, any_graph& g) {
  return g.add_edge(u, v);
};

void remove_edge(any_graph::vertex_descriptor u, any_graph::vertex_descriptor v, any_graph& g) {
  std::pair<any_graph::edge_descriptor, bool> e = g.edge(u, v);
  if(e.second)
    g.remove_edge(e.first);
};

void remove_edge(any_graph::edge_descriptor e, any_graph& g) {
  g.remove_edge(e);
};

void remove_edge(any_graph::edge_iterator e_iter, any_graph& g) {
  g.remove_edge(*e_iter);
};





namespace detail {
  
  
  
  
  template <typename ErasedDesc, typename Iterator>
  class type_erased_graph_iterator : 
    public boost::iterator_facade< 
      type_erased_graph_iterator<ErasedDesc, Iterator>,
      ErasedDesc, 
      typename std::iterator_traits<ContIterator>::iterator_category,
      ErasedDesc
    > {
    public:
      typedef type_erased_graph_iterator<ErasedDesc, Iterator> self;
      
      explicit type_erased_graph_iterator(Iterator aIt = Iterator()) : it(aIt) { };
      
    public: // private:
      friend class boost::iterator_core_access;
      
      void increment() { ++it; };
      void decrement() { --it; };
      bool equal(const self& rhs) const { return this->it == rhs.it; };
      ErasedDesc dereference() const { return ErasedDesc(boost::any(*it)); };
      
      void advance(std::ptrdiff_t i) { it += i; };
      std::ptrdiff_t distance_to(const self& rhs) const { return rhs.it - this->it; }; 
      
      Iterator it;
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_edge_list_graph<Graph>,
  std::size_t >::type try_get_num_edges(const Graph& g) {
    return num_edges(g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_edge_list_graph<Graph>,
  std::size_t >::type try_get_num_edges(const Graph&) {
    return 0;
  };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_edge_list_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_edges(const Graph& g) {
    typedef boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > ResultType;
    typedef typename boost::graph_traits< Graph >::edge_iterator EdgeIter;
    typedef type_erased_graph_iterator<any_graph::edge_descriptor, EdgeIter> ErasedEdgeIter;
    std::pair< EdgeIter, EdgeIter > real_range = edges(g);
    return ResultType(ErasedEdgeIter(real_range.first), ErasedEdgeIter(real_range.second));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_edge_list_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_edges(const Graph&) {
    return boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t >();
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::mpl::or_< boost::is_edge_list_graph<Graph>, boost::is_incidence_graph<Graph> >,
  any_graph::vertex_descriptor >::type try_get_source(const any_graph::edge_descriptor& e, const Graph& g) {
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    return any_graph::vertex_descriptor(boost::any(source(boost::any_cast<Edge>(e), g)));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::mpl::or_< boost::is_edge_list_graph<Graph>, boost::is_incidence_graph<Graph> >,
  any_graph::vertex_descriptor >::type try_get_source(const any_graph::edge_descriptor&, const Graph&) {
    return any_graph::vertex_descriptor(boost::any(boost::graph_traits< Graph >::null_vertex()));
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::mpl::or_< boost::is_edge_list_graph<Graph>, boost::is_incidence_graph<Graph> >,
  any_graph::vertex_descriptor >::type try_get_target(const any_graph::edge_descriptor& e, const Graph& g) {
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    return any_graph::vertex_descriptor(boost::any(target(boost::any_cast<Edge>(e), g)));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::mpl::or_< boost::is_edge_list_graph<Graph>, boost::is_incidence_graph<Graph> >,
  any_graph::vertex_descriptor >::type try_get_target(const any_graph::edge_descriptor&, const Graph&) {
    return any_graph::vertex_descriptor(boost::any(boost::graph_traits< Graph >::null_vertex()));
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_adjacency_matrix<Graph>,
  std::pair< any_graph::edge_descriptor, bool> >::type try_get_edge(const any_graph::vertex_descriptor& u, const any_graph::vertex_descriptor& v, const Graph& g) {
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    std::pair< Edge, bool > result = edge(boost::any_cast<Vertex>(u), boost::any_cast<Vertex>(v), g);
    return std::pair< any_graph::edge_descriptor, bool>(any_graph::edge_descriptor(boost::any(result.first)), result.second);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_adjacency_matrix<Graph>,
  std::pair< any_graph::edge_descriptor, bool> >::type try_get_edge(const any_graph::vertex_descriptor&, const any_graph::vertex_descriptor&, const Graph&) {
    return std::pair< any_graph::edge_descriptor, bool>(any_graph::edge_descriptor(), false);
  };
  
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_incidence_graph<Graph>,
  std::size_t >::type try_get_out_degree(const any_graph::vertex_descriptor& u, const Graph& g) {
    return out_degree(u,g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_incidence_graph<Graph>,
  std::size_t >::type try_get_out_degree(const any_graph::vertex_descriptor&, const Graph&) {
    return 0;
  };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_incidence_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_out_edges(const any_graph::vertex_descriptor& u, const Graph& g) {
    typedef boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > ResultType;
    typedef typename boost::graph_traits< Graph >::out_edge_iterator EdgeIter;
    typedef type_erased_graph_iterator<any_graph::edge_descriptor, EdgeIter> ErasedEdgeIter;
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    std::pair< EdgeIter, EdgeIter > real_range = out_edges(boost::any_cast<Vertex>(u), g);
    return ResultType(ErasedEdgeIter(real_range.first), ErasedEdgeIter(real_range.second));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_incidence_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_out_edges(const any_graph::vertex_descriptor&, const Graph&) {
    return boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t >();
  };
  
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_bidirectional_graph<Graph>,
  std::size_t >::type try_get_in_degree(const any_graph::vertex_descriptor& u, const Graph& g) {
    return in_degree(u,g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_bidirectional_graph<Graph>,
  std::size_t >::type try_get_in_degree(const any_graph::vertex_descriptor&, const Graph&) {
    return 0;
  };
  
  template <typename Graph>
  typename boost::enable_if< boost::is_bidirectional_graph<Graph>,
  std::size_t >::type try_get_degree(const any_graph::edge_descriptor& e, const Graph& g) {
    return degree(e,g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_bidirectional_graph<Graph>,
  std::size_t >::type try_get_degree(const any_graph::edge_descriptor&, const Graph&) {
    return 0;
  };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_bidirectional_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_in_edges(const any_graph::vertex_descriptor& u, const Graph& g) {
    typedef boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > ResultType;
    typedef typename boost::graph_traits< Graph >::in_edge_iterator EdgeIter;
    typedef type_erased_graph_iterator<any_graph::edge_descriptor, EdgeIter> ErasedEdgeIter;
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    std::pair< EdgeIter, EdgeIter > real_range = in_edges(boost::any_cast<Vertex>(u), g);
    return ResultType(ErasedEdgeIter(real_range.first), ErasedEdgeIter(real_range.second));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_bidirectional_graph<Graph>,
  boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t > >::type 
    try_get_in_edges(const any_graph::vertex_descriptor&, const Graph&) {
    return boost::any_range< any_graph::edge_descriptor, boost::forward_pass_traversal_tag, any_graph::edge_descriptor, std::ptrdiff_t >();
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_vertex_list_graph<Graph>,
  std::size_t >::type try_get_num_vertices(const Graph& g) {
    return num_vertices(g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_vertex_list_graph<Graph>,
  std::size_t >::type try_get_num_vertices(const Graph&) {
    return 0;
  };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_vertex_list_graph<Graph>,
  boost::any_range< any_graph::vertex_descriptor, boost::forward_pass_traversal_tag, any_graph::vertex_descriptor, std::ptrdiff_t > >::type 
    try_get_vertices(const Graph& g) {
    typedef boost::any_range< any_graph::vertex_descriptor, boost::forward_pass_traversal_tag, any_graph::vertex_descriptor, std::ptrdiff_t > ResultType;
    typedef typename boost::graph_traits< Graph >::vertex_iterator VertexIter;
    typedef type_erased_graph_iterator<any_graph::vertex_descriptor, VertexIter> ErasedVertexIter;
    std::pair< VertexIter, VertexIter > real_range = vertices(g);
    return ResultType(ErasedVertexIter(real_range.first), ErasedVertexIter(real_range.second));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_vertex_list_graph<Graph>,
  boost::any_range< any_graph::vertex_descriptor, boost::forward_pass_traversal_tag, any_graph::vertex_descriptor, std::ptrdiff_t > >::type 
    try_get_vertices(const Graph&) {
    return boost::any_range< any_graph::vertex_descriptor, boost::forward_pass_traversal_tag, any_graph::vertex_descriptor, std::ptrdiff_t >();
  };
  
  
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_mutable_graph<Graph>,
  any_graph::vertex_descriptor >::type try_add_vertex(Graph& g) {
    return any_graph::vertex_descriptor(boost::any(add_vertex(g)));
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_mutable_graph<Graph>,
  std::size_t >::type try_add_vertex(Graph&) {
    return any_graph::vertex_descriptor(boost::any(boost::graph_traits<Graph>::null_vertex()));
  };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_mutable_graph<Graph>,
  any_graph::vertex_descriptor >::type try_clear_vertex(const any_graph::vertex_descriptor& u, Graph& g) {
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    clear_vertex(boost::any_cast<Vertex>(u), g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_mutable_graph<Graph>,
  void >::type try_clear_vertex(const any_graph::vertex_descriptor&, Graph&) { };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_mutable_graph<Graph>,
  void >::type try_remove_vertex(const any_graph::vertex_descriptor& u, Graph& g) {
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    remove_vertex(boost::any_cast<Vertex>(u), g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_mutable_graph<Graph>,
  void >::type try_remove_vertex(const any_graph::vertex_descriptor&, Graph&) { };
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_mutable_graph<Graph>,
  std::pair<any_graph::edge_descriptor, bool> >::type try_add_edge(const any_graph::vertex_descriptor& u, const any_graph::vertex_descriptor& v, Graph& g) const {
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    std::pair<edge, bool> result = add_edge(boost::any_cast<Vertex>(u), boost::any_cast<Vertex>(v), g);
    return std::pair<any_graph::edge_descriptor, bool>(any_graph::edge_descriptor(boost::any(result.first)), result.second);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_mutable_graph<Graph>,
  std::pair<any_graph::edge_descriptor, bool> >::type try_add_edge(const any_graph::vertex_descriptor&, const any_graph::vertex_descriptor&, Graph&) const {
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    return std::pair<any_graph::edge_descriptor, bool>(edge_descriptor(boost::any(Edge())), false);
  };
  
  
  
  template <typename Graph>
  typename boost::enable_if< boost::is_mutable_graph<Graph>,
  void >::type try_remove_edge(const any_graph::edge_descriptor& e, Graph& g) {
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    remove_edge(boost::any_cast<Edge>(e), g);
  };
  
  template <typename Graph>
  typename boost::disable_if< boost::is_mutable_graph<Graph>,
  void >::type try_remove_edge(const any_graph::edge_descriptor&, Graph&) { };
  
  
  
};




template <typename Graph>
class type_erased_graph : public any_graph {
  public:
    typedef any_graph::any_descriptor any_descriptor;
    typedef any_graph::vertex_descriptor vertex_descriptor;
    typedef any_graph::edge_descriptor edge_descriptor;
    
    typedef any_graph::any_bundled any_bundled;
    typedef any_graph::vertex_bundled vertex_bundled;
    typedef any_graph::edge_bundled edge_bundled;
    
    typedef any_graph::edge_range edge_range;
    typedef any_graph::out_edge_range out_edge_range;
    typedef any_graph::in_edge_range in_edge_range;
    
    typedef any_graph::vertex_range vertex_range;
    
    typedef any_graph::degree_size_type degree_size_type;
    typedef any_graph::vertices_size_type vertices_size_type;
    typedef any_graph::edges_size_type edges_size_type;
    
  protected:
    
    typedef typename boost::graph_traits<Graph>::vertex_descriptor real_vertex_desc;
    typedef typename boost::graph_traits<Graph>::edge_descriptor real_edge_desc;
    
    Graph* p_graph;
    
    virtual vertex_bundled get_vertex_bundled(const vertex_descriptor& aVertex) const {
      return vertex_bundled((*p_graph)[boost::any_cast<real_vertex_desc>(aVertex)]);
    };
    
    virtual edge_bundled get_edge_bundled(const edge_descriptor& aEdge) const {
      return edge_bundled((*p_graph)[boost::any_cast<real_edge_desc>(aEdge)]);
    };
    
  public:
    
    type_erased_graph(Graph* aPGraph) : p_graph(aPGraph) { };
    
    virtual edge_range edges() const {
      return detail::try_get_edges(*p_graph);
    };
    virtual edges_size_type num_edges() const {
      return detail::try_get_num_edges(*p_graph);
    };
    
    virtual vertex_descriptor source(const edge_descriptor& aEdge) const {
      return detail::try_get_source(aEdge, *p_graph);
    };
    virtual vertex_descriptor target(const edge_descriptor& aEdge) const {
      return detail::try_get_target(aEdge, *p_graph);
    };
    
    virtual std::pair<edge_descriptor, bool> edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const {
      return detail::try_get_edge(aU, aV, *p_graph);
    };
    
    
    virtual out_edge_range out_edges(const vertex_descriptor& aVertex) const {
      return detail::try_get_out_edges(aVertex, *p_graph);
    };
    
    virtual degree_size_type out_degree(const vertex_descriptor& aVertex) const {
      return detail::try_get_out_degree(aVertex, *p_graph);
    };
    
    
    virtual in_edge_range in_edges(const vertex_descriptor& aVertex) const {
      return detail::try_get_in_edges(aVertex, *p_graph);
    };
    
    virtual degree_size_type in_degree(const vertex_descriptor& aVertex) const {
      return detail::try_get_in_degree(aVertex, *p_graph);
    };
    
    virtual degree_size_type degree(const edge_descriptor& aEdge) const {
      return detail::try_get_degree(aEdge, *p_graph);
    };
    
    
    virtual vertex_range vertices() const {
      return detail::try_get_vertices(*p_graph);
    };
    
    virtual vertices_size_type num_vertices() const {
      return detail::try_get_num_vertices(*p_graph);
    };
    
    
    virtual vertex_descriptor add_vertex() const {
      return detail::try_add_vertex(*p_graph);
    };
    
    virtual void clear_vertex(const vertex_descriptor& aVertex) const {
      detail::try_clear_vertex(aVertex, *p_graph);
    };
    
    virtual void remove_vertex(const vertex_descriptor& aVertex) const {
      detail::try_remove_vertex(aVertex, *p_graph);
    };
    
    
    virtual std::pair<edge_descriptor, bool> add_edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const {
      return detail::try_add_edge(aU, aV, *p_graph);
    };
    
    virtual void remove_edge(const edge_descriptor& aEdge) const {
      detail::try_remove_edge(aEdge, *p_graph);
    };
    
    
};




};

};


#endif



