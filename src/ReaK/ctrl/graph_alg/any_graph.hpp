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
    typedef any_descriptor vertex_descriptor;
    typedef any_descriptor edge_descriptor;
    
    typedef boost::any any_bundled;
    typedef any_bundled vertex_bundled;
    typedef any_bundled edge_bundled;
    
    typedef boost::any_range<edge_descriptor, boost::forward_pass_traversal_tag, edge_descriptor, std::ptrdiff_t > edge_range;
    typedef edge_range out_edge_range;
    typedef edge_range in_edge_range;
    
    typedef boost::any_range<vertex_descriptor, boost::forward_pass_traversal_tag, vertex_descriptor, std::ptrdiff_t > vertex_range;
    
    typedef std::size_t degree_size_type;
    typedef std::size_t vertices_size_type;
    typedef std::size_t edges_size_type;
    
  protected:
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const any_descriptor& aElement) const = 0;
    virtual boost::any get_property_by_any(const std::string& aProperty, const any_descriptor& aElement) const = 0;
    
  public:
    
    virtual any_bundled operator[](const any_descriptor& aElement) const = 0;
    
    virtual edge_range edges() const = 0;
    virtual edges_size_type num_edges() const = 0;
    
    virtual vertex_descriptor source(const edge_descriptor& aEdge) const = 0;
    virtual vertex_descriptor target(const edge_descriptor& aEdge) const = 0;
    
    virtual edge_descriptor edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const = 0;
    
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
    
    virtual edge_descriptor add_edge() const = 0;
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
    
    Graph* p_graph;
    
  public:
    
    type_erased_graph(Graph* aPGraph) : p_graph(aPGraph) { };
    
    virtual any_bundled operator[](const any_descriptor& aElement) const = 0;
    
    virtual edge_range edges() const = 0;
    virtual edges_size_type num_edges() const = 0;
    
    virtual vertex_descriptor source(const edge_descriptor& aEdge) const = 0;
    virtual vertex_descriptor target(const edge_descriptor& aEdge) const = 0;
    
    virtual edge_descriptor edge(const vertex_descriptor& aU, const vertex_descriptor& aV) const = 0;
    
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
    
    virtual edge_descriptor add_edge() const = 0;
    virtual void remove_edge(const edge_descriptor& aEdge) const = 0;
    
    
};




};

};


#endif



