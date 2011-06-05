/**
 * \file probabilistic_roadmap.hpp
 * \author S. Mikael Persson <mikael.s.persson@gmail.com>
 *
  * This library
 *
 *
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

#ifndef PROBABILISTIC_ROADMAP_HPP
#define PROBABILISTIC_ROADMAP_HPP

#include <functional>
#include <vector>
#include <boost/limits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace ReaK {

namespace graph {

  /**
   * The PRMVisitorConcept is used to customize the behaviour of the PRM algorithm by 
   * implementing a set of required callback functions:
   *  vertex_added(Vertex, Graph&): Is called whenever a new vertex has been added to the graph, but not yet connected.
   *  edge_added(Edge, Graph&): Is called whenever a new edge has been created between the last created vertex and its nearest neighbor in the graph.
   *  is_position_free(const PositionValue&): Is called to query whether a particular configuration (position) is free.
   *  joining_vertex_found(Vertex, Graph&): Is called by the bidirectional RRT algorithm when a vertex is found that can reach the other graph (not the one passed as 2nd parameter).
   *  keep_going(): Is called at each attempt to expand the graph to verify that the user still wants more vertices to be generated in the tree. 
   */
  template <typename Visitor, typename Graph>
  struct PRMVisitorConcept {
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
    void constraints() {
      boost::function_requires< boost::CopyConstructibleConcept<Visitor> >();
      vis.vertex_added(u, g); 
      vis.edge_added(e, g);
      std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> v;
      vis.expand_vertex(u, g, v);
      vis.update_density(u, g);
    };
  };

  /**
   * This class is simply a "null" visitor for the PRM algorithm. It is null in the sense that it
   * will do nothing (except return true on the is_position_free(p) and keep_going() callbacks).
   */
  template <typename Topology, typename PositionMap>
  struct default_prm_visitor {
    default_prm_visitor(const Topology& free_space, PositionMap position) : m_free_space(free_space), m_position(position) {};
    default_prm_visitor(const default_prm_visitor<Topology,PositionMap>& aVis) : m_free_space(aVis.m_free_space), m_position(aVis.m_position) { };
    
    template <typename Vertex, typename Graph>
    void vertex_added(Vertex,Graph&) { };
    template <typename Edge, typename Graph>
    void edge_added(Edge,Graph&) { };
    template <typename Vertex, typename Graph>
    void expand_vertex(Vertex,Graph&,std::vector<Vertex>&) { };
    template <typename Vertex, typename Graph>
    void update_density(Vertex, Graph& g) { };
    
    const Topology& m_free_space;
    PositionMap m_position;
  };
  
  
  /**
   * This class is a composite visitor class template. It can be used to glue together a function pointer 
   * for each of the functions of the PRMVisitorConcept so that it can be used as a light-weight,
   * copyable visitor object for the PRM algorithm (it is especially recommend to use the 
   * make_composite_prm_visitor function template).
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename ExpandVertexFunction,
	    typename UpdateDensityFunction>
  struct composite_prm_visitor {
    VertexAddedCallback vertex_added;
    EdgeAddedCallback edge_added;
    ExpandVertexFunction expand_vertex;
    UpdateDensityFunction update_density;
    composite_prm_visitor(VertexAddedCallback aVertexAdded,
                          EdgeAddedCallback aEdgeAdded,
			  ExpandVertexFunction aExpandVertex,
			  UpdateDensityFunction aUpdateDensity) :
                          vertex_added(aVertexAdded), edge_added(aEdgeAdded),
                          expand_vertex(aExpandVertex), update_density(aUpdateDensity) { };
  };

  /**
   * This is a function template that is used to create an object of a class of the composite_prm_visitor
   * class template. This is particularly convenient to avoid explicitely providing the list of template
   * arguments and let the compiler resolved them from the function parameter types.
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename ExpandVertexFunction,
	    typename UpdateDensityFunction>
  inline composite_prm_visitor<VertexAddedCallback, EdgeAddedCallback,ExpandVertexFunction,UpdateDensityFunction>
    make_composite_prm_visitor(VertexAddedCallback aVertexAdded,
                               EdgeAddedCallback aEdgeAdded,
			       ExpandVertexFunction aExpandVertex,
			       UpdateDensityFunction aUpdateDensity) {
    return composite_prm_visitor<VertexAddedCallback,EdgeAddedCallback,ExpandVertexFunction,UpdateDensityFunction>(aVertexAdded,aEdgeAdded,aExpandVertex,aUpdateDensity);
  };
  
  
  
  namespace detail {
    
    template <typename Graph,
	      typename Topology,
	      typename PRMVisitor,
	      typename MutableQueue,
	      typename PositionMap,
	      typename NcSelector>
    inline typename boost::enable_if< boost::is_undirected_graph<Graph> >::type 
      connect_prm_node(Graph& g,
	               const Topology& free_space,
		       PRMVisitor& vis,
		       MutableQueue& Q,
		       PositionMap position,
		       typename boost::graph_traits<Graph>::vertex_descriptor u,
		       const NcSelector& select_neighborhood) {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      PositionValue p = boost::get(position, u);
      std::vector<Vertex> Nc; 
      select_neighborhood(p,Nc,g,free_space,position);
      for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
	if((u != *it) && (free_space.distance(boost::get(position,*it), p) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from *it.
	  std::pair<Edge, bool> ep = boost::add_edge(*it,u,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g); 
	    Q.update(*it); 
	  };
	};
      }; 
      vis.update_density(u, g); 
      Q.push(u); 
    };
  
    template <typename Graph,
	      typename Topology,
	      typename PRMVisitor,
	      typename MutableQueue,
	      typename PositionMap,
	      typename NcSelector>
    inline typename boost::enable_if< boost::is_directed_graph<Graph> >::type 
      connect_prm_node(Graph& g,
		       const Topology& free_space,
		       PRMVisitor& vis,
		       MutableQueue& Q,
		       PositionMap position,
		       typename boost::graph_traits<Graph>::vertex_descriptor u,
		       const NcSelector& select_neighborhood) {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      PositionValue p = boost::get(position, u);
      std::vector<Vertex> Pred;
      std::vector<Vertex> Succ;
      select_neighborhood(p,Pred,Succ,g,free_space,position);
      for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
	if((u != *it) && (free_space.distance(boost::get(position,*it), p) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from *it.
	  std::pair<Edge, bool> ep = boost::add_edge(*it,u,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g); 
	    Q.update(*it); 
	  };
	};
      }; 
      for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
        PositionValue p_succ = boost::get(position,*it);
	if((u != *it) && (free_space.distance(p, p_succ) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from *it.
	  std::pair<Edge, bool> ep = boost::add_edge(u,*it,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g);
	    Q.update(*it);
	  };
	};
      }; 
      vis.update_density(u, g); 
      Q.push(u);
    };
    
  
  }; //end of detail namespace
  
  
  
  /**
   * This function is the basic PRM algorithm (as of "Geraerts and Overmars, 2002"). Note that 
   * all the design decisions (how to generate vertices, how to select the neighborhood, how to 
   * check the collision) have been externalized via the provided free_space topology, visitor 
   * object, and neighborhood selection functor.
   * \param g A mutable graph that should initially store the starting 
   *          vertex (if not it will be randomly generated) and will store 
   *          the generated graph once the algorithm has finished.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *                   that it is required to generate only random points in 
   *                   the free-space and to only allow interpolation within the free-space.
   * \param vis A PRM visitor implementing the PRMVisitorConcept. This is the 
   *            main point of customization and recording of results that the 
   *            user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *                 the value_type of this map should be the same type as the topology's 
   *                 value_type.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *                            vertices of the graph that ought to be connected to a new 
   *                            vertex. The list should be sorted in order of increasing "distance".
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *                         should stop regardless of whether the resulting graph is satisfactory
   *                         or not.
   * 
   */
  template <typename Graph,
	    typename Topology,
	    typename PRMVisitor,
	    typename PositionMap,
	    typename DensityMap,
	    typename NcSelector,
	    typename RunningPredicate,
	    typename CompareFunction>
  inline void generate_prm(Graph& g,
			   const Topology& free_space,
			   PRMVisitor vis,
			   PositionMap position,
			   DensityMap density,
			   const NcSelector& select_neighborhood,
			   unsigned int max_vertex_count,
			   unsigned int num_constructed_vertices, 
			   unsigned int num_expanded_vertices,
			   RunningPredicate keep_going,
			   CompareFunction compare) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    boost::function_requires< PRMVisitorConcept<PRMVisitor,Graph> >();
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    if(boost::num_vertices(g) == 0) {
      Vertex u = boost::add_vertex(g);
      PositionValue p = free_space.random_point();
      put(position, u, p);
      vis.vertex_added(u, g);
    };
    
    typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
    IndexInHeapMap index_in_heap;
    
    typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, DensityMap, CompareFunction> MutableQueue;
    MutableQueue Q(density, index_in_heap, compare); //priority queue holding the "expandable" vertices.
    
    {
      typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
      for (boost::tie(ui, ui_end) = boost::vertices(g); ui != ui_end; ++ui) {
	vis.update_density(*ui,g);
	Q.push(*ui);
      };
    };

    while(boost::num_vertices(g) < max_vertex_count) {
      //Graph Construction phase:
      unsigned int i = 0;
      while((i < num_constructed_vertices) && (boost::num_vertices(g) < max_vertex_count) && (keep_going())) {
        PositionValue p_rnd = free_space.random_point();
        Vertex u = boost::add_vertex(g); 
        put(position, u, p_rnd);             
        vis.vertex_added(u,g); 
        detail::connect_prm_node(g,free_space,vis,Q,position,u,select_neighborhood);
	++i;
      };
      
      //Node Expansion phase:
      unsigned int j = 0;
      while((j < num_expanded_vertices) && (boost::num_vertices(g) < max_vertex_count) && (keep_going())) {
	//use the priority queue to get the vertices that need expansion.
	Vertex v = Q.top();
	std::vector<Vertex> u_list;
	vis.expand_vertex(v, g, u_list);
	for(typename std::vector<Vertex>::iterator ui = u_list.begin(); ui != u_list.end(); ++ui) {
  	  Vertex u = *ui;     
          vis.vertex_added(u,g);  
          detail::connect_prm_node(g,free_space,vis,Q,position,u,select_neighborhood);
	};
	vis.update_density(v, g); 
	Q.update(v);
	++j;
      };
    };

  };
  


};


};

#endif










