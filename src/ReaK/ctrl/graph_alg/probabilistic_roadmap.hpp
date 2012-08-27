/**
 * \file probabilistic_roadmap.hpp
 *
 * This library provides function templates and concepts that implement the probabilistic roadmap (PRM) 
 * algorithm (as of "Geraerts and Overmars, 2002"). A PRM is generated in two phases. First, a number of 
 * vertices are generated at random from the free-space (configuration-space which is not occupied by 
 * obstacles) and bi-directional (or undirected) connections between those vertices are established, as
 * much as possible (if they can be directly connected through free-space). Second, difficult areas 
 * of the roadmap are determined using some density function (heuristic) and a priority-queue of vertices
 * with high density (whatever that means for the heuristic chosen). Then, high-density vertices are 
 * selected for expansion, in which case a number of neighboring vertices are generated coming out of 
 * the selected vertex into free-space. The first phase is called the construction phase, and the 
 * second is called the expansion phase. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, how to measure density, how to compare 
 * densities, when to stop the algorithm, what proportion of constructed vs. expanded vertices should 
 * be generated, etc. All these customization points are left to the user to implement, some are 
 * defined by the PRMVisitorConcept (expand-vertex and compute-density) while others are provided 
 * as functors to the function template that generates the PRM (generate_prm).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_PROBABILISTIC_ROADMAP_HPP
#define REAK_PROBABILISTIC_ROADMAP_HPP

#include <functional>
#include <vector>
#include <boost/limits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"


namespace ReaK {

namespace graph {

  /**
   * This concept class defines what is required of a class to serve as a visitor to the PRM algorithm.
   * 
   * Required concepts:
   * 
   * the visitor should model the boost::CopyConstructibleConcept.
   * 
   * Valid expressions:
   * 
   * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not yet connected.
   * 
   * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created vertex and its neighbor in the graph (g).
   *
   * boost::tie(pt,b) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph.
   * 
   * vis.update_density(u, g);  This function is called to evaluate the density-measure of the graph (g) around the given vertex (u). This value is used to prioritize the generation of new vertices.
   * 
   * \tparam Visitor The visitor class to be checked for modeling this concept.
   * \tparam Graph The graph on which the visitor class is required to work with.
   */
  template <typename Visitor, typename Graph, typename Topology>
  struct PRMVisitorConcept {
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
    typename ReaK::pp::topology_traits<Topology>::point_type pt;
    bool b;
    
    BOOST_CONCEPT_USAGE(PRMVisitorConcept) {
      BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
      vis.vertex_added(u, g); 
      vis.edge_added(e, g);
      std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> v;
      boost::tie(pt,b) = vis.random_walk(u, g);
      vis.update_density(u, g);
      b = vis.keep_going();
    };
  };

  /**
   * This class is simply a "null" visitor for the PRM algorithm. It is null in the sense that it
   * will do nothing on all accounts.
   * \tparam Topology The topology type that represents the free-space.
   * \tparam PositionMap The property-map type which can store the position associated to each vertex.
   */
  template <typename Topology, typename PositionMap>
  struct default_prm_visitor {
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PointType;
    
    default_prm_visitor(const Topology& free_space, PositionMap position) : m_free_space(free_space), m_position(position) {};
    default_prm_visitor(const default_prm_visitor<Topology,PositionMap>& aVis) : m_free_space(aVis.m_free_space), m_position(aVis.m_position) { };
    
    template <typename Vertex, typename Graph>
    void vertex_added(Vertex,Graph&) { };
    template <typename Edge, typename Graph>
    void edge_added(Edge,Graph&) { };
    template <typename Vertex, typename Graph>
    std::pair<PointType,bool> random_walk(Vertex,Graph&) { return std::make_pair(PointType(),false); };
    template <typename Vertex, typename Graph>
    void update_density(Vertex, Graph& g) { };
    bool keep_going() { return true; };
    
    const Topology& m_free_space;
    PositionMap m_position;
  };
  
  
  /**
   * This class is a composite visitor class template. It can be used to glue together a function pointer (or functor) 
   * for each of the functions of the PRMVisitorConcept so that it can be used as a light-weight,
   * copyable visitor object for the PRM algorithm (it is especially recommended to use the 
   * make_composite_prm_visitor function template).
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename RandomWalkerFunction,
	    typename UpdateDensityFunction,
            typename KeepGoingFunction>
  struct composite_prm_visitor {
    VertexAddedCallback vertex_added;
    EdgeAddedCallback edge_added;
    RandomWalkerFunction random_walk;
    UpdateDensityFunction update_density;
    KeepGoingFunction keep_going;
    composite_prm_visitor(VertexAddedCallback aVertexAdded,
                          EdgeAddedCallback aEdgeAdded,
			  RandomWalkerFunction aRandomWalk,
			  UpdateDensityFunction aUpdateDensity,
                          KeepGoingFunction aKeepGoing) :
                          vertex_added(aVertexAdded), edge_added(aEdgeAdded),
                          random_walk(aRandomWalk), update_density(aUpdateDensity), 
                          keep_going(aKeepGoing) { };
  };

  /**
   * This is a function template that is used to create an object of a class of the composite_prm_visitor
   * class template. This is particularly convenient to avoid explicitely providing the list of template
   * arguments and let the compiler resolved them from the function parameter types.
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename RandomWalkerFunction,
	    typename UpdateDensityFunction,
            typename KeepGoingFunction>
  inline composite_prm_visitor<VertexAddedCallback, EdgeAddedCallback,RandomWalkerFunction,UpdateDensityFunction,KeepGoingFunction>
    make_composite_prm_visitor(VertexAddedCallback aVertexAdded,
                               EdgeAddedCallback aEdgeAdded,
			       RandomWalkerFunction aRandomWalk,
			       UpdateDensityFunction aUpdateDensity,
                               KeepGoingFunction aKeepGoing) {
    return composite_prm_visitor<VertexAddedCallback,EdgeAddedCallback,RandomWalkerFunction,UpdateDensityFunction,KeepGoingFunction>(aVertexAdded,aEdgeAdded,aRandomWalk,aUpdateDensity,aKeepGoing);
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
                       const typename boost::property_traits<PositionMap>::value_type& p,
		       const NcSelector& select_neighborhood) {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      typedef typename Graph::vertex_bundled VertexProp;
      
      typedef boost::composite_property_map< 
        PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
      GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p,std::back_inserter(Nc),g,free_space,g_position);
      
      VertexProp up;
      put(position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      vis.vertex_added(u,g); 
      
      for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
	if((u != *it) && (get(ReaK::pp::distance_metric, free_space)(get(position,g[*it]), p, free_space) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from *it.
	  std::pair<Edge, bool> ep = add_edge(*it,u,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g); 
	    Q.push_or_update(*it); 
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
                       const typename boost::property_traits<PositionMap>::value_type& p,
		       const NcSelector& select_neighborhood) {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      typedef typename Graph::vertex_bundled VertexProp;
      
      typedef boost::composite_property_map< 
        PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
      GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
      
      std::vector<Vertex> Pred;
      std::vector<Vertex> Succ;
      select_neighborhood(p,std::back_inserter(Pred),std::back_inserter(Succ),g,free_space,g_position);
      
      VertexProp up;
      put(position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      vis.vertex_added(u,g); 
      
      for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
	if((u != *it) && (get(ReaK::pp::distance_metric, free_space)(get(position,g[*it]), p, free_space) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from *it.
	  std::pair<Edge, bool> ep = add_edge(*it,u,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g); 
	    Q.push_or_update(*it); 
	  };
	};
      }; 
      for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
	if((u != *it) && (get(ReaK::pp::distance_metric, free_space)(p, get(position,g[*it]), free_space) != std::numeric_limits<double>::infinity())) {
	  //this means that u is reachable from p.
	  std::pair<Edge, bool> ep = add_edge(u,*it,g); 
	  if(ep.second) { 
	    vis.edge_added(ep.first, g); 
	    vis.update_density(*it, g);
	    Q.push_or_update(*it);
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
   * 
   * \tparam Graph A mutable graph type that can store the roadmap, should model boost::MutableGraphConcept 
   *         and boost::VertexListGraphConcept (either bidirectional or undirected graph, the algorithm 
   *         will deal with either cases as appropriate).
   * \tparam Topology A topology type on which the vertex positions lie, should model the TopologyConcept,
   *         the MetricSpaceConcept (has an attached distance-metric) and the PointDistributionConcept (has 
   *         an attached random-sampler functor).
   * \tparam PRMVisitor A PRM visitor type, should model the PRMVisitorConcept.
   * \tparam PositionMap A property-map type that can store the position of each vertex. 
   * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
   *         See classes in the topological_search.hpp header-file.
   * \tparam CompareFunction A functor type that can be used to compare density values (strict weak-ordering).
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated graph once the algorithm has finished.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *        that it is required to generate only random points in 
   *        the free-space and to only allow interpolation within the free-space.
   * \param vis A PRM visitor implementing the PRMVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param density A property-map that provides the density values assiciated to each vertex.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting graph is satisfactory or not.
   * \param num_constructed_vertices The number of vertices to generate in the PRM construction
   *        phase (i.e. vertices are randomly generated all over the free-space).
   * \param num_expanded_vertices The number of vertices to generate in the PRM expansion phase 
   *        where the priority-queue (sorted by density and compare) is used to select vertices for 
   *        expansion (the expansion is done with vis.expand_vertex(u,g,v), see PRMVisitorConcept).
   * \param compare A functor used to compare density values (strict weak-ordering) in the priority-queue 
   *        for expansion of the vertices.
   */
  template <typename Graph,
	    typename Topology,
	    typename PRMVisitor,
	    typename PositionMap,
            typename RandomSampler,
	    typename DensityMap,
	    typename NcSelector,
	    typename CompareFunction>
  inline void generate_prm(Graph& g,
			   const Topology& free_space,
			   PRMVisitor vis,
			   PositionMap position,
                           RandomSampler get_sample,
			   DensityMap density,
			   const NcSelector& select_neighborhood,
			   unsigned int max_vertex_count,
			   unsigned int num_constructed_vertices, 
			   unsigned int num_expanded_vertices,
			   CompareFunction compare) {
    BOOST_CONCEPT_ASSERT((PRMVisitorConcept<PRMVisitor,Graph,Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>)); // for the distance-metric.
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>)); // for the random-sampler.
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    //BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename Graph::vertex_bundled VertexProp;
    
    typedef boost::composite_property_map< 
      DensityMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphDensityMap;
    GraphDensityMap g_density = GraphDensityMap(density, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
    
    if(num_vertices(g) == 0) {
      VertexProp up;
      PositionValue p = get_sample(free_space);
      put(position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      vis.vertex_added(u, g);
    };
    
    typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
    IndexInHeapMap index_in_heap;
    
    typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, GraphDensityMap, CompareFunction> MutableQueue;
    MutableQueue Q(g_density, index_in_heap, compare); //priority queue holding the "expandable" vertices.
    
    {
      typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
      for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
	vis.update_density(*ui,g);
	Q.push(*ui);
      };
    };

    while((num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
      //Graph Construction phase:
      unsigned int i = 0;
      while((i < num_constructed_vertices) && (num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
        PositionValue p_rnd = get_sample(free_space);
        
        detail::connect_prm_node(g,free_space,vis,Q,position,p_rnd,select_neighborhood);
        
	++i;
      };
      
      //Node Expansion phase:
      unsigned int j = 0;
      while((j < num_expanded_vertices) && (num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
	//use the priority queue to get the vertices that need expansion.
	Vertex v = Q.top();
        PositionValue p_rnd; bool expanding_worked;
        boost::tie(p_rnd, expanding_worked) = vis.random_walk(v, g);
        
        if(expanding_worked) {
          
          detail::connect_prm_node(g,free_space,vis,Q,position,p_rnd,select_neighborhood);
          
          vis.update_density(v, g); 
          Q.update(v);
          ++j;
        } else {
          Q.pop(); // if one cannot expand from this vertex, then leave it out of future attempts.
        };
      };
    };

  };
  
  
  


};


};

#endif










