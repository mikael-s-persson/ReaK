/**
 * \file rrt_star.hpp
 * 
 * This library contains the Rapidly-Exploring Random Tree generation algorithm with an
 * A* rewiring policy.
 * This is a method to create a random tree that will span over a non-convex space
 * as rapidly as possible. The method relies on a simple randomized insertion algorithm.
 * At each step, a random point is picked from the underlying topology (i.e. configuration 
 * space in path-planning terms). Then, the points in the current graph that are nearest 
 * to the random point are picked for expansion. Then, edges (of a maximum length) are 
 * added to the nearest vertex towards the random point while it is still possible to
 * add such an edge without leaving the free space (the part of the configuration space which 
 * is not occupied by an obstacle). Finally, the last point of the expansion is used to attempt 
 * the same expansion with the other nearest neighbors and the accumulated cost-to-go of all 
 * these alternate paths are compared and the edge that leads to the shortest accumulated cost-to-go
 * is selected as a new edge in the tree. The algorithm will stop when either the number of vertices
 * in the tree has reached a maximum or when the user callback signals the stop.
 * 
 * This library also provides the bidirectional version of the RRT* algorithm. In this version,
 * two trees are generated. Typically, one tree is initialized with the starting vertex and 
 * the other is initialized with the goal vertex. The algorithm works to try and join the 
 * two graphs as quickly as possible and with the most direct path. The algorithm alternates
 * between the two graphs. It first uses the normal procedure (as in the unidirectional variant)
 * to try and add a vertex to one graph. If it does not succeed (i.e. there was no free-space in 
 * the expanded direction), it will try to expand the other graph. If it does succeed,
 * than the last vertex that was added to the tree will be the point towards which the other tree 
 * will be expanded (if free-space permits). In other words, any successful vertex addition to one
 * tree causes the other to attempt an expansion in that direction, and any failure at adding a vertex
 * to one tree causes the other to attempt to expand towards a random point. This version of the 
 * algorithm will also notify the user (via a visitor's callback) whenever one tree was successfully
 * expanded towards the other tree to the point that they meet at two vertices (one on each graph).
 * The user can thus record successful connections as paths and decide whether it's worth continuing 
 * with the generation the Bi-RRT in the hopes of finding a better path with richer trees.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
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

#ifndef REAK_RRT_STAR_HPP
#define REAK_RRT_STAR_HPP

#include <functional>
#include <vector>
#include <boost/limits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"

#include "rr_graph.hpp"

namespace ReaK {
  
namespace graph {
  
  
  /**
   * This concept class defines what is required of a class to serve as a visitor to the RRT* algorithm.
   * 
   * Required concepts:
   * 
   * the visitor should model the RRGVisitorConcept.
   * 
   * Valid expressions:
   *
   * \tparam Visitor The visitor class to be checked for modeling this concept.
   * \tparam Graph The graph on which the visitor class is required to work with.
   * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
   */
  template <typename Visitor, typename Graph, typename PositionMap>
  struct RRTStarVisitorConcept : RRGVisitorConcept<Visitor,Graph,PositionMap> {
    BOOST_CONCEPT_USAGE(RRTStarVisitorConcept) {
    };
  };

  /**
   * This class is simply a "null" visitor for the RRT* algorithm. It is null in the sense that it
   * will do nothing (always returns true for the can_be_connected function).
   */
  struct default_rrt_star_visitor : default_rrg_visitor {
  };

  /**
   * This class is a composite visitor class template. It can be used to glue together a function pointer (or functor)
   * for each of the functions of the RRTStarVisitorConcept so that it can be used as a light-weight,
   * copyable visitor object for the RRT* algorithm (it is especially recommend to use the 
   * make_composite_rrg_visitor function template).
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
            typename SteerFunction,
            typename IsFreeQuery,
            typename KeepGoingQuery,
            typename CanConnectQuery>
  struct composite_rrt_star_visitor : composite_rrg_visitor<VertexAddedCallback, 
                                                            EdgeAddedCallback, 
                                                            SteerFunction,
                                                            IsFreeQuery,
                                                            KeepGoingQuery,
                                                            CanConnectQuery> {
    CanConnectQuery can_be_connected;
    composite_rrt_star_visitor(VertexAddedCallback aVertexAdded,
                               EdgeAddedCallback aEdgeAdded,
                               SteerFunction aSteerTowardsPosition,
                               IsFreeQuery aIsFree,
                               KeepGoingQuery aKeepGoing,
                               CanConnectQuery aCanBeConnected) :
                               composite_rrg_visitor<
                                 VertexAddedCallback, 
                                 EdgeAddedCallback, 
                                 SteerFunction,
                                 IsFreeQuery,
                                 KeepGoingQuery,
                                 CanConnectQuery>(aVertexAdded, aEdgeAdded, 
                                                  aSteerTowardsPosition, aIsFree, 
                                                  aKeepGoing, aCanBeConnected) { };
  };

  /**
   * This is a function template that is used to create an object of a class of the composite_rrg_visitor
   * class template. This is particularly convenient to avoid explicitely providing the list of template
   * arguments and let the compiler resolved them from the function parameter types.
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
            typename SteerFunction,
            typename IsFreeQuery,
            typename KeepGoingQuery,
            typename CanConnectQuery>
  inline composite_rrt_star_visitor<VertexAddedCallback, 
                                    EdgeAddedCallback, 
                                    SteerFunction, 
                                    IsFreeQuery, 
                                    KeepGoingQuery,
                                    CanConnectQuery>
    make_composite_rrt_star_visitor(VertexAddedCallback aVertexAdded,
                                    EdgeAddedCallback aEdgeAdded,
                                    SteerFunction aSteerTowardsPosition,
                                    IsFreeQuery aIsFree,
                                    KeepGoingQuery aKeepGoing,
                                    CanConnectQuery aCanBeConnected) {
    return composite_rrt_star_visitor<VertexAddedCallback, 
                                      EdgeAddedCallback, 
                                      SteerFunction, 
                                      IsFreeQuery, 
                                      KeepGoingQuery,
                                      CanConnectQuery>(aVertexAdded,
                                                       aEdgeAdded,
                                                       aSteerTowardsPosition,
                                                       aIsFree,
                                                       aKeepGoing,
                                                       aCanBeConnected);
  };


  
namespace detail {
  
  
  
  template <typename Graph,
            typename Topology,
            typename RRTStarVisitor,
            typename PositionMap,
            typename CostMap,
            typename PredecessorMap,
            typename WeightMap,
            typename DistanceMetric>
  inline typename boost::enable_if< boost::is_undirected_graph<Graph> >::type
    rewire_tree_neighborhood(
      typename boost::graph_traits<Graph>::vertex_descriptor x_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Nc,
      Graph& g,
      const Topology& space,
      RRTStarVisitor vis,
      PositionMap position,
      CostMap cost,
      PredecessorMap pred,
      WeightMap weight,
      DistanceMetric distance) {
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter; 
    typedef typename std::vector< Vertex >::const_iterator NcIter;
    typedef typename Graph::vertex_bundled VertexProp;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::back_inserter;

    std::queue< Vertex > incons_set;
    double c_near = get(cost, g[x_new]);
    PositionValue p_new = get(position, g[x_new]);
    for(NcIter it = Nc.begin(); it != Nc.end(); ++it) {
      double c_temp, d_temp;
      if((x_new != get(pred, g[*it])) &&
         (vis.can_be_connected(x_new, *it, g)) &&
         ((c_temp = c_near + (d_temp = distance(p_new, get(position, g[*it]), space))) < get(cost, g[*it]))) {
        
        EdgeProp eprop;
        put(weight, eprop, d_temp);
#ifdef RK_ENABLE_CXX0X_FEATURES
        std::pair<Edge,bool> ep = add_edge(x_new, *it, std::move(eprop), g);
#else
        std::pair<Edge,bool> ep = add_edge(x_new, *it, eprop, g);
#endif
        if(ep.second) {
          put(cost, g[*it], c_temp);
          Vertex x_pred = get(pred, g[*it]);
          put(pred, g[*it], x_new);
          vis.edge_added(ep.first, g);
          if((x_pred != *it) && 
             (x_pred != boost::graph_traits<Graph>::null_vertex()) && 
             (get(pred, g[x_pred]) != *it))
            remove_edge(x_pred, *it, g);
        };
        
        OutEdgeIter ei, ei_end;
        for(boost::tie(ei, ei_end) = out_edges(*it,g); ei != ei_end; ++ei) {
          Vertex x_next = target(*ei,g);
          if(x_next != x_new)
            incons_set.push(x_next);
        };
      };
    };
    
    while(!incons_set.empty()) {
      while(in_degree(incons_set.front(), g) == 0)
        incons_set.pop();
      Vertex x_pred = get(pred, g[incons_set.front()]);
      Edge e; bool e_exists;
      boost::tie(e,e_exists) = edge(x_pred, incons_set.front(), g);
      if((e_exists) &&
         (get(cost, g[x_pred]) + get(weight, g[e]) < get(cost, g[incons_set.front()]))) {
        put(cost, g[incons_set.front()], get(cost, g[x_pred]) + get(weight, g[e]));
        OutEdgeIter ei, ei_end;
        for(boost::tie(ei, ei_end) = out_edges(incons_set.front(),g); ei != ei_end; ++ei) {
          Vertex x_next = target(*ei,g);
          if(x_next != x_pred)
            incons_set.push(x_next);
        };
      };
      incons_set.pop();
    };
  };
  
  
  
  template <typename Graph,
            typename Topology,
            typename RRTStarVisitor,
            typename PositionMap,
            typename CostMap,
            typename WeightMap,
            typename DistanceMetric>
  inline typename boost::enable_if< boost::is_directed_graph<Graph> >::type
    rewire_tree_neighborhood(
      typename boost::graph_traits<Graph>::vertex_descriptor x_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Nc,
      Graph& g,
      const Topology& space,
      RRTStarVisitor vis,
      PositionMap position,
      CostMap cost,
      WeightMap weight,
      DistanceMetric distance) {
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter; 
    typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter; 
    typedef typename std::vector< Vertex >::const_iterator NcIter;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::back_inserter;

    std::queue< Vertex > incons_set;
    double c_near = get(cost, g[x_new]);
    PositionValue p_new = get(position, g[x_new]);
    for(NcIter it = Nc.begin(); it != Nc.end(); ++it) {
      double c_temp, d_temp;
      if((vis.can_be_connected(x_new, *it, g)) &&
         ((c_temp = c_near + (d_temp = distance(p_new, get(position, g[*it]), space))) < get(cost, g[*it]))) {
        
        EdgeProp eprop;
        put(weight, eprop, d_temp);
#ifdef RK_ENABLE_CXX0X_FEATURES
        std::pair<Edge,bool> ep = add_edge(x_new, *it, std::move(eprop), g);
#else
        std::pair<Edge,bool> ep = add_edge(x_new, *it, eprop, g);
#endif
        if(ep.second) {
          put(cost, g[*it], c_temp);
          vis.edge_added(ep.first, g);
          InEdgeIter iei, iei_end;
          for(boost::tie(iei, iei_end) = in_edges(*it,g); iei != iei_end; ++iei) {
            if(*iei != ep.first) {
              remove_edge(*iei, g);
              boost::tie(iei, iei_end) = in_edges(*it,g);
            };
          };
        };
        
        OutEdgeIter ei, ei_end;
        for(boost::tie(ei, ei_end) = out_edges(*it,g); ei != ei_end; ++ei)
          incons_set.push(target(*ei,g));
      };
    };
    
    while(!incons_set.empty()) {
      while(in_degree(incons_set.front(), g) == 0)
        incons_set.pop();
      Edge e = *(in_edges(incons_set.front(),g).first);
      Vertex x_pred = source(e,g);
      if(get(cost, g[x_pred]) + get(weight, g[e]) < get(cost, g[incons_set.front()])) {
        put(cost, g[incons_set.front()], get(cost, g[x_pred]) + get(weight, g[e]));
        OutEdgeIter ei, ei_end;
        for(boost::tie(ei, ei_end) = out_edges(incons_set.front(),g); ei != ei_end; ++ei)
          incons_set.push(target(*ei,g));
      };
      incons_set.pop();
    };
  };

};



  /**
   * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
   * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
   * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
   * \tparam RRTStarVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRTVisitorConcept.
   * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
   * \tparam CostMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or goal).
   * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
   * \tparam DistanceMetric This is a distance-metric over the topology (see pp::DistanceMetricConcept).
   * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
   * \param g A mutable graph that should initially store the starting and goal 
   *        vertex and will store the generated graph once the algorithm has finished.
   * \param space A topology (as defined by the Boost Graph Library). Note 
   *        that it is not required to generate only random points in 
   *        the free-space.
   * \param vis A RRT* visitor implementing the RRTStarVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param cost The property-map which stores the estimated cost-to-go of each vertex to the start (or goal).
   * \param pred The property-map which stores the predecessor of each vertex.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
   * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
   * \param distance A distance metric object to compute the distance between two positions.
   * \param select_neighborhood A callable object (functor) which can perform a 
   *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting tree is satisfactory or not.
   * 
   */
  template <typename Graph,
            typename Topology,
            typename RRTStarVisitor,
            typename PositionMap,
            typename CostMap,
            typename PredecessorMap,
            typename WeightMap,
            typename RandomSampler,
            typename DistanceMetric,
            typename NcSelector>
  inline typename boost::enable_if< boost::is_undirected_graph<Graph> >::type
    generate_rrt_star(Graph& g,
                      const Topology& space,
                      RRTStarVisitor vis,
                      PositionMap position,
                      CostMap cost,
                      PredecessorMap pred,
                      WeightMap weight,
                      RandomSampler get_sample,
                      DistanceMetric distance,
                      NcSelector select_neighborhood,
                      unsigned int max_vertex_count) {
    BOOST_CONCEPT_ASSERT((RRTStarVisitorConcept<RRTStarVisitor,Graph,PositionMap>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::DistanceMetricConcept<DistanceMetric,Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    typedef typename Graph::vertex_bundled VertexProp;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::back_inserter;
    
    typedef boost::composite_property_map< 
      PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
    GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));

    if(num_vertices(g) == 0) {
      PositionValue p = get_sample(space);
      while(!vis.is_position_free(p))
        p = get_sample(space);
      VertexProp up;
      put(position, up, p);
      put(cost, up, 0.0);
      put(pred, up, boost::graph_traits<Graph>::null_vertex());
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up),g);
#else
      Vertex u = add_vertex(up,g);
#endif
      vis.vertex_added(u, g);
    };

    while((num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
      
      PositionValue p_new = get_sample(space);
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p_new, back_inserter(Nc), g, space, g_position);
      //std::cout << "Trying to expand a vertex... " << Nc.size() << std::endl;
      Vertex x_near = boost::graph_traits<Graph>::null_vertex();
      if( ! detail::expand_to_nearest(x_near, p_new, Nc, g, vis) )
        continue;
      //std::cout << "Expanded a vertex." << std::endl;
      
      Nc.clear();
      select_neighborhood(p_new, back_inserter(Nc), g, space, g_position);
      
      VertexProp xp;
      put(position, xp, p_new);
      put(cost, xp, std::numeric_limits<double>::infinity());
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex x_new = add_vertex(std::move(xp),g);
#else
      Vertex x_new = add_vertex(xp,g);
#endif
      put(pred, g[x_new], x_new);
      vis.vertex_added(x_new,g);
      
      // Choose Parent:
      double d_near = distance(get(position, g[x_near]), p_new, space);
      double c_near = get(cost, g[x_near]) + d_near;
      for(typename std::vector<Vertex>::const_iterator it = Nc.begin(); it != Nc.end(); ++it) {
        double c_temp, d_temp;
        if((*it != x_near) && 
           (vis.can_be_connected(*it, x_new, g)) &&
           ((c_temp = get(cost, g[*it]) + (d_temp = distance(get(position, g[*it]), p_new, space))) < c_near)) {
          x_near = *it;
          c_near = c_temp;
          d_near = d_temp;
        };
      };
      
      EdgeProp eprop;
      put(weight, eprop, d_near);
#ifdef RK_ENABLE_CXX0X_FEATURES
      std::pair<Edge, bool> ep = add_edge(x_near, x_new, std::move(eprop), g);
#else
      std::pair<Edge, bool> ep = add_edge(x_near, x_new, eprop, g);
#endif
      if(!ep.second)
        continue;
      put(cost, g[x_new], c_near);
      put(pred, g[x_new], x_near);
      vis.edge_added(ep.first, g);
      
      detail::rewire_tree_neighborhood(x_new, Nc, g, space, vis, position, cost, pred, weight, distance); 
      
    };

  };

  
  
  
  /**
   * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
   * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
   * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
   * \tparam RRTStarVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRTVisitorConcept.
   * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
   * \tparam CostMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or goal).
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
   * \tparam DistanceMetric This is a distance-metric over the topology (see pp::DistanceMetricConcept).
   * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
   * \param g A mutable graph that should initially store the starting and goal 
   *        vertex and will store the generated graph once the algorithm has finished.
   * \param space A topology (as defined by the Boost Graph Library). Note 
   *        that it is not required to generate only random points in 
   *        the free-space.
   * \param vis A RRT* visitor implementing the RRTStarVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param cost The property-map which stores the estimated cost-to-go of each vertex to the start (or goal).
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
   * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
   * \param distance A distance metric object to compute the distance between two positions.
   * \param select_neighborhood A callable object (functor) which can perform a 
   *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting tree is satisfactory or not.
   * 
   */
  template <typename Graph,
            typename Topology,
            typename RRTStarVisitor,
            typename PositionMap,
            typename CostMap,
            typename WeightMap,
            typename RandomSampler,
            typename DistanceMetric,
            typename NcSelector>
  inline typename boost::enable_if< boost::is_directed_graph<Graph> >::type
    generate_rrt_star(Graph& g,
                      const Topology& space,
                      RRTStarVisitor vis,
                      PositionMap position,
                      CostMap cost,
                      WeightMap weight,
                      RandomSampler get_sample,
                      DistanceMetric distance,
                      const NcSelector& select_neighborhood,
                      unsigned int max_vertex_count) {
    BOOST_CONCEPT_ASSERT((RRTStarVisitorConcept<RRTStarVisitor,Graph,PositionMap>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::DistanceMetricConcept<DistanceMetric,Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    typedef typename Graph::vertex_bundled VertexProp;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::back_inserter;
    
    typedef boost::composite_property_map< 
      PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
    GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));

    if(num_vertices(g) == 0) {
      PositionValue p = get_sample(space);
      while(!vis.is_position_free(p))
        p = get_sample(space);
      VertexProp up;
      put(position, up, p);
      put(cost, up, 0.0);
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up),g);
#else
      Vertex u = add_vertex(up,g);
#endif
      vis.vertex_added(u, g);
    };

    while((num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
      
      PositionValue p_new = get_sample(space);
      
      std::vector<Vertex> Pred, Succ; 
      select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g, space, g_position);
      
      Vertex x_near;
      if( ! detail::expand_to_nearest(x_near, p_new, Pred, g, vis) )
        continue;
      
      Pred.clear(); Succ.clear();
      select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g, space, g_position);
      
      VertexProp xp;
      put(position, xp, p_new);
      put(cost, xp, std::numeric_limits<double>::infinity());
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex x_new = add_vertex(std::move(xp),g);
#else
      Vertex x_new = add_vertex(xp,g);
#endif
      vis.vertex_added(x_new,g);
      
      // Choose Parent
      double d_near = distance(get(position, g[x_near]), p_new, space);
      double c_near = get(cost, g[x_near]) + d_near;
      for(typename std::vector<Vertex>::const_iterator it = Pred.begin(); it != Pred.end(); ++it) {
        double c_temp, d_temp;
        if((*it != x_near) && 
           (vis.can_be_connected(*it, x_new, g)) &&
           ((c_temp = get(cost, g[*it]) + (d_temp = distance(get(position, g[*it]), p_new, space))) < c_near)) {
          x_near = *it;
          c_near = c_temp;
          d_near = d_temp;
        };
      };
      
      EdgeProp eprop;
      put(weight, eprop, d_near);
#ifdef RK_ENABLE_CXX0X_FEATURES
      std::pair<Edge, bool> ep = add_edge(x_near, x_new, std::move(eprop), g);
#else
      std::pair<Edge, bool> ep = add_edge(x_near, x_new, eprop, g);
#endif
      if(!ep.second)
        continue;
      put(cost, g[x_new], c_near);
      vis.edge_added(ep.first, g);
      
      detail::rewire_tree_neighborhood(x_new, Succ, g, space, vis, position, cost, weight, distance); 
      
    };

  };
  
  

};

};


#endif

