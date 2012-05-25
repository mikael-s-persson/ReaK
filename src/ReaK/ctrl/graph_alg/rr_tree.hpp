/**
 * \file rr_tree.hpp
 * 
 * This library contains the Rapidly-Exploring Random Tree generation algorithm.
 * This is a method to create a random tree that will span over a non-convex space
 * as rapidly as possible. The method relies on a simple randomized insertion algorithm.
 * At each step, a random point is picked from the underlying topology (i.e. configuration 
 * space in path-planning terms). Then, the point in the current graph that is nearest 
 * to the random point is picked for expansion. Finally, edges (of a maximum length) are 
 * added to the vertex of the graph towards the random point while it is still possible to
 * add such an edge without leaving the free space (the part of the configuration space which 
 * is not occupied by an obstacle). The algorithm will stop when either the number of vertices
 * in the tree has reached a maximum or when the user callback signals the stop.
 * 
 * This library also provides the bidirectional version of the RRT algorithm. In this version,
 * two trees are generated. Typically, one tree is initialized with the starting vertex and 
 * the other is initialized with the goal vertex. The algorithm works to try and join the 
 * two graphs as quickly as possible and with the most direct path. The algorithm alternates
 * between the two graphs. It first uses the normal procedure (as in the unidirectional variant)
 * to try and add a vertex to one graph. If it does not succeed (i.e. there was no free-space in 
 * the expanded direction), it will move to the other graph and try to expand it. If it does succeed,
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
 * \date February 2011
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

#ifndef REAK_RR_TREE_HPP
#define REAK_RR_TREE_HPP

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


namespace ReaK {
  
namespace graph {

  /**
   * This concept class defines what is required of a class to serve as a visitor to the RRT algorithm.
   * 
   * Required concepts:
   * 
   * the visitor should model the boost::CopyConstructibleConcept.
   * 
   * Valid expressions:
   * 
   * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not yet connected.
   * 
   * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created vertex and its nearest neighbor in the graph (g).
   *
   * tie(p,b) = vis.steer_towards_position(p,u,g);  This function is called to attempt to steer from vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether any significant motion occurred (collision-free).
   * 
   * bool b = vis.is_position_free(p);  This function is called to query whether a particular configuration (position, p) is free.
   * 
   * vis.joining_vertex_found(u, g);  This function is called by the bidirectional RRT algorithm when a vertex (u) of graph (g) is found that can reach the other graph (not the one passed as 2nd parameter).
   * 
   * bool b = vis.keep_going();  This function is called at each attempt to expand the graph to verify that the user still wants more vertices to be generated in the tree. 
   *
   * \tparam Visitor The visitor class to be checked for modeling this concept.
   * \tparam Graph The graph on which the visitor class is required to work with.
   * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
   */
  template <typename Visitor, typename Graph, typename PositionMap>
  struct RRTVisitorConcept {
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
    typename boost::property_traits<PositionMap>::value_type p;
    bool b;
    
    BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
    
    BOOST_CONCEPT_USAGE(RRTVisitorConcept) {
      vis.vertex_added(u, g); 
      vis.edge_added(e, g); 
      boost::tie(p,b) = vis.steer_towards_position(p,u,g);
      b = vis.is_position_free(p);
      vis.joining_vertex_found(u, g); 
      b = vis.keep_going(); 
    };
  };

  /**
   * This class is simply a "null" visitor for the RRT algorithm. It is null in the sense that it
   * will do nothing (except return true on the is_position_free(p) and keep_going() callbacks).
   */
  struct default_rrt_visitor {
    template <typename Vertex, typename Graph>
    void vertex_added(Vertex,Graph&) { };
    template <typename Edge, typename Graph>
    void edge_added(Edge,Graph&) { };
    template <typename PositionValue, typename Vertex, typename Graph>
    std::pair<PositionValue,bool> steer_towards_position(const PositionValue& p,Vertex,Graph&) { return std::make_pair(p,true); };
    template <typename PositionValue>
    bool is_position_free(const PositionValue&) { return true; };
    template <typename Vertex, typename Graph>
    void joining_vertex_found(Vertex,Graph&) { };
    bool keep_going() { return true; };
  };

  /**
   * This class is a composite visitor class template. It can be used to glue together a function pointer (or functor)
   * for each of the functions of the RRTVisitorConcept so that it can be used as a light-weight,
   * copyable visitor object for the RRT algorithm (it is especially recommend to use the 
   * make_composite_rrt_visitor function template).
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename SteerFunction,
	    typename IsFreeQuery,
	    typename JoiningVertexFoundCallback,
	    typename KeepGoingQuery>
  struct composite_rrt_visitor {
    VertexAddedCallback vertex_added;
    EdgeAddedCallback edge_added;
    SteerFunction steer_towards_position;
    IsFreeQuery is_position_free;
    JoiningVertexFoundCallback joining_vertex_found;
    KeepGoingQuery keep_going;
    composite_rrt_visitor(VertexAddedCallback aVertexAdded,
                          EdgeAddedCallback aEdgeAdded,
			  SteerFunction aSteerTowardsPosition,
                          IsFreeQuery aIsFree,
                          JoiningVertexFoundCallback aJoiningVertexFound,
                          KeepGoingQuery aKeepGoing) :
                          vertex_added(aVertexAdded), edge_added(aEdgeAdded), steer_towards_position(aSteerTowardsPosition), is_position_free(aIsFree), joining_vertex_found(aJoiningVertexFound), keep_going(aKeepGoing) { };
  };

  /**
   * This is a function template that is used to create an object of a class of the composite_rrt_visitor
   * class template. This is particularly convenient to avoid explicitely providing the list of template
   * arguments and let the compiler resolved them from the function parameter types.
   */
  template <typename VertexAddedCallback,
            typename EdgeAddedCallback,
	    typename SteerFunction,
	    typename IsFreeQuery,
	    typename JoiningVertexFoundCallback,
	    typename KeepGoingQuery>
  inline composite_rrt_visitor<VertexAddedCallback, EdgeAddedCallback, SteerFunction, IsFreeQuery, JoiningVertexFoundCallback, KeepGoingQuery>
    make_composite_rrt_visitor(VertexAddedCallback aVertexAdded,
                               EdgeAddedCallback aEdgeAdded,
			       SteerFunction aSteerTowardsPosition,
                               IsFreeQuery aIsFree,
                               JoiningVertexFoundCallback aJoiningVertexFound,
                               KeepGoingQuery aKeepGoing) {
    return composite_rrt_visitor<VertexAddedCallback, EdgeAddedCallback, SteerFunction, IsFreeQuery, JoiningVertexFoundCallback, KeepGoingQuery>(aVertexAdded,aEdgeAdded,aSteerTowardsPosition,aIsFree,aJoiningVertexFound,aKeepGoing);
  };




namespace detail {

  template <typename Graph,
	    typename Topology,
	    typename RRTVisitor,
	    typename PositionMap>
  inline std::pair< typename boost::graph_traits<Graph>::vertex_descriptor, bool>
    expand_rrt_vertex(Graph& g, const Topology& space, RRTVisitor vis, PositionMap position,
                      typename boost::graph_traits<Graph>::vertex_descriptor u,
                      const typename boost::property_traits<PositionMap>::value_type& p_target) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    
    bool reached_new = false;
    PositionValue p_v;
    boost::tie(p_v, reached_new) = vis.steer_towards_position(p_target,u,g);
    if(reached_new) {  // i.e., a new position was reached, collision-free.
      Vertex v = add_vertex(g);
      put(position, v, p_v);
      vis.vertex_added(v,g);
      std::pair<Edge, bool> ep = add_edge(u,v,g);
      if(ep.second)
        vis.edge_added(ep.first, g);
      u = v;
    };
    return make_pair(u,reached_new);
  };
  
}; //namespace detail



  /**
   * This function template is the unidirectional version of the RRT algorithm (refer to rr_tree.hpp dox).
   * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
   * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
   * \tparam RRTVisitor An RRT visitor type that implements the customizations to this RRT algorithm, should model the RRTVisitorConcept.
   * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
   * \tparam NNFinder A functor type which can perform a nearest-neighbor search of a point to a graph in the topology (see topological_search.hpp).
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated tree once the algorithm has finished.
   * \param space A topology (as defined by the Boost Graph Library). Note 
   *        that it is not required to generate only random points in 
   *        the free-space.
   * \param vis A RRT visitor implementing the RRTVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param find_nearest_neighbor A callable object (functor) which can perform a 
   *        nearest neighbor search of a point to a graph in the topology. (see topological_search.hpp)
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting tree is satisfactory or not.
   * \param max_edge_distance The maximum length (w.r.t. the topology) of an edge.
   * \param min_edge_distance The minimum length (w.r.t. the topology) of an edge. If a free-space 
   *        edge cannot be made longer than this number, it will not be added. This 
   *        parameter also serves as the minimum resolution of the free-query (i.e. 
   *        the collision detection).
   * 
   */
  template <typename Graph,
	    typename Topology,
	    typename RRTVisitor,
	    typename PositionMap,
	    typename NNFinder>
  inline void generate_rrt(Graph& g,
			   const Topology& space,
			   RRTVisitor vis,
			   PositionMap position,
			   NNFinder find_nearest_neighbor,
			   unsigned int max_vertex_count) {
    BOOST_CONCEPT_ASSERT((RRTVisitorConcept<RRTVisitor,Graph,PositionMap>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::LieGroupConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    if(num_vertices(g) == 0) {
      Vertex u = add_vertex(g);
      PositionValue p = get(ReaK::pp::random_sampler, space)(space);
      while(!vis.is_position_free(p))
	p = get(ReaK::pp::random_sampler, space)(space);
      put(position, u, p);
      vis.vertex_added(u, g);
    };

    while((num_vertices(g) < max_vertex_count) && (vis.keep_going())) {
      PositionValue p_rnd = get(ReaK::pp::random_sampler, space)(space);
      Vertex u = find_nearest_neighbor(p_rnd,g,space,position);
      detail::expand_rrt_vertex(g, space, vis, position, u, p_rnd);
    };

  };


  /**
   * This function is the bidirectional version of the RRT algorithm (refer to rr_tree.hpp dox).
   * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
   * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
   * \tparam RRTVisitor An RRT visitor type that implements the customizations to this RRT algorithm, should model the RRTVisitorConcept.
   * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
   * \tparam NNFinder A functor type which can perform a nearest-neighbor search of a point to a graph in the topology (see topological_search.hpp).
   * \param g1 A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated tree once the algorithm has finished.
   * \param g2 A mutable graph that should initially store the goal 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated tree once the algorithm has finished.
   * \param space A topology (as defined by the Boost Graph Library). Note 
   *        that it is not required to generate only random points in 
   *        the free-space.
   * \param vis A RRT visitor implementing the RRTVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position1 A mapping for the first graph that implements the MutablePropertyMap Concept. 
   *        Also, the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param position2 A mapping for the second graph that implements the MutablePropertyMap Concept. 
   *        Also, the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param find_nearest_neighbor A callable object (functor) which can perform a 
   *        nearest neighbor search of a point to a graph in the 
   *        topology. (see topological_search.hpp)
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting tree is satisfactory or not.
   * \param max_edge_distance The maximum length (w.r.t. the topology) of an edge.
   * \param min_edge_distance The minimum length (w.r.t. the topology) of an edge. If a free-space 
   *        edge cannot be made longer than this number, it will not be added. This 
   *        parameter also serves as the minimum resolution of the free-query (i.e. 
   *        the collision detection).
   * 
   */
  template <typename Graph,
	    typename Topology,
	    typename RRTVisitor,
	    typename PositionMap,
	    typename NNFinder>
  inline void generate_bidirectional_rrt(Graph& g1, Graph& g2,
			                 const Topology& space,
                                         RRTVisitor vis,
			                 PositionMap position1, PositionMap position2,
			                 NNFinder find_nearest_neighbor,
			                 unsigned int max_vertex_count,
			                 double max_edge_distance, double min_edge_distance) {
    BOOST_CONCEPT_ASSERT((RRTVisitorConcept<RRTVisitor,Graph,PositionMap>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    PositionValue p_target1; std::pair<Vertex,bool> v_target1;
    PositionValue p_target2; std::pair<Vertex,bool> v_target2;

    if(num_vertices(g1) == 0) {
      Vertex u = add_vertex(g1);
      PositionValue p = get(ReaK::pp::random_sampler, space)(space);
      while(!vis.is_position_free(p))
	p = get(ReaK::pp::random_sampler, space)(space);
      put(position1, u, p);
      p_target2 = p;
      v_target2.first = u; v_target2.second = true;
      vis.vertex_added(u, g1);
    } else {
      Vertex u = vertex(0,g1);
      p_target2 = get(position1, u);
      v_target2.first = u; v_target2.second = true;
    };

    if(num_vertices(g2) == 0) {
      Vertex u = add_vertex(g2);
      PositionValue p = get(ReaK::pp::random_sampler, space)(space);
      while(!vis.is_position_free(p))
	p = get(ReaK::pp::random_sampler, space)(space);
      put(position2, u, p);
      p_target1 = p;
      v_target1.first = u; v_target1.second = true;
      vis.vertex_added(u, g2);
    } else {
      Vertex u = vertex(0,g2);
      p_target1 = get(position2, u);
      v_target1.first = u; v_target1.second = true;
    };

    while((num_vertices(g1) + num_vertices(g2) < max_vertex_count) && (vis.keep_going())) {
      //first, expand the first graph towards its target:
      Vertex u1 = find_nearest_neighbor(p_target1,g1,space,position1);
      std::pair< Vertex, bool> v1 =
        detail::expand_rrt_vertex(g1,space,vis,position1,
                                  u1,p_target1,max_edge_distance,min_edge_distance);
      if((v1.second) && (v_target1.second)) {
        //joining vertex has been reached!
        vis.joining_vertex_found(v1.first, g1);
        vis.joining_vertex_found(v_target1.first, g2);
        p_target2 = get(ReaK::pp::random_sampler, space)(space);
        v_target2.second = false;
      } else {
        if(v1.first == u1) { //we didn't move at all! Unsuccessful expansion.
          p_target2 = get(ReaK::pp::random_sampler, space)(space);
          v_target2.second = false;
        } else {
          p_target2 = get(position1, v1.first);
          v_target2.first = v1.first; v_target2.second = true;
        };
      };

      //then, expand the second graph towards its target:
      Vertex u2 = find_nearest_neighbor(p_target2,g2,space,position2);
      std::pair< Vertex, bool> v2 =
        detail::expand_rrt_vertex(g2,space,vis,position2,
                                  u2,p_target2,max_edge_distance,min_edge_distance);
      if((v2.second) && (v_target2.second)) {
        //joining vertex has been reached!
        vis.joining_vertex_found(v2.first, g2);
        vis.joining_vertex_found(v_target2.first, g1);
        p_target1 = get(ReaK::pp::random_sampler, space)(space);
        v_target1.second = false;
      } else {
        if(v2.first == u2) { //we didn't move at all! Unsuccessful expansion.
          p_target1 = get(ReaK::pp::random_sampler, space)(space);
          v_target1.second = false;
        } else {
          p_target1 = get(position2, v2.first);
          v_target1.first = v2.first; v_target1.second = true;
        };
      };
    };

  };



};

};


#endif

