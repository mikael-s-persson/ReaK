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
	    typename JoiningVertexFoundCallback,
	    typename KeepGoingQuery,
	    typename CanConnectQuery>
  struct composite_rrt_star_visitor : composite_rrt_visitor<VertexAddedCallback, 
                                                            EdgeAddedCallback, 
						            SteerFunction,
						            IsFreeQuery,
						            JoiningVertexFoundCallback,
						            KeepGoingQuery> {
    CanConnectQuery can_be_connected;
    composite_rrt_star_visitor(VertexAddedCallback aVertexAdded,
                               EdgeAddedCallback aEdgeAdded,
			       SteerFunction aSteerTowardsPosition,
                               IsFreeQuery aIsFree,
                               JoiningVertexFoundCallback aJoiningVertexFound,
                               KeepGoingQuery aKeepGoing,
                               CanConnectQuery aCanBeConnected) :
                               vertex_added(aVertexAdded), 
                               edge_added(aEdgeAdded), 
                               steer_towards_position(aSteerTowardsPosition), 
                               is_position_free(aIsFree), 
                               joining_vertex_found(aJoiningVertexFound), 
                               keep_going(aKeepGoing),
                               can_be_connected(aCanBeConnected) { };
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
	    typename JoiningVertexFoundCallback,
	    typename KeepGoingQuery,
	    typename CanConnectQuery>
  inline composite_rrt_star_visitor<VertexAddedCallback, 
                                    EdgeAddedCallback, 
			            SteerFunction, 
			            IsFreeQuery, 
			            JoiningVertexFoundCallback, 
			            KeepGoingQuery,
			            CanConnectQuery>
    make_composite_rrt_star_visitor(VertexAddedCallback aVertexAdded,
                                    EdgeAddedCallback aEdgeAdded,
			            SteerFunction aSteerTowardsPosition,
                                    IsFreeQuery aIsFree,
                                    JoiningVertexFoundCallback aJoiningVertexFound,
                                    KeepGoingQuery aKeepGoing,
                                    CanConnectQuery aCanBeConnected) {
    return composite_rrt_star_visitor<VertexAddedCallback, 
                                      EdgeAddedCallback, 
				      SteerFunction, 
				      IsFreeQuery, 
				      JoiningVertexFoundCallback, 
				      KeepGoingQuery,
			              CanConnectQuery>(aVertexAdded,
						       aEdgeAdded,
			                               aSteerTowardsPosition,
			                               aIsFree,
			                               aJoiningVertexFound,
			                               aKeepGoing,
			                               aCanBeConnected);
  };





  /**
   * This function template is the unidirectional version of the RRT algorithm (refer to rr_tree.hpp dox).
   * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
   * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
   * \tparam RRTVisitor An RRT visitor type that implements the customizations to this RRT algorithm, should model the RRTVisitorConcept.
   * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
   * \tparam CostMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or goal).
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
   * \tparam DistanceMetric This is a distance-metric over the topology (see pp::DistanceMetricConcept);
   * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
   * \param g A mutable graph that should initially store the starting and goal 
   *        vertex and will store the generated graph once the algorithm has finished.
   * \param space A topology (as defined by the Boost Graph Library). Note 
   *        that it is not required to generate only random points in 
   *        the free-space.
   * \param vis A RRT visitor implementing the RRTVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param cost The property-map which stores the estimated cost-to-go of each vertex to the start (or goal).
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
   * \param distance A distance metric object to compute the distance between two positions.
   * \param select_neighborhood A callable object (functor) which can perform a 
   *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
   * \param max_vertex_count The maximum number of vertices beyond which the algorithm 
   *        should stop regardless of whether the resulting tree is satisfactory or not.
   * 
   */
  template <typename Graph,
	    typename Topology,
	    typename RRTVisitor,
	    typename PositionMap,
	    typename CostMap,
	    typename WeightMap,
	    typename DistanceMetric,
	    typename NcSelector>
  inline void generate_rrt_star(Graph& g,
				const Topology& space,
				RRTStarVisitor vis,
				PositionMap position,
				CostMap cost,
				WeightMap weight,
				DistanceMetric distance,
				const NcSelector& select_neighborhood,
				unsigned int max_vertex_count) {
    BOOST_CONCEPT_ASSERT((RRTVisitorConcept<RRTVisitor,Graph,PositionMap>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::LieGroupConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::DistanceMetricConcept<DistanceMetric,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    using std::back_inserter;

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
      std::vector<Vertex> Nc; 
      select_neighborhood(p_rnd, back_inserter(Nc), g, free_space, position);
      if(Nc.empty())
	continue;
      PositionValue p_new; bool expand_succeeded = false;
      std::vector<Vertex>::iterator it = Nc.begin();
      while((!expand_succeeded) && (it != Nc.end())) {
        boost::tie(p_new, expand_succeeded) = detail::expand_rrt_vertex(g, space, vis, position, *it, p_rnd);
	++it;
      };
      if(it == Nc.end())
	continue;
      Vertex x_near = *it;
      
      Nc.clear();
      select_neighborhood(p_new, back_inserter(Nc), g, free_space, position);
      
      Vertex x_new = add_vertex(g);
      put(position, x_new, p_new);
      vis.vertex_added(x_new,g);
        
      for(it = Nc.begin(); it != Nc.end(); ++it) {
	if((*it == x_near) || (vis.can_be_connected(*it, x_new, g))) {
	  std::pair<Edge, bool> ep = add_edge(*it, x_new, g);
          if(ep.second)
            vis.edge_added(ep.first, g);
        };
      };
    };

  };


};

};


#endif

