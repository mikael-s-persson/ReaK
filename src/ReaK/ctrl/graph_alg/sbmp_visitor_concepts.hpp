/**
 * \file sbmp_visitor_concepts.hpp
 * 
 * This library contains the 
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

#ifndef REAK_SBMP_VISITOR_CONCEPTS_HPP
#define REAK_SBMP_VISITOR_CONCEPTS_HPP

#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include <ReaK/ctrl/path_planning/metric_space_concept.hpp>


namespace ReaK {
  
namespace graph {
  
  

/**
  * This concept class defines what is required of a class to serve as a visitor to 
  * a Sampling-based Motion-planning algorithm.
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not yet connected.
  * 
  * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created vertex and its neighbor in the graph (g).
  * 
  * bool b = vis.keep_going();  This function is called at each attempt to expand the graph to verify that the user still wants more vertices to be generated in the motion-graph. 
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph>
struct SBMPVisitorConcept {
  Visitor vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename boost::graph_traits<Graph>::edge_descriptor e;
  bool b;
  
  BOOST_CONCEPT_USAGE(SBMPVisitorConcept) {
    vis.vertex_added(u, g); 
    vis.edge_added(e, g); 
    b = vis.keep_going(); 
  };
};

/**
  * This class is simply an archetype visitor for SBMPVisitorConcept.
  */
struct sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex,Graph&) const { };
  template <typename Edge, typename Graph>
  void edge_added(Edge,Graph&) const { };
  bool keep_going() const { return true; };
};

/**
  * This concept class defines what is required of a class to serve as a visitor to 
  * a Sampling-based Motion-planning algorithm that involves pruning of the motion-graph.
  * 
  * Required concepts:
  * 
  * the visitor should model SBMPVisitorConcept.
  * 
  * Valid expressions:
  * 
  * vis.vertex_to_be_removed(u, g);  This function is called just before a vertex is removed from the motion-graph (i.e., pruned).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph>
struct SBMPPruningVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  BOOST_CONCEPT_USAGE(SBMPPruningVisitorConcept) {
    this->vis.vertex_to_be_removed(this->u, this->g); 
  };
};

/**
  * This class is simply an archetype visitor for SBMPPruningVisitorConcept.
  */
struct sbmp_pruning_visitor_archetype : sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void vertex_to_be_removed(Vertex,Graph&) const { };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to 
  * a Sampling-based Motion-planning algorithm that involves finding joining-points of 
  * multiple motion-graphs (e.g., bi-directional algorithms).
  * 
  * Required concepts:
  * 
  * the visitor should model the SBMPVisitorConcept.
  * 
  * Valid expressions:
  * 
  * vis.joining_vertex_found(u1, u2, g1, g2);  This function is called by the algorithm when two vertices (u1,u2) of two graphs (g1,g2) is found that meet each other.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph>
struct SBMPJoiningVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  BOOST_CONCEPT_USAGE(SBMPJoiningVisitorConcept) {
    this->vis.joining_vertex_found(this->u, this->u, this->g, this->g); 
  };
};

/**
  * This class is simply an archetype visitor for SBMPJoiningVisitorConcept.
  */
struct sbmp_joining_visitor_archetype : sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void joining_vertex_found(Vertex,Vertex,Graph&,Graph&) const { };
};


/**
  * This concept class defines what is required of a class to serve as a node-creator for a motion-planning 
  * algorithm. This is typical used internally (in intermediate visitors) for most motion-planners.
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * u = conn_vis.create_vertex(p, g);  This function is called to request the visitor to create a new vertex (u) in the graph (g), with position (p).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct NodeCreatorVisitorConcept {
  Visitor vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  
  BOOST_CONCEPT_USAGE(NodeCreatorVisitorConcept) {
    u = vis.create_vertex(p, g);
  };
};

/**
  * This class is simply an archetype visitor for NodeCreatorVisitorConcept.
  */
struct node_creator_visitor_archetype {
  template <typename Position, typename Graph>
  typename boost::graph_traits<Graph>::vertex_descriptor create_vertex(const Position&, Graph&) const { 
    return boost::graph_traits<Graph>::null_vertex();
  };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that pulls nodes towards some sample-point (i.e., attempt to steer, at least, partially
  * towards a point). This is typical of RRT-style motion-planners.
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.steer_towards_position(p,u,g);  This function is called to attempt to steer from vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether any significant motion occurred (collision-free).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct NodePullingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_USAGE(NodePullingVisitorConcept) {
    boost::tie(p, b, ep) = vis.steer_towards_position(p, u, g);
  };
};

/**
  * This class is simply an archetype visitor for NodePullingVisitorConcept.
  */
struct node_pulling_visitor_archetype {
  template <typename Position, typename Vertex, typename Graph>
  boost::tuple<Position, bool, typename Graph::edge_bundled> steer_towards_position(const Position& p, Vertex, Graph&) const { 
    typedef typename Graph::edge_bundled EdgeProp;
    return boost::tuple<Position, bool, EdgeProp>(p, false, EdgeProp());
  };
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that pulls nodes backwards towards some sample-point (i.e., attempt to steer (reverse-time), 
  * at least, partially towards a point). This is typical of RRT-style motion-planners.
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.steer_back_to_position(p,u,g);  This function is called to attempt to steer backwards in time from vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether any significant motion occurred (collision-free).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct NodeBackPullingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_USAGE(NodeBackPullingVisitorConcept) {
    boost::tie(p, b, ep) = vis.steer_back_to_position(p, u, g);
  };
};

/**
  * This class is simply an archetype visitor for NodeBackPullingVisitorConcept.
  */
struct node_back_pulling_visitor_archetype {
  template <typename Position, typename Vertex, typename Graph>
  boost::tuple<Position, bool, typename Graph::edge_bundled> steer_back_to_position(const Position& p, Vertex, Graph&) const { 
    typedef typename Graph::edge_bundled EdgeProp;
    return boost::tuple<Position, bool, EdgeProp>(p, false, EdgeProp());
  };
};



/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that expands existing nodes towards some unspecified direction (i.e., attempt to steer, 
  * at least, partially towards some random or guided direction). This is typical of PRM-style 
  * motion-planners and of what one could call "pushed" expansion algorithms (SBA*, FADPRM, etc.).
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph, and with an edge-property object associated with a new edge.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct NodePushingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_USAGE(NodePushingVisitorConcept) {
    boost::tie(p, b, ep) = vis.random_walk(u, g);
  };
};

/**
  * This class is simply an archetype visitor for NodePushingVisitorConcept.
  */
template <typename Topology>
struct node_pushing_visitor_archetype {
  typedef typename ReaK::pp::topology_traits<Topology>::point_type PointType;
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, typename Graph::edge_bundled> random_walk(Vertex, Graph&) const { 
    typedef typename Graph::edge_bundled EdgeProp;
    return boost::tuple<PointType, bool, EdgeProp>(PointType(), false, EdgeProp());
  };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that retracts existing nodes backwards towards some unspecified direction (i.e., attempt to steer, 
  * at least, partially, backwards towards some random or guided direction). This is typical of PRM-style 
  * motion-planners and of what one could call "pushed" expansion algorithms (SBA*, FADPRM, etc.).
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.random_back_walk(u, g);  This function is called to perform the retraction of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph, and with an edge-property object associated with a new edge.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct NodeBackPushingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_USAGE(NodeBackPushingVisitorConcept) {
    boost::tie(p, b, ep) = vis.random_back_walk(u, g);
  };
};

/**
  * This class is simply an archetype visitor for NodeBackPushingVisitorConcept.
  */
template <typename Topology>
struct node_back_pushing_visitor_archetype {
  typedef typename ReaK::pp::topology_traits<Topology>::point_type PointType;
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, typename Graph::edge_bundled> random_back_walk(Vertex, Graph&) const { 
    typedef typename Graph::edge_bundled EdgeProp;
    return boost::tuple<PointType, bool, EdgeProp>(PointType(), false, EdgeProp());
  };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that attempts re-connections of existing nodes (i.e., attempt to steer, 
  * fully, from one node to another). This is typical of PRM-style motion-planners and of 
  * probabilistically optimal connection strategies which must involve re-wirings (e.g., RRT*, SBA*, etc.).
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * tie(b, ep) = vis.can_be_connected(u,v,g);  This function is called to attempt to steer from vertex u to vertex v, it returns true if a local path exists and is collision-free, and it also returns the edge-property of the edge that could connect those two vertices.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph>
struct NodeReConnectVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u, v;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_USAGE(NodeReConnectVisitorConcept) {
    boost::tie(b, ep) = vis.can_be_connected(u, v, g);
  };
};

/**
  * This class is simply an archetype visitor for NodeReConnectVisitorConcept.
  */
struct node_reconnect_visitor_archetype {
  template <typename Vertex, typename Graph>
  std::pair<bool, typename Graph::edge_bundled> can_be_connected(Vertex, Vertex, Graph&) const { 
    typedef typename Graph::edge_bundled EdgeProp;
    return std::pair<bool, EdgeProp>(false, EdgeProp());
  };
};



/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that requires the tracking of the neighborhood of the nodes (e.g., keeping tabs of the
  * density of the neighborhood). This is typical of PRM-style motion-planners and of 
  * other "pushed" expansion algorithms (e.g., SBA*, FADPRM, etc.).
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * vis.travel_explored(u, v, g);  This function is called whenever a source-destination pair of vertices have been matched for a potential edge between them.
  * 
  * vis.travel_succeeded(u, v, g);  This function is called whenever a source-destination pair of vertices led to a successful travel.
  * 
  * vis.travel_failed(u, v, g);  This function is called whenever a source-destination pair of vertices led to a failure to travel (e.g., a collision occurred).
  * 
  * vis.affected_vertex(u,g);  This function is called to notify the visitor that something about the vertex's neighborhood might have changed (possibly nothing). For example, this would be the place where to re-compute a density metric.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph>
struct NeighborhoodTrackingVisitorConcept {
  Visitor vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u, v;
  
  BOOST_CONCEPT_USAGE(NeighborhoodTrackingVisitorConcept) {
    vis.travel_explored(u, v, g);
    vis.travel_succeeded(u, v, g);
    vis.travel_failed(u, v, g);
    vis.affected_vertex(u, g);
  };
};

/**
  * This class is simply an archetype visitor for NeighborhoodTrackingVisitorConcept.
  */
struct neighborhood_tracking_visitor_archetype {
  template <typename Vertex, typename Graph>
  void travel_explored(Vertex, Vertex, Graph&) const { };
  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex, Vertex, Graph&) const { };
  template <typename Vertex, typename Graph>
  void travel_failed(Vertex, Vertex, Graph&) const { };
  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex, Graph&) const { };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to check collision of 
  * a point during a motion-planning algorithm.
  * 
  * Valid expressions:
  * 
  * bool b = vis.is_position_free(p);  This function is called to query whether a particular configuration (position, p) is free.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Topology>
struct CollisionCheckingVisitorConcept {
  Visitor vis;
  bool b;
  typename ReaK::pp::topology_traits<Topology>::point_type p;
  
  BOOST_CONCEPT_USAGE(CollisionCheckingVisitorConcept) {
    b = vis.is_position_free(p);
  };
};

/**
  * This class is simply an archetype visitor for CollisionCheckingVisitorConcept.
  */
struct collision_checking_visitor_archetype {
  template <typename Position>
  bool is_position_free(const Position&) const { return true; };
};



/**
  * This concept class defines the valid expressions required of a class to be used as a 
  * node-exploring visitor. This is used in algorithms that look through a set of vertices 
  * of the motion-graph (either to expand from them or simply doing a shortest-path algorithm).
  * 
  * Required concepts:
  * 
  * Valid expressions:
  * 
  * vis.initialize_vertex(u, g);  A function that gets called whenever a vertex (u) is first initialized before the search.
  * 
  * vis.discover_vertex(u, g);  A function that gets called whenever a vertex (u) is added to the OPEN set (or updated in the OPEN set).
  * 
  * vis.examine_vertex(u, g);  A function that gets called whenever a vertex (u) is taken out of the OPEN set to be examined, this is called before it gets expanded.
  * 
  * vis.examine_edge(e, g);  A function that is called whenever an edge (e) is examined (usually after its source vertex has been examined).
  * 
  * b = vis.has_search_potential(u, g);  This function is called to check if the vertex is deemed minimally useful to examine by some measure.
  * 
  * b = vis.should_close(u, g);  This function is called to check if a vertex should be closed according to some measure, i.e., it should no longer be elligible to be examined in the future.
  * 
  * \tparam Visitor The visitor class to be tested for modeling this visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  */
template <typename Visitor, typename Graph>
struct NodeExploringVisitorConcept {
  
  BOOST_CONCEPT_USAGE(NodeExploringVisitorConcept)
  {
    vis.initialize_vertex(u, g);
    vis.discover_vertex(u, g);
    vis.examine_vertex(u, g);
    vis.examine_edge(e, g);
    b = vis.has_search_potential(u, g); 
    b = vis.should_close(u, g); 
  }
  Visitor vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename boost::graph_traits<Graph>::edge_descriptor e;
  bool b;
};

/**
  * This class is simply an archetype visitor for NodeExploringVisitorConcept.
  */
struct node_exploring_visitor_archetype { 
  
  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex, const Graph&) const { };
  
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex, const Graph&) const { };
  
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex, const Graph&) const { };
  
  template <typename Edge, typename Graph>
  void examine_edge(Edge, const Graph&) const { };
  
  template <typename Vertex, typename Graph>
  bool has_search_potential(Vertex, const Graph&) const { return true; };
  
  template <typename Vertex, typename Graph>
  bool should_close(Vertex, const Graph&) const { return false; };
};




/**
  * This concept class defines what is required of a class to serve as a visitor to the RRT algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model SBMPVisitorConcept, NodePullingVisitorConcept, and CollisionCheckingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct RRTVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_ASSERT((NodePullingVisitorConcept<Visitor, Graph, Topology>));
  BOOST_CONCEPT_ASSERT((CollisionCheckingVisitorConcept<Visitor, Topology>));
  
  BOOST_CONCEPT_USAGE(RRTVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for the RRT algorithm.
  */
struct rrt_visitor_archetype : 
  sbmp_visitor_archetype, 
  node_pulling_visitor_archetype, 
  collision_checking_visitor_archetype { };


/**
  * This concept class defines what is required of a class to serve as a visitor to the Bi-RRT algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model RRTVisitorConcept and SBMPJoiningVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct BiRRTVisitorConcept : 
  RRTVisitorConcept<Visitor, Graph, Topology>, 
  SBMPJoiningVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_USAGE(BiRRTVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for the Bi-RRT algorithm.
  */
struct birrt_visitor_archetype : 
  rrt_visitor_archetype, 
  sbmp_joining_visitor_archetype { };


/**
  * This concept class defines what is required of a class to serve as a visitor to the PRM algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model SBMPVisitorConcept, NodePushingVisitorConcept, 
  * NodeReConnectVisitorConcept, NeighborhoodTrackingVisitorConcept, and CollisionCheckingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct PRMVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_ASSERT((NodePushingVisitorConcept<Visitor, Graph, Topology>));
  BOOST_CONCEPT_ASSERT((NodeReConnectVisitorConcept<Visitor, Graph>));
  BOOST_CONCEPT_ASSERT((NeighborhoodTrackingVisitorConcept<Visitor, Graph>));
  BOOST_CONCEPT_ASSERT((CollisionCheckingVisitorConcept<Visitor, Topology>));
  
  BOOST_CONCEPT_USAGE(PRMVisitorConcept) { };
  
};

/**
  * This class is simply an archetype visitor for the PRM algorithm.
  */
template <typename Topology>
struct prm_visitor_archetype : 
  sbmp_visitor_archetype, 
  node_pushing_visitor_archetype<Topology>, 
  node_reconnect_visitor_archetype, 
  neighborhood_tracking_visitor_archetype, 
  collision_checking_visitor_archetype { };


/**
  * This concept class defines what is required of a class to serve as a visitor to the RRG algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model RRTVisitorConcept, and NodeReConnectVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct RRGVisitorConcept : RRTVisitorConcept<Visitor,Graph,Topology> {
  
  BOOST_CONCEPT_ASSERT((NodeReConnectVisitorConcept<Visitor,Graph>));
  
  BOOST_CONCEPT_USAGE(RRGVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for RRGVisitorConcept.
  */
struct rrg_visitor_archetype : 
  rrt_visitor_archetype, 
  node_reconnect_visitor_archetype { };

/**
  * This concept class defines what is required of a class to serve as a visitor to the bi-directional RRG algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model RRGVisitorConcept, and NodeBackPullingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename Topology>
struct RRGBidirVisitorConcept : RRGVisitorConcept<Visitor,Graph,Topology> {
  
  BOOST_CONCEPT_ASSERT((NodeBackPullingVisitorConcept<Visitor,Graph,Topology>));
  
  BOOST_CONCEPT_USAGE(RRGBidirVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for RRGBidirVisitorConcept.
  */
struct rrg_bidir_visitor_archetype : 
  rrg_visitor_archetype, 
  node_back_pulling_visitor_archetype { };



/**
  * This concept class defines the valid expressions required of a class to be used as a visitor 
  * class for the connection strategies. 
  * 
  * Required concepts:
  * 
  * The visitor class should model SBMPVisitorConcept, NodeCreatorVisitorConcept, NodeReConnectVisitorConcept,
  * and NeighborhoodTrackingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be tested for modeling this visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename ConnectorVisitor, typename Graph, typename Topology>
struct MotionGraphConnectorVisitorConcept {
  
  BOOST_CONCEPT_ASSERT((SBMPVisitorConcept<ConnectorVisitor, Graph>));
  BOOST_CONCEPT_ASSERT((NodeCreatorVisitorConcept<ConnectorVisitor, Graph, Topology>));
  BOOST_CONCEPT_ASSERT((NodeReConnectVisitorConcept<ConnectorVisitor, Graph>));
  BOOST_CONCEPT_ASSERT((NeighborhoodTrackingVisitorConcept<ConnectorVisitor, Graph>));
  
  BOOST_CONCEPT_USAGE(MotionGraphConnectorVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for MotionGraphConnectorVisitorConcept.
  */
struct mg_connector_visitor_archetype : 
  sbmp_visitor_archetype, 
  node_creator_visitor_archetype,
  node_reconnect_visitor_archetype,
  neighborhood_tracking_visitor_archetype { };



/**
  * This concept class defines the valid expressions required of a class to be used as an anytime-heuristic 
  * visitor.
  * 
  * Valid expressions:
  * 
  * d = vis.adjust_relaxation(d, g);  This function should return a new value for the relaxation factor used in the anytime-heuristic algorithm.
  * 
  * \tparam Visitor The visitor class to be tested for modeling this visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  */
template <typename Visitor, typename Graph>
struct AnytimeHeuristicVisitorConcept {
  BOOST_CONCEPT_USAGE(AnytimeHeuristicVisitorConcept)
  {
    d = vis.adjust_relaxation(d, g); 
  }
  Visitor vis;
  Graph g;
  double d;
};

/**
  * This class is simply an archetype visitor for AnytimeHeuristicVisitorConcept.
  */
struct anytime_heuristic_visitor_archetype {
  template <typename Graph>
  double adjust_relaxation(double d, Graph&) const { 
    return d;
  };
};


};

};


#endif

