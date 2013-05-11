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

#include <functional>
#include <vector>
#include <boost/limits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include "path_planning/metric_space_concept.hpp"


namespace ReaK {
  
namespace graph {
  
  

/**
  * This concept class defines what is required of a class to serve as a visitor to 
  * a Sampling-based Motion-planning algorithm.
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
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(SBMPVisitorConcept) {
    vis.vertex_added(u, g); 
    vis.edge_added(e, g); 
    b = vis.keep_going(); 
  };
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
  * This concept class defines what is required of a class to serve as a node-creator for a motion-planning 
  * algorithm. This is typical used internally (in intermediate visitors) for most motion-planners.
  * 
  * Required concepts:
  * 
  * the visitor should model the boost::CopyConstructibleConcept.
  * 
  * Valid expressions:
  * 
  * u = conn_vis.create_vertex(p, g);  This function is called to request the visitor to create a new vertex (u) in the graph (g), with position (p).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct NodeCreatorVisitorConcept {
  Visitor vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename boost::property_traits<PositionMap>::value_type p;
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(NodeCreatorVisitorConcept) {
    u = vis.create_vertex(p, g);
  };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning 
  * algorithm that pulls nodes towards some sample-point (i.e., attempt to steer, at least, partially
  * towards a point). This is typical of RRT-style motion-planners.
  * 
  * Required concepts:
  * 
  * the visitor should model the boost::CopyConstructibleConcept.
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.steer_towards_position(p,u,g);  This function is called to attempt to steer from vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether any significant motion occurred (collision-free).
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct NodePullingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename boost::property_traits<PositionMap>::value_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(NodePullingVisitorConcept) {
    boost::tie(p, b, ep) = vis.steer_towards_position(p, u, g);
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
  * the visitor should model the boost::CopyConstructibleConcept.
  * 
  * Valid expressions:
  * 
  * tie(p,b,ep) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph, and with an edge-property object associated with a new edge.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct NodePushingVisitorConcept {
  Visitor vis;
  Graph g;
  bool b;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  typename boost::property_traits<PositionMap>::value_type p;
  typename Graph::edge_bundled ep;
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(NodePushingVisitorConcept) {
    boost::tie(p, b, ep) = vis.random_walk(u, g);
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
  * the visitor should model the boost::CopyConstructibleConcept.
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
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(NodeReConnectVisitorConcept) {
    boost::tie(b, ep) = vis.can_be_connected(u, v, g);
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
  * the visitor should model the boost::CopyConstructibleConcept.
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
  
  BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
  
  BOOST_CONCEPT_USAGE(NeighborhoodTrackingVisitorConcept) {
    vis.travel_explored(u, v, g);
    vis.travel_succeeded(u, v, g);
    vis.travel_failed(u, v, g);
    vis.affected_vertex(u, g);
  };
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
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename PositionMap>
struct CollisionCheckVisitorConcept {
  Visitor vis;
  bool b;
  typename boost::property_traits<PositionMap>::value_type p;
  
  BOOST_CONCEPT_USAGE(CollisionCheckVisitorConcept) {
    b = vis.is_position_free(p);
  };
};


#if 0
/**
  * This concept class defines what is required of a class to serve as a visitor to the RRT algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model SBMPVisitorConcept, NodePullingVisitorConcept, and CollisionCheckVisitorConcept.
  * 
  * Valid expressions:
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct RRTVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_ASSERT((NodePullingVisitorConcept<Visitor, Graph, PositionMap>));
  BOOST_CONCEPT_ASSERT((CollisionCheckVisitorConcept<Visitor, PositionMap>));
  
  BOOST_CONCEPT_USAGE(RRTVisitorConcept) { };
};


/**
  * This concept class defines what is required of a class to serve as a visitor to the Bi-RRT algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model RRTVisitorConcept and SBMPJoiningVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct BiRRTVisitorConcept : 
  RRTVisitorConcept<Visitor, Graph, PositionMap>, 
  SBMPJoiningVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_USAGE(BiRRTVisitorConcept) { };
};
#endif


/**
  * This concept class defines what is required of a class to serve as a visitor to the PRM algorithm.
  * 
  * Required concepts:
  * 
  * the visitor should model SBMPVisitorConcept, NodePushingVisitorConcept, 
  * NodeReConnectVisitorConcept, NeighborhoodTrackingVisitorConcept, and CollisionCheckVisitorConcept.
  * 
  * Valid expressions:
  * 
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Graph The graph on which the visitor class is required to work with.
  * \tparam PositionMap The position property-map that provides the position descriptors with which the visitor class is required to work.
  */
template <typename Visitor, typename Graph, typename PositionMap>
struct PRMVisitorConcept : SBMPVisitorConcept<Visitor, Graph> {
  
  BOOST_CONCEPT_ASSERT((NodePushingVisitorConcept<Visitor, Graph, PositionMap>));
  BOOST_CONCEPT_ASSERT((NodeReConnectVisitorConcept<Visitor, Graph>));
  BOOST_CONCEPT_ASSERT((NeighborhoodTrackingVisitorConcept<Visitor, Graph>));
  BOOST_CONCEPT_ASSERT((CollisionCheckVisitorConcept<Visitor, PositionMap>));
  
  BOOST_CONCEPT_USAGE(PRMVisitorConcept) { };
  
};





};

};


#endif

