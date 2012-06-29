/**
 * \file tree_concepts.hpp
 * 
 * This library provides concept classes to verify that a data structure has a tree interface, 
 * i.e., models the concepts of trees as used in ReaK::graph. The tree concepts reflect the 
 * same style of interface specification as the Boost.Graph Library, and reuses much of the 
 * general graph traits (see boost::graph_traits). In other words, in a similar way as different
 * graph concepts in the BGL require certain functions on the graph, the tree concepts require 
 * a few additional functions. In practice, one can use a general graph structure to represent 
 * a tree (as a special kind of graph), and thus, using simple wrapper functions like those in 
 * the bgl_tree_adaptor.hpp header, one can make a graph look like a tree. However, 
 * when using a general graph structure as a tree, to be pedantic, one should disable the general
 * functions of the MutableGraphConcept (and MutablePropertyGraphConcept) (e.g., through a tree-view 
 * adaptor on the graph), but in practice, if one uses a general graph in an algorithm that 
 * constructs a tree (via tree-concept functions only), then there shouldn't be any problems.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
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

#ifndef REAK_TREE_CONCEPTS_HPP
#define REAK_TREE_CONCEPTS_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace graph {
  
  
/**
 * This traits class defines a number of nested types associated to a tree structure (TreeConcept).
 */
template <typename TreeType>
struct tree_traits {
  /** This type describes iterators to iterate through child vertices of a vertex. */
  typedef typename TreeType::child_vertex_iterator child_vertex_iterator;
};



/**
 * This concept defines the requirements to fulfill in order to model a tree 
 * as used in ReaK::graph.
 * 
 * Required Concepts:
 * 
 * TreeType should model the GraphConcept and the IncidenceGraphConcept.
 * 
 * Valid expressions (vertex v; vertex_iterator cv_it, cv_it_end; TreeType tree):
 * 
 * v = get_root_vertex(tree);  The root vertex of the tree can be obtained.
 * 
 * tie(cv_it, cv_it_end) = child_vertices(v,tree);  The iterator range for child vertices of a given parent vertex (v) in a tree can be obtained.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct TreeConcept {
  typename boost::graph_traits<TreeType>::vertex_descriptor v;
  typename tree_traits<TreeType>::child_vertex_iterator cv_it, cv_it_end;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(TreeConcept) 
  {
    v = get_root_vertex(tree);
    boost::tie(cv_it, cv_it_end) = child_vertices(v,tree);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a mutable tree 
 * as used in ReaK::graph.
 * 
 * Required concepts:
 * 
 * TreeType should model the TreeConcept.
 * 
 * Valid expressions (vertex u, v; edge e; TreeType tree):
 * 
 * v = create_root(tree);  The root vertex of the tree can be created.
 * 
 * tie(v, e) = add_child_vertex(u, tree);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) in a tree.
 * 
 * remove_branch(v, tree);  The entire branch (or sub-tree) below a given vertex (including the vertex itself) can be removed from the tree.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct MutableTreeConcept {
  typename boost::graph_traits<TreeType>::vertex_descriptor u, v;
  typename boost::graph_traits<TreeType>::edge_descriptor e;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((TreeConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(MutableTreeConcept) 
  {
    v = create_root(tree);
    boost::tie(v,e) = add_child_vertex(u, tree);
    remove_branch(v, tree);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a mutable property-tree 
 * as used in ReaK::graph. A mutable property-tree is essentially a mutable tree whose mutating functions
 * take or deliver the vertex- or edge-property values associated with the vertices or edges.
 * During removal of a branch, all the vertex-properties are collected into an output iterator (e.g., back-inserter).
 * During additions of child nodes, the corresponding vertex-properties can be used to initialize the new 
 * vertex and edge directly. This not only makes such a tree easier to use (not having to manually collect
 * or set vertex properties before or after the mutation), but it can also be necessary in some situations.
 * One typical use-case is when re-balancing a branch of a tree, which often results in collecting properties 
 * of vertices (to preserve them), clearing the branch, re-adding the vertices in a new arrangement, and 
 * restoring their original properties. Now, this use-case is much simpler to do and more efficient if the 
 * mutable property-tree concept is modeled by the tree type in question.
 * Note, this concept is expected to take advantage of C++11 move-semantics on a C++11 capable compiler.
 * 
 * 
 * Required concepts:
 * 
 * TreeType should model the TreeConcept.
 * 
 * Valid expressions (vertex u, v; edge e; TreeType tree):
 * 
 * v = create_root(vp, g);  The root vertex of the tree can be created with a given vertex-property value (vp).
 * 
 * tie(v,e) = add_child_vertex(u, vp, g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) and with a given vertex-property value (vp) in a tree.
 * 
 * tie(v,e) = add_child_vertex(u, vp, ep, g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u), with a given vertex-property value (vp) and  with a given edge-property value (ep) in a tree.
 * 
 * remove_branch(v, back_inserter(vp_vect), g);  The entire branch (or sub-tree) below a given vertex (including the vertex itself) can be removed from the tree, with its vertex-property values put on an output iterator (e.g., a back-inserter).
 * 
 * C++11 only:
 * 
 * v = create_root(std::move(vp), g);  The root vertex of the tree can be created with a given vertex-property value (vp) (to move in).
 *
 * tie(v,e) = add_child_vertex(u, std::move(vp), g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) and with a given vertex-property value (vp) (to move in) in a tree.
 *
 * tie(v,e) = add_child_vertex(u, std::move(vp), std::move(ep), g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u), with a given vertex-property value (vp) (to move in) and  with a given edge-property value (ep) (to move in) in a tree.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct MutablePropertyTreeConcept {
  typename boost::graph_traits<TreeType>::vertex_descriptor u, v;
  typename boost::graph_traits<TreeType>::edge_descriptor e;
  typename TreeType::vertex_property_type vp;
  typename TreeType::edge_property_type ep;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((TreeConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(MutablePropertyTreeConcept) 
  {
    v = create_root(vp, tree);
    boost::tie(v,e) = add_child_vertex(u, vp, tree);
    boost::tie(v,e) = add_child_vertex(u, vp, ep, tree);
    std::vector< typename TreeType::vertex_property_type > vp_vect;
    remove_branch(v, back_inserter(vp_vect), tree);
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    v = create_root(std::move(vp), tree);
    boost::tie(v,e) = add_child_vertex(u, std::move(vp), tree);
    boost::tie(v,e) = add_child_vertex(u, std::move(vp), std::move(ep), tree);
#endif
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a non-compact graph
 * as used in ReaK::graph. A non-compact graph arises from the fact that certain data layouts
 * used to store graphs (and trees) may leave holes (or empty vertices or edges) in order to 
 * preserve the current layout of vertices and edges in memory (e.g., avoid copying and rewiring
 * vertex-pointers, or keeping some regular memory pattern like a breadth-first layout). 
 * Consequently, when iterating through vertices or edges on such a non-compact graph, it is 
 * necessary to test the validity of each vertex or edge, since the iteration might encounter
 * holes (or empty vertices or edges).
 * 
 * Required concepts:
 * 
 * TreeType should model the boost::GraphConcept.
 * 
 * Valid expressions (vertex v; edge e; GraphType g):
 * 
 * b = is_vertex_valid(v,g);  A vertex of the graph can be tested for validity.
 * 
 * b = is_edge_valid(e,g);  An edge of the graph can be tested for validity.
 * 
 * \tparam GraphType The graph type to be checked for this concept.
 */
template <typename GraphType>
struct NonCompactGraphConcept {
  typename boost::graph_traits<GraphType>::vertex_descriptor v;
  typename boost::graph_traits<GraphType>::edge_descriptor e;
  GraphType g;
  
  BOOST_CONCEPT_ASSERT((boost::GraphConcept<GraphType>));
  
  BOOST_CONCEPT_USAGE(NonCompactGraphConcept) 
  {
    bool b = is_vertex_valid(v,g);
    b = is_edge_valid(e,g);
  };
  
};


};

};



#endif


















