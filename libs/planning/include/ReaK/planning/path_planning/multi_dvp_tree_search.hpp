/**
 * \file multi_dvp_tree_search.hpp
 *
 * This library provides a class that implements an adaptor that can act as a
 * neighborhood searching functor by calling upon a dvp-tree space-partitioning
 * of the nodes of an associated graph and over a given topology.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
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

#ifndef REAK_MULTI_DVP_TREE_SEARCH_HPP
#define REAK_MULTI_DVP_TREE_SEARCH_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

#include <map>  // for map in multi_dvp_tree_search

#include <limits>   // for numeric_limits
#include <utility>  // for pair and move

namespace ReaK::pp {

/**
 * This class template is used as a type of kind of associative container for DVP-trees that
 * span distinct graphs. For problems in which multiple graphs exist and nearest-neighbors
 * may be queried from any of these graphs, one can use this class to associate each graph
 * with a dvp-tree. This class is callable as a single query and a KNN / range query.
 * \tparam Graph The graph type which can contain the vertices.
 * \tparam DVPTree The DVP-tree type that is used to perform the nearest-neighbor queries.
 */
template <typename Graph, typename DVPTree, bool IsDirected>
struct multi_dvp_tree_search_base {
  /** This associative container is used to store and access the DVP-trees. */
  std::map<Graph*, DVPTree*> graph_tree_map;

  multi_dvp_tree_search_base() : graph_tree_map() {}

  /**
   * This function computes an approximation of the characteristic size of the vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param g A graph containing the vertices from which to find the nearest-neighbor.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \return The approximation of the characteristic size of the vertices of the given graph.
   */
  template <typename Topology, typename PositionMap>
  double get_characteristic_size(Graph& g, const Topology& space,
                                 PositionMap position) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      return it->second->get_characteristic_size();
    }
    return std::numeric_limits<double>::infinity();
  }

  /**
   * This call-operator finds the nearest vertex of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param p A position in the space, to which the nearest-neighbor is sought.
   * \param g A graph containing the vertices from which to find the nearest-neighbor.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   */
  template <typename Topology, typename PositionMap>
  graph::graph_vertex_t<Graph> operator()(
      const graph::property_value_t<PositionMap>& p, Graph& g,
      const Topology& space, PositionMap position) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      return it->second->find_nearest(p);
    }
    return boost::graph_traits<Graph>::null_vertex();
  }

  /**
   * This call-operator finds the nearest vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions.
   * \tparam PositionMap The property-map type which can store the position associated
   *         with each vertex.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param output_first An iterator to the first place where to put the sorted list of
   *        elements with the smallest distance.
   * \param g A graph containing the vertices from which to find the nearest-neighbors.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \param max_neighbors The maximum number of neighbors to have in the list.
   * \param radius The minimum distance around the position that a vertex should be in to be
   *        considered a neighbor.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename Topology, typename PositionMap, typename OutputIterator>
  OutputIterator operator()(
      const graph::property_value_t<PositionMap>& p,
      OutputIterator output_first, Graph& g, const Topology& space,
      PositionMap position, std::size_t max_neighbors = 1,
      double radius = std::numeric_limits<double>::infinity()) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      return it->second->find_nearest(p, output_first, max_neighbors, radius);
    }
    return output_first;
  }

  /**
   * This is a call-back for when a vertex has been added.
   */
  template <typename Vertex>
  void added_vertex(Vertex v, Graph& g) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      it->second->insert(v);
    }
  }

  /**
   * This is a call-back for when a vertex is about to be removed.
   */
  template <typename Vertex>
  void removed_vertex(Vertex v, Graph& g) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      it->second->erase(v);
    }
  }
};

/**
 * This class template is used as a type of kind of associative container for DVP-trees that
 * span distinct graphs. For problems in which multiple graphs exist and nearest-neighbors
 * may be queried from any of these graphs, one can use this class to associate each graph
 * with a dvp-tree. This class is callable as a single query and a KNN / range query.
 * \tparam Graph The graph type which can contain the vertices.
 * \tparam DVPTree The DVP-tree type that is used to perform the nearest-neighbor queries.
 */
template <typename Graph, typename DVPTree>
struct multi_dvp_tree_search_base<Graph, DVPTree, true> {
  /** This associative container is used to store and access the DVP-trees. */
  std::map<Graph*, DVPTree*> graph_tree_map;

  multi_dvp_tree_search_base() : graph_tree_map() {}

  /**
   * This function computes an approximation of the characteristic size of the vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param g A graph containing the vertices from which to find the nearest-neighbor.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \return The approximation of the characteristic size of the vertices of the given graph.
   */
  template <typename Topology, typename PositionMap>
  double get_characteristic_size(Graph& g, const Topology& space,
                                 PositionMap position) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second != nullptr)) {
      return it->second->get_characteristic_size();
    }
    return std::numeric_limits<double>::infinity();
  }

  /**
   * This call-operator finds the nearest predecesor and successor of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param g A graph containing the vertices from which to find the nearest-neighbor.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \return A pair containing the nearest predecessor and successor vertex.
   */
  template <typename Topology, typename PositionMap>
  std::pair<graph::graph_vertex_t<Graph>, graph::graph_vertex_t<Graph>>
  operator()(const graph::property_value_t<PositionMap>& p, Graph& g,
             const Topology& space, PositionMap position) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second)) {
      return it->second->find_nearest_pred_succ(p);
    }
    return {boost::graph_traits<Graph>::null_vertex(),
            boost::graph_traits<Graph>::null_vertex()};
  }

  /**
   * This call-operator finds the nearest vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions.
   * \tparam PositionMap The property-map type which can store the position associated
   *         with each vertex.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param pred_first An iterator to the first place where to put the sorted list of
   *        the best predecessor elements.
   * \param succ_first An iterator to the first place where to put the sorted list of
   *        the best successor elements.
   * \param g A graph containing the vertices from which to find the nearest-neighbors.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \param max_neighbors The maximum number of neighbors to have in the list.
   * \param radius The minimum distance around the position that a vertex should be in to be
   *        considered a neighbor.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename Topology, typename PositionMap, typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> operator()(
      const graph::property_value_t<PositionMap>& p, OutputIterator pred_first,
      OutputIterator succ_first, Graph& g, const Topology& space,
      PositionMap position, std::size_t max_neighbors = 1,
      double radius = std::numeric_limits<double>::infinity()) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second)) {
      return it->second->find_nearest(p, pred_first, succ_first, max_neighbors,
                                      radius);
    }
    return {pred_first, succ_first};
  }

  /**
   * This is a call-back for when a vertex has been added.
   */
  template <typename Vertex>
  void added_vertex(Vertex v, Graph& g) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second)) {
      it->second->insert(v);
    }
  }

  /**
   * This is a call-back for when a vertex is about to be removed.
   */
  template <typename Vertex>
  void removed_vertex(Vertex v, Graph& g) const {
    auto it = graph_tree_map.find(&g);
    if ((it != graph_tree_map.end()) && (it->second)) {
      it->second->erase(v);
    }
  }
};

template <typename Graph, typename DVPTree>
struct multi_dvp_tree_search
    : multi_dvp_tree_search_base<Graph, DVPTree,
                                 boost::is_directed_graph<Graph>::type::value> {
  multi_dvp_tree_search()
      : multi_dvp_tree_search_base<
            Graph, DVPTree, boost::is_directed_graph<Graph>::type::value>() {}
};

template <typename Graph, typename DVPTree>
struct multi_dvp_tree_pred_succ_search
    : multi_dvp_tree_search_base<Graph, DVPTree, true> {
  multi_dvp_tree_pred_succ_search()
      : multi_dvp_tree_search_base<Graph, DVPTree, true>() {}
};

}  // namespace ReaK::pp

#endif
