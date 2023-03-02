/**
 * \file topological_search.hpp
 *
 * This library contains two simple nearest-neighbor search algorithms implemented
 * as functor templates. This library contains, in fact, three algorithms.
 *
 * First, a simple min_dist_linear_search algorithm is provided which is similar to std::min_element
 * except that it stores and compares the best distance value associated to the best iterator
 * to the current one (it is a simple linear search that avoid recomputation of the distance
 * at every iteration, which would be required if std::min_element was used instead).
 *
 * Second, a linear_neighbor_search algorithm is provided which simply does an exhaustive linear
 * search through all the vertices of a graph to find the nearest one to a given point, in a given
 * topology (as of topologies in the Boost Graph Library). This algorithms simply wraps the
 * min_dist_linear_search with the required distance and comparison function.
 *
 * Third, a best_only_neighbor_search algorithm is provided which is an approximation to an
 * exhaustive linear search by picking a number of random vertices from the graph and performing
 * a best only search down the graph to find the "nearest-neighbor". This is, of course, not going
 * to find the nearest-neighbor, but can significantly cut down on query time if finding the
 * nearest neighbor is not a strict requirement in the algorithm.
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

#ifndef REAK_TOPOLOGICAL_SEARCH_HPP
#define REAK_TOPOLOGICAL_SEARCH_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/topology.hpp>
#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <iterator>
#include <queue>
#include <tuple>
#include <type_traits>
#include <vector>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>

#include <ReaK/core/base/misc_math.hpp>

namespace ReaK::pp {

/**
  * This function template is similar to std::min_element but can be used when the comparison
  * involves computing a derived quantity (a.k.a. distance). This algorithm will search for the
  * the element in the range [first,last) which has the "smallest" distance (of course, both the
  * distance metric and comparison can be overriden to perform something other than the canonical
  * Euclidean distance and less-than comparison, which would yield the element with minimum distance).
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param compare A callable object that returns true if the first element is the preferred one (less-than) of the two.
  * \param inf A DistanceValue which represents infinity (i.e. the very worst value with which to initialize the
 * search).
  * \return The iterator to the best element in the range (best is defined as the one which would compare favorably to
 * all the elements in the range with respect to the distance metric).
  */
template <typename DistanceValue, typename ForwardIterator,
          typename GetDistanceFunction, typename CompareFunction>
inline ForwardIterator min_dist_linear_search(
    ForwardIterator first, ForwardIterator last, GetDistanceFunction distance,
    CompareFunction compare,
    DistanceValue inf = std::numeric_limits<DistanceValue>::infinity()) {
  if (first == last) {
    return last;
  }
  DistanceValue d_best = inf;
  ForwardIterator result = last;
  for (; first != last; ++first) {
    DistanceValue d = distance(*first);
    if (compare(d, d_best)) {
      d_best = d;
      result = first;
    }
  }
  return result;
}

/**
  * This function template is a specialization of min_dist_linear_search for the default comparison
  * function which is the less-than operator.
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param inf A DistanceValue which represents infinity (i.e. the very worst value with which to initialize the
 * search).
  * \return The iterator to the best element in the range (best is defined as the one which would compare favorably to
 * all the elements in the range with respect to the distance metric).
  */
template <typename DistanceValue, typename ForwardIterator,
          typename GetDistanceFunction>
inline ForwardIterator min_dist_linear_search(
    ForwardIterator first, ForwardIterator last, GetDistanceFunction distance,
    DistanceValue inf = std::numeric_limits<DistanceValue>::infinity()) {
  return min_dist_linear_search(first, last, distance, std::less<>(), inf);
}

/**
  * This function template is similar to std::min_element but can be used when the comparison
  * involves computing a derived quantity (a.k.a. distance). This algorithm will search for the
  * the elements in the range [first,last) with the "smallest" distances (of course, both the
  * distance metric and comparison can be overriden to perform something other than the canonical
  * Euclidean distance and less-than comparison, which would yield the element with minimum distance).
  * This function will fill the output container with a number of nearest-neighbors.
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam OutputIterator The forward- output-iterator type which can contain the list of nearest-neighbors.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param output_first An iterator to the first place where to put the sorted list of elements with the smallest
 * distance.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param compare A callable object that returns true if the first element is the preferred one (less-than) of the two.
  * \param max_neighbors The maximum number of elements of smallest distance to output in the sorted list.
  * \param radius The maximum distance value for which an element qualifies to be part of the output list.
  * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
  */
template <typename DistanceValue, typename ForwardIterator,
          typename OutputIterator, typename GetDistanceFunction,
          typename CompareFunction>
inline OutputIterator min_dist_linear_search(
    ForwardIterator first, ForwardIterator last, OutputIterator output_first,
    GetDistanceFunction distance, CompareFunction compare,
    std::size_t max_neighbors = 1,
    DistanceValue radius = std::numeric_limits<DistanceValue>::infinity()) {
  if (first == last) {
    return output_first;
  }
  const auto p_compare = [&](const auto& lhs, const auto& rhs) {
    return compare(lhs.first, rhs.first);
  };
  std::vector<std::pair<DistanceValue, ForwardIterator>> output_queue;
  for (; first != last; ++first) {
    DistanceValue d = distance(*first);
    if (!compare(d, radius)) {
      continue;
    }
    output_queue.emplace_back(d, first);
    std::push_heap(output_queue.begin(), output_queue.end(), p_compare);
    if (output_queue.size() > max_neighbors) {
      std::pop_heap(output_queue.begin(), output_queue.end(), p_compare);
      output_queue.pop_back();
      radius = output_queue.front().first;
    }
  }
  std::sort_heap(output_queue.begin(), output_queue.end(), p_compare);
  for (auto [d, it] : output_queue) {
    RK_UNUSED(d);
    *(output_first++) = *it;
  }
  return output_first;
}

/**
  * This function template is similar to std::min_element but can be used when the comparison
  * involves computing a derived quantity (a.k.a. distance). This algorithm will search for the
  * the element in the range [first,last) which has the "smallest" distance (of course, both the
  * distance metric and comparison can be overriden to perform something other than the canonical
  * Euclidean distance and less-than comparison, which would yield the element with minimum distance).
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam OutputIterator The forward- output-iterator type which can contain the list of nearest-neighbors.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param output_first An iterator to the first place where to put the sorted list of elements with the smallest
 * distance.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param max_neighbors The maximum number of elements of smallest distance to output in the sorted list.
  * \param radius The maximum distance value for which an element qualifies to be part of the output list.
  * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
  */
template <typename DistanceValue, typename ForwardIterator,
          typename OutputIterator, typename GetDistanceFunction>
inline OutputIterator min_dist_linear_search(
    ForwardIterator first, ForwardIterator last, OutputIterator output_first,
    GetDistanceFunction distance, std::size_t max_neighbors = 1,
    DistanceValue radius = std::numeric_limits<DistanceValue>::infinity()) {
  return min_dist_linear_search(first, last, output_first, distance,
                                std::less<>(), max_neighbors, radius);
}

/**
  * This function template is similar to std::min_element but can be used when the comparison
  * involves computing a derived quantity (a.k.a. distance). This algorithm will search for the
  * the elements in the range [first,last) with the "smallest" distances.
  * This version assumes the distance metric is not symmetric and thus, distinguishes between
  * the best predecessors and successors with respect to the query point.
  * This function will fill the output containers with a number of nearest-neighbors.
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam OutputIterator The forward- output-iterator type which can contain the list of nearest-neighbors.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param pred_first An iterator to the first place where to put the sorted list of best predecessor elements.
  * \param succ_first An iterator to the first place where to put the sorted list of best successor elements.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param compare A callable object that returns true if the first element is the preferred one (less-than) of the two.
  * \param max_neighbors The maximum number of elements of smallest distance to output in the sorted list.
  * \param radius The maximum distance value for which an element qualifies to be part of the output list.
  * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
  */
template <typename DistanceValue, typename PointType, typename ForwardIterator,
          typename OutputIterator, typename GetDistanceFunction,
          typename CompareFunction>
inline std::pair<OutputIterator, OutputIterator> min_dist_linear_search(
    const PointType& query_point, ForwardIterator first, ForwardIterator last,
    OutputIterator pred_first, OutputIterator succ_first,
    GetDistanceFunction distance, CompareFunction compare,
    std::size_t max_neighbors = 1,
    DistanceValue radius = std::numeric_limits<DistanceValue>::infinity()) {
  if (first == last) {
    return {pred_first, succ_first};
  }
  const auto p_compare = [&](const auto& lhs, const auto& rhs) {
    return compare(lhs.first, rhs.first);
  };
  DistanceValue radius_pred = radius;
  DistanceValue radius_succ = radius;
  std::vector<std::pair<DistanceValue, ForwardIterator>> pred_queue;
  std::vector<std::pair<DistanceValue, ForwardIterator>> succ_queue;
  for (; first != last; ++first) {
    DistanceValue d_pred = distance(*first, query_point);
    if (compare(d_pred, radius_pred)) {
      pred_queue.emplace_back(d_pred, first);
      std::push_heap(pred_queue.begin(), pred_queue.end(), p_compare);
      if (pred_queue.size() > max_neighbors) {
        std::pop_heap(pred_queue.begin(), pred_queue.end(), p_compare);
        pred_queue.pop_back();
        radius_pred = pred_queue.front().first;
      };
    };
    DistanceValue d_succ = distance(query_point, *first);
    if (compare(d_succ, radius_succ)) {
      succ_queue.emplace_back(d_succ, first);
      std::push_heap(succ_queue.begin(), succ_queue.end(), p_compare);
      if (succ_queue.size() > max_neighbors) {
        std::pop_heap(succ_queue.begin(), succ_queue.end(), p_compare);
        succ_queue.pop_back();
        radius_succ = succ_queue.front().first;
      }
    }
  }
  std::sort_heap(pred_queue.begin(), pred_queue.end(), p_compare);
  for (auto [d, it] : pred_queue) {
    RK_UNUSED(d);
    *(pred_first++) = *it;
  }
  std::sort_heap(succ_queue.begin(), succ_queue.end(), p_compare);
  for (auto [d, it] : succ_queue) {
    RK_UNUSED(d);
    *(succ_first++) = *it;
  }
  return {pred_first, succ_first};
}

/**
  * This function template is similar to std::min_element but can be used when the comparison
  * involves computing a derived quantity (a.k.a. distance). This algorithm will search for the
  * the element in the range [first,last) which has the "smallest" distance (of course, both the
  * distance metric and comparison can be overriden to perform something other than the canonical
  * Euclidean distance and less-than comparison, which would yield the element with minimum distance).
  * \tparam DistanceValue The value-type for the distance measures.
  * \tparam ForwardIterator The forward-iterator type.
  * \tparam OutputIterator The forward- output-iterator type which can contain the list of nearest-neighbors.
  * \tparam GetDistanceFunction The functor type to compute the distance measure.
  * \param first Start of the range in which to search.
  * \param last One element past the last element in the range in which to search.
  * \param pred_first An iterator to the first place where to put the sorted list of best predecessor elements.
  * \param succ_first An iterator to the first place where to put the sorted list of best successor elements.
  * \param distance A callable object that returns a DistanceValue for a given element from the ForwardIterator
 * dereferencing.
  * \param max_neighbors The maximum number of elements of smallest distance to output in the sorted list.
  * \param radius The maximum distance value for which an element qualifies to be part of the output list.
  * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
  */
template <typename DistanceValue, typename PointType, typename ForwardIterator,
          typename OutputIterator, typename GetDistanceFunction>
inline std::pair<OutputIterator, OutputIterator> min_dist_linear_search(
    const PointType& query_point, ForwardIterator first, ForwardIterator last,
    OutputIterator pred_first, OutputIterator succ_first,
    GetDistanceFunction distance, std::size_t max_neighbors = 1,
    DistanceValue radius = std::numeric_limits<DistanceValue>::infinity()) {
  return min_dist_linear_search(query_point, first, last, pred_first,
                                succ_first, distance, std::less<>(),
                                max_neighbors, radius);
}

/**
  * This class is a no-op synchronization class. Some nearest-neighbor methods
  * might require a call-back whenever a vertex has been added or when it is about
  * to be removed in order to update the data-structure maintained to perform the
  * NN queries. This class can be used when no such synchronization is needed.
  */
struct no_NNfinder_synchro {

  /**
    * This is a call-back for when a vertex has been added.
    */
  template <typename Vertex, typename Graph>
  void added_vertex(Vertex, const Graph&) const {}

  /**
    * This is a call-back for when a vertex is about to be removed.
    */
  template <typename Vertex, typename Graph>
  void removed_vertex(Vertex, const Graph&) const {}
};

/**
  * This class is a composite synchronization class. Some nearest-neighbor methods
  * might require a call-back whenever a vertex has been added or when it is about
  * to be removed in order to update the data-structure maintained to perform the
  * NN queries. This class can be used when no such synchronization is needed.
  */
template <typename NNSynchro1, typename NNSynchro2>
struct composite_NNsynchro {
  NNSynchro1 s1;
  NNSynchro2 s2;

  composite_NNsynchro(NNSynchro1 aS1 = NNSynchro1(),
                      NNSynchro2 aS2 = NNSynchro2())
      : s1(aS1), s2(aS2) {}

  /**
    * This is a call-back for when a vertex has been added.
    */
  template <typename Vertex, typename Graph>
  void added_vertex(Vertex u, Graph& g) const {
    s1.added_vertex(u, g);
    s2.added_vertex(u, g);
  }

  /**
    * This is a call-back for when a vertex is about to be removed.
    */
  template <typename Vertex, typename Graph>
  void removed_vertex(Vertex u, Graph& g) const {
    s1.removed_vertex(u, g);
    s2.removed_vertex(u, g);
  }
};

namespace detail {
namespace {

template <typename Topology, typename PositionMap>
struct linear_neighbor_search_distance_functor {
  using PointType = graph::property_value_t<PositionMap>;
  const PointType* p_point;
  const Topology* p_space;
  PositionMap position;

  linear_neighbor_search_distance_functor(const PointType* pPoint,
                                          const Topology* pSpace,
                                          PositionMap aPosition)
      : p_point(pPoint), p_space(pSpace), position(aPosition) {}

  template <typename Vertex>
  double operator()(Vertex u) const {
    return get(distance_metric, *p_space)(*p_point, get(position, u), *p_space);
  }

  template <typename Vertex>
  double operator()(const PointType& p, Vertex u) const {
    return get(distance_metric, *p_space)(p, get(position, u), *p_space);
  }

  template <typename Vertex>
  double operator()(Vertex u, const PointType& p) const {
    return get(distance_metric, *p_space)(get(position, u), p, *p_space);
  }
};

}  // namespace
}  // namespace detail

/**
  * This functor template performs a linear nearest-neighbor search through a graph by invoquing
  * the distance function of an underlying topology. The call operator will return the vertex
  * of the graph whose position value is closest to a given position value.
  * \tparam Graph The graph type which can contain the vertices.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  */
template <typename Graph, typename CompareFunction, bool IsDirected>
struct linear_neighbor_search_base {

  CompareFunction m_compare;
  /**
    * Default constructor.
    * \param compare The comparison functor for ordering the distances (strict weak ordering).
    */
  linear_neighbor_search_base(CompareFunction compare = CompareFunction())
      : m_compare(compare) {}

  /**
    * This function computes an approximation of the characteristic size of the vertices of a graph.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated with each vertex.
    * \param g A graph containing the vertices from which to find the nearest-neighbor.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    */
  template <typename Topology, typename PositionMap>
  double get_characteristic_size(Graph& g, const Topology& space,
                                 PositionMap position) const {
    auto [ui, ui_end] = vertices(g);
    std::size_t N = num_vertices(g);
    std::size_t log_N = math::highest_set_bit(N) + 1;
    std::size_t step = N / log_N;

    auto ui_cur = ui;
    double avg_max_dist = 0.0;
    for (std::size_t i = 0; i < log_N; ++i) {
      // the double half advance is to avoid going over bound at the end.
      std::advance(ui_cur, step / 2);
      double d_best = 0.0;
      for (auto ui_test = ui; ui_test != ui_end; ++ui_test) {
        double d = get(distance_metric, space)(get(position, *ui_cur),
                                               get(position, *ui_test), space);
        if ((d != std::numeric_limits<double>::infinity()) && (d > d_best)) {
          d_best = d;
        }
      }
      if (d_best != 0.0) {
        avg_max_dist += d_best;
      }
      std::advance(ui_cur, step / 2);
    }
    if (avg_max_dist != 0.0) {
      return avg_max_dist / log_N;
    }
    // never found a finite, non-zero distance value.
    return std::numeric_limits<double>::infinity();
  }

  /**
    * This call-operator finds the nearest vertex of a graph, to a given position.
    * \tparam Topology The topology type which contains the positions.
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
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    auto [ui, ui_end] = vertices(g);
    return *(min_dist_linear_search(
        ui, ui_end,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, std::numeric_limits<double>::infinity()));
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
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    auto [ui, ui_end] = vertices(g);
    return min_dist_linear_search(
        ui, ui_end, output_first,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, max_neighbors, radius);
  }

  /**
    * This call-operator finds the nearest vertex of a graph, to a given position.
    * \tparam ForwardIter A forward-iterator type.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated with each vertex.
    * \param p A position in the space, to which the nearest-neighbor is sought.
    * \param first The first of all candidates nearest-neighbors.
    * \param last The last of all candidates nearest-neighbors (one-passed-last, as is usual in C++).
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    */
  template <typename ForwardIter, typename Topology, typename PositionMap>
  ForwardIter operator()(const graph::property_value_t<PositionMap>& p,
                         ForwardIter first, ForwardIter last,
                         const Topology& space, PositionMap position) const {
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    return min_dist_linear_search(
        first, last,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, std::numeric_limits<double>::infinity());
  }

  /**
    * This call-operator finds the nearest vertices of a graph, to a given position.
    * \tparam ForwardIter A forward-iterator type.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated
    *         with each vertex.
    * \tparam OutputIterator The forward- output-iterator type which can contain the
    *         list of nearest-neighbors.
    * \param p A position in the space, to which the nearest-neighbors are sought.
    * \param first The first of all candidates nearest-neighbors.
    * \param last The last of all candidates nearest-neighbors (one-passed-last, as is usual in C++).
    * \param output_first An iterator to the first place where to put the sorted list of
    *        elements with the smallest distance.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    * \param max_neighbors The maximum number of neighbors to have in the list.
    * \param radius The minimum distance around the position that a vertex should be in to be
    *        considered a neighbor.
    * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
    */
  template <typename ForwardIter, typename Topology, typename PositionMap,
            typename OutputIterator>
  OutputIterator operator()(
      const graph::property_value_t<PositionMap>& p, ForwardIter first,
      ForwardIter last, OutputIterator output_first, const Topology& space,
      PositionMap position, std::size_t max_neighbors = 1,
      double radius = std::numeric_limits<double>::infinity()) const {
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    return min_dist_linear_search(
        first, last, output_first,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, max_neighbors, radius);
  }
};

/**
  * This functor template performs a linear nearest-neighbor search through a graph by invoquing
  * the distance function of an underlying topology. The call operator will return the vertex
  * of the graph whose position value is closest to a given position value.
  * \tparam Graph The graph type which can contain the vertices.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  */
template <typename Graph, typename CompareFunction>
struct linear_neighbor_search_base<Graph, CompareFunction, true> {

  CompareFunction m_compare;
  /**
    * Default constructor.
    * \param compare The comparison functor for ordering the distances (strict weak ordering).
    */
  linear_neighbor_search_base(CompareFunction compare = CompareFunction())
      : m_compare(compare) {}

  /**
    * This function computes an approximation of the characteristic size of the vertices of a graph.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated with each vertex.
    * \param g A graph containing the vertices from which to find the nearest-neighbor.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    */
  template <typename Topology, typename PositionMap>
  double get_characteristic_size(Graph& g, const Topology& space,
                                 PositionMap position) const {
    auto [ui, ui_end] = vertices(g);
    std::size_t N = num_vertices(g);
    std::size_t log_N = math::highest_set_bit(N) + 1;
    std::size_t step = N / log_N;

    auto ui_cur = ui;
    double avg_max_dist = 0.0;
    for (std::size_t i = 0; i < log_N; ++i) {
      // the double half advance is to avoid going over bound at the end.
      std::advance(ui_cur, step / 2);
      double d_best = 0.0;
      for (auto ui_test = ui; ui_test != ui_end; ++ui_test) {
        double d1 = get(distance_metric, space)(get(position, *ui_cur),
                                                get(position, *ui_test), space);
        if ((d1 != std::numeric_limits<double>::infinity()) && (d1 > d_best)) {
          d_best = d1;
        }
        double d2 = get(distance_metric, space)(get(position, *ui_test),
                                                get(position, *ui_cur), space);
        if ((d2 != std::numeric_limits<double>::infinity()) && (d2 > d_best)) {
          d_best = d2;
        }
      }
      if (d_best != 0.0) {
        avg_max_dist += d_best;
      }
      std::advance(ui_cur, step / 2);
    }
    if (avg_max_dist != 0.0) {
      return avg_max_dist / log_N;
    }
    // never found a finite, non-zero distance value.
    return std::numeric_limits<double>::infinity();
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
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    auto [ui, ui_end] = vertices(g);
    return min_dist_linear_search(
        p, ui, ui_end, pred_first, succ_first,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, max_neighbors, radius);
  }

  /**
    * This call-operator finds the nearest vertices of a graph, to a given position.
    * \tparam ForwardIter A forward-iterator type.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated
    *         with each vertex.
    * \tparam OutputIterator The forward- output-iterator type which can contain the
    *         list of nearest-neighbors.
    * \param p A position in the space, to which the nearest-neighbors are sought.
    * \param first The first of all candidates nearest-neighbors.
    * \param last The last of all candidates nearest-neighbors (one-passed-last, as is usual in C++).
    * \param pred_first An iterator to the first place where to put the sorted list of
    *        the best predecessor elements.
    * \param succ_first An iterator to the first place where to put the sorted list of
    *        the best successor elements.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    * \param max_neighbors The maximum number of neighbors to have in the list.
    * \param radius The minimum distance around the position that a vertex should be in to be
    *        considered a neighbor.
    * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
    */
  template <typename ForwardIter, typename Topology, typename PositionMap,
            typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> operator()(
      const graph::property_value_t<PositionMap>& p, ForwardIter first,
      ForwardIter last, OutputIterator pred_first, OutputIterator succ_first,
      const Topology& space, PositionMap position,
      std::size_t max_neighbors = 1,
      double radius = std::numeric_limits<double>::infinity()) const {
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    return min_dist_linear_search(
        p, first, last, pred_first, succ_first,
        detail::linear_neighbor_search_distance_functor<Topology, PositionMap>(
            &p, &space, position),
        m_compare, max_neighbors, radius);
  }
};

/**
  * This functor template performs a linear nearest-neighbor search through a graph by invoquing
  * the distance function of an underlying topology. The call operator will return the vertex
  * of the graph whose position value is closest to a given position value.
  * \tparam Graph The graph type which can contain the vertices.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  */
template <typename Graph, typename CompareFunction = std::less<>>
struct linear_neighbor_search
    : linear_neighbor_search_base<
          Graph, CompareFunction,
          boost::is_directed_graph<Graph>::type::value> {
  /**
    * Default constructor.
    * \param compare The comparison functor for ordering the distances (strict weak ordering).
    */
  explicit linear_neighbor_search(CompareFunction compare)
      : linear_neighbor_search_base<
            Graph, CompareFunction,
            boost::is_directed_graph<Graph>::type::value>(compare) {}

  linear_neighbor_search() : linear_neighbor_search(CompareFunction()) {}
};

/**
  * This functor template performs a linear nearest-neighbor search through a graph by invoquing
  * the distance function of an underlying topology. The call operator will return the vertex
  * of the graph whose position value is closest to a given position value.
  * \tparam Graph The graph type which can contain the vertices.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  */
template <typename Graph, typename CompareFunction = std::less<>>
struct linear_pred_succ_search
    : linear_neighbor_search_base<Graph, CompareFunction, true> {
  /**
    * Default constructor.
    * \param compare The comparison functor for ordering the distances (strict weak ordering).
    */
  explicit linear_pred_succ_search(CompareFunction compare)
      : linear_neighbor_search_base<Graph, CompareFunction, true>(compare) {}

  linear_pred_succ_search() : linear_pred_succ_search(CompareFunction()) {}
};

/**
  * This functor template performs a best-only nearest-neighbor search through a tree by invoquing
  * the distance function of an underlying topology. The call operator will return the vertex
  * of the graph whose position value is likely to be closest to a given position value. This
  * algorithm is approximate. It will select a M vertices from the graph from which it starts
  * a best-only search, where M is obtained as M = number_of_vertices / m_vertex_num_divider.
  * \tparam CompareFunction The functor type that can compare two distance measures (strict weak-ordering).
  */
template <typename CompareFunction = std::less<>>
struct best_only_neighbor_search {

  unsigned int m_vertex_num_divider;
  CompareFunction m_compare;
  /**
    * Default constructor.
    * \param aVertexNumDivider The division factor (should be greater than 1) which determines the
    *        fraction of the total number of vertices that is used to stem the best-only searches.
    *        Typical values are between 4 and 10.
    * \param compare The comparison functor for ordering the distances (strict weak ordering).
    */
  explicit best_only_neighbor_search(
      unsigned int aVertexNumDivider,
      CompareFunction compare = CompareFunction())
      : m_vertex_num_divider(aVertexNumDivider), m_compare(compare) {}

  best_only_neighbor_search() : best_only_neighbor_search(10) {}

  /**
    * This function template computes the topological distance between a position and the position of a
    * vertex of a graph. This function is used as a helper to the call-operator overloads.
    * \tparam Vertex The vertex descriptor type.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated with each vertex.
    * \param p A position in the space.
    * \param u A vertex which has a position associated to it, via the position property-map.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    */
  template <typename Vertex, typename Topology, typename PositionMap>
  double distance(const graph::property_value_t<PositionMap>& p, Vertex u,
                  const Topology& space, PositionMap position) const {
    return get(distance_metric, space)(p, get(position, u), space);
  }

  template <typename Graph, typename Topology, typename PositionMap>
  void search(const graph::property_value_t<PositionMap>& p,
              graph::graph_vertex_t<Graph>& u, double& d_min, Graph& g,
              const Topology& space, PositionMap position) {
    using Vertex = graph::graph_vertex_t<Graph>;
    d_min = distance(p, u, space, position);
    while (out_degree(u, g)) {
      Vertex v_min = u;
      for (auto [ei, ei_end] = out_edges(u, g); ei != ei_end; ++ei) {
        Vertex v = target(*ei, g);
        double d_v = distance(p, v, space, position);
        if (m_compare(d_v, d_min)) {
          d_min = d_v;
          v_min = v;
        }
      }
      if (v_min == u) {
        return;
      }
      u = v_min;
    }
    return;
  }

  /**
    * This call-operator finds the nearest vertex of a graph, to a given position.
    * \tparam Graph The graph type which can contain the vertices, should
    *         model boost::VertexListGraphConcept and boost::IncidenceGraphConcept.
    * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
    * \tparam PositionMap The property-map type which can store the position associated with each vertex.
    * \param p A position in the space, to which the nearest-neighbor is sought.
    * \param g A graph containing the vertices from which to find the nearest-neighbor,
    *        should be tree-structured.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    */
  template <typename Graph, typename Topology, typename PositionMap>
  graph::graph_vertex_t<Graph> operator()(
      const graph::property_value_t<PositionMap>& p, Graph& g,
      const Topology& space, PositionMap position) {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    using Vertex = graph::graph_vertex_t<Graph>;
    if (m_vertex_num_divider == 0) {
      m_vertex_num_divider = 1;
    }
    Vertex u_min = vertex(std::rand() % num_vertices(g), g);
    double d_min;
    search(p, u_min, d_min, g, space, position);
    for (int i = 0; i < num_vertices(g) / m_vertex_num_divider; ++i) {
      double d_v;
      Vertex v = vertex(std::rand() % num_vertices(g), g);
      search(p, v, d_v, g, space, position);
      if (m_compare(d_v, d_min)) {
        d_min = d_v;
        u_min = v;
      }
    }
    return u_min;
  }

  template <typename Graph, typename Topology, typename PositionMap,
            typename PriorityQueue, typename PriorityCompare>
  void search(const graph::property_value_t<PositionMap>& p,
              graph::graph_vertex_t<Graph> u, PriorityQueue& output,
              PriorityCompare p_compare, double d_min, Graph& g,
              const Topology& space, PositionMap position,
              std::size_t max_neighbors, double& radius) {
    using Vertex = graph::graph_vertex_t<Graph>;
    if (m_compare(d_min, radius)) {
      output.push_back(std::pair<double, Vertex>(d_min, u));
      std::push_heap(output.begin(), output.end(), p_compare);
      if (output.size() > max_neighbors) {
        std::pop_heap(output.begin(), output.end(), p_compare);
        output.pop_back();
        radius = output.front().first;
      }
    }
    for (auto [ei, ei_end] = out_edges(u, g); ei != ei_end; ++ei) {
      Vertex v = target(*ei, g);
      double d_v = distance(p, v, space, position);
      if (m_compare(d_v, d_min))
        search(p, v, output, p_compare, d_v, g, space, position, max_neighbors,
               radius);
    }
  }

  /**
    * This call-operator finds the nearest vertices of a graph, to a given position.
    * \tparam Graph The graph type which can contain the vertices, should
    *         model boost::VertexListGraphConcept and boost::IncidenceGraphConcept.
    * \tparam Topology The topology type which contains the positions.
    * \tparam PositionMap The property-map type which can store the position associated
    *         with each vertex.
    * \tparam OutputIterator The forward- output-iterator type which can contain the
    *         list of nearest-neighbors.
    * \param p A position in the space, to which the nearest-neighbors are sought.
    * \param output_first An iterator to the first place where to put the sorted list of
    *        elements with the smallest distance.
    * \param g A graph containing the vertices from which to find the nearest-neighbors,
    *        should be tree-structured.
    * \param space The topology objects which define the space in which the positions reside.
    * \param position The property-map which can retrieve the position associated to each vertex.
    * \param max_neighbors The maximum number of neighbors to have in the list.
    * \param radius The minimum distance around the position that a vertex should be in to be
    *        considered a neighbor.
    * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
    */
  template <typename Graph, typename Topology, typename PositionMap,
            typename OutputIterator>
  OutputIterator operator()(
      const graph::property_value_t<PositionMap>& p,
      OutputIterator output_first, Graph& g, const Topology& space,
      PositionMap position, std::size_t max_neighbors = 1,
      double radius = std::numeric_limits<double>::infinity()) {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    using Vertex = graph::graph_vertex_t<Graph>;
    const auto p_compare = [&](const auto& lhs, const auto& rhs) {
      return m_compare(lhs.first, rhs.first);
    };
    std::vector<std::pair<double, Vertex>> output;
    if (m_vertex_num_divider == 0) {
      m_vertex_num_divider = 1;
    }
    for (int i = 0; i < num_vertices(g) / m_vertex_num_divider; ++i) {
      Vertex v = vertex(std::rand() % num_vertices(g), g);
      double d_v = distance(p, v, space, position);
      search(p, v, output, p_compare, d_v, g, space, position, max_neighbors,
             radius);
    }
    std::sort_heap(output.begin(), output.end(), p_compare);
    for (auto [d, v] : output) {
      RK_UNUSED(d);
      *(output_first++) = v;
    }
    return output_first;
  }
};

}  // namespace ReaK::pp

#endif
