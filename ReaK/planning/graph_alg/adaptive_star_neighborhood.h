/**
 * \file adaptive_star_neighborhood.h
 *
 * This library provides a functor to wrap K-NN search functors with an adaptive star-neighborhood strategy
 * to set the number of neighbors and the search radius according to the current cloud of samples in the
 * motion-graph. The adaptive_star_nbhd functor can be used to form probabilistically optimal algorithms
 * like PRM*, RRG, and RRT*.
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

#ifndef REAK_PLANNING_GRAPH_ALG_ADAPTIVE_STAR_NEIGHBORHOOD_H_
#define REAK_PLANNING_GRAPH_ALG_ADAPTIVE_STAR_NEIGHBORHOOD_H_

#include "ReaK/core/base/misc_math.h"

#include <limits>
#include "bagl/property_map.h"

namespace ReaK::graph {

/**
 * This functor template applied the "star" neighborhood function to the nearest neighbor queries.
 * What is called the "star" neighborhood here is the radius and k-nn value that are computed from
 * the formulas proposed by Karaman and Frazzoli (Int. J. of Rob. Research., 2011). We call it the
 * "star" neighborhood because it is the neighborhood definition that transforms the PRM algorithm
 * into the PRM* algorithm, and plays a central part in the RRT* algorithm as well. This version
 * of the star neighborhood implements a scheme which tries to keep track of the characteristic
 * length of the cloud of samples.
 * \tparam NNFinder The functor type of the underlying K-nearest neighbor finder (or predecessor / successor finder).
 */
template <typename NNFinder>
struct adaptive_star_nbhd {

  NNFinder find_neighbors;
  double c_space_dimensions;
  mutable std::size_t next_logN_value = 1;
  mutable double gamma_value = std::numeric_limits<double>::infinity();

  /**
   * Parametrized constructor.
   * \param aFindNeighbors The functor to perform the underlying K-nearest neighbor search (or pred / succ search).
   * \param aCSpaceDimensions The number of dimensions of the free-space (e.g., 3 for R^3).
   */
  adaptive_star_nbhd(NNFinder aFindNeighbors, double aCSpaceDimensions)
      : find_neighbors(aFindNeighbors), c_space_dimensions(aCSpaceDimensions) {}

  /**
   * This function fills the output iterator (like a back-inserter) with the neighborhood of the given
   * position in the given graph, over the given topology and position-map.
   * \tparam OutputIterator A forward- and output-iterator type.
   * \tparam Graph A VertexListGraph (as of BGL).
   * \tparam Topology A topology type (as of topology concepts in ReaK).
   * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the
   * topology.
   * \param p The position with respect to which the nearest-neighbors are sought.
   * \param output_first An iterator to the start of the storage of the resulting neighborhood vertices, should allow
   * for enough room for the resulting neighborhood.
   * \param g The graph on which to look for the nearest neighbors.
   * \param free_space The topology representing the free-space of the positions.
   * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
   */
  template <typename OutputIterator, typename Graph, typename Topology,
            typename PositionMap>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  OutputIterator output_first, Graph& g,
                  const Topology& free_space, PositionMap position) const {
    using std::pow;
    std::size_t N = num_vertices(g);
    std::size_t log_N = math::highest_set_bit(N) + 1;
    if (log_N > next_logN_value) {
      next_logN_value = log_N;
      gamma_value =
          3.0 * find_neighbors.get_characteristic_size(g, free_space, position);
    }
    find_neighbors(
        p, output_first, g, free_space, position, 4 * log_N,
        gamma_value * pow(log_N / double(N), 1.0 / c_space_dimensions));
  }

  /**
   * This function fills the output iterators (like a back-inserters) with the neighborhood of the given
   * position in the given graph, over the given topology and position-map. This overload is expected to
   * fill a neighborhood of possible predecessors and possible successors to a given position. This type
   * of query is relevant for asymmetric problems in which nodes can rarely be connected with bi-directional paths.
   * \tparam PredIterator A forward- and output-iterator type.
   * \tparam SuccIterator A forward- and output-iterator type.
   * \tparam Graph A VertexListGraph (as of BGL).
   * \tparam Topology A topology type (as of topology concepts in ReaK).
   * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the
   * topology.
   * \param p The position with respect to which the nearest-neighbors are sought.
   * \param pred_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices,
   * should allow for enough room for the resulting neighborhood.
   * \param succ_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices,
   * should allow for enough room for the resulting neighborhood.
   * \param g The graph on which to look for the nearest neighbors.
   * \param free_space The topology representing the free-space of the positions.
   * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
   */
  template <typename PredIterator, typename SuccIterator, typename Graph,
            typename Topology, typename PositionMap>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  PredIterator pred_first, SuccIterator succ_first, Graph& g,
                  const Topology& free_space, PositionMap position) const {
    using std::log2;
    using std::pow;
    std::size_t N = num_vertices(g);
    std::size_t log_N = math::highest_set_bit(N) + 1;
    if (log_N > next_logN_value) {
      next_logN_value = log_N;
      gamma_value =
          3.0 * find_neighbors.get_characteristic_size(g, free_space, position);
    }
    find_neighbors(
        p, pred_first, succ_first, g, free_space, position, 4 * log_N,
        gamma_value * pow(log_N / double(N), 1.0 / c_space_dimensions));
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_ADAPTIVE_STAR_NEIGHBORHOOD_H_
