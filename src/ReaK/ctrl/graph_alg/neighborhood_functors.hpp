/**
 * \file neighborhood_functors.hpp
 *
 * This library provides functors to wrap K-NN search functors with a certain strategy to set the 
 * number of neighbors and the search radius. In particular, the star_neighborhood functor can be used 
 * to form probabilistically optimal algorithms like PRM*, RRG, and RRT*. The single_neighbor functor
 * can be used to wrap a general K-NN search functor into a single nearest-neighbor search.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
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

#ifndef REAK_NEIGHBORHOOD_FUNCTORS_HPP
#define REAK_NEIGHBORHOOD_FUNCTORS_HPP


#include <ReaK/core/base/misc_math.hpp>

namespace ReaK {

namespace graph {

  
  
  /**
   * This functor template applied the "star" neighborhood function to the nearest neighbor queries.
   * What is called the "star" neighborhood here is the radius and k-nn value that are computed from
   * the formulas proposed by Karaman and Frazzoli (Int. J. of Rob. Research., 2011). We call it the 
   * "star" neighborhood because it is the neighborhood definition that transforms the PRM algorithm 
   * into the PRM* algorithm, the RRT into the RRG algorithm, and plays a central part in the RRT* 
   * algorithm as well.
   * \tparam NNFinder The functor type of the underlying K-nearest neighbor finder (or predecessor / successor finder).
   */
  template <typename NNFinder>
  struct star_neighborhood {
    
    NNFinder find_neighbors;
    double c_space_dimensions;
    double gamma_value;
    
    /**
     * Parametrized constructor.
     * \param aFindNeighbors The functor to perform the underlying K-nearest neighbor search (or pred / succ search).
     * \param aCSpaceDimensions The number of dimensions of the free-space (e.g., 3 for R^3).
     * \param aGammaValue The gamma factor, should be greater than 2 * (1 + 1/d)^(1/d) * (volume(free_space) / volume(unit_ball))^(1/d).
     *                    A conservative estimate of the bound is 3 * L_c, where L_c is the characteristic length of the configuration 
     *                    space (with respect to the metric used for distances in the NN search, of course).
     */
    star_neighborhood(NNFinder aFindNeighbors, 
                      double aCSpaceDimensions, 
                      double aGammaValue) : 
                      find_neighbors(aFindNeighbors), 
                      c_space_dimensions(aCSpaceDimensions),
                      gamma_value(aGammaValue) { };
    
    /**
     * This function fills the output iterator (like a back-inserter) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map.
     * \tparam OutputIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param output_first An iterator to the start of the storage of the resulting neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename OutputIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    OutputIterator output_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      using std::pow; 
      std::size_t N = num_vertices(g);
      std::size_t log_N = math::highest_set_bit(N) + 1;
      find_neighbors(p, output_first, g, free_space, position, 4 * log_N, gamma_value * pow(log_N / double(N), 1.0 / c_space_dimensions));
    };
    
    /**
     * This function fills the output iterators (like a back-inserters) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map. This overload is expected to 
     * fill a neighborhood of possible predecessors and possible successors to a given position. This type 
     * of query is relevant for asymmetric problems in which nodes can rarely be connected with bi-directional paths.
     * \tparam PredIterator A forward- and output-iterator type.
     * \tparam SuccIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param pred_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param succ_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename PredIterator, 
              typename SuccIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    PredIterator pred_first, SuccIterator succ_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      using std::pow; using std::log2;
      std::size_t N = num_vertices(g);
      std::size_t log_N = math::highest_set_bit(N) + 1;
      find_neighbors(p, pred_first, succ_first, g, free_space, position, 4 * log_N, gamma_value * pow(log_N / double(N), 1.0 / c_space_dimensions));
    };
  };
    
  
  
  
  
  
  /**
   * This functor template applies the single-neighbor neighborhood function to the nearest neighbor queries.
   * This is sort of a trivial functor, but can be used to wrap a k-nn finder.
   * \tparam NNFinder The functor type of the underlying K-nearest neighbor finder (or predecessor / successor finder).
   */
  template <typename NNFinder>
  struct single_neighbor {
    
    NNFinder find_neighbors;
    
    /**
     * Parametrized constructor.
     * \param aFindNeighbors The functor to perform the underlying K-nearest neighbor search (or pred / succ search).
     */
    single_neighbor(NNFinder aFindNeighbors) : find_neighbors(aFindNeighbors) { };
    
    /**
     * This function fills the output iterator (like a back-inserter) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map.
     * \tparam OutputIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param output_first An iterator to the start of the storage of the resulting neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename OutputIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    OutputIterator output_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      find_neighbors(p, output_first, g, free_space, position, 1, std::numeric_limits< double >::infinity());
    };
    
    /**
     * This function fills the output iterators (like a back-inserters) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map. This overload is expected to 
     * fill a neighborhood of possible predecessors and possible successors to a given position. This type 
     * of query is relevant for asymmetric problems in which nodes can rarely be connected with bi-directional paths.
     * \tparam PredIterator A forward- and output-iterator type.
     * \tparam SuccIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param pred_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param succ_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename PredIterator, 
              typename SuccIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    PredIterator pred_first, SuccIterator succ_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      find_neighbors(p, pred_first, succ_first, g, free_space, position, 1, std::numeric_limits< double >::infinity());
    };
  };
    
    
    
  
  
  
  
  /**
   * This functor template applies the single-neighbor neighborhood function to the nearest neighbor queries.
   * This is sort of a trivial functor, but can be used to wrap a k-nn finder.
   * \tparam NNFinder The functor type of the underlying K-nearest neighbor finder (or predecessor / successor finder).
   */
  template <typename NNFinder>
  struct fixed_neighborhood {
    
    NNFinder find_neighbors;
    std::size_t max_neighbor_count;
    double max_radius;
    
    /**
     * Parametrized constructor.
     * \param aFindNeighbors The functor to perform the underlying K-nearest neighbor search (or pred / succ search).
     */
    fixed_neighborhood(NNFinder aFindNeighbors,
                       std::size_t aMaxNeighborCount = 1,
                       double aMaxRadius = std::numeric_limits<double>::infinity()) : 
                       find_neighbors(aFindNeighbors),
                       max_neighbor_count(aMaxNeighborCount),
                       max_radius(aMaxRadius) { };
    
    /**
     * This function fills the output iterator (like a back-inserter) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map.
     * \tparam OutputIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param output_first An iterator to the start of the storage of the resulting neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename OutputIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    OutputIterator output_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      find_neighbors(p, output_first, g, free_space, position, max_neighbor_count, max_radius);
    };
    
    /**
     * This function fills the output iterators (like a back-inserters) with the neighborhood of the given
     * position in the given graph, over the given topology and position-map. This overload is expected to 
     * fill a neighborhood of possible predecessors and possible successors to a given position. This type 
     * of query is relevant for asymmetric problems in which nodes can rarely be connected with bi-directional paths.
     * \tparam PredIterator A forward- and output-iterator type.
     * \tparam SuccIterator A forward- and output-iterator type.
     * \tparam Graph A VertexListGraph (as of BGL).
     * \tparam Topology A topology type (as of topology concepts in ReaK).
     * \tparam PositionMap A property-map type over the vertices of the graph, producing associated positions on the topology.
     * \param p The position with respect to which the nearest-neighbors are sought.
     * \param pred_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param succ_first An iterator to the start of the storage of the resulting predecessor neighborhood vertices, should allow for enough room for the resulting neighborhood.
     * \param g The graph on which to look for the nearest neighbors.
     * \param free_space The topology representing the free-space of the positions.
     * \param position The position-map object which can retrieve the positions associated to each vertex of the graph.
     */
    template <typename PredIterator, 
              typename SuccIterator, 
              typename Graph, 
              typename Topology, 
              typename PositionMap>
    void operator()(const typename boost::property_traits<PositionMap>::value_type& p,
                    PredIterator pred_first, SuccIterator succ_first,
                    Graph& g, const Topology& free_space, PositionMap position) const {
      find_neighbors(p, pred_first, succ_first, g, free_space, position, max_neighbor_count, max_radius);
    };
  };
    
    
    
  
    
  
  
  
  


};


};

#endif










