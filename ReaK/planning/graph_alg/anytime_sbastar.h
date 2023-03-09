/**
 * \file anytime_sbastar.h
 *
 * This library provides function templates and concepts that implement an Anytime Sampling-based A* search
 * algorithm. A ASBA* uses the Anytime A* search algorithm to drive the expansion of a roadmap into the free-space
 * in order to connect a start and goal location. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc.
 * All these customization points are left to the user to implement, some are defined by the
 * ASBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The ASBA* algorithm is a generalization of the Anytime A* algorithm where the neighborhood of a given node of
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an ASBA* algorithm, the same
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the
 * "surprisal") from drawing a new sample in the neighborhood.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_PLANNING_GRAPH_ALG_ANYTIME_SBASTAR_H_
#define REAK_PLANNING_GRAPH_ALG_ANYTIME_SBASTAR_H_

#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include "boost/graph/detail/d_ary_heap.hpp"
#include "boost/graph/graph_concepts.hpp"
#include "boost/graph/properties.hpp"

// BGL-Extra includes:
#include "boost/graph/more_property_maps.hpp"
#include "boost/graph/more_property_tags.hpp"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/sbastar_rrtstar.h"
#include "ReaK/planning/graph_alg/sbastar_search.h"
#include "ReaK/planning/graph_alg/simple_graph_traits.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the ASBA* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the ASBA* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the ASBA* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model SBAStarVisitorConcept, and AnytimeHeuristicVisitorConcept.
  *
  * \tparam Visitor The visitor class to be tested for modeling an ASBA* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct ASBAStarVisitorConcept : SBAStarVisitorConcept<Visitor, Graph, Topology>,
                                AnytimeHeuristicVisitorConcept<Visitor, Graph> {
};

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Topology>
struct asbastar_visitor_archetype : sbastar_visitor_archetype<Topology>,
                                    anytime_heuristic_visitor_archetype {};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the bi-directional ASBA* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the bi-directional ASBA* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the bi-directional ASBA* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model ASBAStarVisitorConcept, and NodeBackPushingVisitorConcept.
  *
  * \tparam Visitor The visitor class to be tested for modeling an ASBA* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct ASBAStarBidirVisitorConcept
    : ASBAStarVisitorConcept<Visitor, Graph, Topology> {

  BOOST_CONCEPT_ASSERT(
      (NodeBackPushingVisitorConcept<Visitor, Graph, Topology>));

  BOOST_CONCEPT_USAGE(ASBAStarBidirVisitorConcept) {}
};

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Topology>
struct asbastar_bidir_visitor_archetype
    : asbastar_visitor_archetype<Topology>,
      node_back_pushing_visitor_archetype<Topology> {};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the ASBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the ASBA*-RRT* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the ASBA*-RRT* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model SBARRTStarVisitorConcept, and AnytimeHeuristicVisitorConcept.
  *
  * \tparam Visitor The visitor class to be tested for modeling an ASBA*-RRT* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct ASBARRTStarVisitorConcept
    : SBARRTStarVisitorConcept<Visitor, Graph, Topology>,
      AnytimeHeuristicVisitorConcept<Visitor, Graph> {};

/**
  * This class is simply an archetype visitor for the ASBA*-RRT* algorithm.
  */
template <typename Topology>
struct asbarrtstar_visitor_archetype : sbarrtstar_visitor_archetype<Topology>,
                                       anytime_heuristic_visitor_archetype {};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the bi-directional ASBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the bi-directional ASBA*-RRT* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the bi-directional ASBA*-RRT* algorithm is in control of execution, but custom
  *behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model ASBARRTStarVisitorConcept, NodeBackPullingVisitorConcept and
  *NodeBackPushingVisitorConcept.
  *
  * \tparam Visitor The visitor class to be tested for modeling an ASBA* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct ASBARRTStarBidirVisitorConcept
    : ASBARRTStarVisitorConcept<Visitor, Graph, Topology> {

  BOOST_CONCEPT_ASSERT(
      (NodeBackPushingVisitorConcept<Visitor, Graph, Topology>));
  BOOST_CONCEPT_ASSERT(
      (NodeBackPullingVisitorConcept<Visitor, Graph, Topology>));

  BOOST_CONCEPT_USAGE(ASBARRTStarBidirVisitorConcept) {}
};

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Topology>
struct asbarrtstar_bidir_visitor_archetype
    : asbarrtstar_visitor_archetype<Topology>,
      node_back_pulling_visitor_archetype,
      node_back_pushing_visitor_archetype<Topology> {};

namespace detail {
namespace {

template <typename Graph, typename UniformCostVisitor, typename UpdatableQueue,
          typename IndexInHeapMap, typename KeyMap, typename PositionMap,
          typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap,
          typename SuccessorMap = null_vertex_prop_map<Graph>>
struct anytime_sbastar_bfs_visitor
    : sbastar_bfs_visitor<Graph, UniformCostVisitor, UpdatableQueue,
                          IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                          DensityMap, ConstrictionMap, DistanceMap,
                          PredecessorMap, FwdDistanceMap, SuccessorMap> {
  using base_type =
      sbastar_bfs_visitor<Graph, UniformCostVisitor, UpdatableQueue,
                          IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                          DensityMap, ConstrictionMap, DistanceMap,
                          PredecessorMap, FwdDistanceMap, SuccessorMap>;

  using Vertex = graph_vertex_t<Graph>;

  anytime_sbastar_bfs_visitor(double current_relax, UniformCostVisitor vis,
                              UpdatableQueue& Q, IndexInHeapMap index_in_heap,
                              KeyMap key, PositionMap pos, WeightMap weight,
                              DensityMap density, ConstrictionMap constriction,
                              DistanceMap dist, PredecessorMap pred,
                              FwdDistanceMap fwd_dist,
                              SuccessorMap succ = SuccessorMap())
      : base_type(vis, Q, index_in_heap, key, pos, weight, density,
                  constriction, dist, pred, fwd_dist, succ),
        m_current_relaxation(current_relax) {}

  void update_key(Vertex u, Graph& g) const {
    this->m_vis.affected_vertex(u, g);
    double g_u = get(this->m_distance, g[u]);
    double h_u = get(this->m_fwd_distance, g[u]);
    // Key-value for the min-heap (priority-queue):
    if (get(this->m_successor, g[u]) ==
        boost::graph_traits<Graph>::null_vertex()) {
      put(this->m_key, u,
          ((g_u + h_u) / (1.0 - get(this->m_constriction, g[u]))) /
                  (1.0 - get(this->m_density, g[u])) +
              m_current_relaxation * h_u);
    } else {
      put(this->m_key, u,
          ((g_u + h_u) / (1.0 - get(this->m_constriction, g[u]))) /
                  (1.0 - get(this->m_density, g[u])) +
              m_current_relaxation * g_u);
    }
  }

  void requeue_vertex(Vertex u, Graph& g) const {
    update_key(u, g);
    if (!this->m_vis.should_close(u, g)) {
      this->m_Q.push_or_update(u);
      this->m_vis.discover_vertex(u, g);
    }
  }
  void affected_vertex(Vertex u, Graph& g) const {
    requeue_vertex(u, g);
  }  // same function, different name.

  void update_relaxation(Graph& g) {
    m_current_relaxation =
        this->m_vis.adjust_relaxation(m_current_relaxation, g);

    for (auto [vi, vi_end] = vertices(g); vi != vi_end; ++vi) {
      requeue_vertex(*vi, g);
    }
  }

  void publish_path(Graph& g) {
    this->m_vis.publish_path(g);
    update_relaxation(g);
  }

  double m_current_relaxation;
};

template <typename Graph, typename Topology, typename ASBAStarVisitor,
          typename NodeConnector, typename KeyMap, typename PositionMap,
          typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename NcSelector>
void generate_anytime_sbastar_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex, const Topology& super_space,
    ASBAStarVisitor vis, NodeConnector connect_vertex, KeyMap key,
    PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    NcSelector select_neighborhood, double init_relaxation) {
  using Vertex = graph_vertex_t<Graph>;

  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap,
                 KeyCompareType());  // priority queue holding the OPEN set.

  anytime_sbastar_bfs_visitor<Graph, ASBAStarVisitor, MutableQueue,
                              IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                              DensityMap, ConstrictionMap, DistanceMap,
                              PredecessorMap, FwdDistanceMap>
      sba_bfs_vis(init_relaxation, vis, Q, index_in_heap, key, position, weight,
                  density, constriction, distance, predecessor, fwd_distance);

  put(distance, g[start_vertex], 0.0);
  put(predecessor, g[start_vertex], start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbastar_search_loop(g, super_space, sba_bfs_vis, connect_vertex,
                      sba_node_generator(), Q, select_neighborhood);
}

template <typename Graph, typename Topology, typename SBAStarVisitor,
          typename NodeConnector, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap, typename KeyMap, typename NcSelector>
void generate_anytime_sbastar_bidir_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex,
    graph_vertex_t<Graph> goal_vertex, const Topology& super_space,
    SBAStarVisitor vis, NodeConnector connect_vertex, KeyMap key,
    PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    SuccessorMap successor, NcSelector select_neighborhood,
    double init_relaxation) {
  using Vertex = graph_vertex_t<Graph>;

  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap,
                 KeyCompareType());  // priority queue holding the OPEN set.

  anytime_sbastar_bfs_visitor<Graph, SBAStarVisitor, MutableQueue,
                              IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                              DensityMap, ConstrictionMap, DistanceMap,
                              PredecessorMap, FwdDistanceMap, SuccessorMap>
      sba_bfs_vis(init_relaxation, vis, Q, index_in_heap, key, position, weight,
                  density, constriction, distance, predecessor, fwd_distance,
                  successor);

  put(distance, g[start_vertex], 0.0);
  put(predecessor, g[start_vertex], start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);
  if (goal_vertex != boost::graph_traits<Graph>::null_vertex()) {
    put(fwd_distance, g[goal_vertex], 0.0);
    put(successor, g[goal_vertex], goal_vertex);
    sba_bfs_vis.requeue_vertex(goal_vertex, g);
  }

  sbastar_bidir_loop(g, super_space, sba_bfs_vis, connect_vertex,
                     sba_bidir_node_generator(), Q, select_neighborhood);
}

template <typename Graph, typename Topology, typename ASBARRTStarVisitor,
          typename NodeConnector, typename KeyMap, typename PositionMap,
          typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename RandomSampler, typename NcSelector>
void generate_anytime_sbarrtstar_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex, const Topology& super_space,
    ASBARRTStarVisitor vis, NodeConnector connect_vertex, KeyMap key,
    PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    RandomSampler get_sample, NcSelector select_neighborhood,
    double init_relaxation, double SA_init_temperature = 0.0) {
  using Vertex = graph_vertex_t<Graph>;

  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap,
                 KeyCompareType());  // priority queue holding the OPEN set.

  anytime_sbastar_bfs_visitor<Graph, ASBARRTStarVisitor, MutableQueue,
                              IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                              DensityMap, ConstrictionMap, DistanceMap,
                              PredecessorMap, FwdDistanceMap>
      sba_bfs_vis(init_relaxation, vis, Q, index_in_heap, key, position, weight,
                  density, constriction, distance, predecessor, fwd_distance);

  put(distance, g[start_vertex], 0.0);
  put(predecessor, g[start_vertex], start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbarrtstar_search_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_node_generator(),
      rrg_node_generator<Topology, RandomSampler, NcSelector>(
          &super_space, get_sample, select_neighborhood),
      Q, select_neighborhood, SA_init_temperature);
}

template <typename Graph, typename Topology, typename ASBARRTStarVisitor,
          typename NodeConnector, typename KeyMap, typename PositionMap,
          typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename SuccessorMap,
          typename RandomSampler, typename NcSelector>
void generate_anytime_sbarrtstar_bidir_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex,
    graph_vertex_t<Graph> goal_vertex, const Topology& super_space,
    ASBARRTStarVisitor vis, NodeConnector connect_vertex, KeyMap key,
    PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    SuccessorMap successor, RandomSampler get_sample,
    NcSelector select_neighborhood, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Vertex = graph_vertex_t<Graph>;

  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap,
                 KeyCompareType());  // priority queue holding the OPEN set.

  anytime_sbastar_bfs_visitor<Graph, ASBARRTStarVisitor, MutableQueue,
                              IndexInHeapMap, KeyMap, PositionMap, WeightMap,
                              DensityMap, ConstrictionMap, DistanceMap,
                              PredecessorMap, FwdDistanceMap, SuccessorMap>
      sba_bfs_vis(init_relaxation, vis, Q, index_in_heap, key, position, weight,
                  density, constriction, distance, predecessor, fwd_distance,
                  successor);

  put(distance, g[start_vertex], 0.0);
  put(predecessor, g[start_vertex], start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);
  if (goal_vertex != boost::graph_traits<Graph>::null_vertex()) {
    put(fwd_distance, g[goal_vertex], 0.0);
    put(successor, g[goal_vertex], goal_vertex);
    sba_bfs_vis.requeue_vertex(goal_vertex, g);
  }

  sbarrtstar_bidir_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_bidir_node_generator(),
      rrg_bidir_generator<Topology, RandomSampler, NcSelector, PredecessorMap,
                          SuccessorMap>(&super_space, get_sample,
                                        select_neighborhood, predecessor,
                                        successor),
      Q, select_neighborhood, SA_init_temperature);
}

}  // namespace
}  // namespace detail

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_sbastar_no_init(const SBAStarBundle& bdl,
                                      double init_relaxation) {
  detail::generate_anytime_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, bdl.m_select_neighborhood, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_sbastar(const SBAStarBundle& bdl,
                              double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_sbastar_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_sbastar_bidir_no_init(const SBAStarBundle& bdl,
                                            double init_relaxation) {
  detail::generate_anytime_sbastar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, pruned_node_connector(), bdl.m_key, bdl.m_position,
      bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
      bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor,
      bdl.m_select_neighborhood, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_sbastar_bidir(const SBAStarBundle& bdl,
                                    double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_sbastar_bidir_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_sbastar_no_init(const SBAStarBundle& bdl,
                                           double init_relaxation) {
  detail::generate_anytime_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, bdl.m_select_neighborhood, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_sbastar(const SBAStarBundle& bdl,
                                   double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_lazy_sbastar_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_sbastar_bidir_no_init(const SBAStarBundle& bdl,
                                                 double init_relaxation) {
  detail::generate_anytime_sbastar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, bdl.m_successor, bdl.m_select_neighborhood,
      init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the
 * search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_sbastar_bidir(const SBAStarBundle& bdl,
                                         double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_lazy_sbastar_bidir_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_bnb_sbastar_no_init(const SBAStarBundle& bdl,
                                               double init_relaxation) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_anytime_sbastar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, bdl.m_select_neighborhood, init_relaxation);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_anytime_sbastar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, bdl.m_select_neighborhood, init_relaxation);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_bnb_sbastar(const SBAStarBundle& bdl,
                                       double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_lazy_bnb_sbastar_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_bnb_sbastar_bidir_no_init(const SBAStarBundle& bdl,
                                                     double init_relaxation) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_anytime_sbastar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor,
        bdl.m_select_neighborhood, init_relaxation);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_anytime_sbastar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor,
        bdl.m_select_neighborhood, init_relaxation);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the
 * search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  */
template <typename SBAStarBundle>
void generate_anytime_lazy_bnb_sbastar_bidir(const SBAStarBundle& bdl,
                                             double init_relaxation) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_lazy_bnb_sbastar_bidir_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_sbarrtstar_no_init(const SBAStarBundle& bdl,
                                         RandomSampler get_sample,
                                         double init_relaxation,
                                         double SA_init_temperature = 0.0) {
  detail::generate_anytime_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
      init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_sbarrtstar(const SBAStarBundle& bdl,
                                 RandomSampler get_sample,
                                 double init_relaxation,
                                 double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                      SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, pruned_node_connector(), bdl.m_key, bdl.m_position,
      bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
      bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
      bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the
 * search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                       RandomSampler get_sample,
                                       double init_relaxation,
                                       double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_sbarrtstar_bidir_no_init(bdl, get_sample, init_relaxation,
                                            SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_sbarrtstar_no_init(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::generate_anytime_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
      init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_sbarrtstar(const SBAStarBundle& bdl,
                                      RandomSampler get_sample,
                                      double init_relaxation,
                                      double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_lazy_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                           SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, bdl.m_successor, get_sample,
      bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start
 * the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                            RandomSampler get_sample,
                                            double init_relaxation,
                                            double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_lazy_sbarrtstar_bidir_no_init(
      bdl, get_sample, init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_bnb_sbarrtstar_no_init(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_anytime_sbarrtstar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
        init_relaxation, SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_anytime_sbarrtstar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
        init_relaxation, SA_init_temperature);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_bnb_sbarrtstar(const SBAStarBundle& bdl,
                                          RandomSampler get_sample,
                                          double init_relaxation,
                                          double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);
  generate_anytime_lazy_bnb_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                               SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* (with branch-and-bound) algorithm, without initialization of the
 * existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_bnb_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
        bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
        bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* (with branch-and-bound) algorithm, with initialization of the
 * existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_anytime_lazy_bnb_sbarrtstar_bidir(
    const SBAStarBundle& bdl, RandomSampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<Visitor, Graph, Space>));
  BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<Visitor, Graph, Space>));

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_anytime_lazy_bnb_sbarrtstar_bidir_no_init(
      bdl, get_sample, init_relaxation, SA_init_temperature);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_ANYTIME_SBASTAR_H_
