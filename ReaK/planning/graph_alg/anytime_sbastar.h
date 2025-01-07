/**
 * \file anytime_sbastar.h
 *
 * This library provides function templates and concepts that implement an Anytime Sampling-based A* search
 * algorithm. A ASBA* uses the Anytime A* search algorithm to drive the expansion of a roadmap into the free-space
 * in order to connect a start and goal location. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc.
 * All these customization points are left to the user to implement, some are defined by the
 * ASBAStarVisitor (random-walk, edge-added, etc.).
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
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/sbastar_rrtstar.h"
#include "ReaK/planning/graph_alg/sbastar_search.h"

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
  * The visitor class should model SBAStarVisitor, and AnytimeHeuristicVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept ASBAStarVisitor = SBAStarVisitor<Visitor, Graph, Space>&&
    AnytimeHeuristicVisitor<Visitor, Graph>;

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Space>
struct asbastar_visitor_archetype : sbastar_visitor_archetype<Space>,
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
  * The visitor class should model ASBAStarVisitor, and NodeBackPushingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept ASBAStarBidirVisitor = ASBAStarVisitor<Visitor, Graph, Space>&&
    NodeBackPushingVisitor<Visitor, Graph, Space>;

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Space>
struct asbastar_bidir_visitor_archetype
    : asbastar_visitor_archetype<Space>,
      node_back_pushing_visitor_archetype<Space> {};

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
  * The visitor class should model SBARRTStarVisitor, and AnytimeHeuristicVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept ASBARRTStarVisitor = SBARRTStarVisitor<Visitor, Graph, Space>&&
    AnytimeHeuristicVisitor<Visitor, Graph>;

/**
  * This class is simply an archetype visitor for the ASBA*-RRT* algorithm.
  */
template <typename Space>
struct asbarrtstar_visitor_archetype : sbarrtstar_visitor_archetype<Space>,
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
  * The visitor class should model ASBARRTStarVisitor, NodeBackPullingVisitor and NodeBackPushingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept ASBARRTStarBidirVisitor =
    ASBARRTStarVisitor<Visitor, Graph, Space>&& NodeBackPushingVisitor<
        Visitor, Graph, Space>&& NodeBackPullingVisitor<Visitor, Graph, Space>;

/**
  * This class is simply an archetype visitor for the ASBA* algorithm.
  */
template <typename Space>
struct asbarrtstar_bidir_visitor_archetype
    : asbarrtstar_visitor_archetype<Space>,
      node_back_pulling_visitor_archetype,
      node_back_pushing_visitor_archetype<Space> {};

namespace sbastar_detail {

template <typename Graph, typename UniformCostVisitor, typename UpdatableQueue,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> IndexInHeapMap,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap =
              null_vertex_prop_map<Graph>>
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

  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  anytime_sbastar_bfs_visitor(double current_relax, UniformCostVisitor vis,
                              UpdatableQueue& Q, IndexInHeapMap index_in_heap,
                              KeyMap key, PositionMap pos, WeightMap weight,
                              DensityMap density, ConstrictionMap constriction,
                              DistanceMap dist, PredecessorMap pred,
                              FwdDistanceMap fwd_dist,
                              SuccessorMap succ = SuccessorMap())
      : base_type(vis, Q, index_in_heap, key, pos, weight, density,
                  constriction, dist, pred, fwd_dist, succ),
        current_relaxation_(current_relax) {}

  void update_key(Vertex u, Graph& g) const {
    this->vis_.affected_vertex(u, g);
    double g_u = get(this->distance_, get_property(g, u));
    double h_u = get(this->fwd_distance_, get_property(g, u));
    // Key-value for the min-heap (priority-queue):
    if (get(this->successor_, get_property(g, u)) ==
        bagl::graph_traits<Graph>::null_vertex()) {
      put(this->key_, u,
          ((g_u + h_u) / (1.0 - get(this->constriction_, get_property(g, u)))) /
                  (1.0 - get(this->density_, get_property(g, u))) +
              current_relaxation_ * h_u);
    } else {
      put(this->key_, u,
          ((g_u + h_u) / (1.0 - get(this->constriction_, get_property(g, u)))) /
                  (1.0 - get(this->density_, get_property(g, u))) +
              current_relaxation_ * g_u);
    }
  }

  void requeue_vertex(Vertex u, Graph& g) const {
    update_key(u, g);
    if (!this->vis_.should_close(u, g)) {
      this->q_.push_or_update(u);
      this->vis_.discover_vertex(u, g);
    }
  }
  void affected_vertex(Vertex u, Graph& g) const {
    requeue_vertex(u, g);
  }  // same function, different name.

  void update_relaxation(Graph& g) {
    current_relaxation_ = this->vis_.adjust_relaxation(current_relaxation_, g);

    for (auto v : vertices(g)) {
      requeue_vertex(v, g);
    }
  }

  void publish_path(Graph& g) {
    this->vis_.publish_path(g);
    update_relaxation(g);
  }

  double current_relaxation_;
};

template <typename Graph, typename V, typename Q, typename... Maps>
anytime_sbastar_bfs_visitor<Graph, V, Q, Maps...>
make_anytime_sbastar_bfs_visitor(double current_relax, V v, Q& q, Maps... m) {
  return {current_relax, v, q, m...};
}

template <typename Graph, pp::MetricSpace Space, typename ASBAStarVisitor,
          typename NodeConnector,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          typename NcSelector>
void generate_anytime_sbastar_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    const Space& super_space, ASBAStarVisitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    NcSelector select_neighborhood, double init_relaxation) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_anytime_sbastar_bfs_visitor<Graph>(
      init_relaxation, vis, q, index_in_heap.ref(), key, position, weight,
      density, constriction, distance, predecessor, fwd_distance);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbastar_search_loop(g, super_space, sba_bfs_vis, connect_vertex,
                      sba_node_generator(), q, select_neighborhood);
}

template <typename Graph, pp::MetricSpace Space, typename SBAStarVisitor,
          typename NodeConnector,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap,
          typename NcSelector>
void generate_anytime_sbastar_bidir_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
    const Space& super_space, SBAStarVisitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    SuccessorMap successor, NcSelector select_neighborhood,
    double init_relaxation) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_anytime_sbastar_bfs_visitor<Graph>(
      init_relaxation, vis, q, index_in_heap.ref(), key, position, weight,
      density, constriction, distance, predecessor, fwd_distance, successor);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);
  if (goal_vertex != bagl::graph_traits<Graph>::null_vertex()) {
    put(fwd_distance, get_property(g, goal_vertex), 0.0);
    put(successor, get_property(g, goal_vertex), goal_vertex);
    sba_bfs_vis.requeue_vertex(goal_vertex, g);
  }

  sbastar_bidir_loop(g, super_space, sba_bfs_vis, connect_vertex,
                     sba_bidir_node_generator(), q, select_neighborhood);
}

template <typename Graph, pp::MetricSpace Space, typename ASBARRTStarVisitor,
          typename NodeConnector,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_anytime_sbarrtstar_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    const Space& super_space, ASBARRTStarVisitor vis,
    NodeConnector connect_vertex, KeyMap key, PositionMap position,
    WeightMap weight, DensityMap density, ConstrictionMap constriction,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, Sampler get_sample,
    NcSelector select_neighborhood, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_anytime_sbastar_bfs_visitor<Graph>(
      init_relaxation, vis, q, index_in_heap.ref(), key, position, weight,
      density, constriction, distance, predecessor, fwd_distance);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbarrtstar_search_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_node_generator(),
      rrg_node_generator(&super_space, get_sample, select_neighborhood), q,
      select_neighborhood, SA_init_temperature);
}

template <typename Graph, pp::MetricSpace Space, typename ASBARRTStarVisitor,
          typename NodeConnector,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_anytime_sbarrtstar_bidir_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
    const Space& super_space, ASBARRTStarVisitor vis,
    NodeConnector connect_vertex, KeyMap key, PositionMap position,
    WeightMap weight, DensityMap density, ConstrictionMap constriction,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, SuccessorMap successor, Sampler get_sample,
    NcSelector select_neighborhood, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_anytime_sbastar_bfs_visitor<Graph>(
      init_relaxation, vis, q, index_in_heap.ref(), key, position, weight,
      density, constriction, distance, predecessor, fwd_distance, successor);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);
  if (goal_vertex != bagl::graph_traits<Graph>::null_vertex()) {
    put(fwd_distance, get_property(g, goal_vertex), 0.0);
    put(successor, get_property(g, goal_vertex), goal_vertex);
    sba_bfs_vis.requeue_vertex(goal_vertex, g);
  }

  sbarrtstar_bidir_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_bidir_node_generator(),
      rrg_bidir_generator(&super_space, get_sample, select_neighborhood,
                          predecessor, successor),
      q, select_neighborhood, SA_init_temperature);
}

}  // namespace sbastar_detail

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
  sbastar_detail::generate_anytime_sbastar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.select_neighborhood_, init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
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
  sbastar_detail::generate_anytime_sbastar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.successor_, bdl.select_neighborhood_,
      init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
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
  sbastar_detail::generate_anytime_sbastar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.select_neighborhood_, init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
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
  sbastar_detail::generate_anytime_sbastar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.successor_, bdl.select_neighborhood_, init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
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
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);

  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_anytime_sbastar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.select_neighborhood_, init_relaxation);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_anytime_sbastar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.select_neighborhood_, init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
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
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_anytime_sbastar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.successor_, bdl.select_neighborhood_,
        init_relaxation);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_anytime_sbastar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_,
        bdl.weight_, bdl.density_, bdl.constriction_, bdl.distance_,
        bdl.predecessor_, bdl.fwd_distance_, bdl.successor_,
        bdl.select_neighborhood_, init_relaxation);
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
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
  generate_anytime_lazy_bnb_sbastar_bidir_no_init(bdl, init_relaxation);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_sbarrtstar_no_init(const SBAStarBundle& bdl,
                                         Sampler get_sample,
                                         double init_relaxation,
                                         double SA_init_temperature = 0.0) {
  sbastar_detail::generate_anytime_sbarrtstar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, get_sample, bdl.select_neighborhood_, init_relaxation,
      SA_init_temperature);
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_sbarrtstar(const SBAStarBundle& bdl, Sampler get_sample,
                                 double init_relaxation,
                                 double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_anytime_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                      SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  sbastar_detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.successor_, get_sample, bdl.select_neighborhood_,
      init_relaxation, SA_init_temperature);
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                       Sampler get_sample,
                                       double init_relaxation,
                                       double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
  generate_anytime_sbarrtstar_bidir_no_init(bdl, get_sample, init_relaxation,
                                            SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_sbarrtstar_no_init(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::generate_anytime_sbarrtstar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      get_sample, bdl.select_neighborhood_, init_relaxation,
      SA_init_temperature);
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_sbarrtstar(const SBAStarBundle& bdl,
                                      Sampler get_sample,
                                      double init_relaxation,
                                      double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_anytime_lazy_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                           SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.successor_, get_sample, bdl.select_neighborhood_, init_relaxation,
      SA_init_temperature);
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                            Sampler get_sample,
                                            double init_relaxation,
                                            double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
  generate_anytime_lazy_sbarrtstar_bidir_no_init(
      bdl, get_sample, init_relaxation, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_bnb_sbarrtstar_no_init(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_anytime_sbarrtstar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, get_sample, bdl.select_neighborhood_,
        init_relaxation, SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_anytime_sbarrtstar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, get_sample, bdl.select_neighborhood_,
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_bnb_sbarrtstar(const SBAStarBundle& bdl,
                                          Sampler get_sample,
                                          double init_relaxation,
                                          double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_anytime_lazy_bnb_sbarrtstar_no_init(bdl, get_sample, init_relaxation,
                                               SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional Anytime-Lazy-SBA*-RRT* (with branch-and-bound) algorithm, without initialization of the
 * existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
  *        Should be greater than 0, the recommeded value is 10.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_bnb_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.successor_, get_sample, bdl.select_neighborhood_,
        init_relaxation, SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_anytime_sbarrtstar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_,
        bdl.weight_, bdl.density_, bdl.constriction_, bdl.distance_,
        bdl.predecessor_, bdl.fwd_distance_, bdl.successor_, get_sample,
        bdl.select_neighborhood_, init_relaxation, SA_init_temperature);
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
template <typename SBAStarBundle, typename Sampler>
void generate_anytime_lazy_bnb_sbarrtstar_bidir(
    const SBAStarBundle& bdl, Sampler get_sample, double init_relaxation,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(ASBAStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
  generate_anytime_lazy_bnb_sbarrtstar_bidir_no_init(
      bdl, get_sample, init_relaxation, SA_init_temperature);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_ANYTIME_SBASTAR_H_
