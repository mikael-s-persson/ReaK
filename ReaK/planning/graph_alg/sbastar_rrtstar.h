/**
 * \file sbastar_rrtstar.h
 *
 * This library provides function templates and concepts that implement a Sampling-based A* search
 * algorithm with an RRT* exploratory phase. A SBA* uses the A* search algorithm to drive the expansion
 * of a roadmap into the free-space in order to connect a start and goal location. When useful nodes
 * are exhausted, an number of RRT* iterations are performed to use the Voronoi bias to generate
 * more useful nodes before continuing the SBA* iterations. This algorithm has many customization points
 * because there are many choices to be made in the method, such as how to find nearest neighbors for
 * attempting to connect them through free-space, how to expand vertices, when to stop the algorithm, etc.
 * All these customization points are left to the user to implement, some are defined by the
 * SBARRTStarVisitorConcept (random-walk, vertex-added, etc.).
 *
 * The SBA* algorithm is a generalization of the A* algorithm where the neighborhood of a given node of
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an SBA* algorithm, the same
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

#ifndef REAK_PLANNING_GRAPH_ALG_SBASTAR_RRTSTAR_H_
#define REAK_PLANNING_GRAPH_ALG_SBASTAR_RRTSTAR_H_

#include <cmath>
#include <limits>
#include <random>
#include <utility>

#include "ReaK/core/base/global_rng.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include <tuple>
#include "adj_list_tree_overlay.h"
#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/node_generators.h"
#include "ReaK/planning/graph_alg/sbastar_search.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the SBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the SBA*-RRT* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the SBA*-RRT* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * the visitor should model SBAStarVisitor and NodePullingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept SBARRTStarVisitor = SBAStarVisitor<Visitor, Graph, Space>&&
    NodePullingVisitor<Visitor, Graph, Space>;

/**
  * This class is simply an archetype visitor for the SBA*-RRT* algorithm.
  */
template <typename Space>
struct sbarrtstar_visitor_archetype : sbastar_visitor_archetype<Space>,
                                      node_pulling_visitor_archetype {};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the bi-directional SBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the bi-directional SBA*-RRT* algorithm.
  * In other words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the bi-directional SBA*-RRT* algorithm is in control of execution,
  * but custom behavior can be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model SBARRTStarVisitor, NodeBackPullingVisitor and NodeBackPushingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept SBARRTStarBidirVisitor =
    SBARRTStarVisitor<Visitor, Graph, Space>&& NodeBackPushingVisitor<
        Visitor, Graph, Space>&& NodeBackPullingVisitor<Visitor, Graph, Space>;

/**
 * This class is simply an archetype visitor for the bi-directional SBA*-RRT* algorithm.
 */
template <typename Space>
struct sbarrtstar_bidir_visitor_archetype
    : sbarrtstar_visitor_archetype<Space>,
      node_back_pulling_visitor_archetype,
      node_back_pushing_visitor_archetype<Space> {};

namespace sbastar_detail {

template <typename Graph, pp::MetricSpace Space, typename Visitor,
          typename MotionGraphConnector, typename SBANodeGenerator,
          typename RRTNodeGenerator, typename MutableQueue, typename NcSelector>
void sbarrtstar_search_loop(Graph& g, const Space& super_space,
                            Visitor& sba_vis,
                            MotionGraphConnector connect_vertex,
                            SBANodeGenerator sba_generate_node,
                            RRTNodeGenerator rrt_generate_node, MutableQueue& q,
                            NcSelector select_neighborhood,
                            double initial_temperature) {
  using std::exp;
  using std::log;
  std::uniform_real_distribution<double> unit_dist{};

  const auto connect_if_needed = [&](auto& x_near, auto& p_new, auto& eprop) {
    if ((x_near != bagl::graph_traits<Graph>::null_vertex()) &&
        (get(sba_vis.predecessor_, get_property(g, x_near)) !=
         bagl::graph_traits<Graph>::null_vertex())) {
      connect_vertex(p_new, x_near, eprop, g, super_space, sba_vis,
                     sba_vis.position_, sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.weight_, select_neighborhood);
    }
  };

  while (!q.empty() && sba_vis.keep_going()) {
    double entropy = 1.0 - exp(-initial_temperature /
                               log(static_cast<double>(num_vertices(g))));
    // generate random-number between 0 and 1.
    double rand_value = unit_dist(ReaK::get_global_rng());
    bool use_sba_sampling = (rand_value > entropy);

    if (use_sba_sampling) {
      auto u = q.top();
      q.pop();

      // stop if the best nodes do not meet the potential threshold.
      while (!sba_vis.has_search_potential(u, g)) {
        if (q.empty()) {
          u = bagl::graph_traits<Graph>::null_vertex();
          break;
        }
        u = q.top();
        q.pop();
      }
      if (u == bagl::graph_traits<Graph>::null_vertex()) {
        break;  // no more nodes with search potential.
      }

      sba_vis.examine_vertex(u, g);

      auto [x_near, p_new, eprop] =
          sba_generate_node(u, g, sba_vis, sba_vis.position_);

      // then push it back on the OPEN queue.
      if ((x_near != bagl::graph_traits<Graph>::null_vertex()) || q.empty()) {
        sba_vis.requeue_vertex(u, g);
      }

      connect_if_needed(x_near, p_new, eprop);
    } else {
      auto [x_near, p_new, eprop] =
          rrt_generate_node(g, sba_vis,
                            bagl::composite_property_map(
                                sba_vis.position_, get(bagl::vertex_all, g)));

      connect_if_needed(x_near, p_new, eprop);
    }
  }  // while
}

template <typename Graph, pp::MetricSpace Space, typename Visitor,
          typename MotionGraphConnector, typename SBANodeGenerator,
          typename RRTNodeGenerator, typename MutableQueue, typename NcSelector>
void sbarrtstar_bidir_loop(Graph& g, const Space& super_space, Visitor& sba_vis,
                           MotionGraphConnector connect_vertex,
                           SBANodeGenerator sba_generate_node,
                           RRTNodeGenerator rrt_generate_node, MutableQueue& q,
                           NcSelector select_neighborhood,
                           double initial_temperature) {
  using std::exp;
  using std::log;
  std::uniform_real_distribution<double> unit_dist{};

  const auto connect_if_needed = [&](auto& x_near_pred, auto& p_new_pred,
                                     auto& ep_pred, auto& x_near_succ,
                                     auto& p_new_succ, auto& ep_succ) {
    if (x_near_pred != bagl::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = bagl::graph_traits<Graph>::null_vertex();
      bagl::edge_bundle_type<Graph> ep_other;
      connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other,
                     g, super_space, sba_vis, sba_vis.position_,
                     sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.fwd_distance_, sba_vis.successor_, sba_vis.weight_,
                     select_neighborhood);
    }
    if (x_near_succ != bagl::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = bagl::graph_traits<Graph>::null_vertex();
      bagl::edge_bundle_type<Graph> ep_other;
      connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ,
                     g, super_space, sba_vis, sba_vis.position_,
                     sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.fwd_distance_, sba_vis.successor_, sba_vis.weight_,
                     select_neighborhood);
    }
  };

  while (!q.empty() && sba_vis.keep_going()) {
    double entropy = 1.0 - exp(-initial_temperature /
                               log(static_cast<double>(num_vertices(g))));
    // generate random-number between 0 and 1.
    double rand_value = unit_dist(ReaK::get_global_rng());
    bool use_sba_sampling = (rand_value > entropy);

    if (use_sba_sampling) {
      auto u = q.top();
      q.pop();

      // stop if the best nodes do not meet the potential threshold.
      while (!sba_vis.has_search_potential(u, g)) {
        if (q.empty()) {
          u = bagl::graph_traits<Graph>::null_vertex();
          break;
        }
        u = q.top();
        q.pop();
      }
      if (u == bagl::graph_traits<Graph>::null_vertex()) {
        break;  // no more nodes with search potential.
      }

      sba_vis.examine_vertex(u, g);

      auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ,
            ep_succ] = sba_generate_node(u, g, sba_vis, sba_vis.position_);

      // then push it back on the OPEN queue.
      if ((x_near_pred != bagl::graph_traits<Graph>::null_vertex()) ||
          (x_near_succ != bagl::graph_traits<Graph>::null_vertex()) ||
          (q.empty())) {
        sba_vis.requeue_vertex(u, g);
      }

      connect_if_needed(x_near_pred, p_new_pred, ep_pred, x_near_succ,
                        p_new_succ, ep_succ);
    } else {
      auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ,
            ep_succ] =
          rrt_generate_node(g, sba_vis,
                            bagl::composite_property_map(
                                sba_vis.position_, get(bagl::vertex_all, g)));

      connect_if_needed(x_near_pred, p_new_pred, ep_pred, x_near_succ,
                        p_new_succ, ep_succ);
    }
  }  // while
}

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBARRTStarVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename RandomSampler, typename NcSelector>
void generate_sbarrtstar_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    const Space& super_space, Visitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    RandomSampler get_sample, NcSelector select_neighborhood,
    double SA_init_temperature = 0.0) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_sbastar_bfs_visitor<Graph>(
      vis, q, index_in_heap.ref(), key, position, weight, density, constriction,
      distance, predecessor, fwd_distance);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbarrtstar_search_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_node_generator(),
      rrg_node_generator(&super_space, get_sample, select_neighborhood), q,
      select_neighborhood, SA_init_temperature);
}

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBARRTStarBidirVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap, typename RandomSampler, typename NcSelector>
void generate_sbarrtstar_bidir_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
    const Space& super_space, Visitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    SuccessorMap successor, RandomSampler get_sample,
    NcSelector select_neighborhood, double SA_init_temperature = 0.0) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_sbastar_bfs_visitor<Graph>(
      vis, q, index_in_heap.ref(), key, position, weight, density, constriction,
      distance, predecessor, fwd_distance, successor);

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
  * using the SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_no_init(const SBAStarBundle& bdl,
                                 RandomSampler get_sample,
                                 double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);

  sbastar_detail::generate_sbarrtstar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, get_sample, bdl.select_neighborhood_,
      SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar(const SBAStarBundle& bdl, RandomSampler get_sample,
                         double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);

  generate_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_bidir_no_init(const SBAStarBundle& bdl,
                                       RandomSampler get_sample,
                                       double SA_init_temperature = 0.0) {
  sbastar_detail::generate_sbarrtstar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.successor_, get_sample, bdl.select_neighborhood_,
      SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_bidir(const SBAStarBundle& bdl,
                               RandomSampler get_sample,
                               double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);

  generate_sbarrtstar_bidir_no_init(bdl, get_sample, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_sbarrtstar_no_init(const SBAStarBundle& bdl,
                                      Sampler get_sample,
                                      double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::generate_sbarrtstar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      get_sample, bdl.select_neighborhood_, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_sbarrtstar(const SBAStarBundle& bdl, Sampler get_sample,
                              double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  using Visitor = std::decay_t<decltype(bdl.vis_)>;
  using Space = std::decay_t<decltype(*bdl.super_space_)>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);

  generate_lazy_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_sbarrtstar_bidir_no_init(const SBAStarBundle& bdl,
                                            Sampler get_sample,
                                            double SA_init_temperature = 0.0) {
  sbastar_detail::generate_sbarrtstar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.successor_, get_sample, bdl.select_neighborhood_,
      SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                    Sampler get_sample,
                                    double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);

  generate_lazy_sbarrtstar_bidir_no_init(bdl, get_sample, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_bnb_sbarrtstar_no_init(const SBAStarBundle& bdl,
                                          Sampler get_sample,
                                          double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_sbarrtstar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, get_sample, bdl.select_neighborhood_,
        SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_sbarrtstar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, get_sample, bdl.select_neighborhood_,
        SA_init_temperature);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_bnb_sbarrtstar(const SBAStarBundle& bdl, Sampler get_sample,
                                  double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);

  generate_lazy_bnb_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_bnb_sbarrtstar_bidir_no_init(
    const SBAStarBundle& bdl, Sampler get_sample,
    double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_sbarrtstar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.successor_, get_sample, bdl.select_neighborhood_,
        SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_sbarrtstar_bidir_no_init_impl(
        *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_,
        bdl.vis_, bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_,
        bdl.weight_, bdl.density_, bdl.constriction_, bdl.distance_,
        bdl.predecessor_, bdl.fwd_distance_, bdl.successor_, get_sample,
        bdl.select_neighborhood_, SA_init_temperature);
  }
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology (see pp::RandomSampler).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches
  *        to SBA*.
  */
template <typename SBAStarBundle, typename Sampler>
void generate_lazy_bnb_sbarrtstar_bidir(const SBAStarBundle& bdl,
                                        Sampler get_sample,
                                        double SA_init_temperature = 0.0) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);

  generate_lazy_bnb_sbarrtstar_bidir_no_init(bdl, get_sample,
                                             SA_init_temperature);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_SBASTAR_RRTSTAR_H_
