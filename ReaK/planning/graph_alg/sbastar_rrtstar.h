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
#include <utility>

#include "ReaK/core/base/global_rng.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include <tuple>
#include "boost/graph/detail/d_ary_heap.hpp"
#include "boost/graph/graph_concepts.hpp"
#include "boost/graph/properties.hpp"

// BGL-Extra includes:
#include "boost/graph/more_property_maps.hpp"
#include "boost/graph/more_property_tags.hpp"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/node_generators.h"
#include "ReaK/planning/graph_alg/sbastar_search.h"
#include "ReaK/planning/graph_alg/simple_graph_traits.h"

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

namespace detail {
namespace {

template <typename Graph, pp::MetricSpace Space, typename Visitor,
          typename MotionGraphConnector, typename SBANodeGenerator,
          typename RRTNodeGenerator, typename MutableQueue, typename NcSelector>
void sbarrtstar_search_loop(Graph& g, const Space& super_space,
                            Visitor& sba_vis,
                            MotionGraphConnector connect_vertex,
                            SBANodeGenerator sba_generate_node,
                            RRTNodeGenerator rrt_generate_node, MutableQueue& Q,
                            NcSelector select_neighborhood,
                            double initial_temperature) {
  using std::exp;
  using std::log;
  std::size_t num_rrt_vertices = 0;
  std::size_t num_sba_vertices = 0;
  std::uniform_real_distribution<double> unit_dist{};

  const auto connect_if_needed = [&](auto& x_near, auto& p_new, auto& eprop) {
    if ((x_near != boost::graph_traits<Graph>::null_vertex()) &&
        (get(sba_vis.m_predecessor, g[x_near]) !=
         boost::graph_traits<Graph>::null_vertex())) {
      connect_vertex(p_new, x_near, eprop, g, super_space, sba_vis,
                     sba_vis.m_position, sba_vis.m_distance,
                     sba_vis.m_predecessor, sba_vis.m_weight,
                     select_neighborhood);
    }
  };

  while (!Q.empty() && sba_vis.keep_going()) {
    double entropy =
        1.0 - exp(-initial_temperature / log(double(num_vertices(g))));
    // generate random-number between 0 and 1.
    double rand_value = unit_dist(ReaK::get_global_rng());
    bool use_sba_sampling = (rand_value > entropy);

    if (use_sba_sampling) {
      auto u = Q.top();
      Q.pop();

      // stop if the best nodes do not meet the potential threshold.
      while (!sba_vis.has_search_potential(u, g)) {
        if (Q.empty()) {
          u = boost::graph_traits<Graph>::null_vertex();
          break;
        }
        u = Q.top();
        Q.pop();
      }
      if (u == boost::graph_traits<Graph>::null_vertex()) {
        break;  // no more nodes with search potential.
      }

      sba_vis.examine_vertex(u, g);

      auto [x_near, p_new, eprop] =
          sba_generate_node(u, g, sba_vis, sba_vis.m_position);
      ++num_sba_vertices;

      // then push it back on the OPEN queue.
      if ((x_near != boost::graph_traits<Graph>::null_vertex()) || Q.empty()) {
        sba_vis.requeue_vertex(u, g);
      }

      connect_if_needed(x_near, p_new, eprop);
    } else {
      auto [x_near, p_new, eprop] = rrt_generate_node(
          g, sba_vis, boost::bundle_prop_to_vertex_prop(sba_vis.m_position, g));
      ++num_rrt_vertices;

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
                           RRTNodeGenerator rrt_generate_node, MutableQueue& Q,
                           NcSelector select_neighborhood,
                           double initial_temperature) {
  using std::exp;
  using std::log;
  std::size_t num_rrt_vertices = 0;
  std::size_t num_sba_vertices = 0;
  std::uniform_real_distribution<double> unit_dist{};

  const auto connect_if_needed = [&](auto& x_near_pred, auto& p_new_pred,
                                     auto& ep_pred, auto& x_near_succ,
                                     auto& p_new_succ, auto& ep_succ) {
    if (x_near_pred != boost::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = boost::graph_traits<Graph>::null_vertex();
      graph_edge_bundle_t<Graph> ep_other;
      connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other,
                     g, super_space, sba_vis, sba_vis.m_position,
                     sba_vis.m_distance, sba_vis.m_predecessor,
                     sba_vis.m_fwd_distance, sba_vis.m_successor,
                     sba_vis.m_weight, select_neighborhood);
    }
    if (x_near_succ != boost::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = boost::graph_traits<Graph>::null_vertex();
      graph_edge_bundle_t<Graph> ep_other;
      connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ,
                     g, super_space, sba_vis, sba_vis.m_position,
                     sba_vis.m_distance, sba_vis.m_predecessor,
                     sba_vis.m_fwd_distance, sba_vis.m_successor,
                     sba_vis.m_weight, select_neighborhood);
    }
  };

  while (!Q.empty() && sba_vis.keep_going()) {
    double entropy =
        1.0 - exp(-initial_temperature / log(double(num_vertices(g))));
    // generate random-number between 0 and 1.
    double rand_value = unit_dist(ReaK::get_global_rng());
    bool use_sba_sampling = (rand_value > entropy);

    if (use_sba_sampling) {
      auto u = Q.top();
      Q.pop();

      // stop if the best nodes do not meet the potential threshold.
      while (!sba_vis.has_search_potential(u, g)) {
        if (Q.empty()) {
          u = boost::graph_traits<Graph>::null_vertex();
          break;
        }
        u = Q.top();
        Q.pop();
      }
      if (u == boost::graph_traits<Graph>::null_vertex()) {
        break;  // no more nodes with search potential.
      }

      sba_vis.examine_vertex(u, g);

      auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ,
            ep_succ] = sba_generate_node(u, g, sba_vis, sba_vis.m_position);
      ++num_sba_vertices;

      // then push it back on the OPEN queue.
      if ((x_near_pred != boost::graph_traits<Graph>::null_vertex()) ||
          (x_near_succ != boost::graph_traits<Graph>::null_vertex()) ||
          (Q.empty())) {
        sba_vis.requeue_vertex(u, g);
      }

      connect_if_needed(x_near_pred, p_new_pred, ep_pred, x_near_succ,
                        p_new_succ, ep_succ);
    } else {
      auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ,
            ep_succ] = rrt_generate_node(g, sba_vis,
                                         boost::bundle_prop_to_vertex_prop(
                                             sba_vis.m_position, g));

      ++num_rrt_vertices;

      connect_if_needed(x_near_pred, p_new_pred, ep_pred, x_near_succ,
                        p_new_succ, ep_succ);
    }
  }  // while
}

template <typename Graph, pp::MetricSpace Space,
          SBARRTStarVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename RandomSampler, typename NcSelector>
void generate_sbarrtstar_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex, const Space& super_space,
    Visitor vis, NodeConnector connect_vertex, KeyMap key, PositionMap position,
    WeightMap weight, DensityMap density, ConstrictionMap constriction,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, RandomSampler get_sample,
    NcSelector select_neighborhood, double SA_init_temperature = 0.0) {
  using Vertex = graph_vertex_t<Graph>;

  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap, KeyCompareType());

  sbastar_bfs_visitor<Graph, Visitor, MutableQueue, IndexInHeapMap, KeyMap,
                      PositionMap, WeightMap, DensityMap, ConstrictionMap,
                      DistanceMap, PredecessorMap, FwdDistanceMap>
      sba_bfs_vis(vis, Q, index_in_heap, key, position, weight, density,
                  constriction, distance, predecessor, fwd_distance);

  put(distance, g[start_vertex], 0.0);
  put(predecessor, g[start_vertex], start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbarrtstar_search_loop(g, super_space, sba_bfs_vis, connect_vertex,
                         sba_node_generator(),
                         rrg_node_generator<Space, RandomSampler, NcSelector>(
                             &super_space, get_sample, select_neighborhood),
                         Q, select_neighborhood, SA_init_temperature);
}

template <typename Graph, pp::MetricSpace Space,
          SBARRTStarBidirVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap, typename RandomSampler, typename NcSelector>
void generate_sbarrtstar_bidir_no_init_impl(
    Graph& g, graph_vertex_t<Graph> start_vertex,
    graph_vertex_t<Graph> goal_vertex, const Space& super_space, Visitor vis,
    NodeConnector connect_vertex, KeyMap key, PositionMap position,
    WeightMap weight, DensityMap density, ConstrictionMap constriction,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, SuccessorMap successor,
    RandomSampler get_sample, NcSelector select_neighborhood,
    double SA_init_temperature = 0.0) {
  using Vertex = graph_vertex_t<Graph>;

  using IndexInHeapMap = boost::vector_property_map<std::size_t>;
  IndexInHeapMap index_in_heap;
  for (auto [ui, ui_end] = vertices(g); ui != ui_end; ++ui) {
    put(index_in_heap, *ui, static_cast<std::size_t>(-1));
  }

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  using MutableQueue = boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap,
                                                  KeyMap, KeyCompareType>;
  MutableQueue Q(key, index_in_heap, KeyCompareType());

  sbastar_bfs_visitor<Graph, Visitor, MutableQueue, IndexInHeapMap, KeyMap,
                      PositionMap, WeightMap, DensityMap, ConstrictionMap,
                      DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>
      sba_bfs_vis(vis, Q, index_in_heap, key, position, weight, density,
                  constriction, distance, predecessor, fwd_distance, successor);

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
      rrg_bidir_generator(&super_space, get_sample, select_neighborhood,
                          predecessor, successor),
      Q, select_neighborhood, SA_init_temperature);
}

}  // namespace
}  // namespace detail

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
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);

  detail::generate_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
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
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);

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
  detail::generate_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, pruned_node_connector(), bdl.m_key, bdl.m_position,
      bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
      bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
      bdl.m_select_neighborhood, SA_init_temperature);
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
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);

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
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  detail::generate_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
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
void generate_lazy_sbarrtstar(const SBAStarBundle& bdl, Sampler get_sample,
                              double SA_init_temperature = 0.0) {
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  using Visitor = std::decay_t<decltype(bdl.m_vis)>;
  using Space = std::decay_t<decltype(*(bdl.m_super_space))>;
  static_assert(SBARRTStarVisitor<Visitor, Graph, Space>);
  static_assert(pp::RandomSampler<Sampler, Space>);

  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);

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
  detail::generate_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
      bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
      bdl.m_fwd_distance, bdl.m_successor, get_sample,
      bdl.m_select_neighborhood, SA_init_temperature);
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
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);

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
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_sbarrtstar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
        SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_sbarrtstar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis,
        bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position, bdl.m_weight,
        bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor,
        bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood,
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
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor);

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
  using Graph = std::decay_t<decltype(*(bdl.m_g))>;
  if (bdl.m_goal_vertex == boost::graph_traits<Graph>::null_vertex()) {
    detail::generate_sbarrtstar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, lazy_node_connector(), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
        bdl.m_select_neighborhood, SA_init_temperature);
  } else {
    bnb_ordering_data<Graph> bnb_data(*(bdl.m_g), bdl.m_start_vertex,
                                      bdl.m_goal_vertex);
    detail::generate_sbarrtstar_bidir_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space),
        bdl.m_vis, bnb_connector<Graph>(bnb_data), bdl.m_key, bdl.m_position,
        bdl.m_weight, bdl.m_density, bdl.m_constriction, bdl.m_distance,
        bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, get_sample,
        bdl.m_select_neighborhood, SA_init_temperature);
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
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key,
                                   bdl.m_distance, bdl.m_predecessor,
                                   bdl.m_fwd_distance, bdl.m_successor);

  generate_lazy_bnb_sbarrtstar_bidir_no_init(bdl, get_sample,
                                             SA_init_temperature);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_SBASTAR_RRTSTAR_H_
