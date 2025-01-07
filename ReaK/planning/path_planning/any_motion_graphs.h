/**
 * \file any_motion_graphs.h
 *
 * This library provides a type-erased classes for motion graphs.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_PLANNING_PATH_PLANNING_ANY_MOTION_GRAPHS_H_
#define REAK_PLANNING_PATH_PLANNING_ANY_MOTION_GRAPHS_H_

#include "bagl/dynamic_graph.h"

#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"

#include <iomanip>
#include <iostream>
#include <type_traits>

namespace ReaK::pp {

/**
 * This struct contains the data required on a per-vertex basis for any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct mg_vertex_data {
  /// The position associated to the vertex.
  topology_point_type_t<Topology> position;
};

/**
 * This struct contains the data required on a per-edge basis for any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology, typename = void>
struct mg_edge_data {};

template <typename Topology>
struct mg_edge_data<
    Topology, std::void_t<std::enable_if_t<is_steerable_space_v<Topology>>>> {
  using steer_record_type = steerable_space_steer_record_t<Topology>;
  steer_record_type steer_record;

  mg_edge_data() : steer_record() {}
  explicit mg_edge_data(const steer_record_type& aRec) : steer_record(aRec) {}
  explicit mg_edge_data(steer_record_type&& aRec)
      : steer_record(std::move(aRec)) {}
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const mg_vertex_data<Topology>& vp) {
  using ReaK::to_vect;
  auto v_pos = to_vect<double>(vp.position);
  for (double x : v_pos) {
    out << " " << std::setw(10) << x;
  }
}

// Make a type-erased encapsulation of a graph used by any path-planning algorithm.
template <typename Topology, typename Graph>
void add_motion_graph_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_position", get(&VBundled::position, g));
  if constexpr (is_steerable_space_v<Topology>) {
    using EBundled = bagl::edge_bundle_type<Graph>;
    dp.property<bagl::graph_edge_descriptor_t<Graph>>(
        "edge_steer_record", get(&EBundled::steer_record, g));
  }
}

/**
 * This struct contains the data required on a per-vertex basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct optimal_mg_vertex : mg_vertex_data<Topology> {
  /// The travel-distance accumulated in the vertex, i.e., the travel-distance from the root vertex to this vertex.
  double distance_accum;
  /// The predecessor associated to the vertex, e.g., following the predecessor links starting at the goal node yields a
  /// backward trace of the optimal path.
  std::size_t predecessor;
};

/**
 * This struct contains the data required on a per-edge basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct optimal_mg_edge : mg_edge_data<Topology> {
  /// The travel-distance associated to the edge (from source to target).
  double weight;

  explicit optimal_mg_edge(double aWeight = 0.0)
      : mg_edge_data<Topology>(), weight(aWeight) {}

  template <typename SteerRec>
  optimal_mg_edge(double aWeight, const SteerRec& aRec)
      : mg_edge_data<Topology>(aRec), weight(aWeight) {}

  template <typename SteerRec>
  optimal_mg_edge(double aWeight, SteerRec&& aRec)
      : mg_edge_data<Topology>(std::move(aRec)), weight(aWeight) {}
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const optimal_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  auto v_pos = to_vect<double>(vp.position);
  for (double x : v_pos) {
    out << " " << std::setw(10) << x;
  }
  out << " " << std::setw(10) << vp.distance_accum;
}

// Make a type-erased encapsulation of a graph used by any optimal path-planning algorithm.
template <typename Topology, typename Graph>
void add_optimal_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_distance_accum", get(&VBundled::distance_accum, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_predecessor", get(&VBundled::predecessor, g));
  using EBundled = bagl::edge_bundle_type<Graph>;
  dp.property<bagl::graph_edge_descriptor_t<Graph>>("edge_weight",
                                                    get(&EBundled::weight, g));
}

/**
 * This struct contains the data required on a per-vertex basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct bidir_optimal_mg_vertex : optimal_mg_vertex<Topology> {
  /// The foward travel-distance accumulated in the vertex, i.e., the travel-distance to the goal vertex from this
  /// vertex.
  double fwd_distance_accum;
  /// The successor associated to the vertex, e.g., following the successor links starting at the start node yields a
  /// trace of the optimal path.
  std::size_t successor;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out,
                     const bidir_optimal_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  auto v_pos = to_vect<double>(vp.position);
  for (double x : v_pos) {
    out << " " << std::setw(10) << x;
  }
  out << " " << std::setw(10) << vp.distance_accum;
  out << " " << std::setw(10) << vp.fwd_distance_accum;
}

// Make a type-erased encapsulation of a graph used by any bidirectional optimal path-planning algorithm.
template <typename Topology, typename Graph>
void add_bidir_optimal_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_fwd_distance_accum", get(&VBundled::fwd_distance_accum, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_successor", get(&VBundled::successor, g));
}

/**
 * This struct contains the data required on a per-vertex basis for any A*-like path-planning algorithm (heuristically
 * driven).
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct astar_mg_vertex : optimal_mg_vertex<Topology> {
  /// The heuristic-value associated to the vertex, i.e., the bird-flight distance to the goal.
  double heuristic_value;
  /// The key-value associated to the vertex, computed by the algorithm (usually a combination of accumulated and
  /// heuristic distances).
  double key_value;
  /// The color-value associated to the vertex, computed by the algorithm.
  bagl::default_color_type astar_color;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const astar_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  auto v_pos = to_vect<double>(vp.position);
  for (double x : v_pos) {
    out << " " << std::setw(10) << x;
  }
  out << " " << std::setw(10) << vp.distance_accum << " " << std::setw(10)
      << vp.heuristic_value << " " << std::setw(10) << vp.key_value;
}

// Make a type-erased encapsulation of a graph used by any A*-like path-planning algorithm.
template <typename Topology, typename Graph>
void add_astar_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_heuristic_value", get(&VBundled::heuristic_value, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_key_value", get(&VBundled::key_value, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_astar_color", get(&VBundled::astar_color, g));
}

/**
 * This struct contains the data required on a per-vertex basis for any A*-like path-planning algorithm (heuristically
 * driven).
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct bidir_astar_mg_vertex : bidir_optimal_mg_vertex<Topology> {
  /// The key-value associated to the vertex, computed by the algorithm (usually a combination of accumulated and
  /// heuristic distances).
  double key_value;
  /// The color-value associated to the vertex, computed by the algorithm.
  bagl::default_color_type astar_color;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out,
                     const bidir_astar_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  auto v_pos = to_vect<double>(vp.position);
  for (double x : v_pos) {
    out << " " << std::setw(10) << x;
  }
  out << " " << std::setw(10) << vp.distance_accum << " " << std::setw(10)
      << vp.fwd_distance_accum << " " << std::setw(10) << vp.key_value;
}

// Make a type-erased encapsulation of a graph used by any bidirectional A*-like path-planning algorithm.
template <typename Topology, typename Graph>
void add_bidir_astar_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_key_value", get(&VBundled::key_value, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_astar_color", get(&VBundled::astar_color, g));
}

/**
 * This struct contains the data required on a per-vertex basis for calculating the density
 * from a survey (non-recursive) of the neighborhood.
 * \tparam BaseVertex The type of the underlying vertex.
 */
template <typename BaseVertex>
struct dense_mg_vertex : BaseVertex {
  /// The density-value associated to the vertex, calculated by a survey of the neighborhood.
  double density;
};

template <typename BaseVertex>
void print_mg_vertex(std::ostream& out, const dense_mg_vertex<BaseVertex>& vp) {
  print_mg_vertex(out, static_cast<const BaseVertex&>(vp));
  out << " " << std::setw(10) << vp.density;
}

// Add a density metric to the type-erased motion graph properties.
template <typename Topology, typename Graph>
void add_density_property_maps(Graph& g, bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_density", get(&VBundled::density, g));
}

/**
 * This struct contains the data required on a per-vertex basis for calculating the density
 * from a recursive accumulation of the neighborhood statistics (e.g., recursive KL-divergence).
 * \tparam BaseVertex The type of the underlying vertex.
 */
template <typename BaseVertex>
struct recursive_dense_mg_vertex : BaseVertex {
  /// The density associated to the vertex, which represents the probability that a sample drawn from the neighborhood
  /// of the vertex will not yield any information gain.
  double density;
  /// Keeps track of the number of neighbors of the vertex.
  std::size_t expansion_trials;
  /// The constriction associated to the vertex, which represents the probability that a sample drawn from the
  /// neighborhood of the vertex will be colliding (or unreachable by a collision-free path).
  double constriction;
  /// Keeps track of the number of neighbors of the vertex that could not be connected to it due to a collision.
  std::size_t collision_count;
};

template <typename BaseVertex>
void print_mg_vertex(std::ostream& out,
                     const recursive_dense_mg_vertex<BaseVertex>& vp) {
  print_mg_vertex(out, static_cast<const BaseVertex&>(vp));
  out << " " << std::setw(10) << vp.density << " " << std::setw(10)
      << vp.expansion_trials << " " << std::setw(10) << vp.constriction << " "
      << std::setw(10) << vp.collision_count;
}

// Add a recursive density metric to the type-erased motion graph properties.
template <typename Topology, typename Graph>
void add_recursive_density_property_maps(Graph& g,
                                         bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_density", get(&VBundled::density, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_constriction", get(&VBundled::constriction, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_expansion_trials", get(&VBundled::expansion_trials, g));
  dp.property<bagl::graph_vertex_descriptor_t<Graph>>(
      "vertex_collision_count", get(&VBundled::collision_count, g));
}

namespace any_motion_graphs_detail {

template <typename VBundled>
struct te_mg_dense_maker {
  template <typename T, typename G>
  static void add_prop_maps(G& /*g*/, bagl::dynamic_properties& /*dp*/) {}
};
template <typename BaseVBundled>
struct te_mg_dense_maker<dense_mg_vertex<BaseVBundled>> {
  template <typename T, typename G>
  static void add_prop_maps(G& g, bagl::dynamic_properties& dp) {
    add_density_property_maps<T>(g, dp);
  }
};
template <typename BaseVBundled>
struct te_mg_dense_maker<recursive_dense_mg_vertex<BaseVBundled>> {
  template <typename T, typename G>
  static void add_prop_maps(G& g, bagl::dynamic_properties& dp) {
    add_recursive_density_property_maps<T>(g, dp);
  }
};

}  // namespace any_motion_graphs_detail

// Make a type-erased encapsulation of a graph used by any path-planning algorithm.
template <typename Topology, typename Graph>
void add_all_motion_graph_property_maps(Graph& g,
                                        bagl::dynamic_properties& dp) {
  using VBundled = bagl::vertex_bundle_type<Graph>;
  if constexpr (std::is_convertible_v<VBundled*, mg_vertex_data<Topology>>) {
    add_motion_graph_property_maps<Topology>(g, dp);
  }
  if constexpr (std::is_convertible_v<VBundled*, optimal_mg_vertex<Topology>>) {
    add_optimal_property_maps<Topology>(g, dp);
  }
  if constexpr (std::is_convertible_v<VBundled*,
                                      bidir_optimal_mg_vertex<Topology>>) {
    add_bidir_optimal_property_maps<Topology>(g, dp);
  }
  if constexpr (std::is_convertible_v<VBundled*, astar_mg_vertex<Topology>>) {
    add_astar_property_maps<Topology>(g, dp);
  }
  if constexpr (std::is_convertible_v<VBundled*,
                                      bidir_astar_mg_vertex<Topology>>) {
    add_bidir_astar_property_maps<Topology>(g, dp);
  }
  any_motion_graphs_detail::te_mg_dense_maker<VBundled>::template add_prop_maps<
      Topology, Graph>(g, dp);
}

/**
 * This stateless functor type can be used to print out the information about an A*-like motion-graph vertex.
 * This is a printing policy type for the vlist_sbmp_report class.
 */
struct mg_vertex_printer : serializable {

  /**
   * This call operator prints all the information about a given vertex to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the planning algorithm.
   * \param out The output-stream to which to print the information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  template <typename Vertex, typename Graph>
  void operator()(std::ostream& out, Vertex u, const Graph& g) const {
    print_mg_vertex(out, g[u]);
    out << std::endl;
  }

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {}
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {}

  RK_RTTI_MAKE_ABSTRACT_1BASE(mg_vertex_printer, 0xC2460011, 1,
                              "mg_vertex_printer", serializable)
};

static constexpr std::size_t BASE_MOTION_GRAPH_KIND_MASK = 0x0F;
static constexpr std::size_t BASIC_MOTION_GRAPH_KIND = 0x00;
static constexpr std::size_t OPTIMAL_MOTION_GRAPH_KIND = 0x01;
static constexpr std::size_t ASTAR_MOTION_GRAPH_KIND = 0x03;
static constexpr std::size_t BIDIR_MOTION_GRAPH_KIND = 0x04;
static constexpr std::size_t BIDIR_OPTIMAL_MOTION_GRAPH_KIND = 0x05;
static constexpr std::size_t BIDIR_ASTAR_MOTION_GRAPH_KIND = 0x07;

static constexpr std::size_t DENSITY_MOTION_GRAPH_KIND_MASK = 0xF0;
static constexpr std::size_t DENSE_MOTION_GRAPH_KIND = 0x10;
static constexpr std::size_t RECURSIVE_DENSE_MOTION_GRAPH_KIND = 0x20;

/**
 * This stateless functor type can be used to print out the information about an A*-like motion-graph vertex.
 * This is a printing policy type for the vlist_sbmp_report class.
 */
template <typename Topology>
struct any_mg_vertex_printer : serializable {
  using self = any_mg_vertex_printer<Topology>;
  using PointType = topology_point_type_t<Topology>;

  std::size_t graph_kind;

  /**
   * This call operator prints all the information about a given vertex to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the planning algorithm.
   * \param out The output-stream to which to print the information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  void operator()(std::ostream& out,
                  bagl::dynamic_graph_observer::vertex_descriptor u,
                  const bagl::dynamic_graph_observer& g) const {
    using ReaK::to_vect;

    auto v_pos = to_vect<double>(
        bagl::get<const PointType&>("vertex_position", g.get_properties(), u));
    for (double x : v_pos) {
      out << " " << std::setw(10) << x;
    }

    if ((graph_kind & OPTIMAL_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << bagl::get<const double&>("vertex_distance_accum",
                                      g.get_properties(), u);
      if ((graph_kind & BIDIR_MOTION_GRAPH_KIND) != 0U) {
        out << " " << std::setw(10)
            << bagl::get<const double&>("vertex_fwd_distance_accum",
                                        g.get_properties(), u);
      }
    }

    if ((graph_kind & ASTAR_MOTION_GRAPH_KIND) != 0U) {
      if ((graph_kind & BIDIR_MOTION_GRAPH_KIND) != 0U) {
        out << " " << std::setw(10)
            << bagl::get<const double&>("vertex_key_value", g.get_properties(),
                                        u);
      } else {
        out << " " << std::setw(10)
            << bagl::get<const double&>("vertex_heuristic_value",
                                        g.get_properties(), u)
            << " " << std::setw(10)
            << bagl::get<const double&>("vertex_key_value", g.get_properties(),
                                        u);
      }
    }

    if ((graph_kind & DENSE_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << bagl::get<const double&>("vertex_density", g.get_properties(), u);
    } else if ((graph_kind & RECURSIVE_DENSE_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << bagl::get<const double&>("vertex_density", g.get_properties(), u)
          << " " << std::setw(10)
          << bagl::get<const std::size_t&>("vertex_expansion_trials",
                                           g.get_properties(), u)
          << " " << std::setw(10)
          << bagl::get<const double&>("vertex_constriction", g.get_properties(),
                                      u)
          << " " << std::setw(10)
          << bagl::get<const std::size_t&>("vertex_collision_count",
                                           g.get_properties(), u);
    }

    out << std::endl;
  }

  explicit any_mg_vertex_printer(std::size_t aGraphKind)
      : graph_kind(aGraphKind) {}

  any_mg_vertex_printer() : any_mg_vertex_printer(0) {}

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(graph_kind);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(graph_kind);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2460012, 1, "any_mg_vertex_printer",
                              serializable)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_ANY_MOTION_GRAPHS_H_
