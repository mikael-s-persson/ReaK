/**
 * \file any_motion_graphs.hpp
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
 *    If not, sMOTION_ee <http://www.gnu.org/licenses/>.
 */

#ifndef RK_ANY_MOTION_GRAPHS_HPP
#define RK_ANY_MOTION_GRAPHS_HPP

#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"
#include "ReaK/planning/graph_alg/any_graph.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/steerable_space_concept.hpp"

#include <iomanip>
#include <iostream>
#include <type_traits>

namespace ReaK::pp {

namespace detail {
namespace {

template <typename Topology, typename Graph>
auto& try_get_steer_record(Graph& g, graph::graph_edge_t<Graph> e) {
  if constexpr (is_steerable_space_v<Topology>) {
    return g[e].steer_record;
  } else {
    throw std::invalid_argument(
        "Required property 'edge_steer_record' on a non-steerable space!");
    return g[e];
  }
}

}  // namespace
}  // namespace detail

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

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_motion_graph : public graph::type_erased_graph<Graph> {
 protected:
  using self = any_motion_graph<Topology, Graph>;
  using base_type = graph::type_erased_graph<Graph>;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {

    if (aProperty == "vertex_position") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .position));
    }
    if (aProperty == "edge_steer_record") {
      return static_cast<void*>(&(detail::try_get_steer_record<Topology>(
          *(this->p_graph), std::any_cast<real_edge_desc>(aElement))));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {

    if (aProperty == "vertex_position") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .position);
    }
    if (aProperty == "edge_steer_record") {
      return std::any(detail::try_get_steer_record<Topology>(
          *(this->p_graph), std::any_cast<real_edge_desc>(aElement)));
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any optimal path-planning
 * algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_optimal_motion_graph : public any_motion_graph<Topology, Graph> {
 protected:
  using self = any_optimal_motion_graph<Topology, Graph>;
  using base_type = any_motion_graph<Topology, Graph>;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {

    if (aProperty == "vertex_distance_accum") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .distance_accum));
    }
    if (aProperty == "vertex_predecessor") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .predecessor));
    }
    if (aProperty == "edge_weight") {
      return static_cast<void*>(&(
          (*(this->p_graph))[std::any_cast<real_edge_desc>(aElement)].weight));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {

    if (aProperty == "vertex_distance_accum") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .distance_accum);
    }
    if (aProperty == "vertex_predecessor") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .predecessor);
    }
    if (aProperty == "edge_weight") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_edge_desc>(aElement)].weight);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_optimal_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any optimal path-planning
 * algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_bidir_optimal_motion_graph
    : public any_optimal_motion_graph<Topology, Graph> {
 protected:
  using self = any_bidir_optimal_motion_graph<Topology, Graph>;
  using base_type = any_optimal_motion_graph<Topology, Graph>;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {

    if (aProperty == "vertex_fwd_distance_accum") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .fwd_distance_accum));
    }
    if (aProperty == "vertex_successor") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .successor));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {

    if (aProperty == "vertex_fwd_distance_accum") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .fwd_distance_accum);
    }
    if (aProperty == "vertex_successor") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .successor);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_bidir_optimal_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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
  boost::default_color_type astar_color;
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

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any A*-like path-planning
 * algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_astar_motion_graph
    : public any_optimal_motion_graph<Topology, Graph> {
 protected:
  using self = any_astar_motion_graph<Topology, Graph>;
  using base_type = any_optimal_motion_graph<Topology, Graph>;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {
    if (aProperty == "vertex_heuristic_value") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .heuristic_value));
    }
    if (aProperty == "vertex_key_value") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .key_value));
    }
    if (aProperty == "vertex_astar_color") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .astar_color));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {
    if (aProperty == "vertex_heuristic_value") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .heuristic_value);
    }
    if (aProperty == "vertex_key_value") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .key_value);
    }
    if (aProperty == "vertex_astar_color") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .astar_color);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_astar_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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
  boost::default_color_type astar_color;
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

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any A*-like path-planning
 * algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_bidir_astar_motion_graph
    : public any_bidir_optimal_motion_graph<Topology, Graph> {
 protected:
  using self = any_bidir_astar_motion_graph<Topology, Graph>;
  using base_type = any_bidir_optimal_motion_graph<Topology, Graph>;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {
    if (aProperty == "vertex_key_value") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .key_value));
    }
    if (aProperty == "vertex_astar_color") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .astar_color));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {
    if (aProperty == "vertex_key_value") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .key_value);
    }
    if (aProperty == "vertex_astar_color") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .astar_color);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_bidir_astar_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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

/**
 * This class template can be used as a type-erased encapsulation of a graph used
 * by a path-planning algorithm that uses a (non-recursive) density metric.
 * \tparam BaseMotionGraph The type-erased motion-graph base-type used.
 */
template <typename BaseMotionGraph>
class any_dense_motion_graph : public BaseMotionGraph {
 protected:
  using self = any_dense_motion_graph<BaseMotionGraph>;
  using base_type = BaseMotionGraph;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {
    if (aProperty == "vertex_density") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .density));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {
    if (aProperty == "vertex_density") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .density);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_dense_motion_graph(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

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

/**
 * This class template can be used as a type-erased encapsulation of a graph used
 * by a path-planning algorithm that uses a recursive density metric.
 * \tparam BaseMotionGraph The type-erased motion-graph base-type used.
 */
template <typename BaseMotionGraph>
class any_recursive_dense_mg : public BaseMotionGraph {
 protected:
  using self = any_recursive_dense_mg<BaseMotionGraph>;
  using base_type = BaseMotionGraph;
  using original_graph_type = typename base_type::original_graph_type;
  using real_vertex_desc = typename base_type::real_vertex_desc;
  using real_edge_desc = typename base_type::real_edge_desc;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& aElement) const override {
    if (aProperty == "vertex_density") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .density));
    }
    if (aProperty == "vertex_constriction") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .constriction));
    }
    if (aProperty == "vertex_expansion_trials") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .expansion_trials));
    }
    if (aProperty == "vertex_collision_count") {
      return static_cast<void*>(
          &((*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
                .collision_count));
    }

    return base_type::get_property_by_ptr(aProperty, aElement);
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {
    if (aProperty == "vertex_density") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .density);
    }
    if (aProperty == "vertex_constriction") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .constriction);
    }
    if (aProperty == "vertex_expansion_trials") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .expansion_trials);
    }
    if (aProperty == "vertex_collision_count") {
      return std::any(
          (*(this->p_graph))[std::any_cast<real_vertex_desc>(aElement)]
              .collision_count);
    }

    return base_type::get_property_by_any(aProperty, aElement);
  }

 public:
  explicit any_recursive_dense_mg(original_graph_type* aPGraph)
      : base_type(aPGraph) {}
};

template <typename Topology, typename Graph>
struct te_mg_selector {
  using VertexProp = graph::graph_vertex_bundle_t<Graph>;

  static constexpr bool IsBasicMG =
      std::is_convertible_v<VertexProp*, mg_vertex_data<Topology>*>;
  static constexpr bool IsOptimMG =
      std::is_convertible_v<VertexProp*, optimal_mg_vertex<Topology>*>;
  static constexpr bool IsAStarMG =
      std::is_convertible_v<VertexProp*, astar_mg_vertex<Topology>*>;

  using BaseVertexProp = std::conditional_t<
      IsBasicMG,
      std::conditional_t<
          IsOptimMG,
          std::conditional_t<IsAStarMG, astar_mg_vertex<Topology>,
                             optimal_mg_vertex<Topology>>,
          mg_vertex_data<Topology>>,
      void>;

  using BaseMG = std::conditional_t<
      IsBasicMG,
      std::conditional_t<
          IsOptimMG,
          std::conditional_t<IsAStarMG, any_astar_motion_graph<Topology, Graph>,
                             any_optimal_motion_graph<Topology, Graph>>,
          any_motion_graph<Topology, Graph>>,
      graph::type_erased_graph<Graph>>;

  static constexpr bool IsDenseMG =
      std::is_convertible_v<VertexProp*, dense_mg_vertex<BaseVertexProp>*>;
  static constexpr bool IsRecDenseMG =
      std::is_convertible_v<VertexProp*,
                            recursive_dense_mg_vertex<BaseVertexProp>*>;

  using type = std::conditional_t<
      IsDenseMG, any_dense_motion_graph<BaseMG>,
      std::conditional_t<IsRecDenseMG, any_recursive_dense_mg<BaseMG>, BaseMG>>;
};

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
  void operator()(std::ostream& out, graph::any_graph::vertex_descriptor u,
                  const graph::any_graph& g) const {
    using ReaK::to_vect;
    using ReaK::graph::get_dyn_prop;

    auto v_pos = to_vect<double>(
        get_dyn_prop<const PointType&>("vertex_position", u, g));
    for (double x : v_pos) {
      out << " " << std::setw(10) << x;
    }

    if ((graph_kind & OPTIMAL_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << get_dyn_prop<const double&>("vertex_distance_accum", u, g);
      if ((graph_kind & BIDIR_MOTION_GRAPH_KIND) != 0U) {
        out << " " << std::setw(10)
            << get_dyn_prop<const double&>("vertex_fwd_distance_accum", u, g);
      }
    }

    if ((graph_kind & ASTAR_MOTION_GRAPH_KIND) != 0U) {
      if ((graph_kind & BIDIR_MOTION_GRAPH_KIND) != 0U) {
        out << " " << std::setw(10)
            << get_dyn_prop<const double&>("vertex_key_value", u, g);
      } else {
        out << " " << std::setw(10)
            << get_dyn_prop<const double&>("vertex_heuristic_value", u, g)
            << " " << std::setw(10)
            << get_dyn_prop<const double&>("vertex_key_value", u, g);
      }
    }

    if ((graph_kind & DENSE_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << get_dyn_prop<const double&>("vertex_density", u, g);
    } else if ((graph_kind & RECURSIVE_DENSE_MOTION_GRAPH_KIND) != 0U) {
      out << " " << std::setw(10)
          << get_dyn_prop<const double&>("vertex_density", u, g) << " "
          << std::setw(10)
          << get_dyn_prop<const std::size_t&>("vertex_expansion_trials", u, g)
          << " " << std::setw(10)
          << get_dyn_prop<const double&>("vertex_constriction", u, g) << " "
          << std::setw(10)
          << get_dyn_prop<const std::size_t&>("vertex_collision_count", u, g);
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

#endif
