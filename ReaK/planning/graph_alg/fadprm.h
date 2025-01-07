/**
 * \file fadprm.h
 *
 * This library provides function templates and concepts that implement the flexible anytime dynamic
 * probabilistic roadmap (FADPRM) algorithm (as of "Belghith et al., 2006"). A FADPRM is uses the
 * AD* search algorithm to drive the expansion of a probabilistic roadmap into the free-space
 * in order to connect a start and goal location. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, how to measure density, how to compare
 * densities, when to stop the algorithm, etc. All these customization points are left to the user
 * to implement, some are defined by the FADPRMVisitorConcept (expand-vertex and compute-density)
 * while others are provided as functors to the function template that generates the FADPRM (generate_fadprm).
 *
 * The FADPRM algorithm is so closely related to the AD* algorithm that it is actually implemented
 * using the AD* loop (a function in the detail namespace of the AD* functions) and a customized
 * AD* BFS visitor that wraps the FADPRM customizations along with the FADPRM visitor for the
 * user's customizations.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_PLANNING_GRAPH_ALG_FADPRM_H_
#define REAK_PLANNING_GRAPH_ALG_FADPRM_H_

#include <concepts>
#include <functional>
#include <limits>
#include <type_traits>

#include "ReaK/planning/graph_alg/adstar_search.h"
#include "ReaK/planning/graph_alg/probabilistic_roadmap.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the FADPRM algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the FADPRM algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the FADPRM algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model PRMVisitor.
  *
  * The visitor class should model ADStarVisitor.
  *
  * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  */
template <typename Visitor, typename Graph, typename Space>
concept FADPRMVisitor =
    std::copy_constructible<Visitor>&& PRMVisitor<Visitor, Graph, Space>;

/**
  * This class is simply a "null" visitor for the FADPRM algorithm. It is null in the sense that it
  * will do nothing on all accounts.
  * \tparam Topology The topology type that represents the free-space.
  * \tparam PositionMap The property-map type which can store the position associated to each vertex.
  */
template <pp::MetricSpace Space, typename PositionMap>
class default_fadprm_visitor {
 public:
  using PointType = typename ReaK::pp::topology_traits<Space>::point_type;

  default_fadprm_visitor(const Space& free_space, PositionMap position)
      : free_space_(free_space), position_(position) {}
  default_fadprm_visitor(const default_fadprm_visitor<Space, PositionMap>&) =
      default;

  // AD* visitor functions:
  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Edge, typename Graph>
  void examine_edge(Edge /*unused*/, const Graph& /*unused*/) const {}
  template <typename Graph>
  void publish_path(const Graph& /*unused*/) const {}
  template <typename EdgeIter, typename Graph>
  std::pair<double, EdgeIter> detect_edge_change(
      EdgeIter ei, const Graph& /*unused*/) const {
    return {0.0, ei};
  }
  template <typename Graph>
  double adjust_relaxation(double old_eps, const Graph& /*unused*/) const {
    return old_eps;
  }

  // PRM visitor functions:
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Edge, typename Graph>
  void edge_added(Edge /*unused*/, const Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  std::pair<PointType, bool> random_walk(Vertex /*unused*/,
                                         const Graph& /*unused*/) const {
    return {PointType(), false};
  }
  template <typename Vertex, typename Graph>
  void update_density(Vertex /*unused*/, const Graph& /*unused*/) const {}

  // Common to both PRM and AD* visitors:
  bool keep_going() const { return true; }

  const Space& free_space_;
  PositionMap position_;
};

namespace fadprm_detail {

template <pp::MetricSpace Space, typename AStarHeuristicMap,
          typename UniformCostVisitor, typename UpdatableQueue, typename List,
          typename IndexInHeapMap, typename PredecessorMap, typename KeyMap,
          typename DistanceMap, typename RHSMap, typename WeightMap,
          typename DensityMap, typename PositionMap, typename NcSelector,
          typename ColorMap>
struct fadprm_bfs_visitor {

  using KeyValue = bagl::property_traits_value_t<KeyMap>;
  using ColorValue = bagl::property_traits_value_t<ColorMap>;
  using Color = bagl::color_traits<ColorValue>;
  using DistanceValue = bagl::property_traits_value_t<DistanceMap>;
  using PositionValue = bagl::property_traits_value_t<PositionMap>;

  fadprm_bfs_visitor(const Space& super_space, AStarHeuristicMap h,
                     UniformCostVisitor vis, UpdatableQueue& Q, List& I,
                     IndexInHeapMap index_in_heap, PredecessorMap p, KeyMap k,
                     DistanceMap d, RHSMap rhs, WeightMap w, DensityMap dens,
                     PositionMap pos, NcSelector select_neighborhood,
                     ColorMap col, double& beta)
      : super_space_(super_space),
        h_(h),
        vis_(vis),
        q_(Q),
        incons_(I),
        index_in_heap_(index_in_heap),
        predecessor_(p),
        key_(k),
        distance_(d),
        rhs_(rhs),
        weight_(w),
        density_(dens),
        position_(pos),
        select_neighborhood_(select_neighborhood),
        color_(col),
        beta_(beta) {}

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, Graph& g) const {
    vis_.initialize_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph& g) const {
    vis_.discover_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void inconsistent_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}

  template <class Graph>
  bagl::graph_vertex_descriptor_t<Graph> create_vertex(const PositionValue& p,
                                                       Graph& g) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    using VertexProp = bagl::vertex_property_type<Graph>;

    VertexProp up;
    put(position_, up, p);
    Vertex u = add_vertex(g, std::move(up));
    vis_.vertex_added(u, g);
    put(color_, u, Color::white());
    put(index_in_heap_, u, std::numeric_limits<std::size_t>::max());
    put(distance_, u, std::numeric_limits<double>::infinity());
    put(rhs_, u, std::numeric_limits<double>::infinity());
    put(key_, u,
        KeyValue(std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity()));
    put(predecessor_, u, u);

    return u;
  }

  template <typename Edge, typename Graph>
  void edge_added(Edge e, Graph& g) const {
    vis_.edge_added(e, g);
  }

  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    vis_.vertex_added(u, g);
  }

  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    vis_.examine_vertex(u, g);

    prm_node_connector connect_vertex;

    std::size_t max_node_degree = 10;
    if (out_degree(u, g) >= max_node_degree) {
      return;
    }
    max_node_degree = (max_node_degree - out_degree(u, g)) / 2;
    for (std::size_t i = 0; i < max_node_degree; ++i) {
      auto [p_rnd, expanding_worked, ep] = vis_.random_walk(u, g);
      if (expanding_worked) {
        connect_vertex(p_rnd, u, ep, g, super_space_, *this, position_,
                       select_neighborhood_);
      } else {
        break;
      }
    }
  }

  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) const {
    if (get(weight_, e) < 0.0) {
      throw bagl::negative_edge();
    }
    vis_.examine_edge(e, g);
  }
  template <class Vertex, class Graph>
  void forget_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}
  template <class Vertex, class Graph>
  void finish_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}
  template <class Vertex, class Graph>
  void recycle_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}

  template <typename Graph>
  void publish_path(const Graph& g) const {
    vis_.publish_path(g);
  }

  bool keep_going() const { return vis_.keep_going(); }

  template <typename EdgeIter, typename Graph>
  std::pair<double, EdgeIter> detect_edge_change(EdgeIter ei,
                                                 const Graph& g) const {
    return vis_.detect_edge_change(ei, g);
  }

  template <typename Graph>
  double adjust_epsilon(double old_eps, double /*unused*/,
                        const Graph& g) const {
    return vis_.adjust_relaxation(old_eps, g);
  }

  template <class Vertex, typename Graph>
  void update_key(Vertex u, Graph& g) const {
    vis_.affected_vertex(u, g);
    DistanceValue g_u = get(distance_, u);
    DistanceValue rhs_u = get(rhs_, u);
    DistanceValue h_u = get(h_, u);
    if (rhs_u < g_u) {
      DistanceValue f_u = rhs_u + h_u;  // 0.5 * ( rhs_u + h_u );
      put(key_, u,
          KeyValue((1.0 - beta_) * get(density_, u) + beta_ * f_u, rhs_u));
    } else {
      DistanceValue f_u = g_u + h_u;  // 0.5 * ( g_u + h_u );
      put(key_, u, KeyValue(f_u, rhs_u));
    }
  }

  template <class Vertex, bagl::concepts::BidirectionalGraph Graph>
  void update_vertex(Vertex u, Graph& g) {
    ColorValue col_u = get(color_, u);

    if (col_u == Color::white()) {
      put(distance_, u, std::numeric_limits<double>::infinity());
      col_u = Color::green();
      put(color_, u, col_u);
    }

    DistanceValue g_u = get(distance_, u);
    DistanceValue rhs_u = get(rhs_, u);

    if (rhs_u != 0.0) {
      rhs_u = std::numeric_limits<double>::infinity();
      for (auto e : in_edges(u, g)) {
        DistanceValue rhs_tmp = get(weight_, e) + get(distance_, source(e, g));
        if (rhs_tmp < rhs_u) {
          rhs_u = rhs_tmp;
          put(rhs_, u, rhs_u);
          put(predecessor_, u, source(e, g));
        }
      }
      rhs_u = get(rhs_, u);
    }

    if (rhs_u != g_u) {
      if ((col_u != Color::black()) && (col_u != Color::red())) {
        update_key(u, g);
        q_.push_or_update(u);
        put(color_, u, Color::gray());
        vis_.discover_vertex(u, g);
      } else if (col_u == Color::black()) {
        incons_.push_back(u);
        put(color_, u, Color::red());
      }
    } else if (q_.contains(u)) {
      put(key_, u, KeyValue(0.0, 0.0));
      q_.update(u);
      q_.pop();  // remove from OPEN set
      update_key(u, g);
      put(color_, u, Color::green());
    }
  }

  template <typename Vertex, typename Graph>
  void travel_explored(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_explored(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_succeeded(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void travel_failed(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_failed(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex u, Graph& g) {
    update_key(u, g);
    update_vertex(u, g);
  }

  template <typename Vertex, typename Graph>
  std::pair<bool, bagl::edge_bundle_type<Graph>> can_be_connected(
      Vertex u, Vertex v, Graph& g) const {
    return vis_.can_be_connected(u, v, g);
  }

  template <typename Vertex, typename Graph>
  std::tuple<PositionValue, bool, bagl::edge_bundle_type<Graph>> random_walk(
      Vertex u, Graph& g) const {
    return vis_.random_walk(u, g);
  }

  bool is_position_free(const PositionValue& p) const {
    return vis_.is_position_free(p);
  }

  const Space& super_space_;
  AStarHeuristicMap h_;
  UniformCostVisitor vis_;
  UpdatableQueue& q_;
  List& incons_;
  IndexInHeapMap index_in_heap_;
  PredecessorMap predecessor_;
  KeyMap key_;
  DistanceMap distance_;
  RHSMap rhs_;
  WeightMap weight_;
  DensityMap density_;
  PositionMap position_;
  NcSelector select_neighborhood_;
  ColorMap color_;
  DistanceValue& beta_;
};

}  // namespace fadprm_detail

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the FADPRM algorithm (as of Belghith et al., 2006), without initialization of the
  * existing graph.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         BidirectionalGraphConcept and MutableGraphConcept.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
  *         for each vertex in the graph.
  * \tparam Visitor The type of the FADPRM visitor to be used.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex together with its optimal predecessor.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex
  *         to the goal.
  * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of
  *         each vertex to the goal (internal use to AD*).
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the
  *         graph (cost of travel along an edge).
  * \tparam PositionMap A property-map type that can store the position of each vertex.
  * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
  * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
  *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
  *         See classes in the topological_search.hpp header-file.
  * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors
  *         are used to mark vertices by their status in the AD* algorithm (white = not visited,
  *         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED,
  *         not OPEN), red = inconsistent (in INCONS)).
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param free_space A topology (as defined by the Born Again Graph Library). Note
  *        that it is required to generate only random points in
  *        the free-space and to only allow interpolation within the free-space.
  * \param hval The property-map of A* heuristic function values for each vertex.
  * \param vis A FADPRM visitor implementing the FADPRMVisitor. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
  *        goal (for internal use).
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param density A property-map that provides the density values assiciated to each vertex.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param select_neighborhood A callable object (functor) that can select a list of
  *        vertices of the graph that ought to be connected to a new
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param color The property-map which stores the color-value of the vertices, colors are used to mark
  *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
  *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  */
template <bagl::concepts::VertexListGraph G,  // this is the actual graph.
          pp::Topology Space,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          FADPRMVisitor<G, Space> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadableVertexPropertyMap<G> DensityMap,
          bagl::concepts::ReadWriteVPropMemberMap<G> PositionMap,
          typename NcSelector,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void generate_fadprm_no_init(G& g,
                             bagl::graph_vertex_descriptor_t<G> start_vertex,
                             const Space& free_space, AStarHeuristicMap hval,
                             Visitor vis, PredecessorMap predecessor,
                             DistanceMap distance, RHSMap rhs, WeightMap weight,
                             DensityMap density, PositionMap position,
                             NcSelector select_neighborhood, ColorMap color,
                             double epsilon) {
  using KeyValue = adstar_key_value<double>;
  using KeyCompareType = typename adstar_key_traits<KeyValue>::compare_type;
  auto key_map = bagl::vector_property_map(
      num_vertices(g), get(bagl::vertex_index, g),
      KeyValue(std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity()));
  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  auto q =
      bagl::make_d_ary_heap_indirect<bagl::graph_vertex_descriptor_t<G>, 4>(
          key_map.ref(), index_in_heap.ref(), KeyCompareType());
  // list holding the INCONS set (inconsistent nodes).
  std::vector<bagl::graph_vertex_descriptor_t<G>> incons;

  auto bfs_vis = fadprm_detail::fadprm_bfs_visitor(
      free_space, hval, vis, q, incons, index_in_heap.ref(), predecessor,
      key_map.ref(), distance, rhs, weight, density, position,
      select_neighborhood, color, epsilon);

  adstar_detail::adstar_search_loop(
      g, start_vertex, hval, bfs_vis, predecessor, distance, rhs, key_map.ref(),
      weight, color, index_in_heap.ref(), q, incons, epsilon,
      std::numeric_limits<double>::infinity(), 0.0, std::less<>(),
      std::equal_to<>(), std::plus<>(), std::multiplies<>());
}

/**
* This function template generates a roadmap to connect a goal location to a start location
* using the FADPRM algorithm (as of Belghith et al., 2006), with initialization of the
* existing graph to (re)start the search.
* \tparam Graph The graph type that can store the generated roadmap, should model
*         BidirectionalGraphConcept and MutableGraphConcept.
* \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
* \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
* \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
*         for each vertex in the graph.
* \tparam Visitor The type of the FADPRM visitor to be used.
* \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
*         vertex together with its optimal predecessor.
* \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex
*         to the goal.
* \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of
*         each vertex to the goal (internal use to AD*).
* \tparam WeightMap This property-map type is used to store the weights of the edges of the
*         graph (cost of travel along an edge).
* \tparam PositionMap A property-map type that can store the position of each vertex.
* \tparam DensityMap A property-map type that can store the density-measures for each vertex.
* \tparam NcSelector A functor type that can select a list of vertices of the graph that are
*         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
*         See classes in the topological_search.hpp header-file.
* \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors
*         are used to mark vertices by their status in the AD* algorithm (white = not visited,
*         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED,
*         not OPEN), red = inconsistent (in INCONS)).
*
* \param g A mutable graph that should initially store the starting
*        vertex (if not it will be randomly generated) and will store
*        the generated graph once the algorithm has finished.
* \param start_vertex The starting point of the algorithm, on the graph.
* \param free_space A topology (as defined by the Born Again Graph Library). Note
*        that it is required to generate only random points in
*        the free-space and to only allow interpolation within the free-space.
* \param hval The property-map of A* heuristic function values for each vertex.
* \param vis A FADPRM visitor implementing the FADPRMVisitor. This is the
*        main point of customization and recording of results that the
*        user can implement.
* \param predecessor The property-map which will store the resulting path by connecting
*        vertices together with their optimal predecessor (follow in reverse to discover the
*        complete path).
* \param distance The property-map which stores the estimated distance of each vertex to the goal.
* \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
*        goal (for internal use).
* \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
*        along the edge).
* \param density A property-map that provides the density values assiciated to each vertex.
* \param position A mapping that implements the MutablePropertyMap Concept. Also,
*        the value_type of this map should be the same type as the topology's
*        value_type.
* \param select_neighborhood A callable object (functor) that can select a list of
*        vertices of the graph that ought to be connected to a new
*        vertex. The list should be sorted in order of increasing "distance".
* \param color The property-map which stores the color-value of the vertices, colors are used to mark
*        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
*        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
* \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
*        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
*/
template <bagl::concepts::VertexListGraph G,  // this is the actual graph.
          pp::MetricSpace Space,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          FADPRMVisitor<G, Space> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadableVertexPropertyMap<G> DensityMap,
          bagl::concepts::ReadWriteVPropMemberMap<G> PositionMap,
          typename NcSelector,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap>
requires pp::PointDistribution<Space> inline void generate_fadprm(
    G& g, bagl::graph_vertex_descriptor_t<G> start_vertex,
    const Space& free_space, AStarHeuristicMap hval, Visitor vis,
    PredecessorMap predecessor, DistanceMap distance, RHSMap rhs,
    WeightMap weight, DensityMap density, PositionMap position,
    NcSelector select_neighborhood, ColorMap color, double epsilon) {
  using ColorValue = bagl::property_traits_value_t<ColorMap>;
  using Color = bagl::color_traits<ColorValue>;
  using VertexProp = bagl::vertex_property_type<G>;

  if (num_vertices(g) == 0) {
    VertexProp up;
    auto p = get(ReaK::pp::random_sampler, free_space)(free_space);
    put(position, up, p);
    start_vertex = add_vertex(g, std::move(up));
    vis.vertex_added(start_vertex, g);
  }

  for (auto u : vertices(g)) {
    put(color, u, Color::white());
    put(distance, u, std::numeric_limits<double>::infinity());
    put(rhs, u, std::numeric_limits<double>::infinity());
    put(predecessor, u, u);
    vis.initialize_vertex(u, g);
  }

  generate_fadprm_no_init(g, start_vertex, free_space, hval, vis, predecessor,
                          distance, rhs, weight, density, position,
                          select_neighborhood, color, epsilon);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_FADPRM_H_
