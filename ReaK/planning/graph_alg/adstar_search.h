/**
 * \file adstar_search.h
 *
 * This library implements the Anytime Dynamic A* (AD*) algorithm according to the algorithmic description given
 * in the original article:
 *
 * M. Likhachev, D. Ferguson, G. Gordon, A. Stentz and S. Thrun (2005). "Anytime Dynamic A*: An Anytime, Replanning
 * Algorithm". Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS), June, 2005.
 *
 * The algorithm assumes a fixed starting point for the search (usually the fixed goal of a path-planner) and an A*
 * heuristic function that can be computed for any vertex that is being explored. The underlying graph search is an
 * A* algorithm that was modified for sub-optimal, anytime computation. It will start at a given initial
 * epsilon value and rapidly find a epsilon-optimal path. Afterwards, the espilon can be adjusted between each run
 * of the graph search as a function of the amount of change in the environment (usually decreasing while no important
 * changes are reported, and a reset to the initial value if important changes occur). The AD* will keep in memory the
 * previously computed solution to improve upon it at future graph searches. AD* uses two callback functions: one to
 * publish a path once it has been resolved by the graph search; and another to query for the list of weights of the
 * graph that have undergone some (significant) change and a scalar value that is representative of the amount of
 * change. For all other aspects related to observing the changes in the graph during the graph search, the A* visitor
 * concept is used, as documented in the Boost.Graph library's homepage.
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

#ifndef REAK_PLANNING_GRAPH_ALG_ADSTAR_SEARCH_H_
#define REAK_PLANNING_GRAPH_ALG_ADSTAR_SEARCH_H_

#include <functional>
#include <limits>
#include <vector>
#include "bagl/d_ary_heap.h"
#include "bagl/exception.h"
#include "bagl/graph_concepts.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/visitors.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the AD* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the AD* algorithm (see fadprm.hpp for
  * example). In other words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the AD* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Valid expressions:
  *
  * vis.initialize_vertex(u,g);  A function that gets called whenever a vertex (u) is first initialized before the
  *search.
  *
  * vis.finish_vertex(u,g);  A function that gets called whenever a vertex (u) is put into the CLOSED set (after being
  *explored by the current search pass).
  *
  * vis.recycle_vertex(u,g);  A function that gets called whenever a vertex (u) is taken out of the CLOSED set to be
  *recycled in the search.
  *
  * vis.discover_vertex(u,g);  A function that gets called whenever a vertex (u) is added to the OPEN set (or updated in
  *the OPEN set).
  *
  * vis.examine_vertex(u,g);  A function that gets called whenever a vertex (u) is taken out of the OPEN set to be
  *examined, this is called before it gets expanded.
  *
  * vis.examine_edge(e,g);  A function that gets called whenever an edge (e) is being looked at, as it comes out of the
  *vertex that is currently being examined in the search.
  *
  * vis.edge_relaxed(e,g);  A function that gets called whenever an edge (e) has been newly declared as a better
  *alternative than the current surrounding edges (i.e. the edge is added to the optimal path).
  *
  * vis.forget_vertex(u,g);  A function that gets called whenever a vertex (u) is deemed uninteresting for the search
  *and is taken out of the OPEN set and not examined.
  *
  * vis.inconsistent_vertex(u,g);  A function that gets called whenever a CLOSEd vertex becomes inconsistent (added to
  *the INCONS set) due to changes in the environment or sub-optimal search.
  *
  * vis.publish_path(g);  A function to notify the visitor that at least one A* round has completed and its resulting
  *path (partial or complete) can be published (the path is encoded in the predecessor property-map).
  *
  * bool b = vis.keep_going();  A function to check to see whether the task is finished (return false) or needs to keep
  *going (true).
  *
  * pair<double, EdgeIter> w_change = vis.detect_edge_change(ei,g); A function call upon the detection of edge changes,
  *ei: back-inserter / forward-iterator for an edge-list. Return the cummulative weight-change and the edge-iterator at
  *the end of the edge list populated by this function.
  *
  * double new_eps = vis.adjust_epsilon(old_eps, w_change.first, g);  A function to adjust the value of epsilon for a
  *given old-value and last cummulative weight-change.
  */
template <typename Visitor, typename Graph>
concept ADStarVisitor = std::copy_constructible<Visitor>&& requires(
    Visitor vis, Graph g, bagl::graph_vertex_descriptor_t<Graph> u,
    bagl::graph_edge_descriptor_t<Graph> e) {
  // whenever the vertex is first initialized.
  vis.initialize_vertex(u, g);
  // whenever a vertex is added to the CLOSED set.
  vis.finish_vertex(u, g);
  // whenever a vertex is taken out of the CLOSED set.
  vis.recycle_vertex(u, g);
  // whenever a vertex is added to the OPEN set (or updated in OPEN).
  vis.discover_vertex(u, g);
  // whenever a vertex is taken out of OPEN, before it gets "expanded".
  vis.examine_vertex(u, g);
  // whenever an edge is being looked at (an out_edge of the vertex under examination).
  vis.examine_edge(e, g);
  // whenever it is newly decided that an edge is relaxed (has improved the distance for its target)
  vis.edge_relaxed(e, g);
  // whenever a vertex is deemed uninteresting and is taken out of OPEN, but not yet expanded.
  vis.forget_vertex(u, g);
  // whenever a closed vertex becomes INCONS.
  vis.inconsistent_vertex(u, g);
  // notify the visitor that at least one A* round has completed and its resulting path
  // (partial or complete) can be published (the path is encoded in the predecessor
  // property-map).
  vis.publish_path(g);
  // check to see whether the task is finished (return false) or needs to keep going (true).
  { vis.keep_going() } -> std::convertible_to<bool>;
}
&&requires(
    Visitor vis, Graph g,
    std::vector<bagl::graph_edge_descriptor_t<Graph>> vect, double d,
    std::back_insert_iterator<std::vector<bagl::graph_edge_descriptor_t<Graph>>>
        e_it) {
  // ei: back-inserter / forward-iterator for an
  // edge-list. Return the cummulative weight-change and
  // the edge-iterator at the end of the edge list
  // populated by this function.
  std::tie(d, e_it) = vis.detect_edge_change(std::back_inserter(vect), g);
  // adjust the value of epsilon for a given old-value and last cummulative weight-change
  { vis.adjust_epsilon(d, d, g) } -> std::convertible_to<double>;
};

namespace adstar_detail {

BAGL_GRAPH_HAS_MEMBER_FUNCTION(recycle_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(forget_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(inconsistent_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(publish_path)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(keep_going)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(detect_edge_change)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(adjust_epsilon)

template <typename T>
bool invoke_keep_going_or_true(T& t) {
  if constexpr (has_keep_going_v<T&>) {
    return t.keep_going();
  }
  return true;
}
template <typename T>
bool invoke_keep_going_or_true_on_all(T& t) {
  return std::apply(
      [](auto&&... ts) {
        return (invoke_keep_going_or_true(ts) && ... && true);
      },
      t);
}

template <typename T, typename EdgeIter, typename Graph>
std::pair<double, EdgeIter> invoke_detect_edge_change_or_pass(T& t, EdgeIter ei,
                                                              const Graph& g) {
  if constexpr (has_detect_edge_change_v<T&, EdgeIter, const Graph&>) {
    return t.detect_edge_change(ei, g);
  }
  return {0.0, ei};
}
template <typename T, typename EdgeIter, typename Graph>
std::pair<double, EdgeIter> invoke_detect_edge_change_or_pass_on_all(
    T& t, EdgeIter ei, const Graph& g) {
  struct max_accum_t {
    double d = 0.0;
    max_accum_t& operator=(double w) {
      d = std::max(d, w);
      return *this;
    }
  } max_accum;
  std::apply(
      [&max_accum, &ei, &g](auto&&... ts) {
        ((std::tie(max_accum, ei) =
              invoke_detect_edge_change_or_pass(ts, ei, g)),
         ...);
      },
      t);
  return {max_accum.d, ei};
}

template <typename T, typename Graph>
double invoke_adjust_epsilon_or_pass(T& t, double old_eps, double w_change,
                                     const Graph& g) {
  if constexpr (has_detect_edge_change_v<T&, const Graph&>) {
    return t.detect_edge_change(old_eps, w_change, g);
  }
  return old_eps;
}
template <typename T, typename Graph>
double invoke_adjust_epsilon_or_pass_on_all(T& t, double old_eps,
                                            double w_change, const Graph& g) {
  struct min_accum_t {
    double d = old_eps;
    min_accum_t& operator=(double w) {
      d = std::min(d, w);
      return *this;
    }
  } min_accum;
  std::apply(
      [&min_accum, old_eps, w_change, &g](auto&&... ts) {
        ((min_accum = invoke_adjust_epsilon_or_pass(ts, old_eps, w_change, g)),
         ...);
      },
      t);
  return min_accum.d;
}

}  // namespace adstar_detail

/**
  * This class is the default implementation of an AD* visitor (see ADStarVisitor).
  * Basically, this implementation models the concept required by AD*, but does nothing at all
  * (all functions are empty).
  */
template <typename Visitors = bagl::null_visitors>
class adstar_visitor {
 public:
  adstar_visitor() = default;
  explicit adstar_visitor(Visitors&& vis) : vis_(std::move(vis)) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, const Graph& g) const {
    bagl::visitors_detail::invoke_initialize_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, const Graph& g) const {
    bagl::visitors_detail::invoke_discover_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void inconsistent_vertex(Vertex u, const Graph& g) const {
    adstar_detail::invoke_inconsistent_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex u, const Graph& g) const {
    bagl::visitors_detail::invoke_examine_vertex_on_all(vis_, u, g);
  }
  template <typename Edge, typename Graph>
  void examine_edge(Edge e, const Graph& g) const {
    bagl::visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_relaxed(Edge e, const Graph& g) const {
    bagl::visitors_detail::invoke_edge_relaxed_on_all(vis_, e, g);
  }
  template <typename Vertex, typename Graph>
  void forget_vertex(Vertex u, const Graph& g) const {
    adstar_detail::invoke_forget_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, const Graph& g) const {
    bagl::visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void recycle_vertex(Vertex u, const Graph& g) const {
    adstar_detail::invoke_recycle_vertex_on_all(vis_, u, g);
  }

  template <typename Graph>
  void publish_path(const Graph& g) const {
    adstar_detail::invoke_publish_path_on_all(vis_, g);
  }

  bool keep_going() const {
    return adstar_detail::invoke_keep_going_or_true_on_all(vis_);
  }

  template <typename EdgeIter, typename Graph>
  std::pair<double, EdgeIter> detect_edge_change(EdgeIter ei,
                                                 const Graph& g) const {
    return adstar_detail::invoke_detect_edge_change_or_pass_on_all(vis_, ei, g);
  }

  template <typename Graph>
  double adjust_epsilon(double old_eps, double w_change, const Graph& g) const {
    return adstar_detail::invoke_adjust_epsilon_or_pass_on_all(vis_, old_eps,
                                                               w_change, g);
  }

 protected:
  Visitors vis_;
};

template <typename... Visitors>
auto make_adstar_visitor(Visitors&&... vis) {
  if constexpr (sizeof...(Visitors) == 0) {
    return adstar_visitor<>();
  } else {
    return adstar_visitor(
        std::tuple<std::decay_t<Visitors>...>(std::forward<Visitors>(vis)...));
  }
}
using default_adstar_visitor = adstar_visitor<>;

/**
  * This class template is used by the AD* algorithm to constitute the key-values which
  * drive the ordering in the priority-queue that the AD* algorithm uses to choose the next
  * vertex to examine. This class simply implements the key-value described in the original
  * paper on AD* (ICAPS 2005).
  * \tparam DistanceValueType The scalar type that describes the distance values.
  * \tparam CompareFunction The strict weak-ordering function that can sort the distance values.
  * \tparam EqualCompareFunction The equal-comparison function that can compare two distance values to be equal.
  */
template <typename DistanceValueType, typename CompareFunction = std::less<>,
          typename EqualCompareFunction = std::equal_to<>>
struct adstar_key_value {
  /**
    * Default constructor.
    */
  adstar_key_value() : k1_(0), k2_(0), compare_(), equal_() {}
  /**
    * Parametrized constructor.
    * \param k1 The first value of the key-value.
    * \param k2 The second value of the key-value.
    * \param compare The functor of type CompareFunction.
    * \param equalCompare The functor of type EqualCompareFunction.
    */
  adstar_key_value(DistanceValueType k1, DistanceValueType k2,
                   CompareFunction compare = CompareFunction(),
                   EqualCompareFunction equal_compare = EqualCompareFunction())
      : k1_(k1), k2_(k2), compare_(compare), equal_(equal_compare) {}

  /**
    * The less-than operator for strict weak-ordering.
    */
  bool operator<(const adstar_key_value<DistanceValueType, CompareFunction,
                                        EqualCompareFunction>& k) const {
    return compare_(k1_, k.k1_) || (equal_(k1_, k.k1_) && compare_(k2_, k.k2_));
  }

  DistanceValueType k1_;
  DistanceValueType k2_;
  CompareFunction compare_;
  EqualCompareFunction equal_;
};

/**
  * This traits class defines that traits that an AD* key-value should have.
  * \tparam ADStarKeyType The key-value type of which the traits are sought.
  */
template <typename ADStarKeyType>
struct adstar_key_traits {
  /** The type of comparison to use for strict weak-ordering of the key-values. */
  using compare_type = std::less<ADStarKeyType>;
};

namespace adstar_detail {

template <typename AStarHeuristicMap, typename UniformCostVisitor,
          typename UpdatableQueue, typename List, typename PredecessorMap,
          typename KeyMap, typename DistanceMap, typename RHSMap,
          typename WeightMap, typename ColorMap, typename CompareFunction,
          typename EqualCompareFunction, typename CombineFunction,
          typename ComposeFunction>
struct adstar_bfs_visitor {

  using KeyValue = bagl::property_traits_value_t<KeyMap>;
  using ColorValue = bagl::property_traits_value_t<ColorMap>;
  using Color = bagl::color_traits<ColorValue>;
  using distance_type = bagl::property_traits_value_t<DistanceMap>;
  using weight_type = bagl::property_traits_value_t<WeightMap>;

  adstar_bfs_visitor(AStarHeuristicMap h, UniformCostVisitor vis,
                     UpdatableQueue& q, List& incons, PredecessorMap p,
                     KeyMap k, DistanceMap d, RHSMap rhs, WeightMap w,
                     ColorMap col, distance_type& epsilon,
                     CompareFunction compare,
                     EqualCompareFunction equal_compare,
                     CombineFunction combine, ComposeFunction compose,
                     distance_type inf, distance_type zero)
      : h_(h),
        vis_(vis),
        q_(q),
        incons_(incons),
        predecessor_(p),
        key_(k),
        distance_(d),
        rhs_(rhs),
        weight_(w),
        color_(col),
        epsilon_(epsilon),
        compare_(compare),
        equal_compare_(equal_compare),
        combine_(combine),
        compose_(compose),
        inf_(inf),
        zero_(zero) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) const {
    vis_.initialize_vertex(u, g);
  }
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, Graph& g) const {
    vis_.discover_vertex(u, g);
  }
  template <typename Vertex, typename Graph>
  void inconsistent_vertex(Vertex u, Graph& g) const {
    vis_.inconsistent_vertex(u, g);
  }
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex u, Graph& g) const {
    vis_.examine_vertex(u, g);
  }
  template <typename Edge, typename Graph>
  void examine_edge(Edge e, Graph& g) const {
    if (compare_(get(weight_, e), zero_)) {
      throw bagl::negative_edge();
    }
    vis_.examine_edge(e, g);
  }
  template <typename Vertex, typename Graph>
  void forget_vertex(Vertex u, Graph& g) const {
    vis_.forget_vertex(u, g);
  }
  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, Graph& g) const {
    vis_.finish_vertex(u, g);
  }
  template <typename Vertex, typename Graph>
  void recycle_vertex(Vertex u, Graph& g) const {
    vis_.recycle_vertex(u, g);
  }

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
  double adjust_epsilon(double old_eps, double w_change, const Graph& g) const {
    return vis_.adjust_epsilon(old_eps, w_change, g);
  }

  template <typename Vertex, typename Graph>
  void update_key(Vertex u, Graph& /*unused*/) {
    distance_type g_u = get(distance_, u);
    distance_type rhs_u = get(rhs_, u);
    if (compare_(rhs_u, g_u)) {
      put(key_, u,
          KeyValue(combine_(rhs_u, compose_(epsilon_, get(h_, u))), rhs_u,
                   compare_, equal_compare_));
    } else {
      put(key_, u,
          KeyValue(combine_(g_u, get(h_, u)), g_u, compare_, equal_compare_));
    }
  }

  template <typename Vertex, bagl::concepts::BidirectionalGraph G>
  void update_vertex(Vertex u, G& g) {
    ColorValue col_u = get(color_, u);

    if (col_u == Color::white()) {
      put(distance_, u, inf_);
      col_u = Color::green();
      put(color_, u, col_u);
    }

    distance_type g_u = get(distance_, u);
    distance_type rhs_u = get(rhs_, u);

    if (!equal_compare_(rhs_u, zero_)) {  // if u is not the start node.
      rhs_u = inf_;                       // This was in the original code!
      for (auto e : in_edges(u, g)) {
        distance_type rhs_tmp =
            combine_(get(weight_, e), get(distance_, source(e, g)));
        if (compare_(rhs_tmp, rhs_u)) {
          rhs_u = rhs_tmp;
          put(rhs_, u, rhs_u);  // this was the original code!
          put(predecessor_, u, source(e, g));
        }
      }
      rhs_u = get(rhs_, u);
    }

    if (!equal_compare_(rhs_u, g_u)) {
      // if not in CLOSED set (i.e. either just closed (black) or inconsistent (red)).
      if ((col_u != Color::black()) && (col_u != Color::red())) {
        update_key(u, g);
        q_.push_or_update(u);
        put(color_, u, Color::gray());
        vis_.discover_vertex(u, g);
      } else if (col_u == Color::black()) {
        incons_.push_back(u);
        put(color_, u, Color::red());
        vis_.inconsistent_vertex(u, g);
      }
    } else if (q_.contains(u)) {  // if u is in the OPEN set, then remove it.
      put(key_, u, KeyValue(-inf_, -inf_, compare_, equal_compare_));
      q_.update(u);
      q_.pop();          // remove from OPEN set
      update_key(u, g);  // this was the original code!
      put(color_, u, Color::green());
      vis_.forget_vertex(u, g);
    }
  }

  AStarHeuristicMap h_;
  UniformCostVisitor vis_;

  UpdatableQueue& q_;
  List& incons_;

  PredecessorMap predecessor_;
  KeyMap key_;
  DistanceMap distance_;
  RHSMap rhs_;
  WeightMap weight_;
  ColorMap color_;

  distance_type& epsilon_;
  CompareFunction compare_;
  EqualCompareFunction equal_compare_;
  CombineFunction combine_;
  ComposeFunction compose_;
  distance_type inf_;
  distance_type zero_;
};

template <
    // this is the actual graph, should comply to BidirectionalGraphConcept.
    typename G,
    // this the map of heuristic function value for each vertex.
    typename AStarHeuristicMap,
    // this is a visitor class that can perform special operations at event points.
    typename ADStarBFSVisitor,
    // this is the map that stores the preceeding edge for each vertex.
    typename PredecessorMap,
    // this is the map of distance values associated with each vertex.
    typename DistanceMap, typename RHSMap,
    // this is the map of key values associated to each vertex.
    typename KeyMap,
    // this is the map of edge weight (or cost) associated to each edge of the graph.
    typename WeightMap,
    // this is a color map for each vertex, i.e. white=not visited, gray=discovered, black=expanded.
    typename ColorMap, typename IndexInHeapMap, typename MutableQueue,
    typename InconsList,
    // a binary comparison function object that returns true if the first operand is
    // strictly better (less-than) than the second operand.
    typename CompareFunction /*= std::less<>*/,
    // a binary comparison function object that returns true if both operands are
    // equal to each other.
    typename EqualCompareFunction /*= std::equal_to<>*/,
    // a binary combination function object that returns the sum of its
    // operands (sum in the broad sense).
    typename CombineFunction /*= std::plus<>*/,
    // a binary composition function object that amplifies a heuristic distance
    // metric by a scalar value (i.e. epsilon x h(u)).
    typename ComposeFunction /*= std::multiplies<>*/>
void adstar_search_loop(
    G& g, bagl::graph_vertex_descriptor_t<G> start_vertex,
    AStarHeuristicMap hval, ADStarBFSVisitor& bfs_vis,
    PredecessorMap predecessor, DistanceMap distance, RHSMap rhs, KeyMap key,
    WeightMap weight, ColorMap color, IndexInHeapMap index_in_heap,
    MutableQueue& q, InconsList& incons,
    bagl::property_traits_value_t<DistanceMap>& epsilon,
    bagl::property_traits_value_t<DistanceMap> inf,
    bagl::property_traits_value_t<DistanceMap> zero =
        bagl::property_traits_value_t<DistanceMap>(0),
    CompareFunction compare = CompareFunction(),
    EqualCompareFunction equal_compare = EqualCompareFunction(),
    CombineFunction combine = CombineFunction(),
    ComposeFunction compose = ComposeFunction()) {
  using Edge = bagl::graph_edge_descriptor_t<G>;
  using ColorValue = bagl::property_traits_value_t<ColorMap>;
  using Color = bagl::color_traits<ColorValue>;

  auto s = start_vertex;
  put(distance, s, inf);
  put(rhs, s, zero);
  put(predecessor, s, s);
  bfs_vis.update_key(s, g);
  put(color, s, Color::gray());
  bfs_vis.discover_vertex(s, g);
  q.push(s);

  std::vector<Edge> affected_edges;

  while (bfs_vis.keep_going()) {

    while (!q.empty()) {
      auto u = q.top();
      q.pop();
      bfs_vis.examine_vertex(u, g);
      auto dist_g_u = get(distance, u);
      auto dist_rhs_u = get(rhs, u);
      // if we have a consistent node at the goal
      if (equal_compare(get(hval, u), zero) &&
          equal_compare(dist_g_u, dist_rhs_u)) {
        break;
      }
      // if dist_g_u is greater than dist_rhs_u, then make u consistent and close it.
      if (compare(dist_rhs_u, dist_g_u)) {
        put(distance, u, dist_rhs_u);
        dist_g_u = dist_rhs_u;
        put(color, u, Color::black());
        bfs_vis.finish_vertex(u, g);
      } else {
        put(distance, u, inf);
        bfs_vis.update_vertex(u, g);
      }
      for (auto e : out_edges(u, g)) {
        bfs_vis.examine_edge(e, g);
        bfs_vis.update_vertex(target(e, g), g);
      }
    }

    bfs_vis.publish_path(g);

    affected_edges.clear();
    auto max_w_change =
        bfs_vis.detect_edge_change(std::back_inserter(affected_edges), g).first;

    // update all nodes that were affected.
    for (auto e : affected_edges) {
      if (get(color, source(e, g)) == Color::black()) {
        put(color, source(e, g), Color::green());
      }
      bfs_vis.update_vertex(source(e, g), g);
      if (get(color, target(e, g)) == Color::black()) {
        put(color, target(e, g), Color::green());
      }
      bfs_vis.update_vertex(target(e, g), g);
    }

    epsilon = bfs_vis.adjust_epsilon(epsilon, max_w_change, g);

    // merge the OPEN and INCONS sets
    for (auto v : incons) {
      bfs_vis.update_key(v, g);
      q.push_or_update(v);
      put(color, v, Color::gray());
    }
    incons.clear();

    // update keys for all OPEN nodes, and change all black nodes to green (empty the CLOSED set).
    for (auto u : vertices(g)) {
      ColorValue u_color = get(color, u);
      if (q.contains(u)) {
        bfs_vis.update_key(u, g);
        q.update(u);
        put(color, u, Color::gray());
        bfs_vis.discover_vertex(u, g);
      } else if (u_color == Color::black()) {
        put(color, u, Color::green());
        bfs_vis.recycle_vertex(u, g);
      }
    }
  }
}

}  // namespace adstar_detail

/**
  * This function template performs an AD* search over a graph, without initialization. The AD* search
  * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
  * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then
  * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it
  * uses callback functions to allow the user to check for changes in the environment, triggering an update
  * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
  * \tparam VertexListGraph The type of the graph on which the search is performed, should model
  *         the BGL's BidirectionalGraphConcept.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
  *         for each vertex in the graph.
  * \tparam Visitor The type of the AD* visitor to be used.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex together with its optimal predecessor.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each
  *         vertex to the goal.
  * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of
  *         each vertex to the goal (internal use to AD*).
  * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the
  *         graph (cost of travel along an edge).
  * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors
  *         are used to mark vertices by their status in the AD* algorithm (white = not visited,
  *         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED,
  *         not OPEN), red = inconsistent (in INCONS)).
  * \tparam CompareFunction A binary comparison functor type that returns true if the first operand
  *         is strictly better (less-than) than the second operand.
  * \tparam EqualCompareFunction A binary comparison functor type that returns true if both operands
  *         are equal to each other.
  * \tparam CombineFunction A binary combination functor type that returns the sum of its operands,
  *         semantically-speaking.
  * \tparam ComposeFunction A binary composition functor type that amplifies a heuristic distance
  *         metric by a scalar value (i.e. epsilon x h(u) ).
  *
  * \param g The graph on which to apply the AD* algorithm.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param hval The property-map of A* heuristic function values for each vertex.
  * \param vis The AD* visitor object.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
  *        goal (for internal use).
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param color The property-map which stores the color-value of the vertices, colors are used to mark
  *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
  *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  * \param inf The quantity that represents infinity (either a very large value or the infinity value for
  *        the underlying value-type).
  * \param zero The quantity that represents zero with the given value-type.
  * \param compare A binary comparison functor that returns true if the first operand is strictly
  *        better (less-than) than the second operand.
  * \param equal_compare A binary comparison functor that returns true if both operands are equal
  *        to each other.
  * \param combine A binary combination functor that returns the sum of its operands,
  *        semantically-speaking.
  * \param compose A binary composition functor that amplifies a heuristic distance metric by a
  *        scalar value (i.e. epsilon x h(u) ).
  */
template <bagl::concepts::VertexListGraph G,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          ADStarVisitor<G> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> KeyMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          typename CompareFunction, typename EqualCompareFunction,
          typename CombineFunction, typename ComposeFunction>
void adstar_search_no_init(
    G& g, bagl::graph_vertex_descriptor_t<G> start_vertex,
    AStarHeuristicMap hval, Visitor vis, PredecessorMap predecessor,
    DistanceMap distance, RHSMap rhs, KeyMap key, WeightMap weight,
    ColorMap color, bagl::property_traits_value_t<DistanceMap> epsilon,
    bagl::property_traits_value_t<DistanceMap> inf,
    bagl::property_traits_value_t<DistanceMap> zero =
        bagl::property_traits_value_t<DistanceMap>(0),
    CompareFunction compare = CompareFunction(),
    EqualCompareFunction equal_compare = EqualCompareFunction(),
    CombineFunction combine = CombineFunction(),
    ComposeFunction compose = ComposeFunction()) {
  using Vertex = bagl::graph_vertex_descriptor_t<G>;
  using KeyValue = bagl::property_traits_value_t<KeyMap>;
  using KeyCompareType = typename adstar_key_traits<KeyValue>::compare_type;

  // priority queue holding the OPEN set.
  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());
  auto q = make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                               KeyCompareType());
  // list holding the INCONS set (inconsistent nodes).
  std::vector<Vertex> incons;

  adstar_detail::adstar_bfs_visitor bfs_vis(
      hval, vis, q, incons, predecessor, key, distance, rhs, weight, color,
      epsilon, compare, equal_compare, combine, compose, inf, zero);

  adstar_detail::adstar_search_loop(
      g, start_vertex, hval, bfs_vis, predecessor, distance, rhs, key, weight,
      color, index_in_heap, q, incons, epsilon, inf, zero, compare,
      equal_compare, combine, compose);
}

/**
  * This function template performs an AD* search over a graph, without initialization. The AD* search
  * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
  * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then
  * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it
  * uses callback functions to allow the user to check for changes in the environment, triggering an update
  * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
  * \tparam VertexListGraph The type of the graph on which the search is performed, should model
  *         the BGL's BidirectionalGraphConcept.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
  *         for each vertex in the graph.
  * \tparam Visitor The type of the AD* visitor to be used.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex together with its optimal predecessor.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each
  *         vertex to the goal.
  * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of
  *         each vertex to the goal (internal use to AD*).
  * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the
  *         graph (cost of travel along an edge).
  * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors
  *         are used to mark vertices by their status in the AD* algorithm (white = not visited,
  *         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED,
  *         not OPEN), red = inconsistent (in INCONS)).
  *
  * \param g The graph on which to apply the AD* algorithm.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param hval The property-map of A* heuristic function values for each vertex.
  * \param vis The AD* visitor object.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
  *        goal (for internal use).
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param color The property-map which stores the color-value of the vertices, colors are used to mark
  *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
  *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  */
template <bagl::concepts::VertexListGraph G,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          ADStarVisitor<G> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> KeyMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void adstar_search_no_init(G& g,
                           bagl::graph_vertex_descriptor_t<G> start_vertex,
                           AStarHeuristicMap hval, Visitor vis,
                           PredecessorMap predecessor, DistanceMap distance,
                           RHSMap rhs, KeyMap key, WeightMap weight,
                           ColorMap color, double epsilon) {
  adstar_search_no_init(
      g, start_vertex, hval, vis, predecessor, distance, rhs, key, weight,
      color, epsilon, std::numeric_limits<double>::infinity(), 0.0,
      std::less<>(), std::equal_to<>(), std::plus<>(), std::multiplies<>());
}

/**
  * This function template performs an AD* search over a graph. The AD* search
  * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
  * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then
  * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it
  * uses callback functions to allow the user to check for changes in the environment, triggering an update
  * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
  * \tparam VertexListGraph The type of the graph on which the search is performed, should model the BGL's
  *BidirectionalGraphConcept.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values for each vertex in
  *the graph.
  * \tparam Visitor The type of the AD* visitor to be used.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting vertex together with
  *its optimal predecessor.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex to the goal.
  * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of each vertex to the
  *goal (internal use to AD*).
  * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors are used to mark
  *vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), black = finished (in
  *CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \tparam CompareFunction A binary comparison functor type that returns true if the first operand is strictly better
  *(less-than) than the second operand.
  * \tparam EqualCompareFunction A binary comparison functor type that returns true if both operands are equal to each
  *other.
  * \tparam CombineFunction A binary combination functor type that returns the sum of its operands,
  *semantically-speaking.
  * \tparam ComposeFunction A binary composition functor type that amplifies a heuristic distance metric by a scalar
  *value (i.e. epsilon x h(u) ).
  *
  * \param g The graph on which to apply the AD* algorithm.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param hval The property-map of A* heuristic function values for each vertex.
  * \param vis The AD* visitor object.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
  *        goal (for internal use).
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param color The property-map which stores the color-value of the vertices, colors are used to mark
  *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
  *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  * \param inf The quantity that represents infinity (either a very large value or the infinity value for
  *        the underlying value-type).
  * \param zero The quantity that represents zero with the given value-type.
  * \param compare A binary comparison functor that returns true if the first operand is strictly better (less-than)
  *than the second operand.
  * \param equal_compare A binary comparison functor that returns true if both operands are equal to each other.
  * \param combine A binary combination functor that returns the sum of its operands, semantically-speaking.
  * \param compose A binary composition functor that amplifies a heuristic distance metric by a scalar value (i.e.
  *epsilon x h(u) ).
  */
template <bagl::concepts::VertexListGraph G,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          ADStarVisitor<G> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> KeyMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          typename CompareFunction, typename EqualCompareFunction,
          typename CombineFunction, typename ComposeFunction>
void adstar_search(G& g, bagl::graph_vertex_descriptor_t<G> start_vertex,
                   AStarHeuristicMap hval, Visitor vis,
                   PredecessorMap predecessor, DistanceMap distance, RHSMap rhs,
                   KeyMap key, WeightMap weight, ColorMap color,
                   bagl::property_traits_value_t<DistanceMap> epsilon,
                   bagl::property_traits_value_t<DistanceMap> inf,
                   bagl::property_traits_value_t<DistanceMap> zero =
                       bagl::property_traits_value_t<DistanceMap>(0),
                   CompareFunction compare = CompareFunction(),
                   EqualCompareFunction equal_compare = EqualCompareFunction(),
                   CombineFunction combine = CombineFunction(),
                   ComposeFunction compose = ComposeFunction()) {

  using ColorValue = bagl::property_traits_value_t<ColorMap>;
  using Color = bagl::color_traits<ColorValue>;
  using KeyValue = bagl::property_traits_value_t<KeyMap>;
  for (auto u : vertices(g)) {
    put(color, u, Color::white());
    put(distance, u, inf);
    put(rhs, u, inf);
    put(key, u, KeyValue(inf, inf, compare, equal_compare));
    put(predecessor, u, u);
    vis.initialize_vertex(u, g);
  }

  adstar_search_no_init(g, start_vertex, hval, vis, predecessor, distance, rhs,
                        key, weight, color, epsilon, inf, zero, compare,
                        equal_compare, combine, compose);
}

/**
  * This function template performs an AD* search over a graph. The AD* search
  * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
  * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then
  * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it
  * uses callback functions to allow the user to check for changes in the environment, triggering an update
  * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
  * \tparam VertexListGraph The type of the graph on which the search is performed, should model the BGL's
  *BidirectionalGraphConcept.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values for each vertex in
  *the graph.
  * \tparam Visitor The type of the AD* visitor to be used.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting vertex together with
  *its optimal predecessor.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex to the goal.
  * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of each vertex to the
  *goal (internal use to AD*).
  * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors are used to mark
  *vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), black = finished (in
  *CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  *
  * \param g The graph on which to apply the AD* algorithm.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param hval The property-map of A* heuristic function values for each vertex.
  * \param vis The AD* visitor object.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the
  *        goal (for internal use).
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param color The property-map which stores the color-value of the vertices, colors are used to mark
  *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN),
  *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  */
template <bagl::concepts::VertexListGraph G,
          bagl::concepts::ReadableVertexPropertyMap<G> AStarHeuristicMap,
          ADStarVisitor<G> Visitor,
          bagl::concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> RHSMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> KeyMap,
          bagl::concepts::ReadableEdgePropertyMap<G> WeightMap,
          bagl::concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void adstar_search(G& g, bagl::graph_vertex_descriptor_t<G> start_vertex,
                   AStarHeuristicMap hval, Visitor vis,
                   PredecessorMap predecessor, DistanceMap distance, RHSMap rhs,
                   KeyMap key, WeightMap weight, ColorMap color,
                   double epsilon) {
  adstar_search(g, start_vertex, hval, vis, predecessor, distance, rhs, key,
                weight, color, epsilon, std::numeric_limits<double>::infinity(),
                0.0, std::less<>(), std::equal_to<>(), std::plus<>(),
                std::multiplies<>());
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_ADSTAR_SEARCH_H_
