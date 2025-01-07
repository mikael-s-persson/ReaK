/**
 * \file dvp_tree_detail.h
 *
 * This library implements the details of a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space, with amortized O(logN)
 * insertion-deletion. A DVP-tree is essentially a generalization of a search tree which only
 * requires the space to have a metric which respects the triangular inequality.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PLANNING_PATH_PLANNING_DVP_TREE_DETAIL_H_
#define REAK_PLANNING_PATH_PLANNING_DVP_TREE_DETAIL_H_

#include "bagl/graph_concepts.h"
#include "bagl/property_map.h"
#include "bagl/tree_adaptor.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/proper_metric_concept.h"

#include "ReaK/core/base/global_rng.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <stack>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ReaK::pp {

/**
 * This class is a callable class that can be used to choose the best
 * vantage-point to use out of a set of points. In theory, the best vantage-point
 * is the one which deviates the most from the other points in the set, however,
 * this functor will approximately select that point by searching for it only
 * in a random subset of the given range of points.
 */
class random_best_vp_chooser {
 private:
  unsigned int divider_ = 10;

 public:
  /**
   * Default construction.
   * \param aDivider The divider of the set (determines the fraction of the points to search), default is 10.
   */
  explicit random_best_vp_chooser(unsigned int divider) : divider_(divider) {}

  random_best_vp_chooser() : random_best_vp_chooser(10) {}

  /**
   * This call-operator will choose a vantage-point from within the given range.
   * \tparam RandomAccessIter A random-access iterator type that can describe the point-range.
   * \tparam DistanceMetric The distance-metric type over the topology type.
   * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of
   * the iterators) to a point (position).
   * \param first The start of the range of vertices.
   * \param last The end of the range of vertices (one element past the end).
   * \param distance The distance-value functor.
   * \param position The property-map used to obtain the positions from the vertices.
   * \return A random-access iterator to the chosen vantage-point.
   */
  template <typename RandomAccessIter, typename DistanceMetric,
            typename PositionMap>
  RandomAccessIter operator()(RandomAccessIter first, RandomAccessIter last,
                              DistanceMetric distance,
                              PositionMap position) const {
    RandomAccessIter best_pt = last;
    double best_dev = -1;
    for (int i = 0; i < (last - first) / divider_ + 1; ++i) {
      RandomAccessIter current_pt =
          first + (get_global_rng()() % (last - first));
      double current_mean = 0.0;
      double current_dev = 0.0;
      auto current_vp = get(position, *current_pt);
      for (int j = 0; first + j != last; ++j) {
        double dist = distance(current_vp, get(position, *(first + j)));
        current_mean = (current_mean * j + dist) / (j + 1);
        current_dev = (current_dev * j + dist * dist) / (j + 1);
      }
      double current_var = current_dev - current_mean * current_mean;
      if (current_var < 0) {
        current_var = 0.0;
      }
      current_dev = std::sqrt(current_var);

      if (current_dev > best_dev) {
        best_pt = current_pt;
        best_dev = current_dev;
      }
    }
    return best_pt;
  }
};

/**
 * This class is a callable class that can be used to choose the
 * vantage-point to use out of a set of points. This functor will
 * select a random point from the set.
 */
class random_vp_chooser {
 public:
  /**
   * Default construction.
   */
  random_vp_chooser() = default;

  /**
   * This call-operator will choose a vantage-point from within the given range.
   * \tparam RandomAccessIter A random-access iterator type that can describe the point-range.
   * \tparam DistanceMetric The distance-metric type over the topology type.
   * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of
   * the iterators) to a point (position).
   * \param first The start of the range of vertices.
   * \param last The end of the range of vertices (one element past the end).
   * \param distance The distance-metric over the given topology.
   * \param position The property-map used to obtain the positions from the vertices.
   * \return A random-access iterator to the chosen vantage-point.
   */
  template <typename RandomAccessIter, typename DistanceMetric,
            typename PositionMap>
  RandomAccessIter operator()(RandomAccessIter first, RandomAccessIter last,
                              DistanceMetric /*distance*/,
                              PositionMap /*position*/) const {
    return first + (get_global_rng()() % (last - first));
  }
};

/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space, with amortized O(logN)
 * insertion-deletion. A DVP-tree is essentially a generalization of a search tree which only
 * requires the space to have a metric which respects the triangular inequality.
 * \tparam TreeType The tree type to be used to store the entries of this DVP search tree.
 * \tparam Space The topology type on which the points associated to each entry resides.
 * \tparam VertexKeyMap The property-map type that can map vertex properties of the tree to key-values (that identify
 * the entries).
 * \tparam DistanceMap The property-map type that can map an edge property to its associated distance value (used
 * internally).
 * \tparam PositionMap The property-map type that can map vertex properties of the tree to associated position values
 * (positions in the topology).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename TreeType, MetricSpace Space, typename VertexPropMap,
          typename VertexKeyMap, typename DistanceMap, typename PositionMap,
          unsigned int Arity, typename VPChooser>
class dvp_tree_impl {
 public:
  using self = dvp_tree_impl<TreeType, Space, VertexPropMap, VertexKeyMap,
                             DistanceMap, PositionMap, Arity, VPChooser>;

  /** Type of the points in the Space. */
  using point_type = bagl::property_traits_value_t<PositionMap>;
  /** Type of the distance values. */
  using distance_type = double;

  using distance_metric_type = metric_space_distance_metric_t<Space>;
  using proper_metric_type = get_proper_metric_t<Space>;

  struct dist_metric_only_impl {
    std::shared_ptr<const Space> space_ptr;
    distance_metric_type distance_func;

    double distance(const point_type& a, const point_type& b) const {
      return distance_func(a, b, *space_ptr);
    }
    double proper_distance(const point_type& a, const point_type& b) const {
      return distance_func(a, b, *space_ptr);
    }
    double operator()(const point_type& a, const point_type& b) const {
      return proper_distance(a, b, *space_ptr);
    }

    explicit dist_metric_only_impl(const std::shared_ptr<const Space>& space)
        : space_ptr(space), distance_func(get(distance_metric, *space)) {}
  };

  struct dist_metric_pair_impl {
    std::shared_ptr<const Space> space_ptr;
    distance_metric_type distance_func;
    proper_metric_type proper_distance_func;

    double distance(const point_type& a, const point_type& b) const {
      return distance_func(a, b, *space_ptr);
    }
    double proper_distance(const point_type& a, const point_type& b) const {
      return proper_distance_func(a, b, *space_ptr);
    }
    double operator()(const point_type& a, const point_type& b) const {
      return proper_distance(a, b, *space_ptr);
    }

    explicit dist_metric_pair_impl(const std::shared_ptr<const Space>& space)
        : space_ptr(space),
          distance_func(get(distance_metric, *space)),
          proper_distance_func(get(proper_metric, *space)) {}
  };

  using parting_metrics_type = std::conditional_t<
      std::is_same_v<distance_metric_type, proper_metric_type>,
      dist_metric_only_impl, dist_metric_pair_impl>;

  /** Type of the key-values that identify entries of the DVP tree. */
  using key_type = bagl::property_traits_value_t<VertexKeyMap>;

  static_assert(DistanceMetric<distance_metric_type, Space>);
  static_assert(DistanceMetric<proper_metric_type, Space>);

 private:
  using tree_indexer = TreeType;

  using vertex_type = bagl::graph_vertex_descriptor_t<tree_indexer>;
  using edge_type = bagl::graph_edge_descriptor_t<tree_indexer>;

  using vertex_property = bagl::property_traits_value_t<VertexPropMap>;
  using edge_property = typename tree_indexer::edge_property_type;

  struct priority_compare_type {
    bool operator()(const std::pair<distance_type, vertex_type>& x,
                    const std::pair<distance_type, vertex_type>& y) const {
      return (x.first < y.first);
    }
  };
  using priority_queue_type =
      std::vector<std::pair<distance_type, vertex_type>>;

  tree_indexer* tree_;
  vertex_type root_;

  // A map from a vertex_descriptor to a vertex_property.
  VertexPropMap vprop_;
  // A map from a vertex_property to a key_type value.
  VertexKeyMap key_;
  // A map from an edge-property to a distance value.
  DistanceMap mu_;
  // A map from a vertex_property to a position value (should be Read-Write).
  PositionMap position_;

  /// The distance-metric functor.
  parting_metrics_type distance_;

  /// The vantage-point chooser (functor).
  VPChooser vp_chooser_;

  using prop_vector_iter = typename std::vector<vertex_property>::iterator;

  struct construction_task {
    vertex_type node;
    prop_vector_iter first;
    prop_vector_iter last;
  };

  prop_vector_iter rearrange_with_chosen_vp(prop_vector_iter first,
                                            prop_vector_iter last) const {
    // choose a vantage-point in the interval:
    auto chosen_vp_it = vp_chooser_(first, last, distance_, position_);
    if (chosen_vp_it == last) {
      // no vp to be chosen in this interval (presumably, empty interval).
      return first;
    }

    // place the vantage-point (or pivot) at the start of interval:
    std::iter_swap(chosen_vp_it, first);
    chosen_vp_it = first;  //<-- chosen_vp_it is now at the start of interval.

    return chosen_vp_it;
  }

  /* NOTE Invalidates vertices */
  /* NOTE This is a non-recursive version of the construct-node algorithm */
  /* Does not require persistent vertices */
  /* This is the main tree construction function. It takes the vertices in the iterator range and organizes them
   * as a sub-tree below the aParentNode node (and aEdgeDist is the minimum distance to the parent of any node in the
   * range). */
  // node is a valid, existing node containing the chosen vantage-point.
  void construct_node(vertex_type node, prop_vector_iter first,
                      prop_vector_iter last) {

    std::unordered_map<key_type, distance_type,
                       bagl::graph_descriptor_hash_t<key_type>>
        dist_map;
    std::queue<construction_task> tasks;
    tasks.emplace(node, first, last);

    auto is_closer = [&](const auto& vp1, const auto& vp2) {
      return dist_map[get(key_, vp1)] < dist_map[get(key_, vp2)];
    };

    while (!tasks.empty()) {
      construction_task cur_task = tasks.front();
      tasks.pop();

      // update values in the dist-map with the distances to the new chosen vantage-point:
      const point_type& chosen_vp_pt =
          get(position_, get(vprop_, cur_task.node));
      for (auto it = cur_task.first; it != cur_task.last; ++it) {
        dist_map[get(key_, *it)] =
            distance_.proper_distance(chosen_vp_pt, get(position_, *it));
      }

      // this loop splits up the children into as equal as possible partitions.
      std::size_t total_count = (cur_task.last - cur_task.first);
      std::array<std::size_t, Arity> child_count = {};
      for (std::size_t i = Arity; i > 0; --i) {
        child_count[i - 1] = total_count / i;
        total_count -= child_count[i - 1];
      }
      for (std::size_t i = 0; (i < Arity) && (child_count[i] > 0); ++i) {
        std::nth_element(cur_task.first, cur_task.first + (child_count[i] - 1),
                         cur_task.last, is_closer);
        auto temp = cur_task.first;
        cur_task.first += child_count[i];
        edge_property ep;
        put(mu_, ep, dist_map[get(key_, *(cur_task.first - 1))]);
        rearrange_with_chosen_vp(temp, cur_task.first);
        dist_map.erase(get(key_, *temp));
        auto [new_vp_node, new_ep_node, new_added] =
            add_child(cur_task.node, *tree_, std::move(*temp), std::move(ep));
        ++temp;
        if (temp != cur_task.first) {
          tasks.emplace(new_vp_node, temp, cur_task.first);
        }
      }
    }
  }

  struct nearest_search_result_set {
    priority_queue_type neighbors;
    std::size_t num_neighbors;
    distance_type radius;

    nearest_search_result_set(std::size_t k, distance_type rad)
        : neighbors(), num_neighbors(k), radius(rad) {}

    void register_vantage_point(const point_type& /*unused*/,
                                const point_type& /*unused*/,
                                distance_type current_dist,
                                vertex_type current_vp,
                                const parting_metrics_type& /*unused*/) {
      // is the vantage point within current search bound?
      if (current_dist < radius) {
        // then add the vantage point to the NN list.
        neighbors.emplace_back(current_dist, current_vp);
        std::push_heap(neighbors.begin(), neighbors.end(),
                       priority_compare_type());
        // are there too many nearest neighbors?
        if (neighbors.size() > num_neighbors) {
          std::pop_heap(neighbors.begin(), neighbors.end(),
                        priority_compare_type());
          // delete last element to keep neighbors with num_neighbors elements
          neighbors.pop_back();
          // distance of the last element is now the search bound radius.
          radius = neighbors.front().first;
        }
      }
    }
  };

  struct pred_succ_search_result_set {
    priority_queue_type pred;
    priority_queue_type succ;
    std::size_t num_neighbors;
    distance_type radius;
    distance_type sigma_pred;
    distance_type sigma_succ;

    pred_succ_search_result_set(std::size_t k, distance_type rad)
        : pred(),
          succ(),
          num_neighbors(k),
          radius(rad),
          sigma_pred(rad),
          sigma_succ(rad) {}

    void register_vantage_point(const point_type& query_point,
                                const point_type& current_vp_point,
                                distance_type current_dist,
                                vertex_type current_vp,
                                const parting_metrics_type& dual_distance) {

      // using the assumption that (current_dist <= current_pred_dist)  -->  if (current_pred_dist < sigma_pred) then
      // (current_dist < sigma_pred)
      if (current_dist < sigma_pred) {
        distance_type current_pred_dist =
            dual_distance.distance(current_vp_point, query_point);
        // is the vantage point within current search bound?
        if (current_pred_dist < sigma_pred) {
          // then add the vantage point to the NN list.
          pred.emplace_back(current_pred_dist, current_vp);
          std::push_heap(pred.begin(), pred.end(), priority_compare_type());
          // are there too many nearest neighbors?
          if (pred.size() > num_neighbors) {
            std::pop_heap(pred.begin(), pred.end(), priority_compare_type());
            // delete last element to keep aList with num_neighbors elements
            pred.pop_back();
            // distance of the last element is now the search bound sigma_pred.
            sigma_pred = pred.front().first;
          }
        }
      }

      // using the assumption that (current_dist <= current_succ_dist)  -->  if (current_succ_dist < sigma_succ) then
      // (current_dist < sigma_succ)
      if (current_dist < sigma_succ) {
        distance_type current_succ_dist =
            dual_distance.distance(query_point, current_vp_point);
        // is the vantage point within current search bound?
        if (current_succ_dist < sigma_succ) {
          // then add the vantage point to the NN list.
          succ.emplace_back(current_succ_dist, current_vp);
          std::push_heap(succ.begin(), succ.end(), priority_compare_type());
          // are there too many nearest neighbors?
          if (succ.size() > num_neighbors) {
            std::pop_heap(succ.begin(), succ.end(), priority_compare_type());
            // delete last element to keep aList with num_neighbors elements
            succ.pop_back();
            // distance of the last element is now the search bound sigma_succ.
            sigma_succ = succ.front().first;
          }
        }
      }
      // radius must remain the maximum of both sigma_pred and sigma_succ (most inclusive search).
      radius = ((sigma_succ > sigma_pred) ? sigma_succ : sigma_pred);
    }
  };

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This is the main nearest-neighbor query function. This takes a query point, a maximum
   * neighborhood radius, node to start recursing from, the current max-heap of neighbors,
   * and the maximum number of neighbors. This function can be used for any kind of NN query (single, kNN, or ranged).
   */
  template <typename SearchResultSet>
  void find_nearest_impl(const point_type& pt, SearchResultSet& result) const {

    std::stack<std::pair<vertex_type, distance_type>> tasks;
    tasks.emplace(root_, 0.0);

    while (!tasks.empty()) {
      auto [cur_node, cur_dist] = tasks.top();
      tasks.pop();

      if (cur_dist > result.radius) {
        continue;
      }

      const point_type& current_vp = get(position_, get(vprop_, cur_node));
      distance_type current_dist = distance_.proper_distance(pt, current_vp);

      result.register_vantage_point(pt, current_vp, current_dist, cur_node,
                                    distance_);

      // first, locate the partition in which pt is:
      if (out_degree(cur_node, *tree_) == 0) {
        continue;
      }
      auto e_rg = out_edges(cur_node, *tree_);
      auto ei = e_rg.begin();
      for (; ei != e_rg.end(); ++ei) {
        if (current_dist <= get(mu_, get_property(*tree_, *ei))) {
          break;
        }
      }
      if (ei == e_rg.end()) {
        --ei;  // back-track if the end was reached.
      }

      std::stack<std::pair<vertex_type, distance_type>> temp_invtasks;
      // search in the most likely node.
      temp_invtasks.emplace(target(*ei, *tree_), 0.0);

      auto ei_left = ei;
      auto ei_right = ei;
      ++ei_right;
      // find the bounds again (start and end).
      bool left_stopped = (ei_left == e_rg.begin());
      bool right_stopped = (ei_right == e_rg.end());
      while (true) {
        if (left_stopped) {
          auto ei_rightleft = ei_right;
          --ei_rightleft;
          distance_type temp_dist = 0.0;
          while ((ei_right != e_rg.end()) &&
                 ((temp_dist = get(mu_, get_property(*tree_, *ei_rightleft)) -
                               current_dist) < result.radius)) {
            temp_invtasks.emplace(target(*ei_right, *tree_), temp_dist);
            ++ei_rightleft;
            ++ei_right;
          }
          break;
        }
        if (right_stopped) {
          auto ei_leftleft = ei_left;
          distance_type temp_dist = 0.0;
          while (
              (ei_left != e_rg.begin()) &&
              ((temp_dist = current_dist -
                            get(mu_, get_property(*tree_, *(--ei_leftleft)))) <
               result.radius)) {
            temp_invtasks.emplace(target(*ei_leftleft, *tree_), temp_dist);
            --ei_left;
          }
          break;
        }
        auto ei_leftleft = ei_left;
        --ei_leftleft;
        // greater than 0 if ei_leftleft should be searched.
        distance_type d1 = get(mu_, get_property(*tree_, *ei_leftleft));
        auto ei_rightleft = ei_right;
        --ei_rightleft;
        // less than 0 if ei_right should be searched.
        distance_type d2 = get(mu_, get_property(*tree_, *ei_rightleft));
        if (d1 + d2 > 2.0 * current_dist) {
          // this means that ei_leftleft's boundary is closer to pt.
          if (d1 + result.radius - current_dist > 0) {
            temp_invtasks.emplace(target(*ei_leftleft, *tree_),
                                  current_dist - d1);
            ei_left = ei_leftleft;
            if (d2 - result.radius - current_dist < 0) {
              temp_invtasks.emplace(target(*ei_right, *tree_),
                                    d2 - current_dist);
              ++ei_right;
            } else {
              right_stopped = true;
            }
          } else {
            break;
          }
        } else {
          if (d2 - result.radius - current_dist < 0) {
            temp_invtasks.emplace(target(*ei_right, *tree_), d2 - current_dist);
            ++ei_right;
            if (d1 + result.radius - current_dist > 0) {
              temp_invtasks.emplace(target(*ei_leftleft, *tree_),
                                    current_dist - d1);
              ei_left = ei_leftleft;
            } else {
              left_stopped = true;
            }
          } else {
            break;
          }
        }
        left_stopped = (ei_left == e_rg.begin());
        right_stopped = (ei_right == e_rg.end());
      }

      // reverse the temporary stack into the main stack.
      while (!temp_invtasks.empty()) {
        tasks.emplace(std::move(temp_invtasks.top()));
        temp_invtasks.pop();
      }
    }
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This function does a single nearest-neighbor query to find the tree-leaf that is closest to
   * the given point (starts to recurse from node). */
  vertex_type get_leaf(const point_type& pt, vertex_type node) const {
    while (out_degree(node, *tree_) != 0) {
      // first, locate the partition in which pt is:
      distance_type current_dist =
          distance_.proper_distance(pt, get(position_, get(vprop_, node)));
      vertex_type result = node;
      for (auto e : out_edges(node, *tree_)) {
        result = target(e, *tree_);
        if (current_dist <= get(mu_, get_property(*tree_, e))) {
          break;
        }
      }
      node = result;
    }
    return node;
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This function finds the tree-vertex with the given key-value and position. */
  vertex_type get_vertex(key_type key, const point_type& pt,
                         vertex_type node) const {
    vertex_type other_branch = bagl::graph_traits<tree_indexer>::null_vertex();
    while (get(key_, get(vprop_, node)) != key) {
      distance_type current_dist =
          distance_.proper_distance(pt, get(position_, get(vprop_, node)));
      if (out_degree(node, *tree_) == 0) {
        if (other_branch != bagl::graph_traits<tree_indexer>::null_vertex()) {
          node = other_branch;
          other_branch = bagl::graph_traits<tree_indexer>::null_vertex();
          continue;
        }
        throw int(0);
      }
      vertex_type result = node;
      auto e_rg = out_edges(node, *tree_);
      for (auto ei = e_rg.begin(); ei != e_rg.end(); ++ei) {
        result = target(*ei, *tree_);
        if (current_dist < get(mu_, get_property(*tree_, *ei))) {
          break;
        }
        if (current_dist == get(mu_, get_property(*tree_, *ei))) {
          ++ei;
          if (ei != e_rg.end()) {
            other_branch = target(*ei, *tree_);
          }
          break;
        }
      }
      node = result;
    }
    return node;
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This function updates the edge distance values as a consequence of a new point being added. */
  void update_mu_upwards(const point_type& pt, vertex_type node) {
    while (node != root_) {
      vertex_type parent = source(*in_edges(node, *tree_).begin(), *tree_);
      distance_type dist =
          distance_.proper_distance(pt, get(position_, get(vprop_, parent)));
      if (dist >
          get(mu_, get_property(*tree_, *in_edges(node, *tree_).begin()))) {
        put(mu_, get_property(*tree_, *in_edges(node, *tree_).begin()), dist);
      }
      node = parent;
    }
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* This function determines if a given node has no children or if all its children have no children. */
  bool is_leaf_node(vertex_type node) const {
    if (out_degree(node, *tree_) == 0) {
      return true;
    }
    for (auto e : out_edges(node, *tree_)) {
      if (out_degree(target(e, *tree_), *tree_) != 0) {
        return false;
      }
    }
    return true;
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This function determines if a given node is the root of a balanced (and full) sub-tree. */
  bool is_node_full(vertex_type node, int& depth_limit) const {
    if (depth_limit < 0) {
      return false;
    }
    std::queue<std::pair<vertex_type, int>> tasks;
    tasks.emplace(node, depth_limit);
    while (!tasks.empty()) {
      auto [cur_node, cur_limit] = tasks.front();
      tasks.pop();
      depth_limit = std::min(cur_limit, depth_limit);

      if (out_degree(cur_node, *tree_) == 0 && cur_limit == 0) {
        continue;
      }

      --cur_limit;

      if (((out_degree(cur_node, *tree_) != 0) && (cur_limit < 0)) ||
          (out_degree(cur_node, *tree_) < Arity) ||
          ((cur_limit > 0) && (is_leaf_node(cur_node)))) {
        depth_limit = cur_limit;
        return false;
      }

      for (auto e : out_edges(cur_node, *tree_)) {
        tasks.emplace(target(e, *tree_), cur_limit);
      }
    }
    return (depth_limit == 0);
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* NOTE This is a non-recursive version. */
  /* This function collects the list of vertices that are in the sub-tree rooted at the given node (exclusively). */
  void collect_vertices(std::vector<vertex_type>& vlist,
                        vertex_type node) const {
    std::queue<vertex_type> tasks;
    tasks.push(node);
    while (!tasks.empty()) {
      vertex_type current_node = tasks.front();
      tasks.pop();
      for (auto e : out_edges(current_node, *tree_)) {
        vlist.push_back(target(e, *tree_));
        tasks.push(target(e, *tree_));
      }
    }
  }

  /* Does not invalidate vertices */
  /* Does not require persistent vertices */
  /* This function computes the maximum depth of the tree rooted at the given node (note: this is an
   * expensive operation, but is useful when debugging and testing). */
  /* NOTE: This is a recursive function, but it isn't practical anyways, and kind of annoying to de-recursify. */
  std::size_t get_depth(vertex_type node) const {
    std::size_t max_depth = 0;
    for (auto e : out_edges(node, *tree_)) {
      std::size_t temp = get_depth(target(e, *tree_));
      if (temp > max_depth) {
        max_depth = temp;
      }
    }
    return max_depth + 1;
  }

 public:
  /**
   * Construct the DVP-tree from a graph, position-map, topology, etc..
   * \tparam Graph The graph type on which the vertices are taken from, should model the bagl::concepts::VertexListGraph.
   * \tparam GraphPositionMap The property-map that associates position values to nodes of the graph.
   * \param g The graph from which to take the vertices.
   * \param g_position The property-map that takes a node of the graph and produces (or looks up) a position value.
   * \param tree The tree object that will be used to store the DVP structure.
   * \param space The topology on which the positions of the vertices reside.
   * \param key The key-map to use to obtain and store the key values for a given vertex-property object.
   * \param mu The property-map which associates a distance-value to each edge of the tree.
   * \param pos The property-map that can be used to obtain and store the positions of the vertices
   * (vertex-property objects).
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  template <typename Graph, typename GraphPositionMap>
  dvp_tree_impl(const Graph& g, GraphPositionMap g_position, tree_indexer& tree,
                const std::shared_ptr<const Space>& space, VertexPropMap vprop,
                VertexKeyMap key, DistanceMap mu, PositionMap pos,
                VPChooser vp_chooser)
      : tree_(&tree),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        vprop_(vprop),
        key_(key),
        mu_(mu),
        position_(pos),
        distance_(space),
        vp_chooser_(vp_chooser) {

    if (num_vertices(g) == 0) {
      return;
    }

    std::vector<vertex_property> v_bin;
    v_bin.reserve(num_vertices(g));
    for (auto v : vertices(g)) {
      vertex_property vp;
      put(key_, vp, v);
      put(position_, vp, get(g_position, v));
      v_bin.emplace_back(std::move(vp));
    }

    auto v_first = v_bin.begin();
    auto v_last = v_bin.end();
    rearrange_with_chosen_vp(v_first, v_last);
    root_ = create_root(*tree_, std::move(*v_first));
    construct_node(root_, ++v_first, v_last);
  }

  // non-copyable.
  dvp_tree_impl(const self&) = delete;
  self& operator=(const self&) = delete;

  /**
   * Construct the DVP-tree from a range, position-map, topology, etc..
   * \tparam ForwardIterator The forward-iterator type from which the vertices can be obtained.
   * \tparam ElemPositionMap The property-map that associates position values to vertices in the given range.
   * \param first The start of the range from which to take the vertices.
   * \param last The end of the range from which to take the vertices (one-past-last).
   * \param elem_position The property-map that takes a node in the given range and produces (or looks up) a position
   * value.
   * \param tree The tree object that will be used to store the DVP structure.
   * \param space The topology on which the positions of the vertices reside.
   * \param key The key-map to use to obtain and store the key values for a given vertex-property object.
   * \param mu The property-map which associates a distance-value to each edge of the tree.
   * \param pos The property-map that can be used to obtain and store the positions of the vertices
   * (vertex-property objects).
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  template <typename ForwardIterator, typename ElemPositionMap>
  dvp_tree_impl(ForwardIterator first, ForwardIterator last,
                ElemPositionMap elem_position, tree_indexer& tree,
                const std::shared_ptr<const Space>& space, VertexPropMap vprop,
                VertexKeyMap key, DistanceMap mu, PositionMap pos,
                VPChooser vp_chooser)
      : tree_(&tree),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        vprop_(vprop),
        key_(key),
        mu_(mu),
        position_(pos),
        distance_(space),
        vp_chooser_(vp_chooser) {
    if (first == last) {
      return;
    }

    // Copy the list of vertices to random access memory.
    std::vector<vertex_property> v_bin;
    for (; first != last; ++first) {
      vertex_property vp;
      put(key_, vp, *first);
      put(position_, vp, get(elem_position, *first));
      v_bin.emplace_back(std::move(vp));
    }

    prop_vector_iter v_first = v_bin.begin();
    prop_vector_iter v_last = v_bin.end();
    rearrange_with_chosen_vp(v_first, v_last);
    root_ = create_root(*tree_, std::move(*v_first));
    construct_node(root_, ++v_first, v_last);
  }

  /**
   * Construct an empty DVP-tree from a topology, etc..
   * \param tree The tree object that will be used to store the DVP structure.
   * \param space The topology on which the positions of the vertices reside.
   * \param key The key-map to use to obtain and store the key values for a given vertex-property object.
   * \param mu The property-map which associates a distance-value to each edge of the tree.
   * \param pos The property-map that can be used to obtain and store the positions of the vertices
   * (vertex-property objects).
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  dvp_tree_impl(tree_indexer& tree, const std::shared_ptr<const Space>& space,
                VertexPropMap vprop, VertexKeyMap key, DistanceMap mu,
                PositionMap pos, VPChooser vp_chooser)
      : tree_(&tree),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        vprop_(vprop),
        key_(key),
        mu_(mu),
        position_(pos),
        distance_(space),
        vp_chooser_(vp_chooser) {}

  // sort-of copyable:
  dvp_tree_impl(tree_indexer& tree, VertexPropMap vprop, const self& rhs)
      : tree_(&tree),
        root_(tree_root(tree)),
        vprop_(vprop),
        key_(rhs.key_),
        mu_(rhs.mu_),
        position_(rhs.position_),
        distance_(rhs.distance_),
        vp_chooser_(rhs.vp_chooser_) {}
  void reassign_copied(tree_indexer& tree, VertexPropMap vprop,
                       const self& rhs) noexcept {
    tree_ = &tree;
    root_ = tree_root(tree);
    vprop_ = vprop;
    key_ = rhs.key_;
    mu_ = rhs.mu_;
    position_ = rhs.position_;
    distance_ = rhs.distance_;
    vp_chooser_ = rhs.vp_chooser_;
  }

  // sort-of movable:
  dvp_tree_impl(tree_indexer& tree, VertexPropMap vprop, self&& rhs) noexcept
      : tree_(&tree),
        root_(tree_root(tree)),
        vprop_(vprop),
        key_(std::move(rhs.key_)),
        mu_(std::move(rhs.mu_)),
        position_(std::move(rhs.position_)),
        distance_(std::move(rhs.distance_)),
        vp_chooser_(std::move(rhs.vp_chooser_)) {
    rhs.tree_ = nullptr;
    rhs.root_ = bagl::graph_traits<tree_indexer>::null_vertex();
  }
  void reassign_moved(tree_indexer& tree, VertexPropMap vprop,
                      self&& rhs) noexcept {
    tree_ = &tree;
    root_ = tree_root(tree);
    rhs.tree_ = nullptr;
    rhs.root_ = bagl::graph_traits<tree_indexer>::null_vertex();
    vprop_ = vprop;
    key_ = std::move(rhs.key_);
    mu_ = std::move(rhs.mu_);
    position_ = std::move(rhs.position_);
    distance_ = std::move(rhs.distance_);
    vp_chooser_ = std::move(rhs.vp_chooser_);
  }

  /**
   * Checks if the DVP-tree is empty.
   * \return True if the DVP-tree is empty.
   */
  bool empty() const { return (num_vertices(*tree_) == 0); }
  /**
   * Returns the size of the DVP-tree (the number of vertices it contains.
   * \return The size of the DVP-tree (the number of vertices it contains.
   */
  std::size_t size() const { return num_vertices(*tree_); }

  /**
   * Returns the depth of the tree.
   * \note This operation must recurse through all the branches of the tree (depth-first), and is
   * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack)
   * w.r.t. the depth of tree).
   * \return The depth of the tree.
   */
  std::size_t depth() const { return get_depth(root_); }

  /**
   * This function computes an approximation of the characteristic size of the vertices in the DVP tree.
   * \return The approximation of the characteristic size of the vertices in the DVP tree.
   */
  double get_characteristic_size() const {
    if (num_vertices(*tree_) == 0) {
      return std::numeric_limits<double>::infinity();
    }

    double max_dist = 0.0;
    for (auto e : out_edges(root_, *tree_)) {
      double cur_dist = get(mu_, get_property(*tree_, e));
      if (cur_dist > max_dist) {
        max_dist = cur_dist;
      }
    }
    if (max_dist != 0.0) {
      return max_dist;
    }  // never found a finite, non-zero distance value.
    return std::numeric_limits<double>::infinity();
  }

  /**
   * Inserts a vertex into the tree.
   * \param up The vertex-property to be added to the DVP-tree.
   */
  void insert(vertex_property up) {
    if (num_vertices(*tree_) == 0) {
      root_ = create_root(*tree_, std::move(up));
      return;
    }

    point_type u_pt = get(position_, up);
    // Store the root of subtree to reconstruct.
    vertex_type u_subroot = get_leaf(u_pt, root_);
    // NOTE: if the root is the leaf, it requires special attention since no parent exists.
    if (u_subroot != root_) {
      vertex_type u_leaf = source(*in_edges(u_subroot, *tree_).begin(), *tree_);
      if ((out_degree(u_leaf, *tree_) == Arity) && (is_leaf_node(u_leaf))) {
        // if u_leaf is a full-leaf, then it is balanced but full,
        // we should then find a non-full parent.
        int actual_depth_limit = 1;
        int last_depth_limit = actual_depth_limit;
        while ((u_leaf != root_) && (is_node_full(u_leaf, last_depth_limit))) {
          u_leaf = source(*in_edges(u_leaf, *tree_).begin(), *tree_);
          last_depth_limit = ++actual_depth_limit;
        }
        bool is_p_full = false;
        if (u_leaf == root_) {
          is_p_full = is_node_full(u_leaf, last_depth_limit);
        }
        if ((!is_p_full) && (last_depth_limit >= 0)) {
          // this means that we can add our key to the sub-tree of u_leaf and reconstruct from there.
          u_subroot = u_leaf;
        }
        // else:
        //  this means that either the root node is full or there are
        //  branches of the tree that are deeper than u_subroot,
        //  and thus, in either case, u_subroot should be expanded.
      } else {
        //  leaf node is not full of children, an additional child can be added
        //  (must be reconstructed to keep ordering, but this is a trivial operation O(Arity)).
        // OR
        //  if leaf is not really a leaf, then it means that this sub-tree is definitely
        //  not balanced and not full either,
        //  then all the Keys ought to be collected and u_leaf ought to be reconstructed.
        u_subroot = u_leaf;
      }
    }

    update_mu_upwards(u_pt, u_subroot);
    std::vector<vertex_property> prop_list;
    prop_list.emplace_back(std::move(up));
    while (out_degree(u_subroot, *tree_) > 0) {
      edge_type e = *out_edges(u_subroot, *tree_).begin();
      remove_branch(target(e, *tree_), *tree_, std::back_inserter(prop_list));
    }
    construct_node(u_subroot, prop_list.begin(), prop_list.end());
  }
  /**
   * Inserts a range of vertices.
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
   * \param first The start of the range from which to take the vertices.
   * \param last The end of the range from which to take the vertices (one-past-last).
   */
  template <typename ForwardIterator>
  void insert(ForwardIterator first, ForwardIterator last) {
    for (; first != last; ++first) {
      insert(*first);
    }
  }

  /**
   * Erases the given vertex from the DVP-tree.
   * \param u_node The vertex to be removed from the DVP-tree.
   */
  void erase(vertex_type u_node) {
    if (num_vertices(*tree_) == 0) {
      return;
    }
    std::vector<vertex_property> prop_list;
    if ((u_node == root_) && (num_vertices(*tree_) == 1)) {
      remove_branch(root_, *tree_, std::back_inserter(prop_list));
      root_ = bagl::graph_traits<tree_indexer>::null_vertex();
      return;
    }
    vertex_type u_parent = bagl::graph_traits<tree_indexer>::null_vertex();
    if (u_node != root_) {
      u_parent = source(*in_edges(u_node, *tree_).begin(), *tree_);
    }

    // remove-and-collect all children of u_node:
    while (out_degree(u_node, *tree_) > 0) {
      edge_type e = *out_edges(u_node, *tree_).begin();
      remove_branch(target(e, *tree_), *tree_, std::back_inserter(prop_list));
    }
    // remove-and-discard u_node:
    {
      std::vector<vertex_property> throwaway_list;
      remove_branch(u_node, *tree_, std::back_inserter(throwaway_list));
    }

    if (u_parent != bagl::graph_traits<tree_indexer>::null_vertex()) {
      // remove-and-collect all other children of u_parent:
      while (out_degree(u_parent, *tree_) > 0) {
        edge_type e = *out_edges(u_parent, *tree_).begin();
        remove_branch(target(e, *tree_), *tree_, std::back_inserter(prop_list));
      }
      // reconstruct parent's subtree:
      construct_node(u_parent, prop_list.begin(), prop_list.end());
    } else {
      // need to re-construct the root node (u_node == root_).
      auto v_first = prop_list.begin();
      auto v_last = prop_list.end();
      rearrange_with_chosen_vp(v_first, v_last);
      root_ = create_root(*tree_, std::move(*v_first));
      construct_node(root_, ++v_first, v_last);
    }
  }

  /**
   * Looks up the given key-value and position from the DVP-tree.
   * \param u_key The key-value to be removed from the DVP-tree.
   * \param u_pt The position-value corresponding to the key-value.
   * \return The vertex descriptor (into the tree-storage) for the given key value.
   */
  vertex_type get_vertex(key_type u_key, const point_type& u_pt) const {
    try {
      return get_vertex(u_key, u_pt, root_);
    } catch (int err) {
      return bagl::graph_traits<tree_indexer>::null_vertex();
    }
  }

  /**
   * Erases the given key-value and position from the DVP-tree.
   * Note that this function is not recommended if the vertex descriptor is known, as it will require a key-lookup.
   * \param u_key The key-value to be removed from the DVP-tree.
   * \param u_pt The position-value corresponding to the key-value.
   */
  void erase(key_type u_key, const point_type& u_pt) {
    vertex_type u_node;
    try {
      u_node = get_vertex(u_key, u_pt, root_);
    } catch (int /*err*/) {
      return;
    }
    erase(u_node);
  }

  /**
   * Erases the given vertex-range from the DVP-tree.
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices (by tree vertex
   * descriptors).
   * \param first The start of the range from which to take the vertices to be erased.
   * \param last The end of the range from which to take the vertices to be erased (one-past-last).
   */
  template <typename ForwardIterator>
  void erase(ForwardIterator first, ForwardIterator last) {
    if (num_vertices(*tree_) == 0) {
      return;
    }

    using vertex_listing =
        std::list<std::pair<vertex_type, std::vector<vertex_type>>>;
    vertex_listing
        v_lists;  // will hold a list of unique nodes and all their non-erased
    std::unordered_set<key_type, bagl::graph_descriptor_hash_t<key_type>>
        invalid_keys;

    // First, generate the vertex-listings in preparation for the deletion.
    for (; first != last; ++first) {
      vertex_type removal_trunk = *first;
      // mark as invalid, for deletion.
      invalid_keys.insert(get(key_, get(vprop_, *first)));
      // go up until a valid parent node is found:
      while ((removal_trunk != root_) &&
             (invalid_keys.count(get(key_, get(vprop_, removal_trunk))) != 0)) {
        removal_trunk =
            source(*in_edges(removal_trunk, *tree_).begin(), *tree_);
      }

      bool already_collected = false;
      for (const auto& v : v_lists) {
        // is removal_trunk contained in the *it listing?
        auto it_key =
            std::binary_search(v.second.begin(), v.second.end(), removal_trunk);
        if (it_key != v.second.end()) {
          already_collected = true;
          break;
        }
      }

      if (!already_collected) {
        std::vector<vertex_type> removal_list;
        removal_list.push_back(removal_trunk);
        collect_vertices(removal_list, removal_trunk);
        std::sort(removal_list.begin(), removal_list.end());

        for (auto it = v_lists.begin(); it != v_lists.end();) {
          // is the *it trunk contained in the removal_trunk?
          auto it_key = std::binary_search(removal_list.begin(),
                                           removal_list.end(), it->first);
          if (it_key != removal_list.end()) {
            it = v_lists.erase(it);
          } else {
            ++it;
          }
        }

        v_lists.emplace_back(removal_trunk, std::vector<vertex_type>());
        v_lists.back().second.swap(removal_list);
      }
    }

    auto is_vertex_prop_valid = [&](const vertex_property& vp) {
      return (invalid_keys.count(get(key_, vp)) == 0);
    };

    if ((v_lists.size() == 1) && (v_lists.front() == root_)) {
      // need to re-construct the root node (u_node == root_).
      std::vector<vertex_property> prop_list;
      remove_branch(root_, *tree_, std::back_inserter(prop_list));
      prop_list.erase(
          remove_if(prop_list.begin(), prop_list.end(), is_vertex_prop_valid),
          prop_list.end());
      prop_vector_iter v_first = prop_list.begin();
      prop_vector_iter v_last = prop_list.end();
      rearrange_with_chosen_vp(v_first, v_last);
      root_ = create_root(*tree_, std::move(*v_first));
      construct_node(root_, ++v_first, v_last);
      return;
    }

    for (auto& v : v_lists) {
      std::vector<vertex_property> prop_list;
      // remove-and-collect all children of v.first:
      while (out_degree(v.first, *tree_) > 0) {
        edge_type e = *out_edges(v.first, *tree_).begin();
        remove_branch(target(e, *tree_), *tree_, std::back_inserter(prop_list));
      }
      // erase removed vertices from the prop-list:
      prop_list.erase(
          remove_if(prop_list.begin(), prop_list.end(), is_vertex_prop_valid),
          prop_list.end());
      // reconstruct parent's subtree:
      construct_node(v.first, prop_list.begin(), prop_list.end());
    }
  }

  /**
   * Clears the DVP-tree.
   */
  void clear() {
    if (num_vertices(*tree_) == 0) {
      std::vector<vertex_property> prop_list;
      remove_branch(root_, *tree_, std::back_inserter(prop_list));
      root_ = bagl::graph_traits<tree_indexer>::null_vertex();
    }
  }

  /**
   * Finds the nearest neighbor to a given position.
   * \param pt The position from which to find the nearest-neighbor of.
   * \return The vertex in the DVP-tree that is closest to the given point.
   */
  vertex_type find_nearest(const point_type& pt) const {
    if (num_vertices(*tree_) == 0) {
      return bagl::graph_traits<tree_indexer>::null_vertex();
    }
    nearest_search_result_set result_set(
        1, std::numeric_limits<distance_type>::infinity());
    find_nearest_impl(pt, result_set);
    if (result_set.neighbors.size()) {
      return result_set.neighbors.front().second;
    }
    return bagl::graph_traits<tree_indexer>::null_vertex();
  }

  /**
   * Finds the nearest predecessor and successor to a given position.
   * \param pt The position from which to find the nearest predecessor and successor of.
   * \return The predecessor and successor vertex in the DVP-tree that is closest to the given point.
   */
  std::pair<vertex_type, vertex_type> find_nearest_pred_succ(
      const point_type& pt) const {
    if (num_vertices(*tree_) == 0) {
      return bagl::graph_traits<tree_indexer>::null_vertex();
    }
    pred_succ_search_result_set result_set(
        1, std::numeric_limits<distance_type>::infinity());
    find_nearest_impl(pt, result_set);
    std::pair<vertex_type, vertex_type> result;
    if (result_set.pred.size()) {
      result.first = result_set.pred.front().second;
    } else {
      result.first = bagl::graph_traits<tree_indexer>::null_vertex();
    }
    if (result_set.succ.size()) {
      result.second = result_set.succ.front().second;
    } else {
      result.second = bagl::graph_traits<tree_indexer>::null_vertex();
    }
    return result;
  }

  /**
   * Finds the K nearest-neighbors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors by tree vertex descriptors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param out_iter An iterator to the first place where to put the sorted list of
   *        elements by tree vertex descriptors with the smallest distance.
   * \param num_neighbors The number of nearest-neighbors.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_nearest(
      const point_type& pt, OutputIterator out_iter, std::size_t num_neighbors,
      distance_type radius =
          std::numeric_limits<distance_type>::infinity()) const {
    if (num_vertices(*tree_) == 0) {
      return out_iter;
    }
    nearest_search_result_set result_set(num_neighbors, radius);
    find_nearest_impl(pt, result_set);
    std::sort_heap(result_set.neighbors.begin(), result_set.neighbors.end(),
                   priority_compare_type());
    for (auto& v : result_set.neighbors) {
      *(out_iter++) = v.second;
    }
    return out_iter;
  }

  /**
   * Finds the K nearest predecessors and successors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors by tree vertex descriptors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param pred_iter An iterator to the first place where to put the sorted list of
   *        predecessors by tree vertex descriptors with the smallest distance.
   * \param succ_iter An iterator to the first place where to put the sorted list of
   *        successors by tree vertex descriptors with the smallest distance.
   * \param num_neighbors The number of nearest-neighbors.
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterators to the end of the lists of nearest predecessors and successors (starting from
   * "output_first").
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_nearest(
      const point_type& pt, OutputIterator pred_iter, OutputIterator succ_iter,
      std::size_t num_neighbors,
      distance_type radius =
          std::numeric_limits<distance_type>::infinity()) const {
    if (num_vertices(*tree_) == 0) {
      return {pred_iter, succ_iter};
    }
    pred_succ_search_result_set result_set(num_neighbors, radius);
    find_nearest_impl(pt, result_set);
    std::sort_heap(result_set.pred.begin(), result_set.pred.end(),
                   priority_compare_type());
    std::sort_heap(result_set.succ.begin(), result_set.succ.end(),
                   priority_compare_type());
    for (auto [d, v] : result_set.pred) {
      *(pred_iter++) = v;
    }
    for (auto [d, v] : result_set.succ) {
      *(succ_iter++) = v;
    }
    return {pred_iter, succ_iter};
  }

  /**
   * Finds the nearest-neighbors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors by tree vertex descriptors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param out_iter An iterator to the first place where to put the sorted list of
   *        elements by tree vertex descriptors with the smallest distance.
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_in_range(const point_type& pt, OutputIterator out_iter,
                               distance_type radius) const {
    if (num_vertices(*tree_) == 0) {
      return out_iter;
    }
    nearest_search_result_set result_set(num_vertices(*tree_), radius);
    find_nearest_impl(pt, result_set);
    std::sort_heap(result_set.neighbors.begin(), result_set.neighbors.end(),
                   priority_compare_type());
    for (auto [d, v] : result_set.neighbors) {
      *(out_iter++) = v;
    }
    return out_iter;
  };

  /**
   * Finds the nearest predecessors and successors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors by tree vertex descriptors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param pred_iter An iterator to the first place where to put the sorted list of
   *        predecessors by tree vertex descriptors with the smallest distance.
   * \param succ_iter An iterator to the first place where to put the sorted list of
   *        successors by tree vertex descriptors with the smallest distance.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterators to the end of the lists of nearest predecessors and successors (starting from
   * "output_first").
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_in_range(
      const point_type& pt, OutputIterator pred_iter, OutputIterator succ_iter,
      distance_type radius) const {
    if (num_vertices(*tree_) == 0) {
      return {pred_iter, succ_iter};
    }
    pred_succ_search_result_set result_set(num_vertices(*tree_), radius);
    find_nearest_impl(pt, result_set);
    std::sort_heap(result_set.pred.begin(), result_set.pred.end(),
                   priority_compare_type());
    std::sort_heap(result_set.succ.begin(), result_set.succ.end(),
                   priority_compare_type());
    for (auto [d, v] : result_set.pred) {
      *(pred_iter++) = v;
    }
    for (auto [d, v] : result_set.succ) {
      *(succ_iter++) = v;
    }
    return {pred_iter, succ_iter};
  }

  struct mutation_visitor {
    self* parent;

    explicit mutation_visitor(self* p) : parent(p) {}

    void remove_vertex(vertex_type v, tree_indexer& /*unused*/) const {
      parent->erase(v);
    }

    void add_vertex(tree_indexer& /*unused*/, const vertex_property& vp) const {
      parent->insert(vp);
    }

    void add_vertex(tree_indexer& /*unused*/, vertex_property&& vp) const {
      parent->insert(std::move(vp));
    }
  };
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_DVP_TREE_DETAIL_H_
