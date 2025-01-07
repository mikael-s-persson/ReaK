/**
 * \file metric_space_search.h
 *
 * This library provides a class that implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which
 * respects the triangular inequality.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_PLANNING_PATH_PLANNING_METRIC_SPACE_SEARCH_H_
#define REAK_PLANNING_PATH_PLANNING_METRIC_SPACE_SEARCH_H_

#include "bagl/bfl_d_ary_tree.h"  // for default tree storage.
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

#include <limits>   // for numeric_limits
#include <utility>  // for pair and move
#include <vector>   // for vector

#include "ReaK/planning/path_planning/dvp_tree_detail.h"
#include "ReaK/planning/path_planning/multi_dvp_tree_search.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

namespace ReaK::pp {

/**
 * This class is a position-caching policy class for the space partitioning trees that are indirectly
 * indexing a set of vertices of a graph (or other container). When building an indirect index, it can
 * be advantageous to store (copy) the position values (or vector) in the index itself to increase
 * the locality of reference during the traversal-mutating operations (lookups, insertions, etc.)
 * through the space partitioning tree. This position-caching policy mandates that no such caching
 * is done.
 */
struct no_position_caching_policy {
  template <typename PointType>
  struct vertex_base_property {};  // intentionally empty.

  template <typename Graph, typename PointType>
  struct effector {
    explicit effector(Graph& /*unused*/) {}

    template <typename Vertex, typename Key, typename PositionMap,
              typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k,
                    const PositionMap& /*unused*/, const KeyMap& key) const {
      put(key, v, k);
    }

    template <typename Vertex, typename PositionMap, typename KeyMap>
    PointType get_position(const Vertex& v, const PositionMap& position,
                           const KeyMap& key) const {
      return get(position, get(key, v));
    }
  };

  template <typename PointType, typename KeyMap, typename PositionMap>
  struct position_map : bagl::composite_property_map<PositionMap, KeyMap> {
    using base = bagl::composite_property_map<PositionMap, KeyMap>;
    using self = position_map<PointType, KeyMap, PositionMap>;

    explicit position_map(KeyMap k, PositionMap p = PositionMap())
        : base(p, k) {}

    position_map() : position_map(KeyMap()) {}
  };
};

/**
 * This class is a position-caching policy class for the space partitioning trees that are indirectly
 * indexing a set of vertices of a graph (or other container). When building an indirect index, it can
 * be advantageous to store (copy) the position values (or vector) in the index itself to increase
 * the locality of reference during the traversal-mutating operations (lookups, insertions, etc.)
 * through the space partitioning tree. This position-caching policy mandates such caching of the
 * position values within the index.
 */
struct position_caching_policy {
  template <typename PointType>
  struct vertex_base_property {
    PointType cached_position_value;
  };

  template <typename Graph, typename PointType>
  struct effector {

    bagl::property_map_t<Graph, PointType vertex_base_property<PointType>::*>
        cached_position;

    explicit effector(Graph& g)
        : cached_position(get(
              &vertex_base_property<PointType>::cached_position_value, g)) {}

    template <typename Vertex, typename Key, typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PointType& pt,
                    const KeyMap& key) const {
      put(key, v, k);
      put(cached_position, v, pt);
    }

    template <typename Vertex, typename Key, typename PositionMap,
              typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PositionMap& position,
                    const KeyMap& key) const {
      put(key, v, k);
      put(cached_position, v, get(position, k));
    }

    template <typename Vertex, typename PositionMap, typename KeyMap>
    PointType get_position(const Vertex& v, const PositionMap& /*unused*/,
                           const KeyMap& /*unused*/) const {
      return get(cached_position, v);
    }
  };

  template <typename PointType, typename KeyMap, typename PositionMap>
  struct position_map
      : bagl::data_member_property_map<PointType,
                                       vertex_base_property<PointType>> {
    using base =
        bagl::data_member_property_map<PointType,
                                       vertex_base_property<PointType>>;
    using self = position_map<PointType, KeyMap, PositionMap>;

    explicit position_map(KeyMap k, PositionMap p = PositionMap())
        : base(&vertex_base_property<PointType>::cached_position_value) {}

    position_map() : position_map(KeyMap()) {}
  };
};

/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which
 * respects the triangular inequality.
 * \tparam Key The key type for the tree, essentially the key value is the vertex descriptor type.
 * \tparam Space The topology type on which the points can reside, should model the MetricSpaceConcept.
 * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the
 * iterators) to a point (position).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename Key, MetricSpace Space, typename PositionMap,
          unsigned int Arity = 2, typename VPChooser = random_vp_chooser,
          typename TreeStorageTag = bagl::bfl_d_ary_tree_storage<Arity>,
          typename PositionCachingPolicy = position_caching_policy>
class dvp_tree {
 public:
  using self = dvp_tree<Key, Space, PositionMap, Arity, VPChooser,
                        TreeStorageTag, PositionCachingPolicy>;

  using point_type = bagl::property_traits_value_t<PositionMap>;
  using distance_type = double;

 private:
  struct vertex_properties
      : PositionCachingPolicy::template vertex_base_property<point_type> {
    Key k;
  };

  struct edge_properties {
    distance_type d;
  };

  using tree_indexer =
      typename bagl::tree_storage<vertex_properties, edge_properties,
                                  TreeStorageTag>::type;
  using vp_map_type = bagl::property_map_t<tree_indexer, bagl::vertex_all_t>;
  using key_map_type = bagl::data_member_property_map<Key, vertex_properties>;
  using distance_map_type =
      bagl::data_member_property_map<distance_type, edge_properties>;
  using vertex_position_map =
      typename PositionCachingPolicy::template position_map<
          point_type, key_map_type, PositionMap>;

  using dvp_impl_type =
      dvp_tree_impl<tree_indexer, Space, vp_map_type, key_map_type,
                    distance_map_type, vertex_position_map, Arity, VPChooser>;

  tree_indexer tree_;
  vp_map_type vp_map_;
  PositionMap position_;
  key_map_type vp_key_;
  vertex_position_map vp_pos_;
  dvp_impl_type impl_;

 public:
  /**
   * Construct the DVP-tree from a graph, topology and property-map.
   * \tparam Graph The graph type on which the vertices are taken from.
   * \param g The graph from which to take the vertices.
   * \param space The topology on which the positions of the vertices reside.
   * \param position The property-map that can be used to obtain the positions of the vertices.
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  template <typename Graph>
  dvp_tree(const Graph& g, const std::shared_ptr<const Space>& space,
           PositionMap position, VPChooser vp_chooser = VPChooser())
      : tree_(),
        vp_map_(get(bagl::vertex_all, tree_)),
        position_(position),
        vp_key_(&vertex_properties::k),
        vp_pos_(vp_key_, position_),
        impl_(g, position_, tree_, space, vp_map_, vp_key_,
              distance_map_type(&edge_properties::d), vp_pos_, vp_chooser) {}

  /**
   * Construct the DVP-tree from a range, topology and property-map.
   * \tparam ForwardIterator The forward-iterator type from which the vertices can be obtained.
   * \param first The start of the range from which to take the vertices.
   * \param last The end of the range from which to take the vertices (one-past-last).
   * \param space The topology on which the positions of the vertices reside.
   * \param position The property-map that can be used to obtain the positions of the vertices.
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  template <typename ForwardIterator>
  dvp_tree(ForwardIterator first, ForwardIterator last,
           const std::shared_ptr<const Space>& space, PositionMap position,
           VPChooser vp_chooser = VPChooser())
      : tree_(),
        vp_map_(get(bagl::vertex_all, tree_)),
        position_(position),
        vp_key_(&vertex_properties::k),
        vp_pos_(vp_key_, position_),
        impl_(first, last, position_, tree_, space, vp_map_, vp_key_,
              distance_map_type(&edge_properties::d), vp_pos_, vp_chooser) {}

  dvp_tree(const self& rhs)
      : tree_(rhs.tree_),
        vp_map_(get(bagl::vertex_all, tree_)),
        position_(rhs.position_),
        vp_key_(rhs.vp_key_),
        vp_pos_(rhs.vp_pos_),
        impl_(tree_, vp_map_, rhs.impl_) {}

  self& operator=(const self& rhs) {
    tree_ = rhs.tree_;
    vp_map_ = get(bagl::vertex_all, tree_);
    position_ = rhs.position_;
    vp_key_ = rhs.vp_key_;
    vp_pos_ = rhs.vp_pos_;
    impl_.reassign_copied(tree_, vp_map_, rhs.impl_);
    return *this;
  }

  dvp_tree(self&& rhs) noexcept
      : tree_(std::move(rhs.tree_)),
        vp_map_(get(bagl::vertex_all, tree_)),
        position_(std::move(rhs.position_)),
        vp_key_(std::move(rhs.vp_key_)),
        vp_pos_(std::move(rhs.vp_pos_)),
        impl_(tree_, vp_map_, std::move(rhs.impl_)) {}

  self& operator=(self&& rhs) noexcept {
    tree_ = std::move(rhs.tree_);
    vp_map_ = get(bagl::vertex_all, tree_);
    position_ = std::move(rhs.position_);
    vp_key_ = std::move(rhs.vp_key_);
    vp_pos_ = std::move(rhs.vp_pos_);
    impl_.reassign_moved(tree_, vp_map_, std::move(rhs.impl_));
    return *this;
  }

  /**
   * Checks if the DVP-tree is empty.
   * \return True if the DVP-tree is empty.
   */
  bool empty() const { return impl_.empty(); }

  /**
   * Returns the size of the DVP-tree (the number of vertices it contains.
   * \return The size of the DVP-tree (the number of vertices it contains.
   */
  std::size_t size() const { return impl_.size(); }

  /**
   * Returns the depth of the tree.
   * \note This operation must recurse through all the branches of the tree (depth-first), and is
   * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack)
   * w.r.t. the depth of tree).
   * \return The depth of the tree.
   */
  std::size_t depth() const { return impl_.depth(); }

  /**
   * This function computes an approximation of the characteristic size of the vertices in the DVP tree.
   * \return The approximation of the characteristic size of the vertices in the DVP tree.
   */
  double get_characteristic_size() const {
    return impl_.get_characteristic_size();
  }

  /**
   * Inserts a key-value (vertex).
   * \param u The vertex to be added to the DVP-tree.
   */
  void insert(Key u) {
    vertex_properties vp;
    put(vp_key_, vp, u);
    put(vp_pos_, vp, get(position_, u));
    impl_.insert(std::move(vp));
  }

  /**
   * Inserts a range of key-values (vertices).
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
   * \param first The start of the range from which to take the vertices.
   * \param last The end of the range from which to take the vertices (one-past-last).
   */
  template <typename ForwardIterator>
  void insert(ForwardIterator first, ForwardIterator last) {
    std::vector<vertex_properties> vp_list;
    for (; first != last; ++first) {
      vp_list.emplace_back();
      put(vp_key_, vp_list.back(), *first);
      put(vp_pos_, vp_list.back(), get(position_, *first));
    }
    impl_.insert(vp_list.begin(), vp_list.end());
  }

  /**
   * Erases the given vertex from the DVP-tree.
   * \param u The vertex to be removed from the DVP-tree.
   */
  void erase(Key u) { impl_.erase(u, get(position_, u)); }

  /**
   * Erases the given vertex-range from the DVP-tree.
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
   * \param first The start of the range from which to take the vertices to be erased.
   * \param last The end of the range from which to take the vertices to be erased (one-past-last).
   */
  template <typename ForwardIterator>
  void erase(ForwardIterator first, ForwardIterator last) {
    using TreeVertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    for (; first != last; ++first) {
      TreeVertex u = impl_.get_vertex(*first, get(position_, *first));
      if (u != bagl::graph_traits<tree_indexer>::null_vertex()) {
        v_list.push_back(u);
      }
    }
    impl_.erase(v_list.begin(), v_list.end());
  }

  /**
   * Clears the DVP-tree.
   */
  void clear() { impl_.clear(); }

  /**
   * Finds the nearest neighbor to a given position.
   * \param pt The position from which to find the nearest-neighbor of.
   * \return The vertex in the DVP-tree that is closest to the given point.
   */
  Key find_nearest(const point_type& pt) const {
    auto u = impl_.find_nearest(pt);
    if (u != bagl::graph_traits<tree_indexer>::null_vertex()) {
      return tree_[u].k;
    }
    return Key();
  }

  /**
   * Finds the nearest predecessor and successor to a given position.
   * \param pt The position from which to find the nearest-neighbor of.
   * \return The vertices in the DVP-tree that are nearest predecessor and successor to the given point.
   */
  std::pair<Key, Key> find_nearest_pred_succ(const point_type& pt) const {
    auto [u_pred, u_succ] = impl_.find_nearest_pred_succ(pt);
    std::pair<Key, Key> result{};
    if (u_pred != bagl::graph_traits<tree_indexer>::null_vertex()) {
      result.first = tree_[u_pred].k;
    }
    if (u_succ != bagl::graph_traits<tree_indexer>::null_vertex()) {
      result.second = tree_[u_succ].k;
    }
    return result;
  }

  /**
   * Finds the K nearest-neighbors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param out_iter An iterator to the first place where to put the sorted list of
   *        elements with the smallest distance.
   * \param num_neighbors The number of nearest-neighbors.
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_nearest(
      const point_type& pt, OutputIterator out_iter, std::size_t num_neighbors,
      distance_type radius =
          std::numeric_limits<distance_type>::infinity()) const {
    using TreeVertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    impl_.find_nearest(pt, std::back_inserter(v_list), num_neighbors, radius);
    for (TreeVertex v : v_list) {
      *(out_iter++) = tree_[v].k;
    }
    return out_iter;
  };

  /**
   * Finds the K nearest predecessors and successors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param pred_iter An iterator to the first place where to put the sorted list of
   *        predecessor elements with the smallest distance.
   * \param succ_iter An iterator to the first place where to put the sorted list of
   *        successor elements with the smallest distance.
   * \param num_neighbors The number of nearest-neighbors.
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the two lists of nearest neighbors (predecessors and successors).
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_nearest(
      const point_type& pt, OutputIterator pred_iter, OutputIterator succ_iter,
      std::size_t num_neighbors,
      distance_type radius =
          std::numeric_limits<distance_type>::infinity()) const {
    using TreeVertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
    std::vector<TreeVertex> pred_list;
    std::vector<TreeVertex> succ_list;
    impl_.find_nearest(pt, std::back_inserter(pred_list),
                       std::back_inserter(succ_list), num_neighbors, radius);
    for (TreeVertex v : pred_list) {
      *(pred_iter++) = tree_[v].k;
    }
    for (TreeVertex v : succ_list) {
      *(succ_iter++) = tree_[v].k;
    }
    return {pred_iter, succ_iter};
  }

  /**
   * Finds the nearest-neighbors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param out_iter An iterator to the first place where to put the sorted list of
   *        elements with the smallest distance.
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_in_range(const point_type& pt, OutputIterator out_iter,
                               distance_type radius) const {
    using TreeVertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    impl_.find_in_range(pt, std::back_inserter(v_list), radius);
    for (TreeVertex v : v_list) {
      *(out_iter++) = tree_[v].k;
    }
    return out_iter;
  }

  /**
   * Finds the K nearest predecessors and successors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param pt The position from which to find the nearest-neighbors.
   * \param pred_iter An iterator to the first place where to put the sorted list of
   *        predecessor elements with the smallest distance.
   * \param succ_iter An iterator to the first place where to put the sorted list of
   *        successor elements with the smallest distance.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the two lists of nearest neighbors (predecessors and successors).
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_in_range(
      const point_type& pt, OutputIterator pred_iter, OutputIterator succ_iter,
      distance_type R) const {
    using TreeVertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
    std::vector<TreeVertex> pred_list;
    std::vector<TreeVertex> succ_list;
    impl_.find_in_range(pt, std::back_inserter(pred_list),
                        std::back_inserter(succ_list), R);
    for (TreeVertex v : pred_list) {
      *(pred_iter++) = tree_[v].k;
    }
    for (TreeVertex v : succ_list) {
      *(succ_iter++) = tree_[v].k;
    }
    return {pred_iter, succ_iter};
  }
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_METRIC_SPACE_SEARCH_H_
