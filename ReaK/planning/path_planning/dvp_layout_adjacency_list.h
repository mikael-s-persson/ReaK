/**
 * \file dvp_layout_adjacency_list.h
 *
 * This library provides a class that implements a Dynamic Vantage-Point Tree (DVP-Tree) which
 * is synchronized with an adjacency-list graph and uses the tree-storage as the layout for the
 * vertices common to both graphs (adj-list and tree). DVP-trees
 * allow for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which
 * respects the triangular inequality.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_PLANNING_PATH_PLANNING_DVP_LAYOUT_ADJACENCY_LIST_H_
#define REAK_PLANNING_PATH_PLANNING_DVP_LAYOUT_ADJACENCY_LIST_H_

#include "ReaK/planning/graph_alg/adj_list_tree_overlay.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

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

namespace ReaK::pp {

/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) which
 * is synchronized with a adjacency-list graph and uses the tree-storage as the layout for the
 * vertices common two both graphs (adj-list and tree). The main advantage of this scheme is that
 * because the DVP-tree tends to group vertices that are closer together w.r.t. some distance metric
 * into memory locations that are also close to each other. This means that operations done on the
 * adjacency-list graph on vertices that are neighbouring each other will have better locality of references.
 * DVP-trees allow for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which
 * respects the triangular inequality.
 * \tparam VertexProperty The bundled vertex properties for the adjacency-list.
 * \tparam EdgeProperty The bundled edge properties for the adjacency-list.
 * \tparam Space The topology type on which the points can reside.
 * \tparam PositionMap The property-map type that can map the vertex property to a position in the topology.
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 * \tparam TreeStorageTag A tree-storage tag which specifies the kind of tree structure to use for the DVP tree.
 * \tparam OutEdgeListS The out-edge list container specifier for the adjacency-list (same as OutEdgeListS in
 * bagl::adjacency_list).
 * \tparam DirectedS The edge's directional specifier for the adjacency-list (same as DirectedS in
 * bagl::adjacency_list).
 * \tparam EdgeListS The edge list container specifier for the adjacency-list (same as EdgeListS in
 * bagl::adjacency_list).
 */
template <typename VertexProperty, typename EdgeProperty, MetricSpace Space,
          typename PositionMap, unsigned int Arity = 2,
          typename VPChooser = random_vp_chooser,
          typename TreeStorageTag = bagl::bfl_d_ary_tree_storage<Arity>,
          typename OutEdgeListS = bagl::vec_s,
          typename DirectedS = bagl::directed_s,
          typename EdgeListS = bagl::vec_s>
class dvp_adjacency_list {
 public:
  using self = dvp_adjacency_list<VertexProperty, EdgeProperty, Space,
                                  PositionMap, Arity, VPChooser, TreeStorageTag,
                                  OutEdgeListS, DirectedS, EdgeListS>;

  using point_type = topology_point_type_t<Space>;
  using point_difference_type = topology_point_difference_type_t<Space>;
  using distance_type = double;

 private:
  using alt_tag =
      ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty,
                                        EdgeProperty, EdgeListS,
                                        TreeStorageTag>;

  struct vertex_properties {};

  struct edge_properties {
    distance_type d;
  };

  using alt_type =
      typename alt_tag::template alt<vertex_properties, edge_properties>::type;

  using tree_indexer =
      typename bagl::tree_storage<vertex_properties, edge_properties,
                                  alt_tag>::type;
  using tree_vertex = bagl::graph_vertex_descriptor_t<tree_indexer>;
  using vertex_raw_map_type =
      bagl::property_map_t<tree_indexer, ReaK::graph::vertex_alt_full_prop_t>;
  using vertex_raw_property_type =
      bagl::property_traits_value_t<vertex_raw_map_type>;

  using key_map_type = decltype(tree_indexer::alt_key_map());
  using distance_map_type =
      bagl::data_member_property_map<distance_type, edge_properties>;
  using position_map_type =
      bagl::composite_property_map<PositionMap,
                                   decltype(tree_indexer::alt_adj_data_map())>;

  using dvp_impl_type =
      dvp_tree_impl<tree_indexer, Space, vertex_raw_map_type, key_map_type,
                    distance_map_type, position_map_type, Arity, VPChooser>;

  using dvp_visitor_type = typename dvp_impl_type::mutation_visitor;

  tree_indexer tree_;
  vertex_raw_map_type vp_map_;
  PositionMap position_;
  key_map_type vp_key_;
  position_map_type vp_pos_;
  dvp_impl_type impl_;

 public:
  // non-copyable: (because of shared-state)
  dvp_adjacency_list(const self&) = delete;
  self& operator=(const self&) = delete;

  using adj_list_type = ReaK::graph::alt_graph_view<alt_type>;
  using adj_list_vertex_type = bagl::graph_vertex_descriptor_t<adj_list_type>;

  /**
   * Construct the DVP-tree from a topology and property-map.
   * \param space The topology on which the positions of the vertices reside.
   * \param position The property-map that can be used to obtain the positions of the vertices.
   * \param vp_chooser The vantage-point chooser functor (policy class).
   */
  explicit dvp_adjacency_list(const std::shared_ptr<const Space>& space,
                              PositionMap position,
                              VPChooser vp_chooser = VPChooser())
      : tree_(),
        vp_map_(get(ReaK::graph::vertex_alt_full_prop, tree_)),
        position_(std::move(position)),
        vp_key_(tree_indexer::alt_key_map()),
        vp_pos_(position_, tree_indexer::alt_adj_data_map()),
        impl_(tree_, space, vp_map_, vp_key_,
              distance_map_type(&edge_properties::d), vp_pos_,
              std::move(vp_chooser)) {}

  /**
   * Move constructor.
   * \note This invalidates the adjacency-list obtained from the moved-from object.
   */
  dvp_adjacency_list(self&& rhs) noexcept
      : tree_(std::move(rhs.tree_)),
        vp_map_(get(ReaK::graph::vertex_alt_full_prop, tree_)),
        position_(std::move(rhs.position_)),
        vp_key_(std::move(rhs.vp_key_)),
        vp_pos_(std::move(rhs.vp_pos_)),
        impl_(tree_, vp_map_, std::move(rhs.impl_)) {}

  /**
   * Move assignment.
   * \note This invalidates the adjacency-list obtained from the moved-from object.
   */
  self& operator=(self&& rhs) noexcept {
    tree_ = std::move(rhs.tree_);
    vp_map_ = get(ReaK::graph::vertex_alt_full_prop, tree_);
    position_ = std::move(rhs.position_);
    vp_key_ = std::move(rhs.vp_key_);
    vp_pos_ = std::move(rhs.vp_pos_);
    impl_.reassign_moved(tree_, vp_map_, std::move(rhs.impl_));
    return *this;
  }

  /**
   * Returns a graph object associated to, stored as and synchronized with this DVP tree layout.
   * \return A graph object associated to, stored as and synchronized with this DVP tree layout.
   */
  adj_list_type get_adjacency_list() {
    return adj_list_type(tree_, dvp_visitor_type(&impl_));
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
   * Finds the nearest neighbor to a given position.
   * \param pt The position from which to find the nearest-neighbor of.
   * \return The vertex in the DVP-tree that is closest to the given point.
   */
  adj_list_vertex_type find_nearest(const point_type& pt) const {
    auto u = impl_.find_nearest(pt);
    if (u != bagl::graph_traits<tree_indexer>::null_vertex()) {
      return get(vp_key_, get(vp_map_, u));
    }
    return bagl::graph_traits<adj_list_type>::null_vertex();
  }

  /**
   * Finds the nearest predecessor and successor to a given position.
   * \param pt The position from which to find the nearest-neighbor of.
   * \return The vertices in the DVP-tree that are nearest predecessor and successor to the given point.
   */
  auto find_nearest_pred_succ(const point_type& pt) const {
    auto [u_pred, u_succ] = impl_.find_nearest_pred_succ(pt);
    std::pair<adj_list_vertex_type, adj_list_vertex_type> result;
    if (u_pred != bagl::graph_traits<tree_indexer>::null_vertex()) {
      result.first = get(vp_key_, get(vp_map_, u_pred));
    } else {
      result.first = bagl::graph_traits<adj_list_type>::null_vertex();
    }
    if (u_succ != bagl::graph_traits<tree_indexer>::null_vertex()) {
      result.second = get(vp_key_, get(vp_map_, u_succ));
    } else {
      result.second = bagl::graph_traits<adj_list_type>::null_vertex();
    }
    return result;
  };

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
  auto find_nearest(const point_type& pt, OutputIterator out_iter,
                    std::size_t num_neighbors,
                    distance_type radius =
                        std::numeric_limits<distance_type>::infinity()) const {
    std::vector<tree_vertex> v_list;
    impl_.find_nearest(pt, std::back_inserter(v_list), num_neighbors, radius);
    for (auto v : v_list) {
      *(out_iter++) = get(vp_key_, get(vp_map_, v));
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
  auto find_nearest(const point_type& pt, OutputIterator pred_iter,
                    OutputIterator succ_iter, std::size_t num_neighbors,
                    distance_type radius =
                        std::numeric_limits<distance_type>::infinity()) const {
    std::vector<tree_vertex> pred_list;
    std::vector<tree_vertex> succ_list;
    impl_.find_nearest(pt, std::back_inserter(pred_list),
                       std::back_inserter(succ_list), num_neighbors, radius);
    for (auto u : pred_list) {
      *(pred_iter++) = get(vp_key_, get(vp_map_, u));
    }
    for (auto v : succ_list) {
      *(succ_iter++) = get(vp_key_, get(vp_map_, v));
    }
    return std::pair(pred_iter, succ_iter);
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
  auto find_in_range(const point_type& pt, OutputIterator out_iter,
                     distance_type radius) const {
    std::vector<tree_vertex> v_list;
    impl_.find_in_range(pt, std::back_inserter(v_list), radius);
    for (auto v : v_list) {
      *(out_iter++) = get(vp_key_, get(vp_map_, v));
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
   * \param radius The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the two lists of nearest neighbors (predecessors and successors).
   */
  template <typename OutputIterator>
  auto find_in_range(const point_type& pt, OutputIterator pred_iter,
                     OutputIterator succ_iter, distance_type radius) const {
    std::vector<tree_vertex> pred_list;
    std::vector<tree_vertex> succ_list;
    impl_.find_in_range(pt, std::back_inserter(pred_list),
                        std::back_inserter(succ_list), radius);
    for (auto u : pred_list) {
      *(pred_iter++) = get(vp_key_, get(vp_map_, u));
    }
    for (auto v : succ_list) {
      *(succ_iter++) = get(vp_key_, get(vp_map_, v));
    }
    return std::pair(pred_iter, succ_iter);
  }
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_DVP_LAYOUT_ADJACENCY_LIST_H_
