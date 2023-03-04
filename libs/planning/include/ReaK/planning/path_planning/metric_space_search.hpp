/**
 * \file metric_space_search.hpp
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

#ifndef REAK_METRIC_SPACE_SEARCH_HPP
#define REAK_METRIC_SPACE_SEARCH_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

#include <limits>   // for numeric_limits
#include <utility>  // for pair and move
#include <vector>   // for vector

#include <ReaK/topologies/spaces/metric_space_concept.hpp>

// BGL-Extra includes:
#include <boost/graph/bfl_d_ary_tree.hpp>  // for default tree storage.
#include <boost/graph/more_property_maps.hpp>
#include <boost/graph/more_property_tags.hpp>

// Pending inclusion in BGL-Extra:
#include <ReaK/planning/graph_alg/bgl_raw_property_graph.hpp>

#include "dvp_tree_detail.hpp"
#include "multi_dvp_tree_search.hpp"

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
  struct position_map : boost::composite_property_map<PositionMap, KeyMap> {
    using base = boost::composite_property_map<PositionMap, KeyMap>;
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

    typename boost::property_map<
        Graph, PointType vertex_base_property<PointType>::*>::type
        m_cached_position;

    explicit effector(Graph& g)
        : m_cached_position(get(
              &vertex_base_property<PointType>::cached_position_value, g)) {}

    template <typename Vertex, typename Key, typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PointType& pt,
                    const KeyMap& key) const {
      put(key, v, k);
      put(m_cached_position, v, pt);
    }

    template <typename Vertex, typename Key, typename PositionMap,
              typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PositionMap& position,
                    const KeyMap& key) const {
      put(key, v, k);
      put(m_cached_position, v, get(position, k));
    }

    template <typename Vertex, typename PositionMap, typename KeyMap>
    PointType get_position(const Vertex& v, const PositionMap& /*unused*/,
                           const KeyMap& /*unused*/) const {
      return get(m_cached_position, v);
    }
  };

  template <typename PointType, typename KeyMap, typename PositionMap>
  struct position_map
      : boost::data_member_property_map<PointType,
                                        vertex_base_property<PointType>> {
    using base =
        boost::data_member_property_map<PointType,
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
 * \tparam Topology The topology type on which the points can reside, should model the MetricSpaceConcept.
 * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the
 * iterators) to a point (position).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename Key, typename Topology, typename PositionMap,
          unsigned int Arity = 2, typename VPChooser = random_vp_chooser,
          typename TreeStorageTag = boost::bfl_d_ary_tree_storage<Arity>,
          typename PositionCachingPolicy = position_caching_policy>
class dvp_tree {
 public:
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));

  using self = dvp_tree<Key, Topology, PositionMap, Arity, VPChooser,
                        TreeStorageTag, PositionCachingPolicy>;

  using point_type = typename boost::property_traits<PositionMap>::value_type;
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
      typename boost::tree_storage<vertex_properties, edge_properties,
                                   TreeStorageTag>::type;
  using vertex_r2b_map_type =
      typename boost::raw_vertex_to_bundle_map<tree_indexer>::type;
  using edge_r2b_map_type =
      typename boost::raw_edge_to_bundle_map<tree_indexer>::type;
  using key_map_type = boost::data_member_property_map<Key, vertex_properties>;
  using distance_map_type =
      boost::data_member_property_map<distance_type, edge_properties>;
  using vertex_position_map =
      typename PositionCachingPolicy::template position_map<
          point_type, key_map_type, PositionMap>;

  using vp_to_key_map_type =
      boost::composite_property_map<key_map_type, vertex_r2b_map_type>;
  using ep_to_distance_map_type =
      boost::composite_property_map<distance_map_type, edge_r2b_map_type>;
  using vp_to_pos_map_type =
      boost::composite_property_map<vertex_position_map, vertex_r2b_map_type>;

  using dvp_impl_type =
      dvp_tree_impl<tree_indexer, Topology, vp_to_key_map_type,
                    ep_to_distance_map_type, vp_to_pos_map_type, Arity,
                    VPChooser>;

  tree_indexer m_tree;
  PositionMap m_position;
  vp_to_key_map_type m_vp_key;
  vp_to_pos_map_type m_vp_pos;
  dvp_impl_type m_impl;

  using vertex_raw_property_type =
      typename boost::property_traits<vp_to_key_map_type>::key_type;
  using edge_raw_property_type =
      typename boost::property_traits<ep_to_distance_map_type>::key_type;

 public:
  /**
   * Construct the DVP-tree from a graph, topology and property-map.
   * \tparam Graph The graph type on which the vertices are taken from, should model the boost::VertexListGraphConcept.
   * \param g The graph from which to take the vertices.
   * \param aSpace The topology on which the positions of the vertices reside.
   * \param aPosition The property-map that can be used to obtain the positions of the vertices.
   * \param aVPChooser The vantage-point chooser functor (policy class).
   */
  template <typename Graph>
  dvp_tree(const Graph& g, const std::shared_ptr<const Topology>& aSpace,
           PositionMap aPosition, VPChooser aVPChooser = VPChooser())
      : m_tree(),
        m_position(aPosition),
        m_vp_key(key_map_type(&vertex_properties::k),
                 get_raw_vertex_to_bundle_map(m_tree)),
        m_vp_pos(
            vertex_position_map(key_map_type(&vertex_properties::k), aPosition),
            get_raw_vertex_to_bundle_map(m_tree)),
        m_impl(g, aPosition, m_tree, aSpace, m_vp_key,
               ep_to_distance_map_type(distance_map_type(&edge_properties::d),
                                       get_raw_edge_to_bundle_map(m_tree)),
               m_vp_pos, aVPChooser) {}

  /**
   * Construct the DVP-tree from a range, topology and property-map.
   * \tparam ForwardIterator The forward-iterator type from which the vertices can be obtained.
   * \param aBegin The start of the range from which to take the vertices.
   * \param aEnd The end of the range from which to take the vertices (one-past-last).
   * \param aSpace The topology on which the positions of the vertices reside.
   * \param aPosition The property-map that can be used to obtain the positions of the vertices.
   * \param aVPChooser The vantage-point chooser functor (policy class).
   */
  template <typename ForwardIterator>
  dvp_tree(ForwardIterator aBegin, ForwardIterator aEnd,
           const std::shared_ptr<const Topology>& aSpace, PositionMap aPosition,
           VPChooser aVPChooser = VPChooser())
      : m_tree(),
        m_position(aPosition),
        m_vp_key(key_map_type(&vertex_properties::k),
                 get_raw_vertex_to_bundle_map(m_tree)),
        m_vp_pos(
            vertex_position_map(key_map_type(&vertex_properties::k), aPosition),
            get_raw_vertex_to_bundle_map(m_tree)),
        m_impl(aBegin, aEnd, aPosition, m_tree, aSpace, m_vp_key,
               ep_to_distance_map_type(distance_map_type(&edge_properties::d),
                                       get_raw_edge_to_bundle_map(m_tree)),
               m_vp_pos, aVPChooser) {}

  dvp_tree(const self& rhs)
      : m_tree(rhs.m_tree),
        m_position(rhs.m_position),
        m_vp_key(rhs.m_vp_key),
        m_vp_pos(rhs.m_vp_pos),
        m_impl(m_tree, rhs.m_impl) {}

  self& operator=(const self& rhs) {
    m_tree = rhs.m_tree;
    m_position = rhs.m_position;
    m_vp_key = rhs.m_vp_key;
    m_vp_pos = rhs.m_vp_pos;
    m_impl.reassign_copied(m_tree, rhs.m_impl);
    return *this;
  }

  dvp_tree(self&& rhs) noexcept
      : m_tree(std::move(rhs.m_tree)),
        m_position(std::move(rhs.m_position)),
        m_vp_key(std::move(rhs.m_vp_key)),
        m_vp_pos(std::move(rhs.m_vp_pos)),
        m_impl(m_tree, std::move(rhs.m_impl)) {}

  self& operator=(self&& rhs) noexcept {
    m_tree = std::move(rhs.m_tree);
    m_position = std::move(rhs.m_position);
    m_vp_key = std::move(rhs.m_vp_key);
    m_vp_pos = std::move(rhs.m_vp_pos);
    m_impl.reassign_moved(m_tree, std::move(rhs.m_impl));
    return *this;
  }

  /**
   * Checks if the DVP-tree is empty.
   * \return True if the DVP-tree is empty.
   */
  bool empty() const { return m_impl.empty(); }

  /**
   * Returns the size of the DVP-tree (the number of vertices it contains.
   * \return The size of the DVP-tree (the number of vertices it contains.
   */
  std::size_t size() const { return m_impl.size(); }

  /**
   * Returns the depth of the tree.
   * \note This operation must recurse through all the branches of the tree (depth-first), and is
   * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack)
   * w.r.t. the depth of tree).
   * \return The depth of the tree.
   */
  std::size_t depth() const { return m_impl.depth(); }

  /**
   * This function computes an approximation of the characteristic size of the vertices in the DVP tree.
   * \return The approximation of the characteristic size of the vertices in the DVP tree.
   */
  double get_characteristic_size() const {
    return m_impl.get_characteristic_size();
  }

  /**
   * Inserts a key-value (vertex).
   * \param u The vertex to be added to the DVP-tree.
   */
  void insert(Key u) {
    vertex_raw_property_type vp;
    put(m_vp_key, vp, u);
    put(m_vp_pos, vp, get(m_position, u));
    m_impl.insert(std::move(vp));
  }

  /**
   * Inserts a range of key-values (vertices).
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
   * \param aBegin The start of the range from which to take the vertices.
   * \param aEnd The end of the range from which to take the vertices (one-past-last).
   */
  template <typename ForwardIterator>
  void insert(ForwardIterator aBegin, ForwardIterator aEnd) {
    std::vector<vertex_raw_property_type> vp_list;
    for (; aBegin != aEnd; ++aBegin) {
      vp_list.push_back(vertex_raw_property_type());
      put(m_vp_key, vp_list.back(), *aBegin);
      put(m_vp_pos, vp_list.back(), get(m_position, *aBegin));
    }
    m_impl.insert(vp_list.begin(), vp_list.end());
  }

  /**
   * Erases the given vertex from the DVP-tree.
   * \param u The vertex to be removed from the DVP-tree.
   */
  void erase(Key u) { m_impl.erase(u, get(m_position, u)); }

  /**
   * Erases the given vertex-range from the DVP-tree.
   * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
   * \param aBegin The start of the range from which to take the vertices to be erased.
   * \param aEnd The end of the range from which to take the vertices to be erased (one-past-last).
   */
  template <typename ForwardIterator>
  void erase(ForwardIterator aBegin, ForwardIterator aEnd) {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    for (; aBegin != aEnd; ++aBegin) {
      TreeVertex u = m_impl.get_vertex(*aBegin, get(m_position, *aBegin));
      if (u != boost::graph_traits<tree_indexer>::null_vertex()) {
        v_list.push_back(u);
      }
    }
    m_impl.erase(v_list.begin(), v_list.end());
  }

  /**
   * Clears the DVP-tree.
   */
  void clear() { m_impl.clear(); }

  /**
   * Finds the nearest neighbor to a given position.
   * \param aPoint The position from which to find the nearest-neighbor of.
   * \return The vertex in the DVP-tree that is closest to the given point.
   */
  Key find_nearest(const point_type& aPoint) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    TreeVertex u = m_impl.find_nearest(aPoint);
    if (u != boost::graph_traits<tree_indexer>::null_vertex()) {
      return m_tree[u].k;
    }
    return Key();
  }

  /**
   * Finds the nearest predecessor and successor to a given position.
   * \param aPoint The position from which to find the nearest-neighbor of.
   * \return The vertices in the DVP-tree that are nearest predecessor and successor to the given point.
   */
  std::pair<Key, Key> find_nearest_pred_succ(const point_type& aPoint) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    auto [u_pred, u_succ] = m_impl.find_nearest_pred_succ(aPoint);
    std::pair<Key, Key> result;
    if (u_pred != boost::graph_traits<tree_indexer>::null_vertex()) {
      result.first = m_tree[u_pred].k;
    } else {
      result.first = Key();
    }
    if (u_succ != boost::graph_traits<tree_indexer>::null_vertex()) {
      result.second = m_tree[u_succ].k;
    } else {
      result.second = Key();
    }
    return result;
  }

  /**
   * Finds the K nearest-neighbors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param aPoint The position from which to find the nearest-neighbors.
   * \param aOutputBegin An iterator to the first place where to put the sorted list of
   *        elements with the smallest distance.
   * \param K The number of nearest-neighbors.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_nearest(
      const point_type& aPoint, OutputIterator aOutputBegin, std::size_t K,
      distance_type R = std::numeric_limits<distance_type>::infinity()) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    m_impl.find_nearest(aPoint, back_inserter(v_list), K, R);
    for (TreeVertex v : v_list) {
      *(aOutputBegin++) = m_tree[v].k;
    }
    return aOutputBegin;
  };

  /**
   * Finds the K nearest predecessors and successors to a given position.
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param aPoint The position from which to find the nearest-neighbors.
   * \param aPredBegin An iterator to the first place where to put the sorted list of
   *        predecessor elements with the smallest distance.
   * \param aSuccBegin An iterator to the first place where to put the sorted list of
   *        successor elements with the smallest distance.
   * \param K The number of nearest-neighbors.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the two lists of nearest neighbors (predecessors and successors).
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_nearest(
      const point_type& aPoint, OutputIterator aPredBegin,
      OutputIterator aSuccBegin, std::size_t K,
      distance_type R = std::numeric_limits<distance_type>::infinity()) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    std::vector<TreeVertex> pred_list;
    std::vector<TreeVertex> succ_list;
    m_impl.find_nearest(aPoint, back_inserter(pred_list),
                        back_inserter(succ_list), K, R);
    for (TreeVertex v : pred_list) {
      *(aPredBegin++) = m_tree[v].k;
    }
    for (TreeVertex v : succ_list) {
      *(aSuccBegin++) = m_tree[v].k;
    }
    return {aPredBegin, aSuccBegin};
  }

  /**
   * Finds the nearest-neighbors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param aPoint The position from which to find the nearest-neighbors.
   * \param aOutputBegin An iterator to the first place where to put the sorted list of
   *        elements with the smallest distance.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename OutputIterator>
  OutputIterator find_in_range(const point_type& aPoint,
                               OutputIterator aOutputBegin,
                               distance_type R) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    std::vector<TreeVertex> v_list;
    m_impl.find_in_range(aPoint, back_inserter(v_list), R);
    for (TreeVertex v : v_list) {
      *(aOutputBegin++) = m_tree[v].k;
    }
    return aOutputBegin;
  }

  /**
   * Finds the K nearest predecessors and successors to a given position within a given range (radius).
   * \tparam OutputIterator The forward- output-iterator type which can contain the
   *         list of nearest-neighbors.
   * \param aPoint The position from which to find the nearest-neighbors.
   * \param aPredBegin An iterator to the first place where to put the sorted list of
   *        predecessor elements with the smallest distance.
   * \param aSuccBegin An iterator to the first place where to put the sorted list of
   *        successor elements with the smallest distance.
   * \param R The maximum distance value for the nearest-neighbors.
   * \return The output-iterator to the end of the two lists of nearest neighbors (predecessors and successors).
   */
  template <typename OutputIterator>
  std::pair<OutputIterator, OutputIterator> find_in_range(
      const point_type& aPoint, OutputIterator aPredBegin,
      OutputIterator aSuccBegin, distance_type R) const {
    using TreeVertex = graph::graph_vertex_t<tree_indexer>;
    std::vector<TreeVertex> pred_list;
    std::vector<TreeVertex> succ_list;
    m_impl.find_in_range(aPoint, back_inserter(pred_list),
                         back_inserter(succ_list), R);
    for (TreeVertex v : pred_list) {
      *(aPredBegin++) = m_tree[v].k;
    }
    for (TreeVertex v : succ_list) {
      *(aSuccBegin++) = m_tree[v].k;
    }
    return {aPredBegin, aSuccBegin};
  }
};

}  // namespace ReaK::pp

#endif
