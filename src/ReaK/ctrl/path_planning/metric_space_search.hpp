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

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include <map>
#include <unordered_map>
#include <vector>
#include "metric_space_concept.hpp"
#include "global_rng.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "topological_search.hpp"

#include "dvp_tree_detail.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "graph_alg/bgl_raw_property_graph.hpp"


namespace ReaK {

namespace pp {

/**
 * This class is a callable class that can be used to choose the best 
 * vantage-point to use out of a set of points. In theory, the best vantage-point 
 * is the one which deviates the most from the other points in the set, however, 
 * this functor will approximately select that point by searching for it only 
 * in a random subset of the given range of points.
 */
class random_best_vp_chooser {
  private:
    unsigned int m_divider;
  public:
  
    /**
     * Default construction.
     * \param aDivider The divider of the set (determines the fraction of the points to search), default is 10.
     */
    random_best_vp_chooser(unsigned int aDivider = 10) : m_divider(aDivider) { };
    
    /**
     * This call-operator will choose a vantage-point from within the given range.
     * \tparam RandomAccessIter A random-access iterator type that can describe the point-range.
     * \tparam Topology The topology type on which the points can reside.
     * \tparam DistanceMetric The distance-metric type over the topology type.
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aSpace The topology on which the points reside.
     * \param aDistance The distance-metric over the given topology.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename Topology, typename DistanceMetric, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, const Topology& aSpace, DistanceMetric aDistance, PositionMap aPosition) {
      BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
      typedef typename topology_traits<Topology>::point_type Point;
      RandomAccessIter best_pt = aEnd;
      double best_dev = -1;
      for(unsigned int i=0; i < (aEnd - aBegin) / m_divider + 1;++i) {
	RandomAccessIter current_pt = aBegin + (get_global_rng()() % (aEnd - aBegin));
	double current_mean = 0.0;
	double current_dev = 0.0;
	Point current_vp = get(aPosition, *current_pt);
	for(unsigned int j=0; aBegin + j != aEnd; ++j) {
	  double dist = aDistance(current_vp, get(aPosition, *(aBegin + j)), aSpace);
	  current_mean = (current_mean * j + dist) / (j + 1);
	  current_dev = (current_dev * j + dist * dist) / (j + 1);
	};
	double current_var = current_dev - current_mean * current_mean;
	if(current_var < 0) current_var = 0.0;
	current_dev = std::sqrt(current_var);
	
	if(current_dev > best_dev) {
	  best_pt = current_pt;
	  best_dev = current_dev;
	};
      };
      return best_pt;
    };
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
    random_vp_chooser() { };
    
    /**
     * This call-operator will choose a vantage-point from within the given range.
     * \tparam RandomAccessIter A random-access iterator type that can describe the point-range.
     * \tparam Topology The topology type on which the points can reside, should model the MetricSpaceConcept.
     * \tparam DistanceMetric The distance-metric type over the topology type.
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aSpace The topology on which the points reside.
     * \param aDistance The distance-metric over the given topology.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename Topology, typename DistanceMetric, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, const Topology& aSpace, DistanceMetric aDistance, PositionMap aPosition) {
      return aBegin + (get_global_rng()() % (aEnd - aBegin));
    };
};

/**
 * This class is a position-caching policy class for the space partitioning trees that are indirectly 
 * indexing a set of vertices of a graph (or other container). When building an indirect index, it can
 * be advantageous to store (copy) the position values (or vector) in the index itself to increase
 * the locality of reference during the traversal-mutating operations (lookups, insertions, etc.) 
 * through the space partitioning tree. This position-caching policy mandates that no such caching 
 * is done.
 */
struct no_position_caching_policy {
  template < typename PointType >
  struct vertex_base_property { }; // intentionally empty.
  
  template <typename Graph, typename PointType>
  struct effector {
    effector(Graph&) { };
    
    template <typename Vertex, typename Key, typename PositionMap, typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PositionMap&, const KeyMap& key) const {
      put(key, v, k);
    };
  
    template <typename Vertex, typename PositionMap, typename KeyMap>
    PointType get_position(const Vertex& v, const PositionMap& position, const KeyMap& key) const {
      return get(position, get(key, v));
    };
  };
  
  template <typename PointType, typename KeyMap, typename PositionMap>
  struct position_map : boost::composite_property_map<PositionMap, KeyMap>  {
    typedef boost::composite_property_map<PositionMap, KeyMap> base; 
    typedef position_map<PointType, KeyMap, PositionMap> self;
    
    position_map(KeyMap k = KeyMap(), PositionMap p = PositionMap()) : base(p,k) { };
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
    
    typename boost::property_map< Graph, PointType vertex_base_property<PointType>::* >::type m_cached_position;
    
    effector(Graph& g) : m_cached_position( get(&vertex_base_property<PointType>::cached_position_value, g) ) { };
    
    template <typename Vertex, typename Key, typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PointType& pt, const KeyMap& key) const {
      put(key, v, k);
      put(m_cached_position, v, pt);
    };
    
    template <typename Vertex, typename Key, typename PositionMap, typename KeyMap>
    void put_vertex(const Vertex& v, const Key& k, const PositionMap& position, const KeyMap& key) const {
      put(key, v, k);
      put(m_cached_position, v, get(position, k));
    };
  
    template <typename Vertex, typename PositionMap, typename KeyMap>
    PointType get_position(const Vertex& v, const PositionMap&, const KeyMap&) const {
      return get(m_cached_position, v);
    };
  };
  
  template <typename PointType, typename KeyMap, typename PositionMap>
  struct position_map : boost::data_member_property_map<PointType, vertex_base_property<PointType> >  {
    typedef boost::data_member_property_map<PointType, vertex_base_property<PointType> > base; 
    typedef position_map<PointType,KeyMap,PositionMap> self;
    
    position_map(KeyMap k = KeyMap(), PositionMap p = PositionMap()) : base(&vertex_base_property<PointType>::cached_position_value) { };
  };
};



/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which 
 * respects the triangular inequality. 
 * \tparam Key The key type for the tree, essentially the key value is the vertex descriptor type.
 * \tparam Topology The topology type on which the points can reside, should model the MetricSpaceConcept.
 * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename Key,
          typename Topology,
          typename PositionMap,
          unsigned int Arity = 2,
          typename VPChooser = random_vp_chooser,
          typename TreeStorageTag = ReaK::graph::d_ary_bf_tree_storage<Arity>,
          typename PositionCachingPolicy = position_caching_policy>
class dvp_tree
{
  public:
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef typename boost::mpl::if_<
      is_metric_symmetric<Topology>,
      typename metric_space_traits<Topology>::distance_metric_type,
      symmetrized_metric< typename metric_space_traits<Topology>::distance_metric_type > >::type distance_metric_type;
    typedef typename boost::property_traits<PositionMap>::value_type point_type;
    typedef double distance_type;
    
  private:
    
    struct vertex_properties : PositionCachingPolicy::template vertex_base_property<point_type> {
      Key k;
    };
    
    struct edge_properties {
      distance_type d;
    };
    
    typedef typename ReaK::graph::tree_storage< vertex_properties, edge_properties, TreeStorageTag>::type tree_indexer;
    typedef typename boost::property_map<tree_indexer, boost::vertex_raw_prop_to_bundle_t>::type vertex_r2b_map_type;
    typedef typename boost::property_map<tree_indexer, boost::edge_raw_prop_to_bundle_t>::type edge_r2b_map_type;
    typedef boost::data_member_property_map<Key, vertex_properties> key_map_type;
    typedef boost::data_member_property_map<distance_type, edge_properties> distance_map_type;
    typedef typename PositionCachingPolicy::template position_map<point_type, key_map_type, PositionMap> vertex_position_map;
    
    typedef boost::composite_property_map<key_map_type, vertex_r2b_map_type> vp_to_key_map_type;
    typedef boost::composite_property_map<distance_map_type, edge_r2b_map_type> ep_to_distance_map_type;
    typedef boost::composite_property_map<vertex_position_map, vertex_r2b_map_type> vp_to_pos_map_type;
    
    
    typedef dvp_tree_impl< 
      tree_indexer,
      Topology,
      distance_metric_type,
      vp_to_key_map_type,
      ep_to_distance_map_type,
      vp_to_pos_map_type,
      Arity,
      VPChooser> dvp_impl_type;
    
    tree_indexer m_tree;
    PositionMap m_position;
    vp_to_key_map_type m_vp_key;
    vp_to_pos_map_type m_vp_pos;
    dvp_impl_type m_impl;
    
    typedef typename boost::property_traits< vp_to_key_map_type >::key_type vertex_raw_property_type;
    typedef typename boost::property_traits< ep_to_distance_map_type >::key_type edge_raw_property_type;
    
    
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
    dvp_tree(const Graph& g, 
             const ReaK::shared_ptr<const Topology>& aSpace, 
             PositionMap aPosition,
             VPChooser aVPChooser = VPChooser()) :
             m_tree(),
             m_position(aPosition),
             m_vp_key(
               key_map_type(&vertex_properties::k),
               get(boost::vertex_raw_prop_to_bundle, m_tree)
             ),
             m_vp_pos(
               vertex_position_map(key_map_type(&vertex_properties::k), aPosition),
               get(boost::vertex_raw_prop_to_bundle, m_tree)
             ),
             m_impl(
               g,
               aPosition,
               m_tree,
               aSpace, 
               distance_metric_type(get(ReaK::pp::distance_metric,*aSpace)),
               m_vp_key,
               ep_to_distance_map_type(
                 distance_map_type(&edge_properties::d),
                 get(boost::edge_raw_prop_to_bundle, m_tree)
               ),
               m_vp_pos,
               aVPChooser
             ) { };
    
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
    dvp_tree(ForwardIterator aBegin,
             ForwardIterator aEnd,
             const ReaK::shared_ptr<const Topology>& aSpace, 
             PositionMap aPosition,
             VPChooser aVPChooser = VPChooser()) : 
             m_tree(),
             m_position(aPosition),
             m_vp_key(
               key_map_type(&vertex_properties::k),
               get(boost::vertex_raw_prop_to_bundle, m_tree)
             ),
             m_vp_pos(
               vertex_position_map(key_map_type(&vertex_properties::k), aPosition),
               get(boost::vertex_raw_prop_to_bundle, m_tree)
             ),
             m_impl(
               aBegin,
               aEnd,
               aPosition,
               m_tree,
               aSpace, 
               distance_metric_type(get(ReaK::pp::distance_metric, *aSpace)),
               m_vp_key,
               ep_to_distance_map_type(
                 distance_map_type(&edge_properties::d),
                 get(boost::edge_raw_prop_to_bundle, m_tree)
               ),
               m_vp_pos,
               aVPChooser
             ) { };
  
    
    
    /**
     * Checks if the DVP-tree is empty.
     * \return True if the DVP-tree is empty.
     */
    bool empty() const { return m_impl.empty(); };
    
    /**
     * Returns the size of the DVP-tree (the number of vertices it contains.
     * \return The size of the DVP-tree (the number of vertices it contains.
     */
    std::size_t size() const { return m_impl.size(); };
    
    /**
     * Returns the depth of the tree.
     * \note This operation must recurse through all the branches of the tree (depth-first), and is 
     * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack) 
     * w.r.t. the depth of tree).
     * \return The depth of the tree.
     */
    std::size_t depth() const { return m_impl.depth(); };
    
    /**
     * Inserts a key-value (vertex).
     * \param u The vertex to be added to the DVP-tree.
     */
    void insert(Key u) {
      vertex_raw_property_type vp;
      put(m_vp_key, vp, u);
      put(m_vp_pos, vp, get(m_position, u));
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
      m_impl.insert(vp);
#else
      m_impl.insert(std::move(vp));
#endif
    };
    
    /**
     * Inserts a range of key-values (vertices).
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     */
    template <typename ForwardIterator>
    void insert(ForwardIterator aBegin, ForwardIterator aEnd) { 
      std::vector< vertex_raw_property_type > vp_list;
      for(; aBegin != aEnd; ++aBegin) {
        vp_list.push_back(vertex_raw_property_type());
        put(m_vp_key, vp_list.back(), *aBegin);
        put(m_vp_pos, vp_list.back(), get(m_position, *aBegin));
      };
      m_impl.insert(vp_list.begin(), vp_list.end());
    };
    
    /**
     * Erases the given vertex from the DVP-tree.
     * \param u The vertex to be removed from the DVP-tree.
     */
    void erase(Key u) { 
      m_impl.erase(u, get(m_position, u));
    };
    
    /**
     * Erases the given vertex-range from the DVP-tree.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices to be erased.
     * \param aEnd The end of the range from which to take the vertices to be erased (one-past-last).
     */
    template <typename ForwardIterator>
    void erase(ForwardIterator aBegin, ForwardIterator aEnd) { 
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      std::vector< TreeVertex > v_list;
      for(; aBegin != aEnd; ++aBegin) {
        TreeVertex u = m_impl.get_vertex(*aBegin, get(m_position, *aBegin));
        if(is_vertex_valid(u, m_tree))
          v_list.push_back(u);
      };
      m_impl.erase(v_list.begin(), v_list.end());
    };
    
    /**
     * Clears the DVP-tree. 
     */
    void clear() {
      m_impl.clear();
    };
    
    /**
     * Finds the nearest neighbor to a given position.
     * \param aPoint The position from which to find the nearest-neighbor of.
     * \return The vertex in the DVP-tree that is closest to the given point.
     */
    Key find_nearest(const point_type& aPoint) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      TreeVertex u = m_impl.find_nearest(aPoint);
      return m_tree[u].k;
    };
    
    /**
     * Finds the nearest predecessor and successor to a given position.
     * \note This only works for an unsymmetric metric (symmetrized in this DVP tree, but unsymmetrized for this query).
     * \param aPoint The position from which to find the nearest-neighbor of.
     * \return The vertices in the DVP-tree that are nearest predecessor and successor to the given point.
     */
    std::pair<Key, Key> find_nearest_pred_succ(const point_type& aPoint) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      std::pair<TreeVertex, TreeVertex> u = m_impl.find_nearest_pred_succ(aPoint);
      return std::pair<Key, Key>(m_tree[u.first].k, m_tree[u.second].k);
    };
    
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
    OutputIterator find_nearest(const point_type& aPoint, OutputIterator aOutputBegin, std::size_t K, distance_type R = std::numeric_limits<distance_type>::infinity()) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      std::vector< TreeVertex > v_list;
      m_impl.find_nearest(aPoint, back_inserter(v_list), K, R);
      for(typename std::vector< TreeVertex >::iterator it = v_list.begin(); it != v_list.end(); ++it)
        *(aOutputBegin++) = m_tree[*it].k;
      return aOutputBegin;
    };
    
    /**
     * Finds the K nearest predecessors and successors to a given position.
     * \note This only works for an unsymmetric metric (symmetrized in this DVP tree, but unsymmetrized for this query).
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
    std::pair<OutputIterator, OutputIterator> find_nearest(const point_type& aPoint, OutputIterator aPredBegin, OutputIterator aSuccBegin, std::size_t K, distance_type R = std::numeric_limits<distance_type>::infinity()) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      std::vector< TreeVertex > pred_list;
      std::vector< TreeVertex > succ_list;
      m_impl.find_nearest(aPoint, back_inserter(pred_list), back_inserter(succ_list), K, R);
      for(typename std::vector< TreeVertex >::iterator it = pred_list.begin(); it != pred_list.end(); ++it)
        *(aPredBegin++) = m_tree[*it].k;
      for(typename std::vector< TreeVertex >::iterator it = succ_list.begin(); it != succ_list.end(); ++it)
        *(aSuccBegin++) = m_tree[*it].k;
      return std::pair<OutputIterator, OutputIterator>(aPredBegin,aSuccBegin);
    };
    
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
    OutputIterator find_in_range(const point_type& aPoint, OutputIterator aOutputBegin, distance_type R) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      std::vector< TreeVertex > v_list;
      m_impl.find_in_range(aPoint, back_inserter(v_list), R);
      for(typename std::vector< TreeVertex >::iterator it = v_list.begin(); it != v_list.end(); ++it)
        *(aOutputBegin++) = m_tree[*it].k;
      return aOutputBegin;
    };
    
    
};


/**
 * This class template is used as a type of kind of associative container for DVP-trees that 
 * span distinct graphs. For problems in which multiple graphs exist and nearest-neighbors 
 * may be queried from any of these graphs, one can use this class to associate each graph 
 * with a dvp-tree. This class is callable as a single query and a KNN / range query.
 * \tparam Graph The graph type which can contain the vertices.
 * \tparam DVPTree The DVP-tree type that is used to perform the nearest-neighbor queries.
 */
template <typename Graph, typename DVPTree>
struct multi_dvp_tree_search {
  /** This associative container is used to store and access the DVP-trees. */
  typename std::map<Graph*, DVPTree*> graph_tree_map;
  
  multi_dvp_tree_search() : graph_tree_map() { };
  
  
  /**
   * This call-operator finds the nearest vertex of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param p A position in the space, to which the nearest-neighbor is sought.
   * \param g A graph containing the vertices from which to find the nearest-neighbor, 
   *        should be tree-structured.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   */
  template <typename Topology, typename PositionMap>
  typename boost::graph_traits<Graph>::vertex_descriptor operator()(const typename boost::property_traits<PositionMap>::value_type& p, 
                                                                    Graph& g, const Topology& space, PositionMap position) const {
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      return it->second->find_nearest(p);
    else
      return typename boost::graph_traits<Graph>::vertex_descriptor();
  };
  
  /**
   * This call-operator finds the nearest vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions.
   * \tparam PositionMap The property-map type which can store the position associated 
   *         with each vertex.
   * \tparam OutputIterator The forward- output-iterator type which can contain the 
   *         list of nearest-neighbors.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param output_first An iterator to the first place where to put the sorted list of 
   *        elements with the smallest distance.
   * \param g A graph containing the vertices from which to find the nearest-neighbors, 
   *        should be tree-structured.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \param max_neighbors The maximum number of neighbors to have in the list.
   * \param radius The minimum distance around the position that a vertex should be in to be 
   *        considered a neighbor.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename Topology, typename PositionMap, typename OutputIterator>
  OutputIterator operator()(const typename boost::property_traits<PositionMap>::value_type& p, 
                            OutputIterator output_first, 
                            Graph& g, const Topology& space, PositionMap position, 
                            std::size_t max_neighbors = 1, double radius = std::numeric_limits<double>::infinity()) const {
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      return it->second->find_nearest(p, output_first, max_neighbors, radius);
    else
      return output_first;
  };
  
  
  /**
   * This is a call-back for when a vertex has been added.
   */
  template <typename Vertex>
  void added_vertex(Vertex v, Graph& g) const { 
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      it->second->insert(v);
  };
  
  /**
   * This is a call-back for when a vertex is about to be removed.
   */
  template <typename Vertex>
  void removed_vertex(Vertex v, Graph& g) const { 
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      it->second->erase(v);
  };
  
};



/**
 * This class template is used as a type of kind of associative container for DVP-trees that 
 * span distinct graphs. For problems in which multiple graphs exist and nearest-neighbors 
 * may be queried from any of these graphs, one can use this class to associate each graph 
 * with a dvp-tree. This class is callable as a single query and a KNN / range query.
 * \tparam Graph The graph type which can contain the vertices.
 * \tparam DVPTree The DVP-tree type that is used to perform the nearest-neighbor queries.
 */
template <typename Graph, typename DVPTree>
struct multi_dvp_tree_pred_succ_search {
  /** This associative container is used to store and access the DVP-trees. */
  typename std::map<Graph*, DVPTree*> graph_tree_map;
  
  multi_dvp_tree_pred_succ_search() : graph_tree_map() { };
  
  /**
   * This call-operator finds the nearest predecesor and successor of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions, should model the MetricSpaceConcept.
   * \tparam PositionMap The property-map type which can store the position associated with each vertex.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param g A graph containing the vertices from which to find the nearest-neighbor, 
   *        should be tree-structured.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \return A pair containing the nearest predecessor and successor vertex.
   */
  template <typename Topology, typename PositionMap>
  std::pair< typename boost::graph_traits<Graph>::vertex_descriptor, 
             typename boost::graph_traits<Graph>::vertex_descriptor > 
    operator()(const typename boost::property_traits<PositionMap>::value_type& p, 
               Graph& g, const Topology& space, PositionMap position) const {
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      return it->second->find_nearest_pred_succ(p);
    else
      return std::pair< typename boost::graph_traits<Graph>::vertex_descriptor, 
                        typename boost::graph_traits<Graph>::vertex_descriptor >();
  };
  
  /**
   * This call-operator finds the nearest vertices of a graph, to a given position.
   * \tparam Topology The topology type which contains the positions.
   * \tparam PositionMap The property-map type which can store the position associated 
   *         with each vertex.
   * \tparam OutputIterator The forward- output-iterator type which can contain the 
   *         list of nearest-neighbors.
   * \param p A position in the space, to which the nearest-neighbors are sought.
   * \param pred_first An iterator to the first place where to put the sorted list of 
   *        the best predecessor elements.
   * \param succ_first An iterator to the first place where to put the sorted list of 
   *        the best successor elements.
   * \param g A graph containing the vertices from which to find the nearest-neighbors, 
   *        should be tree-structured.
   * \param space The topology objects which define the space in which the positions reside.
   * \param position The property-map which can retrieve the position associated to each vertex.
   * \param max_neighbors The maximum number of neighbors to have in the list.
   * \param radius The minimum distance around the position that a vertex should be in to be 
   *        considered a neighbor.
   * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
   */
  template <typename Topology, typename PositionMap, typename OutputIterator>
  std::pair< OutputIterator, OutputIterator > 
    operator()(const typename boost::property_traits<PositionMap>::value_type& p, 
               OutputIterator pred_first, 
               OutputIterator succ_first, 
               Graph& g, const Topology& space, PositionMap position, 
               std::size_t max_neighbors = 1, double radius = std::numeric_limits<double>::infinity()) const {
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      return it->second->find_nearest(p, pred_first, succ_first, max_neighbors, radius);
    else
      return std::pair< OutputIterator, OutputIterator >(pred_first, succ_first);
  };
  
  
  /**
   * This is a call-back for when a vertex has been added.
   */
  template <typename Vertex>
  void added_vertex(Vertex v, Graph& g) const { 
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      it->second->insert(v);
  };
  
  /**
   * This is a call-back for when a vertex is about to be removed.
   */
  template <typename Vertex>
  void removed_vertex(Vertex v, Graph& g) const { 
    typename std::map<Graph*,DVPTree*>::const_iterator it = graph_tree_map.find(&g);
    if((it != graph_tree_map.end()) && (it->second))
      it->second->erase(v);
  };
  
};



};

};


#endif


















