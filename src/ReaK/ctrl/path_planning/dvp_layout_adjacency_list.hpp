/**
 * \file dvp_layout_adjacency_list.hpp
 * 
 * This library provides a class that implements a Dynamic Vantage-Point Tree (DVP-Tree) which 
 * is synchronized with a adjacency-list graph and uses the tree-storage as the layout for the 
 * vertices common two both graphs (adj-list and tree). DVP-trees
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

#ifndef REAK_DVP_LAYOUT_ADJACENCY_LIST_HPP
#define REAK_DVP_LAYOUT_ADJACENCY_LIST_HPP


#include "metric_space_search.hpp"

#include "graph_alg/adj_list_tree_overlay.hpp"

namespace ReaK {

namespace pp {


/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space. A DVP-tree is essentially
 * a generalization of a search tree which only requires the space to have a metric which 
 * respects the triangular inequality. 
 * \tparam VertexProperty The bundled vertex properties for the adjacency-list.
 * \tparam EdgeProperty The bundled edge properties for the adjacency-list.
 * \tparam Topology The topology type on which the points can reside, should model the MetricSpaceConcept.
 * \tparam PositionMap The property-map type that can map the bundled vertex property to a position in the topology.
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 * \tparam TreeStorageTag A tree-storage tag which specifies the kind of tree structure to use for the DVP tree.
 * \tparam OutEdgeListS The out-edge list container specifier for the adjacency-list (same as OutEdgeListS in boost::adjacency_list).
 * \tparam DirectedS The edge's directional specifier for the adjacency-list (same as DirectedS in boost::adjacency_list).
 * \tparam EdgeListS The edge list container specifier for the adjacency-list (same as EdgeListS in boost::adjacency_list).
 */
template <typename VertexProperty,
	  typename EdgeProperty,
	  typename Topology,
          typename PositionMap,
	  unsigned int Arity = 2,
	  typename VPChooser = random_vp_chooser,
	  typename TreeStorageTag = ReaK::graph::d_ary_bf_tree_storage<Arity>,
	  typename OutEdgeListS = boost::vecS,
	  typename DirectedS = boost::directedS,
	  typename EdgeListS = boost::vecS >
class dvp_adjacency_list
{
  public:
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef typename metric_space_traits<Topology>::distance_metric_type distance_metric;
    typedef typename topology_traits<Topology>::point_type point_type;
    typedef typename topology_traits<Topology>::point_difference_type point_difference_type;
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
      distance_metric,
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
	       get(ReaK::pp::distance_metric, *aSpace),
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
	       get(ReaK::pp::distance_metric, *aSpace),
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
#ifndef RK_ENABLE_CXX0X_FEATURES
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
  
  
};



};

};


#endif


















