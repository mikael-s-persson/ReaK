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

#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include "graph_alg/adj_list_tree_overlay.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"

#include "metric_space_search.hpp"



namespace ReaK {

namespace pp {


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
    
    typedef ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, EdgeProperty, EdgeListS, TreeStorageTag> alt_tag;
    
    struct vertex_properties { };
    
    struct edge_properties {
      distance_type d;
    };
    
    typedef typename alt_tag::template alt<vertex_properties, edge_properties>::type alt_type;
    
    typedef typename ReaK::graph::tree_storage< vertex_properties, edge_properties, alt_tag>::type tree_indexer;
    typedef typename boost::property_map<tree_indexer, boost::vertex_raw_prop_to_bundle_t>::type vertex_r2b_map_type;
    typedef typename boost::property_map<tree_indexer, boost::edge_raw_prop_to_bundle_t>::type edge_r2b_map_type;
    typedef typename boost::property_traits< vertex_r2b_map_type >::key_type vertex_raw_property_type;
    typedef typename boost::property_traits< edge_r2b_map_type >::key_type edge_raw_property_type;
    
	typedef typename boost::property_map<vertex_raw_property_type, boost::vertex_key_t>::const_type key_map_type; 
    typedef boost::composite_property_map<
      boost::data_member_property_map<distance_type, edge_properties>,
      edge_r2b_map_type> distance_map_type;
    typedef boost::composite_property_map<
      PositionMap,
	  typename boost::property_map<vertex_raw_property_type, boost::vertex_second_bundle_t>::type > position_map_type;
    
    
    typedef dvp_tree_impl< 
      tree_indexer,
      Topology,
      distance_metric,
      key_map_type,
      distance_map_type,
      position_map_type,
      Arity,
      VPChooser> dvp_impl_type;
      
    typedef typename dvp_impl_type::mutation_visitor dvp_visitor_type;
    
    tree_indexer m_tree;
    PositionMap m_position;
    key_map_type m_vp_key;
    position_map_type m_vp_pos;
    dvp_impl_type m_impl;
    
    
  public:
    
    typedef ReaK::graph::alt_graph_view< alt_type > adj_list_type;
    typedef typename boost::graph_traits< adj_list_type >::vertex_descriptor adj_list_vertex_type;
    
    /**
     * Construct the DVP-tree from a topology and property-map.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aPosition The property-map that can be used to obtain the positions of the vertices.
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    dvp_adjacency_list(const ReaK::shared_ptr<const Topology>& aSpace, 
		       PositionMap aPosition,
		       VPChooser aVPChooser = VPChooser()) :
	     m_tree(),
	     m_position(aPosition),
	     m_vp_key(
	       get(boost::vertex_key, vertex_raw_property_type())
	     ),
	     m_vp_pos(
	       aPosition,
	       get(boost::vertex_second_bundle, vertex_raw_property_type())
	     ),
	     m_impl(
	       m_tree,
               aSpace, 
	       get(ReaK::pp::distance_metric,*aSpace),
	       m_vp_key,
	       distance_map_type(
		 boost::data_member_property_map<distance_type, edge_properties>(&edge_properties::d),
		 get(boost::edge_raw_prop_to_bundle, m_tree)
	       ),
	       m_vp_pos,
	       aVPChooser
	     ) { };
    
    /**
     * Returns a graph object associated to, stored as and synchronized with this DVP tree layout.
     * \return A graph object associated to, stored as and synchronized with this DVP tree layout.
     */
    adj_list_type get_adjacency_list() {
      return adj_list_type(m_tree, dvp_visitor_type(&m_impl));
    };
    
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
     * Finds the nearest neighbor to a given position.
     * \param aPoint The position from which to find the nearest-neighbor of.
     * \return The vertex in the DVP-tree that is closest to the given point.
     */
    adj_list_vertex_type find_nearest(const point_type& aPoint) const {
      typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor TreeVertex;
      TreeVertex u = m_impl.find_nearest(aPoint);
      return get(m_vp_key, get(boost::vertex_raw_property, m_tree, u));
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
	*(aOutputBegin++) = get(m_vp_key, get(boost::vertex_raw_property, m_tree, *it));
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
	*(aOutputBegin++) = get(m_vp_key, get(boost::vertex_raw_property, m_tree, *it));
      return aOutputBegin;
    };
    
    
};





};

};


#endif


















