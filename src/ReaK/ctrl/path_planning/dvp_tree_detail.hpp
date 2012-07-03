/**
 * \file dvp_tree_detail.hpp
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

#ifndef REAK_DVP_TREE_DETAIL_HPP
#define REAK_DVP_TREE_DETAIL_HPP

#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include <unordered_map>
#include <vector>

#include "graph_alg/tree_concepts.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_raw_property_graph.hpp"

#include "metric_space_concept.hpp"
#include "topological_search.hpp"



namespace ReaK {

namespace pp {


/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space, with amortized O(logN) 
 * insertion-deletion. A DVP-tree is essentially a generalization of a search tree which only 
 * requires the space to have a metric which respects the triangular inequality.
 * \tparam TreeType The tree type to be used to store the entries of this DVP search tree.
 * \tparam Topology The topology type on which the points associated to each entry resides.
 * \tparam DistanceMetric The distance-metric type that can compute the distance between two points, should model DistanceMetricConcept over the given Topology. 
 * \tparam VertexKeyMap The property-map type that can map vertex properties of the tree to key-values (that identify the entries).
 * \tparam DistanceMap The property-map type that can map an edge property to its associated distance value (used internally).
 * \tparam PositionMap The property-map type that can map vertex properties of the tree to associated position values (positions in the topology).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename TreeType,
          typename Topology,
          typename DistanceMetric,
	  typename VertexKeyMap,
	  typename DistanceMap,
          typename PositionMap,
	  unsigned int Arity,
	  typename VPChooser>
class dvp_tree_impl
{
  public:
    
    /** Type of the points in the topology. */ 
    typedef typename topology_traits<Topology>::point_type point_type;
    /** Type of the point-differences in the topology. */ 
    typedef typename topology_traits<Topology>::point_difference_type point_difference_type;
    /** Type of the distance values. */ 
    typedef double distance_type;
    /** Type of the distance metric. */ 
    typedef DistanceMetric distance_metric;
    
    /** Type of the key-values that identify entries of the DVP tree. */ 
    typedef typename boost::property_traits< VertexKeyMap >::value_type key_type;
    
    typedef dvp_tree_impl<TreeType, Topology, DistanceMetric, VertexKeyMap, DistanceMap, PositionMap, Arity, VPChooser> self;
    
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<distance_metric, Topology>));
    
  private:
    
    typedef TreeType tree_indexer;
    
    typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<tree_indexer>::edge_descriptor edge_type;
    typedef typename boost::graph_traits<tree_indexer>::out_edge_iterator out_edge_iter;
    typedef typename boost::graph_traits<tree_indexer>::in_edge_iterator in_edge_iter;
    
    typedef typename tree_indexer::vertex_property_type vertex_property;
    typedef typename tree_indexer::edge_property_type edge_property;
    
    
    
    typedef detail::compare_pair_first< distance_type, vertex_type, std::less< distance_type > > priority_compare_type;
    typedef std::vector< std::pair< distance_type, vertex_type > > priority_queue_type;
    
    tree_indexer& m_tree;   ///< Tree storage.
    vertex_type m_root;    ///< Root node of the tree.
    
    VertexKeyMap m_key;      ///< A map from a vertex_property to a key_type value.
    DistanceMap m_mu;        ///< A map from an edge-property to a distance value.
    PositionMap m_position;  ///< A map from a vertex_property to a position value (should be Read-Write).
    
    ReaK::shared_ptr<const Topology> m_space;  ///< The pointer to the topology.
    distance_metric m_distance;  ///< The distance-metric functor.
    
    VPChooser m_vp_chooser;  ///< The vantage-point chooser (functor).
    
    //non-copyable.
    dvp_tree_impl(const self&);
    self& operator=(const self&); 
    
    /* Simple comparison function to sort by pre-computed distance values. */
    static bool closer(std::unordered_map<key_type,distance_type>* m, VertexKeyMap key, const vertex_property& k1, const vertex_property& k2) {
      return (*m)[get(key,k1)] < (*m)[get(key,k2)];
    };
    
    /* Simple predicate function to filter out invalid vertices (invalid key-value). */
    static bool is_vertex_prop_valid(VertexKeyMap key, const vertex_property& k1) {
      return (get(key,k1) != reinterpret_cast<key_type>(-1));
    };
    
    /* NOTE Invalidates vertices */
    /* Does not require persistent vertices */
    /* This is the main tree construction function. It takes the vertices in the iterator range and organizes them 
     * as a sub-tree below the aParentNode node (and aEdgeDist is the minimum distance to the parent of any node in the range). */
    void construct_node(vertex_type aParentNode, 
			double aEdgeDist,
			typename std::vector<vertex_property>::iterator aBegin, 
			typename std::vector<vertex_property>::iterator aEnd, 
			std::unordered_map<key_type,distance_type>& aDistMap) {
      typedef typename std::vector<vertex_property>::iterator PropIter;
      using std::swap;
      PropIter vp_ind = m_vp_chooser(aBegin, aEnd, *m_space, m_distance, m_position);
      point_type vp_pt = get(m_position, *vp_ind);
      for(PropIter it = aBegin; it != aEnd; ++it)
	aDistMap[get(m_key,*it)] = m_distance(vp_pt, get(m_position, *it), *m_space);
      swap(*vp_ind, *aBegin);
      
      vertex_type current_node;
      if( aParentNode != boost::graph_traits<tree_indexer>::null_vertex() ) {
        edge_property ep;
        put(m_mu, ep, aEdgeDist);
	edge_type e;
#ifdef RK_ENABLE_CXX0X_FEATURES
	boost::tie(current_node,e) = add_child_vertex(aParentNode, std::move(*aBegin), std::move(ep), m_tree);
#else
	boost::tie(current_node,e) = add_child_vertex(aParentNode, *aBegin, ep, m_tree);
#endif
      } else {
#ifdef RK_ENABLE_CXX0X_FEATURES
	current_node = create_root(std::move(*aBegin), m_tree);
#else
	current_node = create_root(*aBegin, m_tree);
#endif
	m_root = current_node;
      };
      aDistMap.erase(get(m_key,*aBegin));
      aBegin++;
      if((aEnd - aBegin) < static_cast<int>(Arity)) {
	std::sort(aBegin, aEnd, boost::bind(closer,&aDistMap,m_key,_1,_2));
	for(PropIter it = aBegin; it != aEnd; ++it) {
	  edge_property ep;
	  put(m_mu, ep, aDistMap[get(m_key,*it)]); 
#ifdef RK_ENABLE_CXX0X_FEATURES
	  add_child_vertex(current_node, std::move(*it), std::move(ep), m_tree);
#else
	  add_child_vertex(current_node, *it, ep, m_tree);
#endif
	};
      } else {
	for(unsigned int i = Arity; i >= 1; --i) {
	  int num_children = (aEnd - aBegin) / i;
	  std::nth_element(aBegin, aBegin + (num_children-1), aEnd, boost::bind(closer,&aDistMap,m_key,_1,_2));
	  PropIter temp = aBegin; aBegin += num_children;
	  construct_node(current_node, aDistMap[get(m_key,*(aBegin + (num_children-1)))], 
			 temp, aBegin, aDistMap);
	};
      };
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This is the main nearest-neighbor query function. This takes a query point, a maximum 
     * neighborhood radius (aSigma), aNode to start recursing from, the current max-heap of neighbors,
     * and the maximum number of neighbors. This function can be used for any kind of NN query (single, kNN, or ranged). */
    void find_nearest_impl(const point_type& aPoint, distance_type& aSigma, vertex_type aNode, 
			   priority_queue_type& aList, std::size_t K) const {
      distance_type current_dist = m_distance(aPoint, get(m_position, get(boost::vertex_raw_property,m_tree,aNode)), *m_space);
      if(current_dist < aSigma) { //is the vantage point within current search bound? Yes...
        //then add the vantage point to the NN list.
        aList.push_back(std::pair<distance_type, vertex_type>(current_dist, aNode));
	std::push_heap(aList.begin(), aList.end(), priority_compare_type());
        if(aList.size() > K) { //are there too many nearest neighbors? Yes... 
          std::pop_heap(aList.begin(), aList.end(), priority_compare_type());
	  aList.pop_back(); //delete last element to keep aList with K elements
	  aSigma = aList.front().first; //distance of the last element is now the search bound aSigma.
	};
      };
      out_edge_iter ei,ei_end;
      //first, locate the partition in which aPoint is:
      if(out_degree(aNode,m_tree) == 0)
	return;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	if(current_dist < get(m_mu, get(boost::edge_raw_property,m_tree,*ei))) 
	  break;
      };
      if(ei == ei_end) --ei; //back-track if the end was reached.
      find_nearest_impl(aPoint,aSigma,target(*ei,m_tree),aList,K); //search in the most likely node.
      out_edge_iter ei_left = ei;
      out_edge_iter ei_right = ei; ++ei_right;
      boost::tie(ei,ei_end) = out_edges(aNode,m_tree); //find the bounds again (start and end).
      bool left_stopped  = (ei_left == ei);
      bool right_stopped = (ei_right == ei_end);
      while(true) {
        if(left_stopped) {
	  out_edge_iter ei_rightleft = ei_right; --ei_rightleft;
	  while((ei_right != ei_end) && (get(m_mu,get(boost::edge_raw_property,m_tree,*ei_rightleft)) < current_dist + aSigma)) {
	    find_nearest_impl(aPoint,aSigma,target(*ei_right,m_tree),aList,K);
	    ++ei_rightleft; ++ei_right;
	  };
	  break;
	} else if(right_stopped) {
	  out_edge_iter ei_leftleft = ei_left;
	  while((ei_left != ei) && (get(m_mu,get(boost::edge_raw_property,m_tree,*(--ei_leftleft))) > current_dist - aSigma)) {
	    find_nearest_impl(aPoint,aSigma,target(*ei_leftleft,m_tree),aList,K);
	    --ei_left;
	  };
	  break;
	} else {
	  out_edge_iter ei_leftleft = ei_left; --ei_leftleft;
	  distance_type d1 = get(m_mu,get(boost::edge_raw_property,m_tree,*ei_leftleft)); //greater than 0 if ei_leftleft should be searched.
	  out_edge_iter ei_rightleft = ei_right; --ei_rightleft;
	  distance_type d2 = get(m_mu,get(boost::edge_raw_property,m_tree,*ei_rightleft)); //less than 0 if ei_right should be searched.
	  if(d1 + d2 > 2.0 * current_dist) { //this means that ei_leftleft's boundary is closer to aPoint.
            if(d1 + aSigma - current_dist > 0) {
              find_nearest_impl(aPoint,aSigma,target(*ei_leftleft,m_tree),aList,K);
	      ei_left = ei_leftleft;
	      if(d2 - aSigma - current_dist < 0) {
	        find_nearest_impl(aPoint,aSigma,target(*ei_right,m_tree),aList,K);
	        ++ei_right;
	      } else
		right_stopped = true;
	    } else
	      break;
	  } else {
	    if(d2 - aSigma - current_dist < 0) {
	      find_nearest_impl(aPoint,aSigma,target(*ei_right,m_tree),aList,K);
	      ++ei_right;
	      if(d1 + aSigma - current_dist > 0) {
	        find_nearest_impl(aPoint,aSigma,target(*ei_leftleft,m_tree),aList,K);
	        ei_left = ei_leftleft;
	      } else 
		left_stopped = true;
	    } else
	      break;
	  };
	};
        left_stopped  = (ei_left == ei);
        right_stopped = (ei_right == ei_end);
      };
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function does a single nearest-neighbor query to find the tree-leaf that is closest to 
     * the given point (starts to recurse from aNode). */
    vertex_type get_leaf(const point_type& aPoint, vertex_type aNode) const {
      distance_type current_dist = m_distance(aPoint, get(m_position, get(boost::vertex_raw_property,m_tree,aNode)), *m_space);
      //first, locate the partition in which aPoint is:
      if(out_degree(aNode,m_tree) == 0)
	return aNode;
      vertex_type result = aNode;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	result = target(*ei,m_tree);
	if(current_dist < get(m_mu, get(boost::edge_raw_property,m_tree,*ei))) 
	  break;
      };
      return get_leaf(aPoint,result);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function finds the tree-vertex with the given key-value and position. */
    vertex_type get_vertex(key_type aKey, const point_type& aPoint, vertex_type aNode) const {
      if(get(m_key, get(boost::vertex_raw_property,m_tree,aNode)) == aKey) 
	return aNode;
      distance_type current_dist = m_distance(aPoint, get(m_position, get(boost::vertex_raw_property,m_tree,aNode)), *m_space);
      //first, locate the partition in which aPoint is:
      if(out_degree(aNode,m_tree) == 0)
	throw int(0);
      vertex_type result = aNode;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	result = target(*ei,m_tree);
	if(current_dist < get(m_mu, get(boost::edge_raw_property,m_tree,*ei))) 
	  break;
      };
      return get_vertex(aKey,aPoint,result);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function updates the edge distance values as a consequence of a new point being added. */
    void update_mu_upwards(const point_type& aPoint, vertex_type aNode) {
      if(aNode == m_root) return;
      vertex_type parent = source(*(in_edges(aNode,m_tree).first), m_tree);
      distance_type dist = m_distance(aPoint, get(m_position, get(boost::vertex_raw_property,m_tree,parent)), *m_space);
      if(dist > get(m_mu, get(boost::edge_raw_property, m_tree, *(in_edges(aNode,m_tree).first))))
	put(m_mu, get(boost::edge_raw_property, m_tree, *(in_edges(aNode,m_tree).first)), dist);
      update_mu_upwards(aPoint, parent);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function determines if a given node has no children or if all its children have no children. */
    bool is_leaf_node(vertex_type aNode) const {
      if(out_degree(aNode,m_tree) == 0) return true;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	if(out_degree(target(*ei,m_tree),m_tree) != 0)
	  return false;
      };
      return true;
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function determines if a given node is the root of a balanced (and full) sub-tree. */
    bool is_node_full(vertex_type aNode, int& depth_limit) const {
      if(depth_limit < 0)
	return false;
      if((out_degree(aNode,m_tree) == 0) && (depth_limit == 0))
	return true;
      --depth_limit;
      if((out_degree(aNode,m_tree) == 0) || (out_degree(aNode,m_tree) < Arity)) 
	return false;
      out_edge_iter ei,ei_end;
      if(is_leaf_node(aNode)) {
	if(depth_limit == 0)
	  return true;
	else
	  return false;
      };
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	int new_depth_limit = depth_limit;
	if(!is_node_full(target(*ei,m_tree),new_depth_limit)) { 
	  depth_limit = new_depth_limit;
	  return false;
	};
      };
      return true;
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function collects the list of vertices that are in the sub-tree rooted at the given node (exclusively). */
    void collect_vertices(std::vector<vertex_type>& aList, vertex_type aNode) const {
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	aList.push_back(target(*ei,m_tree));
	collect_vertices(aList,target(*ei,m_tree));
      };
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function computes the maximum depth of the tree rooted at the given node (note: this is an 
     * expensive operation, but is useful when debugging and testing). */
    std::size_t get_depth(vertex_type aNode) const {
      std::size_t max_depth = 0;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	std::size_t temp = get_depth(target(*ei,m_tree));
	if(temp > max_depth)
	  max_depth = temp;
      };
      return max_depth + 1;
    };
    
  public:
    
    /**
     * Construct the DVP-tree from a graph, position-map, topology, etc..
     * \tparam Graph The graph type on which the vertices are taken from, should model the boost::VertexListGraphConcept.
     * \tparam GraphPositionMap The property-map that associates position values to nodes of the graph.
     * \param g The graph from which to take the vertices.
     * \param aGraphPosition The property-map that takes a node of the graph and produces (or looks up) a position value.
     * \param aTree The tree object that will be used to store the DVP structure.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aDistance The distance metric to use to compute distances between points of the topology.
     * \param aKey The key-map to use to obtain and store the key values for a given vertex-property object.
     * \param aMu The property-map which associates a distance-value to each edge of the tree.
     * \param aPosition The property-map that can be used to obtain and store the positions of the vertices (vertex-property objects).
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    template <typename Graph, typename GraphPositionMap>
    dvp_tree_impl(const Graph& g, 
		  GraphPositionMap aGraphPosition,
		  tree_indexer& aTree,
		  const ReaK::shared_ptr<const Topology>& aSpace, 
		  distance_metric aDistance,
		  VertexKeyMap aKey,
		  DistanceMap aMu,
		  PositionMap aPosition,
		  VPChooser aVPChooser) : 
		  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
		  m_key(aKey), m_mu(aMu), m_position(aPosition),
		  m_space(aSpace), m_distance(aDistance), 
		  m_vp_chooser(aVPChooser) {
      
      if(num_vertices(g) == 0) return;
      
      typename boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
      boost::tie(vi,vi_end) = vertices(g);
      std::vector<vertex_property> v_bin; //Copy the list of vertices to random access memory.
      v_bin.reserve(num_vertices(g));
      for(;vi != vi_end; ++vi) {
	vertex_property vp;
	put(m_key, vp, *vi);
	put(m_position, vp, get(aGraphPosition, *vi));
#ifdef RK_ENABLE_CXX0X_FEATURES
	v_bin.push_back(std::move(vp));
#else
	v_bin.push_back(vp);
#endif
      };
      
      std::unordered_map<key_type, distance_type> dist_map;
      construct_node(boost::graph_traits<tree_indexer>::null_vertex(), 0.0, 
		     v_bin.begin(), v_bin.end(), dist_map);
    };
    
    
    /**
     * Construct the DVP-tree from a range, position-map, topology, etc..
     * \tparam ForwardIterator The forward-iterator type from which the vertices can be obtained.
     * \tparam ElemPositionMap The property-map that associates position values to vertices in the given range.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     * \param aElemPosition The property-map that takes a node in the given range and produces (or looks up) a position value.
     * \param aTree The tree object that will be used to store the DVP structure.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aDistance The distance metric to use to compute distances between points of the topology.
     * \param aKey The key-map to use to obtain and store the key values for a given vertex-property object.
     * \param aMu The property-map which associates a distance-value to each edge of the tree.
     * \param aPosition The property-map that can be used to obtain and store the positions of the vertices (vertex-property objects).
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    template <typename ForwardIterator, typename ElemPositionMap>
    dvp_tree_impl(ForwardIterator aBegin,
	          ForwardIterator aEnd,
		  ElemPositionMap aElemPosition,
		  tree_indexer& aTree,
		  const ReaK::shared_ptr<const Topology>& aSpace, 
		  distance_metric aDistance,
		  VertexKeyMap aKey,
		  DistanceMap aMu,
		  PositionMap aPosition,
		  VPChooser aVPChooser) : 
		  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
		  m_key(aKey), m_mu(aMu), m_position(aPosition),
		  m_space(aSpace), m_distance(aDistance), 
		  m_vp_chooser(aVPChooser) {
      if(aBegin == aEnd) return;
      
      std::vector<vertex_property> v_bin; //Copy the list of vertices to random access memory.
      for(;aBegin != aEnd; ++aBegin) {
	vertex_property vp;
	put(m_key, vp, *aBegin);
	put(m_position, vp, get(aElemPosition, *aBegin));
#ifdef RK_ENABLE_CXX0X_FEATURES
	v_bin.push_back(std::move(vp));
#else
	v_bin.push_back(vp);
#endif
      };
      
      std::unordered_map<key_type, distance_type> dist_map;
      construct_node(boost::graph_traits<tree_indexer>::null_vertex(), 0.0, 
		     v_bin.begin(), v_bin.end(), dist_map);
    };
    
    /**
     * Construct an empty DVP-tree from a topology, etc..
     * \param aTree The tree object that will be used to store the DVP structure.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aDistance The distance metric to use to compute distances between points of the topology.
     * \param aKey The key-map to use to obtain and store the key values for a given vertex-property object.
     * \param aMu The property-map which associates a distance-value to each edge of the tree.
     * \param aPosition The property-map that can be used to obtain and store the positions of the vertices (vertex-property objects).
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    dvp_tree_impl(tree_indexer& aTree,
		  const ReaK::shared_ptr<const Topology>& aSpace, 
		  distance_metric aDistance,
		  VertexKeyMap aKey,
		  DistanceMap aMu,
		  PositionMap aPosition,
		  VPChooser aVPChooser) : 
		  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
		  m_key(aKey), m_mu(aMu), m_position(aPosition),
		  m_space(aSpace), m_distance(aDistance), 
		  m_vp_chooser(aVPChooser) { };
    
    /**
     * Checks if the DVP-tree is empty.
     * \return True if the DVP-tree is empty.
     */
    bool empty() const { return (num_vertices(m_tree) == 0); };
    /**
     * Returns the size of the DVP-tree (the number of vertices it contains.
     * \return The size of the DVP-tree (the number of vertices it contains.
     */
    std::size_t size() const { return num_vertices(m_tree); };
    
    /**
     * Returns the depth of the tree.
     * \note This operation must recurse through all the branches of the tree (depth-first), and is 
     * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack) 
     * w.r.t. the depth of tree).
     * \return The depth of the tree.
     */
    std::size_t depth() const { return get_depth(m_root); };
    
    /**
     * Inserts a vertex into the tree.
     * \param up The vertex-property to be added to the DVP-tree.
     */
#ifdef RK_ENABLE_CXX0X_FEATURES
    void insert(vertex_property up) {
#else
    void insert(const vertex_property& up) {
#endif
      if(num_vertices(m_tree) == 0) {
#ifdef RK_ENABLE_CXX0X_FEATURES
	m_root = create_root(std::move(up), m_tree); 
#else
	m_root = create_root(up, m_tree); 
#endif
	return;
      };
      point_type u_pt = get(m_position, up); 
      vertex_type u_realleaf = get_leaf(u_pt,m_root);
      if(u_realleaf == m_root) { //if the root is the leaf, it requires special attention since no parent exists.
        std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
	prop_list.push_back(std::move(up));
#else
	prop_list.push_back(up);
#endif
	remove_branch(u_realleaf, back_inserter(prop_list), m_tree);
	m_root = boost::graph_traits<tree_indexer>::null_vertex();
	u_realleaf = m_root;
	std::unordered_map<key_type, distance_type> dist_map;
	construct_node(u_realleaf, 0.0, prop_list.begin(), prop_list.end(), dist_map); 
	return;
      };
      vertex_type u_leaf = source(*(in_edges(u_realleaf,m_tree).first),m_tree);
      if((out_degree(u_leaf,m_tree) < Arity) || (!is_leaf_node(u_leaf))) {
	// leaf node is not full of children, an additional child can be added 
	//  (must be reconstructed to keep ordering, but this is a trivial operation O(Arity)).
	//OR 
	// if leaf is not really a leaf, then it means that this sub-tree is definitely not balanced and not full either,
	//  then all the Keys ought to be collected and u_leaf ought to be reconstructed.
	update_mu_upwards(u_pt,u_leaf);
	distance_type e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(u_leaf,m_tree).first)));
	vertex_type u_leaf_parent;
	if(u_leaf != m_root)
	  u_leaf_parent = source(*(in_edges(u_leaf,m_tree).first),m_tree);
	else
	  u_leaf_parent = boost::graph_traits<tree_indexer>::null_vertex();
	std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
	prop_list.push_back(std::move(up));
#else
	prop_list.push_back(up);
#endif
	remove_branch(u_leaf, back_inserter(prop_list), m_tree);
	std::unordered_map<key_type,distance_type> dist_map; 
	construct_node(u_leaf_parent, e_dist, prop_list.begin(), prop_list.end(), dist_map); 
      } else {
	//if it is a full-leaf, then this is a leaf node, and it is balanced but full, 
	// we should then find a non-full parent.
	vertex_type p = u_leaf;   
	int actual_depth_limit = 1;
	int last_depth_limit = actual_depth_limit;
  	while((p != m_root) && (is_node_full(p,last_depth_limit))) {
	  p = source(*(in_edges(p,m_tree).first),m_tree);
	  last_depth_limit = ++actual_depth_limit;
	};
	bool is_p_full = false; 
	if(p == m_root)
	  is_p_full = is_node_full(p,last_depth_limit);
	if((!is_p_full) && (last_depth_limit >= 0)) {
	  //this means that we can add our key to the sub-tree of p and reconstruct from there.
	  update_mu_upwards(u_pt,p);
	  distance_type e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(p,m_tree).first)));
	  vertex_type p_parent;
          if(p != m_root)
	    p_parent = source(*(in_edges(p,m_tree).first),m_tree);
	  else
	    p_parent = boost::graph_traits<tree_indexer>::null_vertex();
	  
	  std::vector<vertex_property> prop_list;
#ifdef RK_ENABLE_CXX0X_FEATURES
	  prop_list.push_back(std::move(up));
#else
	  prop_list.push_back(up);
#endif
	  remove_branch(p, back_inserter(prop_list), m_tree);
	  std::unordered_map<key_type,distance_type> dist_map;
	  construct_node(p_parent, e_dist, prop_list.begin(), prop_list.end(), dist_map);
	} else {
	  //this means that either the root node is full or there are branches of the tree that are deeper than u_realleaf, 
	  // and thus, in either case, u_realleaf should be expanded.
	  edge_type l_p;
	  edge_property ep;
	  put(m_mu, ep, m_distance(u_pt, get(m_position, get(boost::vertex_raw_property,m_tree,u_realleaf)), *m_space));
#ifdef RK_ENABLE_CXX0X_FEATURES
	  boost::tie(p, l_p) = add_child_vertex(u_realleaf, std::move(up), std::move(ep), m_tree);
#else
	  boost::tie(p, l_p) = add_child_vertex(u_realleaf, up, ep, m_tree);
#endif
	  update_mu_upwards(u_pt,u_realleaf);
	};
      };
    };
    /**
     * Inserts a range of vertices.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     */
    template <typename ForwardIterator>
    void insert(ForwardIterator aBegin, ForwardIterator aEnd) { 
      std::for_each(aBegin,aEnd,boost::bind(&self::insert,this,_1));
      //TODO: There's got to be a better way to insert many elements (most likely a similar strategy to the erase multiple function).
    };
    
    
    /**
     * Erases the given vertex from the DVP-tree.
     * \param u_node The vertex to be removed from the DVP-tree.
     */
    void erase(vertex_type u_node) { 
      if(num_vertices(m_tree) == 0) 
	return;
      if( (u_node == m_root) && (num_vertices(m_tree) == 1) ) {
	clear();
	return;
      };
      distance_type e_dist = 0.0;
      vertex_type u_parent = boost::graph_traits<tree_indexer>::null_vertex();
      if(u_node != m_root) {
        e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(u_node,m_tree).first)));
	u_parent = source(*(in_edges(u_node,m_tree).first), m_tree);
      };
      
      out_edge_iter ei, ei_end;
      std::vector<vertex_property> prop_list;
      if( out_degree(u_node, m_tree) > 0 ) {
	remove_branch(u_node, back_inserter(prop_list), m_tree);
      } else {
	remove_branch(u_node, back_inserter(prop_list), m_tree);
	u_node = u_parent;
	u_parent = source(*(in_edges(u_node,m_tree).first), m_tree);
	remove_branch(u_node, back_inserter(prop_list), m_tree);
      };
      std::unordered_map<key_type,distance_type> dist_map;
      construct_node(u_parent, e_dist, prop_list.begin() + 1 /* skip first node (u_node) */, prop_list.end(), dist_map);
    };
    
    /**
     * Looks up the given key-value and position from the DVP-tree.
     * \param u_key The key-value to be removed from the DVP-tree.
     * \param u_pt The position-value corresponding to the key-value.
     * \return The vertex descriptor (into the tree-storage) for the given key value.
     */
    vertex_type get_vertex(key_type u_key, const point_type& u_pt) const {
      try {
        return get_vertex(u_key, u_pt, m_root);
      } catch (int err) {
        return boost::graph_traits<tree_indexer>::null_vertex();
      };
    };
    
    /**
     * Erases the given key-value and position from the DVP-tree.
     * Note that this function is not recommended if the vertex descriptor is known, as it will require a key-lookup.
     * \param u_key The key-value to be removed from the DVP-tree.
     * \param u_pt The position-value corresponding to the key-value.
     */
    void erase(key_type u_key, const point_type& u_pt) {
      vertex_type u_node;
      try {
        u_node = get_vertex(u_key, u_pt, m_root);
      } catch (int err) {
        return;
      };
      erase(u_node);
    };
    
    /**
     * Erases the given vertex-range from the DVP-tree.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices (by tree vertex descriptors).
     * \param aBegin The start of the range from which to take the vertices to be erased.
     * \param aEnd The end of the range from which to take the vertices to be erased (one-past-last).
     */
    template <typename ForwardIterator>
    void erase(ForwardIterator aBegin, ForwardIterator aEnd) { 
      if(num_vertices(m_tree) == 0) return;
      
      typedef std::list< std::pair< vertex_type,   // the node at the trunk of the sub-tree to be re-balanced.
                                    std::vector<vertex_type>  // the list of nodes below the re-balanced trunk.
				  > > vertex_listing;
      vertex_listing v_lists; //will hold a list of unique nodes and all their non-erased 
      
      // First, generate the vertex-listings in preparation for the deletion.
      for(ForwardIterator first = aBegin; first != aEnd; ++first) {
	vertex_type removal_trunk = *first;
	if(out_degree(*first, m_tree) == 0)
	  removal_trunk = source(*(in_edges(*first,m_tree).first), m_tree);
	put(m_key, get(boost::vertex_raw_property,m_tree,*first), reinterpret_cast<key_type>(-1)); // mark as invalid, for deletion.
	
	bool already_collected = false;
	for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ++it) {
	  // is removal_trunk contained in the *it listing?
	  typename std::vector<vertex_type>::iterator it_key = std::binary_search(it->second.begin(), it->second.end(), removal_trunk);
	  if( it_key != it->second.end() ) {
	    already_collected = true;
	    break;
	  };
	};
	
	if(!already_collected) {
	  std::vector<vertex_type> removal_list;
	  removal_list.push_back(removal_trunk);
	  collect_vertices(removal_list, removal_trunk);
          std::sort(removal_list.begin(), removal_list.end());
	  
	  for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ) {
	    // is the *it trunk contained in the removal_trunk?
	    typename std::vector<vertex_type>::iterator it_key = std::binary_search(removal_list.begin(), removal_list.end(), it->first);
	    if( it_key != removal_list.end() )
	      it = v_lists.erase(it);
	    else
	      ++it;
	  };
	  
	  v_lists.push_back( std::make_pair(removal_trunk, std::vector<vertex_type>()) );
	  v_lists.back().second.swap(removal_list);
	};
	
      };
      
      for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ++it) {
	
        distance_type e_dist = 0.0;
        vertex_type u_parent = boost::graph_traits<tree_indexer>::null_vertex();
        if(it->first != m_root) {
          e_dist = get(m_mu, get(boost::edge_raw_property,m_tree,*(in_edges(it->first,m_tree).first)));
	  u_parent = source(*(in_edges(it->first, m_tree).first), m_tree);
        };
      
        out_edge_iter ei, ei_end;
        std::vector<vertex_property> prop_list;
	remove_branch(it->first, back_inserter(prop_list), m_tree);
	prop_list.erase( remove_if(prop_list.begin(), prop_list.end(), boost::bind(is_vertex_prop_valid, m_key, _1)), prop_list.end());
        std::unordered_map<key_type,distance_type> dist_map;
        construct_node(u_parent, e_dist, prop_list.begin(), prop_list.end(), dist_map);
      };      
    };
    
    /**
     * Clears the DVP-tree. 
     */
    void clear() {
      if( num_vertices(m_tree) == 0 ) {
        remove_branch(m_root,m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
      };
    }; 
    
    /**
     * Finds the nearest neighbor to a given position.
     * \param aPoint The position from which to find the nearest-neighbor of.
     * \return The vertex in the DVP-tree that is closest to the given point.
     */
    vertex_type find_nearest(const point_type& aPoint) const {
      if(num_vertices(m_tree) == 0) 
	return boost::graph_traits<tree_indexer>::null_vertex();
      priority_queue_type Q;
      distance_type sig = std::numeric_limits<distance_type>::infinity();
      find_nearest_impl(aPoint,sig,m_root,Q,1);
      return Q.front().second;
    };
    
    /**
     * Finds the K nearest-neighbors to a given position.
     * \tparam OutputIterator The forward- output-iterator type which can contain the 
     *         list of nearest-neighbors by tree vertex descriptors.
     * \param aPoint The position from which to find the nearest-neighbors.
     * \param aOutputBegin An iterator to the first place where to put the sorted list of 
     *        elements by tree vertex descriptors with the smallest distance.
     * \param K The number of nearest-neighbors.
     * \param R The maximum distance value for the nearest-neighbors.
     * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
     */
    template <typename OutputIterator>
    OutputIterator find_nearest(const point_type& aPoint, OutputIterator aOutputBegin, std::size_t K, distance_type R = std::numeric_limits<distance_type>::infinity()) const {
      if(num_vertices(m_tree) == 0) 
	return aOutputBegin;
      priority_queue_type Q;
      find_nearest_impl(aPoint,R,m_root,Q,K);
      std::sort_heap(Q.begin(), Q.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = Q.begin(); it != Q.end(); ++it)
	*(aOutputBegin++) = it->second;
      return aOutputBegin;
    };
    
    /**
     * Finds the nearest-neighbors to a given position within a given range (radius).
     * \tparam OutputIterator The forward- output-iterator type which can contain the 
     *         list of nearest-neighbors by tree vertex descriptors.
     * \param aPoint The position from which to find the nearest-neighbors.
     * \param aOutputBegin An iterator to the first place where to put the sorted list of 
     *        elements by tree vertex descriptors with the smallest distance.
     * \param R The maximum distance value for the nearest-neighbors.
     * \return The output-iterator to the end of the list of nearest neighbors (starting from "output_first").
     */
    template <typename OutputIterator>
    OutputIterator find_in_range(const point_type& aPoint, OutputIterator aOutputBegin, distance_type R) const {
      if(num_vertices(m_tree) == 0) 
	return aOutputBegin;
      priority_queue_type Q;
      find_nearest_impl(aPoint,R,m_root,Q,num_vertices(m_tree));
      std::sort_heap(Q.begin(), Q.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = Q.begin(); it != Q.end(); ++it)
	*(aOutputBegin++) = it->second;
      return aOutputBegin;
    };
    
    
    
    struct mutation_visitor {
      self* m_parent;
      
      mutation_visitor(self* aParent) : m_parent(aParent) { };
      
      void remove_vertex(vertex_type v, tree_indexer&) const {
	m_parent->erase(v);
      };
      
      void add_vertex(const vertex_property& vp, tree_indexer&) const {
	m_parent->insert(vp);
      };
      
#ifdef RK_ENABLE_CXX0X_FEATURES
      void add_vertex(vertex_property&& vp, tree_indexer&) const {
	m_parent->insert(std::move(vp));
      };
#endif
    };
    
};





};

};


#endif


















