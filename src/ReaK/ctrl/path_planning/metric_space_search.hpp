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
     * \tparam Topology The topology type on which the points can reside, should model the MetricSpaceConcept.
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aSpace The topology on which the points reside.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename Topology, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, const Topology& aSpace, PositionMap aPosition) {
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
	  double dist = get(distance_metric, aSpace)(current_vp, get(aPosition, *(aBegin + j)), aSpace);
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
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aSpace The topology on which the points reside.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename Topology, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, const Topology& aSpace, PositionMap aPosition) {
      return aBegin + (get_global_rng()() % (aEnd - aBegin));
    };
};


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
};


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
    
    typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<tree_indexer>::edge_descriptor edge_type;
    typedef typename boost::graph_traits<tree_indexer>::out_edge_iterator out_edge_iter;
    typedef typename boost::graph_traits<tree_indexer>::in_edge_iterator in_edge_iter;
    
    typedef detail::compare_pair_first< distance_type, Key, std::less< distance_type > > priority_compare_type;
    typedef std::vector< std::pair< distance_type, Key > > priority_queue_type;
    
    tree_indexer m_tree;
    vertex_type m_root;
    typename boost::property_map< tree_indexer, Key vertex_properties::* >::type m_key;
    typename boost::property_map< tree_indexer, distance_type edge_properties::* >::type m_mu;
    typename PositionCachingPolicy::template effector< tree_indexer, point_type > m_caching_effector;
    
    const Topology& m_space;
    typename metric_space_traits<Topology>::distance_metric_type m_distance;
    PositionMap m_position;
    VPChooser m_vp_chooser;
    
    //non-copyable.
    dvp_tree(const dvp_tree<Key,Topology,PositionMap,Arity,VPChooser>&);
    dvp_tree<Key,Topology,PositionMap,Arity,VPChooser>& operator=(const dvp_tree<Key,Topology,PositionMap,Arity,VPChooser>&); 
    
    static bool closer(std::unordered_map<Key,distance_type>* m, const Key& k1, const Key& k2) {
      return (*m)[k1] < (*m)[k2];
    };
    
    /* NOTE Invalidates vertices */
    /* Does not require persistent vertices */
    void construct_node(vertex_type aNode, 
			typename std::vector<Key>::iterator aBegin, 
			typename std::vector<Key>::iterator aEnd, 
			std::unordered_map<Key,distance_type>& aDistMap) {
      typedef typename std::vector<Key>::iterator KeyIter;
      using std::swap;
      KeyIter vp_ind = m_vp_chooser(aBegin, aEnd, m_space, m_position);
      point_type vp_pt = get(m_position, *vp_ind);
      for(KeyIter it = aBegin; it != aEnd; ++it)
	aDistMap[*it] = m_distance(vp_pt, get(m_position, *it), m_space);
      swap(*vp_ind, *aBegin);
      //std::sort(aBegin,aEnd,boost::bind(closer,&aDistMap,_1,_2));
      m_caching_effector.put_vertex(aNode, *aBegin, vp_pt, m_key);
      aDistMap.erase(*aBegin);
      aBegin++;
      if( out_degree(aNode, m_tree) != 0 )
	clear_node(aNode);
      if((aEnd - aBegin) < static_cast<int>(Arity)) {
	for(KeyIter it = aBegin; it != aEnd; ++it) {
	  vertex_type k; edge_type e; 
	  boost::tie(k,e) = add_child_vertex(aNode, m_tree);
	  m_caching_effector.put_vertex(k, *it, m_position, m_key);
	  put(m_mu, e, aDistMap[*it]);
	};
      } else {
	for(unsigned int i=Arity;i>=1;--i) {
	  vertex_type k; edge_type e;
	  boost::tie(k,e) = add_child_vertex(aNode, m_tree);
	  int num_children = (aEnd - aBegin) / i;
	  std::nth_element(aBegin, aBegin + (num_children-1), aEnd, boost::bind(closer,&aDistMap,_1,_2));
	  put(m_mu, e, aDistMap[*(aBegin + (num_children-1))]);
	  KeyIter temp = aBegin; aBegin += num_children;
	  construct_node(k,temp,aBegin,aDistMap);
	};
      };
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    void find_nearest_impl(const point_type& aPoint, distance_type& aSigma, vertex_type aNode, 
			   priority_queue_type& aList, std::size_t K) const {
      typedef typename std::multimap<distance_type, Key>::value_type ListType;
      distance_type current_dist = m_distance(aPoint, m_caching_effector.get_position(aNode, m_position, m_key), m_space);
      if(current_dist < aSigma) { //is the vantage point within current search bound? Yes...
        //then add the vantage point to the NN list.
        aList.push_back(std::pair<distance_type, Key>(current_dist, get(m_key, aNode)));
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
	if(current_dist < get(m_mu, *ei)) 
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
	  while((ei_right != ei_end) && (get(m_mu,*ei_rightleft) < current_dist + aSigma)) {
	    find_nearest_impl(aPoint,aSigma,target(*ei_right,m_tree),aList,K);
	    ++ei_rightleft; ++ei_right;
	  };
	  break;
	} else if(right_stopped) {
	  out_edge_iter ei_leftleft = ei_left;
	  while((ei_left != ei) && (get(m_mu,*(--ei_leftleft)) > current_dist - aSigma)) {
	    find_nearest_impl(aPoint,aSigma,target(*ei_leftleft,m_tree),aList,K);
	    --ei_left;
	  };
	  break;
	} else {
	  out_edge_iter ei_leftleft = ei_left; --ei_leftleft;
	  distance_type d1 = get(m_mu,*ei_leftleft); //greater than 0 if ei_leftleft should be searched.
	  out_edge_iter ei_rightleft = ei_right; --ei_rightleft;
	  distance_type d2 = get(m_mu,*ei_rightleft); //less than 0 if ei_right should be searched.
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
    vertex_type get_leaf(const point_type& aPoint, vertex_type aNode) const {
      distance_type current_dist = m_distance(aPoint, m_caching_effector.get_position(aNode, m_position, m_key), m_space);
      //first, locate the partition in which aPoint is:
      if(out_degree(aNode,m_tree) == 0)
	return aNode;
      vertex_type result = aNode;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	result = target(*ei,m_tree);
	if(current_dist < get(m_mu, *ei)) 
	  break;
      };
      return get_leaf(aPoint,result);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    vertex_type get_key(Key aVertex, const point_type& aPoint, vertex_type aNode) const {
      if(get(m_key, aNode) == aVertex) return aNode;
      distance_type current_dist = m_distance(aPoint, m_caching_effector.get_position(aNode, m_position, m_key), m_space);
      //first, locate the partition in which aPoint is:
      if(out_degree(aNode,m_tree) == 0)
	throw int(0);
      vertex_type result = aNode;
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
	result = target(*ei,m_tree);
	if(current_dist < get(m_mu, *ei)) 
	  break;
      };
      return get_key(aVertex,aPoint,result);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    void update_mu_upwards(const point_type& aPoint, vertex_type aNode) {
      if(aNode == m_root) return;
      vertex_type parent = source(*(in_edges(aNode,m_tree).first), m_tree);
      distance_type dist = m_distance(aPoint, m_caching_effector.get_position(parent, m_position, m_key), m_space);
      if(dist > get(m_mu,*(in_edges(aNode,m_tree).first)))
	put(m_mu,*(in_edges(aNode,m_tree).first),dist);
      update_mu_upwards(aPoint,parent);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
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
    void collect_keys(std::vector<Key>& aList, vertex_type aNode) const {
      aList.push_back(get(m_key, aNode));
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei)
	collect_keys(aList,target(*ei,m_tree));
    };
    
    /* NOTE Invalidates vertices */
    /* NOTE Requires persistent vertices */
    void clear_node(vertex_type aNode) {
      if(out_degree(aNode,m_tree) == 0) return;
      std::vector<vertex_type> children;
      children.reserve(out_degree(aNode,m_tree));
      out_edge_iter ei,ei_end;
      for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei)
	children.push_back(target(*ei,m_tree));
      for(typename std::vector<vertex_type>::iterator it = children.begin(); it != children.end(); ++it)
	remove_branch(*it,m_tree);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
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
     * Construct the DVP-tree from a graph, topology and property-map.
     * \tparam Graph The graph type on which the vertices are taken from, should model the boost::VertexListGraphConcept.
     * \param g The graph from which to take the vertices.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aPosition The property-map that can be used to obtain the positions of the vertices.
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    template <typename Graph>
    dvp_tree(const Graph& g, 
	     const Topology& aSpace, 
	     PositionMap aPosition,
	     VPChooser aVPChooser = VPChooser()) : 
	     m_tree(), m_root(), 
             m_key(get(&vertex_properties::k,m_tree)), 
             m_mu(get(&edge_properties::d,m_tree)), 
             m_caching_effector(m_tree),
             m_space(aSpace), m_distance(get(distance_metric,aSpace)), 
             m_position(aPosition), m_vp_chooser(aVPChooser) {
      if(num_vertices(g) == 0) return;
      
      m_root = create_root(m_tree);
      typename boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
      boost::tie(vi,vi_end) = vertices(g);
      std::vector<Key> v(vi,vi_end); //Copy the list of vertices to random access memory.
      std::unordered_map<Key,distance_type> dist_map;
      construct_node(m_root, v.begin(), v.end(), dist_map);
    };
    
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
	     const Topology& aSpace, 
	     PositionMap aPosition,
	     VPChooser aVPChooser = VPChooser()) : 
	     m_tree(), m_root(), 
             m_key(get(&vertex_properties::k,m_tree)), 
             m_mu(get(&edge_properties::d,m_tree)), 
             m_caching_effector(m_tree),
             m_space(aSpace), m_distance(get(distance_metric,aSpace)), 
             m_position(aPosition), m_vp_chooser(aVPChooser) {
      if(aBegin == aEnd) return;
      
      m_root = create_root(m_tree);
      std::vector<Key> v(aBegin,aEnd); //Copy the list of vertices to random access memory.
      std::unordered_map<Key,distance_type> dist_map;
      construct_node(m_root, v.begin(), v.end(), dist_map);
    };
    
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
    
    std::size_t depth() const { return get_depth(m_root); };
    
    /**
     * Inserts a key-value (vertex).
     * \param u The vertex to be added to the DVP-tree.
     */
    void insert(Key u) { 
      if(num_vertices(m_tree) == 0) {
	m_root = create_root(m_tree); 
	m_caching_effector.put_vertex(m_root, u, m_position, m_key);
	return;
      };
      point_type u_pt = get(m_position, u); 
      vertex_type u_realleaf = get_leaf(u_pt,m_root);
      if(u_realleaf == m_root) { //if the root is the leaf, it requires special attention since no parent exists.
        std::vector<Key> key_list;       
	collect_keys(key_list,u_realleaf);   
	key_list.push_back(u); 
	clear_node(u_realleaf); 
	std::unordered_map<Key,distance_type> dist_map; 
	construct_node(u_realleaf, key_list.begin(), key_list.end(), dist_map); 
	update_mu_upwards(u_pt,u_realleaf); 
	return;
      };
      vertex_type u_leaf = source(*(in_edges(u_realleaf,m_tree).first),m_tree);
      if((out_degree(u_leaf,m_tree) < Arity) || (!is_leaf_node(u_leaf))) {
	// leaf node is not full of children, an additional child can be added 
	//  (must be reconstructed to keep ordering, but this is a trivial operation O(Arity)).
	//OR 
	// if leaf is not really a leaf, then it means that this sub-tree is definitely not balanced and not full either,
	//  then all the Keys ought to be collected and u_leaf ought to be reconstructed.
	std::vector<Key> key_list;       
	collect_keys(key_list,u_leaf);   
	key_list.push_back(u); 
	clear_node(u_leaf);              
	std::unordered_map<Key,distance_type> dist_map; 
	construct_node(u_leaf, key_list.begin(), key_list.end(), dist_map); 
	update_mu_upwards(u_pt,u_leaf);  
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
	  std::vector<Key> key_list;     
	  collect_keys(key_list,p);      
	  key_list.push_back(u);
	  clear_node(p);                 
	  std::unordered_map<Key,distance_type> dist_map;
	  construct_node(p, key_list.begin(), key_list.end(), dist_map);
	  update_mu_upwards(u_pt,p);     
	} else {
	  //this means that either the root node is full or there are branches of the tree that are deeper than u_realleaf, 
	  // and thus, in either case, u_realleaf should be expanded.
	  edge_type l_p;                 
	  boost::tie(p, l_p) = add_child_vertex(u_realleaf, m_tree);
	  m_caching_effector.put_vertex(p, u, u_pt, m_key);
	  put(m_mu, l_p, m_distance(u_pt, m_caching_effector.get_position(u_realleaf, m_position, m_key), m_space));
//	  put(m_mu, l_p, m_distance(u_pt, get(m_position,get(m_key,u_realleaf)), m_space));
	  update_mu_upwards(u_pt,u_realleaf);   
	};
      };
    };
    /**
     * Inserts a range of key-values (vertices).
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     */
    template <typename ForwardIterator>
    void insert(ForwardIterator aBegin, ForwardIterator aEnd) { 
      std::for_each(aBegin,aEnd,boost::bind(&dvp_tree<Key,Topology,PositionMap,Arity,VPChooser>::insert,this,_1));
      //TODO: There's got to be a better way to insert many elements (most likely a similar strategy to the erase multiple function).
    };
    
    /**
     * Erases the given vertex from the DVP-tree.
     * \param u The vertex to be removed from the DVP-tree.
     */
    void erase(Key u) { 
      if(num_vertices(m_tree) == 0) 
	return;
      point_type u_pt = get(m_position, u);
      vertex_type u_node;
      try {
        u_node = get_key(u, u_pt, m_root);
      } catch (int err) {
        return;
      };
      out_edge_iter ei, ei_end;
      std::vector<Key> key_list;
      if( out_degree(u_node, m_tree) > 0 ) {
        for(boost::tie(ei,ei_end) = out_edges(u_node,m_tree); ei != ei_end; ++ei)
	  collect_keys(key_list,target(*ei,m_tree));
      } else {
	if( u_node == m_root ) {
	  clear();
	  return;
	};
	vertex_type u_child = u_node;
	u_node = source(*(in_edges(u_node,m_tree).first), m_tree);
	for(boost::tie(ei,ei_end) = out_edges(u_node,m_tree); ei != ei_end; ++ei) {
	  if( target(*ei, m_tree) != u_child )
	    collect_keys(key_list,target(*ei, m_tree));
	};
      };
      clear_node(u_node);
      std::unordered_map<Key,distance_type> dist_map;
      construct_node(u_node, key_list.begin(), key_list.end(), dist_map);
    };
    
    /**
     * Erases the given vertex-range from the DVP-tree.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices to be erased.
     * \param aEnd The end of the range from which to take the vertices to be erased (one-past-last).
     */
    template <typename ForwardIterator>
    void erase(ForwardIterator aBegin, ForwardIterator aEnd) { 
      if(num_vertices(m_tree) == 0) return;
      typedef std::list< std::pair< std::pair< Key, Key>, std::vector<Key> > > key_listing;
      key_listing key_lists; //will hold a list of unique nodes and all their non-erased 
      
      for(ForwardIterator first = aBegin; first != aEnd; ++first) {
	// check if the key is already in the key-listings.
	bool already_collected = false;
	for(typename key_listing::iterator it = key_lists.begin(); it != key_lists.end(); ++it) {
	  typename std::vector<Key>::iterator it_key = std::binary_search(it->second.begin(), it->second.end(), *first);
	  if( it_key != it->second.end() ) {
	    // the key was found in this key-listing.
	    it->second.erase(it_key);
	    already_collected = true;
	    break;
	  };
	};
	if(already_collected) continue;
	vertex_type u_node;
	try {
	  u_node = get_key(*first, get(m_position, *first), m_root);
	} catch (int err) {
	  continue;
	};
	
	key_lists.push_back( std::make_pair(std::make_pair(*first, *first), std::vector<Key>()) );
	out_edge_iter ei,ei_end;
	if(out_degree(u_node, m_tree) == 0) {
	  if( u_node == m_root ) {
	    clear();
	    return;
	  };
	  vertex_type u_child = u_node;
	  u_node = source(*(in_edges(u_node,m_tree).first), m_tree);
	  for(boost::tie(ei,ei_end) = out_edges(u_node,m_tree); ei != ei_end; ++ei) {
	    if( target(*ei, m_tree) != u_child )
	      collect_keys(key_lists.back().second,target(*ei, m_tree));
	  };
	  key_lists.back().first.second = get(m_key, u_node, m_tree);
	} else {
	  for(boost::tie(ei,ei_end) = out_edges(u_node, m_tree); ei != ei_end; ++ei)
	    collect_keys(key_lists.back().second, target(*ei,m_tree));
	};
	std::sort(key_lists.back().second.begin(), key_lists.back().second.end());
      };
      
      for(typename key_listing::iterator it = key_lists.begin(); it != key_lists.end(); ++it) {
	typename key_listing::iterator it2 = it; ++it2;
	for(; it2 != key_lists.end(); ++it2) {
	  typename std::vector<Key>::iterator it_key = std::binary_search(it2->second.begin(), it2->second.end(), it->first.first);
	  if( it_key != it2->second.end() ) {
	    // the key was found in this key-listing.
	    it2->second.erase(it_key);
	    key_lists.erase(it--);
	    break;
	  };
	};
      };
      
      for(typename key_listing::iterator it = key_lists.begin(); it != key_lists.end(); ++it) {
	try {
	  vertex_type u_node = get_key(it->first.second, get(m_position, it->first.second), m_root);
	  clear_node(u_node);
	  std::unordered_map<Key,distance_type> dist_map;
	  construct_node(u_node,it->second.begin(),it->second.end(),dist_map);
	} catch (int err) { };
      };
    };
    
    /**
     * Clears the DVP-tree. 
     */
    void clear() {
      if( num_vertices(m_tree) == 0 ) {
        remove_branch(m_root,m_tree);
        m_root = vertex_type();
      };
    }; 
    
    /**
     * Finds the nearest neighbor to a given position.
     * \param aPoint The position from which to find the nearest-neighbor of.
     * \return The vertex in the DVP-tree that is closest to the given point.
     */
    Key find_nearest(const point_type& aPoint) const {
      if(num_vertices(m_tree) == 0) return Key();
      priority_queue_type Q;
      distance_type sig = std::numeric_limits<distance_type>::infinity();
      find_nearest_impl(aPoint,sig,m_root,Q,1);
      return Q.front().second;
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
      if(num_vertices(m_tree) == 0) return aOutputBegin;
      priority_queue_type Q;
      find_nearest_impl(aPoint,R,m_root,Q,K);
      std::sort_heap(Q.begin(), Q.end(), priority_compare_type());
      return detail::copy_neighbors_from_queue<Key, distance_type>(Q, aOutputBegin);
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
      if(num_vertices(m_tree) == 0) return aOutputBegin;
      priority_queue_type Q;
      find_nearest_impl(aPoint,R,m_root,Q,num_vertices(m_tree));
      std::sort_heap(Q.begin(), Q.end(), priority_compare_type());
      return detail::copy_neighbors_from_queue<Key, distance_type>(Q, aOutputBegin);
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


















