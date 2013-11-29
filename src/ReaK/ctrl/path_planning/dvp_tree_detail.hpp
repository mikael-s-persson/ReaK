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

#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include "metric_space_concept.hpp"
#include "proper_metric_concept.hpp"

#include "global_rng.hpp"

// BGL-Extra includes:
#include <boost/graph/tree_traits.hpp>
#include <boost/graph/tree_adaptor.hpp>

// Pending inclusion in BGL-Extra:
#include "graph_alg/bgl_raw_property_graph.hpp"


#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <vector>
#include <stack>
#include <queue>
#include <cmath>
#include <utility>
#include <functional>
#include <algorithm>


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
     * \tparam DistanceMetric The distance-metric type over the topology type.
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aDistance The distance-value functor.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename DistanceMetric, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, DistanceMetric aDistance, PositionMap aPosition) const {
      typedef typename boost::property_traits<PositionMap>::value_type Point;
      RandomAccessIter best_pt = aEnd;
      double best_dev = -1;
      for(unsigned int i=0; i < (aEnd - aBegin) / m_divider + 1;++i) {
        RandomAccessIter current_pt = aBegin + (get_global_rng()() % (aEnd - aBegin));
        double current_mean = 0.0;
        double current_dev = 0.0;
        Point current_vp = get(aPosition, *current_pt);
        for(unsigned int j=0; aBegin + j != aEnd; ++j) {
          double dist = aDistance(current_vp, get(aPosition, *(aBegin + j)));
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
     * \tparam DistanceMetric The distance-metric type over the topology type.
     * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the iterators) to a point (position).
     * \param aBegin The start of the range of vertices.
     * \param aEnd The end of the range of vertices (one element past the end).
     * \param aDistance The distance-metric over the given topology.
     * \param aPosition The property-map used to obtain the positions from the vertices.
     * \return A random-access iterator to the chosen vantage-point.
     */
    template <typename RandomAccessIter, typename DistanceMetric, typename PositionMap>
    RandomAccessIter operator() (RandomAccessIter aBegin, RandomAccessIter aEnd, DistanceMetric aDistance, PositionMap aPosition) const { RK_UNUSED(aDistance); RK_UNUSED(aPosition); 
      return aBegin + (get_global_rng()() % (aEnd - aBegin));
    };
};


namespace detail {
  
  struct dvp_tree_key_hasher {
    std::size_t operator()(std::size_t d) const {
      ::boost::hash< std::size_t > hasher;
      return hasher(d);
    };
    template <typename Iter>
    std::size_t operator()(Iter it) const {
      typedef typename std::iterator_traits<Iter>::pointer PtrType;
      ::boost::hash< PtrType > hasher;
      return hasher(&(*it));
    };
  };
  
};



/**
 * This class implements a Dynamic Vantage-Point Tree (DVP-Tree) that
 * allows for O(logN) time nearest-neighbor queries in a metric-space, with amortized O(logN) 
 * insertion-deletion. A DVP-tree is essentially a generalization of a search tree which only 
 * requires the space to have a metric which respects the triangular inequality.
 * \tparam TreeType The tree type to be used to store the entries of this DVP search tree.
 * \tparam Topology The topology type on which the points associated to each entry resides.
 * \tparam VertexKeyMap The property-map type that can map vertex properties of the tree to key-values (that identify the entries).
 * \tparam DistanceMap The property-map type that can map an edge property to its associated distance value (used internally).
 * \tparam PositionMap The property-map type that can map vertex properties of the tree to associated position values (positions in the topology).
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VPChooser The functor type to use to choose the vantage-point out of a set of vertices.
 */
template <typename TreeType,
          typename Topology,
          typename VertexKeyMap,
          typename DistanceMap,
          typename PositionMap,
          unsigned int Arity,
          typename VPChooser>
class dvp_tree_impl
{
  public:
    
    typedef dvp_tree_impl<TreeType, Topology, VertexKeyMap, DistanceMap, PositionMap, Arity, VPChooser> self;
    
    /** Type of the points in the topology. */ 
    typedef typename boost::property_traits<PositionMap>::value_type point_type;
    /** Type of the distance values. */ 
    typedef double distance_type;
    
    typedef typename metric_space_traits<Topology>::distance_metric_type distance_metric_type;
    typedef typename get_proper_metric<Topology>::type proper_metric_type;
    
    struct dist_metric_only_impl {
      shared_ptr<const Topology> p_space;
      distance_metric_type m_distance;
      
      double distance(const point_type& a, const point_type& b) const {
        return m_distance(a, b, *p_space);
      };
      double proper_distance(const point_type& a, const point_type& b) const {
        return m_distance(a, b, *p_space);
      };
      double operator()(const point_type& a, const point_type& b) const {
        return proper_distance(a, b, *p_space);
      };
      
      dist_metric_only_impl(const shared_ptr<const Topology>& aSpace) : p_space(aSpace), m_distance(get(distance_metric, *aSpace)) { };
    };
    
    struct dist_metric_pair_impl {
      shared_ptr<const Topology> p_space;
      distance_metric_type m_distance;
      proper_metric_type m_proper_distance;
      
      double distance(const point_type& a, const point_type& b) const {
        return m_distance(a, b, *p_space);
      };
      double proper_distance(const point_type& a, const point_type& b) const {
        return m_proper_distance(a, b, *p_space);
      };
      double operator()(const point_type& a, const point_type& b) const {
        return proper_distance(a, b, *p_space);
      };
      
      dist_metric_pair_impl(const shared_ptr<const Topology>& aSpace) : p_space(aSpace), 
        m_distance(get(distance_metric, *aSpace)), m_proper_distance(get(proper_metric, *aSpace)) { };
    };
    
    typedef typename boost::mpl::if_<
      boost::is_same<distance_metric_type, proper_metric_type>, 
      dist_metric_only_impl, dist_metric_pair_impl >::type parting_metrics_type;
    
    /** Type of the key-values that identify entries of the DVP tree. */ 
    typedef typename boost::property_traits< VertexKeyMap >::value_type key_type;
    
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<distance_metric_type, Topology>));
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<proper_metric_type, Topology>));
    
  private:
    
    typedef TreeType tree_indexer;
    
    typedef typename boost::graph_traits<tree_indexer>::vertex_descriptor vertex_type;
    typedef typename boost::graph_traits<tree_indexer>::edge_descriptor edge_type;
    typedef typename boost::graph_traits<tree_indexer>::out_edge_iterator out_edge_iter;
    typedef typename boost::graph_traits<tree_indexer>::in_edge_iterator in_edge_iter;
    
    typedef typename tree_indexer::vertex_property_type vertex_property;
    typedef typename tree_indexer::edge_property_type edge_property;
    
    struct priority_compare_type {
      bool operator()(const std::pair<distance_type, vertex_type>& x, const std::pair<distance_type, vertex_type>& y) const {
        return (x.first < y.first);
      };
    };
    typedef std::vector< std::pair< distance_type, vertex_type > > priority_queue_type;
    
    
    
    typedef boost::unordered_map<key_type, distance_type, detail::dvp_tree_key_hasher> temporary_dist_map_type;
    
    /* Simple comparison function to sort by pre-computed distance values. */
    struct closer {
      temporary_dist_map_type* p_m;
      VertexKeyMap m_key;
      closer(temporary_dist_map_type* m, const VertexKeyMap& key) : p_m(m), m_key(key) { };
      
      bool operator()(const vertex_property& vp1, const vertex_property& vp2) const {
        return (*p_m)[get(m_key, vp1)] < (*p_m)[get(m_key, vp2)];
      };
    };
    
    typedef boost::unordered_set<key_type, detail::dvp_tree_key_hasher> temporary_invalid_key_set_type;
    
    /* Simple predicate function to filter out invalid vertices (invalid key-value). */
    struct is_vertex_prop_valid {
      temporary_invalid_key_set_type* m_invalid_keys;
      VertexKeyMap m_key;
      is_vertex_prop_valid(temporary_invalid_key_set_type* invalid_keys, const VertexKeyMap& key) : m_invalid_keys(invalid_keys), m_key(key) { };
      
      bool operator()(const vertex_property& vp) const {
        return (m_invalid_keys->count(get(m_key, vp)) == 0);
      };
    };
    
    
    
    
    tree_indexer& m_tree;   ///< Tree storage.
    vertex_type m_root;    ///< Root node of the tree.
    
    VertexKeyMap m_key;      ///< A map from a vertex_property to a key_type value.
    DistanceMap m_mu;        ///< A map from an edge-property to a distance value.
    PositionMap m_position;  ///< A map from a vertex_property to a position value (should be Read-Write).
    
    parting_metrics_type m_distance;  ///< The distance-metric functor.
    
    VPChooser m_vp_chooser;  ///< The vantage-point chooser (functor).
    
    //non-copyable.
    dvp_tree_impl(const self&);
    self& operator=(const self&); 
    
    
    
    typedef typename std::vector<vertex_property>::iterator prop_vector_iter;
    
    struct construction_task {
      vertex_type node;
      prop_vector_iter first;
      prop_vector_iter last;
      
      construction_task(vertex_type aNode, prop_vector_iter aBegin, prop_vector_iter aEnd) : 
        node(aNode), first(aBegin), last(aEnd) { };
    };
    
    prop_vector_iter rearrange_with_chosen_vp(prop_vector_iter aBegin, prop_vector_iter aEnd) const {
      using std::iter_swap;
      
      // choose a vantage-point in the interval:
      prop_vector_iter chosen_vp_it = m_vp_chooser(aBegin, aEnd, m_distance, m_position);
      if(chosen_vp_it == aEnd)
        return aBegin; // no vp to be chosen in this interval (presumably, empty interval).
      
      // place the vantage-point (or pivot) at the start of interval:
      iter_swap(chosen_vp_it, aBegin);
      chosen_vp_it = aBegin;  //<-- chosen_vp_it is now at the start of interval.
      
      return chosen_vp_it;
    };
    
    
    /* NOTE Invalidates vertices */
    /* NOTE This is a non-recursive version of the construct-node algorithm */
    /* Does not require persistent vertices */
    /* This is the main tree construction function. It takes the vertices in the iterator range and organizes them 
     * as a sub-tree below the aParentNode node (and aEdgeDist is the minimum distance to the parent of any node in the range). */
    // aNode is a valid, existing node containing the chosen vantage-point.
    void construct_node(vertex_type aNode, prop_vector_iter aBegin, prop_vector_iter aEnd) {
      
      temporary_dist_map_type dist_map;
      std::queue<construction_task> tasks;
      tasks.push(construction_task(aNode, aBegin, aEnd));
      
      while(!tasks.empty()) {
        construction_task cur_task = tasks.front(); tasks.pop();
        
        // update values in the dist-map with the distances to the new chosen vantage-point:
        {
          const point_type& chosen_vp_pt = get(m_position, get_raw_vertex_property(m_tree, cur_task.node));
          for(prop_vector_iter it = cur_task.first; it != cur_task.last; ++it)
            dist_map[ get(m_key, *it) ] = m_distance.proper_distance(chosen_vp_pt, get(m_position, *it));
        };
        
        // this loop splits up the children into as equal as possible partitions.
        std::size_t total_count = (cur_task.last - cur_task.first);
        std::size_t child_count[Arity];
        for(std::size_t i = Arity; i > 0; --i) {
          child_count[i-1] = total_count / i;
          total_count -= child_count[i-1];
        };
        for(std::size_t i = 0; (i < Arity) && (child_count[i] > 0); --i) {
          std::nth_element(cur_task.first, cur_task.first + (child_count[i]-1), 
                           cur_task.last, closer(&dist_map,m_key));
          prop_vector_iter temp = cur_task.first; 
          cur_task.first += child_count[i];
          edge_property ep;
          put(m_mu, ep, dist_map[get(m_key,*(cur_task.first-1))]);
          rearrange_with_chosen_vp(temp, cur_task.first);
          dist_map.erase( get(m_key, *temp) );
          vertex_type new_vp_node;
          boost::tie(new_vp_node, boost::tuples::ignore) = 
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
            add_child_vertex(cur_task.node, std::move(*temp), std::move(ep), m_tree);
#else
            add_child_vertex(cur_task.node, *temp, ep, m_tree);
#endif
          ++temp;
          if(temp != cur_task.first)
            tasks.push(construction_task(new_vp_node, temp, cur_task.first));
        };
      };
    };
    
    
    
    
    struct nearest_search_result_set {
      priority_queue_type Neighbors;
      std::size_t K;
      distance_type Radius;
      
      nearest_search_result_set(std::size_t aK, distance_type aRadius) : Neighbors(), K(aK), Radius(aRadius) { };
      
      void register_vantage_point(const point_type&, const point_type&,
                                  distance_type current_dist, vertex_type current_vp, const parting_metrics_type&) {
        
        if(current_dist < Radius) { //is the vantage point within current search bound? Yes...
          //then add the vantage point to the NN list.
          Neighbors.push_back(std::pair<distance_type, vertex_type>(current_dist, current_vp));
          std::push_heap(Neighbors.begin(), Neighbors.end(), priority_compare_type());
          if(Neighbors.size() > K) { //are there too many nearest neighbors? Yes... 
            std::pop_heap(Neighbors.begin(), Neighbors.end(), priority_compare_type());
            Neighbors.pop_back(); //delete last element to keep Neighbors with K elements
            Radius = Neighbors.front().first; //distance of the last element is now the search bound Radius.
          };
        };
        
      };
      
    };
    
    struct pred_succ_search_result_set {
      priority_queue_type Pred;
      priority_queue_type Succ;
      std::size_t K;
      distance_type Radius;
      distance_type sigma_pred;
      distance_type sigma_succ;
      
      pred_succ_search_result_set(std::size_t aK, distance_type aRadius) : 
        Pred(), Succ(), K(aK), Radius(aRadius), sigma_pred(aRadius), sigma_succ(aRadius) { };
      
      void register_vantage_point(const point_type& query_point, const point_type& current_vp_point,
                                  distance_type current_dist, vertex_type current_vp, const parting_metrics_type& dual_distance) {
        
        // using the assumption that (current_dist <= current_pred_dist)  -->  if (current_pred_dist < sigma_pred) then (current_dist < sigma_pred)
        if(current_dist < sigma_pred) {
          distance_type current_pred_dist = dual_distance.distance(current_vp_point, query_point);
          if(current_pred_dist < sigma_pred) { //is the vantage point within current search bound? Yes...
            //then add the vantage point to the NN list.
            Pred.push_back(std::pair<distance_type, vertex_type>(current_pred_dist, current_vp));
            std::push_heap(Pred.begin(), Pred.end(), priority_compare_type());
            if(Pred.size() > K) { //are there too many nearest neighbors? Yes... 
              std::pop_heap(Pred.begin(), Pred.end(), priority_compare_type());
              Pred.pop_back(); //delete last element to keep aList with K elements
              sigma_pred = Pred.front().first; //distance of the last element is now the search bound sigma_pred.
            };
          };
        };
        
        // using the assumption that (current_dist <= current_succ_dist)  -->  if (current_succ_dist < sigma_succ) then (current_dist < sigma_succ)
        if(current_dist < sigma_succ) {
          distance_type current_succ_dist = dual_distance.distance(query_point, current_vp_point);
          if(current_succ_dist < sigma_succ) { //is the vantage point within current search bound? Yes...
            //then add the vantage point to the NN list.
            Succ.push_back(std::pair<distance_type, vertex_type>(current_succ_dist, current_vp));
            std::push_heap(Succ.begin(), Succ.end(), priority_compare_type());
            if(Succ.size() > K) { //are there too many nearest neighbors? Yes... 
              std::pop_heap(Succ.begin(), Succ.end(), priority_compare_type());
              Succ.pop_back(); //delete last element to keep aList with K elements
              sigma_succ = Succ.front().first; //distance of the last element is now the search bound sigma_succ.
            };
          };
        };
        // Radius must remain the maximum of both sigma_pred and sigma_succ (most inclusive search).
        Radius = ( (sigma_succ > sigma_pred) ? sigma_succ : sigma_pred );
        
      };
      
    };
    
    
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* NOTE This is a non-recursive version. */
    /* This is the main nearest-neighbor query function. This takes a query point, a maximum 
     * neighborhood radius (aSigma), aNode to start recursing from, the current max-heap of neighbors,
     * and the maximum number of neighbors. This function can be used for any kind of NN query (single, kNN, or ranged). */
    template <typename SearchResultSet>
    void find_nearest_impl(const point_type& aPoint, SearchResultSet& aResult) const {
      
      std::stack< std::pair<vertex_type, distance_type> > tasks;
      tasks.push(std::pair<vertex_type, distance_type>(m_root,0.0));
      
      while(!tasks.empty()) {
        std::pair<vertex_type, distance_type> cur_node = tasks.top(); tasks.pop();
        
        if( cur_node.second > aResult.Radius )
          continue;
        
        const point_type& current_vp = get(m_position, get_raw_vertex_property(m_tree,cur_node.first));
        distance_type current_dist = m_distance.proper_distance(aPoint, current_vp);
        
        aResult.register_vantage_point(aPoint, current_vp, current_dist, cur_node.first, m_distance);
        
        out_edge_iter ei,ei_end;
        //first, locate the partition in which aPoint is:
        if(out_degree(cur_node.first,m_tree) == 0)
          continue;
        for(boost::tie(ei,ei_end) = out_edges(cur_node.first,m_tree); ei != ei_end; ++ei) {
          if(current_dist <= get(m_mu, get_raw_edge_property(m_tree,*ei))) 
            break;
        };
        if(ei == ei_end) 
          --ei; //back-track if the end was reached.
        
        std::stack< std::pair<vertex_type, distance_type> > temp_invtasks;
        //search in the most likely node.
        temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei,m_tree),0.0));
        
        out_edge_iter ei_left = ei;
        out_edge_iter ei_right = ei; ++ei_right;
        boost::tie(ei,ei_end) = out_edges(cur_node.first,m_tree); //find the bounds again (start and end).
        bool left_stopped  = (ei_left == ei);
        bool right_stopped = (ei_right == ei_end);
        while(true) {
          if(left_stopped) {
            out_edge_iter ei_rightleft = ei_right; --ei_rightleft;
            distance_type temp_dist = 0.0;
            while((ei_right != ei_end) && 
                  ((temp_dist = get(m_mu,get_raw_edge_property(m_tree,*ei_rightleft)) - current_dist) < aResult.Radius)) {
              temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_right,m_tree), temp_dist));
              ++ei_rightleft; ++ei_right;
            };
            break;
          } else if(right_stopped) {
            out_edge_iter ei_leftleft = ei_left;
            distance_type temp_dist = 0.0;
            while((ei_left != ei) && 
                  ((temp_dist = current_dist - get(m_mu,get_raw_edge_property(m_tree,*(--ei_leftleft)))) < aResult.Radius)) {
              temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_leftleft,m_tree), temp_dist));
              --ei_left;
            };
            break;
          } else {
            out_edge_iter ei_leftleft = ei_left; --ei_leftleft;
            distance_type d1 = get(m_mu,get_raw_edge_property(m_tree,*ei_leftleft)); //greater than 0 if ei_leftleft should be searched.
            out_edge_iter ei_rightleft = ei_right; --ei_rightleft;
            distance_type d2 = get(m_mu,get_raw_edge_property(m_tree,*ei_rightleft)); //less than 0 if ei_right should be searched.
            if(d1 + d2 > 2.0 * current_dist) { //this means that ei_leftleft's boundary is closer to aPoint.
              if(d1 + aResult.Radius - current_dist > 0) {
                temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_leftleft,m_tree), current_dist - d1));
                ei_left = ei_leftleft;
                if(d2 - aResult.Radius - current_dist < 0) {
                  temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_right,m_tree), d2 - current_dist));
                  ++ei_right;
                } else
                  right_stopped = true;
              } else
                break;
            } else {
              if(d2 - aResult.Radius - current_dist < 0) {
                temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_right,m_tree), d2 - current_dist));
                ++ei_right;
                if(d1 + aResult.Radius - current_dist > 0) {
                  temp_invtasks.push(std::pair<vertex_type, distance_type>(target(*ei_leftleft,m_tree), current_dist - d1));
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
        
        // reverse the temporary stack into the main stack.
        while(!temp_invtasks.empty()) {
          tasks.push(temp_invtasks.top());
          temp_invtasks.pop();
        };
      };
    };
    
    
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* NOTE This is a non-recursive version. */
    /* This function does a single nearest-neighbor query to find the tree-leaf that is closest to 
     * the given point (starts to recurse from aNode). */
    vertex_type get_leaf(const point_type& aPoint, vertex_type aNode) const {
      while(out_degree(aNode,m_tree) != 0) {
        //first, locate the partition in which aPoint is:
        distance_type current_dist = m_distance.proper_distance(aPoint, get(m_position, get_raw_vertex_property(m_tree,aNode)));
        vertex_type result = aNode;
        out_edge_iter ei,ei_end;
        for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
          result = target(*ei,m_tree);
          if(current_dist <= get(m_mu, get_raw_edge_property(m_tree,*ei))) 
            break;
        };
        aNode = result;
      };
      return aNode;
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* NOTE This is a non-recursive version. */
    /* This function finds the tree-vertex with the given key-value and position. */
    vertex_type get_vertex(key_type aKey, const point_type& aPoint, vertex_type aNode) const {
      vertex_type aAlternateBranch = boost::graph_traits<tree_indexer>::null_vertex();
//       std::cout << "\n ------ Looking for node " << aKey << " at position: " << aPoint << std::endl;
      while(get(m_key, get_raw_vertex_property(m_tree,aNode)) != aKey) { 
        distance_type current_dist = m_distance.proper_distance(aPoint, get(m_position, get_raw_vertex_property(m_tree,aNode)));
//         std::cout << " ---- Looking at node " << get(m_key, get_raw_vertex_property(m_tree,aNode)) 
//                   << " at position: " << get(m_position, get_raw_vertex_property(m_tree,aNode))
//                   << " at distance " << current_dist << std::endl;
        //first, locate the partition in which aPoint is:
        if(out_degree(aNode,m_tree) == 0) {
          if(aAlternateBranch != boost::graph_traits<tree_indexer>::null_vertex()) {
            aNode = aAlternateBranch;
            aAlternateBranch = boost::graph_traits<tree_indexer>::null_vertex();
            continue;
          } else
            throw int(0);
        };
        vertex_type result = aNode;
        out_edge_iter ei,ei_end;
        for(boost::tie(ei,ei_end) = out_edges(aNode,m_tree); ei != ei_end; ++ei) {
          result = target(*ei,m_tree);
//           std::cout << " -- Looking at child " << std::setw(6) << get(m_key, get_raw_vertex_property(m_tree, result)) 
//                     << " at distance " << get(m_mu, get_raw_edge_property(m_tree, *ei)) << std::endl;
          if(current_dist < get(m_mu, get_raw_edge_property(m_tree,*ei))) 
            break;
          if(current_dist == get(m_mu, get_raw_edge_property(m_tree,*ei))) {
            ++ei;
            if(ei != ei_end)
              aAlternateBranch = target(*ei,m_tree);
            break;
          };
        };
        aNode = result;
      };
      return aNode;
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* NOTE This is a non-recursive version. */
    /* This function updates the edge distance values as a consequence of a new point being added. */
    void update_mu_upwards(const point_type& aPoint, vertex_type aNode) {
      while(aNode != m_root) {
        vertex_type parent = source(*(in_edges(aNode,m_tree).first), m_tree);
        distance_type dist = m_distance.proper_distance(aPoint, get(m_position, get_raw_vertex_property(m_tree,parent)));
        if(dist > get(m_mu, get_raw_edge_property(m_tree, *(in_edges(aNode,m_tree).first))))
          put(m_mu, get_raw_edge_property(m_tree, *(in_edges(aNode,m_tree).first)), dist);
        aNode = parent;
      };
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
    /* NOTE This is a non-recursive version. */
    /* This function determines if a given node is the root of a balanced (and full) sub-tree. */
    bool is_node_full(vertex_type aNode, int& depth_limit) const {
      if(depth_limit < 0)
        return false;
      std::queue< std::pair< vertex_type, int > > tasks;
      tasks.push(std::pair< vertex_type, int >(aNode, depth_limit));
      while(!tasks.empty()) {
        std::pair< vertex_type, int > cur_task = tasks.front(); tasks.pop();
        
        if(cur_task.second < depth_limit)
          depth_limit = cur_task.second;
        
        if((out_degree(cur_task.first, m_tree) == 0) && (cur_task.second == 0))
          continue;
        
        --(cur_task.second);
        
        if(((out_degree(cur_task.first, m_tree) != 0) && (cur_task.second < 0)) || 
           (out_degree(cur_task.first, m_tree) < Arity) ||
           ((cur_task.second > 0) && (is_leaf_node(cur_task.first))) ) {
          depth_limit = cur_task.second;
          return false;
        };
        
        out_edge_iter ei,ei_end;
        for(boost::tie(ei,ei_end) = out_edges(cur_task.first,m_tree); ei != ei_end; ++ei)
          tasks.push(std::pair< vertex_type, int >(target(*ei,m_tree), cur_task.second));
        
      };
      return (depth_limit == 0);
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* NOTE This is a non-recursive version. */
    /* This function collects the list of vertices that are in the sub-tree rooted at the given node (exclusively). */
    void collect_vertices(std::vector<vertex_type>& aList, vertex_type aNode) const {
      std::queue< vertex_type > tasks;
      tasks.push(aNode);
      while(!tasks.empty()) {
        vertex_type current_node = tasks.front(); tasks.pop();
        out_edge_iter ei,ei_end;
        for(boost::tie(ei,ei_end) = out_edges(current_node,m_tree); ei != ei_end; ++ei) {
          aList.push_back(target(*ei,m_tree));
          tasks.push(target(*ei,m_tree));
        };
      };
    };
    
    /* Does not invalidate vertices */
    /* Does not require persistent vertices */
    /* This function computes the maximum depth of the tree rooted at the given node (note: this is an 
     * expensive operation, but is useful when debugging and testing). */
    /* NOTE: This is a recursive function, but it isn't practical anyways, and kind of annoying to de-recursify. */
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
     * \param aKey The key-map to use to obtain and store the key values for a given vertex-property object.
     * \param aMu The property-map which associates a distance-value to each edge of the tree.
     * \param aPosition The property-map that can be used to obtain and store the positions of the vertices (vertex-property objects).
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    template <typename Graph, typename GraphPositionMap>
    dvp_tree_impl(const Graph& g, 
                  GraphPositionMap aGraphPosition,
                  tree_indexer& aTree,
                  const shared_ptr<const Topology>& aSpace, 
                  VertexKeyMap aKey,
                  DistanceMap aMu,
                  PositionMap aPosition,
                  VPChooser aVPChooser) : 
                  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
                  m_key(aKey), m_mu(aMu), m_position(aPosition),
                  m_distance(aSpace), 
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
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        v_bin.push_back(std::move(vp));
#else
        v_bin.push_back(vp);
#endif
      };
      
      prop_vector_iter v_first = v_bin.begin();
      prop_vector_iter v_last  = v_bin.end();
      rearrange_with_chosen_vp(v_first, v_last);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      m_root = create_root(std::move(*v_first), m_tree);
#else
      m_root = create_root(*v_first, m_tree);
#endif
      construct_node(m_root, ++v_first, v_last);
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
                  const shared_ptr<const Topology>& aSpace, 
                  VertexKeyMap aKey,
                  DistanceMap aMu,
                  PositionMap aPosition,
                  VPChooser aVPChooser) : 
                  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
                  m_key(aKey), m_mu(aMu), m_position(aPosition),
                  m_distance(aSpace), 
                  m_vp_chooser(aVPChooser) {
      if(aBegin == aEnd) return;
      
      std::vector<vertex_property> v_bin; //Copy the list of vertices to random access memory.
      for(;aBegin != aEnd; ++aBegin) {
        vertex_property vp;
        put(m_key, vp, *aBegin);
        put(m_position, vp, get(aElemPosition, *aBegin));
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        v_bin.push_back(std::move(vp));
#else
        v_bin.push_back(vp);
#endif
      };
      
      prop_vector_iter v_first = v_bin.begin();
      prop_vector_iter v_last  = v_bin.end();
      rearrange_with_chosen_vp(v_first, v_last);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      m_root = create_root(std::move(*v_first), m_tree);
#else
      m_root = create_root(*v_first, m_tree);
#endif
      construct_node(m_root, ++v_first, v_last);
    };
    
    /**
     * Construct an empty DVP-tree from a topology, etc..
     * \param aTree The tree object that will be used to store the DVP structure.
     * \param aSpace The topology on which the positions of the vertices reside.
     * \param aKey The key-map to use to obtain and store the key values for a given vertex-property object.
     * \param aMu The property-map which associates a distance-value to each edge of the tree.
     * \param aPosition The property-map that can be used to obtain and store the positions of the vertices (vertex-property objects).
     * \param aVPChooser The vantage-point chooser functor (policy class).
     */
    dvp_tree_impl(tree_indexer& aTree,
                  const shared_ptr<const Topology>& aSpace, 
                  VertexKeyMap aKey,
                  DistanceMap aMu,
                  PositionMap aPosition,
                  VPChooser aVPChooser) : 
                  m_tree(aTree), m_root(boost::graph_traits<tree_indexer>::null_vertex()), 
                  m_key(aKey), m_mu(aMu), m_position(aPosition),
                  m_distance(aSpace), 
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
     * This function computes an approximation of the characteristic size of the vertices in the DVP tree.
     * \return The approximation of the characteristic size of the vertices in the DVP tree.
     */
    double get_characteristic_size() const { 
      if(num_vertices(m_tree) == 0)
        return std::numeric_limits<double>::infinity();
      
      out_edge_iter ei,ei_end;
      double max_dist = 0.0;
      for(boost::tie(ei,ei_end) = out_edges(m_root,m_tree); ei != ei_end; ++ei) {
        double cur_dist = get(m_mu, get_raw_edge_property(m_tree,*ei));
        if(cur_dist > max_dist)
          max_dist = cur_dist;
      };
      if(max_dist != 0.0)
        return max_dist;
      else  // never found a finite, non-zero distance value.
        return std::numeric_limits<double>::infinity();
    };
    
    /**
     * Inserts a vertex into the tree.
     * \param up The vertex-property to be added to the DVP-tree.
     */
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    void insert(vertex_property up) {
#else
    void insert(const vertex_property& up) {
#endif
      if(num_vertices(m_tree) == 0) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        m_root = create_root(std::move(up), m_tree); 
#else
        m_root = create_root(up, m_tree); 
#endif
        return;
      };
      
      point_type u_pt = get(m_position, up); 
      vertex_type u_subroot = get_leaf(u_pt,m_root); // <-- to store the root of subtree to reconstruct.
      // NOTE: if the root is the leaf, it requires special attention since no parent exists.
      if(u_subroot != m_root) {
        vertex_type u_leaf = source(*(in_edges(u_subroot, m_tree).first),m_tree);
        if((out_degree(u_leaf, m_tree) == Arity) && (is_leaf_node(u_leaf))) {
          //if u_leaf is a full-leaf, then it is balanced but full, 
          // we should then find a non-full parent.
          int actual_depth_limit = 1;
          int last_depth_limit = actual_depth_limit;
          while((u_leaf != m_root) && (is_node_full(u_leaf, last_depth_limit))) {
            u_leaf = source(*(in_edges(u_leaf, m_tree).first), m_tree);
            last_depth_limit = ++actual_depth_limit;
          };
          bool is_p_full = false; 
          if(u_leaf == m_root)
            is_p_full = is_node_full(u_leaf, last_depth_limit);
          if((!is_p_full) && (last_depth_limit >= 0)) {
            // this means that we can add our key to the sub-tree of u_leaf and reconstruct from there.
            u_subroot = u_leaf;
          };
          // else:
          //  this means that either the root node is full or there are 
          //  branches of the tree that are deeper than u_subroot, 
          //  and thus, in either case, u_subroot should be expanded.
        } else {
          //  leaf node is not full of children, an additional child can be added 
          //  (must be reconstructed to keep ordering, but this is a trivial operation O(Arity)).
          //OR 
          //  if leaf is not really a leaf, then it means that this sub-tree is definitely 
          //  not balanced and not full either,
          //  then all the Keys ought to be collected and u_leaf ought to be reconstructed.
          u_subroot = u_leaf;
        };
      };
      
      update_mu_upwards(u_pt, u_subroot);
      std::vector<vertex_property> prop_list;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      prop_list.push_back(std::move(up));
#else
      prop_list.push_back(up);
#endif
      while( out_degree(u_subroot, m_tree) > 0 ) {
        edge_type e = *(out_edges(u_subroot, m_tree).first);
        remove_branch(target(e, m_tree), back_inserter(prop_list), m_tree);
      };
      construct_node(u_subroot, prop_list.begin(), prop_list.end());
      
    };
    /**
     * Inserts a range of vertices.
     * \tparam ForwardIterator A forward-iterator type that can be used to obtain the vertices.
     * \param aBegin The start of the range from which to take the vertices.
     * \param aEnd The end of the range from which to take the vertices (one-past-last).
     */
    template <typename ForwardIterator>
    void insert(ForwardIterator aBegin, ForwardIterator aEnd) { 
      for(; aBegin != aEnd; ++aBegin) 
        insert(*aBegin);
      //TODO: There's got to be a better way to insert many elements (most likely a similar strategy to the erase multiple function).
    };
    
    
    /**
     * Erases the given vertex from the DVP-tree.
     * \param u_node The vertex to be removed from the DVP-tree.
     */
    void erase(vertex_type u_node) { 
      if(num_vertices(m_tree) == 0) 
        return;
      std::vector<vertex_property> prop_list;
      if( (u_node == m_root) && (num_vertices(m_tree) == 1) ) {
        remove_branch(m_root, back_inserter(prop_list), m_tree);
        m_root = boost::graph_traits<tree_indexer>::null_vertex();
        return;
      };
      vertex_type u_parent = boost::graph_traits<tree_indexer>::null_vertex();
      if(u_node != m_root)
        u_parent = source(*(in_edges(u_node,m_tree).first), m_tree);
      
      // remove-and-collect all children of u_node:
      while( out_degree(u_node, m_tree) > 0 ) {
        edge_type e = *(out_edges(u_node, m_tree).first);
        remove_branch(target(e, m_tree), back_inserter(prop_list), m_tree);
      };
      // remove-and-discard u_node:
      {
        std::vector<vertex_property> throwaway_list;
        remove_branch(u_node, back_inserter(throwaway_list), m_tree);
      };
      
      if( u_parent != boost::graph_traits<tree_indexer>::null_vertex() ) {
        // remove-and-collect all other children of u_parent:
        while( out_degree(u_parent, m_tree) > 0 ) {
          edge_type e = *(out_edges(u_parent, m_tree).first);
          remove_branch(target(e, m_tree), back_inserter(prop_list), m_tree);
        };
        // reconstruct parent's subtree:
        construct_node(u_parent, prop_list.begin(), prop_list.end());
      } else {
        // need to re-construct the root node (u_node == m_root).
        prop_vector_iter v_first = prop_list.begin();
        prop_vector_iter v_last  = prop_list.end();
        rearrange_with_chosen_vp(v_first, v_last);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        m_root = create_root(std::move(*v_first), m_tree);
#else
        m_root = create_root(*v_first, m_tree);
#endif
        construct_node(m_root, ++v_first, v_last);
      };
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
        std::cout << " Could not find the node to be removed from the DVP tree!" << std::endl;
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
      temporary_invalid_key_set_type invalid_keys;
      
      // First, generate the vertex-listings in preparation for the deletion.
      for(ForwardIterator first = aBegin; first != aEnd; ++first) {
        vertex_type removal_trunk = *first;
        // mark as invalid, for deletion.
        invalid_keys.insert(get(m_key, get_raw_vertex_property(m_tree, *first))); 
        // go up until a valid parent node is found:
        while( (removal_trunk != m_root) && 
               (invalid_keys.count(get(m_key, get_raw_vertex_property(m_tree, removal_trunk))) != 0) ) {
          removal_trunk = source(*(in_edges(removal_trunk, m_tree).first), m_tree);
        };
        
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
      
      if( (v_lists.size() == 1) && (v_lists.front() == m_root) ) {
        // need to re-construct the root node (u_node == m_root).
        std::vector<vertex_property> prop_list;
        remove_branch(m_root, back_inserter(prop_list), m_tree);
        prop_list.erase( remove_if(prop_list.begin(), prop_list.end(), is_vertex_prop_valid(&invalid_keys, m_key)), prop_list.end());
        prop_vector_iter v_first = prop_list.begin();
        prop_vector_iter v_last  = prop_list.end();
        rearrange_with_chosen_vp(v_first, v_last);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        m_root = create_root(std::move(*v_first), m_tree);
#else
        m_root = create_root(*v_first, m_tree);
#endif
        construct_node(m_root, ++v_first, v_last);
        return;
      };
      
      for(typename vertex_listing::iterator it = v_lists.begin(); it != v_lists.end(); ++it) {
        std::vector<vertex_property> prop_list;
        // remove-and-collect all children of it->first:
        while( out_degree(it->first, m_tree) > 0 ) {
          edge_type e = *(out_edges(it->first, m_tree).first);
          remove_branch(target(e, m_tree), back_inserter(prop_list), m_tree);
        };
        // erase removed vertices from the prop-list:
        prop_list.erase( remove_if(prop_list.begin(), prop_list.end(), is_vertex_prop_valid(&invalid_keys, m_key)), prop_list.end());
        // reconstruct parent's subtree:
        construct_node(it->first, prop_list.begin(), prop_list.end());
      };
    };
    
    /**
     * Clears the DVP-tree. 
     */
    void clear() {
      if( num_vertices(m_tree) == 0 ) {
        std::vector<vertex_property> prop_list;
        remove_branch(m_root, back_inserter(prop_list), m_tree);
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
      nearest_search_result_set result_set(1, std::numeric_limits<distance_type>::infinity());
      find_nearest_impl(aPoint, result_set);
      if(result_set.Neighbors.size())
        return result_set.Neighbors.front().second;
      else
        return boost::graph_traits<tree_indexer>::null_vertex();
    };
    
    /**
     * Finds the nearest predecessor and successor to a given position.
     * \param aPoint The position from which to find the nearest predecessor and successor of.
     * \return The predecessor and successor vertex in the DVP-tree that is closest to the given point.
     */
    std::pair< vertex_type, vertex_type > find_nearest_pred_succ(const point_type& aPoint) const {
      if(num_vertices(m_tree) == 0) 
        return boost::graph_traits<tree_indexer>::null_vertex();
      pred_succ_search_result_set result_set(1, std::numeric_limits<distance_type>::infinity());
      find_nearest_impl(aPoint, result_set);
      std::pair< vertex_type, vertex_type > result;
      if( result_set.Pred.size() )
        result.first = result_set.Pred.front().second;
      else
        result.first = boost::graph_traits<tree_indexer>::null_vertex();
      if( result_set.Succ.size() )
        result.second = result_set.Succ.front().second;
      else
        result.second = boost::graph_traits<tree_indexer>::null_vertex();
      return result;
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
      nearest_search_result_set result_set(K, R);
      find_nearest_impl(aPoint, result_set);
      std::sort_heap(result_set.Neighbors.begin(), result_set.Neighbors.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = result_set.Neighbors.begin(); it != result_set.Neighbors.end(); ++it)
        *(aOutputBegin++) = it->second;
      return aOutputBegin;
    };
    
    /**
     * Finds the K nearest predecessors and successors to a given position.
     * \tparam OutputIterator The forward- output-iterator type which can contain the 
     *         list of nearest-neighbors by tree vertex descriptors.
     * \param aPoint The position from which to find the nearest-neighbors.
     * \param aPredBegin An iterator to the first place where to put the sorted list of 
     *        predecessors by tree vertex descriptors with the smallest distance.
     * \param aSuccBegin An iterator to the first place where to put the sorted list of 
     *        successors by tree vertex descriptors with the smallest distance.
     * \param K The number of nearest-neighbors.
     * \param R The maximum distance value for the nearest-neighbors.
     * \return The output-iterators to the end of the lists of nearest predecessors and successors (starting from "output_first").
     */
    template <typename OutputIterator>
    std::pair< OutputIterator, OutputIterator > find_nearest(
        const point_type& aPoint, OutputIterator aPredBegin, OutputIterator aSuccBegin, 
        std::size_t K, distance_type R = std::numeric_limits<distance_type>::infinity()) const {
      if(num_vertices(m_tree) == 0) 
        return std::pair< OutputIterator, OutputIterator >(aPredBegin, aSuccBegin);
      pred_succ_search_result_set result_set(K, R);
      find_nearest_impl(aPoint, result_set);
      std::sort_heap(result_set.Pred.begin(), result_set.Pred.end(), priority_compare_type());
      std::sort_heap(result_set.Succ.begin(), result_set.Succ.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = result_set.Pred.begin(); it != result_set.Pred.end(); ++it)
        *(aPredBegin++) = it->second;
      for(typename priority_queue_type::const_iterator it = result_set.Succ.begin(); it != result_set.Succ.end(); ++it)
        *(aSuccBegin++) = it->second;
      return std::pair< OutputIterator, OutputIterator >(aPredBegin, aSuccBegin);
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
      nearest_search_result_set result_set(num_vertices(m_tree), R);
      find_nearest_impl(aPoint, result_set);
      std::sort_heap(result_set.Neighbors.begin(), result_set.Neighbors.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = result_set.Neighbors.begin(); it != result_set.Neighbors.end(); ++it)
        *(aOutputBegin++) = it->second;
      return aOutputBegin;
    };
    
    /**
     * Finds the nearest predecessors and successors to a given position within a given range (radius).
     * \tparam OutputIterator The forward- output-iterator type which can contain the 
     *         list of nearest-neighbors by tree vertex descriptors.
     * \param aPoint The position from which to find the nearest-neighbors.
     * \param aPredBegin An iterator to the first place where to put the sorted list of 
     *        predecessors by tree vertex descriptors with the smallest distance.
     * \param aSuccBegin An iterator to the first place where to put the sorted list of 
     *        successors by tree vertex descriptors with the smallest distance.
     * \param R The maximum distance value for the nearest-neighbors.
     * \return The output-iterators to the end of the lists of nearest predecessors and successors (starting from "output_first").
     */
    template <typename OutputIterator>
    std::pair< OutputIterator, OutputIterator > find_in_range(
        const point_type& aPoint, OutputIterator aPredBegin, OutputIterator aSuccBegin, distance_type R) const {
      if(num_vertices(m_tree) == 0) 
        return std::pair< OutputIterator, OutputIterator >(aPredBegin, aSuccBegin);
      pred_succ_search_result_set result_set(num_vertices(m_tree), R);
      find_nearest_impl(aPoint, result_set);
      std::sort_heap(result_set.Pred.begin(), result_set.Pred.end(), priority_compare_type());
      std::sort_heap(result_set.Succ.begin(), result_set.Succ.end(), priority_compare_type());
      for(typename priority_queue_type::const_iterator it = result_set.Pred.begin(); it != result_set.Pred.end(); ++it)
        *(aPredBegin++) = it->second;
      for(typename priority_queue_type::const_iterator it = result_set.Succ.begin(); it != result_set.Succ.end(); ++it)
        *(aSuccBegin++) = it->second;
      return std::pair< OutputIterator, OutputIterator >(aPredBegin, aSuccBegin);
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
      
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      void add_vertex(vertex_property&& vp, tree_indexer&) const {
        m_parent->insert(std::move(vp));
      };
#endif
    };
    
};





};

};


#endif


















