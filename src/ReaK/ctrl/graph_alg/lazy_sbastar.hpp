/**
 * \file lazy_sbastar.hpp
 *
 * This library provides function templates and concepts that implement a Lazy Sampling-based A* search
 * algorithm. A Lazy-SBA* uses the A* search algorithm to drive the expansion of a roadmap into the free-space 
 * in order to connect a start and goal location. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc. 
 * All these customization points are left to the user to implement, some are defined by the 
 * SBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The Lazy-SBA* algorithm is a generalization of the A* algorithm where the neighborhood of a given node of 
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an SBA* algorithm, the same 
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the 
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the 
 * "surprisal") from drawing a new sample in the neighborhood. In this lazy version, the computation of the 
 * edge weights as a cost-to-go through the free-space is tentatively replaced by the cost-to-go in the 
 * configuration space (without obstacles), and collision along the path is only performed once the edge 
 * has relaxed (identified as a segment of the local optimal path).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_LAZY_SBASTAR_HPP
#define REAK_LAZY_SBASTAR_HPP

#include "sbastar_search.hpp"

#include <functional>
#include <boost/utility/enable_if.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/exception.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

#include "bgl_more_property_maps.hpp"
#include "bgl_more_property_tags.hpp"
#include "bgl_raw_property_graph.hpp"

#include <stack>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {
  
  namespace detail {
    
    struct lazy_sbastar_node_connector {
      
      template <typename Graph,
                typename Topology,
                typename SBAVisitor,
                typename PositionMap,
                typename DistanceMap,
                typename PredecessorMap,
                typename WeightMap,
                typename NcSelector>
      typename boost::enable_if< boost::is_undirected_graph<Graph> >::type operator()(
          const typename boost::property_traits<PositionMap>::value_type& p, 
          typename boost::graph_traits<Graph>::vertex_descriptor u, 
          typename Graph::edge_bundled& ep, 
          Graph& g,
          const Topology& super_space,
          const SBAVisitor& sba_vis,
          PositionMap position,
          DistanceMap distance,
          PredecessorMap predecessor,
          WeightMap weight,
          NcSelector select_neighborhood) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Nc;
        select_neighborhood(p, std::back_inserter(Nc), g, super_space, g_position); 
        
        Vertex v = sba_vis.create_vertex(p, g);
        
        std::pair<Edge, bool> e_new = sba_vis.create_edge(u, v, ep, g);
        if( e_new.second ) {
          put(distance, v, get(distance, u) + get(weight, g[e_new.first]));
          put(predecessor, v, u);
          sba_vis.requeue_vertex(u,g);
        };
        
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
          if(*it == u)
            continue;
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(get(position,g[*it]), p, super_space);
          double g_out = tentative_weight + get(distance, *it);
          if(g_out < get(distance, v)) {
            // edge is useful as an in-edge to v.
            e_new = sba_vis.attempt_connecting_edge(*it, v, g);
            if( e_new.second ) {
              put(distance, v, g_out);
              Vertex old_pred = get(predecessor, v);
              put(predecessor, v, *it); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, v, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        };
        
        sba_vis.requeue_vertex(v,g);
        
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(p, get(position,g[*it]), super_space);
          double g_in  = tentative_weight + get(distance, v);
          if(g_in < get(distance, *it)) {
            // edge is useful as an in-edge to (*it).
            e_new = sba_vis.attempt_connecting_edge(v, *it, g);
            if( e_new.second ) {
              put(distance, *it, g_in);
              Vertex old_pred = get(predecessor, *it);
              put(predecessor, *it, v); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, *it, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        }; 
        
        
#if 0
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
          if(*it == u)
            continue;
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(get(position,g[*it]), p, super_space);
          double g_in  = tentative_weight + get(distance, v);
          double g_out = tentative_weight + get(distance, *it);
          if(g_in < get(distance, *it)) {
            // edge is useful as an in-edge to (*it).
            e_new = sba_vis.attempt_connecting_edge(v, *it, g);
            if( e_new.second ) {
              put(distance, *it, g_in);
              Vertex old_pred = get(predecessor, *it);
              put(predecessor, *it, v); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, *it, g);
            };
          } else if(g_out < get(distance, v)) {
            // edge is useful as an in-edge to v.
            e_new = sba_vis.attempt_connecting_edge(*it, v, g);
            if( e_new.second ) {
              put(distance, v, g_out);
              Vertex old_pred = get(predecessor, v);
              put(predecessor, v, *it); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, v, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        }; 
#endif
        
        sba_vis.requeue_vertex(v,g);
        
        // need to update all the children of the v node:
        std::stack<Vertex> incons;
        incons.push(v);
        while(!incons.empty()) {
          Vertex s = incons.top(); incons.pop();
          OutEdgeIter eo, eo_end;
          for(boost::tie(eo,eo_end) = out_edges(s,g); eo != eo_end; ++eo) {
            Vertex t = target(*eo, g);
            if(t == s)
              t = source(*eo, g);
            if(s != get(predecessor, t))
              continue;
            put(distance, t, get(distance, s) + get(weight, g[*eo]));
            
            sba_vis.requeue_vertex(t,g);
            
            incons.push(t);
          };
        };
      };
      
      template <typename Graph,
                typename Topology,
                typename SBAVisitor,
                typename PositionMap,
                typename DistanceMap,
                typename PredecessorMap,
                typename WeightMap,
                typename NcSelector>
      typename boost::enable_if< boost::is_directed_graph<Graph> >::type operator()(
          const typename boost::property_traits<PositionMap>::value_type& p, 
          typename boost::graph_traits<Graph>::vertex_descriptor u, 
          typename Graph::edge_bundled& ep, 
          Graph& g,
          const Topology& super_space,
          const SBAVisitor& sba_vis,
          PositionMap position,
          DistanceMap distance,
          PredecessorMap predecessor,
          WeightMap weight,
          NcSelector select_neighborhood) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Pred, Succ;
        select_neighborhood(p, std::back_inserter(Pred), std::back_inserter(Succ), g, super_space, g_position); 
        
        Vertex v = sba_vis.create_vertex(p, g);
        
        std::pair<Edge, bool> e_new = sba_vis.create_edge(u, v, ep, g);
        if( e_new.second ) {
          put(distance, v, get(distance, u) + get(weight, g[e_new.first]));
          put(predecessor, v, u);
          sba_vis.requeue_vertex(u,g);
        };
        
        for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
          if(*it == u)
            continue;
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(get(position,g[*it]), p, super_space);
          double g_out = tentative_weight + get(distance, *it);
          if(g_out < get(distance, v)) {
            // edge is useful as an in-edge to v.
            e_new = sba_vis.attempt_connecting_edge(*it, v, g);
            if( e_new.second ) {
              put(distance, v, g_out);
              Vertex old_pred = get(predecessor, v);
              put(predecessor, v, *it); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, v, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        };
        
        sba_vis.requeue_vertex(v,g);
        
        for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(p, get(position,g[*it]), super_space);
          double g_in  = tentative_weight + get(distance, v);
          if(g_in < get(distance, *it)) {
            // edge is useful as an in-edge to (*it).
            e_new = sba_vis.attempt_connecting_edge(v, *it, g);
            if( e_new.second ) {
              put(distance, *it, g_in);
              Vertex old_pred = get(predecessor, *it);
              put(predecessor, *it, v); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, *it, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        }; 
        
        sba_vis.requeue_vertex(v,g);
        
        // need to update all the children of the v node:
        std::stack<Vertex> incons;
        incons.push(v);
        while(!incons.empty()) {
          Vertex s = incons.top(); incons.pop();
          OutEdgeIter eo, eo_end;
          for(boost::tie(eo,eo_end) = out_edges(s,g); eo != eo_end; ++eo) {
            Vertex t = target(*eo, g);
            put(distance, t, get(distance, s) + get(weight, g[*eo]));
            
            sba_vis.requeue_vertex(t,g);
            
            incons.push(t);
          };
        };
      };
      
    };
    
  }; //end of detail namespace.

  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Lazy-SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_lazy_sbastar_no_init(const SBAStarBundle& bdl) {
    detail::generate_sbastar_no_init_impl< detail::lazy_sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, bdl.m_select_neighborhood);
  };

  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_lazy_sbastar(const SBAStarBundle& bdl) {
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_lazy_sbastar_no_init(bdl);
    
  };
  
  
  

};

};

#endif
















