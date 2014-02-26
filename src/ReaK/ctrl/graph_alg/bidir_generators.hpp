/**
 * \file bidir_generators.hpp
 * 
 * This library contains a bidirectional node generator that can be used in RRT-like algorithms. 
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_BIDIR_GENERATORS_HPP
#define REAK_BIDIR_GENERATORS_HPP

#include <utility>
#include <vector>
#include <iterator>
#include <boost/tuple/tuple.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"

namespace ReaK {
  
namespace graph {

namespace detail {
  
  template <typename PositionValue, typename Graph, typename RRGVisitor, typename PredecessorMap>
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor, 
                bool, typename Graph::edge_bundled > 
    expand_to_nearest(
      PositionValue& p_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Pred,
      Graph& g, RRGVisitor vis, PredecessorMap predecessor) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, bool, EdgeProp > ResultType;
    
    if(Pred.empty())
      return false;
    
    PositionValue p_tmp; 
    bool expand_succeeded = false;
    EdgeProp ep;
    typename std::vector<Vertex>::const_iterator it = Pred.begin();
    while((!expand_succeeded) && (it != Pred.end())) {
      if( get(predecessor, g[*it]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_tmp, expand_succeeded, ep) = vis.steer_towards_position(p_new, *it, g);
      ++it;
    };
    
    if(!expand_succeeded)
      return ResultType(boost::graph_traits<Graph>::null_vertex(), false,EdgeProp());
    
    p_new = p_tmp;
    return ResultType(*(--it), true, ep);
  };
  
  template <typename PositionValue, typename Graph, typename RRGVisitor, typename SuccessorMap>
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor, 
                bool, typename Graph::edge_bundled > 
    retract_from_nearest(
      PositionValue& p_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Succ,
      Graph& g, RRGVisitor vis, SuccessorMap successor) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, bool, EdgeProp > ResultType;
    
    if(Succ.empty())
      return false;
    
    PositionValue p_tmp;
    bool retract_succeeded = false;
    EdgeProp ep;
    typename std::vector<Vertex>::const_iterator it = Succ.begin();
    while((!retract_succeeded) && (it != Succ.end())) {
      if( get(successor, g[*it]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_tmp, retract_succeeded, ep) = vis.steer_back_to_position(p_new, *it, g);
      ++it;
    };
    
    if(!retract_succeeded)
      return ResultType(boost::graph_traits<Graph>::null_vertex(), false,EdgeProp());
    
    p_new = p_tmp;
    return ResultType(*(--it), true, ep);
  };
  
  
};
  
/**
 * This bidirectional node generator that can be used in RRT-like algorithms. 
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
 * \tparam Topology The topology type on which the planning is performed (i.e., the configuration space type).
 * \tparam RandomSampler The type of the random-sampler that can generate random points in the configuration space, should model ReaK::pp::RandomSamplerConcept.
 * \tparam NcSelector The type of a functor that can be used to perform nearest-neighbor queries.
 */
template <typename Topology, typename RandomSampler, typename NcSelector, 
          typename PredecessorMap, typename SuccessorMap>
struct rrg_bidir_generator {
  
  const Topology* space;
  RandomSampler get_sample;
  NcSelector select_neighborhood;
  PredecessorMap predecessor;
  SuccessorMap successor;
  
  rrg_bidir_generator(const Topology* aSpace,
                      RandomSampler aGetSample,
                      NcSelector aSelectNeighborhood,
                      PredecessorMap aPred, SuccessorMap aSucc) : 
                      space(aSpace),
                      get_sample(aGetSample),
                      select_neighborhood(aSelectNeighborhood),
                      predecessor(aPred), successor(aSucc) { };
  
  template <typename Graph,
            typename RRGVisitor,
            typename PositionMap>
  typename boost::enable_if< boost::is_undirected_graph<Graph>,
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled,
                typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, PositionValue, EdgeProp, Vertex, PositionValue, EdgeProp > ResultType;
    using std::back_inserter;
    
    std::size_t i = 0;
    while(true) {
      PositionValue p_pred = get_sample(*space);
      PositionValue p_succ = p_pred;
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p_pred, back_inserter(Nc), g, *space, g_position);
      
      Vertex x_pred, x_succ; 
      bool was_expanded, was_retracted; 
      EdgeProp ep_pred, ep_succ;
      boost::tie(x_pred, was_expanded, ep_pred)  = detail::expand_to_nearest(p_pred, Nc, g, vis, predecessor);
      boost::tie(x_succ, was_retracted, ep_succ) = detail::retract_from_nearest(p_succ, Nc, g, vis, successor);
      if( was_expanded || was_retracted )
        return ResultType(x_pred, p_pred, ep_pred, x_succ, p_succ, ep_succ);
      if( i >= 10 )
        return ResultType(boost::graph_traits<Graph>::null_vertex(), p_pred, ep_pred, 
                          boost::graph_traits<Graph>::null_vertex(), p_succ, ep_succ);
      ++i;
    };
    
    return ResultType();
  };
  
  template <typename Graph,
            typename RRGVisitor,
            typename PositionMap>
  typename boost::enable_if< boost::is_directed_graph<Graph>,
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled,
                typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, PositionValue, EdgeProp, Vertex, PositionValue, EdgeProp > ResultType;
    using std::back_inserter;
    
    std::size_t i = 0;
    while(true) {
      PositionValue p_pred = get_sample(*space);
      PositionValue p_succ = p_pred;
      
      std::vector<Vertex> Pred, Succ;
      select_neighborhood(p_pred, back_inserter(Pred), back_inserter(Succ), g, *space, g_position);
      
      Vertex x_pred, x_succ; 
      bool was_expanded, was_retracted; 
      EdgeProp ep_pred, ep_succ;
      boost::tie(x_pred, was_expanded, ep_pred)  = detail::expand_to_nearest(p_pred, Pred, g, vis, predecessor);
      boost::tie(x_succ, was_retracted, ep_succ) = detail::retract_from_nearest(p_succ, Succ, g, vis, successor);
      if( was_expanded || was_retracted )
        return ResultType(x_pred, p_pred, ep_pred, x_succ, p_succ, ep_succ);
      if( i >= 10 )
        return ResultType(boost::graph_traits<Graph>::null_vertex(), p_pred, ep_pred, 
                          boost::graph_traits<Graph>::null_vertex(), p_succ, ep_succ);
      ++i;
    };
  };

};
  

};

};


#endif

