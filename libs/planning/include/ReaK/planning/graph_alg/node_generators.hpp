/**
 * \file node_generators.hpp
 * 
 * This library contains a node generator that can be used in RRT-like algorithms. 
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
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

#ifndef REAK_NODE_GENERATORS_HPP
#define REAK_NODE_GENERATORS_HPP

#include <utility>
#include <vector>
#include <iterator>
#include <boost/tuple/tuple.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/random_sampler_concept.hpp>

namespace ReaK {
  
namespace graph {

namespace detail {
  
  template <typename Graph>
  struct rrg_node_puller {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, bool, EdgeProp > ResultType;
    
    template <typename PositionValue, typename RRGVisitor>
    static
    ResultType expand_to_nearest(PositionValue& p_new, const std::vector< Vertex >& Nc,
                                 Graph& g, RRGVisitor vis) {
      PositionValue p_tmp; 
      bool expand_succeeded = false;
      EdgeProp ep;
      for(typename std::vector<Vertex>::const_iterator it = Nc.begin(); it != Nc.end(); ++it) {
        boost::tie(p_tmp, expand_succeeded, ep) = vis.steer_towards_position(p_new, *it, g);
        if(expand_succeeded) {
          p_new = p_tmp;
          return ResultType(*it, true, ep);
        };
      };
      return ResultType(boost::graph_traits<Graph>::null_vertex(), false, EdgeProp());
    };
    
    template <typename PositionValue, typename RRGVisitor, typename PredecessorMap>
    static
    ResultType expand_to_nearest(PositionValue& p_new, const std::vector< Vertex >& Pred,
                                 Graph& g, RRGVisitor vis, PredecessorMap predecessor) {
      PositionValue p_tmp; 
      bool expand_succeeded = false;
      EdgeProp ep;
      for(typename std::vector<Vertex>::const_iterator it = Pred.begin(); it != Pred.end(); ++it) {
        if( get(predecessor, g[*it]) != boost::graph_traits<Graph>::null_vertex() )
          boost::tie(p_tmp, expand_succeeded, ep) = vis.steer_towards_position(p_new, *it, g);
        if(expand_succeeded) {
          p_new = p_tmp;
          return ResultType(*it, true, ep);
        };
      };
      return ResultType(boost::graph_traits<Graph>::null_vertex(), false,EdgeProp());
    };
    
    template <typename PositionValue, typename RRGVisitor, typename SuccessorMap>
    static
    ResultType retract_from_nearest(PositionValue& p_new, const std::vector< Vertex >& Succ,
                                    Graph& g, RRGVisitor vis, SuccessorMap successor) {
      PositionValue p_tmp;
      bool retract_succeeded = false;
      EdgeProp ep;
      for(typename std::vector<Vertex>::const_iterator it = Succ.begin(); it != Succ.end(); ++it) {
        if( get(successor, g[*it]) != boost::graph_traits<Graph>::null_vertex() )
          boost::tie(p_tmp, retract_succeeded, ep) = vis.steer_back_to_position(p_new, *it, g);
        if(retract_succeeded) {
          p_new = p_tmp;
          return ResultType(*it, true, ep);
        };
      };
      return ResultType(boost::graph_traits<Graph>::null_vertex(), false, EdgeProp());
    };
    
  };
  
};
  
/**
 * This node generator that can be used in RRT-like algorithms. 
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
 * \tparam Topology The topology type on which the planning is performed (i.e., the configuration space type).
 * \tparam RandomSampler The type of the random-sampler that can generate random points in the configuration space, should model ReaK::pp::RandomSamplerConcept.
 * \tparam NcSelector The type of a functor that can be used to perform nearest-neighbor queries.
 */
template <typename Topology, typename RandomSampler, typename NcSelector>
struct rrg_node_generator {
  
  const Topology* space;
  RandomSampler get_sample;
  NcSelector select_neighborhood;
  
  rrg_node_generator(const Topology* aSpace,
                     RandomSampler aGetSample,
                     NcSelector aSelectNeighborhood) : 
                     space(aSpace),
                     get_sample(aGetSample),
                     select_neighborhood(aSelectNeighborhood) { };
  
  template <typename Graph,
            typename RRGVisitor,
            typename PositionMap>
  typename boost::enable_if< boost::is_undirected_graph<Graph>,
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, PositionValue, EdgeProp > ResultType;
    typedef detail::rrg_node_puller<Graph> NodePuller;
    using std::back_inserter;
    
    std::size_t i = 0;
    while(true) {
      PositionValue p_new = get_sample(*space);
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p_new, back_inserter(Nc), g, *space, g_position);
      
      Vertex x_near; bool was_expanded; EdgeProp ep;
      boost::tie(x_near, was_expanded, ep) = NodePuller::expand_to_nearest(p_new, Nc, g, vis);
      if( was_expanded )
        return ResultType(x_near, p_new, ep);
      if( i >= 10 )
        return ResultType(boost::graph_traits<Graph>::null_vertex(), p_new, ep);
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
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, PositionValue, EdgeProp > ResultType;
    typedef detail::rrg_node_puller<Graph> NodePuller;
    using std::back_inserter;
    
    std::size_t i = 0;
    while(true) {
      PositionValue p_new = get_sample(*space);
      
      std::vector<Vertex> Pred, Succ;
      select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g, *space, g_position);
      
      Vertex x_near; bool was_expanded; EdgeProp ep;
      boost::tie(x_near, was_expanded, ep) = NodePuller::expand_to_nearest(p_new, Pred, g, vis);
      if( was_expanded )
        return ResultType(x_near, p_new, ep);
      if( i >= 10 )
        return ResultType(boost::graph_traits<Graph>::null_vertex(), p_new, ep);
      ++i;
    };
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
    typedef detail::rrg_node_puller<Graph> NodePuller;
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
      boost::tie(x_pred, was_expanded, ep_pred)  = NodePuller::expand_to_nearest(p_pred, Nc, g, vis, predecessor);
      boost::tie(x_succ, was_retracted, ep_succ) = NodePuller::retract_from_nearest(p_succ, Nc, g, vis, successor);
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
    typedef detail::rrg_node_puller<Graph> NodePuller;
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
      boost::tie(x_pred, was_expanded, ep_pred)  = NodePuller::expand_to_nearest(p_pred, Pred, g, vis, predecessor);
      boost::tie(x_succ, was_retracted, ep_succ) = NodePuller::retract_from_nearest(p_succ, Succ, g, vis, successor);
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

