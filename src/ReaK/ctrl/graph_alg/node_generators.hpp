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

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"

namespace ReaK {
  
namespace graph {

namespace detail {
  
  template <typename PositionValue,
            typename Graph,
            typename RRGVisitor>
  inline 
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor, 
                bool, typename Graph::edge_bundled > 
    expand_to_nearest(
      PositionValue& p_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Nc,
      Graph& g,
      RRGVisitor vis) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< Vertex, bool, EdgeProp > ResultType;
    
    if(Nc.empty())
      return false;
    
    PositionValue p_tmp; 
    bool expand_succeeded = false;
    EdgeProp ep;
    typename std::vector<Vertex>::const_iterator it = Nc.begin();
    while((!expand_succeeded) && (it != Nc.end())) {
      boost::tie(p_tmp, expand_succeeded, ep) = vis.steer_towards_position(p_new, *it, g);
      ++it;
    };
    
    if(!expand_succeeded)
      return ResultType(Vertex(), false,EdgeProp());
    
    p_new = p_tmp;
    return ResultType(*(--it), true, ep);
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
  inline typename boost::enable_if< boost::is_undirected_graph<Graph>,
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                          typename boost::property_traits<PositionMap>::value_type,
                          typename Graph::edge_bundled > ResultType;
    using std::back_inserter;
    
    while(true) {
      PositionValue p_new = get_sample(*space);
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p_new, back_inserter(Nc), g, *space, g_position);
      
      Vertex x_near; bool was_expanded; EdgeProp ep;
      boost::tie(x_near, was_expanded, ep) = detail::expand_to_nearest(p_new, Nc, g, vis);
      if( was_expanded )
        return ResultType(x_near, p_new, ep);
    };
    
    return ResultType();
  };
  
  template <typename Graph,
            typename RRGVisitor,
            typename PositionMap>
  inline typename boost::enable_if< boost::is_directed_graph<Graph>,
  boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                typename boost::property_traits<PositionMap>::value_type,
                typename Graph::edge_bundled > >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) const {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                          typename boost::property_traits<PositionMap>::value_type,
                          typename Graph::edge_bundled > ResultType;
    using std::back_inserter;
    
    while(true) {
      PositionValue p_new = get_sample(*space);
      
      std::vector<Vertex> Pred, Succ;
      select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g, *space, g_position);
      
      Vertex x_near; bool was_expanded; EdgeProp ep;
      boost::tie(x_near, was_expanded, ep) = detail::expand_to_nearest(p_new, Pred, g, vis);
      if( was_expanded )
        return ResultType(x_near, p_new, ep);
    };
    
    return ResultType();
  };

};
  

};

};


#endif

