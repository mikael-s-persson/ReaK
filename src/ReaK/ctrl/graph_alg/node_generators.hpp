/**
 * \file node_generators.hpp
 * 
 * This library contains
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

#include <functional>
#include <vector>
#include <iterator>
#include <boost/limits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"

#include "rr_tree.hpp"

namespace ReaK {
  
namespace graph {


namespace detail {
  
  
  
  template <typename PositionValue,
            typename Graph,
            typename RRGVisitor>
  inline bool expand_to_nearest(
      typename boost::graph_traits<Graph>::vertex_descriptor& x_near,
      PositionValue& p_new,
      const std::vector< typename boost::graph_traits<Graph>::vertex_descriptor >& Nc,
      Graph& g,
      RRGVisitor vis) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; 
    using std::back_inserter;

    if(Nc.empty())
      return false;
    
    PositionValue p_tmp; bool expand_succeeded = false;
    typename std::vector<Vertex>::const_iterator it = Nc.begin();
    while((!expand_succeeded) && (it != Nc.end())) {
      boost::tie(p_tmp, expand_succeeded) = vis.steer_towards_position(p_new, *it, g);
      ++it;
    };
    
    if(!expand_succeeded)
      return false;
    
    x_near = *(--it);
    p_new = p_tmp;
    return true;
  };
  
  
  
  
};
  

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
  typename boost::property_traits<PositionMap>::value_type >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    using std::back_inserter;
    
    while(true) {
      PositionValue p_new = get_sample(space);
      
      std::vector<Vertex> Nc; 
      select_neighborhood(p_new, back_inserter(Nc), g, space, g_position);
      
      Vertex x_near;
      if( detail::expand_to_nearest(x_near, p_new, Nc, g, vis) )
        return p_new;
    };
    
    return PositionValue();
  };
  
  template <typename Graph,
            typename RRGVisitor,
            typename PositionMap>
  inline typename boost::enable_if< boost::is_directed_graph<Graph>,
  typename boost::property_traits<PositionMap>::value_type >::type 
    operator()(Graph& g,
               RRGVisitor vis,
               PositionMap g_position) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    using std::back_inserter;
    
    while(true) {
      PositionValue p_new = get_sample(space);
      
      std::vector<Vertex> Pred;
      std::vector<Vertex> Succ;
      select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g, space, g_position);
      
      Vertex x_near;
      if( detail::expand_to_nearest(x_near, p_new, Pred, g, vis) )
        return p_new;
    };
    
    return PositionValue();
  };

};
  

};

};


#endif

