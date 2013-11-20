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
#include "lazy_connector.hpp"
#include "branch_and_bound_connector.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Lazy-SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_lazy_sbastar_no_init(const SBAStarBundle& bdl) {
    
    detail::generate_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, lazy_node_connector(), 
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
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Lazy-SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_lazy_bnb_sbastar_no_init(const SBAStarBundle& bdl, typename SBAStarBundle::vertex_type goal_vertex) {
    
    if( goal_vertex == boost::graph_traits<typename SBAStarBundle::graph_type>::null_vertex() ) {
      
      detail::generate_sbastar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, lazy_node_connector(), 
        bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
        bdl.m_distance, bdl.m_predecessor, bdl.m_key, bdl.m_select_neighborhood);
      
    } else {
      
      detail::generate_sbastar_no_init_impl(
        *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
        branch_and_bound_connector<typename SBAStarBundle::graph_type>(
          *(bdl.m_g),
          bdl.m_start_vertex, 
          goal_vertex
        ), 
        bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
        bdl.m_distance, bdl.m_predecessor, bdl.m_key, bdl.m_select_neighborhood);
      
    };
  };

  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_lazy_bnb_sbastar(const SBAStarBundle& bdl, typename SBAStarBundle::vertex_type goal_vertex) {
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_lazy_bnb_sbastar_no_init(bdl, goal_vertex);
    
  };
  

};

};

#endif
















