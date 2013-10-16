/**
 * \file solution_path_factories.hpp
 * 
 * This library contains implementations details for path-planning queries to use for constructing 
 * solution paths or trajectories.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_SOLUTION_PATH_FACTORIES_HPP
#define REAK_SOLUTION_PATH_FACTORIES_HPP

#include "base/defs.hpp"

#include "metric_space_concept.hpp"
#include "subspace_concept.hpp"
#include "steerable_space_concept.hpp"

#include "any_motion_graphs.hpp"

#include <boost/utility/enable_if.hpp>

#include <map>

namespace ReaK {
  
namespace pp {



namespace detail {
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_basic_solution_path_impl(FreeSpaceType& space, 
                                      const graph::any_graph& g, 
                                      graph::any_graph::vertex_descriptor start_node, 
                                      graph::any_graph::vertex_descriptor goal_node,
                                      const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                      double goal_distance,
                                      std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    typedef graph::any_graph::edge_descriptor Edge;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >      position     = graph::get_dyn_prop<const point_type&>("vertex_position", g);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g);
    
    double solutions_total_dist = goal_distance;
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(goal_distance > 0.0) {
      std::pair< point_type, SteerRecordType > goal_steer_result = space.steer_position_toward(position[goal_node], 1.0, goal_pos);
      waypoints.push_front(goal_pos);
      for(SteerIter it = goal_steer_result.second.end_fraction_travel(); it != goal_steer_result.second.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
    };
    
    waypoints.push_front(position[goal_node]);
    
    while((in_degree(goal_node, g)) && (!g.equal_descriptors(goal_node, start_node))) {
      Edge e = *(in_edges(goal_node, g).first);
      Vertex u = source(e, g);
      const SteerRecordType& sr = steer_record[e];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position[u], position[goal_node], *sup_space_ptr);
      goal_node = u;
      waypoints.push_front(position[goal_node]);
    };
    
    if( g.equal_descriptors(goal_node, start_node) && ((solutions.empty()) || (solutions_total_dist < solutions.begin()->first))) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_basic_solution_path_impl(FreeSpaceType& space, 
                                      const graph::any_graph& g, 
                                      graph::any_graph::vertex_descriptor start_node, 
                                      graph::any_graph::vertex_descriptor goal_node,
                                      const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                      double goal_distance,
                                      std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type > position = graph::get_dyn_prop<const point_type&>("vertex_position", g);
    
    double solutions_total_dist = goal_distance;
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(goal_distance > 0.0)
      waypoints.push_front(goal_pos);
    
    waypoints.push_front(position[goal_node]);
    
    while(in_degree(goal_node, g) && (!g.equal_descriptors(goal_node, start_node))) {
      Vertex v = source(*(in_edges(goal_node, g).first), g);
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position[v], position[goal_node], *sup_space_ptr);
      waypoints.push_front(position[v]);
      goal_node = v;
    };
    
    if( g.equal_descriptors(goal_node, start_node) && ((solutions.empty()) || (solutions_total_dist < solutions.begin()->first))) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_optimal_solution_path_impl(FreeSpaceType& space, 
                                        const graph::any_graph& g, 
                                        graph::any_graph::vertex_descriptor start_node, 
                                        graph::any_graph::vertex_descriptor goal_node,
                                        const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                        double goal_distance,
                                        std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    typedef graph::any_graph::in_edge_iterator InEdgeIter;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >      position       = graph::get_dyn_prop<const point_type&>("vertex_position", g);
    graph::any_graph::property_map_by_ptr< const std::size_t >     predecessor    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g);
    graph::any_graph::property_map_by_ptr< const double >          distance_accum = graph::get_dyn_prop<const double&>("vertex_distance_accum", g);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record   = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g);
    
    double solutions_total_dist = distance_accum[goal_node] + goal_distance;
    
    if( ! (solutions_total_dist < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (solutions_total_dist >= solutions.begin()->first) ) )
      return SolutionRecPtr();
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(goal_distance > 0.0) {
      std::pair< point_type, SteerRecordType > goal_steer_result = space.steer_position_toward(position[goal_node], 1.0, goal_pos);
      waypoints.push_front(goal_pos);
      for(SteerIter it = goal_steer_result.second.end_fraction_travel(); it != goal_steer_result.second.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
    };
    
    waypoints.push_front(position[goal_node]);
    
    while(!g.equal_descriptors(goal_node, start_node)) {
      std::pair<InEdgeIter,InEdgeIter> er = in_edges(goal_node, g);
      goal_node = Vertex( boost::any( predecessor[goal_node] ) ); 
      while( ( er.first != er.second ) && ( !g.equal_descriptors(goal_node, source(*(er.first), g)) ) )
        ++(er.first);
      if(er.first == er.second)
        break;
      const SteerRecordType& sr = steer_record[*(er.first)];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      waypoints.push_front(position[goal_node]);
    };
    
    if(g.equal_descriptors(goal_node, start_node)) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_optimal_solution_path_impl(FreeSpaceType& space, 
                                        const graph::any_graph& g, 
                                        graph::any_graph::vertex_descriptor start_node, 
                                        graph::any_graph::vertex_descriptor goal_node,
                                        const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                        double goal_distance,
                                        std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >  position       = graph::get_dyn_prop<const point_type&>("vertex_position", g);
    graph::any_graph::property_map_by_ptr< const std::size_t > predecessor    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g);
    graph::any_graph::property_map_by_ptr< const double >      distance_accum = graph::get_dyn_prop<const double&>("vertex_distance_accum", g);
    
    double solutions_total_dist = distance_accum[goal_node] + goal_distance;
    
    if( ! (solutions_total_dist < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (solutions_total_dist >= solutions.begin()->first) ) )
      return SolutionRecPtr();
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(goal_distance > 0.0)
      waypoints.push_front(goal_pos);
    
    waypoints.push_front(position[goal_node]);
    
    while(!g.equal_descriptors(goal_node, start_node)) {
      goal_node = Vertex( boost::any( predecessor[goal_node] ) ); 
      waypoints.push_front( position[goal_node] );
    };
    
    if(g.equal_descriptors(goal_node, start_node)) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
  
  
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_basic_solution_path_impl(FreeSpaceType& space, 
                                      const graph::any_graph& g1, 
                                      const graph::any_graph& g2, 
                                      graph::any_graph::vertex_descriptor start_node, 
                                      graph::any_graph::vertex_descriptor goal_node, 
                                      graph::any_graph::vertex_descriptor join1_node, 
                                      graph::any_graph::vertex_descriptor join2_node,
                                      double joining_distance,
                                      std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    typedef graph::any_graph::edge_descriptor Edge;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >      position1     = graph::get_dyn_prop<const point_type&>("vertex_position", g1);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record1 = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g1);
    graph::any_graph::property_map_by_ptr< const point_type >      position2     = graph::get_dyn_prop<const point_type&>("vertex_position", g2);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record2 = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g2);
    
    double solutions_total_dist = joining_distance;
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(joining_distance > 0.0) {
      std::pair< point_type, SteerRecordType > join_steer_result = space.steer_position_toward(position1[join1_node], 1.0, position2[join2_node]);
      for(SteerIter it = join_steer_result.second.end_fraction_travel(); it != join_steer_result.second.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
    };
    
    waypoints.push_front(position1[join1_node]);
    
    while((in_degree(join1_node, g1)) && (!g1.equal_descriptors(join1_node, start_node))) {
      Edge e = *(in_edges(join1_node, g1).first);
      Vertex u = source(e, g1);
      const SteerRecordType& sr = steer_record1[e];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position1[u], position1[join1_node], *sup_space_ptr);
      join1_node = u;
      waypoints.push_front(position1[join1_node]);
    };
    
    waypoints.push_back(position2[join2_node]);
    
    while((in_degree(join2_node, g2)) && (!g2.equal_descriptors(join2_node, goal_node))) {
      Edge e = *(in_edges(join2_node, g2).first);
      Vertex u = source(e, g2);
      const SteerRecordType& sr = steer_record2[e];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_back( point_type(*it) );
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position2[u], position2[join2_node], *sup_space_ptr);
      join2_node = u;
      waypoints.push_back(position2[join2_node]);
    };
    
    
    if( g1.equal_descriptors(join1_node, start_node) && 
        g2.equal_descriptors(join2_node, goal_node) && 
        ( (solutions.empty()) || (solutions_total_dist < solutions.begin()->first) ) ) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_basic_solution_path_impl(FreeSpaceType& space, 
                                      const graph::any_graph& g1, 
                                      const graph::any_graph& g2, 
                                      graph::any_graph::vertex_descriptor start_node, 
                                      graph::any_graph::vertex_descriptor goal_node, 
                                      graph::any_graph::vertex_descriptor join1_node, 
                                      graph::any_graph::vertex_descriptor join2_node,
                                      double joining_distance,
                                      std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type > position1 = graph::get_dyn_prop<const point_type&>("vertex_position", g1);
    graph::any_graph::property_map_by_ptr< const point_type > position2 = graph::get_dyn_prop<const point_type&>("vertex_position", g2);
    
    double solutions_total_dist = joining_distance;
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    waypoints.push_front(position1[join1_node]);
    
    while(in_degree(join1_node, g1) && (!g1.equal_descriptors(join1_node, start_node))) {
      Vertex v = source(*(in_edges(join1_node, g1).first), g1);
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position1[v], position1[join1_node], *sup_space_ptr);
      waypoints.push_front(position1[v]);
      join1_node = v;
    };
    
    waypoints.push_back(position2[join2_node]);
    
    while(in_degree(join2_node, g2) && (!g2.equal_descriptors(join2_node, goal_node))) {
      Vertex v = source(*(in_edges(join2_node, g2).first), g2);
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(position2[join2_node], position2[v], *sup_space_ptr);
      waypoints.push_back(position2[v]);
      join2_node = v;
    };
    
    if( g1.equal_descriptors(join1_node, start_node) && 
        g2.equal_descriptors(join2_node, goal_node) && 
        ( (solutions.empty()) || (solutions_total_dist < solutions.begin()->first) ) ) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_optimal_solution_path_impl(FreeSpaceType& space, 
                                        const graph::any_graph& g1, 
                                        const graph::any_graph& g2, 
                                        graph::any_graph::vertex_descriptor start_node, 
                                        graph::any_graph::vertex_descriptor goal_node, 
                                        graph::any_graph::vertex_descriptor join1_node, 
                                        graph::any_graph::vertex_descriptor join2_node,
                                        double joining_distance,
                                        std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    typedef graph::any_graph::in_edge_iterator InEdgeIter;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >      position1       = graph::get_dyn_prop<const point_type&>("vertex_position", g1);
    graph::any_graph::property_map_by_ptr< const std::size_t >     predecessor1    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g1);
    graph::any_graph::property_map_by_ptr< const double >          distance_accum1 = graph::get_dyn_prop<const double&>("vertex_distance_accum", g1);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record1   = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g1);
    graph::any_graph::property_map_by_ptr< const point_type >      position2       = graph::get_dyn_prop<const point_type&>("vertex_position", g2);
    graph::any_graph::property_map_by_ptr< const std::size_t >     predecessor2    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g2);
    graph::any_graph::property_map_by_ptr< const double >          distance_accum2 = graph::get_dyn_prop<const double&>("vertex_distance_accum", g2);
    graph::any_graph::property_map_by_ptr< const SteerRecordType > steer_record2   = graph::get_dyn_prop<const SteerRecordType&>("edge_steer_record", g2);
    
    double solutions_total_dist = distance_accum1[join1_node] + distance_accum2[join2_node] + joining_distance;
    
    if( ! (solutions_total_dist < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (solutions_total_dist >= solutions.begin()->first) ) )
      return SolutionRecPtr();
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    if(joining_distance > 0.0) {
      std::pair< point_type, SteerRecordType > join_steer_result = space.steer_position_toward(position1[join1_node], 1.0, position2[join2_node]);
      for(SteerIter it = join_steer_result.second.end_fraction_travel(); it != join_steer_result.second.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
    };
    
    waypoints.push_front(position1[join1_node]);
    
    while(!g1.equal_descriptors(join1_node, start_node)) {
      std::pair<InEdgeIter,InEdgeIter> er = in_edges(join1_node, g1);
      join1_node = Vertex( boost::any( predecessor1[join1_node] ) ); 
      while( ( er.first != er.second ) && ( !g1.equal_descriptors(join1_node, source(*(er.first), g1)) ) )
        ++(er.first);
      if(er.first == er.second)
        break;
      const SteerRecordType& sr = steer_record1[*(er.first)];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      waypoints.push_front(position1[join1_node]);
    };
    
    waypoints.push_back(position2[join2_node]);
    
    while(!g2.equal_descriptors(join2_node, goal_node)) {
      std::pair<InEdgeIter,InEdgeIter> er = in_edges(join2_node, g2);
      join2_node = Vertex( boost::any( predecessor2[join2_node] ) ); 
      while( ( er.first != er.second ) && ( !g2.equal_descriptors(join2_node, source(*(er.first), g2)) ) )
        ++(er.first);
      if(er.first == er.second)
        break;
      const SteerRecordType& sr = steer_record2[*(er.first)];
      for(SteerIter it = sr.end_fraction_travel(); it != sr.begin_fraction_travel(); it -= 1.0)
        waypoints.push_back( point_type(*it) );
      waypoints.push_back(position2[join2_node]);
    };
    
    
    if(g1.equal_descriptors(join1_node, start_node) && g2.equal_descriptors(join2_node, goal_node)) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
  template <typename SolutionWrapperType, typename FreeSpaceType, typename SolutionRecPtr>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  SolutionRecPtr >::type 
    register_optimal_solution_path_impl(FreeSpaceType& space, 
                                        const graph::any_graph& g1, 
                                        const graph::any_graph& g2, 
                                        graph::any_graph::vertex_descriptor start_node, 
                                        graph::any_graph::vertex_descriptor goal_node, 
                                        graph::any_graph::vertex_descriptor join1_node, 
                                        graph::any_graph::vertex_descriptor join2_node,
                                        double joining_distance,
                                        std::map<double, SolutionRecPtr >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef graph::any_graph::vertex_descriptor Vertex;
    
    typedef typename SolutionWrapperType::wrapped_type SolutionPathType;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    graph::any_graph::property_map_by_ptr< const point_type >  position1       = graph::get_dyn_prop<const point_type&>("vertex_position", g1);
    graph::any_graph::property_map_by_ptr< const std::size_t > predecessor1    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g1);
    graph::any_graph::property_map_by_ptr< const double >      distance_accum1 = graph::get_dyn_prop<const double&>("vertex_distance_accum", g1);
    graph::any_graph::property_map_by_ptr< const point_type >  position2       = graph::get_dyn_prop<const point_type&>("vertex_position", g2);
    graph::any_graph::property_map_by_ptr< const std::size_t > predecessor2    = graph::get_dyn_prop<const std::size_t&>("vertex_predecessor", g2);
    graph::any_graph::property_map_by_ptr< const double >      distance_accum2 = graph::get_dyn_prop<const double&>("vertex_distance_accum", g2);
    
    double solutions_total_dist = distance_accum1[join1_node] + distance_accum2[join2_node] + joining_distance;
    
    if( ! (solutions_total_dist < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (solutions_total_dist >= solutions.begin()->first) ) )
      return SolutionRecPtr();
    
    shared_ptr< SolutionWrapperType > new_sol(new SolutionWrapperType("planning_solution", SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr))));
    SolutionPathType& waypoints = new_sol->get_wrapped_object();
    
    waypoints.push_front(position1[join1_node]);
    
    while(!g1.equal_descriptors(join1_node, start_node)) {
      join1_node = Vertex( boost::any( predecessor1[join1_node] ) ); 
      waypoints.push_front( position1[join1_node] );
    };
    
    waypoints.push_back(position2[join2_node]);
    
    while(!g2.equal_descriptors(join2_node, goal_node)) {
      join2_node = Vertex( boost::any( predecessor2[join2_node] ) ); 
      waypoints.push_back( position2[join2_node] );
    };
    
    if(g1.equal_descriptors(join1_node, start_node) && g2.equal_descriptors(join2_node, goal_node)) {
      solutions[solutions_total_dist] = new_sol;
      return new_sol;
    };
    
    return SolutionRecPtr();
  };
  
  
  
};



};

};

#endif

