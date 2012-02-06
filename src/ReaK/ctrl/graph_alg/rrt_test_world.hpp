/**
 * \file rrt_test_world.hpp
 * 
 * This library defines an example class that uses the RRT algorithm on a world map. This class 
 * basically just wraps the algorithms found in "rr_tree.hpp" for a simple 2D path-planning problem.
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free). It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 * See the test_rrt.cpp file for a program that uses this class.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2011
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

#ifndef REAK_RRT_TEST_WORLD_HPP
#define REAK_RRT_TEST_WORLD_HPP

#include <iostream>
#include <iomanip>

#include "rr_tree.hpp"
#include "path_planning/topological_search.hpp"
#include "path_planning/metric_space_search.hpp"
#include <cmath>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>


#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/default_random_sampler.hpp"

#include "topologies/hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"


namespace boost {

  enum vertex_position_t { vertex_position };

  BOOST_INSTALL_PROPERTY(vertex, position);

};


/**
 * This class is used to test the RRT algorithm on a world map. This class basically just 
 * wraps the algorithms found in "rr_tree.hpp" for a simple 2D path-planning problem.
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free). It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 * See the test_rrt.cpp file for a program that uses this class.
 */
class rrt_test_world {
  public:
    typedef ReaK::pp::hyperbox_topology< ReaK::vect<int,2> > space_type;
    typedef ReaK::pp::topology_traits< space_type >::point_type point_type;
    typedef ReaK::pp::topology_traits< space_type >::point_difference_type point_difference_type;
    
    typedef ReaK::pp::default_distance_metric distance_metric_type;
    typedef ReaK::pp::default_random_sampler random_sampler_type;
    
    typedef point_type pixel_coord;
    typedef point_difference_type pixel_difference;

    typedef boost::property< boost::vertex_position_t, pixel_coord,
            boost::property< boost::vertex_distance_t, double, boost::no_property > > WorldGridVertexProperties;

    typedef boost::no_property WorldGridEdgeProperties;

    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                            WorldGridVertexProperties,
	  	            WorldGridEdgeProperties,
		            boost::vecS> WorldGridType;

    typedef boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::bidirectionalS,boost::vecS>::vertex_descriptor VertexType;
    typedef boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::bidirectionalS,boost::vecS>::edge_descriptor EdgeType;
    typedef boost::graph_traits<WorldGridType>::in_edge_iterator EdgeIter;
    
    typedef boost::function< void(const cv::Mat&, unsigned int) > ProgressCallback;
    typedef boost::function< void(const cv::Mat&, unsigned int, unsigned int, double) > PathFoundCallback;

  private:
    cv::Mat world_map_image;
    cv::Mat world_map_output;
    int grid_width;
    int grid_height;
    int bpp;

    WorldGridType grid;
    WorldGridType grid_goal;
    VertexType current_pos;
    VertexType start_node;
    VertexType goal_node;
    pixel_coord goal_pos;
    double goal_distance;
    int best_solution;
    double max_edge_length;
    double robot_radius;
    unsigned int max_vertex_count;
    unsigned int max_num_results;
    bool unidirectional;
    unsigned int nn_search_divider;
    std::vector< boost::tuples::tuple< double, VertexType, VertexType > > solutions;

    space_type m_space;
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position;
    boost::property_map<WorldGridType, boost::vertex_distance_t>::type m_distance;
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position_goal;
    boost::property_map<WorldGridType, boost::vertex_distance_t>::type m_distance_goal;
    
    typedef ReaK::pp::dvp_tree<VertexType, space_type, 
			       boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
			       4> WorldPartition;
    WorldPartition space_part;
    WorldPartition space_part_goal;
    
    ProgressCallback progress_call_back;
    PathFoundCallback path_found_call_back;

    bool is_free(const pixel_coord& p) {
      if((p[0] < 0) || (p[0] >= grid_width) || (p[1] < 0) || (p[1] >= grid_height)) 
	return true;
      uchar* color_bits = world_map_image.ptr(int(p[1]));
      color_bits += bpp * int(p[0]);
      if( (color_bits[0] == color_bits[1]) &&
	  (color_bits[0] == color_bits[2]) &&
	  (color_bits[0] < 250) ) {
	return false;
      } else {
        return true;
      };
    };

    void draw_edge(const pixel_coord& p_u, const pixel_coord& p_v, bool goal_path = false) {
      pixel_difference diff = m_space.difference(p_v, p_u);
      double dist = get(ReaK::pp::distance_metric, m_space)(diff, m_space);
      double d = 0.0;
      while(d <= dist) {
	pixel_coord p = m_space.adjust(p_u, diff * (d / dist));
	if(p[0] < 0) p[0] = 0;
	if(p[1] < 0) p[1] = 0;
	if(p[0] >= grid_width) p[0] = grid_width-1;
	if(p[1] >= grid_height) p[1] = grid_height-1;
	uchar* color_bits = world_map_output.ptr(int(p[1]));
        color_bits += bpp * int(p[0]);
	if(goal_path) {
	  color_bits[2] = 255;
	  color_bits[1] = 0;
	  color_bits[0] = 0; //red color
	} else {
	  color_bits[2] = 255;
	  color_bits[1] = 140;
	  color_bits[0] = 0; //orange color
	};
	d += 1.0;
      };
    };

  public:
    void get_best_solution(std::list< pixel_coord >& path) {
      path.clear();
      if(unidirectional) {
	VertexType v = goal_node;
	path.push_front(get(m_position,v));
	VertexType u = current_pos;
	path.push_front(get(m_position,u));
	while(in_degree(u,grid)) {
	  v = u;
	  u = source(*(in_edges(v,grid).first),grid);
	  path.push_front(get(m_position,u));
	};
      } else {
	if(best_solution < 0)
	  return;
	VertexType u = solutions[best_solution].get<1>();
	path.push_front(get(m_position,u));
	while(in_degree(u,grid)) {
	  u = source(*(in_edges(u,grid).first),grid);
	  path.push_front(get(m_position,u));
	};
	u = solutions[best_solution].get<2>();
	path.push_back(get(m_position_goal,u));
	while(in_degree(u,grid_goal)) {
	  u = source(*(in_edges(u,grid_goal).first),grid_goal);
	  path.push_back(get(m_position_goal,u));
	};
      };
    };

    void vertex_added(VertexType u, WorldGridType& g) {
      if(nn_search_divider == 1) {
	if(&g == &grid)
	  space_part.insert(u);
	else if(&g == &grid_goal)
	  space_part_goal.insert(u);
      };
    };

    void edge_added(EdgeType e, WorldGridType& g) {
      VertexType u = source(e,g);
      VertexType v = target(e,g);
      pixel_coord p_u = ( &g == &grid ? get(m_position, u) : get(m_position_goal,u));
      pixel_coord p_v = ( &g == &grid ? get(m_position, v) : get(m_position_goal,v));
      if(&g == &grid)
        put(m_distance, v, get(m_distance, u) + get(ReaK::pp::distance_metric, m_space)(p_u, p_v, m_space));
      else
        put(m_distance_goal, v, get(m_distance_goal, u) + get(ReaK::pp::distance_metric, m_space)(p_u, p_v, m_space));

      draw_edge(p_u, p_v);

      if((((num_vertices(grid) + num_vertices(grid_goal)) % 100) == 0) && (progress_call_back))
	progress_call_back(world_map_output, num_vertices(grid) + num_vertices(grid_goal));

      //now, check if v is connected with the goal on a straight collision-free path.
      if(&g == &grid) {
      pixel_difference diff = m_space.difference(goal_pos, p_v);
      double dist = get(ReaK::pp::distance_metric, m_space)(diff, m_space);
      if(dist <= max_edge_length) {
	double d = 1.0;
	pixel_coord p_g = m_space.adjust(p_v, diff * (d / dist));
        while((is_free(p_g)) && (d < dist)) {
	  d += 1.0;
	  p_g = m_space.adjust(p_v, diff * (d / dist));
        };
	if(d >= dist) {
	  //we have a connection to the goal!! Is it the best so far?
	  d = get(m_distance, v) + dist;
	  if( d < goal_distance ) {
	    goal_distance = d;
	    current_pos = v;

	    //draw the path to the goal!
	    draw_edge(p_v, p_g, true);

	    while(in_degree(v, g)) {
	      u = source(*(in_edges(v,grid).first),grid);
	      draw_edge(get(m_position,u), get(m_position,v), true);
	      v = u;
	    };
	    if(path_found_call_back)
	      path_found_call_back(world_map_output, num_vertices(grid), 0, goal_distance);

	  };
	};
      };
      };
    };

    void joining_vertex_found(VertexType u, WorldGridType& g) {
      if((solutions.size() == 0) || (solutions.back().get<0>() > 0.0)) {
        if(&g == &grid)
          solutions.push_back(boost::tuples::make_tuple(-get(m_distance,u),u,goal_node));
        else
          solutions.push_back(boost::tuples::make_tuple(-get(m_distance_goal,u),start_node,u));
      } else {
        if(&g == &grid) {
          solutions.back().get<0>() = get(m_distance,u) - solutions.back().get<0>();
          solutions.back().get<1>() = u;
        } else {
          solutions.back().get<0>() = get(m_distance_goal,u) - solutions.back().get<0>();
          solutions.back().get<2>() = u;
        };
        if(solutions.back().get<0>() < goal_distance) {
          goal_distance = solutions.back().get<0>();
          best_solution = solutions.size() - 1;

          //Draw the edges of the current best solution:
          VertexType v = solutions.back().get<1>();
          while(in_degree(v, grid)) {
	    u = source(*(in_edges(v,grid).first),grid);
	    draw_edge(get(m_position,u), get(m_position,v), true);
	    v = u;
          };
          v = solutions.back().get<2>();
          while(in_degree(v, grid_goal)) {
            u = source(*(in_edges(v,grid_goal).first),grid_goal);
            draw_edge(get(m_position_goal,u),get(m_position_goal,v),true);
            v = u;
          };

	  if(path_found_call_back) {
	    path_found_call_back(world_map_output, num_vertices(grid) + num_vertices(grid_goal), best_solution, goal_distance);
  	    //std::stringstream ss;
            //ss << "test_rrt_results/" << "result_" << max_edge_length << "_" << best_solution << "_" << num_vertices(grid) + num_vertices(grid_goal) << "_" << goal_distance << ".bmp";
            //cv::imwrite(ss.str(),world_map_output);
	  };
	  
          //Redraw the edges of the current best solution as orange:
          v = solutions.back().get<1>();
          while(in_degree(v, grid)) {
	    u = source(*(in_edges(v,grid).first),grid);
	    draw_edge(get(m_position,u), get(m_position,v));
	    v = u;
          };
          v = solutions.back().get<2>();
          while(in_degree(v, grid_goal)) {
            u = source(*(in_edges(v,grid_goal).first),grid_goal);
            draw_edge(get(m_position_goal,u),get(m_position_goal,v));
            v = u;
          };
        };
      };
    };

    bool keep_going() {
      if(solutions.size() >= max_num_results)
        return false;
      return true;
    };

    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aWorldMapImage An image which represents the C-free as white (or colored) pixels and the occupied C-space as gray pixels.
     * \param aMaxEdgeLength The maximum length of an added edge, in pixel-units.
     * \param aRobotRadius The radius of the robot (collision radius), in pixel-units.
     * \param aMaxVertexCount The maximum number of vertices to add to the graph, in case no solution is found, this will stop the madness.
     * \param aMaxNumResults The maximum number of solutions to generate (a value of 1 will stop as soon as a solution is found, regardless of how sub-optimal it is).
     * \param aProgressCallback A callable object which is used to report on the progress of the RRT (called once for every 100 vertices generated).
     * \param aPathFoundCallback A callable object which is used to report that a solution which was better than any pervious solution was just found.
     * \param aUnidirectional A flag to signal whether the unidirectional (true) or bidirectional (false, default) version of the RRT generation algorithm should be used.
     * \param aNNSearchDivider The vertex number divider for the nearest-neighor search (see best_only_neighbor_search), a value of 0 will signify that the linear_neighbor_search shall be used instead of the approximate best_only_neighbor_search algorithm.
     */
    rrt_test_world(const cv::Mat& aWorldMapImage, double aMaxEdgeLength, 
		   double aRobotRadius, unsigned int aMaxVertexCount, unsigned int aMaxNumResults, 
		   ProgressCallback aProgressCallback = ProgressCallback(), 
		   PathFoundCallback aPathFoundCallback = PathFoundCallback(), 
		   bool aUnidirectional = false, unsigned int aNNSearchDivider = 10) :
                   world_map_image(aWorldMapImage.clone()), world_map_output(aWorldMapImage.clone()),
                   grid_width(aWorldMapImage.size().width),
                   grid_height(aWorldMapImage.size().height),
                   goal_distance(std::numeric_limits<double>::infinity()),
                   best_solution(-1), max_edge_length(aMaxEdgeLength), 
                   robot_radius(aRobotRadius), max_vertex_count(aMaxVertexCount),
                   max_num_results(aMaxNumResults), unidirectional(aUnidirectional),
		   nn_search_divider(aNNSearchDivider), 
                   m_space("rrt_space",pixel_coord(0,0), pixel_coord(aWorldMapImage.size().width, aWorldMapImage.size().height)),
                   m_position(get(boost::vertex_position, grid)),
                   m_distance(get(boost::vertex_distance, grid)),
                   m_position_goal(get(boost::vertex_position, grid_goal)),
                   m_distance_goal(get(boost::vertex_distance, grid_goal)),
                   space_part(grid,m_space,m_position), space_part_goal(grid_goal,m_space,m_position_goal),
                   progress_call_back(aProgressCallback), path_found_call_back(aPathFoundCallback)
    {
      
      if(world_map_image.empty()) {
	std::cout << __FILE__ << ":" << __LINE__ << " Error: The world image is empty!" << std::endl;
	throw int(0);
      };

      bpp = world_map_image.elemSize();

      for(int y = 0; y < grid_height; ++y) {
	uchar* color_bits = world_map_image.ptr(y);

	for(int x = 0; x < grid_width; ++x) {
	  if( (color_bits[2] == 0) &&
	      (color_bits[0] == 255) &&
	      (color_bits[1] == 0) ) {
	    //this is the start position.
	    start_node = current_pos = add_vertex(grid);
	    pixel_coord p; p[0] = x; p[1] = y;
	    put(m_position, current_pos, p);
	    put(m_distance, current_pos, 0.0);
	    if(nn_search_divider == 1)
              space_part.insert(start_node);
	  } else if( (color_bits[2] == 0) &&
	             (color_bits[0] == 0) &&
	             (color_bits[1] == 255) ) {
	    //this is the goal position.
	    goal_pos[0] = x;
	    goal_pos[1] = y;
	    goal_node = add_vertex(grid_goal);
	    put(m_position_goal, goal_node, goal_pos);
	    put(m_distance_goal, goal_node, 0.0);
	    if(nn_search_divider == 1)
              space_part_goal.insert(goal_node);
	  };
	  color_bits += bpp;
	};
      };
      
      if(int(std::fabs(robot_radius)) > 0) {
        cv::GaussianBlur(world_map_output,world_map_image,
		         cv::Size(int(std::fabs(robot_radius))*2+1,int(std::fabs(robot_radius))*2+1),
		         robot_radius
		        );
      };
    };

    //This is a RAII class, so the destructor is meaningless. Hurray for RAII!!
    ~rrt_test_world() { };

    void set_start_pos(const pixel_coord& aStart) {
      while(num_vertices(grid)) {
	space_part.erase(*(vertices(grid).first));
	remove_vertex(*(vertices(grid).first),grid);
      };
      start_node = current_pos = add_vertex(grid);
      put(m_position, current_pos, aStart);
      put(m_distance, current_pos, 0.0);
      if(nn_search_divider == 1)
        space_part.insert(start_node);
    };
    
    void set_goal_pos(const pixel_coord& aGoal) {
      while(num_vertices(grid_goal)) {
	space_part_goal.erase(*(vertices(grid_goal).first));
	remove_vertex(*(vertices(grid_goal).first),grid_goal);
      };
      goal_pos = aGoal;
      goal_node = add_vertex(grid_goal);
      put(m_position_goal, goal_node, aGoal);
      put(m_distance_goal, goal_node, 0.0);
      if(nn_search_divider == 1)
        space_part_goal.insert(goal_node);
    };
    
    double run() {
      unsigned int m = max_vertex_count;
      while(!(goal_distance < std::numeric_limits<double>::infinity())) {
	if(unidirectional) {
	  if(nn_search_divider == 0) { //calls the unidirectional RRT, with the linear_neighbor_search
	    ReaK::graph::generate_rrt(
	      grid, 
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
	        boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
	      m_position,
 	      ReaK::pp::linear_neighbor_search<>(),
              m, max_edge_length, 1.0);
	  } else if(nn_search_divider == 1) { //calls the unidirectional RRT, with dvp_tree as nearest-neighbor finder.
	    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition> nn_finder;
	    nn_finder.graph_tree_map[&grid] = &space_part;
	    ReaK::graph::generate_rrt(
	      grid,
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
		boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
              m_position,
	      nn_finder,
	      m, max_edge_length, 1.0);
	  } else {                      //calls the unidirectional RRT, with the best_only_neighbor_search
	    ReaK::graph::generate_rrt(
	      grid, 
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
	        boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
	      m_position,
 	      ReaK::pp::best_only_neighbor_search<>(nn_search_divider),
              m, max_edge_length, 1.0);
	  };
	} else {
	  if(nn_search_divider == 0) {    //calls the bidirectional RRT, with the linear_neighbor_search
	    ReaK::graph::generate_bidirectional_rrt(
	      grid, grid_goal,
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
	        boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
	      m_position, m_position_goal,
 	      ReaK::pp::linear_neighbor_search<>(),
              m, max_edge_length, 1.0);
	  } else if(nn_search_divider == 1) {
	    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition> nn_finder;
	    nn_finder.graph_tree_map[&grid] = &space_part;
	    nn_finder.graph_tree_map[&grid_goal] = &space_part_goal;
	    ReaK::graph::generate_bidirectional_rrt(
	      grid, grid_goal,
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
	        boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
	      m_position, m_position_goal,
 	      nn_finder,
              m, max_edge_length, 1.0);
	    
	    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
	    for(unsigned int i=0;i<100000;++i) {
	      nn_finder(m_space.random_point(),grid,m_space,m_position);
	    };
	    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
	    std::cout << "100000 queries of the vp-tree took: " << dt.total_microseconds() << " microsec on " << num_vertices(grid) << " vertices." << std::endl;
	    
	    WorldPartition fresh_partition(grid,m_space,m_position);
	    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition> nn_finder_fresh;
	    nn_finder_fresh.graph_tree_map[&grid] = &fresh_partition;
	    t_start = boost::posix_time::microsec_clock::local_time();
	    for(unsigned int i=0;i<100000;++i) {
	      nn_finder_fresh(m_space.random_point(),grid,m_space,m_position);
	    };
	    dt = boost::posix_time::microsec_clock::local_time() - t_start;
	    std::cout << "100000 queries of a fresh vp-tree took: " << dt.total_microseconds() << " microsec on " << num_vertices(grid) << " vertices." << std::endl;
	    
	    ReaK::pp::linear_neighbor_search<> lnn_finder;
	    t_start = boost::posix_time::microsec_clock::local_time();
	    for(unsigned int i=0;i<100000;++i) {
	      lnn_finder(m_space.random_point(),grid,m_space,m_position);
	    };
	    dt = boost::posix_time::microsec_clock::local_time() - t_start;
	    std::cout << "100000 queries of the linear search took: " << dt.total_microseconds() << " microsec on " << num_vertices(grid) << " vertices." << std::endl;
	    
	  } else {                         //calls the bidirectional RRT, with the best_only_neighbor_search
            ReaK::graph::generate_bidirectional_rrt(
	      grid, grid_goal,
	      m_space,
	      ReaK::graph::make_composite_rrt_visitor(
	        boost::bind(&rrt_test_world::vertex_added,this,_1,_2),
                boost::bind(&rrt_test_world::edge_added,this,_1,_2),
	        boost::bind(&rrt_test_world::is_free,this,_1),
                boost::bind(&rrt_test_world::joining_vertex_found,this,_1,_2),
                boost::bind(&rrt_test_world::keep_going,this)),
	      m_position, m_position_goal,
 	      ReaK::pp::best_only_neighbor_search<>(nn_search_divider),
              m, max_edge_length, 1.0);
	  };
	};
	m += max_vertex_count;
      };
      return goal_distance;
    };
};

#endif
