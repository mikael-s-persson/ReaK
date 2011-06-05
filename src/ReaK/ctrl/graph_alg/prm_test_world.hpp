
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

#ifndef PRM_TEST_WORLD_HPP
#define PRM_TEST_WORLD_HPP

#include <iostream>
#include <iomanip>

#include "probabilistic_roadmap.hpp"
#include "path_planning/topological_search.hpp"
#include "path_planning/metric_space_search.hpp"
#include <cmath>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/astar_search.hpp>


#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <boost/date_time/posix_time/posix_time.hpp>


namespace boost {

  enum vertex_position_t { vertex_position };
  enum vertex_density_t { vertex_density };

  BOOST_INSTALL_PROPERTY(vertex, position);  
  BOOST_INSTALL_PROPERTY(vertex, density);

};

namespace boost {

  enum vertex_heuristic_t { vertex_heuristic };
  enum vertex_rhs_t { vertex_rhs };

  BOOST_INSTALL_PROPERTY(vertex, heuristic);
  BOOST_INSTALL_PROPERTY(vertex, rhs);

};


/**
 * This class is used to test the PRM algorithm on a world map. This class basically just 
 * wraps the algorithms found in "probabilistic_roadmap.hpp" for a simple 2D path-planning problem.
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free). It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 */
class prm_test_world {
  public:
    typedef boost::rectangle_topology<>::point_type point_type;
    typedef boost::rectangle_topology<>::point_difference_type point_difference_type;

    typedef boost::property< boost::vertex_position_t, point_type, //for PRM
	    boost::property< boost::vertex_rhs_t, double,       //for A*
	    boost::property< boost::vertex_distance_t, double,  //for A*
	    boost::property< boost::vertex_density_t, double,   //for PRM
	    boost::property< boost::vertex_color_t, boost::default_color_type, //for A*
	    boost::property< boost::vertex_predecessor_t, boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::undirectedS,boost::vecS>::vertex_descriptor, //for A*
	    boost::no_property > > > > > > WorldGridVertexProperties;

    typedef boost::property< boost::edge_weight_t, double, //for A*
            boost::no_property> WorldGridEdgeProperties;

    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                   WorldGridVertexProperties,
	  	                   WorldGridEdgeProperties,
		                   boost::vecS> WorldGridType;

    typedef boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::undirectedS,boost::vecS>::vertex_descriptor VertexType;
    typedef boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::undirectedS,boost::vecS>::edge_descriptor EdgeType;
    typedef boost::graph_traits<WorldGridType>::in_edge_iterator InEdgeIter;
    typedef boost::graph_traits<WorldGridType>::out_edge_iterator OutEdgeIter;
    
    typedef boost::function< void(const cv::Mat&, unsigned int) > ProgressCallback;
    typedef boost::function< void(const cv::Mat&, unsigned int, unsigned int, double) > PathFoundCallback;

  private:
    cv::Mat world_map_image;
    cv::Mat world_map_output;
    int grid_width;
    int grid_height;
    int bpp;

    WorldGridType grid;
    VertexType current_pos;
    VertexType start_node;
    VertexType goal_node;
    point_type goal_pos;
    double goal_distance;
    double max_edge_length;
    double robot_radius;
    unsigned int max_vertex_count;
    unsigned int nn_search_divider;
    unsigned int num_neighbors;
    double max_neighbor_radius;
    unsigned int num_progress_reports;

    boost::minstd_rand m_rng;
    boost::rectangle_topology<boost::minstd_rand> m_space;
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position;
    boost::property_map<WorldGridType, boost::vertex_distance_t>::type m_distance;
    boost::property_map<WorldGridType, boost::vertex_predecessor_t>::type m_pred;
    boost::property_map<WorldGridType, boost::vertex_density_t>::type m_density;
    boost::property_map<WorldGridType, boost::edge_weight_t>::type m_weight;
    
    
    typedef ReaK::pp::dvp_tree<VertexType, 
                               boost::rectangle_topology<>, 
			       boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
			       4> WorldPartition;
    WorldPartition space_part;
    
    ProgressCallback progress_call_back;
    PathFoundCallback path_found_call_back;

    bool is_free(const point_type& p) const {
      if((p[0] < 0) || (p[0] >= grid_width) || (p[1] < 0) || (p[1] >= grid_height)) 
	return false;
      const uchar* color_bits = world_map_image.ptr(int(p[1]));
      color_bits += bpp * int(p[0]);
      if( (color_bits[0] == color_bits[1]) &&
	  (color_bits[0] == color_bits[2]) &&
	  (color_bits[0] < 250) ) {
	return false;
      } else {
        return true;
      };
    };

    void draw_edge(const point_type& p_u, const point_type& p_v, bool goal_path = false) {
      point_difference_type diff = m_space.difference(p_v, p_u);
      double dist = m_space.norm(diff);
      double d = 0.0;
      while(d <= dist) {
	point_type p = m_space.adjust(p_u, diff * (d / dist));
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
    
    //Topology concepts:
    point_type random_point() const {
      point_type result;
      while(!is_free(result = m_space.random_point())) ; //output only free C-space points.
      return result;
    };
    
    double distance(const point_type& p1, const point_type& p2) const {
      if(m_space.distance(p2,move_position_toward(p1,1.0,p2)) < std::numeric_limits< double >::epsilon())
	return m_space.distance(p1,p2); //if p2 is reachable from p1, use Euclidean distance.
      else
        return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1.
    };
    
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      point_difference_type diff = m_space.difference(p2,p1);
      double dist = m_space.norm(diff);
      double d = 1.0;
      while(d < dist * fraction) {
	if (!is_free(m_space.adjust(p1, diff * (d / dist)))) {
	  return m_space.adjust(p1,diff * ((d - 1.0) / dist));
	};
        d += 1.0;
      };
      if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
	return p2;
      else if(fraction == 0.0)
	return p1;
      else 
	return m_space.adjust(p1,diff * fraction);
    };
    
    void get_best_solution(std::list< point_type >& path) {
      using namespace boost;
      path.clear();
      VertexType v = goal_node;
      VertexType u = get(m_pred, v);
      point_type p_u = get(m_position, u);
      std::set<VertexType> path_points;
      path_points.insert(v);
      path.push_front(goal_pos);
      while((u != start_node) && (path_points.insert(u).second)) {
        path.push_front(p_u);
	v = u; 
	u = get(m_pred, v);
	p_u = get(m_position, u);
      };
      if(u == start_node)
	path.push_front(p_u);
    };

    //PRM Visitor concepts:
    void select_neighborhood(const point_type& p, std::vector<VertexType>& Nc, WorldGridType& g, const prm_test_world&, boost::property_map<WorldGridType, boost::vertex_position_t>::type) {
      if(nn_search_divider == 0) {
	ReaK::pp::linear_neighbor_search<>()(p,Nc,g,m_space,m_position,num_neighbors,max_neighbor_radius);
      } else if(nn_search_divider == 1) {
	std::multimap<double,VertexType> neighbors;
	space_part.find_nearest(p,neighbors,num_neighbors);
	std::multimap<double,VertexType>::iterator it_end = neighbors.lower_bound(max_neighbor_radius);
	for(std::multimap<double,VertexType>::iterator it = neighbors.begin(); it != it_end; ++it)
	  Nc.push_back(it->second);
      } else {
	ReaK::pp::best_only_neighbor_search<>(nn_search_divider).operator()(p,Nc,g,m_space,m_position,num_neighbors,max_neighbor_radius);
      };
    };
    
    void vertex_added(VertexType u, WorldGridType& g) {
      if(nn_search_divider == 1) {
	space_part.insert(u);
      };
    };
    
    void update_density(VertexType u, WorldGridType& g) {
      //take the sum of all weights of outgoing edges.
      if(boost::out_degree(u,g) == 0) {
	boost::put(m_density, u, 0.0);
	return;
      };
      double sum = 0.0;
      OutEdgeIter ei, ei_end;
      for(boost::tie(ei,ei_end) = boost::out_edges(u,g); ei != ei_end; ++ei) {
	sum += boost::get(m_weight, *ei);
      };
      sum /= boost::out_degree(u,g) * boost::out_degree(u,g); //this gives the average edge distances divided by the number of adjacent nodes.
      boost::put(m_density, u, std::exp(-sum*sum));
    };
    
    VertexType expand_vertex(VertexType u, WorldGridType& g) {
      point_type p_u = boost::get(m_position, u);
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = m_space.random_point();
        double dist = m_space.distance(p_u,p_rnd);
        p_v = move_position_toward(p_u, max_edge_length / dist, p_rnd);
	++i;
      } while((m_space.distance(p_u, p_v) < 1.0) && (i <= 10));
      if(i > 10) {
	//could not expand vertex u, then just generate a random C-free point.
	p_v = random_point();
      };
      
      VertexType v = boost::add_vertex(g);
      boost::put(m_position, v, p_v);
      return v;
    };

    void edge_added(EdgeType e, WorldGridType& g) {
      using namespace boost;
      VertexType u = source(e,g);
      VertexType v = target(e,g);
      point_type p_u = get(m_position, u);
      point_type p_v = get(m_position, v);
      
      boost::put(m_weight, e, m_space.distance(p_u,p_v));
      
      draw_edge(p_u, p_v);

      if(( num_vertices(grid) > (num_progress_reports + 1) * 100) && (progress_call_back)) {
	++num_progress_reports;
	progress_call_back(world_map_output, num_vertices(grid));
      };

#if 0
      //now, check if v is connected with the goal on a straight collision-free path.
      if(&g == &grid) {
      pixel_difference diff = m_space.difference(goal_pos, p_v);
      double dist = m_space.norm(diff);
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
#endif
    };
    
    void register_solution() {
      using namespace boost;
      
      goal_distance = get(m_distance, goal_node);
      
      if(goal_distance < std::numeric_limits<double>::infinity()) {
        //Draw the edges of the current best solution:
        VertexType v = goal_node;
        point_type p_v = goal_pos;
        VertexType u = get(m_pred, v);
        point_type p_u = get(m_position, u);
        std::set<VertexType> path;
        path.insert(v);
        while((u != start_node) && (path.insert(u).second)) {
  	  draw_edge(p_u,p_v,true);
	  v = u; p_v = p_u;
	  u = get(m_pred, v);
	  p_u = get(m_position, u);
        };
	if(u == start_node)
	  draw_edge(p_u,p_v,true);
        path.clear();
      
        if(path_found_call_back) {
	  path_found_call_back(world_map_output, num_vertices(grid), 0, goal_distance);
        };
      };
    };
    
    bool keep_going() {
      return true;
    };
    
    double heuristic(VertexType u) {
      point_type p_u = boost::get(m_position, u);
      return m_space.distance(p_u, goal_pos);
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
    prm_test_world(const cv::Mat& aWorldMapImage, double aMaxEdgeLength, 
		   double aRobotRadius, unsigned int aMaxVertexCount,  
		   ProgressCallback aProgressCallback = ProgressCallback(), 
		   PathFoundCallback aPathFoundCallback = PathFoundCallback(), 
		   unsigned int aNNSearchDivider = 10, unsigned int aNumNeighbors = 6, 
		   double aMaxNeighborRadius = 30.0) :
                   world_map_image(aWorldMapImage.clone()), world_map_output(aWorldMapImage.clone()),
                   grid_width(aWorldMapImage.size().width),
                   grid_height(aWorldMapImage.size().height),
                   goal_distance(std::numeric_limits<double>::infinity()),
                   max_edge_length(aMaxEdgeLength), 
                   robot_radius(aRobotRadius), max_vertex_count(aMaxVertexCount),
                   nn_search_divider(aNNSearchDivider), num_neighbors(aNumNeighbors),
                   max_neighbor_radius(aMaxNeighborRadius), num_progress_reports(0), m_rng(boost::minstd_rand(std::time(0))),
                   m_space(m_rng,0, aWorldMapImage.size().height, aWorldMapImage.size().width, 0),
                   m_position(boost::get(boost::vertex_position, grid)),
                   m_distance(boost::get(boost::vertex_distance, grid)),
                   m_pred(boost::get(boost::vertex_predecessor, grid)),
                   m_density(boost::get(boost::vertex_density, grid)),
                   m_weight(boost::get(boost::edge_weight, grid)),
                   space_part(grid,m_space,m_position),
                   progress_call_back(aProgressCallback), path_found_call_back(aPathFoundCallback)
    {
      using namespace boost;
      
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
	    point_type p; p[0] = x; p[1] = y;
	    put(m_position, current_pos, p);
	    put(m_pred, current_pos, current_pos);
	    put(m_distance, current_pos, 0.0);
	    if(nn_search_divider == 1)
              space_part.insert(start_node);
	  } else if( (color_bits[2] == 0) &&
	             (color_bits[0] == 0) &&
	             (color_bits[1] == 255) ) {
	    //this is the goal position.
	    goal_pos[0] = x;
	    goal_pos[1] = y;
	    goal_node = add_vertex(grid);
	    put(m_position, goal_node, goal_pos);
	    put(m_pred, goal_node, goal_node);
	    put(m_distance, goal_node, std::numeric_limits< double >::infinity());
	    if(nn_search_divider == 1)
              space_part.insert(goal_node);
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
    ~prm_test_world() { };

    void set_start_pos(const point_type& aStart) {
      while(boost::num_vertices(grid)) {
	space_part.erase(*(boost::vertices(grid).first));
	boost::remove_vertex(*(boost::vertices(grid).first),grid);
      };
      start_node = current_pos = boost::add_vertex(grid);
      put(m_position, current_pos, aStart);
      put(m_distance, current_pos, 0.0);
      if(nn_search_divider == 1)
        space_part.insert(start_node);
    };
    
    void set_goal_pos(const point_type& aGoal) {
      goal_pos = aGoal;
      goal_node = boost::add_vertex(grid);
      put(m_position, goal_node, aGoal);
      put(m_distance, goal_node, std::numeric_limits< double >::infinity());
      if(nn_search_divider == 1)
        space_part.insert(goal_node);
    };
    
    double run() {
      using namespace boost;
      unsigned int m = max_vertex_count; 
      num_progress_reports = 0;
      while(!(goal_distance < std::numeric_limits<double>::infinity())) {
	std::cout << "Generating PRM.." << std::endl;
	ReaK::graph::generate_prm(
	      grid, 
	      *this,
	      ReaK::graph::make_composite_prm_visitor(
	        bind(&prm_test_world::vertex_added,this,_1,_2),
                bind(&prm_test_world::edge_added,this,_1,_2),
	        bind(&prm_test_world::expand_vertex,this,_1,_2),
                bind(&prm_test_world::update_density,this,_1,_2)),
	      m_position,
	      m_density,
 	      bind(&prm_test_world::select_neighborhood,this,_1,_2,_3,_4,_5),
              m, max_vertex_count / 10, max_vertex_count / 50,
              bind(&prm_test_world::keep_going,this),
	      std::less<double>()); 
	std::cout << "Solving shortest-distance.." << std::endl;
	astar_search(grid,
	             start_node,
		     bind(&prm_test_world::heuristic,this,_1),
		     default_astar_visitor(),
		     m_pred,
		     get(vertex_rhs, grid),
		     m_distance,
		     m_weight,
		     double(0.0),
		     get(vertex_color,grid),
		     std::less<double>(), std::plus<double>(),
		     std::numeric_limits< double >::infinity(),
		     double(0.0)); 
	std::cout << "Printing results..." << std::endl;
	register_solution();
	m += max_vertex_count;
      };
      return goal_distance;
    };
};

#endif

