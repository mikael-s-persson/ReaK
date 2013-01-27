/**
 * \file fadprm_test_world.hpp
 * 
 * This library defines an example class that uses the FADPRM algorithm on a world map. This class 
 * basically just wraps the algorithms found in "fadprm.hpp" for a simple 2D path-planning problem.
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free), and where the gray-value determines the time at which 
 * the obstacle is discovered. It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 * See the test_fadprm.cpp file for a program that uses this class.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_FADPRM_TEST_WORLD_HPP
#define REAK_FADPRM_TEST_WORLD_HPP


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>
#include <boost/bind.hpp>


#include "fadprm.hpp"
#include "path_planning/topological_search.hpp"
#include <cmath>

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/default_random_sampler.hpp"


#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>


namespace boost {

enum vertex_density_t { vertex_density };

BOOST_INSTALL_PROPERTY(vertex,density);

};


namespace ReaK {

namespace graph {
  
/**
 * This class uses the FADPRM algorithm on a world map. This class 
 * basically just wraps the algorithms found in "fadprm.hpp" for a simple 2D path-planning problem.
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free), and where the gray-value determines the time at which 
 * the obstacle is discovered. It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 * See the test_fadprm.cpp file for a program that uses this class.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
 */
class fadprm_test_world {
  public:
    
    typedef boost::rectangle_topology<>::point_type point_type;
    typedef boost::rectangle_topology<>::point_difference_type point_difference_type;
    
    typedef ReaK::pp::default_distance_metric distance_metric_type;
    typedef ReaK::pp::default_random_sampler random_sampler_type;

    typedef boost::property< boost::vertex_position_t, point_type,
            boost::property< boost::vertex_density_t, double,
            boost::property< boost::vertex_heuristic_t, double,
	    boost::property< boost::vertex_rhs_t, double,
	    boost::property< boost::vertex_key_t, ReaK::graph::adstar_key_value<double>,
	    boost::property< boost::vertex_distance_t, double,
	    boost::property< boost::vertex_color_t, boost::default_color_type,
	    boost::property< boost::vertex_predecessor_t, boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::undirectedS,boost::vecS>::vertex_descriptor, boost::no_property > > > > > > > > WorldGridVertexProperties;

    typedef boost::property< boost::edge_weight_t, double, boost::no_property > WorldGridEdgeProperties;

    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                   WorldGridVertexProperties,
				   WorldGridEdgeProperties,
				   boost::vecS > WorldGridType;

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
    VertexType goal_pos;
    int current_time;
    double initial_beta;
    double goal_distance;
    double max_edge_length;
    double robot_radius;
    unsigned int max_vertex_count;
    unsigned int current_max_vertex_count;
    unsigned int nn_search_divider;
    unsigned int num_neighbors;
    double max_neighbor_radius;
    unsigned int max_node_degree;
    unsigned int num_progress_reports;
    
    boost::minstd_rand m_rng;
    boost::rectangle_topology<boost::minstd_rand> m_space;
    boost::property_map<WorldGridType, boost::vertex_predecessor_t>::type m_pred;
    boost::property_map<WorldGridType, boost::vertex_position_t>::type    m_position;
    boost::property_map<WorldGridType, boost::vertex_density_t>::type     m_density;
    boost::property_map<WorldGridType, boost::vertex_heuristic_t>::type   m_heuristic;
    boost::property_map<WorldGridType, boost::vertex_distance_t>::type    m_distance;
    boost::property_map<WorldGridType, boost::vertex_color_t>::type       m_color;
    boost::property_map<WorldGridType, boost::edge_weight_t>::type        m_weight;
    
    ProgressCallback progress_call_back;
    PathFoundCallback path_found_call_back;

    bool is_free(const point_type& p) const {
      if((p[0] < 0) || (p[0] >= grid_width) || (p[1] < 0) || (p[1] >= grid_height)) 
	return false;
      const uchar* color_bits = world_map_image.ptr(int(p[1]));
      color_bits += bpp * int(p[0]);
      if( (color_bits[0] == color_bits[1]) &&
	  (color_bits[0] == color_bits[2]) &&
	  (color_bits[0] < std::min(255,current_time)) ) {
	return false;
      } else {
        return true;
      };
    };
    
    void draw_edge(const point_type& p_u, const point_type& p_v, bool goal_path = false) {
      point_difference_type diff = m_space.difference(p_v, p_u); 
      double dist = m_space.norm(diff); 
      double d = 0.0; 
      if(dist < std::numeric_limits<double>::epsilon()) return;
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
    
    double updateEdgeWeight(EdgeType e) {
      point_type source_pos = boost::get(m_position, source(e, grid));
      point_type target_pos = boost::get(m_position, target(e, grid));
      if(m_space.distance(target_pos,move_position_toward(source_pos,1.0,target_pos)) < std::numeric_limits< double >::epsilon()) {
	boost::put(m_weight, e, m_space.distance(source_pos,target_pos));
	return 0.0; //if edge is connected in C_free, then no change in weight.
      } else {
	double old_w = boost::get(m_weight, e);
	boost::put(m_weight, e, 2550.0);
        return 2550.0 - old_w; //edge is not good, infinite weight.
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
    
    double norm(const point_difference_type& dp) const {
      return m_space.norm(dp);
    };
    
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1,p2);
    };
    
    point_type origin() const {
      return m_space.origin();
    };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return m_space.adjust(p, dp);
    };
    
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      point_difference_type diff = m_space.difference(p2,p1);
      double dist = m_space.norm(diff);
      double d = 1.0;
      while(d < dist * fraction) {
	if (!is_free(m_space.adjust(p1, diff * (d / dist)))) {
	  return m_space.adjust(p1,diff * (0.95*(d - 1.0) / dist));
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
      VertexType v = current_pos;
      VertexType u = get(m_pred, v);
      point_type p_u = get(m_position, u);
      std::set<VertexType> path_points;
      path_points.insert(v);
      path.push_front(get(m_position, v));
      while((u != goal_pos) && (path_points.insert(u).second)) {
        path.push_front(p_u);
	v = u; 
	u = get(m_pred, v);
	p_u = get(m_position, u);
      };
      if(u == goal_pos)
	path.push_front(p_u);
    };
    
    
    //PRM Visitor concepts:
    void select_neighborhood(const point_type& p, std::vector<VertexType>& Nc, WorldGridType& g, const fadprm_test_world&, boost::property_map<WorldGridType, boost::vertex_position_t>::type) {
      Nc.resize(num_neighbors);
      std::vector<VertexType>::iterator last;
      if(nn_search_divider == 0) {
	last = ReaK::pp::linear_neighbor_search<>()(p,Nc.begin(),g,*this,m_position,num_neighbors,max_neighbor_radius);
      } else {
	last = ReaK::pp::best_only_neighbor_search<>(nn_search_divider).operator()(p,Nc.begin(),g,*this,m_position,num_neighbors,max_neighbor_radius);
      };
      Nc.erase(last, Nc.end());
    };
    
    void vertex_added(VertexType u, WorldGridType& g) {
      //compute heuristic:
      point_type current_coord = get(m_position, current_pos);
      point_type p_u = get(m_position, u);
      
      put(m_heuristic, u, m_space.distance(current_coord,p_u));
    };
    
    void update_density(VertexType u, WorldGridType& g) {
      //take the sum of all weights of outgoing edges.
      unsigned int deg_u = out_degree(u,g);
      if(deg_u == 0) {
	put(m_density, u, 0.0);
	return;
      };
      if(deg_u >= max_node_degree) {
	//no point trying to expand this vertex, at least, from a density stand-point.
	put(m_density, u, (grid_height+grid_width));
	return;
      };
      double sum = 0.0;
      OutEdgeIter ei, ei_end;
      for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei) {
	sum += get(m_weight, *ei);
      };
      sum /= deg_u * (max_node_degree - deg_u); //this gives the average edge distances divided by the number of adjacent nodes.
      put(m_density, u, (grid_height+grid_width)*std::exp(-sum*sum));
    };
    
    void expand_vertex(VertexType u, WorldGridType& g, std::vector<VertexType>& v_list) {
      if(num_vertices(g) > current_max_vertex_count) return;
      point_type p_u = get(m_position, u);
      while(out_degree(u,g) + 2*v_list.size() < max_node_degree) {
        point_type p_rnd, p_v;
        unsigned int i = 0;
        do {
          p_rnd = m_space.random_point();
          double dist = m_space.distance(p_u,p_rnd);
	  double fraction = max_edge_length * 0.001 * (m_rng() % 1000);
          p_v = move_position_toward(p_u, fraction / dist, p_rnd);
	  ++i;
        } while((m_space.distance(p_u, p_v) < 1.0) && (i <= 10));
        if(i <= 100) {
  	  VertexType v = add_vertex(g);
	  put(m_position, v, p_v);
	  v_list.push_back(v);
        };
      };
      if(m_space.distance(p_u,get(m_position,current_pos)) < max_edge_length)
        v_list.push_back(current_pos);
    };

    void edge_added(EdgeType e, WorldGridType& g) {
      VertexType u = source(e,g); 
      VertexType v = target(e,g); 
      point_type p_u = get(m_position, u); 
      point_type p_v = get(m_position, v); 
      
      put(m_weight, e, m_space.distance(p_u,p_v)); 
      
      draw_edge(p_u, p_v); 

      if(( num_vertices(grid) > (num_progress_reports + 1) * 100) && (progress_call_back)) {
	++num_progress_reports; 
	progress_call_back(world_map_output, num_vertices(grid)); 
      };
    };
    
    void register_solution() {
      std::cout << "Publishing path.." << std::endl;
      goal_distance = get(m_distance, current_pos);
      
      if(goal_distance < std::numeric_limits<double>::infinity()) {
        //Draw the edges of the current best solution:
        VertexType v = current_pos;
        point_type p_v = get(m_position, current_pos);
        VertexType u = get(m_pred, v);
        point_type p_u = get(m_position, u);
        std::set<VertexType> path;
        path.insert(v);
	goal_distance = 0.0;
        while((u != goal_pos) && (path.insert(u).second)) {
  	  draw_edge(p_u,p_v,true);
 	  goal_distance += m_space.distance(p_v,p_u);
	  v = u; p_v = p_u;
	  u = get(m_pred, v);
	  p_u = get(m_position, u);
        };
	if(u == goal_pos) {
	  draw_edge(p_u,p_v,true);
	  goal_distance += m_space.distance(p_v,p_u);
	} else {
	  goal_distance = std::numeric_limits< double >::infinity();
	};
        path.clear();
      
        if(path_found_call_back) {
	  path_found_call_back(world_map_output, num_vertices(grid), 0, goal_distance);
        };
	
	//Draw the edges of the current best solution back to orange:
        v = current_pos;
        p_v = get(m_position, current_pos);
        u = get(m_pred, v);
        p_u = get(m_position, u);
        path.insert(v);
	while((u != goal_pos) && (path.insert(u).second)) {
  	  draw_edge(p_u,p_v,false);
 	  v = u; p_v = p_u;
	  u = get(m_pred, v);
	  p_u = get(m_position, u);
        };
	if(u == goal_pos)
	  draw_edge(p_u,p_v,false);
	path.clear();
	
	current_max_vertex_count += max_vertex_count / 100;
      } else if(num_vertices(grid) >= current_max_vertex_count) {
	current_max_vertex_count += max_vertex_count;
      };
      
      current_time++;
      
    };
    
    
    
    

    VertexType getStartNode() {
      return goal_pos;
    };

    double adjustEpsilon(double aOldEpsilon, double aMaxWeightChange) {
      if(aMaxWeightChange > 5)
	return initial_beta;
      else
        return (1.0 - aOldEpsilon) * 0.5 + aOldEpsilon;
    };

    bool isGoalNotReached() {
      if(current_pos == goal_pos)
	return false;
      else
	return true;
    };

    double checkChanges(std::vector< EdgeType >& aList) {
      double max_change = 0.0;
      boost::graph_traits<WorldGridType>::edge_iterator ei, ei_end;
      for( boost::tie(ei,ei_end) = edges(grid); ei != ei_end; ++ei) {
	double ei_change = updateEdgeWeight(*ei);
	if(std::fabs(ei_change) > 1E-3) {
	  aList.push_back(*ei);
	  if(std::fabs(ei_change) > max_change)
	    max_change = std::fabs(ei_change);
	};
      };
      return max_change;
    };

#if 0
    void updatePath() {

      //now update the colors of the world_image and save it.
      for(int y = 0; y < grid_height; ++y) {
	BYTE* color_bits = FreeImage_GetSccurrent_time++;
      std::cout << "\rCurrently at: " << get(m_position, u) << " time: " << current_time << "                ";
      std::cout.flush();

      point_type current_coord = get(m_position, current_pos);
      boost::graph_traits<WorldGridType>::vertex_iterator ui, ui_end;
      for( tie(ui,ui_end) = vertices(grid); ui != ui_end; ++ui) {
	//compute the heuristic value for each node.
	point_type pos = get(m_position, *ui);
	if( *ui != current_pos ) {
	  put(m_heuristic, *ui, std::sqrt((pos[0] - current_coord[0])*(pos[0] - current_coord[0])
	                                + (pos[1] - current_coord[1])*(pos[1] - current_coord[1])));
	} else {
	  put(m_heuristic, *ui, 0.0);
	};
      };anLine(world_map_output, y);
	BYTE* color_bits_orig = FreeImage_GetScanLine(world_map_image, y);

	for(int x = 0; x < grid_width; ++x) {
	  VertexType current_node = vertex(y * grid_width + x, grid);
	  default_color_type col = get(m_color, current_node);
	  if( (color_bits_orig[FI_RGBA_RED] == color_bits_orig[FI_RGBA_GREEN]) &&
	      (color_bits_orig[FI_RGBA_RED] == color_bits_orig[FI_RGBA_BLUE]) &&
	      (((current_time < 255) && (color_bits_orig[FI_RGBA_RED] <= current_time)) ||
               ((current_time >= 255) && (colo
                     r_bits_orig[FI_RGBA_RED] < 255)))) {
	    color_bits[FI_RGBA_RED] = color_bits_orig[FI_RGBA_RED];
	    color_bits[FI_RGBA_GREEN] = color_bits_orig[FI_RGBA_GREEN];
	    color_bits[FI_RGBA_BLUE] = color_bits_orig[FI_RGBA_BLUE];
	  } else {
	    if( col == white_color ) {
	      color_bits[FI_RGBA_RED] = 255;
	      color_bits[FI_RGBA_GREEN] = 255;
	      color_bits[FI_RGBA_BLUE] = 255;
	    } else if( col == gray_color ) {
	      color_bits[FI_RGBA_RED] = 255;
	      color_bits[FI_RGBA_GREEN] = 255;
	      color_bits[FI_RGBA_BLUE] = 80; //light yellow color
	    } else if( col == black_color ) {
	      color_bits[FI_RGBA_RED] = 255;
	      color_bits[FI_RGBA_GREEN] = 140;
	      color_bits[FI_RGBA_BLUE] = 30; //orange color
	    } else if( col == green_color ) {
	      color_bits[FI_RGBA_RED] = 255;
	      color_bits[FI_RGBA_GREEN] = 255;
	      color_bits[FI_RGBA_BLUE] = 80; //light yellow color
	    } else {
	      color_bits[FI_RGBA_RED] = 255;
	      color_bits[FI_RGBA_GREEN] = 140;
	      color_bits[FI_RGBA_BLUE] = 30; //orange color
	    };
	    if(current_node == goal_pos) {
	      color_bits[FI_RGBA_RED] = 0;
	      color_bits[FI_RGBA_GREEN] = 255;
	      color_bits[FI_RGBA_BLUE] = 0; //green color
	    } else if(current_node == current_pos) {
	      color_bits[FI_RGBA_RED] = 0;
	      color_bits[FI_RGBA_GREEN] = 0;
	      color_bits[FI_RGBA_BLUE] = 255; //blue color
	    };
	  };
	  color_bits += bpp;
	  color_bits_orig += bpp;
	};
      };

      VertexType v = current_pos;
      VertexType u = get(m_pred, v);
      std::set<VertexType> path;
      path.insert(v);
      double total_distance = 0;
      while((u != goal_pos) && (path.insert(u).second)) {
	ReaK::vect<int,3> p = get(m_position, u);
	BYTE* color_bits = FreeImage_GetScanLine(world_map_output, p[1]);
	color_bits += bpp*p[0];
	color_bits[FI_RGBA_RED] = 255;
	color_bits[FI_RGBA_GREEN] = 0;
	color_bits[FI_RGBA_BLUE] = 0; //set this color to indicate the planned path.
	total_distance += get(m_weight, edge(u,v,grid).first);
	v = u;
	u = get(m_pred, v);

      };
      total_distance += get(m_weight, edge(u,v,grid).first);
      path.clear();
      std::stringstream ss;
      ss << "test_adstar_results/" << initial_epsilon << "_" << std::setfill('0') << std::setw(5) << current_time << "_" << total_distance << ".bmp";
      FreeImage_Save(FIF_BMP,world_map_output,ss.str().c_str(),BMP_DEFAULT);

      //now v stores the next position, so lets move there:
      u = get(m_pred, current_pos);
      if(get(m_position,u)[2] == 255) current_pos = u;
      current_time++;
      std::cout << "\rCurrently at: " << get(m_position, u) << " time: " << current_time << "                ";
      std::cout.flush();

      point_type current_coord = get(m_position, current_pos);
      boost::graph_traits<WorldGridType>::vertex_iterator ui, ui_end;
      for( tie(ui,ui_end) = vertices(grid); ui != ui_end; ++ui) {
	//compute the heuristic value for each node.
	point_type pos = get(m_position, *ui);
	if( *ui != current_pos ) {
	  put(m_heuristic, *ui, std::sqrt((pos[0] - current_coord[0])*(pos[0] - current_coord[0])
	                                + (pos[1] - current_coord[1])*(pos[1] - current_coord[1])));
	} else {
	  put(m_heuristic, *ui, 0.0);
	};
      };

      return;
    };
#endif
    
    fadprm_test_world(const cv::Mat& aWorldMapImage, double aInitialBeta, double aMaxEdgeLength, 
		      double aRobotRadius, unsigned int aMaxVertexCount,  
		      ProgressCallback aProgressCallback = ProgressCallback(), 
		      PathFoundCallback aPathFoundCallback = PathFoundCallback(), 
		      unsigned int aNNSearchDivider = 10, unsigned int aNumNeighbors = 6, 
		      double aMaxNeighborRadius = 30.0, unsigned int aMaxNodeDegree = 4) :
                      world_map_image(aWorldMapImage.clone()), world_map_output(aWorldMapImage.clone()),
                      grid_width(aWorldMapImage.size().width),
                      grid_height(aWorldMapImage.size().height),
                      current_time(10), initial_beta(aInitialBeta),
                      goal_distance(std::numeric_limits<double>::infinity()),
                      max_edge_length(aMaxEdgeLength), robot_radius(aRobotRadius), 
                      max_vertex_count(aMaxVertexCount), current_max_vertex_count(aMaxVertexCount), nn_search_divider(aNNSearchDivider), 
                      num_neighbors(aNumNeighbors), max_neighbor_radius(aMaxNeighborRadius), 
                      max_node_degree(aMaxNodeDegree), num_progress_reports(0), m_rng(boost::minstd_rand(std::time(0))),
                      m_space(m_rng,0, aWorldMapImage.size().height, aWorldMapImage.size().width, 0),
                      m_pred(get(boost::vertex_predecessor, grid)),
                      m_position(get(boost::vertex_position, grid)),
                      m_density(get(boost::vertex_density, grid)),
                      m_heuristic(get(boost::vertex_heuristic, grid)),
                      m_distance(get(boost::vertex_distance, grid)),
                      m_color(get(boost::vertex_color, grid)),
                      m_weight(get(boost::edge_weight, grid)),
                      progress_call_back(aProgressCallback), path_found_call_back(aPathFoundCallback)
                      
    {
      using namespace boost;
      
      if(world_map_image.empty()) {
	std::cout << __FILE__ << ":" << __LINE__ << " Error: The world image is empty!" << std::endl;
	throw int(0);
      };

      bpp = world_map_image.elemSize();

      int goal_start_count = 0;
      for(int y = 0; y < grid_height; ++y) {
	uchar* color_bits = world_map_image.ptr(y);

	for(int x = 0; x < grid_width; ++x) {
	  if( (color_bits[2] == 0) &&
	      (color_bits[0] == 255) &&
	      (color_bits[1] == 0) ) {
	    //this is the start position.
	    current_pos = add_vertex(grid);
	    point_type p; p[0] = x; p[1] = y;
	    put(m_position, current_pos, p);
	    goal_start_count++;
	  } else if( (color_bits[2] == 0) &&
	             (color_bits[0] == 0) &&
	             (color_bits[1] == 255) ) {
	    //this is the goal position.
	    point_type p; p[0] = x; p[1] = y;
	    goal_pos = add_vertex(grid);
	    put(m_position, goal_pos, p);
	    goal_start_count++;
	  };
	  color_bits += bpp;
	};
      };
      
      if(goal_start_count != 2) {
	std::cout << __FILE__ << ":" << __LINE__ << " Error: The world image must include one start and one goal pixel, in blue and green respectively!" << std::endl;
	throw int(0);
      };
      
      put(m_heuristic, current_pos, 0.0);
      put(m_heuristic, goal_pos, m_space.distance(get(m_position,goal_pos),get(m_position,current_pos)));
      
    };

    ~fadprm_test_world() {
      
    };
    
    double getHeuristicValue(VertexType u) {
      return get(m_heuristic,u);
    };

    void run() {
      
      generate_fadprm(grid, *this,
		      m_heuristic,
		      make_composite_fadprm_visitor(
			default_adstar_visitor(),
			make_composite_prm_visitor(
			  boost::bind(&fadprm_test_world::vertex_added,this,_1,_2),
		          boost::bind(&fadprm_test_world::edge_added,this,_1,_2),
                          boost::bind(&fadprm_test_world::expand_vertex,this,_1,_2,_3),
                          boost::bind(&fadprm_test_world::update_density,this,_1,_2))),
		      m_pred,
		      m_distance,
		      get(boost::vertex_rhs, grid),
		      get(boost::vertex_key, grid),
		      m_weight,
		      m_density,
		      m_position,
		      boost::bind(&fadprm_test_world::select_neighborhood,this,_1,_2,_3,_4,_5),
		      m_color,
		      boost::bind(&fadprm_test_world::adjustEpsilon,this,_1,_2),
		      boost::bind(&fadprm_test_world::isGoalNotReached, this),
		      boost::bind(&fadprm_test_world::getStartNode,this),
		      boost::bind(&fadprm_test_world::checkChanges,this,_1),
		      boost::bind(&fadprm_test_world::register_solution,this),
		      initial_beta, std::numeric_limits< double >::infinity(),
		      double(0.0), std::less<double>(), std::equal_to<double>(), std::plus<double>(), std::multiplies<double>()
 		     );
    };



};

};

};


#endif












