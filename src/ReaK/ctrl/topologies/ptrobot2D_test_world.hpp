/**
 * \file ptrobot2D_test_world.hpp
 * 
 * This library defines a class for path-planning problems on a world map represented by an image of 
 * occupied versus free pixels for a point-robot also the size of one pixel. 
 * It takes a world map (as an OpenCV image) where any non-white gray-scaled pixel 
 * is considered occupied (not C-free). It parses the image for a blue pixel and a green 
 * pixel which each represent the start and goal positions. Alternatively, the start 
 * and goal position can be set via set_start_pos and set_goal_pos functions. The class 
 * also allows for many parameters and callbacks, see the constructor's documentation for details.
 * See the test_prm.cpp file for a program that uses this class.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PTROBOT2D_TEST_WORLD_HPP
#define REAK_PTROBOT2D_TEST_WORLD_HPP


#include <opencv/cv.h>
#include <opencv/highgui.h>

#include "path_planning/random_sampler_concept.hpp"
#include "path_planning/metric_space_concept.hpp"

#include "basic_distance_metrics.hpp"
#include "default_random_sampler.hpp"

#include "hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"

namespace ReaK {

namespace pp {


/**
 * This class is used to represent the free-space consisting of the white pixels of an image. The configuration
 * space is thus the space of all pixels in the image (rectangular), and this class restricts the points to only
 * those which correspond to a white pixel in the given image. This class will also scan the given image to find 
 * a blue and green (pure blue, pure green) pixel which each represent the starting location and goal location in 
 * a path-planning problem, respectively.
 */
class ptrobot2D_test_world {
  public:
    typedef hyperbox_topology< ReaK::vect<int,2> > super_space_type;
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    cv::Mat world_map_image;
    cv::Mat world_map_output;
    int grid_width;
    int grid_height;
    int bpp;
    
    point_type start_pos;
    point_type goal_pos;
    double max_edge_length;
    
    super_space_type m_space;
    typename metric_space_traits<super_space_type>::distance_metric_type m_distance;
    typename point_distribution_traits<super_space_type>::random_sampler_type m_rand_sampler;
    

  public:
    
    super_space_type& get_super_space() { return m_space; };
    const super_space_type& get_super_space() const { return m_space; };
    
    
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
      double dist = m_distance(p_v, p_u, m_space);
      double d = 0.0;
      while(d <= dist) {
	point_type p = m_space.move_position_toward(p_u, (d / dist), p_v);
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
    
    //Topology concepts:
    point_type random_point() const {
      point_type result;
      while(!is_free(result = m_rand_sampler(m_space))) ; //output only free C-space points.
      return result;
    };
    
    double distance(const point_type& p1, const point_type& p2) const {
      if(m_distance(p2,move_position_toward(p1,1.0,p2), m_space) < std::numeric_limits< double >::epsilon())
	return m_distance(p1, p2, m_space); //if p2 is reachable from p1, use Euclidean distance.
      else
        return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1.
    };
    
    double norm(const point_difference_type& dp) const {
      return m_distance(dp, m_space);
    };
    
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1,p2);
    };
    
    point_type origin() const {
      return m_space.origin();
    };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return move_position_toward(p, 1.0, m_space.adjust(p, dp));
    };
    
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      point_difference_type diff = m_space.difference(p2,p1);
      double dist = m_distance(p1, p2, m_space);
      double d = 1.0;
      while(d < dist * fraction) {
	if (!is_free(m_space.move_position_toward(p1, (d / dist), p2))) {
	  return m_space.move_position_toward(p1,((d - 1.0) / dist), p2);
	};
        d += 1.0;
      };
      if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
	return p2;
      else if(fraction == 0.0)
	return p1;
      else 
	return m_space.move_position_toward(p1, fraction, p2);
    };
    
    point_type random_walk(const point_type& p_u) {
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = m_rand_sampler(m_space);
        double dist = m_distance(p_u, p_rnd, m_space);
        p_v = move_position_toward(p_u, max_edge_length / dist, p_rnd);
	++i;
      } while((m_distance(p_u, p_v, m_space) < 1.0) && (i <= 10));
      if(i > 10) {
	//could not expand vertex u, then just generate a random C-free point.
	p_v = random_point();
      };
      return p_v;
    };
    
    double bird_fly_to_goal(const point_type& p_u) {
      return m_distance(p_u, goal_pos, m_space);
    };
    
    double bird_fly_to_start(const point_type& p_u) {
      return m_distance(p_u, start_pos, m_space);
    };
    
    const point_type& get_start_pos() const {
      start_node = aStart;
    };
    
    const point_type& get_goal_pos() const {
      goal_pos = aGoal;
    };
    
    void set_start_pos(const point_type& aStart) {
      start_node = aStart;
    };
    
    void set_goal_pos(const point_type& aGoal) {
      goal_pos = aGoal;
    };
    
    

    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aWorldMapImage An image which represents the C-free as white (or colored) pixels and the occupied C-space as gray pixels.
     * \param aMaxEdgeLength The maximum length of an added edge, in pixel-units.
     * \param aRobotRadius The radius of the robot (collision radius), in pixel-units.
     */
    ptrobot2D_test_world(const cv::Mat& aWorldMapImage, 
			 double aMaxEdgeLength, 
			 double aRobotRadius) :
			 world_map_image(aWorldMapImage.clone()), 
			 world_map_output(aWorldMapImage.clone()),
			 grid_width(aWorldMapImage.size().width),
			 grid_height(aWorldMapImage.size().height),
			 max_edge_length(aMaxEdgeLength), 
			 m_space("ptrobot2D_space", point_type(0,0), point_type(aWorldMapImage.size().width,aWorldMapImage.size().height)),
			 m_distance(get(distance_metric, m_space)),
			 m_rand_sampler(get(random_sampler, m_space))
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
	    start_pos[0] = x; 
	    start_pos[1] = y;
	  } else if( (color_bits[2] == 0) &&
	             (color_bits[0] == 0) &&
	             (color_bits[1] == 255) ) {
	    //this is the goal position.
	    goal_pos[0] = x;
	    goal_pos[1] = y;
	  };
	  color_bits += bpp;
	};
      };
      
      if(int(std::fabs(aRobotRadius)) > 0) {
        cv::GaussianBlur(world_map_output,world_map_image,
		         cv::Size(int(std::fabs(aRobotRadius))*2+1,int(std::fabs(aRobotRadius))*2+1),
		         aRobotRadius
		        );
      };
    };

    //This is a RAII class, so the destructor is meaningless. Hurray for RAII!!
    ~ptrobot2D_test_world() { };

};


};

};

#endif

