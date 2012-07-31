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


#include "path_planning/random_sampler_concept.hpp"
#include "path_planning/metric_space_concept.hpp"

#include "basic_distance_metrics.hpp"
#include "default_random_sampler.hpp"

#include "hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"

namespace ReaK {

namespace pp {


class ptrobot2D_test_world_impl; // forward-declaration, a faint cat.

/**
 * This class is used to represent the free-space consisting of the white pixels of an image. The configuration
 * space is thus the space of all pixels in the image (rectangular), and this class restricts the points to only
 * those which correspond to a white pixel in the given image. This class will also scan the given image to find 
 * a blue and green (pure blue, pure green) pixel which each represent the starting location and goal location in 
 * a path-planning problem, respectively.
 */
class ptrobot2D_test_world : public named_object {
  public:
    typedef hyperbox_topology< ReaK::vect<double,2> > super_space_type;
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    ptrobot2D_test_world_impl* pimpl;
    std::string world_map_file_name;
    double robot_radius;
    double max_edge_length;
    
    super_space_type m_space;
    typename metric_space_traits<super_space_type>::distance_metric_type m_distance;
    typename point_distribution_traits<super_space_type>::random_sampler_type m_rand_sampler;
    
    ptrobot2D_test_world();

  public:
    
    /**
     * Returns a reference to the super-space in which this test-world is embedded.
     * \return A reference to the super-space in which this test-world is embedded.
     */
    super_space_type& get_super_space() { return m_space; };
    
    /**
     * Returns a const-reference to the super-space in which this test-world is embedded.
     * \return A const-reference to the super-space in which this test-world is embedded.
     */
    const super_space_type& get_super_space() const { return m_space; };
    
    
    /**
     * Checks if the given point is within the free-space.
     * \param p The point to be checked for being collision-free.
     * \return True if p is collision-free.
     */
    bool is_free(const point_type& p) const;
    
    /**
     * Resets the output image used to draw the edges of the motion graph.
     */
    void reset_output() const;
    
    /**
     * Saves the output image to a given filename.
     */
    void save_output(const std::string& aFilename) const;
    
    /**
     * Draws the given edge to the output image.
     * \param p_u The start point of the edge.
     * \param p_v The end point of the edge.
     * \param goal_path True if the edge is part of the solution path.
     */
    void draw_edge(const point_type& p_u, const point_type& p_v, bool goal_path = false) const;
    
    //Topology concepts:
    
    /**
     * Produces a random, collision-free point.
     * \return A random, collision-free point.
     */
    point_type random_point() const;
    
    /**
     * Computes the distance between two points. If there is no collision-free line between
     * the two points, the distance is infinite.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return The collision-free distance between the two given points.
     */
    double distance(const point_type& p1, const point_type& p2) const;
    
    /**
     * Computes the norm of the difference between two points. 
     * \param dp The point difference.
     * \return The norm of the difference between the two points.
     */
    double norm(const point_difference_type& dp) const;
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const;
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type origin() const;
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& p, const point_difference_type& dp) const;
    
    /**
     * Returns a point which is at a fraction between two points a to b, or as 
     * far as it can get before a collision.
     */
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const;
    
    /**
     * Returns a random point fairly near to the given point.
     */
    std::pair<point_type, bool> random_walk(const point_type& p_u) const;
    
    double bird_fly_to_goal(const point_type& p_u) const;
    
    double bird_fly_to_start(const point_type& p_u) const;
    
    const point_type& get_start_pos() const;
    
    const point_type& get_goal_pos() const;
    
    void set_start_pos(const point_type& aStart);
    
    void set_goal_pos(const point_type& aGoal);
    
    
    
    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aWorldMapImage The filename of the image which represents the C-free as white (or colored) pixels and the occupied C-space as gray pixels.
     * \param aMaxEdgeLength The maximum length of an added edge, in pixel-units.
     * \param aRobotRadius The radius of the robot (collision radius), in pixel-units.
     */
    ptrobot2D_test_world(const std::string& aWorldMapImage, 
			 double aMaxEdgeLength, 
			 double aRobotRadius);
    
    ptrobot2D_test_world(const ptrobot2D_test_world& rhs);
    
    ptrobot2D_test_world& operator=(const ptrobot2D_test_world& rhs);
    
    ~ptrobot2D_test_world();
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(ptrobot2D_test_world,0xC2400020,1,"ptrobot2D_test_world",named_object)
    
    
};


};

};

#endif

