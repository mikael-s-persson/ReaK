/**
 * \file time_topology.hpp
 * 
 * This library provides classes that define a time-topology. A time-topology is 
 * a simple metric-space where the points are real values (doubles) along a 1D 
 * space. However, because time is unlimited, this topology, although modeling the 
 * MetricSpaceConcept, is not strictly a metric-space since random-points cannot be 
 * generated.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_TIME_TOPOLOGY_HPP
#define REAK_TIME_TOPOLOGY_HPP


#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT
#include <boost/shared_ptr.hpp>

#include <cmath>

namespace ReaK {

namespace pp {

/**
 * This class implements an infinite time-topology. Since the space is 
 * infinite, there is no way to generate random points from it, and thus, 
 * this class does not strictly model the topology concepts, but defines all 
 * the functions required to provide the full model of a MetricSpaceConcept.
 */
class time_topology 
{
  public:
    typedef double point_type;
    typedef double point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 1);
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      return std::fabs(b - a);
    };

    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return a + (b - a) * fraction;
    };

    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      return a - b;
    };

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return a + delta;
    };

    /**
     * Returns the norm of the difference between two points.
     */
    point_type pointwise_min(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MIN();
      return min BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    };

    /**
     * Returns the norm of the difference between two points.
     */
    point_type pointwise_max(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MAX();
      return max BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    };

    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      return std::fabs(delta);
    };

    /**
     * Returns the volume of the difference between two points.
     */
    double volume(const point_difference_type& delta) const {
      return std::fabs(delta);
    };
    
    /**
     * Generates a random point in the space, uniformly distributed.
     * \note This function actually returns the origin of the space.
     */
    point_type random_point() const {
      return 0.0;
    };
    
    /**
     * Returns the origin of the space.
     */
    point_type origin() const {
      return 0.0;
    };

};


};

};

#endif








