/**
 * \file line_topology.hpp
 * 
 * This library provides classes that define a line-topology. A line-topology is 
 * a simple metric-space where the points are real values (doubles) along a 1D 
 * space (line-segment).
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

#ifndef REAK_LINE_TOPOLOGY_HPP
#define REAK_LINE_TOPOLOGY_HPP


#include "base/defs.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include <cmath>

namespace ReaK {

namespace pp {

/**
 * This class implements an infinite line topology. Since the space is 
 * infinite, there is no way to generate random points from it, and thus, 
 * this class does not model the topology concepts, but defines a number 
 * of functions useful to a derived class that can provide the full 
 * model of a topology.
 * \tparam T The value-type for the topology (should be an arithmetic type that is implicitly convertable to double).
 */
template <typename T = double>
class line_topology 
{
  public:
    typedef T point_type;
    typedef T point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 1);
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      using std::fabs;
      return fabs(b - a);
    }

    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return a + (b - a) * fraction;
    }

    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      return a - b;
    }

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return a + delta;
    }

    /**
     * Returns the norm of the difference between two points.
     */
    point_type pointwise_min(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MIN();
      return min BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    }

    /**
     * Returns the norm of the difference between two points.
     */
    point_type pointwise_max(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MAX();
      return max BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    }

    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      using std::fabs;
      return fabs(delta);
    }

    /**
     * Returns the volume of the difference between two points.
     */
    double volume(const point_difference_type& delta) const {
      using std::fabs;
      return fabs(delta);
    }

};

/**
 * This class implements a line-segment topology. The space extends from the origin up to some 
 * maximum value.
 * \tparam RandomNumberGenerator A random number generator functor type.
 */
template<typename T = double, typename RandomNumberGenerator = boost::minstd_rand>
class line_segment_topology : public line_topology<T>
{
  typedef boost::uniform_01<RandomNumberGenerator, T> rand_t;

  public:
    typedef line_topology::point_type point_type;
    typedef line_topology::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = line_topology::dimensions);

    /**
     * Default constructor.
     * \param aScaling The overall span of the line-segment.
     * \param aOrigin The minimum bound of the line-segment.
     */
    explicit line_segment_topology(point_type aScaling = point_type(1.0), point_type aOrigin = point_type(0.0)) 
      : gen_ptr(new RandomNumberGenerator), rand(new rand_t(*gen_ptr)), 
        scaling(scaling), origin(aOrigin) { };

    /**
     * Parametrized constructor.
     * \param aGen A random-number generator to use.
     * \param aScaling The overall span of the line-segment.
     * \param aOrigin The minimum bound of the line-segment.
     */
    explicit line_segment_topology(RandomNumberGenerator& aGen, point_type aScaling = point_type(1.0), point_type aOrigin = point_type(0.0)) 
      : gen_ptr(), rand(new rand_t(aGen)), scaling(aScaling), origin(aOrigin) { };
       
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      return (*rand)() * scaling + origin;
    };

    /**
     * Takes a point and clips it to within this line-segment space.
     */
    point_type bound(point_type a) const {
      if(scaling > 0.0) {
        if(a > origin + scaling)
  	  return origin + scaling;
        else if(a < origin) 
 	  return origin;
        else
	  return a;
      } else {
        if(a < origin + scaling)
	  return origin + scaling;
        else if(a > origin)
	  return origin;
        else
	  return a;
      };
    };

    /**
     * Returns the distance to the boundary of the space.
     */
    double distance_from_boundary(point_type a) const {
      using std::fabs;
      double dist = fabs(scaling - a + origin);
      if(fabs(a - origin) < dist)
        return fabs(a - origin);
      else
        return dist;
    };

    /**
     * Returns the center of the space.
     */
    point_type center() const {
      return origin + scaling * 0.5;
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return origin;
    };

    /**
     * Returns the extent of the space (the upper-limit).
     */
    point_difference_type extent() const {
      return origin + scaling;
    };

  private:
    typename shared_pointer<RandomNumberGenerator>::type gen_ptr;
    typename shared_pointer<rand_t>::type rand;
    point_difference_type scaling;
    point_type origin;
};



};

};

#endif








