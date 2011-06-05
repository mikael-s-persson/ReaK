
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

#ifndef LINE_TOPOLOGY_HPP
#define LINE_TOPOLOGY_HPP


#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT
#include <boost/shared_ptr.hpp>

#include <cmath>

namespace ReaK {

namespace pp {

class line_topology 
{
  public:
    typedef double point_type;
    typedef double point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 1);
    
    double distance(const point_type& a, const point_type& b) const 
    {
      return std::fabs(b - a);
    }

    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return a + (b - a) * fraction;
    }

    point_difference_type difference(const point_type& a, const point_type& b) const {
      return a - b;
    }

    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return a + delta;
    }

    point_type pointwise_min(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MIN();
      return min BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    }

    point_type pointwise_max(const point_type& a, const point_type& b) const {
      BOOST_USING_STD_MAX();
      return max BOOST_PREVENT_MACRO_SUBSTITUTION (a, b);
    }

    double norm(const point_difference_type& delta) const {
      return std::fabs(delta);
    }

    double volume(const point_difference_type& delta) const {
      return delta;
    }

};

template<typename RandomNumberGenerator = boost::minstd_rand>
class line_segment_topology : public line_topology
{
  typedef boost::uniform_01<RandomNumberGenerator, double> rand_t;

  public:
    typedef line_topology::point_type point_type;
    typedef line_topology::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = line_topology::dimensions);

    explicit line_segment_topology(double aScaling = 1.0, double aOrigin = 0.0) 
      : gen_ptr(new RandomNumberGenerator), rand(new rand_t(*gen_ptr)), 
        scaling(scaling), origin(aOrigin) { };

    line_segment_topology(RandomNumberGenerator& aGen, double aScaling = 1.0, double aOrigin = 0.0) 
      : gen_ptr(), rand(new rand_t(aGen)), scaling(aScaling), origin(aOrigin) { };
                     
    point_type random_point() const {
      return (*rand)() * scaling + origin;
    };

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

    double distance_from_boundary(point_type a) const {
      double dist = std::fabs(scaling - a + origin);
      if(std::fabs(a - origin) < dist)
        return std::fabs(a - origin);
      else
        return dist;
    };

    point_type center() const {
      return origin + scaling * 0.5;
    };

    point_type origin() const {
      return origin;
    };

    point_difference_type extent() const {
      return origin + scaling;
    };

   private:
    boost::shared_ptr<RandomNumberGenerator> gen_ptr;
    boost::shared_ptr<rand_t> rand;
    point_difference_type scaling;
    point_type origin;
};



};

};

#endif








