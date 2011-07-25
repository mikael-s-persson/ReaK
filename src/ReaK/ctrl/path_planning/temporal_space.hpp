
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

#ifndef TEMPORAL_SPACE_HPP
#define TEMPORAL_SPACE_HPP

#include <boost/config.hpp>
#include <cmath>
#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/shared_ptr.hpp>

#include "line_topology.hpp"
#include "temporal_space_concept.hpp"

namespace ReaK {

namespace pp {


template <typename Topology, typename DistanceMetric = spatial_distance_only, typename RandomNumberGenerator = boost::minstd_rand>
class temporal_space {
  public:
    typedef line_segment_topology<RandomNumberGenerator> time_topology;
    typedef Topology space_topology;
    typedef DistanceMetric distance_metric;
    
    struct point {
      time_topology::point_type time;
      space_topology::point_type pt;
    
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = space_topology::point_type::dimensions);
      point() { };
      point(const time_topology::point_type& aTime, const space_topology::point_type& aPt) : time(aTime), pt(aPt) { };
      double& operator[](std::size_t i) { return pt[i]; };
      const double& operator[](std::size_t i) const { return pt[i]; };
    };
  
    struct point_difference {
      time_topology::point_difference_type time;
      space_topology::point_difference_type pt;
      
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = space_topology::point_difference_type::dimensions);
      point_difference() : time(0.0) { };
      double& operator[](std::size_t i) { return pt[i]; };
      const double& operator[](std::size_t i) const { return pt[i]; };

      friend point_difference operator+(const point_difference& a, const point_difference& b) {
        point_difference result;
	result.time = a.time + b.time;
	result.pt = a.pt + b.pt;
        return result;
      };

      point_difference& operator+=(const point_difference& b) {
        time += b.time;
	pt += b.pt;
        return *this;
      };

      friend point_difference operator-(const point_difference& a) {
        point_difference result;
        result.time = -a.time;
	result.pt = -a.pt;
        return result;
      };

      friend point_difference operator-(const point_difference& a, const point_difference& b) {
        point_difference result;
	result.time = a.time - b.time;
	result.pt = a.pt - b.pt;
        return result;
      };

      point_difference& operator-=(const point_difference& b) {
        time -= b.time;
	pt -= b.pt;
        return *this;
      };

      friend point_difference operator*(const point_difference& a, double b) {
        point_difference result;
	result.time = a.time * b;
	result.pt = a.pt * b;
        return result;
      };

      friend point_difference operator*(double a, const point_difference& b) {
        point_difference result;
	result.time = a * b.time;
	result.pt = a * b.pt;
        return result;
      };

      friend double dot(const point_difference& a, const point_difference& b) {
        return dot(a.pt, b.pt);
      };

    };
  protected:
    const space_topology& space;
    time_topology time;
    distance_metric dist;
        
    temporal_space(const temporal_space&); //non-copyable.
    temporal_space& operator=(const temporal_space&);
    
  public:
    explicit temporal_space(const space_topology& aSpace, double aMaxTime = 1.0) :
                            space(aSpace), time(aMaxTime,0.0) { };
   
    temporal_space(const space_topology& aSpace, RandomNumberGenerator& aGen, double aMaxTime = 1.0) :
                   space(aSpace), time(aGen,aMaxTime,0.0) { };
    
    typedef point point_type;
    typedef point_difference point_difference_type;

    point random_point() const 
    {
      point p;
      p.time = time.random_point();
      p.pt = space.random_point();
      return p;
    }
    
    double distance(const point& a, const point& b) const {
      return dist(a, b, time, space);
    };

    point move_position_toward(const point& a, double fraction, const point& b) const {
      point result;
      result.time = time.move_position_toward(a.time, fraction, b.time);
      result.pt = space.move_position_toward(a.pt, fraction, b.pt);
      return result;
    };

    point_difference difference(const point& a, const point& b) const {
      point_difference result;
      result.time = time.difference(a.time, b.time);
      result.pt = space.difference(a.pt, b.pt);
      return result;
    };

    point adjust(const point& a, const point_difference& delta) const {
      point result;
      result.time = time.adjust(a.time, delta.time);
      result.pt = space.adjust(a.pt, delta.pt);
      return result;
    };
  
    point origin() const {
      return point(time.origin(), space.origin());
    };

    double norm(const point_difference& a) const {
      return d(a, time, space);
    };
    
};



};

};

#endif
















