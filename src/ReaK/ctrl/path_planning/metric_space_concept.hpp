
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

#ifndef METRIC_SPACE_CONCEPT_HPP
#define METRIC_SPACE_CONCEPT_HPP


#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {
  
  

template <typename Topology>
struct metric_topology_traits {
  typedef Topology::point_type point_type;
  typedef Topology::point_difference_type point_difference_type;
    
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = point_type::dimensions);
  
};


template <typename DistanceMetric, typename Topology>
struct DistanceMetricConcept {
  DistanceMetric d;
  Topology s;
  Topology::point_type p1, p2;
  Topology::point_difference_type pd;
  double dist;
  void constraints() {
    dist = d(p1, p2, s);
    dist = d(pd, s);
  };
  
};


template <typename Topology>
struct MetricSpaceConcept {
  typename metric_topology_traits<Topology>::point_type p1, p2;
  typename metric_topology_traits<Topology>::point_difference_type pd;
  Topology space;
  double d;
  void constraints() {
    d  = space.distance(p1, p2);
    d  = space.norm(pd);
    p1 = space.random_point();
    pd = space.difference(p1,p2);
    p1 = space.move_position_toward(p1,d,p2);
    p1 = space.origin();
    p1 = space.adjust(p1,d * pd);
    pd = -pd;
  };
  
};


struct default_distance_metric {
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    return s.distance(a, b);
  };
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    return s.norm(a);
  };
};


};

};


#endif


