
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

#ifndef TEMPORAL_SPACES_HPP
#define TEMPORAL_SPACES_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

namespace ReaK {

namespace pp {


template <typename TemporalTopology>
struct temporal_topology_traits {
  typedef TemporalTopology::point_type point_type;
  typedef TemporalTopology::point_difference_type point_difference_type;
  
  typedef TemporalTopology::time_topology time_topology;
  typedef TemporalTopology::space_topology space_topology;
  
  BOOST_STATIC_CONSTANT(std::size_t, time_dimensions = time_topology::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, space_dimensions = space_topology::point_type::dimensions);
  
};
  
  
template <typename DistanceMetric, typename TemporalTopology>
struct TemporalDistMetricConcept {
  DistanceMetric d;
  typename temporal_topology_traits<TemporalTopology>::space_topology s;
  typename temporal_topology_traits<TemporalTopology>::time_topology t;
  typename temporal_topology_traits<TemporalTopology>::point_type p1, p2;
  typename temporal_topology_traits<TemporalTopology>::point_difference_type pd;
  double dist;
  void constraints() {
    dist = d(p1, p2, t, s);
    dist = d(pd, t, s);
  };
  
};


template <typename Topology>
struct TemporalSpaceConcept {
  typename temporal_topology_traits<Topology>::point_type p1, p2;
  typename temporal_topology_traits<Topology>::point_difference_type pd;
  Topology spacetime;
  double d;
  void constraints() {
    boost::function_requires< MetricSpaceConcept< typename temporal_topology_traits<Topology>::space_topology > >();
    boost::function_requires< MetricSpaceConcept< typename temporal_topology_traits<Topology>::time_topology > >();
    d  = spacetime.distance(p1, p2);
    d  = spacetime.norm(pd);
    p1 = spacetime.random_point();
    pd = spacetime.difference(p1,p2);
    p1 = spacetime.move_position_toward(p1,d,p2);
    p1 = spacetime.origin();
    p1 = spacetime.adjust(p1,d * pd);
    pd = -pd;
  };
  
};


struct spatial_distance_only {
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology&, const SpaceTopology& s) const {
    return s.distance(a.pt, b.pt);
  };
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology&, const SpaceTopology& s) const {
    return s.norm(a.pt);
  };
};


struct time_distance_only {
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology& t, const SpaceTopology&) const {
    return t.distance(a.time, b.time);
  };
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology&) const {
    return t.norm(a.time);
  };
};


};

};

#endif


