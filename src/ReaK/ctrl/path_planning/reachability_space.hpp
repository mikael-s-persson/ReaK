
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

#ifndef REACHABILITY_SPACE_HPP
#define REACHABILITY_SPACE_HPP

#include "temporal_space.hpp"
#include "reachability_space_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {


struct backward_reachable_norm {
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    return a.time - s.norm(a.pt);
  };
  template <typename TimeDiff, typename SpaceDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const TimeDiff& dt, const SpaceDiff& dp, const TimeTopology& t, const SpaceTopology& s) const {
    return dt - s.norm(dp);
  };
};


struct forward_reachable_norm {
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    return a.time + s.norm(a.pt);
  };
  template <typename TimeDiff, typename SpaceDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const TimeDiff& dt, const SpaceDiff& dp, const TimeTopology& t, const SpaceTopology& s) const {
    return dt + s.norm(dp);
  };
};



struct reachable_distance {
  forward_reachable_norm  forward;
  backward_reachable_norm backward;
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology& t, const SpaceTopology& s) const {
    TimeTopology::point_difference_type time_diff = t.difference(b.time, a.time);
    SpaceTopology::point_difference_type point_diff = s.difference(b.pt, a.pt);
    if(backward(time_diff, point_diff, t, s) >= 0.0)
      return forward(time_diff, point_diff, t, s);
    if(backward(-time_diff, -point_diff, t, s) >= 0.0)
      return forward(-time_diff, -point_diff, t, s);
    return std::numeric_limits<double>::infinity();
  };
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    if(backward(a, t, s) >= 0)
      return forward(a, t, s);
    if(backward(-a, t, s) >= 0)
      return forward(-a, t, s);
    return std::numeric_limits<double>::infinity();
  };
};



template <typename Topology, typename RandomNumberGenerator = boost::minstd_rand>
class reachability_space : public temporal_space<Topology, reachable_distance, RandomNumberGenerator> {
  public:
    typedef temporal_space<Topology, reachable_distance, RandomNumberGenerator> temporal_space_type;
    typedef temporal_space_type::time_topology time_topology;
    typedef temporal_space_type::space_topology space_topology;
    typedef temporal_space_type::distance_metric distance_metric;
    
    typedef temporal_space_type::point_type point_type;
    typedef temporal_space_type::point_difference_type point_difference_type;
    
    
    double forward_reach(const point_type& p) const {
      return forward_reachable_norm()(this->difference(p, this->origin()), this->time, this->space);
    };
    
    double backward_reach(const point_type& p) const {
      return backward_reachable_norm()(this->difference(p, this->origin()), this->time, this->space);
    };
    
    double forward_norm(const point_difference_type& p) const {
      return forward_reachable_norm()(p, this->time, this->space);
    };
    
    double backward_norm(const point_difference_type& p) const {
      return backward_reachable_norm()(p, this->time, this->space);
    };
    
};






};

};

#endif





















