/**
 * \file reachability_space.hpp
 * 
 * This library defines a number of classes to realize a reachability space (see ReachabilitySpaceConcept)
 * from a spatial topology with a temporal distance metric (that is, a norm or distance metric which 
 * represents the minimum travel time between two points).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REACHABILITY_SPACE_HPP
#define REACHABILITY_SPACE_HPP

#include "temporal_space.hpp"
#include "reachability_space_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {

/**
 * This functor class implements a backward reachable norm for a given point-difference in 
 * a temporal-space.
 */
struct backward_reachable_norm {
  /**
   * Computes the backward reachable norm for a given point-difference in a temporal-space.
   * \tparam PointDiff The point-difference type from the temporal-space.
   * \tparam TimeTopology The time-topology type.
   * \tparam SpaceTopology The space-topology type, should model MetricSpaceConcept.
   * \param a The point-difference lying on the temporal-space.
   * \param t The time-topology of the temporal-space.
   * \param s The space-topology of the temporal-space.
   * \return The backward reachable norm for the given point-difference.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    return a.time - s.norm(a.pt);
  };
  /**
   * Computes the backward reachable norm for a given point-difference in a temporal-space.
   * \tparam TimeDiff The time-difference type from the temporal-space.
   * \tparam SpaceDiff The space-difference type from the temporal-space.
   * \tparam TimeTopology The time-topology type.
   * \tparam SpaceTopology The space-topology type, should model MetricSpaceConcept.
   * \param dt The time-difference.
   * \param dp The space-difference.
   * \param t The time-topology of the temporal-space.
   * \param s The space-topology of the temporal-space.
   * \return The backward reachable norm for the given point-difference.
   */
  template <typename TimeDiff, typename SpaceDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const TimeDiff& dt, const SpaceDiff& dp, const TimeTopology& t, const SpaceTopology& s) const {
    return dt - s.norm(dp);
  };
};


/**
 * This functor class implements a forward reachable norm for a given point-difference in 
 * a temporal-space.
 */
struct forward_reachable_norm {
  /**
   * Computes the forward reachable norm for a given point-difference in a temporal-space.
   * \tparam PointDiff The point-difference type from the temporal-space.
   * \tparam TimeTopology The time-topology type.
   * \tparam SpaceTopology The space-topology type, should model MetricSpaceConcept.
   * \param a The point-difference lying on the temporal-space.
   * \param t The time-topology of the temporal-space.
   * \param s The space-topology of the temporal-space.
   * \return The forward reachable norm for the given point-difference.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    return a.time + s.norm(a.pt);
  };
  /**
   * Computes the forward reachable norm for a given point-difference in a temporal-space.
   * \tparam TimeDiff The time-difference type from the temporal-space.
   * \tparam SpaceDiff The space-difference type from the temporal-space.
   * \tparam TimeTopology The time-topology type.
   * \tparam SpaceTopology The space-topology type, should model MetricSpaceConcept.
   * \param dt The time-difference.
   * \param dp The space-difference.
   * \param t The time-topology of the temporal-space.
   * \param s The space-topology of the temporal-space.
   * \return The forward reachable norm for the given point-difference.
   */
  template <typename TimeDiff, typename SpaceDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const TimeDiff& dt, const SpaceDiff& dp, const TimeTopology& t, const SpaceTopology& s) const {
    return dt + s.norm(dp);
  };
};


/**
 * This functor class models the TemporalDistMetricConcept via the reachability norms on a temporal
 * space. This gives a metric that satisfies the triangular inequality.
 */
struct reachable_distance {
  forward_reachable_norm  forward;
  backward_reachable_norm backward;
  /**
   * Computes the reachability distance between two points of a temporal space.
   * \tparam Point The point type of the temporal space.
   * \tparam TimeTopology The time-topology type of the temporal space.
   * \tparam SpaceTopology The space-topology type of the temporal space.
   * \param a The first point.
   * \param b The second point.
   * \param t The time-topology.
   * \param s The space-topology.
   * \return The reachability distance between the two points of the temporal space (t,s).
   */
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
  /**
   * Computes the reachability distance of a point-difference in a temporal space.
   * \tparam PointDiff The point-difference type of the temporal space.
   * \tparam TimeTopology The time-topology type of the temporal space.
   * \tparam SpaceTopology The space-topology type of the temporal space.
   * \param a The point-difference.
   * \param t The time-topology.
   * \param s The space-topology.
   * \return The reachability distance of the point-difference in the temporal space (t,s).
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology& s) const {
    if(backward(a, t, s) >= 0)
      return forward(a, t, s);
    if(backward(-a, t, s) >= 0)
      return forward(-a, t, s);
    return std::numeric_limits<double>::infinity();
  };
};


/**
 * This class realizes a reachability space (see ReachabilitySpaceConcept)
 * from a spatial topology with a temporal distance metric (that is, a norm or distance metric which 
 * represents the minimum travel time between two points). This class extends a temporal_space with
 * the additional functions required to model the ReachabilitySpaceConcept.
 * \tparam Topology A spatial topology whose distance metric represents the minimum travel time between two points, should model the MetricSpaceConcept.
 * \tparam RandomNumberGenerator The random number generator functor type that can introduce the randomness needed for generating samples of the space.
 */
template <typename Topology, typename RandomNumberGenerator = boost::minstd_rand>
class reachability_space : public temporal_space<Topology, reachable_distance, RandomNumberGenerator> {
  public:
    typedef temporal_space<Topology, reachable_distance, RandomNumberGenerator> temporal_space_type;
    typedef temporal_space_type::time_topology time_topology;
    typedef temporal_space_type::space_topology space_topology;
    typedef temporal_space_type::distance_metric distance_metric;
    
    typedef temporal_space_type::point_type point_type;
    typedef temporal_space_type::point_difference_type point_difference_type;
    
    /**
     * Returns the forward reach of a given point in the temporal-space.
     * \param p The point for which the forward-reach is requested.
     * \return The forward reach of p.
     */
    double forward_reach(const point_type& p) const {
      return forward_reachable_norm()(this->difference(p, this->origin()), this->time, this->space);
    };
    
    /**
     * Returns the backward reach of a given point in the temporal-space.
     * \param p The point for which the backward-reach is requested.
     * \return The backward reach of p.
     */
    double backward_reach(const point_type& p) const {
      return backward_reachable_norm()(this->difference(p, this->origin()), this->time, this->space);
    };
    
    /**
     * Returns the forward reachable norm of a given point-difference in the temporal-space.
     * \param p The point-difference for which the forward reachable norm is requested.
     * \return The forward reachable norm of p.
     */
    double forward_norm(const point_difference_type& p) const {
      return forward_reachable_norm()(p, this->time, this->space);
    };
    
    /**
     * Returns the backward reachable norm of a given point-difference in the temporal-space.
     * \param p The point-difference for which the backward reachable norm is requested.
     * \return The backward reachable norm of p.
     */
    double backward_norm(const point_difference_type& p) const {
      return backward_reachable_norm()(p, this->time, this->space);
    };
    
};






};

};

#endif





















