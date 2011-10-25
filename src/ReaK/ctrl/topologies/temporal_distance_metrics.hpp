/**
 * \file temporal_distance_metrics.hpp
 * 
 * This library defines basic temporal distance metrics to use on temporal 
 * spaces (see TemporalSpaceConcept and TemporalDistMetricConcept).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_TEMPORAL_DISTANCE_METRICS_HPP
#define REAK_TEMPORAL_DISTANCE_METRICS_HPP

#include "base/serializable.hpp"

namespace ReaK {

namespace pp {



/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the 
 * distance based only on the distance in the spatial dimensions (space-topology).
 */
struct spatial_distance_only : public serialization::serializable {
  
  spatial_distance_only() { };
  
  /**
   * Computes the distance by calling the distance-function of the space-topology (s) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \param s The space-topology.
   * \return the spatial-distance between the two points.
   */
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology&, const SpaceTopology& s) const {
    return s.distance(a.pt, b.pt);
  };
  /**
   * Computes the norm by calling the norm-function of the space-topology (s) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The point-difference.
   * \param s The space-topology.
   * \return The spatial-norm of the difference between the two points.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology&, const SpaceTopology& s) const {
    return s.norm(a.pt);
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(spatial_distance_only,0xC2410010,1,"spatial_distance_only",serialization::serializable)

};


/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the 
 * distance based only on the distance in the temporal dimensions (time-topology).
 */
struct time_distance_only : public serialization::serializable {
  
  time_distance_only() { };
  
  /**
   * Computes the distance by calling the distance-function of the time-topology (t) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \param t The time-topology.
   * \return the temporal-distance between the two points.
   */
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology& t, const SpaceTopology&) const {
    return t.distance(a.time, b.time);
  };
  /**
   * Computes the norm by calling the norm-function of the time-topology (t) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The point-difference.
   * \param t The time-topology.
   * \return The temporal-norm of the difference between the two points.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology&) const {
    return t.norm(a.time);
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(spatial_distance_only,0xC2410011,1,"spatial_distance_only",serialization::serializable)
  
};


};

};

#endif


