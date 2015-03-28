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

#include <ReaK/core/base/serializable.hpp>

#include "temporal_space_concept.hpp"

namespace ReaK {

namespace pp {


/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the
 * distance based only on the distance in the spatial dimensions (space-topology).
 */
struct spatial_distance_only : public serializable {

  spatial_distance_only(){};

  template < typename TemporalTopology >
  spatial_distance_only( const TemporalTopology& ){};

  /**
   * Computes the distance by calling the distance-function of the space-topology (s) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TemporalTopology The temporal-space type associated to the metric, should model TemporalSpaceConcept.
   * \param a The first point.
   * \param b The second point.
   * \param s The temporal-space.
   * \return the spatial-distance between the two points.
   */
  template < typename Point, typename TemporalTopology >
  double operator()( const Point& a, const Point& b, const TemporalTopology& s ) const {
    BOOST_CONCEPT_ASSERT( (TemporalSpaceConcept< TemporalTopology >));
    return get( distance_metric, s.get_space_topology() )( a.pt, b.pt, s.get_space_topology() );
  };
  /**
   * Computes the norm by calling the norm-function of the space-topology (s) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TemporalTopology The temporal-space type associated to the metric, should model TemporalSpaceConcept.
   * \param a The point-difference.
   * \param s The temporal-space.
   * \return The spatial-norm of the difference between the two points.
   */
  template < typename PointDiff, typename TemporalTopology >
  double operator()( const PointDiff& a, const TemporalTopology& s ) const {
    BOOST_CONCEPT_ASSERT( (TemporalSpaceConcept< TemporalTopology >));
    return get( distance_metric, s.get_space_topology() )( a.pt, s.get_space_topology() );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( spatial_distance_only, 0xC2410010, 1, "spatial_distance_only", serializable )
};


/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the
 * distance based only on the distance in the temporal dimensions (time-topology).
 */
struct time_distance_only : public serializable {

  time_distance_only(){};

  template < typename TemporalTopology >
  time_distance_only( const TemporalTopology& ){};

  /**
   * Computes the distance by calling the distance-function of the time-topology (t) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TemporalTopology The temporal-space type associated to the metric, should model TemporalSpaceConcept.
   * \param a The first point.
   * \param b The second point.
   * \param s The temporal-space.
   * \return the temporal-distance between the two points.
   */
  template < typename Point, typename TemporalTopology >
  double operator()( const Point& a, const Point& b, const TemporalTopology& s ) const {
    return get( distance_metric, s.get_time_topology() )( a.time, b.time, s.get_time_topology() );
  };
  /**
   * Computes the norm by calling the norm-function of the time-topology (t) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TemporalTopology The temporal-space type associated to the metric, should model TemporalSpaceConcept.
   * \param a The point-difference.
   * \param s The temporal-space.
   * \return The temporal-norm of the difference between the two points.
   */
  template < typename PointDiff, typename TemporalTopology >
  double operator()( const PointDiff& a, const TemporalTopology& s ) const {
    return get( distance_metric, s.get_time_topology() )( a.time, s.get_time_topology() );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( time_distance_only, 0xC2410011, 1, "time_distance_only", serializable )
};
};
};

#endif
