/**
 * \file basic_distance_metrics.hpp
 * 
 * This library defines the basic distance metrics classes that work on points of a topology.
 * All the classes satisfy the DistanceMetricConcept.
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

#ifndef REAK_BASIC_DISTANCE_METRICS_HPP
#define REAK_BASIC_DISTANCE_METRICS_HPP

#include "base/serializable.hpp"

#include "path_planning/metric_space_concept.hpp"

namespace ReaK {

namespace pp {
  
/**
 * This class is the default distance metric functor which models the DistanceMetricConcept.
 * This class will simply rely on the distance and norm functions included in the 
 * given topology (assuming it models the MetricSpaceConcept).
 * \note Do not use this distance metric to define a topology, because it will be cyclic (infinite recursion).
 */
struct default_distance_metric : public serialization::serializable {
  
  default_distance_metric() { };
  
  /** 
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    return s.distance(a, b);
  };
  /** 
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    return s.norm(a);
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(default_distance_metric,0xC2410000,1,"default_distance_metric",serialization::serializable)
};


template <typename MetricSpace>
typename boost::enable_if< 
  boost::is_same< typename metric_space_traits<MetricSpace>::distance_metric_type, 
                  default_distance_metric>,
default_distance_metric >::type get(distance_metric_t, const MetricSpace&) {
  return default_distance_metric();
};




/**
 * This class is the default distance metric functor which models the DistanceMetricConcept.
 * This class will simply rely on the distance and norm functions included in the 
 * given topology (assuming it models the MetricSpaceConcept).
 * \note Do not use this distance metric to define a topology, because it will be cyclic (infinite recursion).
 */
template <typename DistanceMetric>
struct symmetrized_metric : public serialization::serializable {
  typedef symmetrized_metric<DistanceMetric> self;
  
  DistanceMetric unsym_distance;
  
  symmetrized_metric(DistanceMetric aUnsymDistance = DistanceMetric()) : unsym_distance(aUnsymDistance) { };
  
  /** 
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    double d_left = unsym_distance(a, b, s);
    double d_right = unsym_distance(b, a, s);
    return (d_left < d_right ? d_left : d_right);
  };
  /** 
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    return unsym_distance(a, s);
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(unsym_distance);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(unsym_distance);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC241000B,1,"symmetrized_metric",serialization::serializable)
};


enum unsymmetrized_metric_t { unsymmetrized_metric };


template <typename DistanceMetric>
struct unsymmetrize {
  typedef DistanceMetric type;
};

template <typename DistanceMetric>
struct unsymmetrize< symmetrized_metric<DistanceMetric> > {
  typedef DistanceMetric type;
};


template <typename DistanceMetric>
const DistanceMetric& get(unsymmetrized_metric_t, const DistanceMetric& d) {
  return d;
};

template <typename DistanceMetric>
DistanceMetric get(unsymmetrized_metric_t, const symmetrized_metric< DistanceMetric >& d) {
  return d.unsym_distance;
};








};

};


#endif


