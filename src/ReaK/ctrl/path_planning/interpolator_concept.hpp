/**
 * \file interpolator_concept.hpp
 * 
 * This library defines the traits and concepts related to an interpolator and interpolator factories. An 
 * interpolator is simply an object that stores some representation of a trajectory
 * segment between two points. The basic scheme used in ReaK is to have one interpolator factory object 
 * that serves to generate fresh interpolators for a given trajectory segment. In other words,
 * from a user's perspective, the factory can be used to store certain data that all interpolators will 
 * need (e.g. information about motion limits or general parameters of the interpolation scheme). The 
 * intended use from within the library is to store one interpolator factory object (set by the user) 
 * and use that factory object to create interpolators for each segment of a trajectory.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_INTERPOLATOR_CONCEPT_HPP
#define REAK_INTERPOLATOR_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include "temporal_space_concept.hpp"

#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {

  
/**
 * This traits class defines the traits that characterize an interpolator factory
 * class.
 * \tparam InterpolatorFactory The trajectory for which the traits are sought.
 */
template <typename InterpolatorFactory>
struct interpolator_factory_traits {
  /** This type describes a point in the temporal space or topology. */
  typedef typename InterpolatorFactory::point_type point_type;
  
  /** This type is the interpolator type that is generated from the factory. */
  typedef typename InterpolatorFactory::interpolator_type interpolator_type;
};


/**
 * This concept class defines the requirements for a class to model an interpolator
 * concept as used in ReaK::pp. An interpolator is simply an object that stores some 
 * representation of an interpolated segment between two points.
 * 
 * Required concepts:
 * 
 * The topology should model the TemporalSpaceConcept.
 * 
 * The distance-metric should model the TemporalDistMetricConcept.
 * 
 * Valid expressions:
 * 
 * interp.set_segment(&pt, &pt);  The start and end point of the interpolated segment can be set (as a const pointers to a point).
 * 
 * const point_type* ppt = interp.get_start_point();  The starting point of the interpolated segment can be obtained as a const pointer.
 * 
 * const point_type* ppt = interp.get_end_point();  The ending point of the interpolated segment can be obtained as a const pointer.
 * 
 * d = interp.travel_distance_to(pt);  The travel distance, along the interpolated segment (interp), between the starting point and a given point (pt), can be obtained.
 * 
 * d = interp.travel_distance_from(pt);  The travel distance, along the interpolated segment (interp), between a given point (pt) and the ending point, can be obtained.
 * 
 * pt = interp.get_point_at_time(t);  The point, along the interpolated segment (interp), at a given time (t) can be obtained.
 * 
 * \tparam Interpolator The trajectory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the trajectory should be able to exist.
 * \tparam DistanceMetric The distance metric type to be used on the topology and along the interpolated segments.
 */
template <typename Interpolator, typename Topology, typename DistanceMetric>
struct InterpolatorConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((DistanceMetricConcept< DistanceMetric, Topology >));
  
  Interpolator interp;
  typename topology_traits<Topology>::point_type pt;
  const typename topology_traits<Topology>::point_type* ppt;
  typedef typename temporal_space_traits<Topology>::time_topology time_topology;
  typename topology_traits<time_topology>::point_type t;
  double d;
  DistanceMetric dist;
  
  BOOST_CONCEPT_USAGE(InterpolatorConcept)
  {
    interp.set_segment(&pt,&pt);
    ppt = interp.get_start_point();
    ppt = interp.get_end_point();
    d = interp.travel_distance_to(pt, dist);
    d = interp.travel_distance_from(pt, dist);
    pt = interp.get_point_at_time(t);
  };
  
};



/**
 * This concept class defines the requirements for a class to model a limited interpolator
 * concept as used in ReaK::pp. A limited interpolator is simply an object that stores some 
 * representation of an interpolated segment between two points and also provides functions 
 * to characterize the feasibility of the interpolated segment with respect to the limits 
 * imposed by the topology on the interpolation (the nature of those limits are internal 
 * to the implementation, but should translate into some feasibility test and some minimum
 * travel time required to make the interpolation feasible, in theory).
 * 
 * Required concepts:
 * 
 * The interpolator should model the InterpolatorConcept with the given topology and distance-metric.
 * 
 * Valid expressions:
 * 
 * dt = interp.get_minimum_travel_time();  The minimum travel time required to reach the end-point can be obtained as a time-difference (interval).
 * 
 * bool b = interp.is_segment_feasible();  The segment can be tested for feasibility, i.e., if it is possible to reach the end-point at its associated time while respecting the limits.
 * 
 * \tparam Interpolator The trajectory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the trajectory should be able to exist.
 * \tparam DistanceMetric The distance metric type to be used on the topology and along the interpolated segments.
 */
template <typename LimitedInterpolator, typename Topology, typename DistanceMetric>
struct LimitedInterpolatorConcept : public InterpolatorConcept<LimitedInterpolator,Topology,DistanceMetric> {
  
  typename topology_traits< typename temporal_space_traits<Topology>::time_topology >::point_difference_type dt;
  
  BOOST_CONCEPT_USAGE(LimitedInterpolatorConcept)
  {
    dt = this->interp.get_minimum_travel_time();
    bool b = this->interp.is_segment_feasible(); RK_UNUSED(b);
  };
  
};



/**
 * This concept class defines the requirements for a class to model an interpolator
 * concept as used in ReaK::pp. An interpolator is simply an object that stores some 
 * representation of an interpolated segment between two points. The basic scheme used 
 * in ReaK is to have one interpolator factory object that serves to generate fresh interpolators 
 * for a given trajectory segment. In other words, from a user's perspective, the factory can 
 * be used to store certain data that all interpolators will need (e.g. information about motion 
 * limits or general parameters of the interpolation scheme). The intended use from within the 
 * library is to store one interpolator factory object (set by the user) and use that factory 
 * object to create interpolators for each segment of a trajectory.
 * 
 * Required concepts:
 * 
 * The topology should model the TemporalSpaceConcept.
 * 
 * The interpolator type should model the InterpolatorConcept on the given topology.
 * 
 * Valid expressions:
 * 
 * interp = interp_fact.create_interpolator(&pt, &pt);  An interpolator object (interp) can be created from the factory object (interp_fact) given a start and end point as pointers (&pt,&pt).
 * 
 * interp_fact.set_temporal_space(pspace);  The temporal space object used by the interpolators can be set as a const shared-pointer to a topology.
 * 
 * \tparam InterpolatorFactory The interpolator factory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the interpolated segments should exist.
 * \tparam DistanceMetric The distance metric type to be used on the topology and along the interpolated segments.
 */
template <typename InterpolatorFactory, typename Topology, typename DistanceMetric>
struct InterpolatorFactoryConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((InterpolatorConcept< typename interpolator_factory_traits<InterpolatorFactory>::interpolator_type, Topology, DistanceMetric>));
  
  InterpolatorFactory interp_fact;
  typename interpolator_factory_traits<InterpolatorFactory>::interpolator_type interp;
  typename interpolator_factory_traits<InterpolatorFactory>::point_type pt;
  shared_ptr<Topology> pspace;
  
  BOOST_CONCEPT_USAGE(InterpolatorFactoryConcept)
  {
    interp = interp_fact.create_interpolator(&pt, &pt);
    interp_fact.set_temporal_space(pspace);
  };
  
};


};

};

#endif














