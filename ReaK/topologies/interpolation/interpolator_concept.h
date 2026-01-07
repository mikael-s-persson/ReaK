/**
 * \file interpolator_concept.h
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATOR_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATOR_CONCEPT_H_


#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include <concepts>

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize an interpolator factory
 * class.
 * \tparam InterpolatorFactory The trajectory for which the traits are sought.
 */
template <typename InterpolatorFactory>
struct interpolator_factory_traits {
  /** This type describes a point in the temporal space or topology. */
  using point_type = typename InterpolatorFactory::point_type;

  /** This type is the interpolator type that is generated from the factory. */
  using interpolator_type = typename InterpolatorFactory::interpolator_type;
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
 * interp.set_segment(&pt, &pt);  The start and end point of the interpolated segment can be set (as a const pointers to
 *a point).
 *
 * const point_type* ppt = interp.get_start_point();  The starting point of the interpolated segment can be obtained as
 *a const pointer.
 *
 * const point_type* ppt = interp.get_end_point();  The ending point of the interpolated segment can be obtained as a
 *const pointer.
 *
 * d = interp.travel_distance_to(pt);  The travel distance, along the interpolated segment (interp), between the
 *starting point and a given point (pt), can be obtained.
 *
 * d = interp.travel_distance_from(pt);  The travel distance, along the interpolated segment (interp), between a given
 *point (pt) and the ending point, can be obtained.
 *
 * pt = interp.get_point_at_time(t);  The point, along the interpolated segment (interp), at a given time (t) can be
 *obtained.
 */
template <typename Interp, typename Space, typename Metric>
concept Interpolator = TemporalSpace<Space>&& DistanceMetric<Metric, Space>&&
requires(Interp& interp, const topology_point_type_t<Space>& pt,
         const Metric& dist,
         const topology_point_type_t<
             typename temporal_space_traits<Space>::time_topology>& t) {
  interp.set_segment(&pt, &pt);
  {
    interp.get_start_point()
    } -> std::convertible_to<const topology_point_type_t<Space>*>;
  {
    interp.get_end_point()
    } -> std::convertible_to<const topology_point_type_t<Space>*>;
  { interp.travel_distance_to(pt, dist) } -> std::convertible_to<double>;
  { interp.travel_distance_from(pt, dist) } -> std::convertible_to<double>;
  {
    interp.get_point_at_time(t)
    } -> std::convertible_to<topology_point_type_t<Space>>;
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
 * The interpolator should model Interpolator with the given topology and distance-metric.
 *
 * Valid expressions:
 *
 * dt = interp.get_minimum_travel_time();  The minimum travel time required to reach the end-point can be obtained as a
 *time-difference (interval).
 *
 * bool b = interp.is_segment_feasible();  The segment can be tested for feasibility, i.e., if it is possible to reach
 *the end-point at its associated time while respecting the limits.
 */
template <typename Interp, typename Space, typename Metric>
concept LimitedInterpolator =
    Interpolator<Interp, Space, Metric>&& requires(const Interp& interp) {
  {
    interp.get_minimum_travel_time()
    } -> std::convertible_to<topology_point_difference_type_t<
        typename temporal_space_traits<Space>::time_topology>>;
  { interp.is_segment_feasible() } -> std::convertible_to<bool>;
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
 * interp = interp_fact.create_interpolator(&pt, &pt);  An interpolator object (interp) can be created from the factory
 *object (interp_fact) given a start and end point as pointers (&pt,&pt).
 *
 * interp_fact.set_temporal_space(pspace);  The temporal space object used by the interpolators can be set as a const
 *shared-pointer to a topology.
 */
template <typename Factory, typename Space, typename Metric>
concept InterpolatorFactory = Interpolator<
    typename interpolator_factory_traits<Factory>::interpolator_type, Space,
    Metric>&&
requires(Factory& factory, std::shared_ptr<Space> pspace,
         const typename interpolator_factory_traits<Factory>::point_type& pt) {
  {
    factory.create_interpolator(&pt, &pt)
    } -> std::convertible_to<
        typename interpolator_factory_traits<Factory>::interpolator_type>;
  factory.set_temporal_space(pspace);
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATOR_CONCEPT_H_
