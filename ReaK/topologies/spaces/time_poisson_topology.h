/**
 * \file time_poisson_topology.h
 *
 * This library provides class that define a time-topology with a Poisson distribution.
 * A time-topology is a simple metric-space where the points are real values (doubles).
 * However, because time is unlimited, this topology uses a Poisson distribution to
 * generate random-points at discrete intervals.
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

#ifndef REAK_TOPOLOGIES_SPACES_TIME_POISSON_TOPOLOGY_H_
#define REAK_TOPOLOGIES_SPACES_TIME_POISSON_TOPOLOGY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/global_rng.h"

#include <cmath>
#include <random>

#include "ReaK/topologies/spaces/time_topology.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"

namespace ReaK::pp {

/**
 * This class implements a time-topology with a Poisson distribution. A time-topology is a
 * simple metric-space where the points are real values (doubles).
 * However, because time is unlimited, this topology uses a Poisson distribution to
 * generate random-points at discrete intervals. This class models the MetricSpaceConcept.
 */
class time_poisson_topology : public time_topology {
 public:
  using point_type = double;
  using point_difference_type = double;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = 1;

  double time_step;
  double mean_discrete_time;
  double time_delay;

  explicit time_poisson_topology(
      const std::string& aName = "time_poisson_topology",
      double aTimeStep = 1.0, double aMeanDiscreteTime = 10.0,
      double aTimeDelay = 0.0)
      : time_topology(aName),
        time_step(aTimeStep),
        mean_discrete_time(aMeanDiscreteTime),
        time_delay(aTimeDelay) {}

  /**
   * Generates a random point in the space, uniformly distributed.
   * \note This function actually returns the origin of the space.
   */
  point_type random_point() const {
    std::lognormal_distribution<double> rand_dist(0.0, 1.0);
    return time_step * int(0.5 * mean_discrete_time / time_step *
                               point_type(rand_dist(get_global_rng())) +
                           time_delay / time_step);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int aVersion) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(time_step) &
        RK_SERIAL_SAVE_WITH_NAME(mean_discrete_time);
    if (aVersion > 1) {
      A& RK_SERIAL_SAVE_WITH_NAME(time_delay);
    }
  }

  void load(serialization::iarchive& A, unsigned int aVersion) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(time_step) &
        RK_SERIAL_LOAD_WITH_NAME(mean_discrete_time);
    if (aVersion > 1) {
      A& RK_SERIAL_LOAD_WITH_NAME(time_delay);
    } else {
      time_delay = 0.0;
    }
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(time_poisson_topology, 0xC240000B, 2,
                              "time_poisson_topology", time_topology)
};

template <>
struct is_metric_space<time_poisson_topology> : std::true_type {};

template <>
struct is_reversible_space<time_poisson_topology> : std::true_type {};

template <>
struct is_point_distribution<time_poisson_topology> : std::true_type {};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_TIME_POISSON_TOPOLOGY_H_
