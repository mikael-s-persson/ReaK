/**
 * \file svp_Ndof_samplers.h
 *
 * This library provides an implementation of a random-sampler within a differentiable N-dof topology
 * which respects the boundaries of the topology within the context of sustained velocity motions
 * that can stop the motions before reaching the boundaries of the topology.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SVP_NDOF_SAMPLERS_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SVP_NDOF_SAMPLERS_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/prob_distribution_concept.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/time_topology.h"

#include "ReaK/topologies/interpolation/sustained_velocity_pulse_ndof.h"

#include <cmath>
#include <utility>

namespace ReaK::pp {

/**
 * This functor class is a random-sampler based on the rate-limited motions of a SVP interpolation
 * between points within a bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <Topology TimeSpaceType = time_topology>
struct svp_Ndof_rate_limited_sampler : public serializable {
  using self = svp_Ndof_rate_limited_sampler<TimeSpaceType>;

  std::shared_ptr<TimeSpaceType> t_space;

  explicit svp_Ndof_rate_limited_sampler(
      std::shared_ptr<TimeSpaceType> aTimeSpace)
      : t_space(std::move(aTimeSpace)) {}

  svp_Ndof_rate_limited_sampler()
      : svp_Ndof_rate_limited_sampler(std::make_shared<TimeSpaceType>()) {}

  /**
   * This function returns a random sample-point on a topology.
   * \param s The topology or space on which the sample-point lies.
   * \return A random sample-point on the topology.
   */
  template <PointDistribution Space>
  topology_point_type_t<Space> operator()(const Space& s) const {
    static_assert(TangentBundle<Space, TimeSpaceType, 0, 1>);

    const auto& get_sample = get(random_sampler, s);

    while (true) {
      auto pt = get_sample(s);

      if (svp_Ndof_is_in_bounds(pt, s, *t_space)) {
        return pt;
      }
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(t_space);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(t_space);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2450003, 1,
                              "svp_Ndof_rate_limited_sampler", serializable)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SVP_NDOF_SAMPLERS_H_
