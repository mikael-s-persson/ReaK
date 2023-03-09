/**
 * \file sap_samplers.h
 *
 *
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SAP_SAMPLERS_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SAP_SAMPLERS_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/prob_distribution_concept.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/time_topology.h"

#include "ReaK/topologies/interpolation/sustained_acceleration_pulse.h"

#include "boost/concept_check.hpp"

namespace ReaK::pp {

/**
 * This functor class is a random-sampler based on the rate-limited motions of a SAP interpolation
 * between points within a bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology>
struct sap_rate_limited_sampler : public serializable {

  using self = sap_rate_limited_sampler<TimeSpaceType>;

  std::shared_ptr<TimeSpaceType> t_space;

  sap_rate_limited_sampler(const std::shared_ptr<TimeSpaceType>& aTimeSpace =
                               std::make_shared<TimeSpaceType>())
      : t_space(aTimeSpace) {}

  /**
   * This function returns a random sample-point on a topology.
   * \tparam Topology The topology.
   * \param s The topology or space on which the sample-point lies.
   * \return A random sample-point on the topology.
   */
  template <typename Topology>
  topology_point_type_t<Topology> operator()(const Topology& s) const {
    BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 2, TimeSpaceType>));

    using PointType = topology_point_type_t<Topology>;

    const auto& get_sample = get(random_sampler, s);

    while (true) {
      PointType pt = get_sample(s);
      get<2>(pt) =
          get_space<2>(s, *t_space)
              .origin();  // the acceleration value should always be 0 in SAP interpolation end-points.

      if (sap_is_in_bounds<Topology, TimeSpaceType>(pt, s, *t_space))
        return pt;
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(t_space);
  }

  void load(serialization::iarchive& A, unsigned int) override {
    A& RK_SERIAL_LOAD_WITH_NAME(t_space);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2450002, 1, "sap_rate_limited_sampler",
                              serializable)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SAP_SAMPLERS_H_
