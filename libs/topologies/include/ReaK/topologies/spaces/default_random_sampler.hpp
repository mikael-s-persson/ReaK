/**
 * \file default_random_sampler.hpp
 *
 * This library defines the default random-sampler class that work on points of a topology.
 * All the classes satisfy the RandomSamplerConcept.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2012
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

#ifndef REAK_DEFAULT_RANDOM_SAMPLER_HPP
#define REAK_DEFAULT_RANDOM_SAMPLER_HPP

#include <ReaK/core/base/serializable.hpp>

#include "random_sampler_concept.hpp"

#include <type_traits>

namespace ReaK::pp {

/**
 * This class is the default random-sampler functor which models the RandomSamplerConcept.
 * This class will simply rely on the random_point() function included in the
 * given topology.
 * \note Do not use this random-sampler to define a topology, because it will be cyclic (infinite recursion).
 */
struct default_random_sampler : public serializable {

  default_random_sampler() = default;

  template <typename Topology>
  explicit default_random_sampler(const Topology& /*unused*/) {}

  /**
   * This function returns a random sample point on a topology.
   * \tparam Topology The topology.
   * \param s The topology or space on which the points lie.
   * \return A random sample point on the given topology.
   */
  template <typename Topology>
  typename topology_traits<Topology>::point_type operator()(
      const Topology& s) const {
    return s.random_point();
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {}

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {}

  RK_RTTI_MAKE_ABSTRACT_1BASE(default_random_sampler, 0xC2450000, 1,
                              "default_random_sampler", serializable)
};

template <typename PointDistribution>
auto get(random_sampler_t /*unused*/, const PointDistribution& s) {
  static_assert(is_point_distribution_v<PointDistribution>);
  return point_distribution_random_sampler_t<PointDistribution>{s};
}

}  // namespace ReaK::pp

#endif
