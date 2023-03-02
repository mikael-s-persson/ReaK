/**
 * \file prob_distribution_concept.hpp
 *
 * This library defines the traits and concepts that pertain to a probability distribution over
 * a topology, as used in ReaK::pp. Probability distributions are used in ReaK to assess the
 * probablitity of obtaining a given sample-point from a probability distribution within a topology.
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

#ifndef REAK_PROB_DISTRIBUTION_CONCEPT_HPP
#define REAK_PROB_DISTRIBUTION_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This concept defines the requirements to fulfill in order to model a probability distribution function
 * as used in ReaK::pp. A probability distribution function is essentially a callable type that can compute
 * the probability of a given sample-point from a distribution within a given topology.
 *
 * Required concepts:
 *
 * Topology should model the TopologyConcept.
 *
 * Valid expressions:
 *
 * d = pdf_function(p, space);  A probablitity value (d) can be obtained from a sample-point (p) by calling the
 *pdf-function (pdf_function) by providing a const-ref to the topology (space).
 *
 * \tparam ProbDistFunction The probability distribution function type to be checked for this concept.
 * \tparam Topology The topology on which the probability distribution function should apply.
 */
template <typename ProbDistFunction, typename Topology>
struct ProbDistFunctionConcept {
  ProbDistFunction pdf_function;
  Topology space;
  topology_point_type_t<Topology> p;

  BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));

  BOOST_CONCEPT_USAGE(ProbDistFunctionConcept) {
    double d = pdf_function(p, space);
    RK_UNUSED(d);
  }
};

/**
 * This tag-type is used to identify (during a "get" call) that the probability distribution function object is
 * to be fetched, if one is associated to a given topology (usually the default one).
 */
enum prob_dist_function_t { prob_dist_function };

/**
 * This traits class defines the types and constants associated to a probability distribution.
 * \tparam ProbabilityDistribution The topology type for which the probability distribution traits are sought.
 */
template <typename ProbabilityDistribution>
struct probability_distribution_traits {
  /** The type that describes the probability distribution function type for the distribution. */
  using prob_dist_function_type =
      typename ProbabilityDistribution::prob_dist_function_type;
};

/**
 * This concept defines the requirements to fulfill in order to model a point distribution
 * as used in ReaK::pp. A point distribution is a special kind of topology which has a
 * random-sampler over it.
 *
 * Required concepts:
 *
 * The PointDistribution should model the TopologyConcept.
 *
 * The random-sampler should model the RandomSamplerConcept.
 *
 * Valid expressions:
 *
 * rand_sampler = get(random_sampler, space);  The random-sampler can be obtained by a tagged "get" call on the
 *point-distribution.
 *
 * \tparam ProbabilityDistribution The topology type to be checked for this concept.
 */
template <typename ProbabilityDistribution>
struct ProbabilityDistributionConcept {
  topology_point_type_t<ProbabilityDistribution> p1, p2;
  typename probability_distribution_traits<
      ProbabilityDistribution>::prob_dist_function_type pdf_function;
  ProbabilityDistribution space;

  BOOST_CONCEPT_ASSERT((TopologyConcept<ProbabilityDistribution>));
  BOOST_CONCEPT_ASSERT((ProbDistFunctionConcept<
                        typename probability_distribution_traits<
                            ProbabilityDistribution>::prob_dist_function_type,
                        ProbabilityDistribution>));

  BOOST_CONCEPT_USAGE(ProbabilityDistributionConcept) {
    pdf_function = get(prob_dist_function, space);
  }
};

}  // namespace ReaK::pp

#endif
