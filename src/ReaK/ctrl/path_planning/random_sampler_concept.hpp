/**
 * \file random_sampler_concept.hpp
 * 
 * This library defines the traits and concepts that pertain to a random-sampler over 
 * a topology, as used in ReaK::pp. Random-samplers are used in ReaK to generate 
 * points within a topology, often with the hopes of connecting them into a planned path
 * or to propagate or evaluate them in estimation or optimization methods.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2012
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

#ifndef REAK_RANDOM_SAMPLER_CONCEPT_HPP
#define REAK_RANDOM_SAMPLER_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {
  
  
/**
 * This concept defines the requirements to fulfill in order to model a random-sampler 
 * as used in ReaK::pp. A random-sampler is essentially a callable type that can produce 
 * a point within a given topology.
 * 
 * Required concepts:
 * 
 * Topology should model the TopologyConcept.
 * 
 * Valid expressions:
 * 
 * p = rand_sampler(space);  A random point (p) can be obtained by calling the random-sampler (rand_sampler) by providing a const-ref to the topology (space).
 * 
 * \tparam RandomSampler The random-sampler type to be checked for this concept.
 * \tparam Topology The topology to which the random sampler should apply.
 */
template <typename RandomSampler, typename Topology>
struct RandomSamplerConcept {
  RandomSampler rand_sampler;
  Topology space;
  typename topology_traits<Topology>::point_type p;

  BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));
  
  BOOST_CONCEPT_USAGE(RandomSamplerConcept) 
  {
    p = rand_sampler(space);
  };
  
};

/**
 * This tag-type is used to identify (during a "get" call) that the random-sampler object is 
 * to be fetched, if one is associated to a given topology (usually the default one).
 */
enum random_sampler_t { random_sampler };

/**
 * This traits class defines the types and constants associated to a point distribution.
 * \tparam PointDistribution The topology type for which the point distribution traits are sought.
 */
template <typename PointDistribution>
struct point_distribution_traits {
  /** The type that describes the random-sampler type for the distribution. */
  typedef typename PointDistribution::random_sampler_type random_sampler_type;
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
 * rand_sampler = get(random_sampler,space);  The random-sampler can be obtained by a tagged "get" call on the point-distribution.
 * 
 * \tparam PointDistribution The topology type to be checked for this concept.
 */
template <typename PointDistribution>
struct PointDistributionConcept {
  typename topology_traits<PointDistribution>::point_type p1, p2;
  typename point_distribution_traits<PointDistribution>::random_sampler_type rand_sampler;
  PointDistribution space;
  
  BOOST_CONCEPT_ASSERT((TopologyConcept<PointDistribution>));
  BOOST_CONCEPT_ASSERT((RandomSamplerConcept<typename point_distribution_traits<PointDistribution>::random_sampler_type, PointDistribution>));
  
  BOOST_CONCEPT_USAGE(PointDistributionConcept) 
  {
    rand_sampler = get(random_sampler,space);
  };
  
};


template <typename T>
struct is_point_distribution : boost::mpl::false_ { };



};

};


#endif


