/**
 * \file svp_samplers.hpp
 * 
 * This library provides an implementation of a random-sampler within a differentiable topology
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

#ifndef REAK_SVP_SAMPLERS_HPP
#define REAK_SVP_SAMPLERS_HPP

#include "base/defs.hpp"

#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/rate_limited_spaces.hpp"
#include "topologies/time_topology.hpp"

#include "sustained_velocity_pulse.hpp"

namespace ReaK {

namespace pp {


/**
 * This functor class is a random-sampler based on the rate-limited motions of a SVP interpolation 
 * between points within a bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = ReaK::pp::time_topology>
struct svp_rate_limited_sampler : public serialization::serializable {
  
  typedef svp_rate_limited_sampler<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  svp_rate_limited_sampler(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
                           t_space(aTimeSpace) { };
  
  /** 
   * This function returns a random sample-point on a topology.
   * \tparam Topology The topology.
   * \param s The topology or space on which the sample-point lies.
   * \return A random sample-point on the topology.
   */
  template <typename Topology>
  typename topology_traits<Topology>::point_type operator()(const Topology& s) const {
    BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 1, TimeSpaceType>));
    
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 0>::type Space0;
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 1>::type Space1;
    typedef typename topology_traits<Topology>::point_type PointType;
    typedef typename topology_traits<Space0>::point_type Point0;
    typedef typename topology_traits<Space1>::point_type Point1;
    typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
    
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< Space1 >));
    
    const typename point_distribution_traits<Topology>::random_sampler_type& get_sample = get(random_sampler,s);
//    const Space0& s0 = get_space<0>(s, *t_space);
    const Space1& s1 = get_space<1>(s, *t_space);
    const typename metric_space_traits< Space1 >::distance_metric_type& get_vel_dist = get(distance_metric,s1);
    
    while(true) {
      PointType pt = get_sample(s);
      
      PointDiff1 dp1 = s1.difference(s1.origin(), get<1>(pt));
      double dt = get_vel_dist(s1.origin(), get<1>(pt), s1);
      
      // Check if we can stop the motion before the boundary.
      PointType result = pt;
      detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result, dp1, s, *t_space, dt);
      if( !s.is_in_bounds(result) )
        continue; //reject the sample.
      
      // Check if we could have initiated the motion from within the boundary.
      result = pt;
      detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result, -dp1, s, *t_space, -dt);
      if( !s.is_in_bounds(result) )
        continue; //reject the sample.
      
      // If this point is reached, it means that the sample is acceptable:
      return pt;
    };
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(t_space);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(t_space);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(svp_rate_limited_sampler,0xC2450001,1,"svp_rate_limited_sampler",serialization::serializable)
};



};

};



#endif









