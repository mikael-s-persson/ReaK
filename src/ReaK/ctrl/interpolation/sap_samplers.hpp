/**
 * \file sap_samplers.hpp
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

#ifndef REAK_SAP_SAMPLERS_HPP
#define REAK_SAP_SAMPLERS_HPP

#include "base/defs.hpp"

#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/rate_limited_spaces.hpp"

#include "sustained_acceleration_pulse.hpp"
#include "topologies/time_topology.hpp"

namespace ReaK {

namespace pp {


/**
 * This functor class is a random-sampler based on the rate-limited motions of a SAP interpolation 
 * between points within a bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology>
struct sap_rate_limited_sampler : public serialization::serializable {
  
  typedef sap_rate_limited_sampler<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  sap_rate_limited_sampler(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
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
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 2, TimeSpaceType>));
    
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 0>::type Space0;
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 1>::type Space1;
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 2>::type Space2;
    typedef typename topology_traits<Topology>::point_type PointType;
    typedef typename topology_traits<Space0>::point_type Point0;
    typedef typename topology_traits<Space1>::point_type Point1;
    typedef typename topology_traits<Space2>::point_type Point2;
    typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
    typedef typename topology_traits<Space2>::point_difference_type PointDiff2;
    
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< Space1 >));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< Space2 >));
    
    const typename point_distribution_traits<Topology>::random_sampler_type& get_sample = get(random_sampler,s);
    //const Space0& s0 = get_space<0>(s, *t_space);
    const Space1& s1 = get_space<1>(s, *t_space);
    const Space2& s2 = get_space<2>(s, *t_space);
    const typename metric_space_traits< Space1 >::distance_metric_type& get_dist1 = get(distance_metric,s1);
    const typename metric_space_traits< Space2 >::distance_metric_type& get_dist2 = get(distance_metric,s2);
    
    while(true) {
      PointType pt = get_sample(s);
      get<2>(pt) = s2.origin();   // the acceleration value should always be 0 in SAP interpolation end-points.
      
      // Get the descended acceleration to the point of zero velocity (origin).
      PointDiff1 dp1 = s1.difference(s1.origin(), get<1>(pt));
      double dt1 = get_dist1(s1.origin(), get<1>(pt), s1);
      // Get the corresponding acceleration.
      Point2 p2 = lift_to_space<2>(dp1, dt1, s, *t_space);
      // Get the descended jerk to reach that acceleration.
      PointDiff2 dp2 = s2.difference(p2, s2.origin());
      double dt2 = get_dist2(p2, s2.origin(), s2);
      
      // Check if we can safely ramp-up to that acceleration:
      PointType result_a = pt;
      detail::sap_constant_jerk_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result_a, dp2, s, *t_space, dt2);
      if( !s.is_in_bounds(result_a) )
        continue; //reject the sample.
      
      if( dt1 > get_dist1(get<1>(result_a), get<1>(pt), s1) ) {
        // This means, we didn't cross the zero-velocity during the jerk-down.
        
        //Get the new descended acceleration to the point of zero velocity (origin).
        dp1 = s1.difference(s1.origin(), get<1>(result_a));
        double dt1a = get_dist1(s1.origin(), get<1>(result_a), s1);
      
        // Check if we can safely stop before the boundary:
        detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result_a, dp1, s, *t_space, dt1a);
        if( !s.is_in_bounds(result_a) )
          continue; //reject the sample.
      };
      
      // Check if we could have ramped-down from that inverse acceleration:
      PointType result_b = pt;
      detail::sap_constant_jerk_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result_b, -dp2, s, *t_space, -dt2);
      if( !s.is_in_bounds(result_b) )
        continue; //reject the sample.
      
      if( dt1 > get_dist1(get<1>(pt), get<1>(result_b), s1) ) {
        // This means, the zero-velocity point is not in the wake of the last jerk-down.
        
        //Get the new descended acceleration to the point of zero velocity (origin).
        dp1 = s1.difference(s1.origin(), get<1>(result_b));
        double dt1b = get_dist1(s1.origin(), get<1>(result_b), s1);
      
        // Check if we could have ramped up from within the boundary:
        detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result_b, -dp1, s, *t_space, -dt1b);
        if( !s.is_in_bounds(result_b) )
          continue; //reject the sample.
      };
      
      // if this point is reached it means that the sample is acceptable:
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2450002,1,"sap_rate_limited_sampler",serialization::serializable)
};




};

};



#endif









