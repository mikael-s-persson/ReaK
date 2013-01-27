/**
 * \file sap_Ndof_samplers.hpp
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

#ifndef REAK_SAP_NDOF_SAMPLERS_HPP
#define REAK_SAP_NDOF_SAMPLERS_HPP

#include "base/defs.hpp"

#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/rate_limited_spaces.hpp"

#include "sustained_acceleration_pulse_Ndof.hpp"
#include "topologies/time_topology.hpp"

#include <cmath>

namespace ReaK {

namespace pp {


/**
 * This functor class is a random-sampler based on the rate-limited motions of a SAP interpolation 
 * between points within a N-dof bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology>
struct sap_Ndof_rate_limited_sampler : public serialization::serializable {
  
  typedef sap_rate_limited_sampler<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  sap_Ndof_rate_limited_sampler(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
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
    
    using std::fabs;
    
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
    BOOST_CONCEPT_ASSERT((BoundedSpaceConcept< Space0 >));
    BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept< Space1 >));
    BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept< Space2 >));
    BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point0>));
    BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point1>));
    BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point2>));
    
    const typename point_distribution_traits<Topology>::random_sampler_type& get_sample = get(random_sampler,s);
    const Space0& s0 = get_space<0>(s, *t_space);
    const Space1& s1 = get_space<1>(s, *t_space);
    const Space2& s2 = get_space<2>(s, *t_space);
    
    Point1 max_velocity     = s1.get_upper_corner();
    Point2 max_acceleration = s2.get_upper_corner();
    
    while(true) {
      PointType pt = get_sample(s);
      get<2>(pt) = s2.origin();   // the acceleration value should always be 0 in SAP interpolation end-points.
      
      Point0 stopping_point = get<0>(pt);
      Point0 starting_point = get<0>(pt);
      
      
      for(std::size_t i = 0; i < max_velocity.size(); ++i) {
        double dt_vp1_1st = fabs(get<1>(pt)[i]);
        if(dt_vp1_1st <= std::numeric_limits<double>::epsilon()) // no motion to do (v == 0).
          continue;
        // we know that dt_vp_2nd = dt_vp_1st + dt_amax
        double dt_vp1 = dt_vp1_1st - max_acceleration[i];
        if( dt_vp1 < 0.0 ) {
          //means that we don't have time to reach the maximum acceleration:
          double integ = get<1>(pt)[i] * dt_vp1_1st / max_velocity[i];
          stopping_point[i] += integ;
          starting_point[i] -= integ;
        } else {
          double integ = 0.5 * get<1>(pt)[i] * (dt_vp1_1st + max_acceleration[i]) / max_velocity[i];
          stopping_point[i] += integ;
          starting_point[i] -= integ;
        };
      };
      
      // Check if we could have initiated the motion from within the boundary or if we can stop the motion before the boundary.
      if( !s0.is_in_bounds(stopping_point) || !s0.is_in_bounds(starting_point) )
        continue; //reject the sample.
      
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2450004,1,"sap_Ndof_rate_limited_sampler",serialization::serializable)
};




};

};



#endif









