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

#include <ReaK/core/base/defs.hpp>

#include <ReaK/ctrl/topologies/tangent_bundle_concept.hpp>
#include <ReaK/ctrl/topologies/bounded_space_concept.hpp>
#include <ReaK/ctrl/topologies/prob_distribution_concept.hpp>
#include <ReaK/ctrl/topologies/rate_limited_spaces.hpp>
#include <ReaK/ctrl/topologies/time_topology.hpp>

#include "sustained_acceleration_pulse_Ndof.hpp"

#include <boost/concept_check.hpp>

#include <cmath>

namespace ReaK {

namespace pp {


/**
 * This functor class is a random-sampler based on the rate-limited motions of a SAP interpolation 
 * between points within a N-dof bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology>
struct sap_Ndof_rate_limited_sampler : public serializable {
  
  typedef sap_Ndof_rate_limited_sampler<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  sap_Ndof_rate_limited_sampler(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
                                t_space(aTimeSpace) { };
  
  
  template <typename Topology>
  bool is_in_bounds(const Topology& s, const typename topology_traits<Topology>::point_type& pt) const {
    return sap_Ndof_is_in_bounds(pt, s, *t_space);
  };
  
  
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
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 2, TimeSpaceType>));
    
    typedef typename topology_traits<Topology>::point_type PointType;
    
    const typename point_distribution_traits<Topology>::random_sampler_type& get_sample = get(random_sampler,s);
    
    while(true) {
      PointType pt = get_sample(s);
      // the acceleration value should always be 0 in SAP interpolation end-points.
      get<2>(pt) = get_space<2>(s, *t_space).origin();
      
      if( sap_Ndof_is_in_bounds(pt, s, *t_space) )
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2450004,1,"sap_Ndof_rate_limited_sampler",serializable)
};




};

};



#endif









