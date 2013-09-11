/**
 * \file svp_Ndof_reach_topologies.hpp
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

#ifndef REAK_SVP_NDOF_REACH_TOPOLOGIES_HPP
#define REAK_SVP_NDOF_REACH_TOPOLOGIES_HPP

#include "base/defs.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"
#include "path_planning/tangent_bundle_concept.hpp"

#include "interpolated_trajectory.hpp"
#include "generic_interpolator_factory.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "topologies/generic_sampler_factory.hpp"
#include "topologies/rate_limited_spaces.hpp"

#include "sustained_velocity_pulse_Ndof.hpp"
#include "svp_Ndof_metrics.hpp"
#include "svp_Ndof_samplers.hpp"


namespace ReaK {

namespace pp {



/**
 * This class wraps a reach-time topology with SVP-based distance metric and a sampler.
 * \tparam BaseTopology The topology underlying this space, should express values as reach-time values and metrics (distance), and should model TopologyConcept, PointDistributionConcept, BoundedSpaceConcept and TangentBundleConcept for time_topology and up to 2nd order (acceleration).
 */
template <typename BaseTopology>
class svp_Ndof_reach_topology : public BaseTopology
{
  public:
    BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<BaseTopology>));
    
    typedef svp_Ndof_reach_topology<BaseTopology> self;
    
    typedef typename topology_traits< BaseTopology >::point_type point_type;
    typedef typename topology_traits< BaseTopology >::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef BaseTopology super_space_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< BaseTopology >::dimensions);
    
  protected:
    
    shared_ptr<time_topology> t_space;
    svp_Ndof_reach_time_metric<time_topology> rt_dist;
    generic_sampler<svp_Ndof_rate_limited_sampler<time_topology>, BaseTopology> rl_sampler;
    
    
  public:
    
    const svp_Ndof_reach_time_metric<time_topology>& get_pseudo_factory() const { return rt_dist; };
    
    svp_Ndof_reach_topology(const BaseTopology& aTopo) : 
                            BaseTopology(aTopo),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    svp_Ndof_reach_topology(Args&&... args) : 
                            BaseTopology(std::forward<Args>(args)...),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
#else
    svp_Ndof_reach_topology() : 
                            BaseTopology(),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1>
    svp_Ndof_reach_topology(const A1& a1) : 
                            BaseTopology(a1),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2>
    svp_Ndof_reach_topology(const A1& a1, const A2& a2) : 
                            BaseTopology(a1, a2),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3>
    svp_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3) : 
                            BaseTopology(a1, a2, a3),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    svp_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                            BaseTopology(a1, a2, a3, a4),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    svp_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                            BaseTopology(a1, a2, a3, a4, a5),
                            t_space(new time_topology),
                            rt_dist(t_space), 
                            rl_sampler(svp_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
#endif
                       
   /**
    * Returns a const-reference to the super-space of this topology.
    * \note This function returns a const-reference to itself since the super-space is also 
    *       the base-class of this topology. The base class is not polymorphic, meaning that its
    *       distance metric and random-sampler are not overridden (non-virtual calls).
    */
   const super_space_type& get_super_space() const { return *this; };
    
   /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      return rt_dist(a, b, static_cast<const BaseTopology&>(*this));
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      return rt_dist(delta, static_cast<const BaseTopology&>(*this));
    };
    
    /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      return svp_Ndof_is_in_bounds(a, static_cast<const BaseTopology&>(*this), *t_space);
    };
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed within the reachable space.
     */
    point_type random_point() const {
      return rl_sampler(static_cast<const BaseTopology&>(*this));
    };

   /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      detail::generic_interpolator_impl<svp_Ndof_interpolator,BaseTopology,time_topology> interp;
      interp.initialize(a, b, 0.0, static_cast<const BaseTopology&>(*this), *t_space, rt_dist);
      double dt_min = interp.get_minimum_travel_time();
      double dt = dt_min * fraction;
      point_type result = a;
      interp.compute_point(result, a, b, static_cast<const BaseTopology&>(*this), *t_space, dt, dt_min, rt_dist);
      return result;
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      BaseTopology::save(A,BaseTopology::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      BaseTopology::load(A,BaseTopology::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400030,1,"svp_Ndof_reach_topology",BaseTopology)
    
};


template <typename BaseTopology>
struct is_metric_space< svp_Ndof_reach_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_point_distribution< svp_Ndof_reach_topology<BaseTopology> > : boost::mpl::true_ { };



template <typename BaseTopology>
struct get_rate_illimited_space< svp_Ndof_reach_topology<BaseTopology> > : 
  get_rate_illimited_space< BaseTopology > { };



template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator< svp_Ndof_interpolation_tag, SpaceType, TimeTopology> {
  typedef detail::generic_interpolator_impl<svp_Ndof_interpolator, SpaceType, TimeTopology> type; 
  typedef svp_Ndof_reach_time_metric<TimeTopology> pseudo_factory_type;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator< svp_Ndof_interpolation_tag, TemporalSpaceType> {
  typedef generic_interpolator<svp_Ndof_interp_factory<TemporalSpaceType>, svp_Ndof_interpolator> type; 
};





};


};



namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename BaseTopology>
  struct arithmetic_tuple_size< pp::svp_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_size< BaseTopology > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, pp::svp_Ndof_reach_topology<BaseTopology> > :
    arithmetic_tuple_element< Idx, BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const pp::svp_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, volatile pp::svp_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, volatile BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const volatile pp::svp_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const volatile BaseTopology > { };
  
};





#endif









