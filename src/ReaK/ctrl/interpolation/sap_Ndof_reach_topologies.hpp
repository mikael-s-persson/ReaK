/**
 * \file sap_Ndof_reach_topologies.hpp
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

#ifndef REAK_SAP_NDOF_REACH_TOPOLOGIES_HPP
#define REAK_SAP_NDOF_REACH_TOPOLOGIES_HPP

#include "base/defs.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"
#include "path_planning/tangent_bundle_concept.hpp"

#include "interpolated_topologies.hpp"
#include "generic_interpolator_factory.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "topologies/basic_distance_metrics.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "topologies/generic_sampler_factory.hpp"
#include "topologies/rate_limited_spaces.hpp"

#include "sustained_acceleration_pulse_Ndof.hpp"
#include "sap_Ndof_metrics.hpp"
#include "sap_Ndof_samplers.hpp"

namespace ReaK {

namespace pp {



/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class interpolated_topology<BaseTopology, sap_Ndof_interpolation_tag> : public interpolated_topology_base<BaseTopology> {
  public:
    
    typedef interpolated_topology_base<BaseTopology> base_type;
    typedef interpolated_topology<BaseTopology, sap_Ndof_interpolation_tag> self;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef typename base_type::distance_metric_type distance_metric_type;
    typedef typename base_type::random_sampler_type random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
    
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
    typedef boost::function< bool(const point_type&) > validity_predicate_type;
#else
    typedef std::function< bool(const point_type&) > validity_predicate_type;
#endif
    
  protected:
    
    shared_ptr<time_topology> t_space;
    sap_Ndof_reach_time_metric<time_topology> rt_dist;
    generic_sampler<sap_Ndof_rate_limited_sampler<time_topology>, BaseTopology> rl_sampler;
    
    virtual point_type interp_topo_move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      try {
        detail::generic_interpolator_impl<sap_Ndof_interpolator,BaseTopology,time_topology> interp;
        interp.initialize(a, b, 0.0, b_space, *t_space, rt_dist);
        double dt_min = interp.get_minimum_travel_time();
        double dt = dt_min * fraction;
        point_type result = a;
        interp.compute_point(result, a, b, b_space, *t_space, dt, dt_min, rt_dist);
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    virtual double interp_topo_get_distance(const point_type& a, const point_type& b) const {
      return rt_dist(a, b, static_cast<const BaseTopology&>(*this));
    };
    
    virtual point_type interp_topo_move_position_toward_pred(const point_type& a, double fraction, const point_type& b,
                                                        double min_dist_interval, validity_predicate_type predicate) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      try {
        detail::generic_interpolator_impl<sap_Ndof_interpolator,BaseTopology,time_topology> interp;
        interp.initialize(a, b, 0.0, b_space, *t_space, rt_dist);
        double dt_min = interp.get_minimum_travel_time();
        double dt = dt_min * fraction;
        double d = min_dist_interval;
        point_type result = a;
        point_type last_result = a;
        while(d < dt) {
          interp.compute_point(result, a, b, b_space, *t_space, d, dt_min, rt_dist);
          if(!predicate(result))
            return last_result;
          d += min_dist_interval;
          last_result = result;
        };
        if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
          return b;
        if(fraction == 0.0)
          return a;
        interp.compute_point(result, a, b, b_space, *t_space, dt, dt_min, rt_dist);
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    virtual double interp_topo_get_norm(const point_difference_type& dp) const {
      return rt_dist(dp, static_cast<const BaseTopology&>(*this));
    };
    
    virtual bool interp_topo_is_in_bounds(const point_type& a) const {
      return sap_Ndof_is_in_bounds(a, static_cast<const BaseTopology&>(*this), *t_space);
    };
    
    virtual point_type interp_topo_random_point() const {
      return rl_sampler(static_cast<const BaseTopology&>(*this));
    };
    
    
  public:
    
    const sap_Ndof_reach_time_metric<time_topology>& get_pseudo_factory() const { return rt_dist; };
    
    
    interpolated_topology(const BaseTopology& aTopo) : base_type(aTopo),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    interpolated_topology(Args&&... args) : base_type(std::forward<Args>(args)...),
                                            t_space(new time_topology), rt_dist(t_space), 
                                            rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
#else
    interpolated_topology() : base_type(),
                              t_space(new time_topology), rt_dist(t_space), 
                              rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1>
    interpolated_topology(const A1& a1) : 
                          base_type(a1),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2>
    interpolated_topology(const A1& a1, const A2& a2) : 
                          base_type(a1, a2),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3) : 
                          base_type(a1, a2, a3),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                          base_type(a1, a2, a3, a4),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                          base_type(a1, a2, a3, a4, a5),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                          base_type(a1, a2, a3, a4, a5, a6),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7, a8),
                          t_space(new time_topology), rt_dist(t_space), 
                          rl_sampler(sap_Ndof_rate_limited_sampler<time_topology>(t_space)) { };
#endif
   
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC240003A,1,"interpolated_topology",base_type)
    
};


  
/**
 * This class wraps a reach-time topology with SAP-based distance metric and a sampler.
 * \tparam BaseTopology The topology underlying this space, should express values as reach-time values and metrics (distance), and should model TopologyConcept, PointDistributionConcept, BoundedSpaceConcept and TangentBundleConcept for time_topology and up to 2nd order (acceleration).
 */
template <typename BaseTopology>
class sap_Ndof_reach_topology : public interpolated_topology<BaseTopology, sap_Ndof_interpolation_tag>
{
  public:
    BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<BaseTopology>));
    
    typedef interpolated_topology<BaseTopology, sap_Ndof_interpolation_tag> base_type;
    typedef sap_Ndof_reach_topology<BaseTopology> self;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef typename base_type::distance_metric_type distance_metric_type;
    typedef typename base_type::random_sampler_type random_sampler_type;
    
    typedef typename base_type::super_space_type super_space_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
  public:
    
    sap_Ndof_reach_topology(const BaseTopology& aTopo) : base_type(aTopo) { };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    sap_Ndof_reach_topology(Args&&... args) : base_type(std::forward<Args>(args)...) { };
#else
    sap_Ndof_reach_topology() : base_type() { };
    
    template <typename A1>
    sap_Ndof_reach_topology(const A1& a1) : 
                            base_type(a1) { };
    
    template <typename A1, typename A2>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2) : 
                            base_type(a1, a2) { };
    
    template <typename A1, typename A2, typename A3>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3) : 
                            base_type(a1, a2, a3) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                            base_type(a1, a2, a3, a4) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                            base_type(a1, a2, a3, a4, a5) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                            base_type(a1, a2, a3, a4, a5, a6) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                            base_type(a1, a2, a3, a4, a5, a6, a7) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    sap_Ndof_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
                            base_type(a1, a2, a3, a4, a5, a6, a7, a8) { };
#endif
   
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC240002F,1,"sap_Ndof_reach_topology",base_type)
    
};

template <typename BaseTopology>
struct is_metric_space< sap_Ndof_reach_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_point_distribution< sap_Ndof_reach_topology<BaseTopology> > : boost::mpl::true_ { };


template <typename BaseTopology>
struct get_rate_illimited_space< sap_Ndof_reach_topology<BaseTopology> > : 
  get_rate_illimited_space< BaseTopology > { };




template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator< sap_Ndof_interpolation_tag, SpaceType, TimeTopology> {
  typedef detail::generic_interpolator_impl<sap_Ndof_interpolator, SpaceType, TimeTopology> type; 
  typedef sap_Ndof_reach_time_metric<TimeTopology> pseudo_factory_type;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator< sap_Ndof_interpolation_tag, TemporalSpaceType> {
  typedef generic_interpolator<sap_Ndof_interp_factory<TemporalSpaceType>, sap_Ndof_interpolator> type; 
};



};

};


namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename BaseTopology>
  struct arithmetic_tuple_size< pp::sap_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_size< BaseTopology > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, pp::sap_Ndof_reach_topology<BaseTopology> > :
    arithmetic_tuple_element< Idx, BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const pp::sap_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, volatile pp::sap_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, volatile BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const volatile pp::sap_Ndof_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const volatile BaseTopology > { };
  
};





#endif









