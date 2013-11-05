/**
 * \file sap_reach_topologies.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a rate-limited sustained acceleration pulse (SAP) interpolation.
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

#ifndef REAK_SAP_REACH_TOPOLOGIES_HPP
#define REAK_SAP_REACH_TOPOLOGIES_HPP

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

#include "sustained_acceleration_pulse.hpp"
#include "sap_metrics.hpp"
#include "sap_samplers.hpp"

#include "optimization/optim_exceptions.hpp"

namespace ReaK {

namespace pp {


namespace detail {
  
  template <typename BaseTopology, bool IsTemporal>
  struct sap_reach_topo_impl {
    
    typedef typename topology_traits<BaseTopology>::point_type point_type;
    typedef typename topology_traits<BaseTopology>::point_difference_type point_difference_type;
    
    typedef sap_reach_time_metric<time_topology> rt_metric_type;
    typedef generic_sampler<sap_rate_limited_sampler<time_topology>, BaseTopology> sampler_type;
    
    static rt_metric_type make_rt_metric(const BaseTopology&) {
      return rt_metric_type();
    };
    
    static sampler_type make_sampler(const BaseTopology&) {
      return sampler_type();
    };
    
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
    typedef boost::function< bool(const point_type&) > validity_predicate_type;
#else
    typedef std::function< bool(const point_type&) > validity_predicate_type;
#endif
    
    static 
    point_type move_pt_toward(const BaseTopology& b_space, 
                              const point_type& a, double fraction, const point_type& b) const {
      try {
        generic_interpolator_impl<sap_interpolator,BaseTopology,time_topology> interp;
        interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
        double dt_min = interp.get_minimum_travel_time();
        double dt = dt_min * fraction;
        point_type result = a;
        interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min, rt_dist);
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    static 
    point_type move_pt_toward(const BaseTopology& b_space, 
                              const point_type& a, double fraction, const point_type& b,
                              double min_dist_interval, validity_predicate_type predicate) const {
      try {
        generic_interpolator_impl<sap_interpolator,BaseTopology,time_topology> interp;
        interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
        double dt_min = interp.get_minimum_travel_time();
        double dt = dt_min * fraction;
        double d = min_dist_interval;
        point_type result = a;
        point_type last_result = a;
        while(d < dt) {
          interp.compute_point(result, a, b, b_space, time_topology(), d, dt_min, rt_dist);
          if(!predicate(result))
            return last_result;
          d += min_dist_interval;
          last_result = result;
        };
        if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
          return b;
        if(fraction == 0.0)
          return a;
        interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min, rt_dist);
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    static 
    double get_distance(const BaseTopology& b_space, const rt_metric_type& rt_dist,
                        const point_type& a, const point_type& b) const {
      return rt_dist(a, b, b_space);
    };
    
    static 
    double get_norm(const BaseTopology& b_space, const rt_metric_type& rt_dist,
                    const point_difference_type& dp) const {
      return rt_dist(dp, b_space);
    };
    
    static 
    bool is_in_bounds(const BaseTopology& b_space,
                      const point_type& a) const {
      return sap_is_in_bounds(a, b_space, time_topology());
    };
    
    static 
    point_type random_point(const BaseTopology& b_space, const sampler_type& rl_sampler) const {
      return rl_sampler(b_space);
    };
    
    
  };
  
  
  
  // Implementation for the temporal spaces:
  template <typename BaseTopology>
  struct sap_reach_topo_impl<BaseTopology, true> {
    
    typedef typename topology_traits<BaseTopology>::point_type point_type;
    typedef typename topology_traits<BaseTopology>::point_difference_type point_difference_type;
    
    typedef typename temporal_space_traits<BaseTopology>::time_topology base_time_topo;
    typedef typename temporal_space_traits<BaseTopology>::space_topology base_space_topo;
    
    typedef sap_reach_time_metric< base_time_topo > rt_metric_type;
    typedef generic_sampler<sap_rate_limited_sampler< base_time_topo >, base_space_topo > sampler_type;
    
    static rt_metric_type make_rt_metric(const BaseTopology& b_space) {
      return rt_metric_type(shared_ptr<const base_time_topo>(&(b_space.get_time_topology()), null_deleter()));
    };
    
    static sampler_type make_sampler(const BaseTopology& b_space) {
      return sampler_type(shared_ptr<const base_time_topo>(&(b_space.get_time_topology()), null_deleter()));
    };
    
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
    typedef boost::function< bool(const point_type&) > validity_predicate_type;
#else
    typedef std::function< bool(const point_type&) > validity_predicate_type;
#endif
    
    point_type move_pt_toward(const BaseTopology& b_space, 
                              const point_type& a, double fraction, const point_type& b) const {
      
      if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
        return a; //b is not reachable from a.
      
      try {
        generic_interpolator_impl<sap_interpolator, base_space_topo, base_time_topo> interp;
        double dt_total = (b.time - a.time);  // the free time that I have along the path.
        interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), rt_dist);
        double dt = dt_total * fraction;
        point_type result = a;
        interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, rt_dist);
        result.time += dt;
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    point_type move_pt_toward(const BaseTopology& b_space, 
                              const point_type& a, double fraction, const point_type& b,
                              double min_dist_interval, validity_predicate_type predicate) const {
      
      if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
        return a; //b is not reachable from a.
      
      try {
        double dt_total = (b.time - a.time);  // the free time that I have along the path.
        if(dt_total < min_dist_interval)
          return move_pt_toward(b_space, a, fraction, b);
        
        generic_interpolator_impl<sap_interpolator, base_space_topo, base_time_topo> interp;
        interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), rt_dist);
        double dt = dt_total * fraction;
        double d = min_dist_interval;
        point_type result = a;
        point_type last_result = a;
        while(d < dt) {
          interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), d, dt_total, rt_dist);
          result.time = a.time + d;
          if(!predicate(result))
            return last_result;
          d += min_dist_interval;
          last_result = result;
        };
        if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
          return b;
        if(fraction == 0.0)
          return a;
        interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, rt_dist);
        result.time = a.time + dt;
        return result;
      } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        return a;
      };
    };
    
    double get_distance(const BaseTopology& b_space, const rt_metric_type& rt_dist,
                        const point_type& a, const point_type& b) const {
      return rt_dist(a.pt, b.pt, b_space.get_space_topology());
    };
    
    double get_norm(const BaseTopology& b_space, const rt_metric_type& rt_dist,
                    const point_difference_type& dp) const {
      return rt_dist(dp.pt, b_space.get_space_topology());
    };
    
    bool is_in_bounds(const BaseTopology& b_space,
                      const point_type& a) const {
      return sap_is_in_bounds(a.pt, b_space.get_space_topology(), b_space.get_time_topology());
    };
    
    point_type random_point(const BaseTopology& b_space, const sampler_type& rl_sampler) const {
      return point_type(get(random_sampler, b_space.get_time_topology())(b_space.get_time_topology()), 
                        rl_sampler(b_space.get_space_topology()));
    };
    
    
  };
  
  template <typename BaseTopology>
  struct sap_reach_topo_selector {
    typedef sap_reach_topo_impl<BaseTopology, is_temporal_space<BaseTopology>::type::value > type;
  };
  
  
};




/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class interpolated_topology<BaseTopology, sap_interpolation_tag> : public interpolated_topology_base<BaseTopology> {
  public:
    
    typedef interpolated_topology_base<BaseTopology> base_type;
    typedef interpolated_topology<BaseTopology, sap_interpolation_tag> self;
    
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
    
    typedef typename detail::sap_reach_topo_selector<BaseTopology>::type Impl;
    typedef typename Impl::rt_metric_type rt_metric_type;
    typedef typename Impl::sampler_type sampler_type;
    
    rt_metric_type rt_dist;
    sampler_type rl_sampler;
    
    virtual point_type interp_topo_move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return Impl::move_pt_toward(*this, a, fraction, b);
    };
    
    virtual point_type interp_topo_move_position_toward_pred(const point_type& a, double fraction, const point_type& b,
                                                             double min_dist_interval, validity_predicate_type predicate) const {
      return Impl::move_pt_toward(*this, a, fraction, b, min_dist_interval, predicate);
    };
    
    virtual double interp_topo_get_distance(const point_type& a, const point_type& b) const {
      return Impl::get_distance(*this, rt_dist, a, b);
    };
    virtual double interp_topo_get_norm(const point_difference_type& dp) const { 
      return Impl::get_norm(*this, rt_dist, dp);
    };
    virtual bool interp_topo_is_in_bounds(const point_type& a) const { 
      return Impl::is_in_bounds(*this, a); 
    };
    virtual point_type interp_topo_random_point() const { 
      return Impl::random_point(*this, rl_sampler); 
    };
    
    
  public:
    
    const rt_metric_type& get_pseudo_factory() const { return rt_dist; };
    
    
    interpolated_topology(const BaseTopology& aTopo) : base_type(aTopo),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    interpolated_topology(Args&&... args) : base_type(std::forward<Args>(args)...),
                                            rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
#else
    interpolated_topology() : base_type(),
                              rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1>
    interpolated_topology(const A1& a1) : 
                          base_type(a1),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2>
    interpolated_topology(const A1& a1, const A2& a2) : 
                          base_type(a1, a2),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3) : 
                          base_type(a1, a2, a3),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                          base_type(a1, a2, a3, a4),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                          base_type(a1, a2, a3, a4, a5),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                          base_type(a1, a2, a3, a4, a5, a6),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7, a8),
                          rt_dist(Impl::make_rt_metric(*this)), rl_sampler(Impl::make_sampler(*this)) { };
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


#if 0

/**
 * This class wraps a reach-time topology with SAP-based distance metric and a sampler.
 * \tparam BaseTopology The topology underlying this space, should express values as reach-time values and metrics (distance), and should model TopologyConcept, PointDistributionConcept, BoundedSpaceConcept and TangentBundleConcept for time_topology and up to 2nd order (acceleration).
 */
template <typename BaseTopology>
class sap_reach_topology : public interpolated_topology<BaseTopology, sap_interpolation_tag>
{
  public:
    BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<BaseTopology>));
    
    typedef interpolated_topology<BaseTopology, sap_interpolation_tag> base_type;
    typedef sap_reach_topology<BaseTopology> self;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef typename base_type::distance_metric_type distance_metric_type;
    typedef typename base_type::random_sampler_type random_sampler_type;
    
    typedef typename base_type::super_space_type super_space_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
  public:
    
    sap_reach_topology(const BaseTopology& aTopo) : base_type(aTopo) { };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    sap_reach_topology(Args&&... args) : base_type(std::forward<Args>(args)...) { };
#else
    sap_reach_topology() : base_type() { };
    
    template <typename A1>
    sap_reach_topology(const A1& a1) : 
                       base_type(a1) { };
    
    template <typename A1, typename A2>
    sap_reach_topology(const A1& a1, const A2& a2) : 
                       base_type(a1, a2) { };
    
    template <typename A1, typename A2, typename A3>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3) : 
                       base_type(a1, a2, a3) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                       base_type(a1, a2, a3, a4) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                       base_type(a1, a2, a3, a4, a5) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                       base_type(a1, a2, a3, a4, a5, a6) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                       base_type(a1, a2, a3, a4, a5, a6, a7) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    sap_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
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

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400022,1,"sap_reach_topology",base_type)
    
};

template <typename BaseTopology>
struct is_metric_space< sap_reach_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_point_distribution< sap_reach_topology<BaseTopology> > : boost::mpl::true_ { };


template <typename BaseTopology>
struct get_rate_illimited_space< sap_reach_topology<BaseTopology> > : 
  get_rate_illimited_space< BaseTopology > { };

#endif


template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator< sap_interpolation_tag, SpaceType, TimeTopology> {
  typedef detail::generic_interpolator_impl<sap_interpolator, SpaceType, TimeTopology> type; 
  typedef sap_reach_time_metric<TimeTopology> pseudo_factory_type;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator< sap_interpolation_tag, TemporalSpaceType> {
  typedef generic_interpolator<sap_interpolator_factory<TemporalSpaceType>, sap_interpolator> type; 
};



};

};

#if 0

namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename BaseTopology>
  struct arithmetic_tuple_size< pp::sap_reach_topology<BaseTopology> > : 
    arithmetic_tuple_size< BaseTopology > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, pp::sap_reach_topology<BaseTopology> > :
    arithmetic_tuple_element< Idx, BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const pp::sap_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, volatile pp::sap_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, volatile BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const volatile pp::sap_reach_topology<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const volatile BaseTopology > { };
  
};

#endif




#endif









