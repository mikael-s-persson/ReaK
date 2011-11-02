/**
 * \file quintic_hermite_interp.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a quintic Hermite interpolation (quintic hermite spline, or qspline).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_QUINTIC_HERMITE_INTERP_HPP
#define REAK_QUINTIC_HERMITE_INTERP_HPP

#include "path_planning/spatial_trajectory_concept.hpp"

#include "path_planning/differentiable_space_concept.hpp"

#include "interpolated_trajectory.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal_to.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include "lin_alg/mat_num_exceptions.hpp"

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type quintic_hermite_interpolate_HOT_impl(PointType&, const PointDiff2&,
						    const PointDiff2&, const PointDiff2&, 
                                                    const DiffSpace&, const TimeSpace&,
					            double, double) {
    /* nothing to do, but this function overload should be available, i.e., it is a no-op. */
  };
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type quintic_hermite_interpolate_HOT_impl(PointType& result, const PointDiff2& da1a0,
						    const PointDiff2& da_term1, const PointDiff2& da_term2, 
                                                    const DiffSpace& space, const TimeSpace& t_space,
					            double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
        
   // lift( 
    get<3>(result) = lift_to_space<3>(
   //   (5 - 30 t + 30 t^2) * ( diff( lift( 6 * diff( lift( diff( p1, p0 ) ), v0 ) ), a0 ) + diff( a1, lift( 6 * diff( v1, lift( diff( p1, p0 ) ) ) ) ) ) 
      (5.0 - (30.0 - 30.0 * t_normal) * t_normal) * da_term1
   //   + (3 - 6 t) * ( diff( lift( diff( v1, v0 ) ), a0 ) - diff( a1, lift( diff( v1, v0 ) ) ) ) + diff(a1, a0)
      + (3.0 - 6.0 * t_normal) * da_term2 + da1a0, t_factor, space, t_space);

  };
  
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<4> 
    >,
  void >::type quintic_hermite_interpolate_HOT_impl(PointType& result, const PointDiff2& da1a0,
						    const PointDiff2& da_term1, const PointDiff2& da_term2, 
                                                    const DiffSpace& space, const TimeSpace& t_space,
					            double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    quintic_hermite_interpolate_HOT_impl< boost::mpl::size_t<3>, PointType, PointDiff2, DiffSpace, TimeSpace >(result,da1a0,da_term1,da_term2,space,t_space,t_factor,t_normal);
      
   // lift( diff(
    get<4>(result) = lift_to_space<4>( get_space<3>(space,t_space).difference(
   //   (-30 + 60 t) * (
      lift_to_space<3>((60.0 * t_normal - 30.0) * da_term1, t_factor, space, t_space),
   //   , 6 * (
      lift_to_space<3>(6.0 * da_term2, t_factor, space, t_space)
    ), t_factor, space, t_space);
    
  };
  
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<5> 
    >,
  void >::type quintic_hermite_interpolate_HOT_impl(PointType& result, const PointDiff2& da1a0,
						    const PointDiff2& da_term1, const PointDiff2& da_term2, 
                                                    const DiffSpace& space, const TimeSpace& t_space,
					            double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    quintic_hermite_interpolate_HOT_impl< boost::mpl::size_t<4>, PointType, PointDiff2, DiffSpace, TimeSpace >(result,da1a0,da_term1,da_term2,space,t_space,t_factor,t_normal);
    
   // lift( diff(
    get<5>(result) = lift_to_space<5>( get_space<4>(space,t_space).difference(
   //   lift( 60 * diff(
       lift_to_space<4>( 60.0 * get_space<3>(space,t_space).difference(   
   //     lift( diff( lift( 6 * diff( lift( diff( p1, p0 ) ), v0 ) ), a0 )
        lift_to_space<3>( da_term1, t_factor, space, t_space), 
   //     , lift( diff( lift( 6 * diff( v1, lift( diff( p1, p0 ) ) ) ), a1 ) )
        get_space<3>(space,t_space).origin() ), t_factor, space, t_space),
   //   , origin<4>
      get_space<4>(space,t_space).origin()
    ), t_factor, space, t_space);
    
  };
  
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointDiff1, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<6> 
    >,
  void >::type quintic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
						const PointDiff0& dp1p0, const PointDiff1& dv1v0,
						const PointDiff1& d_ldp1p0_v0, const PointDiff1& d_v1_ldp1p0,
						const PointDiff1& i_a0, const PointDiff1& i_a1,
						const PointDiff2& da1a0, const PointDiff2& da_term1, const PointDiff2& da_term2,
                                                const DiffSpace& space, const TimeSpace& t_space,
					        double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    /*
    typedef typename derived_N_order_space<DiffSpace,TimeSpace,0>::type Space0;
    typedef typename derived_N_order_space<DiffSpace,TimeSpace,1>::type Space1;
    typedef typename derived_N_order_space<DiffSpace,TimeSpace,2>::type Space2;
    
    typedef typename metric_topology_traits<Space0>::point_type PointType0;
    typedef typename metric_topology_traits<Space1>::point_type PointType1;
    typedef typename metric_topology_traits<Space2>::point_type PointType2;
    
    typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff1;
    typedef typename metric_topology_traits<Space2>::point_difference_type PointDiff2;
    
    PointDiff0 dp1p0 = get_space<0>(space,t_space).difference( get<0>(b), get<0>(a) );
    PointDiff1 dv1v0 = get_space<1>(space,t_space).difference( get<1>(b), get<1>(a) );
    PointType1 ldp1p0 = lift_to_space<1>(dp1p0, t_factor, space, t_space);
    PointDiff1 d_ldp1p0_v0 = get_space<1>(space,t_space).difference( ldp1p0, get<1>(a) );
    PointDiff1 d_v1_ldp1p0 = get_space<1>(space,t_space).difference( get<1>(b), ldp1p0 );
    PointDiff1 i_a0 = descend_to_space<1>(get<2>(a), t_factor, space, t_space);
    PointDiff1 i_a1 = descend_to_space<1>(get<2>(b), t_factor, space, t_space);
    */
    
    double t2 = t_normal * t_normal;
    double t3 = t_normal * t2;
    double t4 = t2 * t2;
    double t5 = t4 * t_normal;
   
   // p0 +
    get<0>(result) = get_space<0>(space,t_space).adjust(get<0>(a), 
   //   (10 t^3 - 15 t^4 + 6 t^5) * (diff(p1,p0) - 0.5 (desc(v0) + desc(v1)) )
      (10.0 * t3 - 15.0 * t4 + 6.0 * t5) * ( dp1p0 - 0.5 * ( descend_to_space<0>(get<1>(a),t_factor, space, t_space) 
	                                                   + descend_to_space<0>(get<1>(b),t_factor, space, t_space) ) )
   //   + t desc(v0 +
      + t_normal * descend_to_space<0>( get_space<1>(space,t_space).adjust(get<1>(a),
   //           (t^2 - 0.5 t^3) diff(v1,v0) + (0.5 t - 1.5 t^2 + 1.5 t^3 - 0.5 t^4) desc(a0) + (0.5 t^2 - t^3 + 0.5 t^4) desc(a1) )
        (t2 - 0.5 * t3) * dv1v0 + (0.5 * (t_normal - 3.0 * t2 + 3.0 * t3 - t4)) * i_a0 + (0.5 * (t2 - 2.0 * t3 + t4)) * i_a1
      ), t_factor, space, t_space));
    
   // v0 + 
    get<1>(result) = get_space<1>(space,t_space).adjust(get<1>(a),
   //   (15 t^2 - 30 t^3 + 15 t^4) * ( diff( lift( diff(p1,p0) ) , v0 ) - diff( v1, lift( diff(p1,p0) ) ) )
      (15.0 * (t2 - 2.0 * t3 + t4)) * (d_ldp1p0_v0 - d_v1_ldp1p0)
   //   (3 t^2 - 2 t^3) diff(v1,v0) + (t - 4.5 t^2 + 6 t^3 - 2.5 t^4) desc(a0) + (1.5 t^2 - 4 t^3 + 2.5 t^4) desc(a1)
      + (3.0 * t2 - 2.0 * t3) * dv1v0 + (t_normal - 4.5 * t2 + 6.0 * t3 - 2.5 * t4) * i_a0 + (1.5 * t2 - 4.0 * t3 + 2.5 * t4) * i_a1 
    );
    /*
    PointType2 l6d_ldp1p0_v0 = lift_to_space<2>( 6.0 * d_ldp1p0_v0, t_factor, space, t_space);
    PointType2 l6d_v1_ldp1p0 = lift_to_space<2>( 6.0 * d_v1_ldp1p0, t_factor, space, t_space);
    PointDiff2 d_l6d_ldp1p0_v0_a0 = get_space<2>(space,t_space).difference( l6d_ldp1p0_v0, get<2>(a) );
    PointDiff2 d_a1_l6d_v1_ldp1p0 = get_space<2>(space,t_space).difference( get<2>(b), l6d_v1_ldp1p0 );
    PointType2 ldv1v0 = lift_to_space<2>( dv1v0, t_factor, space, t_space);
    PointDiff2 d_ldv1v0_a0 = get_space<2>(space,t_space).difference( ldv1v0, get<2>(a) );
    PointDiff2 d_a1_ldv1v0 = get_space<2>(space,t_space).difference( get<2>(b), ldv1v0 );
    PointDiff2 da1a0 = get_space<2>(space,t_space).difference( get<2>(b), get<2>(a) );
    
    PointDiff2 da_term1 = d_l6d_ldp1p0_v0_a0 + d_a1_l6d_v1_ldp1p0;
    PointDiff2 da_term2 = d_ldv1v0_a0 - d_a1_ldv1v0;
    */
   // a0 + 
    get<2>(result) = get_space<2>(space,t_space).adjust( get<2>(a), 
   //   (5 t - 15 t^2 + 10 t^3) * (diff( lift( 6 * diff( lift( diff( p1, p0 ) ), v0 ) ), a0 ) + diff( a1, lift( 6 * diff( v1, lift( diff( p1, p0 ) ) ) ) ) )
      (5.0 * (t_normal - 3.0 * t2 + 2.0 * t3)) * da_term1
   //   + (3 t - 3 t^2) * ( diff( lift( diff( v1, v0 ) ), a0 ) - diff( a1, lift( diff( v1, v0 ) ) ) ) + t diff(a1, a0)
      + (3.0 * (t_normal - t2)) * da_term2 + t_normal * da1a0
    );
    
    quintic_hermite_interpolate_HOT_impl< Idx, PointType, PointDiff2, DiffSpace, TimeSpace >(result,da1a0,da_term1,da_term2,space,t_space,t_factor,t_normal);
    
  };
  
  
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointDiff1, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<5> 
    >,
  void >::type quintic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
						const PointDiff0& dp1p0, const PointDiff1& dv1v0,
						const PointDiff1& d_ldp1p0_v0, const PointDiff1& d_v1_ldp1p0,
						const PointDiff1& i_a0, const PointDiff1& i_a1,
						const PointDiff2& da1a0, const PointDiff2& da_term1, const PointDiff2& da_term2,
                                                const DiffSpace& space, const TimeSpace& t_space,
					        double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    quintic_hermite_interpolate_impl< typename boost::mpl::prior<Idx>::type, PointType, PointDiff0, PointDiff1, PointDiff2, DiffSpace, TimeSpace >(result,a,b,dp1p0,dv1v0,d_ldp1p0_v0,d_v1_ldp1p0,i_a0,i_a1,da1a0,da_term1,da_term2,space,t_space,t_factor,t_normal);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
};



/**
 * This function template computes a quintic Hermite interpolation between two points in a 
 * temporal and twice-differentiable topology.
 * \tparam PointType The point type on the temporal and twice-differentiable topology.
 * \tparam Topology The temporal and twice-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType quintic_hermite_interpolate(const PointType& a, const PointType& b, double t, const Topology& space) {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  typedef typename temporal_topology_traits< Topology >::space_topology SpaceType;
  typedef typename temporal_topology_traits< Topology >::time_topology TimeSpaceType;
  BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< SpaceType, 2, TimeSpaceType >));
  
  double t_factor = b.time - a.time;
  if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
    throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
  double t_normal = (t - a.time) / (b.time - a.time);
      
  PointType result;
  result.time = t;
  
  typedef typename derived_N_order_space<SpaceType,TimeSpaceType,0>::type Space0;
  typedef typename derived_N_order_space<SpaceType,TimeSpaceType,1>::type Space1;
  typedef typename derived_N_order_space<SpaceType,TimeSpaceType,2>::type Space2;
    
  typedef typename metric_topology_traits<Space0>::point_type PointType0;
  typedef typename metric_topology_traits<Space1>::point_type PointType1;
  typedef typename metric_topology_traits<Space2>::point_type PointType2;
    
  typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
  typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff1;
  typedef typename metric_topology_traits<Space2>::point_difference_type PointDiff2;
    
  PointDiff0 dp1p0 = get_space<0>(space.get_space_topology(),space.get_time_topology()).difference( get<0>(b.pt), get<0>(a.pt) );
  PointDiff1 dv1v0 = get_space<1>(space.get_space_topology(),space.get_time_topology()).difference( get<1>(b.pt), get<1>(a.pt) );
  PointType1 ldp1p0 = lift_to_space<1>(dp1p0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff1 d_ldp1p0_v0 = get_space<1>(space.get_space_topology(),space.get_time_topology()).difference( ldp1p0, get<1>(a.pt) );
  PointDiff1 d_v1_ldp1p0 = get_space<1>(space.get_space_topology(),space.get_time_topology()).difference( get<1>(b.pt), ldp1p0 );
  PointDiff1 i_a0 = descend_to_space<1>(get<2>(a.pt), t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff1 i_a1 = descend_to_space<1>(get<2>(b.pt), t_factor, space.get_space_topology(), space.get_time_topology());
  
  PointType2 l6d_ldp1p0_v0 = lift_to_space<2>( 6.0 * d_ldp1p0_v0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointType2 l6d_v1_ldp1p0 = lift_to_space<2>( 6.0 * d_v1_ldp1p0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff2 d_l6d_ldp1p0_v0_a0 = get_space<2>(space.get_space_topology(),space.get_time_topology()).difference( l6d_ldp1p0_v0, get<2>(a.pt) );
  PointDiff2 d_a1_l6d_v1_ldp1p0 = get_space<2>(space.get_space_topology(),space.get_time_topology()).difference( get<2>(b.pt), l6d_v1_ldp1p0 );
  PointType2 ldv1v0 = lift_to_space<2>( dv1v0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff2 d_ldv1v0_a0 = get_space<2>(space.get_space_topology(),space.get_time_topology()).difference( ldv1v0, get<2>(a.pt) );
  PointDiff2 d_a1_ldv1v0 = get_space<2>(space.get_space_topology(),space.get_time_topology()).difference( get<2>(b.pt), ldv1v0 );
  PointDiff2 da1a0 = get_space<2>(space.get_space_topology(),space.get_time_topology()).difference( get<2>(b.pt), get<2>(a.pt) );
    
  PointDiff2 da_term1 = d_l6d_ldp1p0_v0_a0 + d_a1_l6d_v1_ldp1p0;
  PointDiff2 da_term2 = d_ldv1v0_a0 - d_a1_ldv1v0;
  
  detail::quintic_hermite_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result.pt, a.pt, b.pt, dp1p0, dv1v0, d_ldp1p0_v0, d_v1_ldp1p0, i_a0, i_a1, da1a0, da_term1, da_term2, space.get_space_topology(), space.get_time_topology(), t_factor, t_normal);
      
  return result;      
};






/**
 * This functor class implements a quintic Hermite interpolation in a temporal and twice-differentiable 
 * topology.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done.
 */
template <typename Factory>
class quintic_hermite_interpolator {
  public:
    typedef cubic_hermite_interpolator<Factory> self;
    typedef typename Factory::point_type point_type;
    typedef typename Factory::topology topology;
  
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,0>::type Space0;
    typedef typename metric_topology_traits<Space0>::point_type PointType0;
    typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,1>::type Space1;
    typedef typename metric_topology_traits<Space1>::point_type PointType1;
    typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff1;
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,2>::type Space2;
    typedef typename metric_topology_traits<Space1>::point_type PointType2;
    typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff2;
  
  private:
    const Factory* parent;
    const point_type* start_point;
    const point_type* end_point;
    PointDiff0 delta_first_order;
    PointDiff1 delta_second_order;
    PointDiff1 delta_lifted_first_and_second;
    PointDiff1 delta_second_and_lifted_first;
    PointDiff1 delta_integral_third_0;
    PointDiff1 delta_integral_third_1;
    PointDiff2 delta_third_order;
    PointDiff2 third_order_term1;
    PointDiff2 third_order_term2;
    
    void update_delta_value() {
      if(parent && start_point && end_point) {
	double t_factor = end_point->time - start_point->time;
        
	delta_first_order = get_space<0>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<0>(end_point->pt), get<0>(start_point->pt) );
	delta_second_order = get_space<1>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<1>(end_point->pt), get<1>(start_point->pt) );
	
	PointType1 ldp1p0 = lift_to_space<1>(delta_first_order, t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	delta_lifted_first_and_second = get_space<1>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( ldp1p0, get<1>(start_point->pt));
	delta_second_and_lifted_first = get_space<1>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<1>(end_point->pt), ldp1p0 );
	
        delta_integral_third_0 = descend_to_space<1>(get<2>(start_point->pt), t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
        delta_integral_third_1 = descend_to_space<1>(get<2>(end_point->pt), t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
  
        PointType2 l6d_ldp1p0_v0 = lift_to_space<2>( 6.0 * delta_lifted_first_and_second, t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
        PointType2 l6d_v1_ldp1p0 = lift_to_space<2>( 6.0 * delta_second_and_lifted_first, t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
        PointType2 ldv1v0 = lift_to_space<2>( delta_second_order, t_factor, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
        
	delta_third_order = get_space<2>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<2>(end_point->pt), get<2>(start_point->pt) );
        third_order_term1 = get_space<2>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( l6d_ldp1p0_v0, get<2>(start_point->pt) )
                          + get_space<2>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<2>(end_point->pt), l6d_v1_ldp1p0 );
        third_order_term2 = get_space<2>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( ldv1v0, get<2>(start_point->pt) )
                          - get_space<2>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()).difference( get<2>(end_point->pt), ldv1v0 );
      };
    };
  
  public:
    
    
    /**
     * Default constructor.
     */
    quintic_hermite_interpolator(const Factory* aParent = NULL, const point_type* aStart = NULL, const point_type* aEnd = NULL) :
                                 parent(aParent), start_point(aStart), end_point(aEnd) {
      update_delta_value();
    };
    
    void set_segment(const point_type* aStart, const point_type* aEnd) {
      start_point = aStart;
      end_point = aEnd;
      update_delta_value();
    };
    
    const point_type* get_start_point() const { return start_point; };
    const point_type* get_end_point() const { return end_point; };
    
    template <typename DistanceMetric>
    double travel_distance_to(const point_type& pt, const DistanceMetric& dist) const {
      BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,topology>));
      if(parent && start_point)
	return dist(pt, *start_point, *(parent->get_temporal_space()));
      else
	return 0.0;
    };
    
    template <typename DistanceMetric>
    double travel_distance_from(const point_type& pt, const DistanceMetric& dist) const {
      BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,topology>));
      if(parent && end_point)
	return dist(*end_point, pt, *(parent->get_temporal_space()));
      else
	return 0.0;
    };
    
    point_type get_point_at_time(double t) const {
      if(!parent || !start_point || !end_point)
	return point_type();
      double t_factor = end_point->time - start_point->time;
      if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
        throw singularity_error("Normalizing factor in quintic Hermite spline is zero!");
      double t_normal = (t - start_point->time) / t_factor;
      
      point_type result;
      result.time = t;
      
      detail::quintic_hermite_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< typename temporal_topology_traits<topology>::space_topology >::order> >(result.pt, start_point->pt, end_point->pt, delta_first_order, delta_second_order, delta_lifted_first_and_second, delta_second_and_lifted_first, delta_integral_third_0, delta_integral_third_1, delta_third_order, third_order_term1, third_order_term2, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology(), t_factor, t_normal);
  
      return result;   
    };
    
};



/**
 * This class is a factory class for quintic Hermite interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template <typename TemporalTopology>
class quintic_hermite_interp_factory : public serialization::serializable {
  public:
    typedef quintic_hermite_interp_factory<TemporalTopology> self;
    typedef TemporalTopology topology;
    typedef typename temporal_topology_traits<TemporalTopology>::point_type point_type;
    typedef quintic_hermite_interpolator<self> interpolator_type;
  
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<TemporalTopology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<TemporalTopology>::space_topology, 2, typename temporal_topology_traits<TemporalTopology>::time_topology >));
  private:
    const topology* p_space;
  public:
    quintic_hermite_interp_factory(const topology* aPSpace = NULL) : p_space(aPSpace) { };
  
    void set_temporal_space(const topology* aPSpace) { p_space = aPSpace; };
    const topology* get_temporal_space() const { return p_space; };
  
    interpolator_type create_interpolator(const point_type* pp1, const point_type* pp2) const {
      return interpolator_type(this, pp1, pp2);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2420003,1,"quintic_hermite_interp_factory",serialization::serializable)
};



  
/**
 * This class implements a trajectory in a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a quintic Hermite interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class quintic_hermite_interp_traj : public interpolated_trajectory<Topology,quintic_hermite_interp_factory<Topology>,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 2, typename temporal_topology_traits<Topology>::time_topology >));
    
    typedef quintic_hermite_interp_traj<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,quintic_hermite_interp_factory<Topology>,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    typedef typename base_class_type::point_type point_type;
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit quintic_hermite_interp_traj(const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                         base_class_type(aSpace, aDist, quintic_hermite_interp_factory<Topology>(&aSpace)) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    quintic_hermite_interp_traj(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                                base_class_type(aSpace, aStart, aEnd, aDist, quintic_hermite_interp_factory<Topology>(&aSpace)) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    quintic_hermite_interp_traj(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                base_class_type(aBegin, aEnd, aSpace, aDist, quintic_hermite_interp_factory<Topology>(&aSpace)) { };
    
};



};

};

#endif









