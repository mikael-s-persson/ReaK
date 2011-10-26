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
  
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<6> 
    >,
  void >::type quintic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                                const DiffSpace& space, const TimeSpace& t_space,
					        double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
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
    
   // a0 + 
    get<2>(result) = get_space<2>(space,t_space).adjust( get<2>(a), 
   //   (5 t - 15 t^2 + 10 t^3) * (diff( lift( 6 * diff( lift( diff( p1, p0 ) ), v0 ) ), a0 ) + diff( a1, lift( 6 * diff( v1, lift( diff( p1, p0 ) ) ) ) ) )
      (5.0 * (t_normal - 3.0 * t2 + 2.0 * t3)) * da_term1
   //   + (3 t - 3 t^2) * ( diff( lift( diff( v1, v0 ) ), a0 ) - diff( a1, lift( diff( v1, v0 ) ) ) ) + t diff(a1, a0)
      + (3.0 * (t_normal - t2)) * da_term2 + t_normal * da1a0
    );
    
    quintic_hermite_interpolate_HOT_impl< Idx, PointType, PointDiff2, DiffSpace, TimeSpace >(result,da1a0,da_term1,da_term2,space,t_space,t_factor,t_normal);
    
  };
  
  
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<5> 
    >,
  void >::type quintic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                                const DiffSpace& space, const TimeSpace& t_space,
					        double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    quintic_hermite_interpolate_impl< typename boost::mpl::prior<Idx>::type, PointType, DiffSpace, TimeSpace >(result,a,b,space,t_space,t_factor,t_normal);
    
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
  BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 2, typename temporal_topology_traits<Topology>::time_topology >));
  typedef typename temporal_topology_traits< Topology >::space_topology SpaceType;

  double t_factor = b.time - a.time;
  if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
    throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
  double t_normal = (t - a.time) / (b.time - a.time);
      
  PointType result;
  result.time = t;
  detail::quintic_hermite_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result.pt, a.pt, b.pt, space.get_space_topology(), space.get_time_topology(), t_factor, t_normal);
      
  return result;      
};


/**
 * This functor class implements a quintic Hermite interpolation in a temporal and twice-differentiable 
 * topology.
 */
struct quintic_hermite_interpolator {
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
  PointType operator()(const PointType& a, const PointType& b, double t, const Topology& space) const {
    return quintic_hermite_interpolate(a,b,t,space);
  };
};




  
/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class quintic_hermite_interp : public interpolated_trajectory<Topology,quintic_hermite_interpolator,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 2, typename temporal_topology_traits<Topology>::time_topology >));
    
    typedef quintic_hermite_interp<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,quintic_hermite_interpolator,DistanceMetric> base_class_type;
    
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
    explicit quintic_hermite_interp(const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                    base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    quintic_hermite_interp(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                           base_class_type(aSpace, aStart, aEnd, aDist) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    quintic_hermite_interp(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                           base_class_type(aBegin, aEnd, aSpace, aDist) { };
    
};



};

};

#endif









