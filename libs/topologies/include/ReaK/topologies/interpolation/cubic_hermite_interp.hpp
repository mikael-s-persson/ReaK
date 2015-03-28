/**
 * \file cubic_hermite_interp.hpp
 *
 * This library provides an implementation of a trajectory within a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a cubic Hermite interpolation (cubic hermite spline, or cspline).
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

#ifndef REAK_CUBIC_HERMITE_INTERP_HPP
#define REAK_CUBIC_HERMITE_INTERP_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include "spatial_trajectory_concept.hpp"
#include <ReaK/topologies/spaces/tangent_bundle_concept.hpp>

#include "interpolated_trajectory.hpp"
#include "generic_interpolator_factory.hpp"

#include <boost/concept_check.hpp>

#include <cmath>
#include <list>
#include <map>
#include <limits>

namespace ReaK {

namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct cubic_hermite_interpolation_tag {};


namespace detail {

template < typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 2 > >, void >::type
  cubic_hermite_interpolate_HOT_impl( PointType& result, const PointDiff1& dv1v0, const PointDiff1& d_ldp1p0_v0,
                                      const DiffSpace& space, const TimeSpace& t_space, double t_factor,
                                      double t_normal ){/* nothing to do. */
  };

template < typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::equal_to< Idx, boost::mpl::size_t< 2 > >, void >::type
  cubic_hermite_interpolate_HOT_impl( PointType& result, const PointDiff1& dv1v0, const PointDiff1& d_ldp1p0_v0,
                                      const DiffSpace& space, const TimeSpace& t_space, double t_factor,
                                      double t_normal ) {
  get< 2 >( result ) = get_space< 2 >( space, t_space ).adjust(
    lift_to_space< 2 >( dv1v0, t_factor, space, t_space ),
    ( 6.0 - 12.0 * t_normal )
    * get_space< 2 >( space, t_space ).difference( lift_to_space< 2 >( d_ldp1p0_v0, t_factor, space, t_space ),
                                                   lift_to_space< 2 >( 0.5 * dv1v0, t_factor, space, t_space ) ) );
};

template < typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace >
inline typename boost::enable_if< boost::mpl::equal_to< Idx, boost::mpl::size_t< 3 > >, void >::type
  cubic_hermite_interpolate_HOT_impl( PointType& result, const PointDiff1& dv1v0, const PointDiff1& d_ldp1p0_v0,
                                      const DiffSpace& space, const TimeSpace& t_space, double t_factor,
                                      double t_normal ) {
  cubic_hermite_interpolate_HOT_impl< boost::mpl::size_t< 2 >, PointType, PointDiff1, DiffSpace, TimeSpace >(
    result, dv1v0, d_ldp1p0_v0, space, t_space, t_factor, t_normal );

  get< 3 >( result ) = lift_to_space< 3 >(
    -12.0
    * get_space< 2 >( space, t_space ).difference( lift_to_space< 2 >( d_ldp1p0_v0, t_factor, space, t_space ),
                                                   lift_to_space< 2 >( 0.5 * dv1v0, t_factor, space, t_space ) ),
    t_factor, space, t_space );
};


template < typename Idx, typename PointType, typename PointDiff0, typename PointDiff1, typename DiffSpace,
           typename TimeSpace >
inline typename boost::enable_if< boost::mpl::less< Idx, boost::mpl::size_t< 4 > >, void >::type
  cubic_hermite_interpolate_impl( PointType& result, const PointType& a, const PointType& b, const PointDiff0& dp1p0,
                                  const PointDiff1& dv1v0, const PointDiff1 d_ldp1p0_v0, const DiffSpace& space,
                                  const TimeSpace& t_space, double t_factor, double t_normal ) {

  double t2 = t_normal * t_normal;
  double t3 = t_normal * t2;

  get< 0 >( result ) = get_space< 0 >( space, t_space ).adjust(
    get< 0 >( a ), ( 3.0 * t2 - 2.0 * t3 ) * dp1p0
                   + ( t_normal - t2 * 2.0 + t3 ) * descend_to_space< 0 >( get< 1 >( a ), t_factor, space, t_space )
                   + ( t3 - t2 ) * descend_to_space< 0 >( get< 1 >( b ), t_factor, space, t_space ) );

  get< 1 >( result ) = get_space< 1 >( space, t_space ).adjust(
    get< 1 >( a ), ( ( t_normal - t2 ) * 6.0 ) * d_ldp1p0_v0 - ( 2.0 * t_normal - 3.0 * t2 ) * dv1v0 );

  cubic_hermite_interpolate_HOT_impl< Idx, PointType, PointDiff1, DiffSpace, TimeSpace >(
    result, dv1v0, d_ldp1p0_v0, space, t_space, t_factor, t_normal );
};

template < typename Idx, typename PointType, typename PointDiff0, typename PointDiff1, typename DiffSpace,
           typename TimeSpace >
inline typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 3 > >, void >::type
  cubic_hermite_interpolate_impl( PointType& result, const PointType& a, const PointType& b, const PointDiff0& dp1p0,
                                  const PointDiff1& dv1v0, const PointDiff1 d_ldp1p0_v0, const DiffSpace& space,
                                  const TimeSpace& t_space, double t_factor, double t_normal ) {
  cubic_hermite_interpolate_impl< typename boost::mpl::prior< Idx >::type, PointType, DiffSpace, TimeSpace >(
    result, a, b, dp1p0, dv1v0, d_ldp1p0_v0, space, t_space, t_factor, t_normal );

  get< Idx::type::value >( result ) = get_space< Idx::type::value >( space, t_space ).origin();
};
};


/**
 * This function template computes a cubic Hermite interpolation between two points in a
 * temporal and once-differentiable topology.
 * \tparam PointType The point type on the temporal and once-differentiable topology.
 * \tparam Topology The temporal and once-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template < typename PointType, typename Topology >
PointType cubic_hermite_interpolate( const PointType& a, const PointType& b, double t, const Topology& space ) {
  typedef typename temporal_space_traits< Topology >::space_topology SpaceType;
  typedef typename temporal_space_traits< Topology >::time_topology TimeSpaceType;
  BOOST_CONCEPT_ASSERT( (TemporalSpaceConcept< Topology >));
  BOOST_CONCEPT_ASSERT( (TangentBundleConcept< SpaceType, 1, TimeSpaceType >));

  double t_factor = b.time - a.time;
  if( std::fabs( t_factor ) < std::numeric_limits< double >::epsilon() )
    throw singularity_error( "Normalizing factor in cubic Hermite spline is zero!" );
  double t_normal = ( t - a.time ) / ( b.time - a.time );

  PointType result;
  result.time = t;

  typedef typename derived_N_order_space< SpaceType, TimeSpaceType, 0 >::type Space0;
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType, 1 >::type Space1;

  BOOST_CONCEPT_ASSERT( (LieGroupConcept< Space0 >));
  BOOST_CONCEPT_ASSERT( (LieGroupConcept< Space1 >));

  typedef typename topology_traits< Space0 >::point_difference_type PointDiff0;
  typedef typename topology_traits< Space1 >::point_difference_type PointDiff1;

  PointDiff0 dp1p0 = get_space< 0 >( space.get_space_topology(), space.get_time_topology() )
                       .difference( get< 0 >( b.pt ), get< 0 >( a.pt ) );
  PointDiff1 dv1v0 = get_space< 1 >( space.get_space_topology(), space.get_time_topology() )
                       .difference( get< 1 >( b.pt ), get< 1 >( a.pt ) );
  PointDiff1 d_ldp1p0_v0 = get_space< 1 >( space.get_space_topology(), space.get_time_topology() ).difference(
    lift_to_space< 1 >( dp1p0, t_factor, space.get_space_topology(), space.get_time_topology() ), get< 1 >( a.pt ) );

  detail::cubic_hermite_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(
    result.pt, a.pt, b.pt, dp1p0, dv1v0, d_ldp1p0_v0, space.get_space_topology(), space.get_time_topology(), t_factor,
    t_normal );

  return result;
};


/**
 * This functor class implements a cubic Hermite interpolation in a temporal and once-differentiable
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and
 * DifferentiableSpaceConcept once against time.
 * \tparam TimeSpaceType The time topology.
 */
template < typename SpaceType, typename TimeSpaceType >
class cubic_hermite_interpolator {
public:
  typedef cubic_hermite_interpolator< SpaceType, TimeSpaceType > self;
  typedef typename topology_traits< SpaceType >::point_type point_type;

  typedef typename derived_N_order_space< SpaceType, TimeSpaceType, 0 >::type Space0;
  typedef typename topology_traits< Space0 >::point_type PointType0;
  typedef typename topology_traits< Space0 >::point_difference_type PointDiff0;
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType, 1 >::type Space1;
  typedef typename topology_traits< Space1 >::point_type PointType1;
  typedef typename topology_traits< Space1 >::point_difference_type PointDiff1;

  BOOST_CONCEPT_ASSERT( ( TopologyConcept< SpaceType > ) );
  BOOST_CONCEPT_ASSERT( ( LieGroupConcept< Space0 > ) );
  BOOST_CONCEPT_ASSERT( ( LieGroupConcept< Space1 > ) );
  BOOST_CONCEPT_ASSERT( ( TangentBundleConcept< SpaceType, 1, TimeSpaceType > ) );

private:
  PointDiff0 delta_first_order;
  PointDiff1 delta_second_order;
  PointDiff1 delta_lifted_first_and_second;

public:
  /**
   * Default constructor.
   */
  cubic_hermite_interpolator(){};

  /**
   * Constructs the interpolator with its start and end points.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template < typename Factory >
  cubic_hermite_interpolator( const point_type& start_point, const point_type& end_point, double dt,
                              const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory ) {
    initialize( start_point, end_point, dt, space, t_space, factory );
  };

  /**
   * Initializes the interpolator with its start and end points.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param dt The time difference between the start point to the end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template < typename Factory >
  void initialize( const point_type& start_point, const point_type& end_point, double dt, const SpaceType& space,
                   const TimeSpaceType& t_space, const Factory& factory ) {
    delta_first_order = get_space< 0 >( space, t_space ).difference( get< 0 >( end_point ), get< 0 >( start_point ) );
    delta_second_order = get_space< 1 >( space, t_space ).difference( get< 1 >( end_point ), get< 1 >( start_point ) );
    delta_lifted_first_and_second = get_space< 1 >( space, t_space ).difference(
      lift_to_space< 1 >( delta_first_order, dt, space, t_space ), get< 1 >( start_point ) );
  };

  /**
   * Computes the point at a given delta-time from the start-point.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param result The result point of the interpolation.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param dt The time difference from the start-point to the resulting interpolated point.
   * \param dt_total The time difference from the start-point to the end point.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template < typename Factory >
  void compute_point( point_type& result, const point_type& start_point, const point_type& end_point,
                      const SpaceType& space, const TimeSpaceType& t_space, double dt, double dt_total,
                      const Factory& factory ) const {
    if( std::fabs( dt_total ) < std::numeric_limits< double >::epsilon() )
      throw singularity_error( "Normalizing factor in cubic Hermite spline is zero!" );
    double t_normal = dt / dt_total;

    detail::cubic_hermite_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(
      result, start_point, end_point, delta_first_order, delta_second_order, delta_lifted_first_and_second, space,
      t_space, dt_total, t_normal );
  };

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return 0.0; };
};


/**
 * This class is a factory class for cubic Hermite interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template < typename TemporalTopology >
class cubic_hermite_interp_factory : public serializable {
public:
  typedef cubic_hermite_interp_factory< TemporalTopology > self;
  typedef TemporalTopology topology;
  typedef typename topology_traits< TemporalTopology >::point_type point_type;
  typedef generic_interpolator< self, cubic_hermite_interpolator > interpolator_type;

  BOOST_CONCEPT_ASSERT( ( TemporalSpaceConcept< TemporalTopology > ) );

private:
  shared_ptr< topology > space;

public:
  cubic_hermite_interp_factory( const shared_ptr< topology >& aSpace = shared_ptr< topology >() ) : space( aSpace ){};

  void set_temporal_space( const shared_ptr< topology >& aSpace ) { space = aSpace; };
  const shared_ptr< topology >& get_temporal_space() const { return space; };

  interpolator_type create_interpolator( const point_type* pp1, const point_type* pp2 ) const {
    return interpolator_type( this, pp1, pp2 );
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const { A& RK_SERIAL_SAVE_WITH_NAME( space ); };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) { A& RK_SERIAL_LOAD_WITH_NAME( space ); };

  RK_RTTI_MAKE_ABSTRACT_1BASE( self, 0xC2430002, 1, "cubic_hermite_interp_factory", serializable )
};


template < typename SpaceType, typename TimeTopology >
struct get_tagged_spatial_interpolator< cubic_hermite_interpolation_tag, SpaceType, TimeTopology > {
  typedef detail::generic_interpolator_impl< cubic_hermite_interpolator, SpaceType, TimeTopology > type;
  typedef void* pseudo_factory_type;
};

template < typename TemporalSpaceType >
struct get_tagged_temporal_interpolator< cubic_hermite_interpolation_tag, TemporalSpaceType > {
  typedef generic_interpolator< cubic_hermite_interp_factory< TemporalSpaceType >, cubic_hermite_interpolator > type;
};


/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a cubic Hermite interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept
 * and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template < typename Topology, typename DistanceMetric = typename metric_space_traits< Topology >::distance_metric_type >
class cubic_hermite_interp_traj
  : public interpolated_trajectory< Topology, cubic_hermite_interp_factory< Topology >, DistanceMetric > {
public:
  BOOST_CONCEPT_ASSERT( ( TemporalSpaceConcept< Topology > ) );
  BOOST_CONCEPT_ASSERT( ( TangentBundleConcept< typename temporal_space_traits< Topology >::space_topology, 1,
                                                typename temporal_space_traits< Topology >::time_topology > ) );

  typedef cubic_hermite_interp_traj< Topology, DistanceMetric > self;
  typedef interpolated_trajectory< Topology, cubic_hermite_interp_factory< Topology >, DistanceMetric > base_class_type;

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
  explicit cubic_hermite_interp_traj( const shared_ptr< topology >& aSpace = shared_ptr< topology >( new topology() ),
                                      const distance_metric& aDist = distance_metric() )
      : base_class_type( aSpace, aDist, cubic_hermite_interp_factory< Topology >( aSpace ) ){};

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  cubic_hermite_interp_traj( const shared_ptr< topology >& aSpace, const point_type& aStart, const point_type& aEnd,
                             const distance_metric& aDist = distance_metric() )
      : base_class_type( aSpace, aStart, aEnd, aDist, cubic_hermite_interp_factory< Topology >( aSpace ) ){};

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template < typename ForwardIter >
  cubic_hermite_interp_traj( ForwardIter aBegin, ForwardIter aEnd, const shared_ptr< topology >& aSpace,
                             const distance_metric& aDist = distance_metric() )
      : base_class_type( aBegin, aEnd, aSpace, aDist, cubic_hermite_interp_factory< Topology >( aSpace ) ){};


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    base_class_type::save( A, base_class_type::getStaticObjectType()->TypeVersion() );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    base_class_type::load( A, base_class_type::getStaticObjectType()->TypeVersion() );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2440004, 1, "cubic_hermite_interp_traj", base_class_type )
};
};


namespace rtti {

template <>
struct get_type_id< pp::cubic_hermite_interpolation_tag > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 2 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "cubic_hermite_interpolation_tag" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "cubic_hermite_interpolation_tag"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};
};
};

#endif
