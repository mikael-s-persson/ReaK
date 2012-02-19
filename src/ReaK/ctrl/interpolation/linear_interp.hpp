/**
 * \file linear_interp.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation.
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

#ifndef REAK_LINEAR_INTERP_HPP
#define REAK_LINEAR_INTERP_HPP

#include "path_planning/spatial_trajectory_concept.hpp"

#include "path_planning/tangent_bundle_concept.hpp"

#include "interpolated_trajectory.hpp"
#include "generic_interpolator_factory.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include "lin_alg/mat_num_exceptions.hpp"
#include "topologies/basic_distance_metrics.hpp"

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  template <typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type linear_interpolate_HOT_impl(PointType& result, const PointDiff1& dv1v0, const PointDiff1& d_ldp1p0_v0,
                                                  const DiffSpace& space, const TimeSpace& t_space,
				 	          double t_factor, double t_normal) {
    /* nothing to do. */
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type linear_interpolate_HOT_impl(PointType& result, const PointDiff0& dp1p0,
                                           const DiffSpace& space, const TimeSpace& t_space,
				 	   double t_factor, double t_normal) {
    get<1>(result) = lift_to_space<1>(dp1p0, t_factor, space, t_space);
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type linear_interpolate_impl(PointType& result, const PointType& a, const PointDiff0& dp1p0,
                                       const DiffSpace& space, const TimeSpace& t_space,
				       double t_factor, double t_normal) {
    
    get<0>(result) = get_space<0>(space,t_space).adjust(get<0>(a), t_normal * dp1p0);
    
    linear_interpolate_HOT_impl< Idx, PointType, PointDiff0, DiffSpace, TimeSpace >(result, dp1p0, space, t_space, t_factor, t_normal);
    
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type linear_interpolate_impl(PointType& result, const PointType& a, const PointDiff0 dp1p0, 
				       const DiffSpace& space, const TimeSpace& t_space,
				       double t_factor, double t_normal) {
    linear_interpolate_impl< typename boost::mpl::prior<Idx>::type, PointType, PointDiff0, DiffSpace, TimeSpace >(result,a,dp1p0,space,t_space,t_factor,t_normal);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
};



/**
 * This function template computes a linear interpolation between two points in a 
 * temporal and zero-differentiable topology.
 * \tparam PointType The point type on the temporal and zero-differentiable topology.
 * \tparam Topology The temporal and zero-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType linear_interpolate(const PointType& a, const PointType& b, double t, const Topology& space) {
  typedef typename temporal_space_traits<Topology>::space_topology SpaceType;
  typedef typename temporal_space_traits<Topology>::time_topology TimeSpaceType;
  
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept< SpaceType, 0, TimeSpaceType>));
  
  double t_factor = b.time - a.time;
  if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
    throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
  double t_normal = (t - a.time) / (b.time - a.time);
      
  PointType result;
  result.time = t;
  
  typedef typename derived_N_order_space< SpaceType ,TimeSpaceType,0>::type Space0;
  
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  
  typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
    
  PointDiff0 dp1p0 = get_space<0>(space.get_space_topology(),space.get_time_topology()).difference( get<0>(b.pt), get<0>(a.pt) );
  
  detail::linear_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(result.pt, a.pt, dp1p0, space.get_space_topology(), space.get_time_topology(), t_factor, t_normal);
  
  return result;
};



/**
 * This functor class implements a linear interpolation in a temporal and zero-differentiable 
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and DifferentiableSpaceConcept zero times against time.
 * \tparam TimeSpaceType The time topology.
 */

template <typename SpaceType, typename TimeSpaceType>
class linear_interpolator {
  public:
    typedef linear_interpolator<SpaceType,TimeSpaceType> self;
    typedef typename topology_traits<SpaceType>::point_type point_type;
    
    typedef typename derived_N_order_space< SpaceType,TimeSpaceType,0>::type Space0;
    typedef typename topology_traits<Space0>::point_type PointType0;
    typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
  
    BOOST_CONCEPT_ASSERT((TopologyConcept<SpaceType>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept< SpaceType, 0, TimeSpaceType >));
    
  private:
    PointDiff0 delta_first_order;
    
  public:
    
    /**
     * Default constructor.
     */
    linear_interpolator() { };
    
    /**
     * Constructs the interpolator with its start and end points.
     * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
     * \param start_point The start point of the interpolation.
     * \param end_point The end point of the interpolation.
     * \param space The metric space on which the interpolation resides.
     * \param t_space The time-space against which the interpolation is done.
     * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
     */
    template <typename Factory>
    linear_interpolator(const point_type& start_point, const point_type& end_point, double dt,
		        const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory) {
      initialize(start_point,end_point,dt,space,t_space,factory);
    };
    
    /**
     * Initializes the interpolator with its start and end points.
     * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
     * \param start_point The start point of the interpolation.
     * \param end_point The end point of the interpolation.
     * \param space The metric space on which the interpolation resides.
     * \param t_space The time-space against which the interpolation is done.
     * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
     */
    template <typename Factory>
    void initialize(const point_type& start_point, const point_type& end_point, double,
		    const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory) {
      delta_first_order = get_space<0>(space,t_space).difference( get<0>(end_point), get<0>(start_point) );
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
    template <typename Factory>
    void compute_point(point_type& result, const point_type& start_point, const point_type& end_point,
		       const SpaceType& space, const TimeSpaceType& t_space, 
		       double dt, double dt_total, const Factory& factory) const {
      if(std::fabs(dt_total) < std::numeric_limits<double>::epsilon())
        throw singularity_error("Normalizing factor in linear interpolation is zero!");
      double t_normal = dt / dt_total;
      
      detail::linear_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(result, start_point, delta_first_order, space, t_space, dt_total, t_normal);
    };
    
    /**
     * Returns the minimum travel time between the initialized start and end points.
     * \return The minimum travel time between the initialized start and end points.
     */
    double get_minimum_travel_time() const {
      return 0.0;
    };
    
};


/**
 * This class is a factory class for linear interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template <typename TemporalTopology>
class linear_interpolator_factory : public serialization::serializable {
  public:
    typedef linear_interpolator_factory<TemporalTopology> self;
    typedef TemporalTopology topology;
    typedef typename topology_traits<TemporalTopology>::point_type point_type;
    typedef generic_interpolator<self,linear_interpolator> interpolator_type;
  
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
  private:
    typename shared_pointer<topology>::type space;
  public:
    linear_interpolator_factory(const typename shared_pointer<topology>::type& aSpace = typename shared_pointer<topology>::type()) : space(aSpace) { };
  
    void set_temporal_space(const typename shared_pointer<topology>::type& aSpace) { space = aSpace; };
    const typename shared_pointer<topology>::type& get_temporal_space() const { return space; };
  
    interpolator_type create_interpolator(const point_type* pp1, const point_type* pp2) const {
      return interpolator_type(this, pp1, pp2);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      A & RK_SERIAL_LOAD_WITH_NAME(space);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2430001,1,"linear_interpolator_factory",serialization::serializable)
};




  
/**
 * This class implements a trajectory in a temporal and zero-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class linear_interp_traj : public interpolated_trajectory<Topology,linear_interpolator_factory<Topology>,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept< typename temporal_space_traits<Topology>::space_topology, 1, typename temporal_space_traits<Topology>::time_topology >));
    
    typedef linear_interp_traj<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,linear_interpolator_factory<Topology>,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric_type distance_metric_type;
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit linear_interp_traj(const typename shared_pointer<topology>::type& aSpace = typename shared_pointer<topology>::type(new topology()), const distance_metric_type& aDist = distance_metric_type()) : 
                                base_class_type(aSpace, aDist, linear_interpolator_factory<Topology>(aSpace)) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    linear_interp_traj(const typename shared_pointer<topology>::type& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric_type& aDist = distance_metric_type()) :
                       base_class_type(aSpace, aStart, aEnd, aDist, linear_interpolator_factory<Topology>(aSpace)) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    linear_interp_traj(ForwardIter aBegin, ForwardIter aEnd, const typename shared_pointer<topology>::type& aSpace, const distance_metric_type& aDist = distance_metric_type()) : 
                       base_class_type(aBegin, aEnd, aSpace, aDist, linear_interpolator_factory<Topology>(aSpace)) { };
		       
		       
		       
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440003,1,"linear_interp_traj",base_class_type)
    
};



};

};

#endif









