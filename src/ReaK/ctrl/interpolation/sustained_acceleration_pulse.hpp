/**
 * \file sustained_acceleration_pulse.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a rate-limited sustained acceleration pulse (SAP) interpolation.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_SUSTAINED_ACCELERATION_PULSE_HPP
#define REAK_SUSTAINED_ACCELERATION_PULSE_HPP

#include "base/defs.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"

#include "path_planning/differentiable_space_concept.hpp"

#include "interpolated_trajectory.hpp"
#include "generic_interpolator_factory.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include "topologies/basic_distance_metrics.hpp"
#include "path_planning/bounded_space_concept.hpp"

#include "sustained_acceleration_pulse_detail.hpp"

namespace ReaK {

namespace pp {


/**
 * This function template computes a Sustained Acceleration Pulse (SAP) interpolation between two points in a 
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
PointType sap_interpolate(const PointType& a, const PointType& b, double t, const Topology& space) {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  typedef typename temporal_topology_traits<Topology>::space_topology SpaceType;
  typedef typename temporal_topology_traits<Topology>::time_topology TimeSpaceType;
  BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< SpaceType, 2, TimeSpaceType>));
  BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 1>::type >));
  BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 2>::type >));
  
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType,0>::type Space0;
  typedef typename metric_topology_traits<Space0>::point_type PointType0;
  typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
  
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType,1>::type Space1;
  typedef typename metric_topology_traits<Space1>::point_type PointType1;
  typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff1;

  typedef typename derived_N_order_space< SpaceType, TimeSpaceType,2>::type Space2;
  typedef typename metric_topology_traits<Space2>::point_type PointType2;
  typedef typename metric_topology_traits<Space2>::point_difference_type PointDiff2;
  
  if(t <= a.time)
    return a;
  if(t >= b.time)
    return b;
  
  PointDiff0 delta_first_order;
  PointType1 peak_velocity;  
  double delta_time = b.time - a.time;
  
  double min_delta_time = detail::sap_compute_interpolation_data_impl(a.pt, b.pt,
                                                                  delta_first_order, peak_velocity,
								  space.get_space_topology(),
								  space.get_time_topology(),
								  delta_time, NULL,
								  1e-6, 60);
  
  if(min_delta_time > delta_time)
    delta_time = min_delta_time;
  double dt = t - a.time;
      
  PointType result;
  result.time = t;
      
  detail::sap_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result.pt, a.pt, b.pt, delta_first_order, peak_velocity, space.get_space_topology(), space.get_time_topology(), dt, delta_time);
  
  return result;    
};





/**
 * This functor class implements a sustained acceleration pulse (SAP) interpolation in a temporal 
 * and twice-differentiable topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and DifferentiableSpaceConcept once against time.
 * \tparam TimeSpaceType The time topology.
 */
template <typename SpaceType, typename TimeSpaceType>
class sap_interpolator {
  public:
    typedef sap_interpolator<SpaceType,TimeSpaceType> self;
    typedef typename metric_topology_traits<SpaceType>::point_type point_type;
  
    typedef typename derived_N_order_space< SpaceType,TimeSpaceType,0>::type Space0;
    typedef typename metric_topology_traits<Space0>::point_type PointType0;
    typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename derived_N_order_space< SpaceType,TimeSpaceType,1>::type Space1;
    typedef typename metric_topology_traits<Space1>::point_type PointType1;
    typedef typename metric_topology_traits<Space1>::point_difference_type PointDiff1;
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<SpaceType>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< SpaceType, 2, TimeSpaceType >));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 1>::type >));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 2>::type >));
  
  private:
    PointDiff0 delta_first_order;
    PointType1 peak_velocity;
    double min_delta_time;
    PointType1 best_peak_velocity;
  
  public:
    
    
    /**
     * Default constructor.
     */
    sap_interpolator() : min_delta_time(std::numeric_limits<double>::infinity()) { };
    
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
    sap_interpolator(const point_type& start_point, const point_type& end_point, double dt,
		     const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory) {
      initialize(start_point,end_point,dt,space,t_space,factory);
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
    template <typename Factory>
    void initialize(const point_type& start_point, const point_type& end_point, double dt,
		    const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory) {
      
      min_delta_time = detail::sap_compute_interpolation_data_impl(start_point, end_point,
                                                                   delta_first_order, peak_velocity,
								   space, t_space,
								   dt, &best_peak_velocity,
								   factory.tolerance, factory.maximum_iterations);
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
      if(dt <= 0.0) {
	result = start_point;
	return;
      };
      if(dt >= dt_total) {
	result = end_point;
	return;
      };
      
      detail::sap_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result, start_point, end_point, delta_first_order, peak_velocity, space, t_space, dt, dt_total);
    };
    
    /**
     * Returns the minimum travel time between the initialized start and end points.
     * \return The minimum travel time between the initialized start and end points.
     */
    double get_minimum_travel_time() const {
      return min_delta_time;
    };
    
};



/**
 * This class is a factory class for sustained acceleration pulse (SAP) interpolators on a temporal 
 * differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept, 
 *                          with a spatial topology that is twice-differentiable (see DifferentiableSpaceConcept) and 
 *                          whose 1-order and 2-order derivative space has a spherical bound (see SphereBoundedSpaceConcept).
 */
template <typename TemporalTopology>
class sap_interpolator_factory : public serialization::serializable {
  public:
    typedef sap_interpolator_factory<TemporalTopology> self;
    typedef TemporalTopology topology;
    typedef typename temporal_topology_traits<TemporalTopology>::point_type point_type;
  
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
    
    typedef generic_interpolator<self,sap_interpolator> interpolator_type;
    
  private:
    typename shared_pointer<topology>::type space;
  public:
    double tolerance;
    unsigned int maximum_iterations;
    
    sap_interpolator_factory(const typename shared_pointer<topology>::type& aSpace = typename shared_pointer<topology>::type(), 
			     double aTolerance = 1e-6, 
			     unsigned int aMaxIter = 60) : 
			     space(aSpace),
			     tolerance(aTolerance),
			     maximum_iterations(aMaxIter)  { };
  
    void set_temporal_space(const typename shared_pointer<topology>::type& aSpace) { space = aSpace; };
    const typename shared_pointer<topology>::type& get_temporal_space() const { return space; };
  
    interpolator_type create_interpolator(const point_type* pp1, const point_type* pp2) const {
      return interpolator_type(this, pp1, pp2);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { 
      A & RK_SERIAL_SAVE_WITH_NAME(space)
        & RK_SERIAL_SAVE_WITH_NAME(tolerance)
        & RK_SERIAL_SAVE_WITH_NAME(maximum_iterations);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      A & RK_SERIAL_LOAD_WITH_NAME(space)
        & RK_SERIAL_LOAD_WITH_NAME(tolerance)
        & RK_SERIAL_LOAD_WITH_NAME(maximum_iterations);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2430005,1,"sap_interpolator_factory",serialization::serializable)
};



/**
 * This functor class is a distance metric based on the reach-time of a SAP interpolation between
 * two points in a differentiable space.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType>
struct sap_reach_time_metric : public serialization::serializable {
  
  typedef sap_reach_time_metric<TimeSpaceType> self;
  
  typename shared_pointer<TimeSpaceType>::type t_space;
  double tolerance;
  unsigned int maximum_iterations;
  
  sap_reach_time_metric(const typename shared_pointer<TimeSpaceType>::type& aTimeSpace = typename shared_pointer<TimeSpaceType>::type(new TimeSpaceType()),
                        double aTolerance = 1e-6, 
			unsigned int aMaxIter = 60) : 
			t_space(aTimeSpace),
			tolerance(aTolerance),
			maximum_iterations(aMaxIter) { };
  
  /** 
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    detail::generic_interpolator_impl<sap_interpolator,Topology,TimeSpaceType> interp;
    interp.initialize(a, b, 0.0, s, *t_space, *this);
    return interp.get_minimum_travel_time();
  };
  
  /** 
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    detail::generic_interpolator_impl<sap_interpolator,Topology,TimeSpaceType> interp;
    interp.initialize(s.origin(), s.adjust(s.origin(),a), 0.0, s, *t_space, *this);
    return interp.get_minimum_travel_time();
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(t_space)
      & RK_SERIAL_SAVE_WITH_NAME(tolerance)
      & RK_SERIAL_SAVE_WITH_NAME(maximum_iterations);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(t_space)
      & RK_SERIAL_LOAD_WITH_NAME(tolerance)
      & RK_SERIAL_LOAD_WITH_NAME(maximum_iterations);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(sap_reach_time_metric,0xC241000A,1,"sap_reach_time_metric",serialization::serializable)
};


  
/**
 * This class implements a trajectory in a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a sustained acceleration pulse interpolation (limited). This class models 
 * the SpatialTrajectoryConcept.
 * \tparam Topology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept, 
 *                  with a spatial topology that is twice-differentiable (see DifferentiableSpaceConcept) and 
 *                  whose 1-order and 2-order derivative spaces have a spherical bound (see SphereBoundedSpaceConcept).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class sap_interp_traj : public interpolated_trajectory<Topology,sap_interpolator_factory<Topology>,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef sap_interp_traj<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,sap_interpolator_factory<Topology>,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit sap_interp_traj(const typename shared_pointer<topology>::type& aSpace = typename shared_pointer<topology>::type(new topology()), const distance_metric& aDist = distance_metric()) : 
                             base_class_type(aSpace, aDist, sap_interpolator_factory<Topology>(aSpace)) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    sap_interp_traj(const typename shared_pointer<topology>::type& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                    base_class_type(aSpace, aStart, aEnd, aDist, sap_interpolator_factory<Topology>(aSpace)) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    sap_interp_traj(ForwardIter aBegin, ForwardIter aEnd, const typename shared_pointer<topology>::type& aSpace, const distance_metric& aDist = distance_metric()) : 
                    base_class_type(aBegin, aEnd, aSpace, aDist, sap_interpolator_factory<Topology>(aSpace)) { };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440007,1,"sap_interp_traj",base_class_type)
};



};

};

#endif









