/**
 * \file sustained_velocity_pulse.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points 
 * are computed with a rate-limited sustained velocity pulse (SVP) interpolation.
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

#ifndef REAK_SUSTAINED_VELOCITY_PULSE_HPP
#define REAK_SUSTAINED_VELOCITY_PULSE_HPP

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
#include "topologies/basic_distance_metrics.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "topologies/generic_sampler_factory.hpp"

#include "sustained_velocity_pulse_detail.hpp"

namespace ReaK {

namespace pp {


/**
 * This function template computes a Sustained Velocity Pulse (SVP) interpolation between two points in a 
 * temporal and once-differentiable topology.
 * \tparam PointType The point type on the temporal and once-differentiable topology.
 * \tparam Topology The temporal and once-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType svp_interpolate(const PointType& a, const PointType& b, double t, const Topology& space, 
			  double tolerance = 1e-6, unsigned int maximum_iterations = 60) {
  typedef typename temporal_space_traits<Topology>::space_topology SpaceType;
  typedef typename temporal_space_traits<Topology>::time_topology TimeSpaceType;
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept< SpaceType, 1, TimeSpaceType>));
  BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 1>::type >));
  
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType,0>::type Space0;
  typedef typename topology_traits<Space0>::point_type PointType0;
  typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
  
  typedef typename derived_N_order_space< SpaceType, TimeSpaceType,1>::type Space1;
  typedef typename topology_traits<Space1>::point_type PointType1;
  typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
  
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space0>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  
  if(t <= a.time)
    return a;
  if(t >= b.time)
    return b;
  
  PointDiff0 delta_first_order;
  PointType1 peak_velocity;
  
  double delta_time = b.time - a.time;
  
  double slack = detail::svp_compute_interpolation_data_impl(a.pt, b.pt,
                                                             delta_first_order, peak_velocity,
					                     space.get_space_topology(),
						             space.get_time_topology(),
						             delta_time, NULL,
							     tolerance, maximum_iterations);
  
  if(slack < 0.0)
    delta_time -= slack;
  double dt = t - a.time;
      
  PointType result;
  result.time = t;
      
  detail::svp_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(result.pt, a.pt, b.pt, delta_first_order, peak_velocity, space.get_space_topology(), space.get_time_topology(), dt, delta_time);
  
  return result;    
};







/**
 * This functor class implements a sustained velocity pulse (SVP) interpolation in a temporal and once-differentiable 
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and DifferentiableSpaceConcept once against time.
 * \tparam TimeSpaceType The time topology.
 */
template <typename SpaceType, typename TimeSpaceType>
class svp_interpolator {
  public:
    typedef svp_interpolator<SpaceType,TimeSpaceType> self;
    typedef typename topology_traits<SpaceType>::point_type point_type;
  
    typedef typename derived_N_order_space< SpaceType,TimeSpaceType,0>::type Space0;
    typedef typename topology_traits<Space0>::point_type PointType0;
    typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename derived_N_order_space< SpaceType,TimeSpaceType,1>::type Space1;
    typedef typename topology_traits<Space1>::point_type PointType1;
    typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
    
    BOOST_CONCEPT_ASSERT((TopologyConcept<SpaceType>));
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space0>));
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space1>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept< SpaceType, 1, TimeSpaceType >));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< typename derived_N_order_space<SpaceType, TimeSpaceType, 1>::type >));
    
  private:
    PointDiff0 delta_first_order;
    PointType1 peak_velocity;
    double min_delta_time;
    PointType1 best_peak_velocity;
  
  public:
    
    
    /**
     * Default constructor.
     */
    svp_interpolator() : min_delta_time(std::numeric_limits<double>::infinity()) { };
    
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
    svp_interpolator(const point_type& start_point, const point_type& end_point, double dt,
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
      
      min_delta_time = detail::svp_compute_interpolation_data_impl(start_point, end_point,
	                                                           delta_first_order, peak_velocity,
						                   space, t_space, dt, &best_peak_velocity,
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
      
      detail::svp_interpolate_impl< max_derivation_order< SpaceType, TimeSpaceType > >(result, start_point, end_point, delta_first_order, peak_velocity, space, t_space, dt, dt_total);
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
 * This class is a factory class for sustained velocity pulse (SVP) interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept, 
 *                          with a spatial topology that is once-differentiable (see DifferentiableSpaceConcept) and 
 *                          whose 1-order derivative space has a spherical bound (see SphereBoundedSpaceConcept).
 */
template <typename TemporalTopology>
class svp_interpolator_factory : public serialization::serializable {
  public:
    typedef svp_interpolator_factory<TemporalTopology> self;
    typedef TemporalTopology topology;
    typedef typename topology_traits<TemporalTopology>::point_type point_type;
  
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
    
    typedef generic_interpolator<self, svp_interpolator> interpolator_type;
    
  private:
    shared_ptr<topology> space;
  public:
    double tolerance;
    unsigned int maximum_iterations;
    
    svp_interpolator_factory(const shared_ptr<topology>& aSpace = shared_ptr<topology>(), 
			     double aTolerance = 1e-6, 
			     unsigned int aMaxIter = 60) : 
			     space(aSpace),
			     tolerance(aTolerance),
			     maximum_iterations(aMaxIter) { };
  
    void set_temporal_space(const shared_ptr<topology>& aSpace) { space = aSpace; };
    const shared_ptr<topology>& get_temporal_space() const { return space; };
  
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

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2430004,1,"svp_interpolator_factory",serialization::serializable)
};


/**
 * This functor class is a distance metric based on the reach-time of a SVP interpolation between
 * two points in a differentiable space.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType>
struct svp_reach_time_metric : public serialization::serializable {
  
  typedef svp_reach_time_metric<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  double tolerance;
  unsigned int maximum_iterations;
  
  svp_reach_time_metric(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType()),
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
    detail::generic_interpolator_impl<svp_interpolator,Topology,TimeSpaceType> interp;
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
    detail::generic_interpolator_impl<svp_interpolator,Topology,TimeSpaceType> interp;
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(svp_reach_time_metric,0xC2410009,1,"svp_reach_time_metric",serialization::serializable)
};



/**
 * This functor class is a random-sampler based on the rate-limited motions of a SVP interpolation 
 * between points within a bounded tangent-bundle.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType>
struct svp_rate_limited_sampler : public serialization::serializable {
  
  typedef svp_rate_limited_sampler<TimeSpaceType> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  svp_rate_limited_sampler(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
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
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 1, TimeSpaceType>));
    
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 0>::type Space0;
    typedef typename derived_N_order_space<Topology, TimeSpaceType, 1>::type Space1;
    typedef typename topology_traits<Topology>::point_type PointType;
    typedef typename topology_traits<Space0>::point_type Point0;
    typedef typename topology_traits<Space1>::point_type Point1;
    typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
    typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
    
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
    BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< Space1 >));
    
    const typename point_distribution_traits<Topology>::random_sampler_type& get_sample = get(random_sampler, s);
    const Space0& s0 = get_space<0>(s, *t_space);
    const Space1& s1 = get_space<1>(s, *t_space);
    const typename metric_space_traits< Space1 >::distance_metric_type& get_vel_dist = get(distance_metric, s1);
    
    while(true) {
      PointType pt = get_sample(s);
      
      PointDiff1 dp1 = s1.difference(s1.origin(), get<1>(pt));
      double dt = get_vel_dist(s1.origin(), get<1>(pt), s1);
      
      // Check if we can stop the motion before the boundary.
      PointType result = pt;
      detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result, dp1, s, *t_space, dt);
      if( !s.is_in_bounds(result) )
	continue; //reject the sample.
      
      // Check if we could have initiated the motion from within the boundary.
      result = pt;
      detail::svp_constant_accel_motion_impl< max_derivation_order< Topology, TimeSpaceType > >(result, -dp1, s, *t_space, -dt);
      if( !s.is_in_bounds(result) )
	continue; //reject the sample.
      
      // If this point is reached, it means that the sample is acceptable:
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(svp_rate_limited_sampler,0xC2450001,1,"svp_rate_limited_sampler",serialization::serializable)
};





/**
 * This class wraps a reach-time topology with SVP-based distance metric and a sampler.
 * \tparam BaseTopology The topology underlying this space, should express values as reach-time values and metrics (distance), and should model TopologyConcept, PointDistributionConcept, BoundedSpaceConcept and TangentBundleConcept for time_topology and up to 2nd order (acceleration).
 */
template <typename BaseTopology>
class svp_reach_topology : public BaseTopology
{
  public:
    BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((PointDistributionConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<BaseTopology>));
    BOOST_CONCEPT_ASSERT((TangentBundleConcept<BaseTopology, 2, time_topology>));
    
    typedef svp_reach_topology<BaseTopology> self;
    
    typedef typename topology_traits< BaseTopology >::point_type point_type;
    typedef typename topology_traits< BaseTopology >::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef BaseTopology super_space_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< BaseTopology >::dimensions);
    
  protected:
    
    shared_ptr<time_topology> t_space;
    svp_reach_time_metric<time_topology> rt_dist;
    generic_sampler<svp_rate_limited_sampler<time_topology>, BaseTopology> rl_sampler;
    
    
  public:
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    svp_reach_topology(Args&&... args) : 
                       BaseTopology(std::forward<Args>(args)...),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
#else
    svp_reach_topology() : 
                       BaseTopology(),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1>
    svp_reach_topology(const A1& a1) : 
                       BaseTopology(a1),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2>
    svp_reach_topology(const A1& a1, const A2& a2) : 
                       BaseTopology(a1, a2),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1, typename A2, typename A3>
    svp_reach_topology(const A1& a1, const A2& a2, const A3& a3) : 
                       BaseTopology(a1, a2, a3),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1>
    svp_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                       BaseTopology(a1, a2, a3, a4),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
    
    template <typename A1>
    svp_reach_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                       BaseTopology(a1, a2, a3, a4, a5),
                       t_space(new time_topology),
                       rt_dist(t_space), rl_sampler(svp_rate_limited_sampler<time_topology>(t_space)) { };
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
    }
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      return rt_dist(delta, static_cast<const BaseTopology&>(*this));
    }
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed within the reachable space.
     */
    point_type random_point() const {
      return rl_sampler(static_cast<const BaseTopology&>(*this));
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400023,1,"svp_reach_topology",BaseTopology)
    
};






  
/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a sustained velocity pulse interpolation (limited). This class models 
 * the SpatialTrajectoryConcept.
 * \tparam Topology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept, 
 *                  with a spatial topology that is once-differentiable (see DifferentiableSpaceConcept) and 
 *                  whose 1-order derivative space has a spherical bound (see SphereBoundedSpaceConcept).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type>
class svp_interp_traj : public interpolated_trajectory<Topology,svp_interpolator_factory<Topology>,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef svp_interp_traj<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,svp_interpolator_factory<Topology>,DistanceMetric> base_class_type;
    
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
    explicit svp_interp_traj(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), const distance_metric& aDist = distance_metric()) : 
                             base_class_type(aSpace, aDist, svp_interpolator_factory<Topology>(aSpace)) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    svp_interp_traj(const shared_ptr<topology>& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                    base_class_type(aSpace, aStart, aEnd, aDist, svp_interpolator_factory<Topology>(aSpace)) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    svp_interp_traj(ForwardIter aBegin, ForwardIter aEnd, const shared_ptr<topology>& aSpace, const distance_metric& aDist = distance_metric()) : 
                    base_class_type(aBegin, aEnd, aSpace, aDist, svp_interpolator_factory<Topology>(aSpace)) { };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440006,1,"svp_interp_traj",base_class_type)
};



};

};

#endif









