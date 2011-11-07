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

#include "path_planning/differentiable_space_concept.hpp"

#include "interpolated_trajectory.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include "lin_alg/mat_num_exceptions.hpp"
#include "topologies/basic_distance_metrics.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "root_finders/secant_method.hpp"
#include "root_finders/newton_raphson_method.hpp"
#include <root_finders/bisection_method.hpp>
#include <boost/bind.hpp>

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
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
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
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    
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
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    linear_interpolate_impl< typename boost::mpl::prior<Idx>::type, PointType, PointDiff0, DiffSpace, TimeSpace >(result,a,dp1p0,space,t_space,t_factor,t_normal);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
  inline
  double svp_compute_slack_time(double beta, double dt, double norm_delta, const vect<double,5>& coefs) {
    using std::sqrt;
    return beta * ( dt - sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4])) 
                       - sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4])) )
           - norm_delta;
  };
    
  inline
  double svp_compute_derivative_slack_time(double beta, double dt, double norm_delta, const vect<double,5>& coefs) {
    using std::sqrt;
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    return dt - term1 - term2 + 2.0 * beta * coefs[4] * ((coefs[0] * coefs[1] - coefs[4] * beta) / term1 + (coefs[2] * coefs[3] - coefs[4] * beta) / term2);
  };
  
  
  inline
  double svp_compute_travel_time(double beta, double norm_delta, const vect<double,5>& coefs) {
    using std::sqrt;
    return sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4])) 
           + sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4])) 
           + norm_delta / beta;
  };
    
  inline
  double svp_compute_derivative_travel_time(double beta, double norm_delta, const vect<double,5>& coefs) {
    using std::sqrt;
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    return 2.0 * coefs[4] * ((coefs[4] * beta - coefs[0] * coefs[1]) / term1 + (coefs[4] * beta - coefs[2] * coefs[3]) / term2) 
           - norm_delta / (beta * beta);
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
PointType svp_interpolate(const PointType& a, const PointType& b, double t, const Topology& space) {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  typedef typename temporal_topology_traits<Topology>::space_topology SpaceType;
  typedef typename temporal_topology_traits<Topology>::time_topology TimeSpaceType;
  BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< SpaceType, 1, TimeSpaceType>));
  BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< derived_N_order_space<SpaceType, TimeSpaceType, 1> >));
  
  double t_factor = b.time - a.time;
  if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
    throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
  double t_normal = (t - a.time) / (b.time - a.time);
      
  PointType result;
  result.time = t;
  
  typedef typename derived_N_order_space< SpaceType ,TimeSpaceType,0>::type Space0;
  typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
    
  PointDiff0 dp1p0 = get_space<0>(space.get_space_topology(),space.get_time_topology()).difference( get<0>(b.pt), get<0>(a.pt) );
  
  detail::linear_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result.pt, a.pt, dp1p0, space.get_space_topology(), space.get_time_topology(), t_factor, t_normal);
  
  return result;      
};

/**
 * This functor class implements a sustained velocity pulse (SVP) interpolation in a temporal and once-differentiable 
 * topology.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done.
 */
template <typename Factory, 
          std::size_t Dimensions >
class svp_interpolator {
  public:
    typedef svp_interpolator<Factory> self;
    typedef typename Factory::point_type point_type;
    typedef typename Factory::topology topology;
  
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,0>::type Space0;
    typedef typename metric_topology_traits<Space0>::point_type PointType0;
    typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff0;
  
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,1>::type Space1;
    typedef typename metric_topology_traits<Space0>::point_type PointType1;
    typedef typename metric_topology_traits<Space0>::point_difference_type PointDiff1;
  
  private:
    const Factory* parent;
    const point_type* start_point;
    const point_type* end_point;
    PointDiff0 delta_first_order;
    PointType1 peak_velocity;
    double min_delta_time;
    PointType1 best_peak_velocity;
    
    void update_delta_value() {
      if(parent && start_point && end_point) {
	const Space0& space_0 = get_space<0>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology());
        const Space1& space_1 = get_space<1>(parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology());
        delta_first_order = space_0.difference( get<0>(end_point->pt), get<0>(start_point->pt) );
	double norm_delta = space_0.norm( delta_first_order );
	double beta = 1.0;
	peak_velocity = lift_to_space<1>(delta_first_order, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	vect<double,5> coefs( space_1.distance( get<1>(start_point->pt), space_1.origin() ),
	                      0.0,
		              space_1.distance( get<1>(end_point->pt), space_1.origin() ),
	                      0.0,
		              space_1.get_radius());
	for(unsigned int i = 0; i < 20; ++i) {
	  PointType1 prev_peak_velocity = peak_velocity;
	  
	  peak_velocity = lift_to_space<1>(delta_first_order, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  double temp = space_1.distance(get<1>(start_point->pt), peak_velocity);
	  coefs[1] = (coefs[0] * coefs[0] + coefs[4] * coefs[4] - temp * temp) / (2.0 * coefs[0] * coefs[4]);
	  temp = space_1.distance(get<1>(end_point->pt), peak_velocity);
	  coefs[3] = (coefs[2] * coefs[2] + coefs[4] * coefs[4] - temp * temp) / (2.0 * coefs[2] * coefs[4]);
	  
	  if(detail::svp_compute_derivative_travel_time(1.0,norm_delta,coefs) > 0.0) {
	    double upper = 1.0;
	    double lower = 0.5;
	    while(detail::svp_compute_derivative_travel_time(lower,norm_delta,coefs) > 0.0) {
	      upper = lower;
	      lower *= 0.5;
	    };
	    bisection_method(lower, upper, boost::bind(detail::svp_compute_derivative_travel_time,_1,boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	    beta = upper;
	  } else {
	    beta = 1.0;
	  };
	  
	  min_delta_time = detail::svp_compute_travel_time(beta,norm_delta,coefs);
	  peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  
	  PointDiff0 descended_peak_velocity = descend_to_space<0>(peak_velocity,1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology());
	  delta_first_order = space_0.difference(
	    space_0.adjust(get<0>(end_point->pt), (-0.5 * space_1.distance(peak_velocity, get<1>(end_point->pt))) * ( descended_peak_velocity
	                                                                                                             + descend_to_space<0>(get<1>(end_point->pt),1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()) )),
	    space_0.adjust(get<0>(start_point->pt), (0.5 * space_1.distance(peak_velocity, get<1>(start_point->pt))) * (descended_peak_velocity
	                                                                                                                + descend_to_space<0>(get<1>(start_point->pt),1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()) ))
	  );
	  norm_delta = space_0.norm(delta_first_order);
	  peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  
	  if( space_1.distance(prev_peak_velocity, peak_velocity) < 1e-6 * space_1.get_radius() ) {
	    break;
	  };
	  
	};
	best_peak_velocity = peak_velocity;
	double delta_time = start_point->time - end_point->time;
	if(min_delta_time > delta_time)
	  return;
	for(unsigned int i = 0; i < 20; ++i) {
	  PointType1 prev_peak_velocity = peak_velocity;
	  
	  peak_velocity = lift_to_space<1>(delta_first_order, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  double temp = space_1.distance(get<1>(start_point->pt), peak_velocity);
	  coefs[1] = (coefs[0] * coefs[0] + coefs[4] * coefs[4] - temp * temp) / (2.0 * coefs[0] * coefs[4]);
	  temp = space_1.distance(get<1>(end_point->pt), peak_velocity);
	  coefs[3] = (coefs[2] * coefs[2] + coefs[4] * coefs[4] - temp * temp) / (2.0 * coefs[2] * coefs[4]);
	  double beta_peak1 = coefs[0] * coefs[1] / (coefs[4] * coefs[4]);
	  double beta_peak2 = coefs[2] * coefs[3] / (coefs[4] * coefs[4]);
	  if(beta_peak1 > beta_peak2)
	    std::swap(beta_peak1, beta_peak2);
	  if(beta_peak2 > 1.0)
	    beta_peak2 = 1.0;
	  if((beta_peak1 > 0.0) && (detail::svp_compute_slack_time(beta_peak1,delta_time,norm_delta,coefs) > 0.0)) {
	    //there must be a root between peak1 and 0.
	    double beta_low = 0.0;
	    double beta_temp = beta_peak1;
	    bisection_method(beta_low, beta_peak1, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	    
	    if(detail::svp_compute_slack_time(beta_peak2,delta_time,norm_delta,coefs) < 0.0) {
	      //there must be a root between peak1 and peak2.
	      bisection_method(beta_temp, beta_peak2, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	    } else if(detail::svp_compute_slack_time(1.0,delta_time,norm_delta,coefs) < 0.0) {
	      //there must be a root between peak2 and 1.0.
	      beta_temp = beta_peak2;
	      beta_peak2 = 1.0;
	      bisection_method(beta_temp, beta_peak2, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	    } else {
	      beta_peak2 = 5.0;
	    };
	    
	  } else {
	    double beta_low = 0.0;
	    if(beta_peak1 > 0.0) {
	      beta_low = beta_peak1;
	    };
	    if((beta_peak2 > 0.0) && (detail::svp_compute_slack_time(beta_peak2,delta_time,norm_delta,coefs) > 0.0)) {
	      //there must be a root between beta_low and beta_peak2.
	      beta_peak1 = beta_peak2;
	      bisection_method(beta_low, beta_peak1, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	      if(detail::svp_compute_slack_time(1.0,delta_time,norm_delta,coefs) < 0.0) {
		//there must be a root between beta_peak2 and 1.0.
		beta_low = beta_peak2;
		beta_peak2 = 1.0;
		bisection_method(beta_low,beta_peak2, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	      } else {
		beta_peak2 = 5.0;
	      };
	    } else {
	      if(beta_peak2 > beta_low)
		beta_low = beta_peak2;
	      beta_peak1 = 1.0;
	      if(detail::svp_compute_slack_time(1.0,delta_time,norm_delta,coefs) > 0.0) {
		//there must be a root between beta_peak2 and 1.0.
		bisection_method(beta_low, beta_peak1, boost::bind(detail::svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(norm_delta),boost::cref(coefs)), 1e-6);
	      };
	      beta_peak2 = 5.0;
	    };
	  };
	  
	  if( std::fabs(beta - beta_peak1) < std::fabs(beta - beta_peak2) ) {
	    beta = beta_peak1;
	  } else {
	    beta = beta_peak2;
	  };
	  
	  peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  
	  PointDiff0 descended_peak_velocity = descend_to_space<0>(peak_velocity,1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology());
	  delta_first_order = space_0.difference(
	    space_0.adjust(get<0>(end_point->pt), (-0.5 * space_1.distance(peak_velocity, get<1>(end_point->pt))) * ( descended_peak_velocity
	                                                                                                             + descend_to_space<0>(get<1>(end_point->pt),1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()) )),
	    space_0.adjust(get<0>(start_point->pt), (0.5 * space_1.distance(peak_velocity, get<1>(start_point->pt))) * (descended_peak_velocity
	                                                                                                                + descend_to_space<0>(get<1>(start_point->pt),1.0,parent->get_temporal_space()->get_space_topology(),parent->get_temporal_space()->get_time_topology()) ))
	  );
	  norm_delta = space_0.norm(delta_first_order);
	  peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology());
	  
	  if( space_1.distance(prev_peak_velocity, peak_velocity) < 1e-6 * space_1.get_radius() ) {
	    break;
	  };
	};
	
      };
    };
  
  public:
    
    
    /**
     * Default constructor.
     */
    svp_interpolator(const Factory* aParent = NULL, const point_type* aStart = NULL, const point_type* aEnd = NULL) :
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
        throw singularity_error("Normalizing factor in linear interpolation is zero!");
      double t_normal = (t - start_point->time) / t_factor;
      
      point_type result;
      result.time = t;
      
      detail::linear_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< typename temporal_topology_traits<topology>::space_topology >::order> >(result.pt, start_point->pt, delta_first_order, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology(), t_factor, t_normal);
  
      return result;   
    };
    
};

/**
 * This class is a factory class for sustained velocity pulse (SVP) interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template <typename TemporalTopology>
class svp_interpolator_factory : public serialization::serializable {
  public:
    typedef svp_interpolator_factory<TemporalTopology> self;
    typedef TemporalTopology topology;
    typedef typename temporal_topology_traits<TemporalTopology>::point_type point_type;
  
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<TemporalTopology>::space_topology, 1, typename temporal_topology_traits<TemporalTopology>::time_topology >));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept< derived_N_order_space<typename temporal_topology_traits<TemporalTopology>::space_topology, typename temporal_topology_traits<TemporalTopology>::time_topology, 1> >));
    
    typedef typename derived_N_order_space< typename temporal_topology_traits<topology>::space_topology,
                                            typename temporal_topology_traits<topology>::time_topology,0>::type Space0;
    
    typedef svp_interpolator<self, metric_topology_traits< Space0 >::dimensions> interpolator_type;
    
  private:
    const topology* p_space;
  public:
    linear_interpolator_factory(const topology* aPSpace = NULL) : p_space(aPSpace) { };
  
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

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2420001,1,"linear_interpolator_factory",serialization::serializable)
};




  
/**
 * This class implements a trajectory in a temporal and zero-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class linear_interp_traj : public interpolated_trajectory<Topology,linear_interpolator_factory<Topology>,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 1, typename temporal_topology_traits<Topology>::time_topology >));
    
    typedef linear_interp_traj<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,linear_interpolator_factory<Topology>,DistanceMetric> base_class_type;
    
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
    explicit linear_interp_traj(const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                base_class_type(aSpace, aDist, linear_interpolator_factory<Topology>(&aSpace)) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    linear_interp_traj(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                       base_class_type(aSpace, aStart, aEnd, aDist, linear_interpolator_factory<Topology>(&aSpace)) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    linear_interp_traj(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                       base_class_type(aBegin, aEnd, aSpace, aDist, linear_interpolator_factory<Topology>(&aSpace)) { };
    
};



};

};

#endif









