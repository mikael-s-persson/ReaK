/**
 * \file sustained_acceleration_pulse_detail.hpp
 * 
 * This library contains the implementation details of the rate-limited sustained acceleration 
 * pulse (SAP) interpolation.
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

#ifndef REAK_SUSTAINED_ACCELERATION_PULSE_DETAIL_HPP
#define REAK_SUSTAINED_ACCELERATION_PULSE_DETAIL_HPP

#include "lin_alg/arithmetic_tuple.hpp"

#include "lin_alg/mat_num_exceptions.hpp"
#include "lin_alg/vect_alg.hpp"
#include "path_planning/bounded_space_concept.hpp"
#include "path_planning/differentiable_space_concept.hpp"

#include "root_finders/bisection_method.hpp"

#include "sustained_velocity_pulse_detail.hpp"

#include <boost/bind.hpp>

#include <cmath>

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  template <typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type sap_interpolate_HOT_impl(PointType& result, const PointDiff1& delta_second_order,
                                        const DiffSpace& space, const TimeSpace& t_space,
				 	double t_factor) {
    /* nothing to do. */
  };
  
  template <typename Idx, typename PointType, typename PointDiff1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type sap_interpolate_HOT_impl(PointType& result, const PointDiff1& delta_second_order,
                                        const DiffSpace& space, const TimeSpace& t_space,
					double t_factor) {
    get<2>(result) = lift_to_space<2>(delta_second_order, t_factor, space, t_space);
  };
  
  
  
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_constant_jerk_motion_HOT_impl(PointType& result, const PointDiff2& descended_jerk,
                                                 const DiffSpace& space, const TimeSpace& t_space) {
    /* Nothing to do. */ 
  };
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::equal_to<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_constant_jerk_motion_HOT_impl(PointType& result, const PointDiff2& descended_jerk,
                                                 const DiffSpace& space, const TimeSpace& t_space) {
    get<3>(result) = lift_to_space<3>(descended_jerk,1.0,space,t_space);
  };
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::greater<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_constant_jerk_motion_HOT_impl(PointType& result, const PointDiff2& descended_jerk,
                                                 const DiffSpace& space, const TimeSpace& t_space) {
    sap_constant_jerk_motion_HOT_impl< typename boost::mpl::prior<Idx>::type >(result, descended_jerk, space, t_space);
    
    get<Idx::type::value>(result) = get_space<Idx::type::value>(space,t_space).origin();
  };
  
  template <typename Idx, typename PointType, typename PointDiff2, typename DiffSpace, typename TimeSpace>
  inline
  void sap_constant_jerk_motion_impl(PointType& result, const PointDiff2& descended_jerk,
                                             const DiffSpace& space, const TimeSpace& t_space,
					     double dt) {
      
    get<0>(result) = 
      get_space<0>(space,t_space).adjust(
	get<0>(result),
        descend_to_space<0>( get_space<1>(space,t_space).adjust( get<1>(result), 
          descend_to_space<1>( get_space<2>(space,t_space).adjust( get<2>(result), (1.0 / 6.0) * dt * descended_jerk ), dt, space, t_space) ), dt, space, t_space)
      );
      
    
    get<1>(result) = 
      get_space<1>(space,t_space).adjust( get<1>(result), 
        descend_to_space<1>( get_space<2>(space,t_space).adjust( get<2>(result), 0.5 * dt * descended_jerk ), dt, space, t_space) );

    get<2>(result) = get_space<2>(space,t_space).adjust( get<2>(result), dt * descended_jerk );
    
    if(Idx::type::value > 2) 
      sap_constant_jerk_motion_HOT_impl(result,descended_jerk,space,t_space);
    
    return;
  };
  
  
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<4> 
    >,
  void >::type sap_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
				    const PointDiff0& delta_first_order, const PointType1& peak_velocity,
                                    const DiffSpace& space, const TimeSpace& t_space,
				    double dt, double dt_total) {
    using std::sqrt;
    
    result = start_point;
    
    double dt_amax = get_space<2>(space,t_space).get_radius();
    
    double dt_vp_1st = get_space<1>(space,t_space).distance(get<1>(start_point), peak_velocity);
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp = dt_vp_1st - dt_amax;
    double dt_ap = dt_amax;
    if( dt_vp < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp = 0.0;
      dt_ap = sqrt(4.0 * dt_amax * dt_vp_1st) * 0.5;
    };
    
    double dt_vp2_1st = get_space<1>(space,t_space).distance(peak_velocity, get<1>(end_point));
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - dt_amax;
    double dt_ap2 = dt_amax;
    if( dt_vp2 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt(4.0 * dt_amax * dt_vp2_1st) * 0.5;
    };
    dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp + 2.0 * dt_ap;
    
    //Phase 1: in the jerk-up phase of velocity ramp-up.
    if( dt < dt_ap)
      dt_ap = dt;
    
    sap_constant_jerk_motion_impl<Idx>(result,
      get_space<2>(space,t_space).difference( 
        lift_to_space<2>(get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point)),dt_vp_1st * dt_amax,space,t_space),
	get_space<2>(space,t_space).origin()
      ),
      space, t_space, dt_ap
    );
    dt -= dt_ap;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    //Phase 2: in the constant accel phase of velocity ramp-up.
    if( dt < dt_vp )
      dt_vp = dt;
    
    svp_constant_accel_motion_impl<Idx>(result, descend_to_space<1>(get<2>(result),1.0,space,t_space), space, t_space, dt_vp);
    dt -= dt_vp;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    //Phase 3: in the jerk-down phase of velocity ramp-up.
    if( dt < dt_ap )
      dt_ap = dt;
    
    sap_constant_jerk_motion_impl<Idx>(result,
      get_space<2>(space,t_space).difference( 
        get_space<2>(space,t_space).origin(),
	lift_to_space<2>(get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point)),dt_vp_1st * dt_amax,space,t_space)
      ), space, t_space, dt_ap);
    dt -= dt_ap;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    
    //TODO case 4: in the cruise phase.
    if( dt < dt_total )
      dt_total = dt;
    
    svp_constant_vel_motion_impl<Idx>(result, descend_to_space<0>(get<1>(result),1.0,space,t_space), space, t_space, dt_total);
    dt -= dt_total;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    //TODO case 5: in the jerk-up phase of velocity ramp-down.
    if( dt < dt_ap2)
      dt_ap2 = dt;
    
    sap_constant_jerk_motion_impl<Idx>(result,
      get_space<2>(space,t_space).difference( 
        lift_to_space<2>(get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity),dt_vp2_1st * dt_amax,space,t_space),
	get_space<2>(space,t_space).origin()
      ),
      space, t_space, dt_ap2
    );
    dt -= dt_ap2;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    //TODO case 6: in the constant accel phase of velocity ramp-down.
    if( dt < dt_vp2 )
      dt_vp2 = dt;
    
    svp_constant_accel_motion_impl<Idx>(result, descend_to_space<1>(get<2>(result),1.0,space,t_space), space, t_space, dt_vp2);
    dt -= dt_vp2;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    //TODO case 7: in the jerk-down phase of velocity ramp-down.
    if( dt < dt_ap2 )
      dt_ap2 = dt;
    
    sap_constant_jerk_motion_impl<Idx>(result,
      get_space<2>(space,t_space).difference( 
        get_space<2>(space,t_space).origin(),
	lift_to_space<2>(get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity),dt_vp2_1st * dt_amax,space,t_space)
      ), space, t_space, dt_ap2);
    dt -= dt_ap2;
    
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type sap_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
				    const PointDiff0& delta_first_order, const PointType1& peak_velocity, 
				    const DiffSpace& space, const TimeSpace& t_space,
				    double dt, double dt_total) {
    sap_interpolate_impl< typename boost::mpl::prior<Idx>::type >(result,start_point,end_point,delta_first_order,peak_velocity,space,t_space,dt,dt_total);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
  inline
  double sap_compute_slack_time(double beta, double dt, double beta_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    
    if(term1 < dt_amax)
      term1 = sqrt(4.0 * term1 * dt_amax);
    else
      term1 += dt_amax;
    if(term2 < dt_amax)
      term2 = sqrt(4.0 * term2 * dt_amax);
    else
      term2 += dt_amax;
    
    double c = (norm_delta + coefs[4] * (beta_0 * beta_0 - beta * beta) );
    
    return dt - term1 - term2 - fabs(c) / beta;
  };
    
  inline
  double sap_compute_derivative_slack_time(double beta, double dt, double beta_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double c = (norm_delta + coefs[4] * (beta_0 * beta_0 - beta * beta) );
    
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    
    double dterm1 = coefs[4] * (coefs[4] * beta - coefs[0] * coefs[1]) / term1;
    double dterm2 = coefs[4] * (coefs[4] * beta - coefs[2] * coefs[3]) / term2;
    if(term1 < dt_amax)
      dterm1 *= sqrt(dt_amax / term1);
    if(term2 < dt_amax)
      dterm2 *= sqrt(dt_amax / term2);
    
    return fabs(c) / (beta * beta) - ( c > 0.0 ? -2.0 : 2.0) * coefs[4] - dterm1 - dterm2;
  };
  
  
  inline
  double sap_compute_travel_time(double beta, double beta_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double c = (norm_delta + coefs[4] * (beta_0 * beta_0 - beta * beta) );
    
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    
    if(term1 < dt_amax)
      term1 = sqrt(4.0 * term1 * dt_amax);
    else
      term1 += dt_amax;
    if(term2 < dt_amax)
      term2 = sqrt(4.0 * term2 * dt_amax);
    else
      term2 += dt_amax;
    
    return term1 + term2 + fabs(c) / beta;
  };
    
  inline
  double sap_compute_derivative_travel_time(double beta, double beta_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double c = (norm_delta + coefs[4] * (beta_0 * beta_0 - beta * beta) );
    
    double term1 = sqrt(coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[0] * coefs[1] - beta * coefs[4]));
    double term2 = sqrt(coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[2] * coefs[3] - beta * coefs[4]));
    
    double dterm1 = coefs[4] * (coefs[4] * beta - coefs[0] * coefs[1]) / term1;
    double dterm2 = coefs[4] * (coefs[4] * beta - coefs[2] * coefs[3]) / term2;
    if(term1 < dt_amax)
      dterm1 *= sqrt(dt_amax / term1);
    if(term2 < dt_amax)
      dterm2 *= sqrt(dt_amax / term2);
    
    return dterm1 + dterm2 - fabs(c) / (beta * beta) + ( c > 0.0 ? -2.0 : 2.0) * coefs[4];
  };
  
  struct sap_update_delta_first_order {
    template <typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
    void operator()(const PointType& start_point, const PointType& end_point,
		    PointDiff0& delta_first_order, PointType1& peak_velocity,
		    double& norm_delta, const double& beta,
		    const DiffSpace& space, const TimeSpace& t_space,
		    double num_tol = 1E-6) {
    
      double dt_amax = get_space<2>(space,t_space).get_radius();
      
      double dt_vp1_1st = get_space<1>(space,t_space).distance(get<1>(start_point), peak_velocity);
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp1 = dt_vp1_1st - dt_amax;
      double dt_ap1 = dt_amax;
      if( dt_vp1 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp1 = 0.0;
        dt_ap1 = sqrt(4.0 * dt_amax * dt_vp1_1st) * 0.5;
      };
      
      double dt_vp2_1st = get_space<1>(space,t_space).distance(peak_velocity, get<1>(end_point));
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp2 = dt_vp2_1st - dt_amax;
      double dt_ap2 = dt_amax;
      if( dt_vp2 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp2 = 0.0;
        dt_ap2 = sqrt(4.0 * dt_amax * dt_vp2_1st) * 0.5;
      };
      
      
      PointDiff0 start_to_peak = get_space<0>(space,t_space).difference(get<0>(start_point),get<0>(start_point));
      
      if( dt_vp1_1st > num_tol * get_space<1>(space,t_space).get_radius() ) {
        start_to_peak =
          descend_to_space<0>( get_space<1>(space,t_space).adjust( get<1>(start_point),
	    (0.5 * (dt_ap1 + dt_vp1 + dt_vp1 * dt_ap1 / (dt_vp1 + dt_ap1)) * dt_ap1 / (dt_vp1_1st * dt_amax)) * 
	    get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point))), dt_vp1 + dt_ap1, space, t_space);
      };
      
      PointDiff0 peak_to_end = get_space<0>(space,t_space).difference(get<0>(end_point),get<0>(end_point));
      
      if( dt_vp2_1st > num_tol * get_space<0>(space,t_space).get_radius() ) {
        peak_to_end =
          descend_to_space<0>( get_space<1>(space,t_space).adjust( peak_velocity,
	    (0.5 * (dt_ap2 + dt_vp2 + dt_vp2 * dt_ap2 / (dt_vp2 + dt_ap2)) * dt_ap2 / (dt_vp2_1st * dt_amax)) * 
	    get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity)), dt_vp2 + dt_ap2, space, t_space);
      };
      
      delta_first_order = get_space<0>(space,t_space).difference(
        get_space<0>(space,t_space).adjust(get<0>(end_point), - peak_to_end ),
        get_space<0>(space,t_space).adjust(get<0>(start_point), start_to_peak)
      );
      
      norm_delta = get_space<0>(space,t_space).norm(delta_first_order);
      if(norm_delta > num_tol)
        peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, space, t_space);
      
      PointDiff0 descended_peak_velocity = descend_to_space<0>(peak_velocity,1.0,space,t_space);
      if( get_space<0>(space,t_space).norm(descended_peak_velocity - delta_first_order) > get_space<0>(space,t_space).norm(descended_peak_velocity) )
        norm_delta = -norm_delta;
    };
  };
  
  double sap_solve_for_min_dt_beta(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double dt_amax) {
    if(sap_compute_derivative_travel_time(1.0,beta,norm_delta,coefs,dt_amax) > 0.0) {
      double upper = 1.0;
      double lower = 0.5;
      while((lower < 0.99) && (sap_compute_derivative_travel_time(lower,beta,norm_delta,coefs,dt_amax) > 0.0)) {
        lower += 0.5 * (upper - lower);
      };
      if(lower < 0.99) {
        bisection_method(lower, upper, boost::bind(sap_compute_derivative_travel_time,_1,boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
      } else {
        upper = 1.0;
        lower = 0.5;
        while(sap_compute_derivative_travel_time(lower,beta,norm_delta,coefs) > 0.0) {
          upper = lower;
          lower *= 0.5;
        };
        bisection_method(lower, upper, boost::bind(sap_compute_derivative_travel_time,_1,boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
      };
      //make sure that the second root does not cause a reversal of the travel direction:
      if( (norm_delta > 0.0) && ( norm_delta + coefs[4] * (beta * beta - upper * upper) < 0.0 )) {
        upper = std::sqrt(beta * beta + norm_delta / coefs[4]);
      };
      beta = upper;
    } else {
      beta = 1.0;
    };
    return beta;
  };
  
  bool sap_min_dt_predicate(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double& result, double dt_amax) {
    result = sap_compute_travel_time(beta,beta,norm_delta,coefs,dt_amax);
    RK_NOTICE(1,"   gives min-dt = " << result);
    return true;
  };
  
  double sap_solve_for_no_slack_beta(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double delta_time, double dt_amax) {
    double beta_peak1 = 1.0;
    double beta_peak2 = 5.0;
      
    if(svp_compute_slack_time(1.0,delta_time,beta,norm_delta,coefs,dt_amax) > 0.0) {
      //means I have a single root in the interval, so I can solve for it:
      double beta_low = 0.5;
      while(svp_compute_slack_time(beta_low,delta_time,beta,norm_delta,coefs,dt_amax) > 0.0) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method(beta_low, beta_peak1, boost::bind(svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
    } else {
      //This means that I must have either a parabola-looking curve, or it never goes positive.
      // so, find the maximum in the interval, by finding the zero of the derivative:
      double beta_low = 0.5;
      while(svp_compute_derivative_slack_time(beta_low,delta_time,beta,norm_delta,coefs,dt_amax) < 0.0) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method(beta_low,beta_peak1, boost::bind(svp_compute_derivative_slack_time,_1,boost::cref(delta_time),boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
      if( svp_compute_slack_time(beta_peak1,delta_time,beta,norm_delta,coefs,dt_amax) > 0.0 ) {
        //this means the maximum slack-time is actually positive, meaning there must be a root on either side.
        beta_peak2 = beta_peak1;
        beta_low = 0.5 * beta_peak1;
        while(svp_compute_slack_time(beta_low,delta_time,beta,delta_time,coefs,dt_amax) > 0.0) {
          beta_peak1 = beta_low;
          beta_low *= 0.5;
        };
        bisection_method(beta_low, beta_peak1, boost::bind(svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
        beta_low = beta_peak2;
        beta_peak2 = 1.0;
        bisection_method(beta_low, beta_peak2, boost::bind(svp_compute_slack_time,_1,boost::cref(delta_time),boost::cref(beta),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
          
        //make sure that the second root does not cause a reversal of the travel direction:
        if( norm_delta + coefs[4] * (beta * beta - beta_peak2 * beta_peak2) < 0.0 ) {
          beta_peak2 = 5.0;
        };
      };
    };
      
    RK_NOTICE(1,"   delta-time = " << delta_time);
    RK_NOTICE(1,"   beta-peaks = (" << beta_peak1 << " " << beta_peak2 << ")");
    
    if( std::fabs(beta - beta_peak1) < std::fabs(beta - beta_peak2) ) {
      beta = beta_peak1;
    } else {
      beta = beta_peak2;
    };
    if(beta <= num_tol)
      beta = num_tol;
    return beta;
  };
  
  bool sap_no_slack_predicate(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double& slack, double delta_time, double dt_amax) {
    slack = svp_compute_slack_time(beta,delta_time,beta,norm_delta,coefs,dt_amax);
    RK_NOTICE(1,"   gives slack-time = " << slack);
    return (std::fabs(slack) < 100.0 * num_tol * delta_time );
  };
  
  
  
  
  template <typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  double sap_compute_min_delta_time(const PointType& start_point, const PointType& end_point,
                                    PointDiff0& delta_first_order, PointType1& peak_velocity,
				    double& norm_delta, double& beta,
				    const DiffSpace& space, const TimeSpace& t_space,
				    double num_tol = 1E-6, unsigned int max_iter = 20) {
    
    double result = 0.0;
    RK_NOTICE(1,"*********  Starting Min-dt Iterations    *****************");
    
    svp_peak_velocity_iteration(
      start_point, end_point,
      delta_first_order, peak_velocity,
      norm_delta, beta,
      space, t_space,
      sap_update_delta_first_order(),
      boost::bind(sap_solve_for_min_dt_beta,_1,_2,_3,_4,boost::cref(get_space<2>(space,t_space).get_radius())),
      boost::bind(sap_min_dt_predicate,_1,_2,_3,_4,boost::ref(result),boost::cref(get_space<2>(space,t_space).get_radius())),
      num_tol, max_iter
    );
    
    return result;
  };
  
  template <typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  double sap_compute_peak_velocity(const PointType& start_point, const PointType& end_point,
                                   PointDiff0& delta_first_order, PointType1& peak_velocity,
				   double& norm_delta, double& beta, double delta_time,
				   const DiffSpace& space, const TimeSpace& t_space,
				   double num_tol = 1E-6, unsigned int max_iter = 20) {
    
    double slack = 0.0;
    RK_NOTICE(1,"*********  Starting No-slack Iterations  *****************");
    
    svp_peak_velocity_iteration(
      start_point, end_point,
      delta_first_order, peak_velocity,
      norm_delta, beta,
      space, t_space,
      sap_update_delta_first_order(),
      boost::bind(sap_solve_for_no_slack_beta,_1,_2,_3,_4,boost::cref(delta_time),boost::cref(get_space<2>(space,t_space).get_radius())),
      boost::bind(sap_no_slack_predicate,_1,_2,_3,_4,boost::ref(slack),boost::cref(delta_time),boost::cref(get_space<2>(space,t_space).get_radius())),
      num_tol, max_iter
    );
    
    return slack;
  };
  
};


};

};

#endif









