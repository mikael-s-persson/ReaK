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
#include "path_planning/tangent_bundle_concept.hpp"

#include "root_finders/bisection_method.hpp"
#include "root_finders/secant_method.hpp"

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
      sap_constant_jerk_motion_HOT_impl<Idx>(result,descended_jerk,space,t_space);
    
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
    
    double dt_vp_1st = get(distance_metric, get_space<1>(space,t_space))(get<1>(start_point), peak_velocity, get_space<1>(space,t_space));
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp = dt_vp_1st - dt_amax;
    double dt_ap = dt_amax;
    if( dt_vp < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp = 0.0;
      dt_ap = sqrt(dt_amax * dt_vp_1st);
    };
    
    double dt_vp2_1st = get(distance_metric, get_space<1>(space,t_space))(peak_velocity, get<1>(end_point), get_space<1>(space,t_space));
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - dt_amax;
    double dt_ap2 = dt_amax;
    if( dt_vp2 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt(dt_amax * dt_vp2_1st);
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
  vect<double,2> sap_compute_projected_deltas(double beta, const vect<double,5>& coefs, double dt_amax, double& A_0, double& A_1) {
    using std::sqrt;
    vect<double,2> result;
    A_0 = sqrt( coefs[0] * coefs[0] - beta * coefs[4] * (2.0 * coefs[1] - beta * coefs[4]) );
    double t_0 = coefs[1] / coefs[4];
    if(A_0 < dt_amax)
      result[0] = 0.5 * sqrt( dt_amax * A_0 ) * ( t_0 + 3.0 * beta );
    else
      result[0] = 0.5 * ((dt_amax + A_0) * (beta + t_0) + dt_amax * dt_amax / A_0 * (beta - t_0));
    A_1 = sqrt( coefs[2] * coefs[2] - beta * coefs[4] * (2.0 * coefs[3] - beta * coefs[4]) );
    double t_1 = coefs[3] / coefs[4];
    if(A_1 < dt_amax)
      result[1] = 0.5 * sqrt( dt_amax * A_1 ) * (t_1 + 3.0 * beta);
    else
      result[1] = 0.5 * ((dt_amax + A_1) * (beta + t_1) + dt_amax * dt_amax / A_1 * (beta - t_1));
    return result;
  };

  inline
  vect<double,2> sap_compute_projected_deltas(double beta, const vect<double,5>& coefs, double dt_amax) {
    double A_0, A_1;
    return sap_compute_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
  };
  
  inline
  vect<double,2> sap_compute_derivative_projected_deltas(double beta, const vect<double,5>& coefs, double dt_amax, double A_0, double A_1) {
    using std::sqrt;
    using std::fabs;
    vect<double,2> result;
    if(A_0 < dt_amax) {
//       if((fabs( coefs[1] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[0] - coefs[1]) / coefs[4]) < 1e-6))
//         result[0] = sqrt(dt_amax * A_0);
//       else
//         result[0] = sqrt(dt_amax * A_0) * (1.5 + 0.25 * (3.0 * beta * coefs[4] + coefs[1]) * 
//                                                         (beta * coefs[4] - coefs[1]) / (A_0 * A_0));
      result[0] = beta * coefs[4] * (beta * coefs[4] - coefs[1]) / dt_amax 
                + (coefs[0] * coefs[0] - coefs[1] * coefs[1]) / dt_amax
                + 0.5 * dt_amax;
    } else {
      result[0] = beta * coefs[4] * (beta * coefs[4] - coefs[1]) / A_0 
                + 0.5 * ((1.0 + dt_amax * dt_amax / (A_0 * A_0)) * (coefs[0] * coefs[0] - coefs[1] * coefs[1]) / A_0
                         + dt_amax);
    };
    if(A_1 < dt_amax) {
//       if((fabs( coefs[3] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[2] - coefs[3]) / coefs[4]) < 1e-6))
//         result[1] = sqrt(dt_amax * A_1);
//       else
//         result[1] = sqrt(dt_amax * A_1) * (1.5 + 0.25 * (3.0 * beta * coefs[4] + coefs[3]) * 
//                                                         (beta * coefs[4] - coefs[3]) / (A_1 * A_1));
      result[1] = beta * coefs[4] * (beta * coefs[4] - coefs[3]) / dt_amax 
                + (coefs[2] * coefs[2] - coefs[3] * coefs[3]) / dt_amax
                + 0.5 * dt_amax;
    } else {
      result[1] = beta * coefs[4] * (beta * coefs[4] - coefs[3]) / A_1 
                + 0.5 * ((1.0 + dt_amax * dt_amax / (A_1 * A_1)) * (coefs[2] * coefs[2] - coefs[3] * coefs[3]) / A_1
                         + dt_amax);
    };
    return result;
  };
  
  inline
  double sap_compute_slack_time(double beta, double dt, const vect<double,2>& deltas_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double A_0, A_1;
    vect<double,2> deltas_1 = sap_compute_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    if(A_0 < dt_amax)
      A_0 = sqrt(4.0 * A_0 * dt_amax);
    else
      A_0 += dt_amax;
    if(A_1 < dt_amax)
      A_1 = sqrt(4.0 * A_1 * dt_amax);
    else
      A_1 += dt_amax;
    
    double c = (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]);
    
    return dt - A_0 - A_1 - fabs(c) / beta;
  };
    
  inline
  double sap_compute_derivative_slack_time(double beta, double dt, const vect<double,2>& deltas_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double A_0, A_1;
    vect<double,2> deltas_1 = sap_compute_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    double c = (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]);
    
    vect<double,2> deltas_dot_1 = sap_compute_derivative_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    double c_dot = -deltas_dot_1[0] - deltas_dot_1[1];
    
    double dt0, dt1;
    if(A_0 < dt_amax) {
//       if((fabs( coefs[1] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[0] - coefs[1]) / coefs[4]) < 1e-6))
//         dt0 = -0.5 * sqrt(dt_amax * A_0) / beta;
//       else
//         dt0 = sqrt(dt_amax / A_0) * coefs[4] * (coefs[4] * beta - coefs[1]) / A_0;
      dt0 = coefs[4] * (coefs[4] * beta - coefs[1]) / dt_amax;
    } else {
      dt0 = coefs[4] * (coefs[4] * beta - coefs[1]) / A_0;
    };
    if(A_1 < dt_amax) {
//       if((fabs( coefs[3] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[2] - coefs[3]) / coefs[4]) < 1e-6))
//         dt1 = -0.5 * sqrt(dt_amax * A_1) / beta;
//       else
//         dt1 = sqrt(dt_amax / A_1) * coefs[4] * (coefs[4] * beta - coefs[3]) / A_1;
      dt1 = coefs[4] * (coefs[4] * beta - coefs[3]) / dt_amax;
    } else {
      dt1 = coefs[4] * (coefs[4] * beta - coefs[3]) / A_1;
    };
    
    return fabs(c) / (beta * beta) - ( c > 0.0 ? c_dot : -c_dot) / beta - dt0 - dt1;
  };
  
  
  inline
  double sap_compute_travel_time(double beta, const vect<double,2>& deltas_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
    double A_0, A_1;
    vect<double,2> deltas_1 = sap_compute_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    if(A_0 < dt_amax)
      A_0 = sqrt(4.0 * A_0 * dt_amax);
    else
      A_0 += dt_amax;
    if(A_1 < dt_amax)
      A_1 = sqrt(4.0 * A_1 * dt_amax);
    else
      A_1 += dt_amax;
    
    double c = (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]);
    
    return A_0 + A_1 + fabs(c) / beta;
  };
    
  inline
  double sap_compute_derivative_travel_time(double beta, const vect<double,2>& deltas_0, double norm_delta, const vect<double,5>& coefs, double dt_amax) {
    using std::sqrt;
    using std::fabs;
    
        
    double A_0, A_1;
    vect<double,2> deltas_1 = sap_compute_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    double c = (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]);
    
    //RK_NOTICE(1,"  dt-der: norm-dp = " << norm_delta << " deltas_1 = " << deltas_1 << " deltas_0 = " << deltas_0 << " c = " << c);
    
    vect<double,2> deltas_dot_1 = sap_compute_derivative_projected_deltas(beta,coefs,dt_amax,A_0,A_1);
    
    double c_dot = -deltas_dot_1[0] - deltas_dot_1[1];
    
    double dt0, dt1;
    if(A_0 < dt_amax) {
//       if((fabs( coefs[1] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[0] - coefs[1]) / coefs[4]) < 1e-6))
//         dt0 = -0.5 * sqrt(dt_amax * A_0) / beta;
//       else
//         dt0 = sqrt(dt_amax / A_0) * coefs[4] * (coefs[4] * beta - coefs[1]) / A_0;
      dt0 = coefs[4] * (coefs[4] * beta - coefs[1]) / dt_amax;
    } else {
      dt0 = coefs[4] * (coefs[4] * beta - coefs[1]) / A_0;
    };
    if(A_1 < dt_amax) {
//       if((fabs( coefs[3] / coefs[4] - beta ) < 1e-6) && (fabs((coefs[2] - coefs[3]) / coefs[4]) < 1e-6))
//         dt1 = -0.5 * sqrt(dt_amax * A_1) / beta;
//       else
//         dt1 = sqrt(dt_amax / A_1) * coefs[4] * (coefs[4] * beta - coefs[3]) / A_1;
      dt1 = coefs[4] * (coefs[4] * beta - coefs[3]) / dt_amax;
    } else {
      dt1 = coefs[4] * (coefs[4] * beta - coefs[3]) / A_1;
    };
    
    //RK_NOTICE(1,"   dt-der iter");
    
    return dt0 + dt1 - fabs(c) / (beta * beta) + ( c > 0.0 ? c_dot : -c_dot) / beta;
  };
  
  struct sap_update_delta_first_order {
    template <typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
    void operator()(const PointType& start_point, const PointType& end_point,
		    PointDiff0& delta_first_order, PointType1& peak_velocity,
		    double& norm_delta, const double& beta,
		    const DiffSpace& space, const TimeSpace& t_space,
		    double num_tol = 1E-6) {
      using std::fabs;
      using std::sqrt;
      
      double dt_amax = get_space<2>(space,t_space).get_radius();
      
      double dt_vp1_1st = get(distance_metric, get_space<1>(space,t_space))(get<1>(start_point), peak_velocity, get_space<1>(space,t_space));
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp1 = dt_vp1_1st - dt_amax;
      double dt_ap1 = dt_amax;
      if( dt_vp1 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp1 = 0.0;
        dt_ap1 = sqrt(dt_amax * dt_vp1_1st);
      };
      
      double dt_vp2_1st = get(distance_metric, get_space<1>(space,t_space))(peak_velocity, get<1>(end_point), get_space<1>(space,t_space));
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp2 = dt_vp2_1st - dt_amax;
      double dt_ap2 = dt_amax;
      if( dt_vp2 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp2 = 0.0;
        dt_ap2 = sqrt(dt_amax * dt_vp2_1st);
      };
      
      
      PointDiff0 start_to_peak = get_space<0>(space,t_space).difference(get<0>(start_point),get<0>(start_point));
      
      if( dt_vp1_1st > num_tol * get_space<1>(space,t_space).get_radius() ) {
	start_to_peak = 
	  descend_to_space<0>( get<1>(start_point), dt_vp1, space, t_space)
	  + descend_to_space<0>( get_space<1>(space,t_space).adjust( get<1>(start_point),
	      ((0.75 * dt_ap1 * (dt_ap1 + dt_vp1) + 0.25 * dt_vp1 * dt_vp1) / (dt_amax * dt_vp1_1st)) * 
	      get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point))
            ),
	    2.0 * dt_ap1, space, t_space);
	
      };
      
      PointDiff0 peak_to_end = get_space<0>(space,t_space).difference(get<0>(end_point),get<0>(end_point));
      
      if( dt_vp2_1st > num_tol * get_space<1>(space,t_space).get_radius() ) {
	peak_to_end = 
	  descend_to_space<0>( peak_velocity, dt_vp2, space, t_space)
	  + descend_to_space<0>( get_space<1>(space,t_space).adjust( peak_velocity,
	      ((0.75 * dt_ap2 * (dt_ap2 + dt_vp2) + 0.25 * dt_vp2 * dt_vp2) / (dt_amax * dt_vp2_1st)) * 
	      get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity)
            ),
	    2.0 * dt_ap2, space, t_space);

      };
      
      delta_first_order = get_space<0>(space,t_space).difference(
        get_space<0>(space,t_space).adjust(get<0>(end_point), - peak_to_end ),
        get_space<0>(space,t_space).adjust(get<0>(start_point), start_to_peak)
      );
      
      norm_delta = get(distance_metric, get_space<0>(space,t_space))(delta_first_order, get_space<0>(space,t_space));
      PointDiff0 descended_peak_velocity = descend_to_space<0>(peak_velocity,1.0,space,t_space);
      double normA = get(distance_metric, get_space<0>(space,t_space))(descended_peak_velocity, get_space<0>(space,t_space));
      double normC = get(distance_metric, get_space<0>(space,t_space))(descended_peak_velocity - delta_first_order, get_space<0>(space,t_space));
      if( normC * normC > normA * normA + norm_delta * norm_delta )
        norm_delta = -norm_delta;
      if(fabs(norm_delta) > num_tol)
        peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta, space, t_space);
      
    };
  };
  
  inline 
  double sap_solve_for_min_dt_beta(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double dt_amax) {
    vect<double,2> deltas_0 = sap_compute_projected_deltas(beta,coefs,dt_amax);
    if(sap_compute_derivative_travel_time(1.0,deltas_0,norm_delta,coefs,dt_amax) > 0.0) {
      //RK_NOTICE(1,"   dt-rate at edge is " << sap_compute_derivative_travel_time(1.0,deltas_0,norm_delta,coefs,dt_amax));
      double upper = 1.0;
      double lower = 0.1;
      while((lower < 0.99) && (sap_compute_derivative_travel_time(lower,deltas_0,norm_delta,coefs,dt_amax) > 0.0)) {
        lower += 0.5 * (upper - lower);
      };
      if(lower < 0.99) {
	//RK_NOTICE(1,"  starting dt-der iterations....");
        brent_method(lower, upper, boost::bind(sap_compute_derivative_travel_time,_1,boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
	//RK_NOTICE(1,"  done.");
      } else {
        upper = 1.0;
        lower = 0.9;
        while(sap_compute_derivative_travel_time(lower,deltas_0,norm_delta,coefs,dt_amax) > 0.0) {
          upper = lower;
          lower *= 0.5;
        };
	//RK_NOTICE(1,"  starting dt-der iterations....");
        brent_method(lower, upper, boost::bind(sap_compute_derivative_travel_time,_1,boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
	//RK_NOTICE(1,"  done.");
      };
      //make sure that the second root does not cause a reversal of the travel direction:
      vect<double,2> deltas_1 = sap_compute_projected_deltas(upper,coefs,dt_amax);
      if( (norm_delta > 0.0) && ( (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]) < 0.0 )) {
        upper = std::sqrt(beta * beta + norm_delta / coefs[4]);
      };
      beta = upper;
    } else {
      beta = 1.0;
    };
    return beta;
  };
  
  inline 
  bool sap_min_dt_predicate(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double& result, double dt_amax) {
    result = sap_compute_travel_time(beta,sap_compute_projected_deltas(beta,coefs,dt_amax),norm_delta,coefs,dt_amax);
    //RK_NOTICE(1,"   gives min-dt = " << result);
    return true;
  };
  
  inline 
  double sap_solve_for_no_slack_beta(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double delta_time, double dt_amax) {
    double beta_peak1 = 1.0;
    double beta_peak2 = 5.0;
      
    vect<double,2> deltas_0 = sap_compute_projected_deltas(beta,coefs,dt_amax);
    
    if(sap_compute_slack_time(1.0,delta_time,deltas_0,norm_delta,coefs,dt_amax) > 0.0) {
      //means I have a single root in the interval, so I can solve for it:
      double beta_low = 0.5;
      while(sap_compute_slack_time(beta_low,delta_time,deltas_0,norm_delta,coefs,dt_amax) > 0.0) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method(beta_low, beta_peak1, boost::bind(sap_compute_slack_time,_1,boost::cref(delta_time),boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
    } else {
      //This means that I must have either a parabola-looking curve, or it never goes positive.
      // so, find the maximum in the interval, by finding the zero of the derivative:
      double beta_low = 0.5;
      while(sap_compute_derivative_slack_time(beta_low,delta_time,deltas_0,norm_delta,coefs,dt_amax) < 0.0) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method(beta_low,beta_peak1, boost::bind(sap_compute_derivative_slack_time,_1,boost::cref(delta_time),boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
      if( sap_compute_slack_time(beta_peak1,delta_time,deltas_0,norm_delta,coefs,dt_amax) > 0.0 ) {
        //this means the maximum slack-time is actually positive, meaning there must be a root on either side.
        beta_peak2 = beta_peak1;
        beta_low = 0.5 * beta_peak1;
        while(sap_compute_slack_time(beta_low,delta_time,deltas_0,delta_time,coefs,dt_amax) > 0.0) {
          beta_peak1 = beta_low;
          beta_low *= 0.5;
        };
        bisection_method(beta_low, beta_peak1, boost::bind(sap_compute_slack_time,_1,boost::cref(delta_time),boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
        beta_low = beta_peak2;
        beta_peak2 = 1.0;
        bisection_method(beta_low, beta_peak2, boost::bind(sap_compute_slack_time,_1,boost::cref(delta_time),boost::cref(deltas_0),boost::cref(norm_delta),boost::cref(coefs),boost::cref(dt_amax)), num_tol);
          
        //make sure that the second root does not cause a reversal of the travel direction:
        vect<double,2> deltas_1 = sap_compute_projected_deltas(beta_peak2,coefs,dt_amax);
        if( (norm_delta > 0.0) && ( (norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1]) < 0.0 )) {
          beta_peak2 = 5.0;
        };
      };
    };
      
    //RK_NOTICE(1,"   delta-time = " << delta_time);
    //RK_NOTICE(1,"   beta-peaks = (" << beta_peak1 << " " << beta_peak2 << ")");
    
    if( std::fabs(beta - beta_peak1) < std::fabs(beta - beta_peak2) ) {
      beta = beta_peak1;
    } else {
      beta = beta_peak2;
    };
    if(beta <= num_tol)
      beta = num_tol;
    return beta;
  };
  
  inline 
  bool sap_no_slack_predicate(double beta, double norm_delta, const vect<double,5>& coefs, double num_tol, double& slack, double delta_time, double dt_amax) {
    slack = sap_compute_slack_time(beta,delta_time,sap_compute_projected_deltas(beta,coefs,dt_amax),norm_delta,coefs,dt_amax);
    //RK_NOTICE(1,"   gives slack-time = " << slack);
    return (std::fabs(slack) < 100.0 * num_tol * delta_time );
  };
  
  
  
  
  template <typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  double sap_compute_min_delta_time(const PointType& start_point, const PointType& end_point,
                                    PointDiff0& delta_first_order, PointType1& peak_velocity,
				    double& norm_delta, double& beta,
				    const DiffSpace& space, const TimeSpace& t_space,
				    double num_tol = 1E-6, unsigned int max_iter = 20) {
    
    double result = 0.0;
    //RK_NOTICE(1,"*********  Starting Min-dt Iterations    *****************");
    
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
    //RK_NOTICE(1,"*********  Starting No-slack Iterations  *****************");
    
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
  
  template <typename PointType, typename DiffSpace, typename TimeSpace>
  double sap_compute_interpolation_data_impl(const PointType& start_point, const PointType& end_point,
					     typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,0>::type >::point_difference_type& delta_first_order,
					     typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type& peak_velocity,
					     const DiffSpace& space, const TimeSpace& t_space,
					     double delta_time = 0.0,
					     typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type* best_peak_velocity = NULL, 
					     double num_tol = 1e-6, unsigned int max_iter = 20) {
    
    delta_first_order = get_space<0>(space,t_space).difference( get<0>(end_point), get<0>(start_point) );
    double norm_delta = get(distance_metric, get_space<0>(space,t_space))( delta_first_order, get_space<0>(space,t_space) );
    double beta = 0.0;
    peak_velocity = get_space<1>(space,t_space).origin();
    
    double min_delta_time = sap_compute_min_delta_time(start_point, end_point, 
	                                               delta_first_order, peak_velocity,
						       norm_delta, beta,
						       space, t_space, num_tol, max_iter);
    
    if(best_peak_velocity)
      *best_peak_velocity = peak_velocity;
	
    if(min_delta_time > delta_time)
      return min_delta_time;
    
    beta = beta * min_delta_time / delta_time;
    peak_velocity = get_space<1>(space,t_space).adjust(
      get_space<1>(space,t_space).origin(),
	(min_delta_time / delta_time) *
	get_space<1>(space,t_space).difference(
	  peak_velocity,
	  get_space<1>(space,t_space).origin()
	)
      );
	
    sap_compute_peak_velocity(start_point, end_point, 
	                      delta_first_order, peak_velocity,
			      norm_delta, beta, delta_time,
			      space, t_space, num_tol, max_iter);
    
    return min_delta_time;
  };
  
  
  
};


};

};

#endif









