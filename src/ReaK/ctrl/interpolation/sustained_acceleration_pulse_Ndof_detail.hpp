/**
 * \file sustained_acceleration_pulse_Ndof_detail.hpp
 * 
 * This library contains the implementation details of the rate-limited sustained acceleration 
 * pulse (SAP) interpolation.
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

#ifndef REAK_SUSTAINED_ACCELERATION_PULSE_NDOF_DETAIL_HPP
#define REAK_SUSTAINED_ACCELERATION_PULSE_NDOF_DETAIL_HPP

#include "lin_alg/arithmetic_tuple.hpp"

#include "lin_alg/vect_alg.hpp"
#include "path_planning/tangent_bundle_concept.hpp"

#include "sustained_velocity_pulse_Ndof_detail.hpp"

#include <cmath>

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_Ndof_constant_jerk_motion_HOT_impl(PointType&, double, std::size_t,
                                                 const DiffSpace&, const TimeSpace&) {
    /* Nothing to do. */ 
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::equal_to<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_Ndof_constant_jerk_motion_HOT_impl(PointType& result, double descended_jerk, std::size_t i,
                                                 const DiffSpace& space, const TimeSpace& t_space) {
    typedef typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,2>::type >::point_type PointType2;
    const PointType2& max_jerk = get_space<3>(space,t_space).get_upper_corner();
    get<3>(result)[i] = descended_jerk * max_jerk[i];
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::greater<
      Idx,
      boost::mpl::size_t<3>
    >,
  void >::type sap_Ndof_constant_jerk_motion_HOT_impl(PointType& result, double descended_jerk, std::size_t i,
                                                      const DiffSpace& space, const TimeSpace& t_space) {
    sap_Ndof_constant_jerk_motion_HOT_impl< typename boost::mpl::prior<Idx>::type >(result, descended_jerk, i, space, t_space);
    
    get<Idx::type::value>(result)[i] = 0.0;
  };
  
  
  
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<4> 
    >,
  void >::type sap_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                         const PointDiff0& delta_first_order, const PointType1& peak_velocity,
                                         const DiffSpace& space, const TimeSpace& t_space,
                                         double dt, double dt_total) {
    using std::sqrt;
    typedef typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1>::type >::point_type PointDiff1;
    
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1>::type >::point_type max_velocity = get_space<1>(space,t_space).get_upper_corner();
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 2>::type >::point_type max_acceleration = get_space<2>(space,t_space).get_upper_corner();
    
    PointDiff1 dv1 = get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point));
    PointDiff1 dv2 = get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity);
    result = start_point;
    
    //double dt_amax = get_space<2>(space,t_space).get_radius();
    
    for(std::size_t i = 0; i < dv1.size(); ++i) {
      
      double dt_vp1_1st = fabs(dv1[i]);
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp1 = dt_vp1_1st - max_acceleration[i];
      double dt_ap1 = max_acceleration[i];
      if( dt_vp1 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp1 = 0.0;
        dt_ap1 = sqrt(max_acceleration[i] * dt_vp1_1st);
      };
      
      double dt_vp2_1st = fabs(dv2[i]);
      // we know that dt_vp_2nd = dt_vp_1st + dt_amax
      double dt_vp2 = dt_vp2_1st - max_acceleration[i];
      double dt_ap2 = max_acceleration[i];
      if( dt_vp2 < 0.0 ) {
        //means that we don't have time to reach the maximum acceleration:
        dt_vp2 = 0.0;
        dt_ap2 = sqrt(max_acceleration[i] * dt_vp2_1st);
      };
      
      double dt_tmp = dt;
      double dt_total_tmp = dt_total - dt_vp2 - 2.0 * dt_ap2 - dt_vp1 - 2.0 * dt_ap1;
      
      //Phase 1: in the jerk-up phase of velocity ramp-up.
      double descended_jerk = dv1[i] / dt_vp1_1st; // NOTE
      if( dt_tmp < dt_ap1 )
        dt_ap1 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + ( get<2>(result)[i] + (1.0 / 6.0) * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration[i] ) * dt_ap1 / max_velocity[i];
      get<1>(result)[i] += ( get<2>(result)[i] + 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration[i];
      get<2>(result)[i] += dt_ap1 * descended_jerk;
      
      if(Idx::type::value > 2) 
        sap_Ndof_constant_jerk_motion_HOT_impl<Idx>(result,descended_jerk,i,space,t_space);
      
      dt_tmp -= dt_ap1;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 2: in the constant accel phase of velocity ramp-up.
      double descended_accel = get<2>(result)[i] / max_acceleration[i];
      if( dt_tmp < dt_vp1 )
        dt_vp1 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + 0.5 * dt_vp1 * descended_accel ) * dt_vp1 / max_velocity[i];
      get<1>(result)[i] += dt_vp1 * descended_accel;
      
      if(Idx::type::value > 1) 
        svp_Ndof_constant_accel_motion_HOT_impl<Idx>(result, descended_accel, i, space, t_space);
      
      dt_tmp -= dt_vp1;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 3: in the jerk-down phase of velocity ramp-up.
      descended_jerk = -dv1[i] / dt_vp1_1st; // NOTE
      if( dt_tmp < dt_ap1 )
        dt_ap1 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + ( get<2>(result)[i] + (1.0 / 6.0) * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration[i] ) * dt_ap1 / max_velocity[i];
      get<1>(result)[i] += ( get<2>(result)[i] + 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration[i];
      get<2>(result)[i] += dt_ap1 * descended_jerk;
      
      if(Idx::type::value > 2) 
        sap_Ndof_constant_jerk_motion_HOT_impl<Idx>(result,descended_jerk,i,space,t_space);
      
      dt_tmp -= dt_ap1;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 4: in the cruise phase.
      double descended_vel = get<1>(result)[i] / max_velocity[i];
      if( dt_tmp < dt_total_tmp )
        dt_total_tmp = dt_tmp;
      
      get<0>(result)[i] += dt_total_tmp * descended_vel;
      //get<1>(result)[i] = constant.
      
      if(Idx::type::value > 1) 
        svp_Ndof_constant_vel_motion_HOT_impl<Idx>(result,i);
      
      dt_tmp -= dt_total_tmp;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 5: in the jerk-up phase of velocity ramp-down.
      descended_jerk = dv2[i] / dt_vp2_1st; // NOTE
      if( dt_tmp < dt_ap2 )
        dt_ap2 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + ( get<2>(result)[i] + (1.0 / 6.0) * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration[i] ) * dt_ap2 / max_velocity[i];
      get<1>(result)[i] += ( get<2>(result)[i] + 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration[i];
      get<2>(result)[i] += dt_ap2 * descended_jerk;
      
      if(Idx::type::value > 2) 
        sap_Ndof_constant_jerk_motion_HOT_impl<Idx>(result,descended_jerk,i,space,t_space);
      
      dt_tmp -= dt_ap2;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 6: in the constant accel phase of velocity ramp-down.
      descended_accel = get<2>(result)[i] / max_acceleration[i];
      if( dt_tmp < dt_vp2 )
        dt_vp2 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + 0.5 * dt_vp2 * descended_accel ) * dt_vp2 / max_velocity[i];
      get<1>(result)[i] += dt_vp2 * descended_accel;
      
      if(Idx::type::value > 1) 
        svp_Ndof_constant_accel_motion_HOT_impl<Idx>(result, descended_accel, i, space, t_space);
      
      dt_tmp -= dt_vp2;
      if(dt_tmp <= std::numeric_limits<double>::epsilon())
        continue;
      
      //Phase 7: in the jerk-down phase of velocity ramp-down.
      descended_jerk = -dv2[i] / dt_vp2_1st; // NOTE
      if( dt_tmp < dt_ap2 )
        dt_ap2 = dt_tmp;
      
      get<0>(result)[i] += ( get<1>(result)[i] + ( get<2>(result)[i] + (1.0 / 6.0) * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration[i] ) * dt_ap2 / max_velocity[i];
      get<1>(result)[i] += ( get<2>(result)[i] + 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration[i];
      get<2>(result)[i] += dt_ap2 * descended_jerk;
      
      if(Idx::type::value > 2) 
        sap_Ndof_constant_jerk_motion_HOT_impl<Idx>(result,descended_jerk,i,space,t_space);
      
    };
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type sap_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                         const PointDiff0& delta_first_order, const PointType1& peak_velocity, 
                                         const DiffSpace& space, const TimeSpace& t_space,
                                         double dt, double dt_total) {
    sap_Ndof_interpolate_impl< typename boost::mpl::prior<Idx>::type >(result,start_point,end_point,delta_first_order,peak_velocity,space,t_space,dt,dt_total);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
  
  
  
  inline
  double sap_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& delta_first_order, double& peak_velocity, 
                                         double max_velocity, double max_acceleration, double& norm_delta) {
    using std::fabs;
    using std::sqrt;
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    double sa_v0_v1_2 = sign_p1_p0 * max_acceleration * (start_velocity + end_velocity) * 0.5;
    double vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    double vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    
    // try to assume that a0 and a1 are a_max
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //           + vp/v_max (dt - fabs(vp - v0) - fabs(vp - v1) - 2*a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    delta_first_order = end_position - start_position
      - 0.5 * (fabs(peak_velocity -   end_velocity) + max_acceleration) * (peak_velocity +   end_velocity) / max_velocity
      - 0.5 * (fabs(peak_velocity - start_velocity) + max_acceleration) * (peak_velocity + start_velocity) / max_velocity;
    norm_delta = fabs(delta_first_order);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return norm_delta + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + 2.0 * max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    delta_first_order = 0.0; norm_delta = 0.0;
    //   this assumes that we have 
    //    p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    double descrim = max_acceleration * max_acceleration - 4.0 * (sa_v0_v1_2 - vm_p1_p0 - vsqr_avg);
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration - 4.0 * (vm_p1_p0 - sa_v0_v1_2 - vsqr_avg);
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
    };
    
    
    // if this point is reached, we have to assume that at least one of the acceleration pulse is too short.
    
    // try to assume that fabs(a0) = fabs(vp - v0) and fabs(a1) = a_max.
    // giving:
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max 
    //           + vp/v_max (dt - 2*fabs(vp - v0) - fabs(vp - v1) - a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    delta_first_order = end_position - start_position
      - 0.5 * (fabs(peak_velocity -   end_velocity) + max_acceleration) * (peak_velocity +   end_velocity) / max_velocity
      - fabs(peak_velocity - start_velocity) * (peak_velocity + start_velocity) / max_velocity;
    norm_delta = fabs(delta_first_order);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return norm_delta + fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    delta_first_order = 0.0; norm_delta = 0.0;
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * max_acceleration * end_velocity - 2.0 * vm_p1_p0 - 2.0 * start_velocity * start_velocity - end_velocity * end_velocity);
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * ( 2.0 * vm_p1_p0 + 2.0 * start_velocity * start_velocity + end_velocity * end_velocity - sign_p1_p0 * max_acceleration * end_velocity );
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    
    // try to assume that fabs(a1) = fabs(vp - v1) and fabs(a0) = a_max.
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max
    //           + vp/v_max (dt - fabs(vp - v0) - 2*fabs(vp - v1) - a_max)
    //           + fabs(vp - v1) (vp + v1) / v_max 
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    delta_first_order = end_position - start_position
      - fabs(peak_velocity -   end_velocity) * (peak_velocity +   end_velocity) / max_velocity
      - 0.5 * (fabs(peak_velocity -   start_velocity) + max_acceleration) * (peak_velocity + start_velocity) / max_velocity;
    norm_delta = fabs(delta_first_order);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return norm_delta + 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    delta_first_order = 0.0; norm_delta = 0.0;
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v1) (vp + v1) / v_max + (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * max_acceleration * start_velocity - 2.0 * vm_p1_p0 - start_velocity * start_velocity - 2.0 * end_velocity * end_velocity);
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        // this means that this root works for this case.
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * ( 2.0 * vm_p1_p0 + start_velocity * start_velocity + 2.0 * end_velocity * end_velocity - sign_p1_p0 * max_acceleration * start_velocity );
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    
    // try to assume that fabs(a0) = fabs(vp - v0) and fabs(a1) = fabs(vp - v1).
    // giving:
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max 
    //           + vp/v_max (dt - 2*fabs(vp - v0) - 2*fabs(vp - v1))
    //           + fabs(vp - v1) (vp + v1) / v_max 
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    delta_first_order = end_position - start_position
      - fabs(peak_velocity -   end_velocity) * (peak_velocity +   end_velocity) / max_velocity
      - fabs(peak_velocity -   start_velocity) * (peak_velocity + start_velocity) / max_velocity;
    norm_delta = fabs(delta_first_order);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return norm_delta + 2.0 * fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity);
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    delta_first_order = 0.0; norm_delta = 0.0;
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max + fabs(vp - v1) (vp + v1) / v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    peak_velocity = sqrt(0.5 * vm_p1_p0 + vsqr_avg);
    if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= sign_p1_p0 * start_velocity) && (peak_velocity >= sign_p1_p0 * end_velocity)) {
      // this means that this root works for this case.
      peak_velocity *= sign_p1_p0;
      return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = vsqr_avg - 0.5 * vm_p1_p0;
    if( descrim > 0.0 ) {
      // this means there exists a root for vp:
      peak_velocity = sqrt(descrim);
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity <= sign_p1_p0 * start_velocity) && (peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        peak_velocity *= sign_p1_p0;
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      if((fabs(peak_velocity) <= max_velocity) && (-peak_velocity <= sign_p1_p0 * start_velocity) && (-peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        peak_velocity *= -sign_p1_p0;
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    
    // What the fuck!! This point should never be reached, unless the motion is completely impossible:
    peak_velocity = 0.0;
    return std::numeric_limits<double>::infinity();
  };
  
  inline
  void sap_Ndof_compute_peak_velocity(double start_position, double end_position,
                                      double start_velocity, double end_velocity,
                                      double& peak_velocity, double max_velocity, 
                                      double max_acceleration, double delta_time) {
    // NOTE: Assume that delta-time is larger than minimum reachable delta-time 
    //       (avoid checking for that, even if it could mean that vp is higher than maximum)
    
    using std::fabs;
    using std::sqrt;
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    
    
    // try to assume that a0 and a1 are a_max
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //           + vp/v_max (dt - fabs(vp - v0) - fabs(vp - v1) - 2*a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    // try to assume that vp more in the direction (p1-p0) than both v1 and v0 (i.e., ramp-up and ramp-down).
    double v0_v1_ts_sa = start_velocity + end_velocity + (delta_time - max_acceleration) * sign_p1_p0;
    double vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    double a_v0_v1_2 = max_acceleration * (start_velocity + end_velocity) * 0.5;
    double vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    double descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // try to assume that vp is somewhere between v1 and v0 (i.e., ramp-up ramp-up or ramp-down ramp-down).
    if(end_velocity > start_velocity) {
      // ramp-up ramp-up
      v0_v1_ts_sa = start_velocity - end_velocity + delta_time - max_acceleration;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = 0.5 * (start_velocity * start_velocity - end_velocity * end_velocity);
      if(fabs(v0_v1_ts_sa) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg - a_v0_v1_2) / v0_v1_ts_sa;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= start_velocity) && (peak_velocity <= end_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    } else {
      // ramp-down ramp-down
      v0_v1_ts_sa = end_velocity - start_velocity + delta_time - max_acceleration;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = 0.5 * (end_velocity * end_velocity - start_velocity * start_velocity);
      if(fabs(v0_v1_ts_sa) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg - a_v0_v1_2) / v0_v1_ts_sa;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    };
    
    // try to assume (the last case) that vp is less in the direction (p1-p0) and both v0 and v1.
    v0_v1_ts_sa = start_velocity + end_velocity + (max_acceleration - delta_time) * sign_p1_p0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity * sign_p1_p0 <= start_velocity * sign_p1_p0) && (peak_velocity * sign_p1_p0 <= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity * sign_p1_p0 <= start_velocity * sign_p1_p0) && (peak_velocity * sign_p1_p0 <= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    
    // if this point is reached, we have to assume that at least one of the acceleration pulse is too short.
    
    // try to assume that fabs(a0) = fabs(vp - v0) and fabs(a1) = a_max.
    // giving:
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max 
    //           + vp/v_max (dt - 2*fabs(vp - v0) - fabs(vp - v1) - a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    // try to assume that vp more in the direction (p1-p0) than both v1 and v0 (i.e., ramp-up and ramp-down).
    v0_v1_ts_sa = (2.0 * start_velocity + end_velocity + (delta_time - 0.5 * max_acceleration) * sign_p1_p0) * 2.0 / 3.0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position) * 2.0 / 3.0;
    a_v0_v1_2 = max_acceleration * end_velocity * sign_p1_p0 / 3.0;
    vsqr_avg = (2.0 * start_velocity * start_velocity + end_velocity * end_velocity) / 3.0;
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - a_v0_v1_2 + vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // try to assume that vp is somewhere between v1 and v0 (i.e., ramp-up ramp-up or ramp-down ramp-down).
    if(end_velocity > start_velocity) {
      // ramp-up ramp-up
      v0_v1_ts_sa = 2.0 * (2.0 * start_velocity - end_velocity + delta_time - 0.5 * max_acceleration);
      vm_p1_p0 = 2.0 * max_velocity * (end_position - start_position);
      a_v0_v1_2 = max_acceleration * end_velocity;
      vsqr_avg = 2.0 * start_velocity * start_velocity - end_velocity * end_velocity;
      descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - a_v0_v1_2 + vm_p1_p0);
      if(descrim > 0.0) {
        peak_velocity = (v0_v1_ts_sa + sqrt(descrim)) * 0.5;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= start_velocity) && (peak_velocity <= end_velocity)) {
          return;
        };
        peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= start_velocity) && (peak_velocity <= end_velocity)) {
          return;
        };
        // then, the solution doesn't fit the assumption.
      };
    } else {
      // ramp-down ramp-down
      v0_v1_ts_sa = 2.0 * (end_velocity - 2.0 * start_velocity + delta_time - 0.5 * max_acceleration);
      vm_p1_p0 = 2.0 * max_velocity * (end_position - start_position);
      a_v0_v1_2 = max_acceleration * end_velocity;
      vsqr_avg = 2.0 * start_velocity * start_velocity - end_velocity * end_velocity;
      descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + a_v0_v1_2 - vm_p1_p0);
      if(descrim > 0.0) {
        peak_velocity = 0.5 * (-v0_v1_ts_sa + sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity <= start_velocity) && (peak_velocity >= end_velocity)) {
          return;
        };
        peak_velocity = 0.5 * (-v0_v1_ts_sa - sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity <= start_velocity) && (peak_velocity >= end_velocity)) {
          return;
        };
        // then, the solution doesn't fit the assumption.
      };
    };
    
    // try to assume (the last case) that vp is less in the direction (p1-p0) and both v0 and v1.
    v0_v1_ts_sa = (-2.0 * start_velocity - end_velocity + (delta_time - 0.5 * max_acceleration) * sign_p1_p0) * 2.0 / 3.0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position) * 2.0 / 3.0;
    a_v0_v1_2 = max_acceleration * end_velocity * sign_p1_p0 / 3.0;
    vsqr_avg = (2.0 * start_velocity * start_velocity + end_velocity * end_velocity) / 3.0;
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + a_v0_v1_2 - vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (-v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity <= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (-v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity <= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    
    // try to assume that fabs(a1) = fabs(vp - v1) and fabs(a0) = a_max.
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max
    //           + vp/v_max (dt - fabs(vp - v0) - 2*fabs(vp - v1) - a_max)
    //           + fabs(vp - v1) (vp + v1) / v_max 
    
    // try to assume that vp more in the direction (p1-p0) than both v1 and v0 (i.e., ramp-up and ramp-down).
    v0_v1_ts_sa = (2.0 * end_velocity + start_velocity + (delta_time - 0.5 * max_acceleration) * sign_p1_p0) * 2.0 / 3.0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position) * 2.0 / 3.0;
    a_v0_v1_2 = max_acceleration * start_velocity * sign_p1_p0 / 3.0;
    vsqr_avg = (2.0 * end_velocity * end_velocity + start_velocity * start_velocity) / 3.0;
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - a_v0_v1_2 + vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // try to assume that vp is somewhere between v1 and v0 (i.e., ramp-up ramp-up or ramp-down ramp-down).
    if(end_velocity < start_velocity) {
      // ramp-down ramp-down
      v0_v1_ts_sa = 2.0 * (2.0 * end_velocity - start_velocity + delta_time - 0.5 * max_acceleration);
      vm_p1_p0 = 2.0 * max_velocity * (end_position - start_position);
      a_v0_v1_2 = max_acceleration * start_velocity;
      vsqr_avg = 2.0 * end_velocity * end_velocity - start_velocity * start_velocity;
      descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - a_v0_v1_2 + vm_p1_p0);
      if(descrim > 0.0) {
        peak_velocity = (v0_v1_ts_sa + sqrt(descrim)) * 0.5;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity)) {
          return;
        };
        peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity)) {
          return;
        };
        // then, the solution doesn't fit the assumption.
      };
    } else {
      // ramp-up ramp-up
      v0_v1_ts_sa = 2.0 * (start_velocity - 2.0 * end_velocity + delta_time - 0.5 * max_acceleration);
      vm_p1_p0 = 2.0 * max_velocity * (end_position - start_position);
      a_v0_v1_2 = max_acceleration * start_velocity;
      vsqr_avg = 2.0 * end_velocity * end_velocity - start_velocity * start_velocity;
      descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + a_v0_v1_2 - vm_p1_p0);
      if(descrim > 0.0) {
        peak_velocity = 0.5 * (-v0_v1_ts_sa + sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity)) {
          return;
        };
        peak_velocity = 0.5 * (-v0_v1_ts_sa - sqrt(descrim));
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity)) {
          return;
        };
        // then, the solution doesn't fit the assumption.
      };
    };
    
    // try to assume (the last case) that vp is less in the direction (p1-p0) and both v0 and v1.
    v0_v1_ts_sa = (-2.0 * end_velocity - start_velocity + (delta_time - 0.5 * max_acceleration) * sign_p1_p0) * 2.0 / 3.0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position) * 2.0 / 3.0;
    a_v0_v1_2 = max_acceleration * start_velocity * sign_p1_p0 / 3.0;
    vsqr_avg = (2.0 * end_velocity * end_velocity + start_velocity * start_velocity) / 3.0;
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + a_v0_v1_2 - vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (-v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity <= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (-v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity <= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    
    // try to assume that fabs(a0) = fabs(vp - v0) and fabs(a1) = fabs(vp - v1).
    // giving:
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max 
    //           + vp/v_max (dt - 2*fabs(vp - v0) - 2*fabs(vp - v1))
    //           + fabs(vp - v1) (vp + v1) / v_max 
    
    // try to assume that vp more in the direction (p1-p0) than both v1 and v0 (i.e., ramp-up and ramp-down).
    v0_v1_ts_sa = start_velocity + end_velocity + 0.5 * delta_time * sign_p1_p0;
    vm_p1_p0 = 0.5 * max_velocity * fabs(end_position - start_position);
    vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg + vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= start_velocity * sign_p1_p0) && (sign_p1_p0 * peak_velocity >= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // try to assume that vp is somewhere between v1 and v0 (i.e., ramp-up ramp-up or ramp-down ramp-down).
    if(end_velocity > start_velocity) {
      // ramp-up ramp-up
      v0_v1_ts_sa = 2.0 * start_velocity - 2.0 * end_velocity + delta_time;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = start_velocity * start_velocity - end_velocity * end_velocity;
      if(fabs(v0_v1_ts_sa) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg) / v0_v1_ts_sa;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= start_velocity) && (peak_velocity <= end_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    } else {
      // ramp-down ramp-down
      v0_v1_ts_sa = 2.0 * end_velocity - 2.0 * start_velocity + delta_time;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = end_velocity * end_velocity - start_velocity * start_velocity;
      if(fabs(v0_v1_ts_sa) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg) / v0_v1_ts_sa;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    };
    
    // try to assume (the last case) that vp is less in the direction (p1-p0) and both v0 and v1.
    v0_v1_ts_sa = start_velocity + end_velocity - 0.5 * delta_time * sign_p1_p0;
    vm_p1_p0 = 0.5 * max_velocity * fabs(end_position - start_position);
    vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    descrim = v0_v1_ts_sa * v0_v1_ts_sa - 4.0 * (vsqr_avg - vm_p1_p0);
    if(descrim > 0.0) {
      // then there is a real root to the quadratic equation.
      peak_velocity = 0.5 * (v0_v1_ts_sa + sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity * sign_p1_p0 <= start_velocity * sign_p1_p0) && (peak_velocity * sign_p1_p0 <= end_velocity * sign_p1_p0)) {
        return;
      };
      peak_velocity = 0.5 * (v0_v1_ts_sa - sqrt(descrim));
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity * sign_p1_p0 <= start_velocity * sign_p1_p0) && (peak_velocity * sign_p1_p0 <= end_velocity * sign_p1_p0)) {
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // What the fuck!! This point should never be reached, unless the motion is completely impossible:
    peak_velocity = 0.0;
    return;
  };  
  
  template <typename PointType, typename DiffSpace, typename TimeSpace>
  double sap_compute_Ndof_interpolation_data_impl(const PointType& start_point, const PointType& end_point,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,0>::type >::point_difference_type& delta_first_order,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type& peak_velocity,
                                                  const DiffSpace& space, const TimeSpace& t_space,
                                                  double delta_time = 0.0,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type* best_peak_velocity = NULL) {
    using std::fabs;
    
    delta_first_order = get_space<0>(space,t_space).difference( get<0>(end_point), get<0>(start_point) );
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 1>::type >::point_type max_velocity = get_space<1>(space,t_space).get_upper_corner();
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace, 2>::type >::point_type max_acceleration = get_space<2>(space,t_space).get_upper_corner();
    peak_velocity = max_velocity;
    double min_dt_final = 0.0;
    
    for(std::size_t i = 0; i < delta_first_order.size(); ++i) {
      
      double norm_delta = fabs(delta_first_order[i]);
      double dp = 0.0;
      double vp = 0.0;
      
      double min_delta_time = sap_Ndof_compute_min_delta_time(
        get<0>(start_point)[i], get<0>(end_point)[i], 
        get<1>(start_point)[i], get<1>(end_point)[i], 
        dp, vp, max_velocity[i], max_acceleration[i], norm_delta);
      delta_first_order[i] = dp;
      peak_velocity[i]     = vp;
      
      if(min_dt_final < min_delta_time)
        min_dt_final = min_delta_time;
    };

    if(best_peak_velocity)
      *best_peak_velocity = peak_velocity;
    
    if(min_dt_final > delta_time)
      return min_dt_final;
    
    for(std::size_t i = 0; i < peak_velocity.size(); ++i) {
      double vp = 0.0;
      sap_Ndof_compute_peak_velocity(
        get<0>(start_point)[i], get<0>(end_point)[i], 
        get<1>(start_point)[i], get<1>(end_point)[i], 
        vp, max_velocity[i], max_acceleration[i], delta_time);
      peak_velocity[i] = vp;
      
    };
    
    return min_dt_final;
  };
  
  
  
  
  
};


};

};

#endif









