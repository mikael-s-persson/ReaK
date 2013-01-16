/**
 * \file sustained_velocity_pulse_Ndof_detail.hpp
 * 
 * This library contains the implementation details of the rate-limited sustained velocity 
 * pulse (SVP) interpolation.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_SUSTAINED_VELOCITY_PULSE_NDOF_DETAIL_HPP
#define REAK_SUSTAINED_VELOCITY_PULSE_NDOF_DETAIL_HPP

#include "lin_alg/arithmetic_tuple.hpp"

#include "lin_alg/vect_alg.hpp"
#include "path_planning/tangent_bundle_concept.hpp"

#include <cmath>

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<2>
    >,
  void >::type svp_Ndof_constant_accel_motion_HOT_impl(PointType&, double, std::size_t, const DiffSpace&, const TimeSpace& t_space) {
    /* Nothing to do. */ 
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline
  typename boost::enable_if<
    boost::mpl::equal_to<
      Idx,
      boost::mpl::size_t<2>
    >,
  void >::type svp_Ndof_constant_accel_motion_HOT_impl(PointType& result, double descended_accel, std::size_t i, const DiffSpace& space, const TimeSpace& t_space) {
    typedef typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type PointType2;
    const PointType2& max_acceleration = get_space<2>(space,t_space).get_upper_corner();
    get<2>(result)[i] = descended_accel * max_acceleration[i];
  };
  
  
  template <typename Idx, typename PointType>
  inline
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<2>
    >,
  void >::type svp_Ndof_constant_vel_motion_HOT_impl(PointType&, std::size_t) {
    /* Nothing to do. */ 
  };
  
  template <typename Idx, typename PointType>
  inline
  typename boost::enable_if<
    boost::mpl::equal_to<
      Idx,
      boost::mpl::size_t<2>
    >,
  void >::type svp_Ndof_constant_vel_motion_HOT_impl(PointType& result, std::size_t i) {
    get<2>(result)[i] = 0.0;
  };
  
  
  
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type svp_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                    const PointDiff0& delta_first_order, const PointType1& peak_velocity,
                                    const DiffSpace& space, const TimeSpace& t_space,
                                    double dt, double dt_total) {
    using std::fabs;
    typedef typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_difference_type PointDiff1;
    
    PointType1 max_velocity = get_space<1>(space,t_space).get_upper_corner();
    
    PointDiff1 dv1 = get_space<1>(space,t_space).difference(peak_velocity, get<1>(start_point));
    PointDiff1 dv2 = get_space<1>(space,t_space).difference(get<1>(end_point), peak_velocity);
    result = start_point;
    
    for(std::size_t i = 0; i < dv1.size(); ++i) {
      
      double dt1 = fabs(dv1[i]);
      double dt2 = fabs(dv2[i]);
      double dt_total_tmp = dt_total - dt1 - dt2;
      double dt_tmp = dt;
      
      //Phase 1: constant acceleration to the peak-velocity:
      
      if(dt1 > std::numeric_limits<double>::epsilon()) {
        if(dt_tmp > dt1)
          dt_tmp = dt1;
        
        get<0>(result)[i] += ( get<1>(result)[i] + 0.5 * dt_tmp / dt1 * dv1[i] ) * dt_tmp / max_velocity[i];
        
        get<1>(result)[i] += dt_tmp / dt1 * dv1[i];
        
        if(Idx::type::value > 1) 
          svp_Ndof_constant_accel_motion_HOT_impl<Idx>(result, dv1[i] / dt1, i, space, t_space);
        
      };
      dt_tmp = dt - dt1;
      if(dt_tmp < 0.0)
        continue;
      
      //Phase 2: constant velocity (or cruise phase):
      if(dt_tmp > dt_total_tmp)
        dt_tmp = dt_total_tmp;
      
      get<0>(result)[i] += dt_tmp * peak_velocity[i] / max_velocity[i];
      get<1>(result)[i] = peak_velocity[i];
      
      if(Idx::type::value > 1) 
        svp_Ndof_constant_vel_motion_HOT_impl<Idx>(result,i);
      
      dt_tmp = dt - dt1 - dt_total_tmp;
      if(dt_tmp < 0.0)
        continue;
      
      //Phase 3: constant acceleration to end-velocity:
      
      if(dt2 > std::numeric_limits<double>::epsilon()) {
        if(dt_tmp > dt2)
          dt_tmp = dt2;
        
        get<0>(result)[i] += ( get<1>(result)[i] + 0.5 * dt_tmp / dt2 * dv2[i] ) * dt_tmp / max_velocity[i];
        
        get<1>(result)[i] += dt_tmp / dt2 * dv2[i];
        
        if(Idx::type::value > 1) 
          svp_Ndof_constant_accel_motion_HOT_impl<Idx>(result, dv2[i] / dt2, i, space, t_space);
        
      };
    };
    
  };
  
  template <typename Idx, typename PointType, typename PointDiff0, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type svp_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                    const PointDiff0& delta_first_order, const PointType1& peak_velocity, 
                                    const DiffSpace& space, const TimeSpace& t_space,
                                    double dt, double dt_total) {
    svp_Ndof_interpolate_impl< typename boost::mpl::prior<Idx>::type >(result,start_point,end_point,delta_first_order,peak_velocity,space,t_space,dt,dt_total);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
  
  
  
  inline
  double svp_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& delta_first_order, double& peak_velocity, 
                                         double max_velocity, double& norm_delta) {
    using std::fabs;
    using std::sqrt;
    
    // try to assume that peak_velocity = sign(p1 - p0) * max_velocity
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    peak_velocity = sign_p1_p0 * max_velocity;
    double descended_peak_velocity = peak_velocity / max_velocity;
    delta_first_order = end_position - start_position
      - (0.5 * fabs(peak_velocity -   end_velocity)) * (descended_peak_velocity +   end_velocity / max_velocity)
      - (0.5 * fabs(peak_velocity - start_velocity)) * (descended_peak_velocity + start_velocity / max_velocity);
    norm_delta = fabs(delta_first_order);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-position):
      return norm_delta + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity);
    };
    // if not, then can try to see if we simply can quite reach max velocity before having to ramp-down:
    delta_first_order = 0.0; norm_delta = 0.0;
    // this assumes that we have p1 - p0 == 0.5 / vm * ( fabs(vp - v1) * (vp + v1) + fabs(vp - v0) * (vp + v0) )
    // first try if vp is more in the direction (p1-p0) than both v1 and v0:
    peak_velocity = sqrt(max_velocity * fabs(end_position - start_position) + 0.5 * start_velocity * start_velocity + 0.5 * end_velocity * end_velocity);
    if( ( peak_velocity > sign_p1_p0 * start_velocity ) && ( peak_velocity > sign_p1_p0 * end_velocity ) ) {
      // this means that the vp solution is consistent with the assumption:
      peak_velocity *= sign_p1_p0;
      return fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity);
    };
    // else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible):
    if( max_velocity * fabs(end_position - start_position) < 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity) ) {
      // this means there exists a solution for the magnitude of vp:
      peak_velocity = sqrt(0.5 * start_velocity * start_velocity + 0.5 * end_velocity * end_velocity - max_velocity * fabs(end_position - start_position));
      if( ( peak_velocity < sign_p1_p0 * start_velocity ) && ( peak_velocity < sign_p1_p0 * end_velocity ) ) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        peak_velocity *= sign_p1_p0;
      } else {
        // this means there must be a valid solution for when vp is in the opposite direction of (p1-p0):
        peak_velocity *= -sign_p1_p0;
      };
      return fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity);
    };
    // What the fuck!! This point should never be reached, unless the motion is completely impossible:
    peak_velocity = 0.0;
    return std::numeric_limits<double>::infinity();
  };
  
  inline
  void svp_Ndof_compute_peak_velocity(double start_position, double end_position,
                                      double start_velocity, double end_velocity,
                                      double& peak_velocity, double max_velocity, double delta_time) {
    // NOTE: Assume that delta-time is larger than minimum reachable delta-time 
    //       (avoid checking for that, even if it could mean that vp is higher than maximum)
    
    using std::fabs;
    using std::sqrt;
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    
    // try to assume that vp more in the direction (p1-p0) than both v1 and v0 (i.e., ramp-up and ramp-down).
    double v0_v1_ts = start_velocity + end_velocity + delta_time * sign_p1_p0;
    double vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    double vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    if(v0_v1_ts * v0_v1_ts > 4.0 * (vm_p1_p0 + vsqr_avg)) {
      // then there is a real root to the quadratic equation.
      double r1 = 0.5 * (v0_v1_ts + sqrt(v0_v1_ts * v0_v1_ts - 4.0 * (vm_p1_p0 + vsqr_avg))) * sign_p1_p0;
      if((fabs(r1) <= max_velocity) && (r1 >= start_velocity * sign_p1_p0) && (r1 >= end_velocity * sign_p1_p0)) {
        peak_velocity = sign_p1_p0 * r1;
        return;
      };
      r1 = 0.5 * (v0_v1_ts - sqrt(v0_v1_ts * v0_v1_ts - 4.0 * (vm_p1_p0 + vsqr_avg))) * sign_p1_p0;
      if((fabs(r1) <= max_velocity) && (r1 >= start_velocity * sign_p1_p0) && (r1 >= end_velocity * sign_p1_p0)) {
        peak_velocity = sign_p1_p0 * r1;
        return;
      };
      // then, there was no suitable root in this region.
    };
    
    // try to assume that vp is somewhere between v1 and v0 (i.e., ramp-up ramp-up or ramp-down ramp-down).
    if(end_velocity > start_velocity) {
      // ramp-up ramp-up
      v0_v1_ts = start_velocity - end_velocity + delta_time;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = 0.5 * (start_velocity * start_velocity - end_velocity * end_velocity);
      if(fabs(v0_v1_ts) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg) / v0_v1_ts;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= start_velocity) && (peak_velocity <= end_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    } else {
      // ramp-down ramp-down
      v0_v1_ts = end_velocity - start_velocity + delta_time;
      vm_p1_p0 = max_velocity * (end_position - start_position);
      vsqr_avg = 0.5 * (end_velocity * end_velocity - start_velocity * start_velocity);
      if(fabs(v0_v1_ts) > 1e-6 * max_velocity) {
        peak_velocity = (vm_p1_p0 + vsqr_avg) / v0_v1_ts;
        if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= end_velocity) && (peak_velocity <= start_velocity))
          return;
        // then, the solution doesn't fit the assumption.
      };
    };
    
    // try to assume (the last case) that vp is less in the direction (p1-p0) and both v0 and v1.
    v0_v1_ts = start_velocity + end_velocity - delta_time * sign_p1_p0;
    vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    if(v0_v1_ts * v0_v1_ts > 4.0 * (vsqr_avg - vm_p1_p0)) {
      // then there is a real root to the quadratic equation.
      double r1 = 0.5 * (v0_v1_ts + sqrt(v0_v1_ts * v0_v1_ts - 4.0 * (vsqr_avg - vm_p1_p0))) * sign_p1_p0;
      if((fabs(r1) <= max_velocity) && (r1 <= start_velocity * sign_p1_p0) && (r1 <= end_velocity * sign_p1_p0)) {
        peak_velocity = sign_p1_p0 * r1;
        return;
      };
      r1 = 0.5 * (v0_v1_ts - sqrt(v0_v1_ts * v0_v1_ts - 4.0 * (vsqr_avg - vm_p1_p0))) * sign_p1_p0;
      if((fabs(r1) <= max_velocity) && (r1 <= start_velocity * sign_p1_p0) && (r1 <= end_velocity * sign_p1_p0)) {
        peak_velocity = sign_p1_p0 * r1;
        return;
      };
      // then, there was no suitable root in this region.
    };
    // What the fuck!! This point should never be reached, unless the motion is completely impossible:
    peak_velocity = 0.0;
    return;
  };  
  
  template <typename PointType, typename DiffSpace, typename TimeSpace>
  double svp_compute_Ndof_interpolation_data_impl(const PointType& start_point, const PointType& end_point,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,0>::type >::point_difference_type& delta_first_order,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type& peak_velocity,
                                                  const DiffSpace& space, const TimeSpace& t_space,
                                                  double delta_time = 0.0,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type* best_peak_velocity = NULL) {
    using std::fabs;
    
    delta_first_order = get_space<0>(space,t_space).difference( get<0>(end_point), get<0>(start_point) );
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type max_velocity = get_space<1>(space,t_space).get_upper_corner();
    peak_velocity = max_velocity;
    double min_dt_final = 0.0;
    
    for(std::size_t i = 0; i < delta_first_order.size(); ++i) {
      
      double norm_delta = fabs(delta_first_order[i]);
      double dp = 0.0;
      double vp = 0.0;
      
      double min_delta_time = svp_Ndof_compute_min_delta_time(
        get<0>(start_point)[i], get<0>(end_point)[i], 
        get<1>(start_point)[i], get<1>(end_point)[i], 
        dp, vp, max_velocity[i], norm_delta);
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
      svp_Ndof_compute_peak_velocity(
        get<0>(start_point)[i], get<0>(end_point)[i], 
        get<1>(start_point)[i], get<1>(end_point)[i], 
        vp, max_velocity[i], delta_time);
      peak_velocity[i] = vp;
      
    };
    
    return min_dt_final;
  };
  
  
  
  
  
};


};

};

#endif









