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

#include "path_planning/tangent_bundle_concept.hpp"

#include <cmath>

namespace ReaK {

namespace pp {
  
namespace detail {
  
  
  double svp_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& peak_velocity, double max_velocity);
  
  void svp_Ndof_compute_peak_velocity(double start_position, double end_position,
                                      double start_velocity, double end_velocity,
                                      double& peak_velocity, double max_velocity, double delta_time);  
  
  
  
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
  
  
  
  template <typename Idx, typename PointType, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type svp_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                         const PointType1& peak_velocity,
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
  
  template <typename Idx, typename PointType, typename PointType1, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type svp_Ndof_interpolate_impl(PointType& result, const PointType& start_point, const PointType& end_point, 
                                         const PointType1& peak_velocity, 
                                         const DiffSpace& space, const TimeSpace& t_space,
                                         double dt, double dt_total) {
    svp_Ndof_interpolate_impl< typename boost::mpl::prior<Idx>::type >(result,start_point,end_point,peak_velocity,space,t_space,dt,dt_total);
    
    get< Idx::type::value >(result) = get_space< Idx::type::value >(space,t_space).origin();
  };
  
  
  template <typename PointType, typename DiffSpace, typename TimeSpace>
  double svp_compute_Ndof_interpolation_data_impl(const PointType& start_point, const PointType& end_point,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type& peak_velocity,
                                                  const DiffSpace& space, const TimeSpace& t_space,
                                                  double delta_time = 0.0,
                                                  typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type* best_peak_velocity = NULL) {
    typename topology_traits< typename derived_N_order_space< DiffSpace, TimeSpace,1>::type >::point_type max_velocity = get_space<1>(space,t_space).get_upper_corner();
    peak_velocity = max_velocity;
    double min_dt_final = 0.0;
    
    for(std::size_t i = 0; i < peak_velocity.size(); ++i) {
      double vp = 0.0;
      double min_delta_time = svp_Ndof_compute_min_delta_time(
        get<0>(start_point)[i], get<0>(end_point)[i], 
        get<1>(start_point)[i], get<1>(end_point)[i], 
        vp, max_velocity[i]);
      peak_velocity[i] = vp;
      
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









