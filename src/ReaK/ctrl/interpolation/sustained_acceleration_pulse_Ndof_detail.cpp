
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

#include "sustained_acceleration_pulse_Ndof_detail.hpp"

#include <cmath>


#include "root_finders/bisection_method.hpp"
#include "root_finders/secant_method.hpp"


// #define RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION


#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION

#include <fstream>

#endif

namespace ReaK {

namespace pp {
  
  
namespace detail {

  
  
  
  static void sap_Ndof_compute_interpolated_values_closedform(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double peak_velocity, double max_velocity, double max_acceleration,
    double dt, double dt_total,
    double& result_pos, double& result_vel, 
    double& result_acc, double& result_desc_jerk) {
    
    using std::fabs;
    using std::sqrt;
    
    double dv1 = peak_velocity - start_velocity;
    double dv2 = end_velocity - peak_velocity;
    result_pos = start_position;
    result_vel = start_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    double dt_vp1_1st = fabs(dv1);
    double sgn_vp1 = 1.0;
    if( dv1 < 0.0 ) sgn_vp1 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp1 = dt_vp1_1st - max_acceleration;
    double dt_ap1 = max_acceleration;
    if( dt_vp1 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp1 = 0.0;
      dt_ap1 = sqrt(max_acceleration * dt_vp1_1st);
    };
    
    double dt_vp2_1st = fabs(dv2);
    double sgn_vp2 = 1.0;
    if( dv2 < 0.0 ) sgn_vp2 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - max_acceleration;
    double dt_ap2 = max_acceleration;
    if( dt_vp2 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt(max_acceleration * dt_vp2_1st);
    };
    
    
    dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp1 + 2.0 * dt_ap1;
    
    
    if( dt < (2.0 * dt_ap1 + dt_vp1) ) {
      
      if( dt < dt_ap1 ) {
        
        //Segment 1: in the jerk-up phase of velocity ramp-up.
        // assume result_acc == 0
        result_pos += ( result_vel + ( dt * sgn_vp1 / 6.0 ) * dt / max_acceleration ) * dt / max_velocity;
        result_vel += ( 0.5 * dt * sgn_vp1 ) * dt / max_acceleration;
        result_acc  = dt * sgn_vp1;
        result_desc_jerk = sgn_vp1;
        
      } else if( dt < dt_ap1 + dt_vp1 ) {
        //Segment 2: in the constant accel phase of velocity ramp-up.
        
        dt -= dt_ap1;
        
        result_pos += ( ( result_vel + 0.5 * dt * sgn_vp1) * ( max_acceleration + dt ) 
                      + max_acceleration * max_acceleration * sgn_vp1 / 6.0 ) / max_velocity;
        result_vel += (dt + 0.5 * max_acceleration) * sgn_vp1;
        result_acc  = max_acceleration * sgn_vp1;
        result_desc_jerk = 0.0;
        
      } else {
        //Segment 3: in the jerk-down phase of velocity ramp-up.
        
        dt -= dt_ap1 + dt_vp1;
        
        result_pos += ( result_vel * ( dt_ap1 + dt_vp1 + dt )
          + ( dt_ap1 * dt * ( 0.5 * dt + dt_vp1 + 0.5 * dt_ap1 )
            + 0.5 * dt_ap1 * ( dt_ap1 * dt_ap1 / 3.0 + dt_vp1 * dt_ap1 + dt_vp1 * dt_vp1 ) 
            - dt * dt * dt / 6.0 ) * sgn_vp1 / max_acceleration ) / max_velocity;
        result_vel += ( dt_ap1 * ( dt + dt_vp1 + 0.5 * dt_ap1 )
                      - 0.5 * dt * dt ) * sgn_vp1 / max_acceleration;
        result_acc  = ( dt_ap1 - dt ) * sgn_vp1;
        result_desc_jerk = -sgn_vp1;
        
        
        // // alternative calculation (from the end of segment, backwards)
        // double mdt -= 2.0 * dt_ap1 + dt_vp1 - dt;
        
        // result_pos += ( peak_velocity * ( 2.0 * dt_ap1 + dt_vp1 - mdt )
        //   + ( mdt * mdt * mdt / 6.0 - dt_ap1 * ( dt_ap1 + dt_vp1 ) * ( dt_ap1 + 0.5 * dt_vp1 ) ) * sgn_vp1 / max_acceleration ) / max_velocity;
        // // peak_velocity == result_vel + ( dt_ap1 * dt_vp1 + dt_ap1 * dt_ap1 ) * sgn_vp1 / max_acceleration
        // result_vel  = peak_velocity - ( 0.5 * mdt * sgn_vp1 ) * mdt / max_acceleration;
        // result_acc  = mdt * sgn_vp1;
        // result_desc_jerk = -sgn_vp1;
        
      };
      
      return;
      
    };
    
    result_pos += ( peak_velocity * ( 2.0 * dt_ap1 + dt_vp1 )
        - dt_ap1 * ( dt_ap1 + dt_vp1 ) * ( dt_ap1 + 0.5 * dt_vp1 ) * sgn_vp1 / max_acceleration ) / max_velocity;
    result_vel  = peak_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    if( dt < (2.0 * dt_ap1 + dt_vp1 + dt_total) ) {
      //Segment 4: in the cruise phase.
      
      dt -= 2.0 * dt_ap1 + dt_vp1;
      
      result_pos += dt * peak_velocity / max_velocity;
      
      return;
    };
    
    result_pos += dt_total * peak_velocity / max_velocity;
    
    
    if( dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + 2.0 * dt_ap2 + dt_vp2) ) {
      
      if( dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2) ) {
        
        //Segment 5: in the jerk-up phase of velocity ramp-down.
        
        dt -= 2.0 * dt_ap1 + dt_vp1 + dt_total;
        
        result_pos += ( result_vel + ( dt * sgn_vp2 / 6.0 ) * dt / max_acceleration ) * dt / max_velocity;
        result_vel += ( 0.5 * dt * sgn_vp2 ) * dt / max_acceleration;
        result_acc  = dt * sgn_vp2;
        result_desc_jerk = sgn_vp2;
        
      } else if( dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2 + dt_vp2) ) {
        
        //Segment 6: in the constant accel phase of velocity ramp-down.
        
        dt -= 2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2;
        
        result_pos += ( result_vel * ( max_acceleration + dt )
          + ( max_acceleration * max_acceleration / 6.0 + 0.5 * max_acceleration * dt + 0.5 * dt * dt ) * sgn_vp2 ) / max_velocity;
        result_vel += ( dt + 0.5 * max_acceleration ) * sgn_vp2;
        result_acc  = max_acceleration * sgn_vp2;
        result_desc_jerk = 0.0;
        
      } else {
        
        //Segment 7: in the jerk-down phase of velocity ramp-down.
        
        dt -= 2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2 + dt_vp2;
        
        result_pos += ( result_vel * ( dt_ap2 + dt_vp2 + dt )
          + ( dt_ap2 * dt_ap2 * dt_ap2 / 6.0 
            + 0.5 * dt_vp2 * dt_ap2 * dt_ap2 
            + 0.5 * dt_vp2 * dt_vp2 * dt_ap2 
            + dt * dt_vp2 * dt_ap2 
            + 0.5 * dt * dt_ap2 * dt_ap2 
            + 0.5 * dt * dt * dt_ap2 
            - dt * dt * dt / 6.0 ) * sgn_vp2 / max_acceleration ) / max_velocity;
        result_vel += ( dt_ap2 * sgn_vp2 - 0.5 * dt * sgn_vp2 ) * dt / max_acceleration + ( dt_vp2 + 0.5 * dt_ap2 ) * dt_ap2 * sgn_vp2 / max_acceleration;
        result_acc  = ( dt_ap2 - dt ) * sgn_vp2;
        result_desc_jerk = -sgn_vp2;
        
        // // alternative calculation (from the end of segment, backwards)
        // double mdt -= 2.0 * dt_ap1 + dt_vp1 + dt_total + 2.0 * dt_ap1 + dt_vp1 - dt;
        
        // result_pos += ( end_velocity * ( 2.0 * dt_ap2 + dt_vp2 - mdt )
        //   + ( mdt * mdt * mdt / 6.0 - dt_ap2 * ( dt_ap2 + dt_vp2 ) * ( dt_ap2 + 0.5 * dt_vp2 ) ) * sgn_vp2 / max_acceleration ) / max_velocity;
        // // end_velocity == result_vel + ( dt_ap2 * dt_vp2 + dt_ap2 * dt_ap2 ) * sgn_vp2 / max_acceleration
        // result_vel  = end_velocity - ( 0.5 * mdt * sgn_vp2 ) * mdt / max_acceleration;
        // result_acc  = mdt * sgn_vp2;
        // result_desc_jerk = -sgn_vp2;
        
      };
      
      return;
    };
    
    
//     result_pos += ( result_vel + 0.5 * ( dt_vp2 + dt_ap2 ) * dt_ap2 * sgn_vp2 / max_acceleration ) * ( 2.0 * dt_ap2 + dt_vp2 ) / max_velocity;
    result_vel += ( dt_ap2 + dt_vp2 ) * dt_ap2 * sgn_vp2 / max_acceleration;
    result_acc  = ( dt_ap2 - dt_ap2 ) * sgn_vp2;
    result_desc_jerk = -sgn_vp2;
    
//     result_pos = end_position;
//     result_vel = end_velocity;
//     result_acc = 0.0;
//     result_desc_jerk = 0.0;
    
    
    // resulting position formula:
    
    
//     result_pos += ( 
//         peak_velocity * ( dt_total + 2.0 * dt_ap1 + dt_vp1 + 2.0 * dt_ap2 + dt_vp2 )
//       - 0.5 * ( dt_ap1 + dt_vp1 ) * ( 2.0 * dt_ap1 + dt_vp1 ) * dt_ap1 * sgn_vp1 / max_acceleration
//       + 0.5 * ( dt_vp2 + dt_ap2 ) * ( 2.0 * dt_ap2 + dt_vp2 ) * dt_ap2 * sgn_vp2 / max_acceleration ) / max_velocity;
    
    result_pos = start_position + ( peak_velocity * dt_total 
      + 0.5 * ( start_velocity + peak_velocity ) * ( 2.0 * dt_ap1 + dt_vp1 ) 
      + 0.5 * ( peak_velocity  + end_velocity  ) * ( 2.0 * dt_ap2 + dt_vp2 ) ) / max_velocity;
    
      
  };
  
  
  
  
#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION
  
  
  
  static void sap_Ndof_compute_interpolated_values_incremental_logged_version(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double peak_velocity, double max_velocity, double max_acceleration,
    double dt, double dt_total,
    double& result_pos, double& result_vel, 
    double& result_acc, double& result_desc_jerk,
    std::ostream& log_output) {
    
    using std::fabs;
    using std::sqrt;
    
    double original_dt = dt;
    double original_dt_total = dt_total;
    double current_dt = 0.0;
    
    double dv1 = peak_velocity - start_velocity;
    double dv2 = end_velocity - peak_velocity;
    result_pos = start_position;
    result_vel = start_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    double dt_vp1_1st = fabs(dv1);
    double sgn_vp1 = 1.0;
    if( dv1 < 0.0 ) sgn_vp1 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp1 = dt_vp1_1st - max_acceleration;
    double dt_ap1 = max_acceleration;
    if( dt_vp1 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp1 = 0.0;
      dt_ap1 = sqrt(max_acceleration * dt_vp1_1st);
    };
    
    double dt_vp2_1st = fabs(dv2);
    double sgn_vp2 = 1.0;
    if( dv2 < 0.0 ) sgn_vp2 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - max_acceleration;
    double dt_ap2 = max_acceleration;
    if( dt_vp2 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt(max_acceleration * dt_vp2_1st);
    };
    
    
    dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp1 + 2.0 * dt_ap1;
    
    
    log_output << "-----------------------------------------------------------\n"
               << " dt_vp1 = " << dt_vp1 << "\n"
               << " dt_ap1 = " << dt_ap1 << "\n"
               << " dt_vp2 = " << dt_vp2 << "\n"
               << " dt_ap2 = " << dt_ap2 << "\n"
               << " dt_total = " << dt_total << std::endl;
    
    
    if(dt_ap1 > std::numeric_limits<double>::epsilon()) {
      //Phase 1: in the jerk-up phase of velocity ramp-up.
      double descended_jerk = sgn_vp1;
      if( dt < dt_ap1 )
        dt_ap1 = dt;
      
      log_output << "-----------------------------------------------------------\n"
                 << " Phase 1: in the jerk-up phase of velocity ramp-up.\n"
                 << " dt = " << dt_ap1 << "\n"
                 << " delta_pos = " << ( ( result_vel + ( dt_ap1 * descended_jerk / 6.0 ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity ) << "\n"
                 << " delta_vel = " << ( ( 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration ) << "\n"
                 << " delta_acc = " << ( dt_ap1 * descended_jerk ) << std::endl;
      
      {
        double exp_delta_pos = 0.0;
        for(double exp_t = 0.0; exp_t < 1.00005 * dt_ap1; exp_t += 0.0001 * dt_ap1) {
          double exp_dv = ( 0.5 * exp_t * descended_jerk ) * exp_t / max_acceleration;
          exp_delta_pos += 0.0001 * dt_ap1 * ( result_vel + exp_dv ) / max_velocity;
        };
        log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
      };
      
      current_dt += dt_ap1;
      {
        double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
        
        sap_Ndof_compute_interpolated_values_closedform(
          start_position, end_position, start_velocity, end_velocity,
          peak_velocity, max_velocity, max_acceleration,
          current_dt - 1e-4, original_dt_total,
          cf_pos, cf_vel, cf_acc, cf_desc_jerk);
        
        log_output << " closed-form pos = " << cf_pos << "\n"
                   << " closed-form vel = " << cf_vel << "\n"
                   << " closed-form acc = " << cf_acc << std::endl;
        
      };
      
      
      // assume result_acc == 0
      result_pos += ( result_vel + ( dt_ap1 * descended_jerk / 6.0 ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity;
      result_vel += ( 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration;
      result_acc  = dt_ap1 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
      log_output << " pos = " << result_pos << "\n"
                 << " vel = " << result_vel << "\n"
                 << " acc = " << result_acc << std::endl;
      
      dt -= dt_ap1;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
      
      //Phase 2: in the constant accel phase of velocity ramp-up.
      if(dt_vp1 > std::numeric_limits<double>::epsilon()) {
        double descended_accel = result_acc / max_acceleration;
        if( dt < dt_vp1 )
          dt_vp1 = dt;
        
        log_output << "-----------------------------------------------------------\n"
                   << " Phase 2: in the constant accel phase of velocity ramp-up.\n"
                   << " dt = " << dt_vp1 << "\n"
                   << " delta_pos = " << ( ( result_vel + 0.5 * dt_vp1 * descended_accel ) * dt_vp1 / max_velocity ) << "\n"
                   << " delta_vel = " << ( dt_vp1 * descended_accel ) << std::endl;
        
        {
          double exp_delta_pos = 0.0;
          for(double exp_t = 0.0; exp_t < 1.00005 * dt_vp1; exp_t += 0.0001 * dt_vp1) {
            double exp_dv = ( exp_t * descended_accel );
            exp_delta_pos += 0.0001 * dt_vp1 * ( result_vel + exp_dv ) / max_velocity;
          };
          log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
        };
        
        current_dt += dt_vp1;
        {
          double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
          
          sap_Ndof_compute_interpolated_values_closedform(
            start_position, end_position, start_velocity, end_velocity,
            peak_velocity, max_velocity, max_acceleration,
            current_dt - 1e-4, original_dt_total,
            cf_pos, cf_vel, cf_acc, cf_desc_jerk);
          
          log_output << " closed-form pos = " << cf_pos << "\n"
                     << " closed-form vel = " << cf_vel << "\n"
                     << " closed-form acc = " << cf_acc << std::endl;
          
        };
        
        
        result_pos += ( result_vel + 0.5 * dt_vp1 * descended_accel ) * dt_vp1 / max_velocity;
        result_vel += dt_vp1 * descended_accel;
        result_desc_jerk = 0.0;
        
        log_output << " pos = " << result_pos << "\n"
                   << " vel = " << result_vel << std::endl;
        
        dt -= dt_vp1;
        if(dt <= std::numeric_limits<double>::epsilon())
          return;
      };
      
      //Phase 3: in the jerk-down phase of velocity ramp-up.
      descended_jerk = -sgn_vp1;
      if( dt < dt_ap1 )
        dt_ap1 = dt;
      
      log_output << "-----------------------------------------------------------\n"
                 << " Phase 3: in the jerk-down phase of velocity ramp-up.\n"
                 << " dt = " << dt_ap1 << "\n"
                 << " delta_pos = " << ( ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity ) << "\n"
                 << " delta_vel = " << ( ( result_acc + 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration ) << "\n"
                 << " delta_acc = " << ( dt_ap1 * descended_jerk ) << std::endl;
      
      {
        double exp_delta_pos = 0.0;
        for(double exp_t = 0.0; exp_t < 1.00005 * dt_ap1; exp_t += 0.0001 * dt_ap1) {
          double exp_dv = ( result_acc + 0.5 * exp_t * descended_jerk ) * exp_t / max_acceleration;
          exp_delta_pos += 0.0001 * dt_ap1 * ( result_vel + exp_dv ) / max_velocity;
        };
        log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
      };
      
      current_dt += dt_ap1;
      {
        double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
        
        sap_Ndof_compute_interpolated_values_closedform(
          start_position, end_position, start_velocity, end_velocity,
          peak_velocity, max_velocity, max_acceleration,
          current_dt - 1e-4, original_dt_total,
          cf_pos, cf_vel, cf_acc, cf_desc_jerk);
        
        log_output << " closed-form pos = " << cf_pos << "\n"
                   << " closed-form vel = " << cf_vel << "\n"
                   << " closed-form acc = " << cf_acc << std::endl;
        
      };
      
      result_pos += ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity;
      result_vel += ( result_acc + 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration;
      result_acc += dt_ap1 * descended_jerk;
      
      result_desc_jerk = descended_jerk;
      
      log_output << " pos = " << result_pos << "\n"
                 << " vel = " << result_vel << "\n"
                 << " acc = " << result_acc << std::endl;
      
      dt -= dt_ap1;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
    };
    
    //Phase 4: in the cruise phase.
    if( dt < dt_total )
      dt_total = dt;
    
    log_output << "-----------------------------------------------------------\n"
               << " Phase 4: in the cruise phase.\n"
               << " dt = " << dt_total << "\n"
               << " delta_pos = " << ( dt_total * peak_velocity / max_velocity ) << std::endl;
    
    {
      double exp_delta_pos = 0.0;
      for(double exp_t = 0.0; exp_t < 1.00005 * dt_total; exp_t += 0.0001 * dt_total) {
        exp_delta_pos += 0.0001 * dt_total * result_vel / max_velocity;
      };
      log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
    };
    
    current_dt += dt_total;
    {
      double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
      
      sap_Ndof_compute_interpolated_values_closedform(
        start_position, end_position, start_velocity, end_velocity,
        peak_velocity, max_velocity, max_acceleration,
        current_dt - 1e-4, original_dt_total,
        cf_pos, cf_vel, cf_acc, cf_desc_jerk);
      
      log_output << " closed-form pos = " << cf_pos << "\n"
                 << " closed-form vel = " << cf_vel << "\n"
                 << " closed-form acc = " << cf_acc << std::endl;
      
    };
    
    result_pos += dt_total * peak_velocity / max_velocity;
    result_vel = peak_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    log_output << " pos = " << result_pos << std::endl;
    
    dt -= dt_total;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    if(dt_ap2 > std::numeric_limits<double>::epsilon()) {
      //Phase 5: in the jerk-up phase of velocity ramp-down.
      double descended_jerk = sgn_vp2;
      if( dt < dt_ap2 )
        dt_ap2 = dt;
      
      log_output << "-----------------------------------------------------------\n"
                 << " Phase 5: in the jerk-up phase of velocity ramp-down.\n"
                 << " dt = " << dt_ap2 << "\n"
                 << " delta_pos = " << ( ( result_vel + ( dt_ap2 * descended_jerk / 6.0 ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity ) << "\n"
                 << " delta_vel = " << ( ( 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration ) << "\n"
                 << " delta_acc = " << ( dt_ap2 * descended_jerk ) << std::endl;
      
      {
        double exp_delta_pos = 0.0;
        for(double exp_t = 0.0; exp_t < 1.00005 * dt_ap2; exp_t += 0.0001 * dt_ap2) {
          double exp_dv = ( 0.5 * exp_t * descended_jerk ) * exp_t / max_acceleration;
          exp_delta_pos += 0.0001 * dt_ap2 * ( result_vel + exp_dv ) / max_velocity;
        };
        log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
      };
      
      current_dt += dt_ap2;
      {
        double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
        
        sap_Ndof_compute_interpolated_values_closedform(
          start_position, end_position, start_velocity, end_velocity,
          peak_velocity, max_velocity, max_acceleration,
          current_dt - 1e-4, original_dt_total,
          cf_pos, cf_vel, cf_acc, cf_desc_jerk);
        
        log_output << " closed-form pos = " << cf_pos << "\n"
                   << " closed-form vel = " << cf_vel << "\n"
                   << " closed-form acc = " << cf_acc << std::endl;
        
      };
      
      result_pos += ( result_vel + ( dt_ap2 * descended_jerk / 6.0 ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity;
      result_vel += ( 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration;
      result_acc  = dt_ap2 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
      log_output << " pos = " << result_pos << "\n"
                 << " vel = " << result_vel << "\n"
                 << " acc = " << result_acc << std::endl;
      
      dt -= dt_ap2;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
      
      //Phase 6: in the constant accel phase of velocity ramp-down.
      if(dt_vp2 > std::numeric_limits<double>::epsilon()) {
        double descended_accel = result_acc / max_acceleration;
        if( dt < dt_vp2 )
          dt_vp2 = dt;
        
        log_output << "-----------------------------------------------------------\n"
                   << " Phase 6: in the constant accel phase of velocity ramp-down.\n"
                   << " dt = " << dt_vp2 << "\n"
                   << " delta_pos = " << ( ( result_vel + 0.5 * dt_vp2 * descended_accel ) * dt_vp2 / max_velocity ) << "\n"
                   << " delta_vel = " << ( dt_vp2 * descended_accel ) << std::endl;
        
        {
          double exp_delta_pos = 0.0;
          for(double exp_t = 0.0; exp_t < 1.00005 * dt_vp2; exp_t += 0.0001 * dt_vp2) {
            double exp_dv = exp_t * descended_accel;
            exp_delta_pos += 0.0001 * dt_vp2 * ( result_vel + exp_dv ) / max_velocity;
          };
          log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
        };
        
        current_dt += dt_vp2;
        {
          double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
          
          sap_Ndof_compute_interpolated_values_closedform(
            start_position, end_position, start_velocity, end_velocity,
            peak_velocity, max_velocity, max_acceleration,
            current_dt - 1e-4, original_dt_total,
            cf_pos, cf_vel, cf_acc, cf_desc_jerk);
          
          log_output << " closed-form pos = " << cf_pos << "\n"
                     << " closed-form vel = " << cf_vel << "\n"
                     << " closed-form acc = " << cf_acc << std::endl;
          
        };
        
        result_pos += ( result_vel + 0.5 * dt_vp2 * descended_accel ) * dt_vp2 / max_velocity;
        result_vel += dt_vp2 * descended_accel;
        result_desc_jerk = 0.0;
        
        log_output << " pos = " << result_pos << "\n"
                   << " vel = " << result_vel << std::endl;
        
        dt -= dt_vp2;
        if(dt <= std::numeric_limits<double>::epsilon())
          return;
      };
      
      //Phase 7: in the jerk-down phase of velocity ramp-down.
      descended_jerk = -sgn_vp2;
      if( dt < dt_ap2 )
        dt_ap2 = dt;
      
      log_output << "-----------------------------------------------------------\n"
                 << " Phase 7: in the jerk-down phase of velocity ramp-down.\n"
                 << " dt = " << dt_ap2 << "\n"
                 << " delta_pos = " << ( ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity ) << "\n"
                 << " delta_vel = " << ( ( result_acc + 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration ) << "\n"
                 << " delta_acc = " << ( dt_ap2 * descended_jerk ) << std::endl;
      
      {
        double exp_delta_pos = 0.0;
        for(double exp_t = 0.0; exp_t < 1.00005 * dt_ap2; exp_t += 0.0001 * dt_ap2) {
          double exp_dv = ( result_acc + 0.5 * exp_t * descended_jerk ) * exp_t / max_acceleration;
          exp_delta_pos += 0.0001 * dt_ap2 * ( result_vel + exp_dv ) / max_velocity;
        };
        log_output << " num-integrated delta_pos = " << exp_delta_pos << std::endl;
      };
      
      current_dt += dt_ap2;
      {
        double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
        
        sap_Ndof_compute_interpolated_values_closedform(
          start_position, end_position, start_velocity, end_velocity,
          peak_velocity, max_velocity, max_acceleration,
          current_dt - 1e-4, original_dt_total,
          cf_pos, cf_vel, cf_acc, cf_desc_jerk);
        
        log_output << " closed-form pos = " << cf_pos << "\n"
                   << " closed-form vel = " << cf_vel << "\n"
                   << " closed-form acc = " << cf_acc << std::endl;
        
      };
      
      result_pos += ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity;
      result_vel += ( result_acc + 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration;
      result_acc += dt_ap2 * descended_jerk;
      
      result_desc_jerk = descended_jerk;
      
      log_output << " pos = " << result_pos << "\n"
                 << " vel = " << result_vel << "\n"
                 << " acc = " << result_acc << std::endl;
      
    };
  };
  
  
  
  
  
  
  void sap_Ndof_compute_interpolated_values(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double peak_velocity, double max_velocity, double max_acceleration,
    double dt, double dt_total,
    double& result_pos, double& result_vel, 
    double& result_acc, double& result_desc_jerk) {
    
    
    static std::ofstream log_output("sap_interp_log.txt");
    
    log_output << "\n\n0000000000000000000000000000000000000000000000000000000000000000000\n\n"
               << "Starting to compute the interpolation for....\n"
               << "Positions: ( " << start_position << " , " << end_position << " )\n"
               << "Velocities: ( " << start_velocity << " , " << end_velocity << " )\n"
               << "Max Velocity: " << max_velocity << "\n"
               << "Max Acceleration: " << max_acceleration << "\n"
               << "Peak Velocity: " << peak_velocity << "\n"
               << "Requested time: " << dt << "\n"
               << "over Total time: " << dt_total << std::endl;
    
    sap_Ndof_compute_interpolated_values_incremental_logged_version(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration, dt, dt_total,
      result_pos, result_vel, result_acc, result_desc_jerk, log_output);
    
    log_output << "The calculation of the minimum delta-time has finished and got the following:\n";
    
  };
  
  
#else
  
  
  static void sap_Ndof_compute_interpolated_values_incremental(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double peak_velocity, double max_velocity, double max_acceleration,
    double dt, double dt_total,
    double& result_pos, double& result_vel, 
    double& result_acc, double& result_desc_jerk) {
    
    using std::fabs;
    using std::sqrt;
    
    double dv1 = peak_velocity - start_velocity;
    double dv2 = end_velocity - peak_velocity;
    result_pos = start_position;
    result_vel = start_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    double dt_vp1_1st = fabs(dv1);
    double sgn_vp1 = 1.0;
    if( dv1 < 0.0 ) sgn_vp1 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp1 = dt_vp1_1st - max_acceleration;
    double dt_ap1 = max_acceleration;
    if( dt_vp1 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp1 = 0.0;
      dt_ap1 = sqrt(max_acceleration * dt_vp1_1st);
    };
    
    double dt_vp2_1st = fabs(dv2);
    double sgn_vp2 = 1.0;
    if( dv2 < 0.0 ) sgn_vp2 = -1.0;
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp2 = dt_vp2_1st - max_acceleration;
    double dt_ap2 = max_acceleration;
    if( dt_vp2 < 0.0 ) {
      //means that we don't have time to reach the maximum acceleration:
      dt_vp2 = 0.0;
      dt_ap2 = sqrt(max_acceleration * dt_vp2_1st);
    };
    
    
    dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp1 + 2.0 * dt_ap1;
    
    
    if(dt_ap1 > std::numeric_limits<double>::epsilon()) {
      //Phase 1: in the jerk-up phase of velocity ramp-up.
      double descended_jerk = sgn_vp1;
      if( dt < dt_ap1 )
        dt_ap1 = dt;
      
      // assume result_acc == 0
      result_pos += ( result_vel + ( dt_ap1 * descended_jerk / 6.0 ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity;
      result_vel += ( 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration;
      result_acc  = dt_ap1 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
      dt -= dt_ap1;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
      
      //Phase 2: in the constant accel phase of velocity ramp-up.
      if(dt_vp1 > std::numeric_limits<double>::epsilon()) {
        double descended_accel = result_acc / max_acceleration;
        if( dt < dt_vp1 )
          dt_vp1 = dt;
        
        result_pos += ( result_vel + 0.5 * dt_vp1 * descended_accel ) * dt_vp1 / max_velocity;
        result_vel += dt_vp1 * descended_accel;
        result_desc_jerk = 0.0;
        
        dt -= dt_vp1;
        if(dt <= std::numeric_limits<double>::epsilon())
          return;
      };
      
      //Phase 3: in the jerk-down phase of velocity ramp-up.
      descended_jerk = -sgn_vp1;
      if( dt < dt_ap1 )
        dt_ap1 = dt;
      
      result_pos += ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration ) * dt_ap1 / max_velocity;
      result_vel += ( result_acc + 0.5 * dt_ap1 * descended_jerk ) * dt_ap1 / max_acceleration;
      result_acc += dt_ap1 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
      dt -= dt_ap1;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
    };
    
    //Phase 4: in the cruise phase.
    if( dt < dt_total )
      dt_total = dt;
    
    result_pos += dt_total * peak_velocity / max_velocity;
    result_vel = peak_velocity;
    result_acc = 0.0;
    result_desc_jerk = 0.0;
    
    dt -= dt_total;
    if(dt <= std::numeric_limits<double>::epsilon())
      return;
    
    if(dt_ap2 > std::numeric_limits<double>::epsilon()) {
      //Phase 5: in the jerk-up phase of velocity ramp-down.
      double descended_jerk = sgn_vp2;
      if( dt < dt_ap2 )
        dt_ap2 = dt;
      
      result_pos += ( result_vel + ( dt_ap2 * descended_jerk / 6.0 ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity;
      result_vel += ( 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration;
      result_acc  = dt_ap2 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
      dt -= dt_ap2;
      if(dt <= std::numeric_limits<double>::epsilon())
        return;
      
      //Phase 6: in the constant accel phase of velocity ramp-down.
      if(dt_vp2 > std::numeric_limits<double>::epsilon()) {
        double descended_accel = result_acc / max_acceleration;
        if( dt < dt_vp2 )
          dt_vp2 = dt;
        
        result_pos += ( result_vel + 0.5 * dt_vp2 * descended_accel ) * dt_vp2 / max_velocity;
        result_vel += dt_vp2 * descended_accel;
        result_desc_jerk = 0.0;
        
        dt -= dt_vp2;
        if(dt <= std::numeric_limits<double>::epsilon())
          return;
      };
      
      //Phase 7: in the jerk-down phase of velocity ramp-down.
      descended_jerk = -sgn_vp2;
      if( dt < dt_ap2 )
        dt_ap2 = dt;
      
      result_pos += ( result_vel + ( 0.5 * result_acc + (1.0 / 6.0) * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration ) * dt_ap2 / max_velocity;
      result_vel += ( result_acc + 0.5 * dt_ap2 * descended_jerk ) * dt_ap2 / max_acceleration;
      result_acc += dt_ap2 * descended_jerk;
      result_desc_jerk = descended_jerk;
      
    };
  };
  
  
  
  void sap_Ndof_compute_interpolated_values(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double peak_velocity, double max_velocity, double max_acceleration,
    double dt, double dt_total,
    double& result_pos, double& result_vel, 
    double& result_acc, double& result_desc_jerk) {
    
    sap_Ndof_compute_interpolated_values_closedform(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration, dt, dt_total,
      result_pos, result_vel, result_acc, result_desc_jerk);
    
  };
  
  
  
  
#endif
  
  
  
  
  inline void sap_Ndof_compute_ramp_dist_and_time(
    double v1, double v2, double vmax, double amax, 
    double& d_pos, double& dt) {
    
    using std::fabs;
    using std::sqrt;
    
    if( fabs(v2 - v1) >= amax ) {
      dt = fabs(v2 - v1) + amax;
      d_pos = 0.5 * dt * (v1 + v2) / vmax;
    } else {
      dt = 2.0 * sqrt(amax * fabs(v2 - v1));
      d_pos = 0.5 * dt * (v1 + v2) / vmax;
    };
    
  };
  
  struct sap_Ndof_no_cruise_calculator {
    double dp, v1, v2, vmax, amax;
    
    sap_Ndof_no_cruise_calculator(
      double a_dp, double a_v1, double a_v2, double a_vmax, double a_amax) : 
      dp(a_dp), v1(a_v1), v2(a_v2), vmax(a_vmax), amax(a_amax) { };
    
    double operator()(double vp) {
      double dp1, dt1, dp2, dt2;
      sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
      sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);
      
      if( dp < 0.0 )
        return dp1 + dp2 - dp;
      else 
        return dp - dp1 - dp2;
    };
    
  };
  
  
  static double sap_Ndof_compute_min_delta_time_bisection(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double& peak_velocity, double max_velocity, double max_acceleration) {
    using std::fabs;
    using std::sqrt;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      return 0.0;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      return std::numeric_limits<double>::infinity();
    };
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    
    double dt_ramp1 = 0.0, dp_ramp1 = 0.0, dt_ramp2 = 0.0, dp_ramp2 = 0.0;
    
    sap_Ndof_no_cruise_calculator nc_calc(end_position - start_position, start_velocity, end_velocity, max_velocity, max_acceleration);
    double peak_vel_low = -sign_p1_p0 * max_velocity;
    double peak_vel_hi  =  sign_p1_p0 * max_velocity;
    
    double delta_first_order = nc_calc(peak_vel_hi);
    if( delta_first_order > 0.0 ) {
      peak_velocity = peak_vel_hi;
      sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
      sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
      return delta_first_order + dt_ramp1 + dt_ramp2;
    };
    
    delta_first_order = nc_calc(peak_vel_low);
    if( delta_first_order < 0.0 ) {
      peak_velocity = peak_vel_low;
      sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
      sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
      return -delta_first_order + dt_ramp1 + dt_ramp2;
    };
    
    bisection_method(peak_vel_low, peak_vel_hi, nc_calc, 1e-6 * max_velocity);
    
    peak_velocity = peak_vel_hi;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
    return fabs(end_position - start_position - dp_ramp1 - dp_ramp2) + dt_ramp1 + dt_ramp2;
  };
  
  
  
  
  
#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION
  
  
  static double sap_Ndof_compute_min_delta_time_logged_version(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double& peak_velocity, 
    double max_velocity, double max_acceleration, std::ostream& log_output) {
    using std::fabs;
    using std::sqrt;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      log_output << " the two points are the same, within a tolerance margin.\n"
                 << " peak_velocity = " << peak_velocity << "\n"
                 << " min-dt = 0.0" << std::endl;
      return 0.0;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      log_output << " one of the end condition is infeasible!\n"
                 << " peak_velocity = 0.0\n"
                 << " min-dt = inf" << std::endl;
      return std::numeric_limits<double>::infinity();
    };
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    double sa_v0_v1_2 = sign_p1_p0 * max_acceleration * (start_velocity + end_velocity) * 0.5;
    double vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    double vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    
    double dt_ramp1 = 0.0, dp_ramp1 = 0.0, dt_ramp2 = 0.0, dp_ramp2 = 0.0;
    
    log_output << "-----------------------------------------------------------\n"
               << " sign_p1_p0 = sign(end_position - start_position) = " << sign_p1_p0 << "\n"
               << " sa_v0_v1_2 = sign_p1_p0 * max_acceleration * (start_velocity + end_velocity) * 0.5 = " << sa_v0_v1_2 << "\n"
               << " vm_p1_p0 = max_velocity * fabs(end_position - start_position) = " << vm_p1_p0 << "\n"
               << " vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity) = " << vsqr_avg << std::endl;
    
    // try to assume that a0 and a1 are a_max
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //           + vp/v_max (dt - fabs(vp - v0) - fabs(vp - v1) - 2*a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
    
    double delta_first_order = end_position - start_position - dp_ramp1 - dp_ramp2;
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - dp_ramp1 - dp_ramp2 = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + dt_ramp1 + dt_ramp2) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + dt_ramp1 + dt_ramp2;
    };
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = -sign_p1_p0 * max_velocity;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
    
    delta_first_order = end_position - start_position - dp_ramp1 - dp_ramp2;
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = -sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - dp_ramp1 - dp_ramp2 = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + dt_ramp1 + dt_ramp2) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + dt_ramp1 + dt_ramp2;
    };
    
    log_output << "-----------------------------------------------------------\n"
               << " looping on peak_velocity values ..." << std::endl;
    for(double cur_vp = sign_p1_p0 * max_velocity; cur_vp * sign_p1_p0 > -1.05 * max_velocity; cur_vp -= 0.1 * sign_p1_p0 * max_velocity) {
      double cur_dt_ramp1 = 0.0, cur_dp_ramp1 = 0.0, cur_dt_ramp2 = 0.0, cur_dp_ramp2 = 0.0;
      
      sap_Ndof_compute_ramp_dist_and_time(start_velocity, cur_vp, max_velocity, max_acceleration, cur_dp_ramp1, cur_dt_ramp1);
      sap_Ndof_compute_ramp_dist_and_time(cur_vp, end_velocity, max_velocity, max_acceleration, cur_dp_ramp2, cur_dt_ramp2);
      
      log_output << "  vp = " << cur_vp << " dt = " << (cur_dt_ramp1 + cur_dt_ramp2)
                 << " gives dp = " << sign_p1_p0 * (end_position - start_position - cur_dp_ramp1 - cur_dp_ramp2) << std::endl; 
      
    };
    
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //    p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //             + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max 
    //  or
    //    p1 - p0 == sqrt(fabs(vp - v0) * a_max) (vp + v0) / v_max 
    //             + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max 
    //  or
    //    p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //             + sqrt(fabs(vp - v1) * a_max) (vp + v1) / v_max 
    //  or
    //    p1 - p0 == sqrt(fabs(vp - v0) * a_max) (vp + v0) / v_max 
    //             + sqrt(fabs(vp - v1) * a_max) (vp + v1) / v_max 
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    double descrim = max_acceleration * max_acceleration - 4.0 * (sa_v0_v1_2 - vm_p1_p0 - vsqr_avg);
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am - 4.0 * (sa_v0_v1_2 - vm_p1_p0 - vsqr_avg) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (-sign_p1_p0 * am + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration) << std::endl;
        // this means that this root works for this case.
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (-sign_p1_p0 * am - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration) << std::endl;
        // this means that this root works for this case.
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration - 4.0 * (vm_p1_p0 - sa_v0_v1_2 - vsqr_avg);
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am - 4.0 * (vm_p1_p0 - sa_v0_v1_2 - vsqr_avg) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (sign_p1_p0 * am + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration) << std::endl;
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (sign_p1_p0 * am - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - start_velocity) + fabs(peak_velocity - end_velocity) + 2.0 * max_acceleration) << std::endl;
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
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - 0.5 * (fabs(vp - ve) + am) * (vp + ve) / vm - fabs(vp - vs) * (vp + vs) / vm = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * max_acceleration * end_velocity - 2.0 * vm_p1_p0 - 2.0 * start_velocity * start_velocity - end_velocity * end_velocity);
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * am * ve - 2.0 * vm_p1_p0 - 2.0 * vs * vs - ve * ve) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (-sign_p1_p0 * am / 3.0 + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means that this root works for this case.
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (-sign_p1_p0 * am / 3.0 - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means that this root works for this case.
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * ( 2.0 * vm_p1_p0 + 2.0 * start_velocity * start_velocity + end_velocity * end_velocity - sign_p1_p0 * max_acceleration * end_velocity );
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am / 9.0 - 4.0 / 3.0 * (2.0 * vm_p1_p0 + 2.0 * vs * vs + ve * ve - sign_p1_p0 * am * ve) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (sign_p1_p0 * am / 3.0 + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (sign_p1_p0 * am / 3.0 - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
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
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - fabs(vp - ve) * (vp + ve) / vm - 0.5 * (fabs(vp - vs) + am) * (vp + vs) / vm = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v1) (vp + v1) / v_max + (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * max_acceleration * start_velocity - 2.0 * vm_p1_p0 - start_velocity * start_velocity - 2.0 * end_velocity * end_velocity);
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am / 9.0 - 4.0 / 3.0 * (sign_p1_p0 * am * vs - 2.0 * vm_p1_p0 - vs * vs - 2.0 * ve * ve) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (-sign_p1_p0 * am / 3.0 + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means that this root works for this case.
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (-sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (-sign_p1_p0 * am / 3.0 - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp >= sign_p1_p0 * vs) && (sign_p1_p0 * vp >= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity >= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means that this root works for this case.
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = max_acceleration * max_acceleration / 9.0 - 4.0 / 3.0 * ( 2.0 * vm_p1_p0 + start_velocity * start_velocity + 2.0 * end_velocity * end_velocity - sign_p1_p0 * max_acceleration * start_velocity );
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = am * am / 9.0 - 4.0 / 3.0 * ( 2.0 * vm_p1_p0 + vs * vs + 2.0 * ve * ve - sign_p1_p0 * am * vs) = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a quadratic root for vp:
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 + sqrt(descrim)) * 0.5;
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = (sign_p1_p0 * am / 3.0 + sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        return 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      peak_velocity = (sign_p1_p0 * max_acceleration / 3.0 - sqrt(descrim)) * 0.5;
      
      log_output << "  peak_velocity = (sign_p1_p0 * am / 3.0 - sqrt(descrim)) * 0.5 = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (sign_p1_p0 * vp <= sign_p1_p0 * vs) && (sign_p1_p0 * vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * start_velocity) && (sign_p1_p0 * peak_velocity <= sign_p1_p0 * end_velocity)) {
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
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
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - fabs(vp - ve) * (vp + ve) / vm - fabs(vp - vs) * (vp + vs) / vm = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity)) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity);
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //  p1 - p0 == fabs(vp - v0) (vp + v0) / v_max + fabs(vp - v1) (vp + v1) / v_max
    
    //   first try if vp is more in the direction (p1-p0) than both v1 and v0:
    peak_velocity = sqrt(0.5 * vm_p1_p0 + vsqr_avg);
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sqrt(0.5 * vm_p1_p0 + vsqr_avg) = " << peak_velocity << "\n"
               << " 'if((fabs(vp) <= vm) && (vp >= sign_p1_p0 * vs) && (vp >= sign_p1_p0 * ve))'" << std::endl;
    
    if((fabs(peak_velocity) <= max_velocity) && (peak_velocity >= sign_p1_p0 * start_velocity) && (peak_velocity >= sign_p1_p0 * end_velocity)) {
      // this means that this root works for this case.
      peak_velocity *= sign_p1_p0;
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = " << peak_velocity << "\n"
                 << "  min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
      return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   else, try if vp is less in the direction (p1-p0) than both v0 and v1 (because in-between is impossible for min-dt):
    descrim = vsqr_avg - 0.5 * vm_p1_p0;
    
    log_output << "-----------------------------------------------------------\n"
               << " descrim = vsqr_avg - 0.5 * vm_p1_p0 = " << descrim << "\n"
               << " 'if(descrim > 0.0)'" << std::endl;
    
    if( descrim > 0.0 ) {
      // this means there exists a root for vp:
      peak_velocity = sqrt(descrim);
      
      log_output << "  condition passed!" << "\n"
                 << "  peak_velocity = sqrt(descrim) = " << peak_velocity << "\n"
                 << "  'if((fabs(vp) <= vm) && (vp <= sign_p1_p0 * vs) && (vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (peak_velocity <= sign_p1_p0 * start_velocity) && (peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        peak_velocity *= sign_p1_p0;
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
      
      log_output << "  'if((fabs(vp) <= vm) && (-vp <= sign_p1_p0 * vs) && (-vp <= sign_p1_p0 * ve))'" << std::endl;
      
      if((fabs(peak_velocity) <= max_velocity) && (-peak_velocity <= sign_p1_p0 * start_velocity) && (-peak_velocity <= sign_p1_p0 * end_velocity)) {
        // this means there is a valid solution for which vp is still in the direction of (p1-p0):
        peak_velocity *= -sign_p1_p0;
        log_output << "   condition passed!" << "\n"
                   << "   min-dt = " << (fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration) << std::endl;
        return fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
      };
    };
    
    log_output << "-----------------------------------------------------------\n"
               << "Reached the end of the calculation and did not find a suitable root! This should never happen!" << std::endl;
    
    // What the fuck!! This point should never be reached, unless the motion is completely impossible:
    peak_velocity = 0.0;
    return std::numeric_limits<double>::infinity();
  };
  
  
  
  double sap_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& peak_velocity, 
                                         double max_velocity, double max_acceleration) {
    
    static std::ofstream log_output("sap_mindt_log.txt");
    
    log_output << "\n\n0000000000000000000000000000000000000000000000000000000000000000000\n\n"
               << "Starting to compute the mininum delta-time for....\n"
               << "Positions: ( " << start_position << " , " << end_position << " )\n"
               << "Velocities: ( " << start_velocity << " , " << end_velocity << " )\n"
               << "Max Velocity: " << max_velocity << "\n"
               << "Max Acceleration: " << max_acceleration << std::endl;
    
    /*
    double min_delta_time = sap_Ndof_compute_min_delta_time_logged_version(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration, log_output);
    */
    
    double min_delta_time = sap_Ndof_compute_min_delta_time_bisection(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration);
    
    log_output << "The calculation of the minimum delta-time has finished and got the following:\n"
               << "Minimum Time: " << min_delta_time << "\n"
               << "Peak Velocity: " << peak_velocity << std::endl;
    
    log_output << "**********************************************\n"
               << "Checking the interpolation...." << std::endl;
    
    double result_pos, result_vel, result_acc, result_desc_jerk;
    
    sap_Ndof_compute_interpolated_values_incremental_logged_version(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration, min_delta_time, min_delta_time,
      result_pos, result_vel, result_acc, result_desc_jerk, log_output);
    
    log_output << "The interpolation resulted in the following state:\n"
               << "Position: " << result_pos << "\n"
               << "Velocity: " << result_vel << "\n"
               << "Acceleration: " << result_acc << "\n"
               << "Descended-Jerk: " << result_desc_jerk << std::endl;
    
    
    return min_delta_time;
  };
  
  
#else
  
  
  static double sap_Ndof_compute_min_delta_time_closedform(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double& peak_velocity, double max_velocity, double max_acceleration) {
    using std::fabs;
    using std::sqrt;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      return 0.0;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      return std::numeric_limits<double>::infinity();
    };
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    double sa_v0_v1_2 = sign_p1_p0 * max_acceleration * (start_velocity + end_velocity) * 0.5;
    double vm_p1_p0 = max_velocity * fabs(end_position - start_position);
    double vsqr_avg = 0.5 * (start_velocity * start_velocity + end_velocity * end_velocity);
    
    double dt_ramp1 = 0.0, dp_ramp1 = 0.0, dt_ramp2 = 0.0, dp_ramp2 = 0.0;
    
    // try to assume that a0 and a1 are a_max
    // giving:
    //  p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //           + vp/v_max (dt - fabs(vp - v0) - fabs(vp - v1) - 2*a_max)
    //           + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = sign_p1_p0 * max_velocity;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
    
    double delta_first_order = end_position - start_position - dp_ramp1 - dp_ramp2;
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + dt_ramp1 + dt_ramp2;
    };
    
    //   try to assume that vp = sign(p1-p0) * v_max
    peak_velocity = -sign_p1_p0 * max_velocity;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp_ramp2, dt_ramp2);
    
    delta_first_order = end_position - start_position - dp_ramp1 - dp_ramp2;
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + dt_ramp1 + dt_ramp2;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //    p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max 
    //             + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max 
    
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
    
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
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
    
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
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
    
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + 2.0 * fabs(peak_velocity - end_velocity) + 2.0 * fabs(peak_velocity - start_velocity);
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
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
  
  
  
  
  double sap_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& peak_velocity, 
                                         double max_velocity, double max_acceleration) {
    
    /*
    double min_delta_time = sap_Ndof_compute_min_delta_time_closedform(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration);
    */
    
    double min_delta_time = sap_Ndof_compute_min_delta_time_bisection(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration);
    
    return min_delta_time;
  };
  
  
  
#endif
  
  
  
  
  
  
  
  struct sap_Ndof_pos_diff_calculator {
    double dp, v1, v2, vmax, amax, dt;
    
    sap_Ndof_pos_diff_calculator(
      double a_dp, double a_v1, double a_v2, double a_vmax, double a_amax, double a_dt) : 
      dp(a_dp), v1(a_v1), v2(a_v2), vmax(a_vmax), amax(a_amax), dt(a_dt) { };
    
    double operator()(double vp) {
      using std::fabs;
      
      double dp1, dt1, dp2, dt2;
      sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
      sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);
//       if(dt + 1e-6 >= dt1 + dt2)
        return dp - dp1 - dp2 - vp / vmax * (dt - dt1 - dt2);
//       else {
//         double cruise_gap = dp - dp1 - dp2;
//         if( fabs(cruise_gap) < 1e-6 ) {
//           return dp - vp / vmax * dt + dt1 * (vp / vmax - dp1 / dt1) + dt2 * (vp / vmax - dp2 / dt2);
//         } else if( cruise_gap * vp > 0.0 )
//           return cruise_gap - vp / vmax * (dt - dt1 - dt2);
//         else
//           return cruise_gap + vp / vmax * (dt - dt1 - dt2);
//       };
    };
    
    double get_delta_time_diff(double vp) {
      double dp1, dt1, dp2, dt2;
      sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
      sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);
      return dt - dt1 - dt2;
    };
    
  };
  
  
#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION
  
  static void sap_Ndof_compute_peak_velocity_bisection(
    double start_position, double end_position,
    double start_velocity, double end_velocity, double& peak_velocity, 
    double max_velocity, double max_acceleration, double delta_time) {
    
    static std::ofstream log_output("sap_compute_vp_log.txt");
    
    log_output << "\n\n0000000000000000000000000000000000000000000000000000000000000000000\n\n"
               << "Starting to compute the peak-velocity for....\n"
               << "Positions: ( " << start_position << " , " << end_position << " )\n"
               << "Velocities: ( " << start_velocity << " , " << end_velocity << " )\n"
               << "Max Velocity: " << max_velocity << "\n"
               << "Max Acceleration: " << max_acceleration << "\n"
               << "Delta-Time: " << delta_time << std::endl;
    
    using std::fabs;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      return;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      RK_NOTICE(1," Warning: violation of the velocity bounds was detected on SAP interpolations!");
      return;
    };
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    
    sap_Ndof_pos_diff_calculator pd_calc(
      end_position - start_position, start_velocity, end_velocity, 
      max_velocity, max_acceleration, delta_time);
    
    if( fabs(pd_calc(sign_p1_p0 * max_velocity)) < 1e-6 * max_velocity ) {
      peak_velocity = sign_p1_p0 * max_velocity;
      log_output << "---------------------------------------------\n"
                 << " found max-velocity as peak = " << peak_velocity << std::endl;
      return;
    };
    
    log_output << "---------------------------------------------\n"
               << " looping on the peak velocities... " << std::endl;
    for(double cur_vp = sign_p1_p0 * max_velocity; cur_vp * sign_p1_p0 > -1.05 * max_velocity; cur_vp -= 0.1 * sign_p1_p0 * max_velocity) {
      log_output << " peak-velocity = " << cur_vp
                 << " residual dp = " << pd_calc(cur_vp) << std::endl;
    };
    
    double prev_vp = sign_p1_p0 * max_velocity;
    double prev_pd = pd_calc(prev_vp);
    for(double cur_vp = prev_vp - 0.1 * sign_p1_p0 * max_velocity; cur_vp * sign_p1_p0 > -1.05 * max_velocity; prev_vp = cur_vp, cur_vp -= 0.1 * sign_p1_p0 * max_velocity) {
      double cur_pd = pd_calc(cur_vp);
      if(cur_pd * prev_pd < 0.0) {
        bisection_method(prev_vp, cur_vp, pd_calc, 1e-6 * max_velocity);
        peak_velocity = cur_vp;
        return;
      };
    };
    
    peak_velocity = sign_p1_p0 * max_velocity;
    return;
  };
  
  
#else
  
  
  static void sap_Ndof_compute_peak_velocity_bisection(
    double start_position, double end_position,
    double start_velocity, double end_velocity, double& peak_velocity, 
    double max_velocity, double max_acceleration, double delta_time) {
    
    using std::fabs;
    
//     if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
//         ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
//       peak_velocity = start_velocity;
//       return;
//     };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      RK_NOTICE(1," Warning: violation of the velocity bounds was detected on SAP interpolations!");
      return;
    };
    
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    
    sap_Ndof_pos_diff_calculator pd_calc(
      end_position - start_position, start_velocity, end_velocity, 
      max_velocity, max_acceleration, delta_time);
    
//     if( pd_calc(1.02 * sign_p1_p0 * max_velocity) * pd_calc(0.98 * sign_p1_p0 * max_velocity) < 0.0 ) {
//       double peak_vel_low = 0.98 * sign_p1_p0 * max_velocity;
//       peak_velocity = 1.02 * sign_p1_p0 * max_velocity;
//       bisection_method(peak_vel_low, peak_velocity, pd_calc, 1e-6 * max_velocity);
//       return;
//     };
    
    double prev_vp = 1.03 * sign_p1_p0 * max_velocity;
    double prev_pd = pd_calc(prev_vp);
    for(double cur_vp = prev_vp - 0.02 * sign_p1_p0 * max_velocity; cur_vp * sign_p1_p0 > -1.04 * max_velocity; cur_vp -= 0.02 * sign_p1_p0 * max_velocity) {
      double cur_pd = pd_calc(cur_vp);
      if( cur_pd * prev_pd < 0.0 ) {
        double orig_cur_vp = cur_vp;
        ford3_method(prev_vp, cur_vp, pd_calc, 1e-7);
//         bisection_method(prev_vp, cur_vp, pd_calc, 1e-12 * max_velocity);
        cur_vp = (prev_vp + cur_vp) * 0.5;
        cur_pd = pd_calc(cur_vp);
        if( (pd_calc.get_delta_time_diff(cur_vp) >= -1e-3 * max_velocity) && ( fabs(cur_pd) < 1e-3 * max_velocity ) ) {
          peak_velocity = cur_vp;
          return;
        } else {
//           std::cout << " Could not really finish the root-finding for the peak-velocity:\n"
//                     << "  got a delta-time error of " << pd_calc.get_delta_time_diff(cur_vp) << "\n"
//                     << "  got vp-interval: " << prev_vp << " -- mid: " << cur_vp << " .. diff = " << (cur_vp - prev_vp) << "\n"
//                     << "  corresponding to solutions: " << pd_calc(prev_vp) << " -- mid: " << cur_pd << std::endl;
          cur_vp = orig_cur_vp;
          cur_pd = pd_calc(cur_vp);
        };
      };
      prev_vp = cur_vp;
      prev_pd = cur_pd;
    };
    
    RK_NOTICE(1," Warning: There was no solution to the peak-velocity for the given delta-time!");
    peak_velocity = -sign_p1_p0 * max_velocity;
    return;
  };
  
#endif
  
  
  
  static void sap_Ndof_compute_peak_velocity_closedform(
    double start_position, double end_position,
    double start_velocity, double end_velocity, double& peak_velocity, 
    double max_velocity, double max_acceleration, double delta_time) {
    // NOTE: Assume that delta-time is larger than minimum reachable delta-time 
    //       (avoid checking for that, even if it could mean that vp is higher than maximum)
    
    using std::fabs;
    using std::sqrt;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      return;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      RK_NOTICE(1," Warning: violation of the velocity bounds was detected on SAP interpolations!");
      return;
    };
    
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
  
  
  
  
  
  
  void sap_Ndof_compute_peak_velocity(double start_position, double end_position,
                                      double start_velocity, double end_velocity,
                                      double& peak_velocity, double max_velocity, 
                                      double max_acceleration, double delta_time) {
    
    
//     sap_Ndof_compute_min_delta_time(
//       start_position, end_position, 
//       start_velocity, end_velocity, 
//       peak_velocity, max_velocity, max_acceleration);
    
    sap_Ndof_compute_peak_velocity_bisection(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration, delta_time);
    
    
    
    double cf_pos, cf_vel, cf_acc, cf_desc_jerk;
    
    sap_Ndof_compute_interpolated_values_closedform(
      start_position, end_position, start_velocity, end_velocity,
      peak_velocity, max_velocity, max_acceleration,
      delta_time, delta_time,
      cf_pos, cf_vel, cf_acc, cf_desc_jerk);
    
    if( std::fabs(cf_pos - end_position) > 1e-3 ) {
      double dt1, dt2, dp1, dp2;
      
      sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity, max_velocity, max_acceleration, dp1, dt1);
      sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity, max_acceleration, dp2, dt2);
      
      sap_Ndof_pos_diff_calculator pd_calc(
        end_position - start_position, start_velocity, end_velocity, 
        max_velocity, max_acceleration, delta_time);
      
      std::cout << "The calculation of the peak velocity yielded a bad interpolated path!\n"
                << " Start position = " << start_position << "\n"
                << " End position   = " << end_position << "\n"
                << " Start velocity = " << start_velocity << "\n"
                << " End velocity   = " << end_velocity << "\n"
                << " Peak velocity  = " << peak_velocity << "\n"
                << " Delta-time     = " << delta_time << "\n"
                << " Delta-time-1   = " << dt1 << "\n"
                << " Delta-time-2   = " << dt2 << "\n"
                << " Delta-pos-1    = " << dp1 << "\n"
                << " Delta-pos-2    = " << dp2 << "\n"
                << " Calculated EDP = " << pd_calc(peak_velocity) << "\n"
                << " Actual EDP     = " << (cf_pos - end_position) << std::endl;
    };
    
    
  };
  
  
  
};


};

};






