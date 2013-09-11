
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


#define RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION


#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION

#include <fstream>

namespace ReaK {

namespace pp {
  
  
namespace detail {

  static double sap_Ndof_compute_min_delta_time_logged_version(
    double start_position, double end_position,
    double start_velocity, double end_velocity,
    double& peak_velocity, 
    double max_velocity, double max_acceleration) {
    using std::fabs;
    using std::sqrt;
    
    static std::ofstream log_output("sap_mindt_log.txt");
    
    log_output << "\n\n0000000000000000000000000000000000000000000000000000000000000000000\n\n"
               << "Starting to compute the mininum delta-time for....\n"
               << "Positions: ( " << start_position << " , " << end_position << " )\n"
               << "Velocities: ( " << start_velocity << " , " << end_velocity << " )\n"
               << "Max Velocity: " << max_velocity << "\n"
               << "Max Acceleration: " << max_acceleration << std::endl;
    
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
    double delta_first_order = end_position - start_position
      - 0.5 * (fabs(peak_velocity -   end_velocity) + max_acceleration) * (peak_velocity +   end_velocity) / max_velocity
      - 0.5 * (fabs(peak_velocity - start_velocity) + max_acceleration) * (peak_velocity + start_velocity) / max_velocity;
    
    log_output << "-----------------------------------------------------------\n"
               << " peak_velocity = sign_p1_p0 * max_velocity = " << peak_velocity << "\n"
               << " delta_first_order = pe - ps - 0.5 * (fabs(vp - ve) + am) * (vp + ve) / vm - 0.5 * (fabs(vp - vs) + am) * (vp + vs) / vm = " << delta_first_order << "\n"
               << " 'if(delta_first_order * peak_velocity > 0.0)' with delta_first_order * peak_velocity = " << (delta_first_order * peak_velocity) << std::endl;
    
    if(delta_first_order * peak_velocity > 0.0) {
      log_output << "  condition passed!" << "\n"
                 << "  min-dt = " << (fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + 2.0 * max_acceleration) << std::endl;
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + 2.0 * max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
    //   this assumes that we have 
    //    p1 - p0 == (fabs(vp - v0) + a_max) (vp + v0) / 2*v_max + (fabs(vp - v1) + a_max) (vp + v1) / 2*v_max
    
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

};

};

};

#endif




namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  
  double sap_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& peak_velocity, 
                                         double max_velocity, double max_acceleration) {
    using std::fabs;
    using std::sqrt;
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_acceleration ) ) {
      peak_velocity = start_velocity;
      return 0.0;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      RK_NOTICE(1," Warning: violation of the velocity bounds was detected on SAP interpolations!");
      return std::numeric_limits<double>::infinity();
    };
    
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
    double delta_first_order = end_position - start_position
      - 0.5 * (fabs(peak_velocity -   end_velocity) + max_acceleration) * (peak_velocity +   end_velocity) / max_velocity
      - 0.5 * (fabs(peak_velocity - start_velocity) + max_acceleration) * (peak_velocity + start_velocity) / max_velocity;
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-point):
      return fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity) + 2.0 * max_acceleration;
    };
    
    //   if not, then can try to see if we simply can't quite reach max velocity before having to ramp-down:
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
    
#ifdef RK_SAP_DETAIL_IMPLEMENTATION_USE_LOGGED_VERSION
    sap_Ndof_compute_min_delta_time_logged_version(start_position, end_position, start_velocity, end_velocity,
                                                   peak_velocity, max_velocity, max_acceleration);
#endif
    
    return std::numeric_limits<double>::infinity();
  };
  
  void sap_Ndof_compute_peak_velocity(double start_position, double end_position,
                                      double start_velocity, double end_velocity,
                                      double& peak_velocity, double max_velocity, 
                                      double max_acceleration, double delta_time) {
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
  
  
  
};


};

};






