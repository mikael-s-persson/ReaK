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

#include "sustained_velocity_pulse_Ndof_detail.hpp"

#include <cmath>

namespace ReaK {

namespace pp {
  
namespace detail {
  
  
  double svp_Ndof_compute_min_delta_time(double start_position, double end_position,
                                         double start_velocity, double end_velocity,
                                         double& peak_velocity, double max_velocity) {
    using std::fabs;
    using std::sqrt;
    
    // try to assume that peak_velocity = sign(p1 - p0) * max_velocity
    double sign_p1_p0 = 1.0;
    if(start_position > end_position)
      sign_p1_p0 = -1.0;
    peak_velocity = sign_p1_p0 * max_velocity;
    double descended_peak_velocity = peak_velocity / max_velocity;
    double delta_first_order = end_position - start_position
      - (0.5 * fabs(peak_velocity -   end_velocity)) * (descended_peak_velocity +   end_velocity / max_velocity)
      - (0.5 * fabs(peak_velocity - start_velocity)) * (descended_peak_velocity + start_velocity / max_velocity);
    if(delta_first_order * peak_velocity > 0.0) {
      // this means that we guessed correctly (we can reach max cruise speed in the direction of the end-position):
      return fabs(delta_first_order) + fabs(peak_velocity - end_velocity) + fabs(peak_velocity - start_velocity);
    };
    // if not, then can try to see if we simply can quite reach max velocity before having to ramp-down:
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
  
};


};

};









