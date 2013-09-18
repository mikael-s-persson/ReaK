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
    
    if( ( fabs(end_position - start_position) < 1e-6 * max_velocity ) &&
        ( fabs(end_velocity - start_velocity) < 1e-6 * max_velocity ) ) {
      peak_velocity = start_velocity;
      return 0.0;
    };
    
    if( ( fabs(start_velocity) > max_velocity ) || ( fabs(end_velocity) > max_velocity ) ) {
      peak_velocity = 0.0;
      RK_NOTICE(1," Warning: violation of the velocity bounds was detected on SVP interpolations!");
      return std::numeric_limits<double>::infinity();
    };
    
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
  
  
};


};

};









