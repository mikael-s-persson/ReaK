
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

#include <ReaK/topologies/interpolation/sustained_velocity_pulse_detail.hpp>

#include <ReaK/math/root_finders/bisection_method.hpp>

#include <ReaK/math/lin_alg/vect_alg.hpp>

#include <boost/bind.hpp>
#include <cmath>

namespace ReaK {

namespace pp {


namespace detail {

static double svp_compute_slack_time( double beta, double dt, double beta_0, double norm_delta,
                                      const double ( &coefs )[5] ) {
  using std::sqrt;
  using std::fabs;

  double c = ( norm_delta + coefs[4] * ( beta_0 * beta_0 - beta * beta ) );

  return dt - sqrt( coefs[0] * coefs[0] - beta * coefs[4] * ( 2.0 * coefs[1] - beta * coefs[4] ) )
         - sqrt( coefs[2] * coefs[2] - beta * coefs[4] * ( 2.0 * coefs[3] - beta * coefs[4] ) ) - fabs( c ) / beta;
};

static double svp_compute_derivative_slack_time( double beta, double dt, double beta_0, double norm_delta,
                                                 const double ( &coefs )[5] ) {
  using std::sqrt;
  using std::fabs;

  double c = ( norm_delta + coefs[4] * ( beta_0 * beta_0 - beta * beta ) );

  double term1 = sqrt( coefs[0] * coefs[0] - beta * coefs[4] * ( 2.0 * coefs[1] - beta * coefs[4] ) );
  double term2 = sqrt( coefs[2] * coefs[2] - beta * coefs[4] * ( 2.0 * coefs[3] - beta * coefs[4] ) );
  return fabs( c ) / ( beta * beta ) - ( c > 0.0 ? -2.0 : 2.0 ) * coefs[4]
         - coefs[4] * ( ( coefs[4] * beta - coefs[1] ) / term1 + ( coefs[4] * beta - coefs[3] ) / term2 );
};

static double svp_compute_travel_time( double beta, double beta_0, double norm_delta, const double ( &coefs )[5] ) {
  using std::sqrt;
  using std::fabs;

  double c = ( norm_delta + coefs[4] * ( beta_0 * beta_0 - beta * beta ) );

  return sqrt( coefs[0] * coefs[0] - beta * coefs[4] * ( 2.0 * coefs[1] - beta * coefs[4] ) )
         + sqrt( coefs[2] * coefs[2] - beta * coefs[4] * ( 2.0 * coefs[3] - beta * coefs[4] ) ) + fabs( c ) / beta;
};

static double svp_compute_derivative_travel_time( double beta, double beta_0, double norm_delta,
                                                  const double ( &coefs )[5] ) {
  using std::sqrt;
  using std::fabs;

  double c = ( norm_delta + coefs[4] * ( beta_0 * beta_0 - beta * beta ) );

  double term1 = sqrt( coefs[0] * coefs[0] - beta * coefs[4] * ( 2.0 * coefs[1] - beta * coefs[4] ) );
  double term2 = sqrt( coefs[2] * coefs[2] - beta * coefs[4] * ( 2.0 * coefs[3] - beta * coefs[4] ) );
  return coefs[4] * ( ( coefs[4] * beta - coefs[1] ) / term1 + ( coefs[4] * beta - coefs[3] ) / term2 )
         - fabs( c ) / ( beta * beta ) + ( c > 0.0 ? -2.0 : 2.0 ) * coefs[4];
};

double svp_solve_for_min_dt_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol ) {
  if( svp_compute_derivative_travel_time( 1.0, beta, norm_delta, coefs ) > 0.0 ) {
    double upper = 1.0;
    double lower = 0.5;
    while( ( lower < 0.99 ) && ( svp_compute_derivative_travel_time( lower, beta, norm_delta, coefs ) > 0.0 ) ) {
      lower += 0.5 * ( upper - lower );
    };
    if( lower < 0.99 ) {
      bisection_method( lower, upper, boost::bind( svp_compute_derivative_travel_time, _1, boost::cref( beta ),
                                                   boost::cref( norm_delta ), boost::cref( coefs ) ),
                        num_tol );
    } else {
      upper = 1.0;
      lower = 0.5;
      while( svp_compute_derivative_travel_time( lower, beta, norm_delta, coefs ) > 0.0 ) {
        upper = lower;
        lower *= 0.5;
      };
      bisection_method( lower, upper, boost::bind( svp_compute_derivative_travel_time, _1, boost::cref( beta ),
                                                   boost::cref( norm_delta ), boost::cref( coefs ) ),
                        num_tol );
    };
    // make sure that the second root does not cause a reversal of the travel direction:
    if( ( norm_delta > 0.0 ) && ( norm_delta + coefs[4] * ( beta * beta - upper * upper ) < 0.0 ) ) {
      upper = std::sqrt( beta * beta + norm_delta / coefs[4] );
    };
    beta = upper;
  } else {
    beta = 1.0;
  };
  return beta;
};

bool svp_min_dt_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                           double& result ) {
  result = svp_compute_travel_time( beta, beta, norm_delta, coefs );
  return true;
};

double svp_solve_for_no_slack_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                                    double delta_time ) {
  double beta_peak1 = 1.0;
  double beta_peak2 = 5.0;

  if( svp_compute_slack_time( 1.0, delta_time, beta, norm_delta, coefs ) > 0.0 ) {
    // means I have a single root in the interval, so I can solve for it:
    double beta_low = 0.5;
    while( svp_compute_slack_time( beta_low, delta_time, beta, norm_delta, coefs ) > 0.0 ) {
      beta_peak1 = beta_low;
      beta_low *= 0.5;
    };
    bisection_method( beta_low, beta_peak1,
                      boost::bind( svp_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( beta ),
                                   boost::cref( norm_delta ), boost::cref( coefs ) ),
                      num_tol );
  } else {
    // This means that I must have either a parabola-looking curve, or it never goes positive.
    // so, find the maximum in the interval, by finding the zero of the derivative:
    double beta_low = 0.5;
    while( svp_compute_derivative_slack_time( beta_low, delta_time, beta, norm_delta, coefs ) < 0.0 ) {
      beta_peak1 = beta_low;
      beta_low *= 0.5;
    };
    bisection_method( beta_low, beta_peak1,
                      boost::bind( svp_compute_derivative_slack_time, _1, boost::cref( delta_time ),
                                   boost::cref( beta ), boost::cref( norm_delta ), boost::cref( coefs ) ),
                      num_tol );
    if( svp_compute_slack_time( beta_peak1, delta_time, beta, norm_delta, coefs ) > 0.0 ) {
      // this means the maximum slack-time is actually positive, meaning there must be a root on either side.
      beta_peak2 = beta_peak1;
      beta_low = 0.5 * beta_peak1;
      while( svp_compute_slack_time( beta_low, delta_time, beta, delta_time, coefs ) > 0.0 ) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method( beta_low, beta_peak1,
                        boost::bind( svp_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( beta ),
                                     boost::cref( norm_delta ), boost::cref( coefs ) ),
                        num_tol );
      beta_low = beta_peak2;
      beta_peak2 = 1.0;
      bisection_method( beta_low, beta_peak2,
                        boost::bind( svp_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( beta ),
                                     boost::cref( norm_delta ), boost::cref( coefs ) ),
                        num_tol );

      // make sure that the second root does not cause a reversal of the travel direction:
      if( norm_delta + coefs[4] * ( beta * beta - beta_peak2 * beta_peak2 ) < 0.0 ) {
        beta_peak2 = 5.0;
      };
    };
  };

  if( std::fabs( beta - beta_peak1 ) < std::fabs( beta - beta_peak2 ) ) {
    beta = beta_peak1;
  } else {
    beta = beta_peak2;
  };
  if( beta <= num_tol )
    beta = num_tol;
  return beta;
};

bool svp_no_slack_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol, double& slack,
                             double delta_time ) {
  slack = svp_compute_slack_time( beta, delta_time, beta, norm_delta, coefs );
  return ( std::fabs( slack ) < 100.0 * num_tol * delta_time );
};
};
};
};
