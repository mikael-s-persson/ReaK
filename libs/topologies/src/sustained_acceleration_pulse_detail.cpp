
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


#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_detail.hpp>

#include <ReaK/math/lin_alg/vect_alg.hpp>

#include <ReaK/math/root_finders/bisection_method.hpp>
#include <ReaK/math/root_finders/secant_method.hpp>

#include <boost/bind.hpp>
#include <cmath>

namespace ReaK {

namespace pp {

namespace detail {

static vect< double, 2 > sap_compute_projected_deltas( double beta, const double ( &coefs )[5], double dt_amax,
                                                       double& A_0, double& A_1 ) {
  using std::sqrt;
  vect< double, 2 > result;
  A_0 = sqrt( coefs[0] * coefs[0] - beta * coefs[4] * ( 2.0 * coefs[1] - beta * coefs[4] ) );
  double t_0 = coefs[1] / coefs[4];
  if( A_0 < dt_amax )
    result[0] = 0.5 * sqrt( dt_amax * A_0 ) * ( t_0 + 3.0 * beta );
  else
    result[0] = 0.5 * ( ( dt_amax + A_0 ) * ( beta + t_0 ) + dt_amax * dt_amax / A_0 * ( beta - t_0 ) );
  A_1 = sqrt( coefs[2] * coefs[2] - beta * coefs[4] * ( 2.0 * coefs[3] - beta * coefs[4] ) );
  double t_1 = coefs[3] / coefs[4];
  if( A_1 < dt_amax )
    result[1] = 0.5 * sqrt( dt_amax * A_1 ) * ( t_1 + 3.0 * beta );
  else
    result[1] = 0.5 * ( ( dt_amax + A_1 ) * ( beta + t_1 ) + dt_amax * dt_amax / A_1 * ( beta - t_1 ) );
  return result;
};

static vect< double, 2 > sap_compute_projected_deltas( double beta, const double ( &coefs )[5], double dt_amax ) {
  double A_0, A_1;
  return sap_compute_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );
};

static vect< double, 2 > sap_compute_derivative_projected_deltas( double beta, const double ( &coefs )[5],
                                                                  double dt_amax, double A_0, double A_1 ) {
  using std::sqrt;
  using std::fabs;
  vect< double, 2 > result;
  if( A_0 < dt_amax ) {
    result[0] = beta * coefs[4] * ( beta * coefs[4] - coefs[1] ) / dt_amax
                + ( coefs[0] * coefs[0] - coefs[1] * coefs[1] ) / dt_amax + 0.5 * dt_amax;
  } else {
    result[0] = beta * coefs[4] * ( beta * coefs[4] - coefs[1] ) / A_0
                + 0.5 * ( ( 1.0 + dt_amax * dt_amax / ( A_0 * A_0 ) ) * ( coefs[0] * coefs[0] - coefs[1] * coefs[1] )
                          / A_0 + dt_amax );
  };
  if( A_1 < dt_amax ) {
    result[1] = beta * coefs[4] * ( beta * coefs[4] - coefs[3] ) / dt_amax
                + ( coefs[2] * coefs[2] - coefs[3] * coefs[3] ) / dt_amax + 0.5 * dt_amax;
  } else {
    result[1] = beta * coefs[4] * ( beta * coefs[4] - coefs[3] ) / A_1
                + 0.5 * ( ( 1.0 + dt_amax * dt_amax / ( A_1 * A_1 ) ) * ( coefs[2] * coefs[2] - coefs[3] * coefs[3] )
                          / A_1 + dt_amax );
  };
  return result;
};

static double sap_compute_slack_time( double beta, double dt, const vect< double, 2 >& deltas_0, double norm_delta,
                                      const double ( &coefs )[5], double dt_amax ) {
  using std::sqrt;
  using std::fabs;

  double A_0, A_1;
  vect< double, 2 > deltas_1 = sap_compute_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  if( A_0 < dt_amax )
    A_0 = sqrt( 4.0 * A_0 * dt_amax );
  else
    A_0 += dt_amax;
  if( A_1 < dt_amax )
    A_1 = sqrt( 4.0 * A_1 * dt_amax );
  else
    A_1 += dt_amax;

  double c = ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] );

  return dt - A_0 - A_1 - fabs( c ) / beta;
};

static double sap_compute_derivative_slack_time( double beta, double dt, const vect< double, 2 >& deltas_0,
                                                 double norm_delta, const double ( &coefs )[5], double dt_amax ) {
  using std::sqrt;
  using std::fabs;

  double A_0, A_1;
  vect< double, 2 > deltas_1 = sap_compute_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  double c = ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] );

  vect< double, 2 > deltas_dot_1 = sap_compute_derivative_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  double c_dot = -deltas_dot_1[0] - deltas_dot_1[1];

  double dt0, dt1;
  if( A_0 < dt_amax ) {
    dt0 = coefs[4] * ( coefs[4] * beta - coefs[1] ) / dt_amax;
  } else {
    dt0 = coefs[4] * ( coefs[4] * beta - coefs[1] ) / A_0;
  };
  if( A_1 < dt_amax ) {
    dt1 = coefs[4] * ( coefs[4] * beta - coefs[3] ) / dt_amax;
  } else {
    dt1 = coefs[4] * ( coefs[4] * beta - coefs[3] ) / A_1;
  };

  return fabs( c ) / ( beta * beta ) - ( c > 0.0 ? c_dot : -c_dot ) / beta - dt0 - dt1;
};


static double sap_compute_travel_time( double beta, const vect< double, 2 >& deltas_0, double norm_delta,
                                       const double ( &coefs )[5], double dt_amax ) {
  using std::sqrt;
  using std::fabs;

  double A_0, A_1;
  vect< double, 2 > deltas_1 = sap_compute_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  if( A_0 < dt_amax )
    A_0 = sqrt( 4.0 * A_0 * dt_amax );
  else
    A_0 += dt_amax;
  if( A_1 < dt_amax )
    A_1 = sqrt( 4.0 * A_1 * dt_amax );
  else
    A_1 += dt_amax;

  double c = ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] );

  return A_0 + A_1 + fabs( c ) / beta;
};

static double sap_compute_derivative_travel_time( double beta, const vect< double, 2 >& deltas_0, double norm_delta,
                                                  const double ( &coefs )[5], double dt_amax ) {
  using std::sqrt;
  using std::fabs;

  double A_0, A_1;
  vect< double, 2 > deltas_1 = sap_compute_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  double c = ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] );

  vect< double, 2 > deltas_dot_1 = sap_compute_derivative_projected_deltas( beta, coefs, dt_amax, A_0, A_1 );

  double c_dot = -deltas_dot_1[0] - deltas_dot_1[1];

  double dt0, dt1;
  if( A_0 < dt_amax ) {
    dt0 = coefs[4] * ( coefs[4] * beta - coefs[1] ) / dt_amax;
  } else {
    dt0 = coefs[4] * ( coefs[4] * beta - coefs[1] ) / A_0;
  };
  if( A_1 < dt_amax ) {
    dt1 = coefs[4] * ( coefs[4] * beta - coefs[3] ) / dt_amax;
  } else {
    dt1 = coefs[4] * ( coefs[4] * beta - coefs[3] ) / A_1;
  };

  return dt0 + dt1 - fabs( c ) / ( beta * beta ) + ( c > 0.0 ? c_dot : -c_dot ) / beta;
};

double sap_solve_for_min_dt_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                                  double dt_amax ) {
  vect< double, 2 > deltas_0 = sap_compute_projected_deltas( beta, coefs, dt_amax );
  if( sap_compute_derivative_travel_time( 1.0, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) {
    double upper = 1.0;
    double lower = 0.1;
    while( ( lower < 0.99 )
           && ( sap_compute_derivative_travel_time( lower, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) ) {
      lower += 0.5 * ( upper - lower );
    };
    if( lower < 0.99 ) {
      brent_method( lower, upper,
                    boost::bind( sap_compute_derivative_travel_time, _1, boost::cref( deltas_0 ),
                                 boost::cref( norm_delta ), boost::cref( coefs ), boost::cref( dt_amax ) ),
                    num_tol );
    } else {
      upper = 1.0;
      lower = 0.9;
      while( sap_compute_derivative_travel_time( lower, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) {
        upper = lower;
        lower *= 0.5;
      };
      brent_method( lower, upper,
                    boost::bind( sap_compute_derivative_travel_time, _1, boost::cref( deltas_0 ),
                                 boost::cref( norm_delta ), boost::cref( coefs ), boost::cref( dt_amax ) ),
                    num_tol );
    };
    // make sure that the second root does not cause a reversal of the travel direction:
    vect< double, 2 > deltas_1 = sap_compute_projected_deltas( upper, coefs, dt_amax );
    if( ( norm_delta > 0.0 ) && ( ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] ) < 0.0 ) ) {
      upper = std::sqrt( beta * beta + norm_delta / coefs[4] );
    };
    beta = upper;
  } else {
    beta = 1.0;
  };
  return beta;
};

bool sap_min_dt_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol, double& result,
                           double dt_amax ) {
  result
    = sap_compute_travel_time( beta, sap_compute_projected_deltas( beta, coefs, dt_amax ), norm_delta, coefs, dt_amax );
  return true;
};

double sap_solve_for_no_slack_beta( double beta, double norm_delta, const double ( &coefs )[5], double num_tol,
                                    double delta_time, double dt_amax ) {
  double beta_peak1 = 1.0;
  double beta_peak2 = 5.0;

  vect< double, 2 > deltas_0 = sap_compute_projected_deltas( beta, coefs, dt_amax );

  if( sap_compute_slack_time( 1.0, delta_time, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) {
    // means I have a single root in the interval, so I can solve for it:
    double beta_low = 0.5;
    while( sap_compute_slack_time( beta_low, delta_time, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) {
      beta_peak1 = beta_low;
      beta_low *= 0.5;
    };
    bisection_method( beta_low, beta_peak1,
                      boost::bind( sap_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( deltas_0 ),
                                   boost::cref( norm_delta ), boost::cref( coefs ), boost::cref( dt_amax ) ),
                      num_tol );
  } else {
    // This means that I must have either a parabola-looking curve, or it never goes positive.
    // so, find the maximum in the interval, by finding the zero of the derivative:
    double beta_low = 0.5;
    while( sap_compute_derivative_slack_time( beta_low, delta_time, deltas_0, norm_delta, coefs, dt_amax ) < 0.0 ) {
      beta_peak1 = beta_low;
      beta_low *= 0.5;
    };
    bisection_method( beta_low, beta_peak1,
                      boost::bind( sap_compute_derivative_slack_time, _1, boost::cref( delta_time ),
                                   boost::cref( deltas_0 ), boost::cref( norm_delta ), boost::cref( coefs ),
                                   boost::cref( dt_amax ) ),
                      num_tol );
    if( sap_compute_slack_time( beta_peak1, delta_time, deltas_0, norm_delta, coefs, dt_amax ) > 0.0 ) {
      // this means the maximum slack-time is actually positive, meaning there must be a root on either side.
      beta_peak2 = beta_peak1;
      beta_low = 0.5 * beta_peak1;
      while( sap_compute_slack_time( beta_low, delta_time, deltas_0, delta_time, coefs, dt_amax ) > 0.0 ) {
        beta_peak1 = beta_low;
        beta_low *= 0.5;
      };
      bisection_method( beta_low, beta_peak1,
                        boost::bind( sap_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( deltas_0 ),
                                     boost::cref( norm_delta ), boost::cref( coefs ), boost::cref( dt_amax ) ),
                        num_tol );
      beta_low = beta_peak2;
      beta_peak2 = 1.0;
      bisection_method( beta_low, beta_peak2,
                        boost::bind( sap_compute_slack_time, _1, boost::cref( delta_time ), boost::cref( deltas_0 ),
                                     boost::cref( norm_delta ), boost::cref( coefs ), boost::cref( dt_amax ) ),
                        num_tol );

      // make sure that the second root does not cause a reversal of the travel direction:
      vect< double, 2 > deltas_1 = sap_compute_projected_deltas( beta_peak2, coefs, dt_amax );
      if( ( norm_delta > 0.0 ) && ( ( norm_delta - deltas_1[0] - deltas_1[1] + deltas_0[0] + deltas_0[1] ) < 0.0 ) ) {
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

bool sap_no_slack_predicate( double beta, double norm_delta, const double ( &coefs )[5], double num_tol, double& slack,
                             double delta_time, double dt_amax ) {
  slack = sap_compute_slack_time( beta, delta_time, sap_compute_projected_deltas( beta, coefs, dt_amax ), norm_delta,
                                  coefs, dt_amax );
  return ( std::fabs( slack ) < 100.0 * num_tol * delta_time );
};
};
};
};
