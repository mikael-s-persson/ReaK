/**
 * \file nl_interior_points_methods.hpp
 *
 * The following library provides methods to perform constrained non-linear optimization using
 * interior-points methods (or penalty-barrier methods). One option is a trust-region approach
 * that first solves a quadratic program within a trust-region and then solves an
 * equality-constrained quadratic program by the projected conjugate gradient method (this
 * method is essentially an extension of the SQP method to also include inequality constraints).
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_NL_INTERIOR_POINTS_METHODS_HPP
#define REAK_NL_INTERIOR_POINTS_METHODS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/math/lin_alg/mat_svd_method.hpp>

#include "newton_search_directions.hpp"
#include "limit_functions.hpp"
#include "line_search.hpp"
#include "trust_region_search.hpp"
#include "hessian_approx_update.hpp"

#include "newton_methods.hpp"
#include "quadratic_programs.hpp"
#include "augmented_lagrangian_methods.hpp"

namespace ReaK {


namespace optim {


namespace detail {


template < typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename EqFunction,
           typename EqJacFunction, typename IneqFunction, typename IneqJacFunction, typename TrustRegionSolver,
           typename LimitFunction >
void nl_intpoint_method_tr_impl(
  Function f, GradFunction df, HessianFunction fill_hessian, EqFunction g, EqJacFunction fill_g_jac, IneqFunction h,
  IneqJacFunction fill_h_jac, Vector& x, typename vect_traits< Vector >::value_type max_radius,
  typename vect_traits< Vector >::value_type mu, unsigned int max_iter, TrustRegionSolver solve_step,
  LimitFunction impose_limits,
  typename vect_traits< Vector >::value_type abs_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type kappa = typename vect_traits< Vector >::value_type( 1e-4 ),
  typename vect_traits< Vector >::value_type tau = typename vect_traits< Vector >::value_type( 0.995 ) ) {
  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;
  using std::sqrt;
  using std::fabs;
  using std::log;

  SizeType N = x.size();
  Vector ht_value = h( x );
  SizeType K = ht_value.size();
  mat< ValueType, mat_structure::rectangular > Jac_h( K, N );
  fill_h_jac( Jac_h, x, ht_value );

  // compute initial slack vector and roughly adjust x if needed.
  Vector s = ht_value;
  for( SizeType i = 0; i < K; ++i )
    s[i] = ValueType( 1.0 );
  ValueType min_s( 0.0 );
  for( SizeType i = 0; i < K; ++i )
    if( s[i] < min_s )
      min_s = s[i];
  min_s *= ValueType( -0.1 );

  while( min_s > abs_tol ) {
    // RK_NOTICE(1," s = " << s << " x = " << x << " h = " << ht_value << " Jac = " << Jac_h);
    Vector e_s( s );
    for( SizeType i = 0; i < K; ++i )
      e_s[i] = min_s;
    s += e_s;
    Vector dx = x;
    mat_vect_adaptor< Vector > dx_mat( dx );
    minnorm_QR( Jac_h, dx_mat, mat_vect_adaptor< Vector >( e_s ), abs_tol );
    // RK_NOTICE(1," e_s = " << e_s << " dx = " << dx);
    x += dx;
    ht_value = h( x );
    fill_h_jac( Jac_h, x, ht_value );
    s = ht_value;
    min_s = ValueType( 0.0 );
    for( SizeType i = 0; i < K; ++i )
      if( s[i] < min_s )
        min_s = s[i];
    min_s *= ValueType( -0.1 );
  };
  if( min_s > ValueType( 0.0 ) )
    s += vect_scalar< double >( K, ValueType( 100.0 ) * abs_tol );

  Vector gt_value = g( x );
  SizeType M = gt_value.size();
  Vector c;
  c.resize( M + K );
  c[range( 0, M )] = gt_value;

  mat< ValueType, mat_structure::rectangular > Jac_g( M, N );
  fill_g_jac( Jac_g, x, gt_value );

  if( ( M == 0 ) && ( K == 0 ) ) { // this means it is an unconstrained problem. TODO change this to dispatch on the
                                   // type of the fill-hessian functor.
    newton_method_tr_impl( f, df, fill_hessian, x, max_radius, max_iter, solve_step, impose_limits, abs_tol, abs_tol,
                           kappa );
    return;
  };

  mat< ValueType, mat_structure::diagonal > mS( K );
  mat< ValueType, mat_structure::nil > Zero_mk( M, K );

  mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >, mat< ValueType, mat_structure::nil > >
    Jac_upper( Jac_g, Zero_mk );
  mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >, mat< ValueType, mat_structure::diagonal > >
    Jac_lower( Jac_h, mS );

  mat_const_ref_vert_cat< mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >,
                                                   mat< ValueType, mat_structure::nil > >,
                          mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >,
                                                   mat< ValueType, mat_structure::diagonal > > > Jac_aug( Jac_upper,
                                                                                                          Jac_lower );


  mS = mat< ValueType, mat_structure::diagonal >( -s );
  c[range( M, M + K )] = ht_value - s;

  ValueType gt_norm = norm_2( c[range( 0, M )] ) / sqrt( ValueType( M ) );
  ValueType ht_norm = norm_2( c[range( M, M + K )] ) / sqrt( ValueType( K ) );

  ValueType radius = ValueType( 0.5 ) * max_radius;
  ValueType x_value = f( x );
  Vector x_grad = df( x );
  mat< ValueType, mat_structure::symmetric > H( mat< ValueType, mat_structure::identity >( x.size() ) );
  fill_hessian( H, x, x_value, x_grad );


  Vector xt = x;
  Vector p_x = x;
  Vector st = s;
  Vector p_s = s;
  Vector v;
  v.resize( N + K );
  Vector p = v;
  Vector p_grad;
  p_grad.resize( N + K );
  p_grad[range( 0, N )] = x_grad;
  p_grad[range( N, N + K )] = vect_scalar< ValueType >( K, -mu );
  ValueType norm_p = std::numeric_limits< ValueType >::max();
  ValueType norm_v = std::numeric_limits< ValueType >::max();
  Vector r = c;
  // ValueType rho(0.1);
  // ValueType tol_mu = 1e-2 + tol;

  ValueType log_s( 0.0 );
  for( SizeType i = 0; i < K; ++i )
    log_s += log( s[i] );


  ValueType c_norm_star = norm_2( c );
  ValueType nu( std::numeric_limits< ValueType >::min() );

  Vector yz;
  yz.resize( M + K );
  mat_vect_adaptor< Vector > yz_mat( yz );
  vect_ref_view< Vector > y( yz[range( 0, M )] );
  vect_ref_view< Vector > z( yz[range( M, M + K )] );
  try {
    linlsq_QR( transpose_view( Jac_aug ), yz_mat, mat_vect_adaptor< Vector >( p_grad ), abs_tol );
  } catch( singularity_error& ) {
    SVD_linlsqsolver()( transpose_view( Jac_aug ), yz_mat, mat_vect_adaptor< Vector >( p_grad ), abs_tol );
  };
  //    for(SizeType i = 0; i < K; ++i)
  //      if(z[i] < ValueType(0.0))
  //        z[i] = mu;
  for( SizeType i = 0; i < K; ++i )
    if( z[i] < ValueType( 0.0 ) )
      z[i] = mu / s[i];

  ValueType norm_star = norm_2( x_grad ) + norm_2( y * Jac_g ) + norm_2( z * Jac_h );
  Vector l( x_grad - y * Jac_g - z * Jac_h );
  Vector lt = l;
  norm_star = norm_2( l ) / sqrt( ValueType( N ) );

  mat< ValueType, mat_structure::diagonal > SES( K );
  for( SizeType i = 0; i < K; ++i )
    SES( i, i ) = z[i] * s[i];

  mat< ValueType, mat_structure::nil > Zero_nk( N, K );
  mat< ValueType, mat_structure::nil > Zero_kn( K, N );

  mat_const_ref_horiz_cat< mat< ValueType, mat_structure::symmetric >, mat< ValueType, mat_structure::nil > > H_upper(
    H, Zero_nk );
  mat_const_ref_horiz_cat< mat< ValueType, mat_structure::nil >, mat< ValueType, mat_structure::diagonal > > H_lower(
    Zero_kn, SES );

  mat_const_ref_vert_cat< mat_const_ref_horiz_cat< mat< ValueType, mat_structure::symmetric >,
                                                   mat< ValueType, mat_structure::nil > >,
                          mat_const_ref_horiz_cat< mat< ValueType, mat_structure::nil >,
                                                   mat< ValueType, mat_structure::diagonal > > > H_aug( H_upper,
                                                                                                        H_lower );

  ValueType Err_value = norm_2( SES * vect_scalar< ValueType >( K, 1.0 ) ) / sqrt( ValueType( K ) );
  if( Err_value < norm_star )
    Err_value = norm_star;
  if( Err_value < gt_norm )
    Err_value = gt_norm;
  if( Err_value < ht_norm )
    Err_value = ht_norm;

  unsigned int k = 0;

  do {

    // compute error with mu.
    Err_value = norm_2( SES * vect_scalar< ValueType >( K, 1.0 ) + p_grad[range( N, N + K )] ) / sqrt( ValueType( K ) );
    if( Err_value < norm_star )
      Err_value = norm_star;
    if( Err_value < gt_norm )
      Err_value = gt_norm;
    if( Err_value < ht_norm )
      Err_value = ht_norm;

    ValueType abs_tol_mu = mu;
    if( abs_tol_mu < abs_tol )
      abs_tol_mu = abs_tol;

    radius = 0.5 * max_radius;
    ValueType rho = ( 1.0 - mu ) * 0.99;

    while( ( ++k <= max_iter ) && ( Err_value > abs_tol_mu ) ) {
      solve_step( c, Jac_aug, v, norm_v, ValueType( 0.9 ) * radius, abs_tol ); // RK_NOTICE(1," reached");
      for( SizeType i = 0; i < K; ++i )
        if( v[N + i] < -0.5 * tau )
          v[N + i] = -0.5 * tau;
      r = Jac_aug * v;

      try {
        null_space_QP_method( Jac_aug, r, H_aug, p_grad, p, abs_tol, radius ); // RK_NOTICE(1," reached");
      } catch( singularity_error& ) {
        try {
          // RK_NOTICE(4," falling back to PCG method...");
          projected_CG_method( Jac_aug, r, H_aug, p_grad, p, max_iter, abs_tol_mu ); // RK_NOTICE(1," reached");
        } catch( maximum_iteration& ) {
        };
      };

      norm_p = norm_2( p );
      if( norm_p < abs_tol_mu ) {
        if( abs_tol_mu <= abs_tol ) {
          return;
        };
        break;
      };

      if( norm_p > radius ) {
        p *= radius / norm_p;
        norm_p = radius;
      };
      for( SizeType i = 0; i < N; ++i )
        p_x[i] = p[i];
      for( SizeType i = 0; i < K; ++i ) {
        if( p[i + N] < -tau )
          p[i + N] = -tau;
        p_s[i] = s[i] * p[i + N];
      };
      impose_limits( x, p_x );

      // RK_NOTICE(1,"Err_value = " << Err_value << " mu = " << mu << " abs_tol_mu = " << abs_tol_mu << " norm_p = " <<
      // norm_p);
      // RK_NOTICE(1," c = " << c << "\n J v = " << (Jac_aug * v));

      ValueType pHp = p * ( H_aug * p );
      ValueType dm_p = c_norm_star - norm_2( c + Jac_aug * p );
      ValueType dq_p = p_grad * p + ValueType( 0.5 ) * pHp;
      ValueType nu_t;
      if( pHp > ValueType( 0.0 ) )
        nu_t = ( p_grad * p + ValueType( 0.5 ) * pHp ) / ( ( ValueType( 1.0 ) - rho ) * norm_1( c ) );
      else
        nu_t = ( p_grad * p ) / ( ( ValueType( 1.0 ) - rho ) * norm_1( c ) );
      if( nu < nu_t )
        nu = nu_t + abs_tol;
      // else if(nu_t > abs_tol)
      // nu = (nu + nu_t) * ValueType(0.5);
      ValueType pred = nu * dm_p - dq_p;

      xt = x;
      xt += p_x;
      st = s;
      st += p_s;
      ValueType xt_value = f( xt );
      gt_value = g( xt );
      ht_value = h( xt ) - st;
      ValueType log_st( 0.0 );
      for( SizeType i = 0; i < K; ++i )
        log_st += log( st[i] );
      ValueType ct_norm_star = sqrt( norm_2_sqr( gt_value ) + norm_2_sqr( ht_value ) );

      // RK_NOTICE(1," gt_value = " << gt_value);
      // RK_NOTICE(1," x = " << x << " s = " << s << " x_value = " << x_value << " log_s = " << log_s << " c_norm_star =
      // " << c_norm_star << " mu = " << mu);
      // RK_NOTICE(1," xt = " << xt << " st = " << st << " xt_value = " << xt_value << " log_st = " << log_st << "
      // ct_norm_star = " << ct_norm_star);

      ValueType ared = x_value - xt_value - mu * ( log_s - log_st ) + nu * ( c_norm_star - ct_norm_star );

      // RK_NOTICE(1," x_value = " << x_value << " xt_value = " << xt_value << " nu = " << nu << " log_s = " << log_s <<
      // " log_st = " << log_st << " c_norm_star = " << c_norm_star << " ct_norm_star = " << ct_norm_star);
      // RK_NOTICE(1," ared = " << ared << " pred = " << pred << " norm_p = " << norm_2(p_x));

      if( ( pred > ValueType( 0.0 ) ) && ( ared >= kappa * pred ) ) {
        // if(ared >= kappa * pred) {
        // step is accepted.
        x += p_x;
        x_value = xt_value;
        x_grad = df( x );
        p_grad[range( 0, N )] = x_grad;
        s += p_s;
        mS = mat< ValueType, mat_structure::diagonal >( -s );
        log_s = log_st;
        c[range( 0, M )] = gt_value;
        gt_norm = norm_2( gt_value ) / sqrt( ValueType( M ) );
        c[range( M, M + K )] = ht_value;
        ht_norm = norm_2( ht_value ) / sqrt( ValueType( K ) );
        c_norm_star = ct_norm_star;
        fill_g_jac( Jac_g, x, gt_value );
        fill_h_jac( Jac_h, x, ht_value + st );
        if( ( ared >= ValueType( 0.5 ) * pred ) && ( norm_p > ValueType( 0.5 ) * radius ) ) {
          radius *= ValueType( 1.5 );
          if( radius > max_radius )
            radius = max_radius;
        };
        p_grad[range( N, N + K )] = vect_scalar< ValueType >( K, -mu );
        linlsq_QR( transpose_view( Jac_aug ), yz_mat, mat_vect_adaptor< Vector >( p_grad ),
                   abs_tol ); // RK_NOTICE(1," reached");
        //          for(SizeType i = 0; i < K; ++i)
        //            if(z[i] < ValueType(0.0))
        //              z[i] = mu;
        for( SizeType i = 0; i < K; ++i )
          if( z[i] < ValueType( 0.0 ) )
            z[i] = mu / s[i];
        lt = x_grad - y * Jac_g - z * Jac_h;
        // p_grad[range(0,N)] = x_grad + y * Jac_g + z * Jac_h;  // NOTE this is not part of original.
        norm_star = norm_2( lt ) / sqrt( ValueType( N ) );
        fill_hessian( H, x, x_value, x_grad, p_x, lt - l );
        l = lt;
        for( SizeType i = 0; i < K; ++i )
          SES( i, i ) = z[i] * s[i];
        // compute new error with mu.
        Err_value = norm_2( SES * vect_scalar< ValueType >( K, 1.0 ) + p_grad[range( N, N + K )] )
                    / sqrt( ValueType( K ) );
        if( Err_value < norm_star )
          Err_value = norm_star;
        if( Err_value < gt_norm )
          Err_value = gt_norm;
        if( Err_value < ht_norm )
          Err_value = ht_norm;
        // RK_NOTICE(1," Err_value = " << Err_value << " c_norm_star = " << c_norm_star << "\n y = " << y << "\n g = "
        // << gt_value);
      } else {
        if( radius < abs_tol_mu )
          break;
        radius = ValueType( 0.7 ) * norm_p;
      };

      if( norm_2( r ) < abs_tol * sqrt( ValueType( M + K ) ) ) {
        // RK_NOTICE(1," Cannot solve constraints by themselves, trying to get QP to do it..");
        throw infeasible_problem( "Cannot improve on the constraint satisfaction!" );
      };
    };
    if( k > max_iter )
      throw maximum_iteration( max_iter );

    // decrease mu;
    if( K > 1 ) {
      ValueType sz_k = s * z / ValueType( K );
      /*ValueType sigma(2.0);
      ValueType zeta(1.0);
      for(SizeType i = 0; i < K; ++i)
        if( SES(i,i) / sz_k < zeta )
          zeta = SES(i,i) / sz_k;
      if(zeta < 0.5)  //this is just for numerical stability (takes the denominator where it will equalize the
      quantities involved).
        sigma = (ValueType(1.0) - zeta) / (ValueType(20.0) * zeta);
      else
        sigma = ValueType(0.05) * ((ValueType(1.0) - zeta + abs_tol) / zeta);
      if(sigma > ValueType(2.0))
        sigma = ValueType(2.0);
      sigma *= sigma * sigma * ValueType(0.1);
      mu = sigma * sz_k;*/
      mu = sz_k * ValueType( 0.1 );
      // mu *= 0.1;
      // RK_NOTICE(1," inter-step: s = " << s << " z = " << z << " sz_k = " << sz_k << " zeta = " << zeta << " sigma = "
      // << sigma << " mu = " << mu);
    } else {
      mu *= 0.1;
    };

    p_grad[range( N, N + K )] = vect_scalar< ValueType >( K, -mu );
    linlsq_QR( transpose_view( Jac_aug ), yz_mat, mat_vect_adaptor< Vector >( p_grad ), abs_tol );
    //      for(SizeType i = 0; i < K; ++i)
    //        if(z[i] < ValueType(0.0))
    //          z[i] = mu;
    for( SizeType i = 0; i < K; ++i )
      if( z[i] < ValueType( 0.0 ) )
        z[i] = mu / s[i];
    l = x_grad - y * Jac_g - z * Jac_h;
    norm_star = norm_2( l ) / sqrt( ValueType( N ) );
    for( SizeType i = 0; i < K; ++i )
      SES( i, i ) = z[i] * s[i];

    tau = 1.0 - mu;
    // tau = 1.0 - 10.0 * mu;

    // compute error without mu.
    Err_value = norm_2( SES * vect_scalar< ValueType >( K, 1.0 ) ) / sqrt( ValueType( K ) );
    if( Err_value < norm_star )
      Err_value = norm_star;
    if( Err_value < gt_norm )
      Err_value = gt_norm;
    if( Err_value < ht_norm )
      Err_value = ht_norm;
    // RK_NOTICE(1,"Err_value = " << Err_value << " mu = " << mu);
  } while( Err_value > abs_tol );
};


template < typename Function, typename GradFunction, typename Vector, typename EqFunction, typename EqJacFunction,
           typename IneqFunction, typename IneqJacFunction >
struct merit_function_computer {

  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;

  Function f;
  GradFunction df;
  const Vector* x0;
  const Vector* p_x;
  const Vector* s0;
  const Vector* p_s;
  ValueType penalty;
  ValueType mu;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;

  merit_function_computer( Function aF, GradFunction aDF, const Vector& aX0, const Vector& aPX, const Vector& aS0,
                           const Vector& aPS, ValueType aPenalty, ValueType aMu, EqFunction aG, EqJacFunction aFillGJac,
                           IneqFunction aH, IneqJacFunction aFillHJac )
      : f( aF ), df( aDF ), x0( &aX0 ), p_x( &aPX ), s0( &aS0 ), p_s( &aPS ), penalty( aPenalty ), mu( aMu ), g( aG ),
        fill_g_jac( aFillGJac ), h( aH ), fill_h_jac( aFillHJac ){};


  ValueType compute_merit( ValueType alpha_s ) const {
    Vector s = *s0;
    s += alpha_s * ( *p_s );
    Vector x = *x0;
    x += alpha_s * ( *p_x );
    ValueType mulog_s( 0.0 );
    // for(SizeType i = 0; i < s.size(); ++i)
    // mulog_s += mu * log(s[i]);
    return f( x ) - mulog_s + penalty * ( norm_2( g( x ) ) + norm_2( h( x ) - s ) );
  };

  ValueType compute_derivative_merit( ValueType alpha_s ) const {
    Vector s = *s0;
    s += alpha_s * ( *p_s );
    Vector x = *x0;
    x += alpha_s * ( *p_x );
    ValueType muDlog_s( 0.0 );
    // for(SizeType i = 0; i < s.size(); ++i)
    // muDlog_s +=  (mu * (*p_s)[i]) / s[i];
    Vector c_g = g( x );
    mat< ValueType, mat_structure::rectangular > Jac_g( c_g.size(), x.size() );
    fill_g_jac( Jac_g, x, c_g );

    Vector c_h = h( x );
    mat< ValueType, mat_structure::rectangular > Jac_h( c_h.size(), x.size() );
    fill_h_jac( Jac_h, x, c_h );
    c_h -= s;

    ValueType result = ( *p_x ) * df( x );
    if( c_g.size() > 0 ) {
      c_g = unit( c_g );
      result += penalty * ( ( *p_x ) * ( c_g * Jac_g ) );
    };
    if( c_h.size() > 0 ) {
      c_h = unit( c_h );
      result += penalty * ( ( *p_x ) * ( c_h * Jac_h ) - ( *p_s ) * c_h ) - muDlog_s;
    };
    return result;
  };
};


template < typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename EqFunction,
           typename EqJacFunction, typename IneqFunction, typename IneqJacFunction >
void nl_intpoint_method_ls_impl(
  Function f, GradFunction df, HessianFunction fill_hessian, EqFunction g, EqJacFunction fill_g_jac, IneqFunction h,
  IneqJacFunction fill_h_jac, Vector& x,
  typename vect_traits< Vector >::value_type mu = typename vect_traits< Vector >::value_type( 0.1 ),
  unsigned int max_iter = 100,
  typename vect_traits< Vector >::value_type abs_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type kappa = typename vect_traits< Vector >::value_type( 1e-4 ),
  typename vect_traits< Vector >::value_type tau = typename vect_traits< Vector >::value_type( 0.995 ) ) {
  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;
  using std::sqrt;
  using std::fabs;
  using std::log;

  SizeType N = x.size();
  Vector c_h = h( x );
  SizeType K = c_h.size();
  mat< ValueType, mat_structure::rectangular > Jac_h( K, N );
  fill_h_jac( Jac_h, x, c_h );

  // compute initial slack vector and roughly adjust x if needed.
  Vector s = c_h;
  for( SizeType i = 0; i < K; ++i )
    s[i] = ValueType( 1.0 );
  ValueType min_s( 0.0 );
  for( SizeType i = 0; i < K; ++i )
    if( s[i] < min_s )
      min_s = s[i];
  min_s *= ValueType( -0.1 );

  while( min_s > abs_tol ) {
    RK_NOTICE( 1, " s = " << s << " x = " << x << " h = " << c_h << " Jac = " << Jac_h );
    Vector e_s( s );
    for( SizeType i = 0; i < K; ++i )
      e_s[i] = min_s;
    s += e_s;
    Vector dx = x;
    mat_vect_adaptor< Vector > dx_mat( dx );
    minnorm_QR( Jac_h, dx_mat, mat_vect_adaptor< Vector >( e_s ), abs_tol );
    RK_NOTICE( 1, " e_s = " << e_s << " dx = " << dx );
    x += dx;
    c_h = h( x );
    fill_h_jac( Jac_h, x, c_h );
    s = c_h;
    min_s = ValueType( 0.0 );
    for( SizeType i = 0; i < K; ++i )
      if( s[i] < min_s )
        min_s = s[i];
    min_s *= ValueType( -0.1 );
  };
  if( min_s > ValueType( 0.0 ) )
    s += vect_scalar< double >( K, ValueType( 1000.0 ) * abs_tol );

  Vector c_g = g( x );
  SizeType M = c_g.size();
  ValueType g_norm = norm_2( c_g );

  mat< ValueType, mat_structure::rectangular > Jac_g( M, N );
  fill_g_jac( Jac_g, x, c_g );

  if( ( M == 0 ) && ( K == 0 ) ) { // this means it is an unconstrained problem. TODO change this to dispatch on the
                                   // type of the fill-hessian functor.
    newton_method_ls_impl( f, df, fill_hessian, x, max_iter,
                           line_search_expand_and_zoom< ValueType >( kappa, ValueType( 0.9 ) ), no_limit_functor(),
                           newton_directioner(), abs_tol, abs_tol );
    return;
  };

  ValueType h_norm = norm_2( c_h - s );

  ValueType x_value = f( x );
  Vector x_grad = df( x );
  mat< ValueType, mat_structure::symmetric > H( mat< ValueType, mat_structure::identity >( x.size() ) );
  fill_hessian( H, x, x_value, x_grad );

  Vector y = c_g;
  Vector z = c_h;
  {
    mat< ValueType, mat_structure::diagonal > mS( K );
    mat< ValueType, mat_structure::nil > Zero_mk( M, K );

    mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >, mat< ValueType, mat_structure::nil > >
      Jac_upper( Jac_g, Zero_mk );
    mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >, mat< ValueType, mat_structure::diagonal > >
      Jac_lower( Jac_h, mS );

    mat_const_ref_vert_cat< mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >,
                                                     mat< ValueType, mat_structure::nil > >,
                            mat_const_ref_horiz_cat< mat< ValueType, mat_structure::rectangular >,
                                                     mat< ValueType, mat_structure::diagonal > > > Jac_aug( Jac_upper,
                                                                                                            Jac_lower );

    mS = mat< ValueType, mat_structure::diagonal >( -s );

    Vector yz;
    yz.resize( M + K );
    mat_vect_adaptor< Vector > yz_mat( yz );
    Vector p_grad;
    p_grad.resize( N + K );
    p_grad[range( 0, N )] = x_grad;
    p_grad[range( N, N + K )] = vect_scalar< ValueType >( K, -mu );
    linlsq_QR( transpose_view( Jac_aug ), yz_mat, mat_vect_adaptor< Vector >( p_grad ), abs_tol );
    y = yz[range( 0, M )];
    z = yz[range( M, M + K )];
  };


  Vector xt = x;
  Vector st = s;
  Vector p_x = x;
  Vector dx = x;
  ValueType norm_p = std::numeric_limits< ValueType >::max();
  Vector p_s = s;
  Vector p_y = y;
  Vector p_z = z;

  ValueType log_s( 0.0 );
  for( SizeType i = 0; i < K; ++i )
    log_s += log( s[i] );


  // ValueType nu(std::numeric_limits<ValueType>::min());
  ValueType nu( 0.1 );


  for( SizeType i = 0; i < K; ++i )
    if( z[i] < abs_tol )
      z[i] = abs_tol;

  ValueType norm_star = norm_2( x_grad ) + norm_2( y * Jac_g ) + norm_2( z * Jac_h );
  Vector l( x_grad - y * Jac_g - z * Jac_h );
  Vector lt = l;
  norm_star = norm_2( l );

  mat< ValueType, mat_structure::diagonal > muSES_inv( K );
  for( SizeType i = 0; i < K; ++i )
    muSES_inv( i, i ) = mu / ( z[i] * s[i] );

  mat< ValueType, mat_structure::rectangular > ZJac_h( Jac_h );
  mat< ValueType, mat_structure::rectangular > SJac_h( Jac_h );
  for( SizeType i = 0; i < K; ++i )
    for( SizeType j = 0; j < N; ++j ) {
      ZJac_h( i, j ) *= z[i];
      SJac_h( i, j ) /= s[i];
    };

  mat< ValueType, mat_structure::symmetric > qp_G( H + transpose_view( SJac_h ) * ZJac_h );

  Vector qp_c( x_grad );
  Vector c_h_s = c_h;
  for( SizeType i = 0; i < K; ++i )
    c_h_s[i] /= s[i];
  qp_c -= y * Jac_g + z * Jac_h - ( c_h_s - muSES_inv * vect_scalar< ValueType >( K, 1.0 ) ) * ZJac_h;

  typedef merit_function_computer< Function, GradFunction, Vector, EqFunction, EqJacFunction, IneqFunction,
                                   IneqJacFunction > MeritFuncComputer;
  MeritFuncComputer m_func( f, df, x, p_x, s, p_s, nu, mu, g, fill_g_jac, h, fill_h_jac );


  ValueType Err_value = ValueType( 0.0 );
  for( SizeType i = 0; i < K; ++i ) {
    ValueType tmp = z[i] * s[i];
    if( fabs( tmp ) > Err_value )
      Err_value = fabs( tmp );
  };
  if( Err_value < norm_star )
    Err_value = norm_star;
  if( Err_value < g_norm )
    Err_value = g_norm;
  if( Err_value < h_norm )
    Err_value = h_norm;

  unsigned int k = 0;

  do {

    // compute error with mu.
    Err_value = ValueType( 0.0 );
    for( SizeType i = 0; i < K; ++i ) {
      ValueType tmp = z[i] * s[i] - mu;
      if( fabs( tmp ) > Err_value )
        Err_value = fabs( tmp );
    };
    if( Err_value < norm_star )
      Err_value = norm_star;
    if( Err_value < g_norm )
      Err_value = g_norm;
    if( Err_value < h_norm )
      Err_value = h_norm;

    ValueType abs_tol_mu = mu;
    // if(abs_tol_mu < abs_tol)
    // abs_tol_mu = abs_tol;

    ValueType rho = ValueType( 0.1 );

    while( ( ++k <= max_iter ) && ( Err_value > abs_tol_mu ) ) {

      //         Vector c_aug; c_aug.resize(N + K);
      //         c_aug[range(0,N)] = l;
      //         ValueType mu_sqrt = sqrt(mu);
      //         for(SizeType i = 0; i < K; ++i)
      //           c_aug[i+N] = s[i] * z[i] / mu_sqrt - mu_sqrt;
      //         Vector b_aug; b_aug.resize(M + K);
      //         b_aug[range(0,M)] = -c_g;
      //         b_aug[range(M,M+K)] = s - c_h;
      //
      //         mat<ValueType,mat_structure::rectangular> J_aug(M+K,N+K);
      //         sub(J_aug)(range(0,M),range(0,N)) = Jac_g;
      //         sub(J_aug)(range(M,M+K),range(0,N)) = Jac_h;
      //         sub(J_aug)(range(M,M+K),range(N,N+K)) = (ValueType(-1.0) / mu_sqrt) *
      //         mat<ValueType,mat_structure::diagonal>(s);
      //
      //         mat<ValueType,mat_structure::diagonal> E_aug(K);
      //         for(SizeType i = 0; i < K; ++i)
      //           E_aug(i,i) = (z[i] * s[i]) / mu;
      //         mat<ValueType,mat_structure::symmetric> H_aug( ((H + mat<ValueType,mat_structure::scalar>(N,abs_tol)) &
      //         mat<ValueType,mat_structure::nil>(N,K)) | (mat<ValueType,mat_structure::nil>(K,N) &
      //         mat<ValueType,mat_structure::identity>(K) ) );
      //
      //         Vector p_aug; p_aug.resize(N+K);
      //         Vector l_aug; l_aug.resize(M+K);
      //
      //         mat<ValueType,mat_structure::rectangular> rhs_mat(N+M+K+K,1);
      //         for(SizeType i = 0; i < N+K; ++i)
      //           rhs_mat(i,0) = -c_aug[i];
      //         for(SizeType i = 0; i < M+K; ++i)
      //           rhs_mat(i+N+K,0) = b_aug[i];
      //         mat<ValueType,mat_structure::rectangular> lhs_mat(N+M+K+K,1);
      //         mat<ValueType,mat_structure::rectangular> H_aug2( (H_aug & transpose_view(J_aug)) | (J_aug &
      //         mat<ValueType,mat_structure::scalar>(M+K,-abs_tol) ) );
      //         RK_NOTICE(1," reached");
      //         QR_linlsqsolver()(H_aug2,lhs_mat,rhs_mat,abs_tol);
      //         RK_NOTICE(1," reached");
      //         for(SizeType i = 0; i < N; ++i)
      //           p_x[i] = lhs_mat(i,0);
      //         for(SizeType i = 0; i < K; ++i)
      //           p_s[i] = lhs_mat(i+N,0);
      //         for(SizeType i = 0; i < M; ++i)
      //           p_y[i] = -lhs_mat(i+N+K,0);
      //         for(SizeType i = 0; i < K; ++i)
      //           p_z[i] = -lhs_mat(i+N+K+M,0);

      //         try {
      //           null_space_QP_method(J_aug, b_aug, H_aug, c_aug, p_aug, abs_tol,
      //           std::numeric_limits<ValueType>::infinity(), &l_aug);
      //         } catch(singularity_error&) {
      //           try {
      //             projected_CG_method(J_aug, b_aug, H_aug, c_aug, p_aug, max_iter, abs_tol_mu, &l_aug);
      //           } catch(maximum_iteration&) { };
      //         };
      //         p_x = p_aug[range(0,N)];
      //         p_s = p_aug[range(N,N+K)];
      //         p_y = l_aug[range(0,M)];
      //         p_z = l_aug[range(M,M+K)];
      //
      //        for(SizeType i = 0; i < K; ++i)
      //          p_s[i] = (mu_sqrt * p_s[i]) / s[i];


      try {
        null_space_QP_method( Jac_g, -c_g, qp_G, qp_c, p_x, abs_tol, std::numeric_limits< ValueType >::infinity(),
                              &p_y );
      } catch( singularity_error& ) {
        try {
          projected_CG_method( Jac_g, -c_g, qp_G, qp_c, p_x, max_iter, abs_tol_mu, &p_y );
        } catch( maximum_iteration& ) {
        };
      };

      //         p_s = Jac_h * p_x + c_h - s;
      ValueType alpha_s_max( 1.0 );
      ValueType dq_p( 0.0 );
      ValueType pHp = p_x * ( H * p_x );
      for( SizeType i = 0; i < K; ++i ) {
        if( alpha_s_max * p_s[i] < -tau * s[i] )
          alpha_s_max = -tau * ( s[i] / p_s[i] );
        dq_p -= mu * p_s[i] - z[i] * ( p_s[i] / s[i] );
        pHp += p_s[i] * z[i] * ( p_s[i] / s[i] );
      };
      dq_p += l * p_x;

      m_func.penalty = nu;
      m_func.mu = mu;
      ValueType alpha_s( 0.0 );

      // RK_NOTICE(1," nu = " << nu << " alpha_max = " << alpha_s_max << " phi(0) = " << m_func.compute_merit(0.0) << "
      // phi'(0) = " << m_func.compute_derivative_merit(0.0) << " phi(a) = " << m_func.compute_merit(alpha_s_max) << "
      // phi'(a) = " << m_func.compute_derivative_merit(alpha_s_max));

      if( m_func.compute_derivative_merit( alpha_s_max ) > ValueType( 0.0 ) ) {

        //           expand_and_zoom_search([&m_func](ValueType alpha_s) -> ValueType { return
        //           m_func.compute_merit(alpha_s); },
        //                                  [&m_func](ValueType alpha_s) -> ValueType { return
        //                                  m_func.compute_derivative_merit(alpha_s); },
        //                                   alpha_s, alpha_s_max, abs_tol_mu, kappa, ValueType(0.1));
        backtracking_search(
          [&m_func]( ValueType alpha_s ) -> ValueType { return m_func.compute_merit( alpha_s ); },
          [&m_func]( ValueType alpha_s ) -> ValueType { return m_func.compute_derivative_merit( alpha_s ); }, alpha_s,
          alpha_s_max, 0.0001, kappa, ValueType( 0.5 ), ValueType( 0.75 ) );
      } else {
        alpha_s = alpha_s_max;
      };

      // RK_NOTICE(1," alpha_s = " << alpha_s << " phi(0) = " << m_func.compute_merit(0.0) << " phi'(0) = " <<
      // m_func.compute_derivative_merit(0.0) << " phi(a) = " << m_func.compute_merit(alpha_s) << " phi'(a) = " <<
      // m_func.compute_derivative_merit(alpha_s));

      dx = alpha_s * p_x;
      norm_p = norm_2( dx );
      x += dx;
      s += alpha_s * p_s;

      pHp *= alpha_s * alpha_s;
      dq_p *= alpha_s;
      ValueType nu_t;
      if( pHp > ValueType( 0.0 ) )
        nu_t = ( dq_p + ValueType( 0.5 ) * pHp )
               / ( ( ValueType( 1.0 ) - rho ) * ( norm_1( c_g ) + norm_1( c_h - s ) ) );
      else
        nu_t = dq_p / ( ( ValueType( 1.0 ) - rho ) * ( norm_1( c_g ) + norm_1( c_h - s ) ) );
      if( nu < nu_t )
        nu = nu_t + abs_tol;
      else
        nu *= 1.1;

      //         p_z = muSES_inv * vect_scalar<ValueType>(K,ValueType(1.0)) - c_h_s - SJac_h * p_x;
      ValueType alpha_z( ( mu * ValueType( K ) - s * z ) / ( s * p_z ) );
      for( SizeType i = 0; i < K; ++i ) {
        if( alpha_z * p_z[i] < -tau )
          alpha_z = -tau / p_z[i];
        p_z[i] = z[i] * p_z[i];
      };
      // RK_NOTICE(1," alpha_s = " << alpha_s << " p_x = " << p_x << " norm_p = " << norm_2(alpha_s * p_x) << " p_s = "
      // << p_s << " p_y = " << p_y << " p_z = " << p_z);

      y += alpha_z * p_y;
      z += alpha_z * p_z;

      c_g = g( x );
      g_norm = norm_2( c_g );
      fill_g_jac( Jac_g, x, c_g );
      c_h = h( x );
      fill_h_jac( Jac_h, x, c_h );
      h_norm = norm_2( c_h - s );
      ZJac_h = Jac_h;
      for( SizeType i = 0; i < K; ++i )
        for( SizeType j = 0; j < N; ++j ) {
          ZJac_h( i, j ) *= z[i];
          SJac_h( i, j ) /= s[i];
        };
      for( SizeType i = 0; i < K; ++i )
        muSES_inv( i, i ) = mu / ( z[i] * s[i] );

      x_value = f( x );
      x_grad = df( x );
      lt = x_grad - y * Jac_g - z * Jac_h;

      for( SizeType i = 0; i < K; ++i )
        c_h_s[i] = c_h[i] / s[i];

      qp_c = lt + ( c_h_s - muSES_inv * vect_scalar< ValueType >( K, 1.0 ) ) * ZJac_h;

      fill_hessian( H, x, x_value, x_grad, dx, lt - l );
      l = lt;
      norm_star = norm_2( l );

      qp_G = H + transpose_view( SJac_h ) * ZJac_h;

      Err_value = ValueType( 0.0 );
      for( SizeType i = 0; i < K; ++i ) {
        ValueType tmp = z[i] * s[i] - mu;
        if( fabs( tmp ) > Err_value )
          Err_value = fabs( tmp );
      };
      if( Err_value < norm_star )
        Err_value = norm_star;
      if( Err_value < g_norm )
        Err_value = g_norm;
      if( Err_value < h_norm )
        Err_value = h_norm;
      // RK_NOTICE(1," Err_value = " << Err_value << " g_norm = " << g_norm << " h_norm = " << h_norm << " s = " << s <<
      // " z = " << z << " y = " << y);

      if( norm_p < abs_tol )
        break;
    };
    if( k > max_iter )
      throw maximum_iteration( max_iter );
    // if(radius < abs_tol)
    // return;

    if( ( abs_tol_mu <= abs_tol ) && ( norm_p < abs_tol ) )
      return;

    // decrease mu;
    if( K > 1 ) {
      ValueType sigma( 2.0 );
      ValueType zeta( 1.0 );
      ValueType sz_k = s * z / ValueType( K );
      for( SizeType i = 0; i < K; ++i )
        if( ( z[i] * s[i] ) < zeta * sz_k )
          zeta = ( z[i] * s[i] ) / sz_k;
      if( zeta < 0.5 ) // this is just for numerical stability (takes the denominator where it will equalize the
                       // quantities involved).
        sigma = ( ValueType( 1.0 ) - zeta ) / ( ValueType( 20.0 ) * zeta );
      else
        sigma = ValueType( 0.05 ) * ( ( ValueType( 1.0 ) - zeta + abs_tol ) / zeta );
      if( sigma > ValueType( 2.0 ) )
        sigma = ValueType( 2.0 );
      sigma *= sigma * sigma * ValueType( 0.1 );
      mu = sigma * sz_k;
      // RK_NOTICE(1," inter-step: s = " << s << " z = " << z << " sz_k = " << sz_k << " zeta = " << zeta << " sigma = "
      // << sigma << " mu = " << mu);
    } else {
      mu *= 0.5;
    };

    for( SizeType i = 0; i < K; ++i )
      muSES_inv( i, i ) = mu / ( z[i] * s[i] );
    qp_c = l + ( c_h_s - muSES_inv * vect_scalar< ValueType >( K, 1.0 ) ) * ZJac_h;

    // tau += ValueType(0.5) * (ValueType(1.0) - tau);

    // compute error without mu.
    Err_value = ValueType( 0.0 );
    for( SizeType i = 0; i < K; ++i ) {
      ValueType tmp = z[i] * s[i];
      if( fabs( tmp ) > Err_value )
        Err_value = fabs( tmp );
    };
    if( Err_value < norm_star )
      Err_value = norm_star;
    if( Err_value < g_norm )
      Err_value = g_norm;
    if( Err_value < h_norm )
      Err_value = h_norm;
    // RK_NOTICE(1,"Err_value = " << Err_value);
  } while( Err_value > abs_tol );
};
};


/**
 * This functor is a factory class to construct a Non-Linear Interior-Point optimizer routine
 * that uses a trust-region approach and incorporates equality and inequality constraints via
 * the sequential quadratic programming method and a Newton direction. Use make_nlip_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see
 * trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 */
template < typename Function, typename GradFunction, typename HessianFunction, typename T,
           typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor,
           typename IneqFunction = no_constraint_functor, typename IneqJacFunction = no_constraint_jac_functor,
           typename TrustRegionSolver = tr_solver_right_pinv_dogleg, typename LimitFunction = no_limit_functor >
struct nlip_newton_tr_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T max_radius;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;

  typedef nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                  IneqJacFunction, TrustRegionSolver, LimitFunction > self;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aG The function that computes the inequality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the inequality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or
   * non-negativity).
   */
  nlip_newton_tr_factory( Function aF, GradFunction aDf, HessianFunction aFillHessian, T aMaxRadius, T aMu,
                          unsigned int aMaxIter, EqFunction aG = EqFunction(),
                          EqJacFunction aFillGJac = EqJacFunction(), IneqFunction aH = EqFunction(),
                          IneqJacFunction aFillHJac = IneqJacFunction(), T aTol = T( 1e-6 ), T aEta = T( 1e-4 ),
                          T aTau = T( 0.995 ), TrustRegionSolver aSolveStep = TrustRegionSolver(),
                          LimitFunction aImposeLimits = LimitFunction() )
      : f( aF ), df( aDf ), fill_hessian( aFillHessian ), max_radius( aMaxRadius ), mu( aMu ), max_iter( aMaxIter ),
        tol( aTol ), eta( aEta ), tau( aTau ), g( aG ), fill_g_jac( aFillGJac ), h( aH ), fill_h_jac( aFillHJac ),
        solve_step( aSolveStep ), impose_limits( aImposeLimits ){};
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template < typename Vector >
  void operator()( Vector& x ) const {
    detail::nl_intpoint_method_tr_impl( f, df, hessian_update_dual_exact< HessianFunction >( fill_hessian ), g,
                                        fill_g_jac, h, fill_h_jac, x, max_radius, mu, max_iter, solve_step,
                                        impose_limits, tol, eta, tau );
  };

  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius( T aMaxRadius ) {
    max_radius = aMaxRadius;
    return *this;
  };
  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"),
   * this parameter is positive and should start with a rather large value (relative to the scale
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param( T aMu ) {
    mu = aMu;
    return *this;
  };
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration( unsigned int aMaxIter ) {
    max_iter = aMaxIter;
    return *this;
  };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance( T aTol ) {
    tol = aTol;
    return *this;
  };
  /**
   * Sets the tolerance on the decrease in order to accept a step in the trust region.
   */
  self& set_decrease_tolerance( T aEta ) {
    eta = aEta;
    return *this;
  };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality
   * constraint (barrier).
   */
  self& set_barrier_step_margin( T aTau ) {
    tau = aTau;
    return *this;
  };

  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param aTau The initial relative damping factor to regularize the Hessian matrix.
   */
  nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                          IneqJacFunction, tr_solver_right_pinv_dogleg_reg< T >, LimitFunction >
    regularize( const T& aTau ) const {
    return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                   IneqJacFunction, tr_solver_right_pinv_dogleg_reg< T >, LimitFunction >(
      f, df, fill_hessian, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau,
      tr_solver_right_pinv_dogleg_reg< T >( aTau ), impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see
   * trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template < typename NewTrustRegionSolver >
  nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                          IneqJacFunction, NewTrustRegionSolver, LimitFunction >
    set_tr_solver( NewTrustRegionSolver new_solver ) const {
    return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                   IneqJacFunction, NewTrustRegionSolver, LimitFunction >(
      f, df, fill_hessian, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, new_solver,
      impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given limit-function for the search domain.
   * Using a limit-function boils down to a gradient projection method, and thus, for more complex
   * constraints (such as non-linear ones or mix of equality-inequality constraints), please use a
   * constraint optimization method instead (see augmented_lagrangian_methods.hpp for example).
   * \tparam NewLimitFunction A new functor type that can impose limits on a proposed solution step (see
   * no_limit_functor or box_limit_function for examples).
   * \param new_limits The functor that can impose simple limits on the search domain (i.e. using this boils down to a
   * gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  template < typename NewLimitFunction >
  nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                          IneqJacFunction, TrustRegionSolver, NewLimitFunction >
    set_limiter( NewLimitFunction new_limits ) const {
    return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                   IneqJacFunction, TrustRegionSolver, NewLimitFunction >(
      f, df, fill_hessian, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, solve_step,
      new_limits );
  };

  /**
   * This function remaps the factory to one which will use the given equality constraints for
   * the search domain. The equality constraints must be formulated at g(x) = 0.
   * \tparam NewEqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewEqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_g The function that computes the equality constraints.
   * \param new_fill_g_jac The function that computes the jacobian matrix of the equality constaints.
   */
  template < typename NewEqFunction, typename NewEqJacFunction >
  nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                          IneqJacFunction, TrustRegionSolver, LimitFunction >
    set_eq_constraints( NewEqFunction new_g, NewEqJacFunction new_fill_g_jac ) const {
    return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, NewEqFunction, NewEqJacFunction,
                                   IneqFunction, IneqJacFunction, TrustRegionSolver, LimitFunction >(
      f, df, fill_hessian, max_radius, mu, max_iter, new_g, new_fill_g_jac, h, fill_h_jac, tol, eta, tau, solve_step,
      impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given inequality constraints for
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template < typename NewIneqFunction, typename NewIneqJacFunction >
  nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                          NewIneqJacFunction, TrustRegionSolver, LimitFunction >
    set_ineq_constraints( NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac ) const {
    return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction,
                                   NewIneqFunction, NewIneqJacFunction, TrustRegionSolver, LimitFunction >(
      f, df, fill_hessian, max_radius, mu, max_iter, g, fill_g_jac, new_h, new_fill_h_jac, tol, eta, tau, solve_step,
      impose_limits );
  };
};

/**
 * This function template creates a factory object to construct a Non-Linear
 * Interior-Point optimizer routine that uses a trust-region approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method and a Newton steps. This algorithm is roughly as described
 * in Nocedal's Numerical Optimization book. This algorithm solves the following
 * problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
 * positive and should start with a rather large value (relative to the scale of the function) and will be progressively
 * decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
 * (barrier).
 */
template < typename Function, typename GradFunction, typename HessianFunction, typename T >
nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T >
  make_nlip_newton_tr( Function f, GradFunction df, HessianFunction fill_hessian, T max_radius, T mu,
                       unsigned int max_iter, T tol = T( 1e-6 ), T eta = T( 1e-4 ), T tau = T( 0.995 ) ) {
  return nlip_newton_tr_factory< Function, GradFunction, HessianFunction, T >(
    f, df, fill_hessian, max_radius, mu, max_iter, no_constraint_functor(), no_constraint_jac_functor(),
    no_constraint_functor(), no_constraint_jac_functor(), tol, eta, tau );
};


/**
 * This functor is a factory class to construct a Non-Linear
 * Interior-Point optimizer routine that uses a trust-region approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update). Use make_bosqp_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \tparam HessianUpdater The functor type to update the Hessian approximation of the function to optimize.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see
 * trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 */
template < typename Function, typename GradFunction, typename T, typename EqFunction = no_constraint_functor,
           typename EqJacFunction = no_constraint_jac_functor, typename IneqFunction = no_constraint_functor,
           typename IneqJacFunction = no_constraint_jac_functor, typename HessianUpdater = hessian_update_sr1,
           typename TrustRegionSolver = tr_solver_right_pinv_dogleg, typename LimitFunction = no_limit_functor >
struct nlip_quasi_newton_tr_factory {
  Function f;
  GradFunction df;
  T max_radius;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  HessianUpdater update_hessian;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;

  typedef nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                        IneqJacFunction, HessianUpdater, TrustRegionSolver, LimitFunction > self;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   * \param aUpdateHessian The functor object that can update the approximate Hessian matrix of the function to be
   * optimized.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or
   * non-negativity).
   */
  nlip_quasi_newton_tr_factory( Function aF, GradFunction aDf, T aMaxRadius, T aMu, unsigned int aMaxIter,
                                EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
                                IneqFunction aH = IneqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
                                T aTol = T( 1e-6 ), T aEta = T( 1e-4 ), T aTau = T( 0.995 ),
                                HessianUpdater aUpdateHessian = HessianUpdater(),
                                TrustRegionSolver aSolveStep = TrustRegionSolver(),
                                LimitFunction aImposeLimits = LimitFunction() )
      : f( aF ), df( aDf ), max_radius( aMaxRadius ), mu( aMu ), max_iter( aMaxIter ), tol( aTol ), eta( aEta ),
        tau( aTau ), g( aG ), fill_g_jac( aFillGJac ), h( aH ), fill_h_jac( aFillHJac ),
        update_hessian( aUpdateHessian ), solve_step( aSolveStep ), impose_limits( aImposeLimits ){};
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template < typename Vector >
  void operator()( Vector& x ) const {
    detail::nl_intpoint_method_tr_impl( f, df, hessian_update_dual_quasi< HessianUpdater >( update_hessian ), g,
                                        fill_g_jac, h, fill_h_jac, x, max_radius, mu, max_iter, solve_step,
                                        impose_limits, tol, eta, tau );
  };

  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius( T aMaxRadius ) {
    max_radius = aMaxRadius;
    return *this;
  };
  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"),
   * this parameter is positive and should start with a rather large value (relative to the scale
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param( T aMu ) {
    mu = aMu;
    return *this;
  };
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration( unsigned int aMaxIter ) {
    max_iter = aMaxIter;
    return *this;
  };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance( T aTol ) {
    tol = aTol;
    return *this;
  };
  /**
   * Sets the tolerance on the decrease in order to accept a step in the trust region.
   */
  self& set_decrease_tolerance( T aEta ) {
    eta = aEta;
    return *this;
  };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality
   * constraint (barrier).
   */
  self& set_barrier_step_margin( T aTau ) {
    tau = aTau;
    return *this;
  };

  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                HessianUpdater, tr_solver_right_pinv_dogleg_reg< T >, LimitFunction >
    regularize( const T& tau ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                         IneqJacFunction, HessianUpdater, tr_solver_right_pinv_dogleg_reg< T >,
                                         LimitFunction >( f, df, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac,
                                                          tol, eta, tau, update_hessian,
                                                          tr_solver_right_pinv_dogleg_reg< T >( tau ), impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see
   * trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template < typename NewTrustRegionSolver >
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                HessianUpdater, NewTrustRegionSolver, LimitFunction >
    set_tr_solver( NewTrustRegionSolver new_solver ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                         IneqJacFunction, HessianUpdater, NewTrustRegionSolver, LimitFunction >(
      f, df, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, update_hessian, new_solver,
      impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given limit-function for the search domain.
   * Using a limit-function boils down to a gradient projection method, and thus, for more complex
   * constraints (such as non-linear ones or mix of equality-inequality constraints), please use a
   * constraint optimization method instead (see augmented_lagrangian_methods.hpp for example).
   * \tparam NewLimitFunction A new functor type that can impose limits on a proposed solution step (see
   * no_limit_functor or box_limit_function for examples).
   * \param new_limits The functor that can impose simple limits on the search domain (i.e. using this boils down to a
   * gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  template < typename NewLimitFunction >
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                HessianUpdater, TrustRegionSolver, NewLimitFunction >
    set_limiter( NewLimitFunction new_limits ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                         IneqJacFunction, HessianUpdater, TrustRegionSolver, NewLimitFunction >(
      f, df, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, update_hessian, solve_step,
      new_limits );
  };

  /**
   * This function remaps the factory to one which will use the given equality constraints for
   * the search domain. The equality constraints must be formulated at g(x) = 0.
   * \tparam NewEqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewEqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_g The function that computes the equality constraints.
   * \param new_fill_g_jac The function that computes the jacobian matrix of the equality constaints.
   */
  template < typename NewEqFunction, typename NewEqJacFunction >
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                                IneqJacFunction, HessianUpdater, TrustRegionSolver, LimitFunction >
    set_eq_constraints( NewEqFunction new_g, NewEqJacFunction new_fill_g_jac ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                                         IneqJacFunction, HessianUpdater, TrustRegionSolver, LimitFunction >(
      f, df, max_radius, mu, max_iter, new_g, new_fill_g_jac, h, fill_h_jac, tol, eta, tau, update_hessian, solve_step,
      impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given inequality constraints for
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template < typename NewIneqFunction, typename NewIneqJacFunction >
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                                NewIneqJacFunction, HessianUpdater, TrustRegionSolver, LimitFunction >
    set_ineq_constraints( NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                                         NewIneqJacFunction, HessianUpdater, TrustRegionSolver, LimitFunction >(
      f, df, max_radius, mu, max_iter, g, fill_g_jac, new_h, new_fill_h_jac, tol, eta, tau, update_hessian, solve_step,
      impose_limits );
  };

  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewHessianUpdater A new functor type to update the Hessian approximation of the function to optimize.
   * \param new_update_hessian The new functor object that can update the approximate Hessian matrix of the function to
   * be optimized.
   */
  template < typename NewHessianUpdater >
  nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                NewHessianUpdater, TrustRegionSolver, LimitFunction >
    set_hessian_updater( NewHessianUpdater new_update_hessian ) const {
    return nlip_quasi_newton_tr_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                         IneqJacFunction, NewHessianUpdater, TrustRegionSolver, LimitFunction >(
      f, df, max_radius, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, new_update_hessian, solve_step,
      impose_limits );
  };
};

/**
 * This function template creates a factory object to construct a Non-Linear
 * Interior-Point optimizer routine that uses a trust-region approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update).
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
 * positive and should start with a rather large value (relative to the scale of the function) and will be progressively
 * decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
 * (barrier).
 */
template < typename Function, typename GradFunction, typename T >
nlip_quasi_newton_tr_factory< Function, GradFunction, T >
  make_nlip_quasi_newton_tr( Function f, GradFunction df, T max_radius, T mu, unsigned int max_iter, T tol = T( 1e-6 ),
                             T eta = T( 1e-4 ), T tau = T( 0.995 ) ) {
  return nlip_quasi_newton_tr_factory< Function, GradFunction, T >(
    f, df, max_radius, mu, max_iter, no_constraint_functor(), no_constraint_jac_functor(), no_constraint_functor(),
    no_constraint_jac_functor(), tol, eta, tau );
};


/**
 * This functor is a factory class to construct a Non-Linear Interior-Point optimizer routine
 * that uses a line-search approach and incorporates equality and inequality constraints via
 * the sequential quadratic programming method and a Newton direction. Use make_nlip_newton_ls to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this. So far, all tests have failed, there must be something really wrong in this
 * implementation, but I cannot pin-point it.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 */
template < typename Function, typename GradFunction, typename HessianFunction, typename T,
           typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor,
           typename IneqFunction = no_constraint_functor, typename IneqJacFunction = no_constraint_jac_functor >
struct nlip_newton_ls_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;

  typedef nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                  IneqJacFunction > self;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aG The function that computes the inequality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the inequality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the sufficient decrease in order to accept a line-search step.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  nlip_newton_ls_factory( Function aF, GradFunction aDf, HessianFunction aFillHessian, T aMu, unsigned int aMaxIter,
                          EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
                          IneqFunction aH = EqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
                          T aTol = T( 1e-6 ), T aEta = T( 1e-4 ), T aTau = T( 0.995 ) )
      : f( aF ), df( aDf ), fill_hessian( aFillHessian ), mu( aMu ), max_iter( aMaxIter ), tol( aTol ), eta( aEta ),
        tau( aTau ), g( aG ), fill_g_jac( aFillGJac ), h( aH ), fill_h_jac( aFillHJac ){};
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template < typename Vector >
  void operator()( Vector& x ) const {
    detail::nl_intpoint_method_ls_impl( f, df, hessian_update_dual_exact< HessianFunction >( fill_hessian ), g,
                                        fill_g_jac, h, fill_h_jac, x, mu, max_iter, tol, eta, tau );
  };

  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"),
   * this parameter is positive and should start with a rather large value (relative to the scale
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param( T aMu ) {
    mu = aMu;
    return *this;
  };
  /**
   * Sets the maximum number of iterations to perform.
   */
  self& set_max_iteration( unsigned int aMaxIter ) {
    max_iter = aMaxIter;
    return *this;
  };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance( T aTol ) {
    tol = aTol;
    return *this;
  };
  /**
   * Sets the tolerance on the sufficient decrease in order to accept a line-search step.
   */
  self& set_decrease_tolerance( T aEta ) {
    eta = aEta;
    return *this;
  };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality
   * constraint (barrier).
   */
  self& set_barrier_step_margin( T aTau ) {
    tau = aTau;
    return *this;
  };


  /**
   * This function remaps the factory to one which will use the given equality constraints for
   * the search domain. The equality constraints must be formulated at g(x) = 0.
   * \tparam NewEqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewEqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_g The function that computes the equality constraints.
   * \param new_fill_g_jac The function that computes the jacobian matrix of the equality constaints.
   */
  template < typename NewEqFunction, typename NewEqJacFunction >
  nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                          IneqJacFunction >
    set_eq_constraints( NewEqFunction new_g, NewEqJacFunction new_fill_g_jac ) const {
    return nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T, NewEqFunction, NewEqJacFunction,
                                   IneqFunction, IneqJacFunction >( f, df, fill_hessian, mu, max_iter, new_g,
                                                                    new_fill_g_jac, h, fill_h_jac, tol, eta, tau );
  };

  /**
   * This function remaps the factory to one which will use the given inequality constraints for
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template < typename NewIneqFunction, typename NewIneqJacFunction >
  nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                          NewIneqJacFunction >
    set_ineq_constraints( NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac ) const {
    return nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T, EqFunction, EqJacFunction,
                                   NewIneqFunction, NewIneqJacFunction >(
      f, df, fill_hessian, mu, max_iter, g, fill_g_jac, new_h, new_fill_h_jac, tol, eta, tau );
  };
};

/**
 * This function template creates a factory object to construct a Non-Linear
 * Interior-Point optimizer routine that uses a line-search approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method and a Newton steps. This algorithm is roughly as described
 * in Nocedal's Numerical Optimization book. This algorithm solves the following
 * problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this. So far, all tests have failed, there must be something really wrong in this
 * implementation, but I cannot pin-point it.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
 * positive and should start with a rather large value (relative to the scale of the function) and will be progressively
 * decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
 * (barrier).
 */
template < typename Function, typename GradFunction, typename HessianFunction, typename T >
nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T >
  make_nlip_newton_ls( Function f, GradFunction df, HessianFunction fill_hessian, T mu, unsigned int max_iter,
                       T tol = T( 1e-6 ), T eta = T( 1e-4 ), T tau = T( 0.995 ) ) {
  return nlip_newton_ls_factory< Function, GradFunction, HessianFunction, T >(
    f, df, fill_hessian, mu, max_iter, no_constraint_functor(), no_constraint_jac_functor(), no_constraint_functor(),
    no_constraint_jac_functor(), tol, eta, tau );
};


/**
 * This functor is a factory class to construct a Non-Linear
 * Interior-Point optimizer routine that uses a line-search approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update). Use make_nlip_quasi_newton_ls to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this. So far, all tests have failed, there must be something really wrong in this
 * implementation, but I cannot pin-point it.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \tparam HessianUpdater The functor type to update the Hessian approximation of the function to optimize.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see
 * trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 */
template < typename Function, typename GradFunction, typename T, typename EqFunction = no_constraint_functor,
           typename EqJacFunction = no_constraint_jac_functor, typename IneqFunction = no_constraint_functor,
           typename IneqJacFunction = no_constraint_jac_functor, typename HessianUpdater = hessian_update_bfgs >
struct nlip_quasi_newton_ls_factory {
  Function f;
  GradFunction df;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  HessianUpdater update_hessian;

  typedef nlip_quasi_newton_ls_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                        IneqJacFunction, HessianUpdater > self;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the sufficient decrease in order to accept a line-search step.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   * \param aUpdateHessian The functor object that can update the approximate Hessian matrix of the function to be
   * optimized.
   */
  nlip_quasi_newton_ls_factory( Function aF, GradFunction aDf, T aMu, unsigned int aMaxIter,
                                EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
                                IneqFunction aH = IneqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
                                T aTol = T( 1e-6 ), T aEta = T( 1e-4 ), T aTau = T( 0.995 ),
                                HessianUpdater aUpdateHessian = HessianUpdater() )
      : f( aF ), df( aDf ), mu( aMu ), max_iter( aMaxIter ), tol( aTol ), eta( aEta ), tau( aTau ), g( aG ),
        fill_g_jac( aFillGJac ), h( aH ), fill_h_jac( aFillHJac ), update_hessian( aUpdateHessian ){};
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template < typename Vector >
  void operator()( Vector& x ) const {
    detail::nl_intpoint_method_ls_impl( f, df, hessian_update_dual_quasi< HessianUpdater >( update_hessian ), g,
                                        fill_g_jac, h, fill_h_jac, x, mu, max_iter, tol, eta, tau );
  };

  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"),
   * this parameter is positive and should start with a rather large value (relative to the scale
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param( T aMu ) {
    mu = aMu;
    return *this;
  };
  /**
   * Sets the maximum number of iterations to perform.
   */
  self& set_max_iteration( unsigned int aMaxIter ) {
    max_iter = aMaxIter;
    return *this;
  };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance( T aTol ) {
    tol = aTol;
    return *this;
  };
  /**
   * Sets the tolerance on the sufficient decrease in order to accept a line-search step.
   */
  self& set_decrease_tolerance( T aEta ) {
    eta = aEta;
    return *this;
  };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality
   * constraint (barrier).
   */
  self& set_barrier_step_margin( T aTau ) {
    tau = aTau;
    return *this;
  };


  /**
   * This function remaps the factory to one which will use the given equality constraints for
   * the search domain. The equality constraints must be formulated at g(x) = 0.
   * \tparam NewEqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewEqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_g The function that computes the equality constraints.
   * \param new_fill_g_jac The function that computes the jacobian matrix of the equality constaints.
   */
  template < typename NewEqFunction, typename NewEqJacFunction >
  nlip_quasi_newton_ls_factory< Function, GradFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                                IneqJacFunction, HessianUpdater >
    set_eq_constraints( NewEqFunction new_g, NewEqJacFunction new_fill_g_jac ) const {
    return nlip_quasi_newton_ls_factory< Function, GradFunction, T, NewEqFunction, NewEqJacFunction, IneqFunction,
                                         IneqJacFunction, HessianUpdater >(
      f, df, mu, max_iter, new_g, new_fill_g_jac, h, fill_h_jac, tol, eta, tau, update_hessian );
  };

  /**
   * This function remaps the factory to one which will use the given inequality constraints for
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template < typename NewIneqFunction, typename NewIneqJacFunction >
  nlip_quasi_newton_ls_factory< Function, GradFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                                NewIneqJacFunction, HessianUpdater >
    set_ineq_constraints( NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac ) const {
    return nlip_quasi_newton_ls_factory< Function, GradFunction, T, EqFunction, EqJacFunction, NewIneqFunction,
                                         NewIneqJacFunction, HessianUpdater >(
      f, df, mu, max_iter, g, fill_g_jac, new_h, new_fill_h_jac, tol, eta, tau, update_hessian );
  };

  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewHessianUpdater A new functor type to update the Hessian approximation of the function to optimize.
   * \param new_update_hessian The new functor object that can update the approximate Hessian matrix of the function to
   * be optimized.
   */
  template < typename NewHessianUpdater >
  nlip_quasi_newton_ls_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                NewHessianUpdater >
    set_hessian_updater( NewHessianUpdater new_update_hessian ) const {
    return nlip_quasi_newton_ls_factory< Function, GradFunction, T, EqFunction, EqJacFunction, IneqFunction,
                                         IneqJacFunction, NewHessianUpdater >(
      f, df, mu, max_iter, g, fill_g_jac, h, fill_h_jac, tol, eta, tau, new_update_hessian );
  };
};

/**
 * This function template creates a factory object to construct a Non-Linear
 * Interior-Point optimizer routine that uses a line-search approach and
 * incorporates equality and inequality constraints via the sequential quadratic
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update).
 * This algorithm is roughly as described in Nocedal's Numerical Optimization
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this. So far, all tests have failed, there must be something really wrong in this
 * implementation, but I cannot pin-point it.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
 * positive and should start with a rather large value (relative to the scale of the function) and will be progressively
 * decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the sufficient decrease in order to accept a line-search step.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
 * (barrier).
 */
template < typename Function, typename GradFunction, typename T >
nlip_quasi_newton_ls_factory< Function, GradFunction, T >
  make_nlip_quasi_newton_ls( Function f, GradFunction df, T mu, unsigned int max_iter, T tol = T( 1e-6 ),
                             T eta = T( 1e-4 ), T tau = T( 0.995 ) ) {
  return nlip_quasi_newton_ls_factory< Function, GradFunction, T >(
    f, df, mu, max_iter, no_constraint_functor(), no_constraint_jac_functor(), no_constraint_functor(),
    no_constraint_jac_functor(), tol, eta, tau );
};
};
};


#endif
