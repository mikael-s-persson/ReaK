/**
 * \file MEAQR_topology.hpp
 *
 * This library provides classes
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_MEAQR_TOPOLOGY_HPP
#define REAK_MEAQR_TOPOLOGY_HPP


#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include <ReaK/topologies/spaces/direct_kinematics_topomap.hpp> // for write_joint_coordinates_impl
#include <ReaK/topologies/spaces/hyperbox_topology.hpp>
#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/steerable_space_concept.hpp>
#include <ReaK/topologies/spaces/proper_metric_concept.hpp>
#include <ReaK/topologies/spaces/tuple_distance_metrics.hpp>

#include <ReaK/topologies/interpolation/discrete_point_path.hpp>

#include <ReaK/control/systems/linear_ss_system_concept.hpp>

#include "IHAQR_topology.hpp"

#include <ReaK/control/integrators/runge_kutta4_integrator_sys.hpp>

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>
#include <ReaK/math/lin_alg/mat_operators.hpp>
#include <ReaK/math/lin_alg/mat_are_solver.hpp>
#include <ReaK/math/lin_alg/mat_exp_methods.hpp>

namespace ReaK {

namespace pp {

namespace detail {


template < typename StateSpace, typename StateSpaceSystem >
class MEAQR_ZIR_system {
public:
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;

  typedef arithmetic_tuple< mat< double, mat_structure::square >, mat< double, mat_structure::square >,
                            vect_n< double >, vect_n< double > > point_type;
  typedef point_type point_difference_type;
  typedef point_type point_derivative_type;

  typedef vector_topology< point_type > state_space_type;

  typedef double time_type;
  typedef double time_difference_type;

  typedef vect_n< double > input_type;
  typedef vect_n< double > output_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, input_dimensions = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, output_dimensions = 0 );

  typedef typename IHAQR_point_type< StateSpace, StateSpaceSystem >::linearization_payload lin_payload;
  typedef typename IHAQR_point_type< StateSpace, StateSpaceSystem >::IHAQR_payload IHAQR_payload;

public:
  lin_payload* lin_data;
  mat< double, mat_structure::square > BRBmatrix;
  vect_n< double > stat_drift;

public:
  std::size_t get_state_dimensions() const {
    return BRBmatrix.get_row_count() * ( 2 * BRBmatrix.get_row_count() + 2 );
  };

  std::size_t get_input_dimensions() const { return 0; };

  std::size_t get_output_dimensions() const { return 0; };

  MEAQR_ZIR_system( lin_payload* aLinData = nullptr,
                    const mat< double, mat_structure::diagonal >& aR = ( mat< double, mat_structure::diagonal >() ) )
      : lin_data( aLinData ) {
    if( !lin_data )
      return;

    BRBmatrix = lin_data->B * invert( aR ) * transpose_view( lin_data->B );
    stat_drift = to_vect< double >( lin_data->c );
  };


  point_derivative_type get_state_derivative( const state_space_type&, const point_type& p, input_type, double ) const {
    if( !lin_data )
      return point_derivative_type();

    const mat< double, mat_structure::square >& L = get< 0 >( p );
    const mat< double, mat_structure::square >& H = get< 1 >( p );
    const vect_n< double >& eta = get< 2 >( p );
    const vect_n< double >& zir = get< 3 >( p );

    // Compute derivative
    vect_n< double > T_rhs = BRBmatrix * eta - stat_drift;
    mat_vect_adaptor< vect_n< double > > T_rhs_m( T_rhs );
    ReaK::detail::backsub_Cholesky_impl( H, T_rhs_m );

    mat< double, mat_structure::square > Temp_L = lin_data->A * L;
    ReaK::detail::forwardsub_L_impl( L, Temp_L, 1e-6 );
    mat< double, mat_structure::square > F = BRBmatrix;
    ReaK::detail::forwardsub_L_impl( L, F, 1e-6 );
    F = transpose( F );
    ReaK::detail::forwardsub_L_impl( L, F, 1e-6 );
    F += Temp_L;
    F += transpose_view( Temp_L );
    for( std::size_t i = 0; i < F.get_row_count(); ++i )
      F( i, i ) *= 0.5;
    // L_dot = L * tril(F);
    ReaK::detail::inplace_lower_multiply_with_fill_impl( L, F );

    Temp_L = lin_data->A * H;
    ReaK::detail::forwardsub_L_impl( H, Temp_L, 1e-6 );
    mat< double, mat_structure::square > G = BRBmatrix;
    ReaK::detail::forwardsub_L_impl( H, G, 1e-6 );
    G = transpose( G );
    ReaK::detail::forwardsub_L_impl( H, G, 1e-6 );
    G -= Temp_L;
    G -= transpose_view( Temp_L );
    for( std::size_t i = 0; i < G.get_row_count(); ++i )
      G( i, i ) *= 0.5;
    // H_dot     = H * tril(G);
    ReaK::detail::inplace_lower_multiply_with_fill_impl( H, G );

    return point_derivative_type( std::move( F ),                              // L_dot
                                  std::move( G ),                              // H_dot
                                  transpose_view( lin_data->A ) * eta - T_rhs, // eta_dot
                                  lin_data->A * zir + stat_drift );            // zir_dot
  };

  output_type get_output( const state_space_type&, const point_type&, input_type, double ) const {
    return output_type();
  };
};
};


template < typename StateSpace, typename StateSpaceSystem >
class MEAQR_point_type : public IHAQR_point_type< StateSpace, StateSpaceSystem > {
public:
  typedef IHAQR_point_type< StateSpace, StateSpaceSystem > base_type;
  typedef MEAQR_point_type< StateSpace, StateSpaceSystem > self;

  typedef typename detail::MEAQR_ZIR_system< StateSpace, StateSpaceSystem >::point_type MEAQR_bundle_type;

  typedef typename base_type::matrixA_type matrixA_type;
  typedef typename base_type::matrixB_type matrixB_type;

  typedef typename base_type::state_type state_type;
  typedef typename base_type::state_difference_type state_difference_type;
  typedef typename base_type::state_derivative_type state_derivative_type;
  typedef typename base_type::system_input_type system_input_type;

  typedef typename base_type::linearization_payload linearization_payload;
  typedef typename base_type::IHAQR_payload IHAQR_payload;

  struct MEAQR_payload {
    std::vector< std::pair< double, MEAQR_bundle_type > > pts;
  };

  mutable shared_ptr< MEAQR_payload > MEAQR_data;

  MEAQR_point_type( const base_type& rhs ) : base_type( rhs ){};
  explicit MEAQR_point_type( const state_type& aX = state_type() ) : base_type( aX ){};

  MEAQR_point_type& operator()( const base_type& rhs ) {
    this->x = rhs.x;
    this->lin_data = rhs.lin_data;
    this->IHAQR_data = rhs.IHAQR_data;
    this->MEAQR_data = shared_ptr< MEAQR_payload >();
    return *this;
  };

  MEAQR_point_type( base_type&& rhs ) : base_type( std::move( rhs ) ){};
  explicit MEAQR_point_type( state_type&& aX ) : base_type( std::move( aX ) ){};

  MEAQR_point_type& operator()( base_type&& rhs ) {
    this->x = std::move( rhs.x );
    this->lin_data = std::move( rhs.lin_data );
    this->IHAQR_data = std::move( rhs.IHAQR_data );
    this->MEAQR_data = shared_ptr< MEAQR_payload >();
    return *this;
  };


  virtual ~MEAQR_point_type(){};

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    base_type::save( A, base_type::getStaticObjectType()->TypeVersion() );
  };
  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    base_type::load( A, base_type::getStaticObjectType()->TypeVersion() );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2400034, 1, "MEAQR_point_type", base_type )
};


class MEAQR_to_state_mapper : public named_object {
public:
  MEAQR_to_state_mapper() : named_object() { setName( "MEAQR_to_state_mapper" ); };

  template < typename SourceSpace, typename StateSpace, typename StateSpaceSystem >
  typename topology_traits< StateSpace >::point_type
    map_to_space( const MEAQR_point_type< StateSpace, StateSpaceSystem >& pt, const SourceSpace&,
                  const StateSpace& ) const {
    return pt.x;
  };

  template < typename DestSpace, typename StateSpace >
  typename topology_traits< DestSpace >::point_type
    map_to_space( const typename topology_traits< StateSpace >::point_type& pt, const StateSpace&,
                  const DestSpace& ) const {
    return typename topology_traits< DestSpace >::point_type( pt );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
  };
  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    named_object::load( A, named_object::getStaticObjectType()->TypeVersion() );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( MEAQR_to_state_mapper, 0xC2400038, 1, "MEAQR_to_state_mapper", named_object )
};


/**
 * This class implements a topology on a system controlled by an minimum-energy affine
 * quadratic regulator (MEAQR). This topology class takes an underlying state-space (topology),
 * a state-space system of that topology, and a sampler that is appropriate for the state-space.
 * One main feature of this topology is that the distance-metric is a reflection of the
 * cost to go from one state to another by a finite-horizon minimum-energy AQR controller,
 * as defined by the R matrix (quadratic input cost) and the ratio of the cost for the
 * input-bias (drift-compensation) and the feedback input (driven by the state-space error).
 * Another main feature is that movements between points are governed by the MEAQR controller
 * applied to the underlying state-space system (dynamics), and numerically integrated. This
 * topology fulfills the SteerableSpaceConcept which means that the record of a steer between
 * two points can be recorded as a path in the state-space.
 * \tparam StateSpace A topology type which represents the space in which the states of the system can exist, should
 * model the TopologyConcept.
 * \tparam StateSpaceSystem A state-space system type, should model the SSSystemConcept and the LinearSSSystemConcept
 * (linearizable system).
 * \tparam StateSpaceSampler A random sampler type for the given state-space, should model the RandomSamplerConcept.
 */
template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
class MEAQR_topology : public named_object {
public:
  typedef MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > self;

  typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixA_type matrixA_type;
  typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixB_type matrixB_type;
  typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixC_type matrixC_type;
  typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixD_type matrixD_type;

  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
  typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;

  typedef StateSpaceSampler state_sampler_type;

  typedef IHAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > IHAQR_space_type;
  typedef typename IHAQR_space_type::point_type IHAQR_point_type;
  typedef typename IHAQR_space_type::point_difference_type IHAQR_point_difference_type;

  typedef typename detail::MEAQR_ZIR_system< StateSpace, StateSpaceSystem >::point_type MEAQR_bundle_type;

  typedef MEAQR_point_type< StateSpace, StateSpaceSystem > point_type;
  typedef typename point_type::MEAQR_payload MEAQR_payload;

  typedef IHAQR_point_difference_type point_difference_type;

  typedef default_distance_metric distance_metric_type;
  typedef default_random_sampler random_sampler_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = 0 );

  typedef discrete_point_path< StateSpace > steer_record_type;

protected:
  shared_ptr< IHAQR_space_type > m_IHAQR_space; ///< Holds the IHAQR_topology for the same system. The IHAQR topology
  /// has many of the same parameters and necessary functions for the
  /// MEAQR, and is thus stored internally in the MEAQR instead of
  /// repeating all that code.
  double m_idle_to_cost_ratio; ///< The ratio of the idle-power cost to the feedback input cost, e.g., a value of 0.1
                               /// would mean that the input-bias cost (idle or hover cost) is penalized ten times less
                               /// than the feedback term.

  /**
   * This function fills the vector of MEAQR data points (P,M,eta,D) for a given point p and up to the maximum
   * time-horizon.
   * \param p The point for which the data points are required.
   */
  void compute_MEAQR_data( const point_type& p ) const {
    m_IHAQR_space->compute_linearization_data( p );
    if( p.MEAQR_data )
      return;
    p.MEAQR_data = shared_ptr< MEAQR_payload >( new MEAQR_payload() );

    detail::MEAQR_ZIR_system< StateSpace, StateSpaceSystem > MEAQR_sys( p.lin_data.get(), m_IHAQR_space->m_R );

    vector_topology< MEAQR_bundle_type > MEAQR_sys_space;

    double init_damping_value = 1e-5 * norm_inf( MEAQR_sys.BRBmatrix );
    std::size_t N = p.lin_data->A.get_row_count();
    MEAQR_bundle_type start_point(
      mat< double, mat_structure::square >( mat< double, mat_structure::scalar >( N, init_damping_value ) ),
      mat< double, mat_structure::square >( mat< double, mat_structure::scalar >( N, 100.0 * init_damping_value ) ),
      vect_n< double >( N, 0.0 ), vect_n< double >( N, 0.0 ) );

    ReaK::detail::decompose_Cholesky_impl( MEAQR_sys.BRBmatrix * m_IHAQR_space->m_time_step + get< 0 >( start_point ),
                                           get< 0 >( start_point ), 1e-8 );
    get< 1 >( start_point ) = get< 0 >( start_point );

    MEAQR_bundle_type current_point = start_point;

    typedef typename detail::MEAQR_ZIR_system< StateSpace, StateSpaceSystem >::input_type InputType;
    constant_trajectory< vector_topology< InputType > > input_traj
      = constant_trajectory< vector_topology< InputType > >( InputType() );

    // while not reached (fly-by) the goal yet:
    double current_time = 0.0;
    p.MEAQR_data->pts.push_back( std::make_pair( current_time, current_point ) );
    while( current_time < m_IHAQR_space->m_max_time_horizon ) {
      // integrate for one time-step.
      ReaK::ctrl::detail::runge_kutta4_integrate_impl(
        MEAQR_sys_space, MEAQR_sys, start_point, current_point, input_traj, current_time,
        current_time + m_IHAQR_space->m_time_step, m_IHAQR_space->m_time_step );
      start_point = current_point;
      current_time += m_IHAQR_space->m_time_step;
      p.MEAQR_data->pts.push_back( std::make_pair( current_time, current_point ) );
    };

    return;
  };

  /**
   * This virtual function evaluates if a given state-space point is within free-space (non-colliding).
   * \param a The state for which collision is checked.
   * \return True if the given state-space point is not colliding (in free-space).
   */
  virtual bool is_free_impl( const state_type& a ) const { return true; };


  /**
   * This function computes the cost-to-go and optimal time-horizon for a potential travel between two points.
   * \param a The starting point of the travel.
   * \param b The end point of the travel.
   * \param K The constant gain matrix at the destination (R^-1 * B^T)
   * \param u0 The input bias or drift conpensation at the destination point.
   * \param T_optimal Stores, as output, the optimal time-horizon for the travel.
   * \param J_optimal Stores, as output, the optimal cost-to-go for the travel.
   * \param it_optimal Stores, as output, an iterator to the MEAQR data point corresponding to where T_optimal and
   * J_optimal are located.
   */
  void find_optimal_cost_and_time(
    const point_type& a, const point_type& b, const mat< double, mat_structure::rectangular >& K,
    const system_input_type& u0, double& T_optimal, double& J_optimal,
    typename std::vector< std::pair< double, MEAQR_bundle_type > >::const_iterator& it_optimal ) const {
    if( !b.lin_data )
      m_IHAQR_space->compute_linearization_data( b );
    if( !b.MEAQR_data )
      compute_MEAQR_data( b );

    J_optimal = std::numeric_limits< double >::infinity();
    T_optimal = std::numeric_limits< double >::infinity();
    it_optimal = b.MEAQR_data->pts.end();

    double rho = m_idle_to_cost_ratio * ( u0 * ( m_IHAQR_space->m_R * u0 ) );

    typedef typename std::vector< std::pair< double, MEAQR_bundle_type > >::const_iterator Iter;

    mat< double, mat_structure::square > eAdt( b.lin_data->A.get_row_count() );
    ReaK::exp_PadeSAS( b.lin_data->A * m_IHAQR_space->m_time_step, eAdt, QR_linlsqsolver(), 1e-6 );
    mat< double, mat_structure::square > eAt = eAdt;

    vect_n< double > a_rel = to_vect< double >( m_IHAQR_space->m_space.difference( a.x, b.x ) );

    for( Iter it = b.MEAQR_data->pts.begin(); it != b.MEAQR_data->pts.end(); ++it ) {

      vect_n< double > s = eAt * a_rel + get< 3 >( it->second ); // ZIR(T)  from Nir's paper.

      // s = invert(L) * s
      mat_vect_adaptor< vect_n< double > > s_m( s );
      ReaK::detail::forwardsub_L_impl( get< 0 >( it->second ), s_m, 1e-6 );

      double J = rho * it->first + 0.5 * ( s * s );

      ReaK::detail::backsub_R_impl( transpose_view( get< 0 >( it->second ) ), s_m, 1e-6 );
      system_input_type u_T = u0 - ( K * transpose_view( eAt ) * s ) - K * get< 2 >( it->second );

      if( ( J < J_optimal ) && ( m_IHAQR_space->m_input_space.is_in_bounds( u_T ) ) ) {
        J_optimal = J;
        T_optimal = it->first;
        it_optimal = it;
      };
      if( J_optimal < rho * it->first )
        break;

      eAt = eAdt * eAt;
    };

    if( it_optimal == b.MEAQR_data->pts.end() - 1 ) {
      T_optimal += ( J_optimal / rho - T_optimal ) * 0.5;
      J_optimal = 0.25 * rho * m_IHAQR_space->m_max_time_horizon + 0.75 * J_optimal;
    };
  };


  /**
   * This function attempts to travel (steer) for some time interval.
   * \param H The current H matrix (lower-triangular part of M^-1).
   * \param K The current feedback gain matrix (R^-1 B^T).
   * \param eta The current drift term (affine term of the MEAQR).
   * \param u0 The current input-bias (or drift compensation) for the destination point.
   * \param u_prev The previous input to the system, used for bandwidth limitations on the input.
   * \param x_current The current state-space point (will be modified to contain the resulting state-space point).
   * \param x_goal The ultimate state-space point to be reached (the steering will stop if it is reached at any time).
   * \param current_time The current time (will be modified to contain the resulting time).
   * \param time_limit The time limit after which to stop this steering segment.
   * \param with_collision_check A flag to tell if collision should be checked or not (through the virtual
   * is_free_impl() function).
   * \param st_rec A pointer to a steer-record object. If nullptr, then the steering doesn't get recorded.
   * \return False if a collision occurred.
   */
  bool steer_with_constant_control( const mat< double, mat_structure::square >& H,
                                    const mat< double, mat_structure::rectangular >& K, const vect_n< double >& eta,
                                    const system_input_type& u0, system_input_type& u_prev, state_type& x_current,
                                    const state_type& x_goal, double& current_time, double time_limit,
                                    bool with_collision_check, steer_record_type* st_rec = nullptr ) const {
    bool was_collision_free = true;
    state_type x_next = x_current;
    // while not reached (fly-by) the goal yet:
    while( ( current_time < time_limit )
           && ( m_IHAQR_space->m_space.distance( x_current, x_goal ) > m_IHAQR_space->m_goal_proximity_threshold ) ) {
      // compute the current MEAQR input
      vect_n< double > HHx = to_vect< double >( m_IHAQR_space->m_space.difference( x_current, x_goal ) );
      mat_vect_adaptor< vect_n< double > > HHx_m( HHx );
      ReaK::detail::backsub_Cholesky_impl( H, HHx_m );


      system_input_type u_current;
      if( current_time < m_IHAQR_space->m_time_step ) {
        u_current = u0 - K * eta - K * HHx; // don't saturate the first time step.
      } else {
        u_current = m_IHAQR_space->get_bounded_input( u_prev, u0 - K * eta, -K * HHx );
      };
      //         system_input_type u_current = u0 - K * eta - K * HHx;
      //         std::cout << " current u = " << u_current << std::endl;

      //         accum_steer_cost += to_vect<double>(u_current) * (m_IHAQR_space->m_R * to_vect<double>(u_current)) *
      //         m_IHAQR_space->m_time_step;

      constant_trajectory< vector_topology< system_input_type > > input_traj( u_current );

      // integrate for one time-step.
      ReaK::ctrl::detail::runge_kutta4_integrate_impl(
        m_IHAQR_space->m_space, *( m_IHAQR_space->m_system ), x_current, x_next, input_traj, current_time,
        current_time + m_IHAQR_space->m_time_step, m_IHAQR_space->m_time_step * 1e-1 );
      //         std::cout << " current x = " << x_current << std::endl;

      if( ( !with_collision_check ) || is_free_impl( x_next ) ) {
        x_current = x_next;
        current_time += m_IHAQR_space->m_time_step;
        u_prev = u_current;
        if( st_rec )
          st_rec->push_back( x_current );
      } else {
        was_collision_free = false;
        break;
      };
    };
    return was_collision_free;
  };


  /**
   * This function attempts to travel (steer) between two points for a given fraction between them.
   * \param a The starting point of the travel.
   * \param fraction The fraction of the complete travel that should be done.
   * \param b The end point of the travel.
   * \param with_collision_check A flag to tell if collision should be checked or not (through the virtual
   * is_free_impl() function).
   * \param st_rec A pointer to a steer-record object. If nullptr, then the steering doesn't get recorded.
   * \return The resulting point after the steering.
   */
  point_type move_position_toward_impl( const point_type& a, double fraction, const point_type& b,
                                        bool with_collision_check, steer_record_type* st_rec = nullptr ) const {
    if( !a.lin_data )
      m_IHAQR_space->compute_linearization_data( a );
    if( !b.lin_data )
      m_IHAQR_space->compute_linearization_data( b );
    if( !b.MEAQR_data )
      compute_MEAQR_data( b );

    typedef typename std::vector< std::pair< double, MEAQR_bundle_type > >::const_iterator Iter;

    mat< double, mat_structure::rectangular > K = invert( m_IHAQR_space->m_R ) * transpose_view( b.lin_data->B );

    // first, find the T-optimal:
    double min_J;
    double T_optimal;
    Iter min_it;
    find_optimal_cost_and_time( a, b, K, b.lin_data->u, T_optimal, min_J, min_it );
    if( min_it == b.MEAQR_data->pts.end() )
      return a; // cost is infinite!

    // scale back the end-time to meet the fractional travel requirement:
    double T_goal = fraction * ( T_optimal + 2.0 * m_IHAQR_space->m_time_step );

    //       std::cout << " t_optimal = " << T_optimal << std::endl;

    // compute basic controller components: gain matrix and previous actuation.
    system_input_type u_prev = a.lin_data->u;
    m_IHAQR_space->m_input_space.bring_point_in_bounds( u_prev );

    double current_time = 0.0;
    state_type x_current = a.x;
    if( st_rec )
      st_rec->push_back( x_current );

    // First, steer up to the time-horizon if the optimal time was beyond it.
    if( min_it == b.MEAQR_data->pts.end() - 1 ) {
      mat< double, mat_structure::square > H = get< 1 >( min_it->second );
      vect_n< double > eta = get< 2 >( min_it->second );

      double time_limit = T_optimal - m_IHAQR_space->m_max_time_horizon;
      if( time_limit > T_goal )
        time_limit = T_goal;
      if( !steer_with_constant_control( get< 1 >( min_it->second ), K,                     // H, K
                                        get< 2 >( min_it->second ), b.lin_data->u, u_prev, // eta, u0, u_previous
                                        x_current, b.x, current_time, time_limit, with_collision_check, st_rec ) )
        return point_type( x_current ); // resulted in a collision, output last non-colliding point.
    };

    // Then, iterate through the MEAQR sequence points.
    while( ( current_time < T_goal ) && ( min_it != b.MEAQR_data->pts.begin() ) ) {
      Iter prev_it = min_it;
      --min_it;

      // integrate for some time steps.
      double time_limit = T_optimal - min_it->first;
      if( time_limit > T_goal )
        time_limit = T_goal;
      if( !steer_with_constant_control( get< 1 >( prev_it->second ), K,                     // H, K
                                        get< 2 >( prev_it->second ), b.lin_data->u, u_prev, // eta, u0, u_previous
                                        x_current, b.x, current_time, time_limit, with_collision_check, st_rec ) )
        return point_type( x_current ); // resulted in a collision, output last non-colliding point.
    };

    // Finally, apply a bit more control, in case the goal state wasn't quite reached yet (mostly due to discrete-time
    // effects).
    steer_with_constant_control( get< 1 >( min_it->second ), K,                     // H, K
                                 get< 2 >( min_it->second ), b.lin_data->u, u_prev, // eta, u0, u_previous
                                 x_current, b.x, current_time, T_goal, with_collision_check, st_rec );

    return point_type( x_current ); // output last non-colliding point, whether collision occurred or not.
  };

  /**
   * This function samples a point from the free-space (if collision check is enabled).
   * \param with_collision_check A flag to tell if collision should be checked or not (through the virtual
   * is_free_impl() function).
   * \return The resulting point of the random sampling.
   */
  point_type random_point_impl( bool with_collision_check ) const {
    state_type result_pt = m_IHAQR_space->m_get_sample( m_IHAQR_space->m_space );
    while( with_collision_check && !is_free_impl( result_pt ) )
      result_pt = m_IHAQR_space->m_get_sample( m_IHAQR_space->m_space );
    point_type result( result_pt );
    m_IHAQR_space->compute_linearization_data( result );
    compute_MEAQR_data( result );
    return result;
  };

public:
  /**
   * This function returns a reference to the underlying state-space used by this topology.
   * \return A reference to the underlying state-space used by this topology.
   */
  StateSpace& get_state_space() { return m_IHAQR_space->get_state_space(); };
  /**
   * This function returns a const-reference to the underlying state-space used by this topology.
   * \return A const-reference to the underlying state-space used by this topology.
   */
  const StateSpace& get_state_space() const { return m_IHAQR_space->get_state_space(); };

  /**
   * This function returns a reference to the underlying IHAQR topology used by this topology.
   * \return A reference to the underlying IHAQR topology used by this topology.
   */
  IHAQR_space_type& get_IHAQR_space() { return *m_IHAQR_space; };
  /**
   * This function returns a const-reference to the underlying IHAQR topology used by this topology.
   * \return A const-reference to the underlying IHAQR topology used by this topology.
   */
  const IHAQR_space_type& get_IHAQR_space() const { return *m_IHAQR_space; };


  /**
   * This function returns a maximum allowable time-horizon.
   * \return A maximum allowable time-horizon.
   */
  double get_max_time_horizon() const { return m_IHAQR_space->m_max_time_horizon; };

  /**
   * This function returns the idle-power (hover-power) cost for a given point.
   * \return The idle-power (hover-power) cost for a given point.
   */
  double get_idle_power_cost( const point_type& b ) const {
    if( !b.lin_data )
      m_IHAQR_space->compute_linearization_data( b );
    return m_idle_to_cost_ratio * ( b.lin_data->u * ( m_IHAQR_space->m_R * b.lin_data->u ) );
  };


  /**
   * Default constructor.
   * \param aName The name of this topology / object.
   * \param aIHAQRSpace A pointer to a IHAQR topology that can be used by this MEAQR topology.
   * \param aIdleToCostRatio The ratio of the idle-power cost to the feedback input cost, e.g.,
   *                         a value of 0.1 would mean that the input-bias cost (idle or hover cost)
   *                         is penalized ten times less than the feedback term.
   */
  MEAQR_topology( const std::string& aName = "MEAQR_topology",
                  const shared_ptr< IHAQR_space_type >& aIHAQRSpace = shared_ptr< IHAQR_space_type >(),
                  double aIdleToCostRatio = 0.01 )
      : named_object(), m_IHAQR_space( aIHAQRSpace ), m_idle_to_cost_ratio( aIdleToCostRatio ) {
    setName( aName );
  };

  virtual ~MEAQR_topology(){};


  /*************************************************************************
  *                             MetricSpaceConcept
  * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance( const point_type& a, const point_type& b ) const {
    if( !b.lin_data )
      m_IHAQR_space->compute_linearization_data( b );
    if( !b.MEAQR_data )
      compute_MEAQR_data( b );

    typedef typename std::vector< std::pair< double, MEAQR_bundle_type > >::const_iterator Iter;

    mat< double, mat_structure::rectangular > K = invert( m_IHAQR_space->m_R ) * transpose_view( b.lin_data->B );

    double min_J;
    double T_optimal;
    Iter min_it;
    find_optimal_cost_and_time( a, b, K, b.lin_data->u, T_optimal, min_J, min_it );

    return min_J;
  };

  /**
   * Returns the norm of the difference between two points.
   */
  double norm( const point_difference_type& delta ) const {
    vect_n< double > dx = to_vect< double >( delta.dx );
    return dx * dx;
  };

  /**
   * Returns the proper distance between two points.
   */
  double proper_distance( const point_type& a, const point_type& b ) const {
    double d1 = distance( a, b );
    double d2 = distance( b, a );
    if( d1 < d2 )
      return d1;
    return d2;
  };

  /**
   * Returns the proper norm of the difference between two points.
   */
  double proper_norm( const point_difference_type& delta ) const { return norm( delta ); };

  /*************************************************************************
   *                         for PointDistributionConcept
   * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const { return random_point_impl( false ); };

  /*************************************************************************
   *                             TopologyConcept
   * **********************************************************************/

  /**
   * Returns the difference between two points (analogous to a - b, but implemented in SO(3) Lie algebra).
   */
  point_difference_type difference( const point_type& a, const point_type& b ) const {
    return point_difference_type( m_IHAQR_space->m_space.difference( b.x, a.x ) );
  };

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type adjust( const point_type& a, const point_difference_type& delta ) const {
    return point_type( m_IHAQR_space->m_space.adjust( a.x, delta.dx ) );
  };

  /**
   * Returns the origin of the space (the lower-limit).
   */
  point_type origin() const { return point_type( m_IHAQR_space->m_space.origin() ); };

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds( const point_type& a ) const { return m_IHAQR_space->m_space.is_in_bounds( a.x ); };

  // NOTE: don't know if I can get rid of this. (only seems useful in bounded interpolators (and samplers)).
  point_difference_type get_diff_to_boundary( const point_type& ) const { return point_difference_type(); };

  /*************************************************************************
  *                             LieGroupConcept
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward( const point_type& a, double fraction, const point_type& b ) const {
    return move_position_toward_impl( a, fraction, b, false );
  };

  /*************************************************************************
  *                             SteerableSpaceConcept
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  std::pair< point_type, steer_record_type > steer_position_toward( const point_type& a, double fraction,
                                                                    const point_type& b ) const {
    std::pair< point_type, steer_record_type > result(
      point_type(), steer_record_type( shared_ptr< StateSpace >( m_IHAQR_space, &( m_IHAQR_space->m_space ) ) ) );
    result.first = move_position_toward_impl( a, fraction, b, false, &result.second );
    return result;
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    ReaK::named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( m_IHAQR_space ) & RK_SERIAL_SAVE_WITH_NAME( m_idle_to_cost_ratio );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    ReaK::named_object::load( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( m_IHAQR_space ) & RK_SERIAL_LOAD_WITH_NAME( m_idle_to_cost_ratio );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2400035, 1, "MEAQR_topology", named_object )
};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_metric_space< MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > > : boost::mpl::true_ {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_point_distribution< MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > >
  : is_point_distribution< StateSpace > {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_metric_symmetric< MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > > : boost::mpl::false_ {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_steerable_space< MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > > : boost::mpl::true_ {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct get_proper_metric< MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > > {
  typedef default_proper_metric type;
};


/**
 * This class is an implementation of the MEAQR_topology with collision detection.
 */
template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
class MEAQR_topology_with_CD : public MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > {
public:
  typedef MEAQR_topology_with_CD< StateSpace, StateSpaceSystem, StateSpaceSampler > self;
  typedef MEAQR_topology< StateSpace, StateSpaceSystem, StateSpaceSampler > base_type;

  typedef base_type super_space_type;

  typedef typename base_type::state_type state_type;
  typedef typename base_type::state_difference_type state_difference_type;
  typedef typename base_type::state_derivative_type state_derivative_type;
  typedef typename base_type::system_input_type system_input_type;

  typedef typename base_type::point_type point_type;
  typedef typename base_type::point_difference_type point_difference_type;
  typedef default_distance_metric distance_metric_type;
  typedef default_random_sampler random_sampler_type;

  typedef typename base_type::IHAQR_space_type IHAQR_space_type;

  typedef typename base_type::steer_record_type steer_record_type;

protected:
  shared_ptr< kte::direct_kinematics_model > m_model;

  virtual bool is_free_impl( const state_type& a ) const {

    if( !this->get_state_space().is_in_bounds( a ) )
      return false;

    detail::write_joint_coordinates_impl( a, this->m_IHAQR_space->get_state_space(), m_model );
    // update the kinematics model with the given joint states.
    m_model->doDirectMotion();

    for( std::vector< shared_ptr< geom::proxy_query_pair_2D > >::const_iterator it = m_proxy_env_2D.begin();
         it != m_proxy_env_2D.end(); ++it ) {
      geom::proximity_record_2D tmp = ( *it )->findMinimumDistance();
      if( tmp.mDistance < 0.0 )
        return false;
    };
    for( std::vector< shared_ptr< geom::proxy_query_pair_3D > >::const_iterator it = m_proxy_env_3D.begin();
         it != m_proxy_env_3D.end(); ++it ) {
      geom::proximity_record_3D tmp = ( *it )->findMinimumDistance();
      if( tmp.mDistance < 0.0 )
        return false;
    };

    return true;
  };

public:
  std::vector< shared_ptr< geom::proxy_query_pair_2D > > m_proxy_env_2D;
  std::vector< shared_ptr< geom::proxy_query_pair_3D > > m_proxy_env_3D;

  /**
   * Default constructor.
   * \param aName The name of this topology / object.
   * \param aIHAQRSpace A pointer to a IHAQR topology that can be used by this MEAQR topology.
   * \param aIdleToCostRatio The ratio of the idle-power cost to the feedback input cost, e.g.,
   *                         a value of 0.1 would mean that the input-bias cost (idle or hover cost)
   *                         is penalized ten times less than the feedback term.
   * \param aModel The direct-kinematics model for the system (this is for the application of
   *               state-space points to synchronize the collision-detection code which rely on KTE-based systems).
   */
  MEAQR_topology_with_CD( const std::string& aName = "MEAQR_topology_with_CD",
                          const shared_ptr< IHAQR_space_type >& aIHAQRSpace = shared_ptr< IHAQR_space_type >(),
                          double aIdleToCostRatio = 0.01, const shared_ptr< kte::direct_kinematics_model >& aModel
                                                          = shared_ptr< kte::direct_kinematics_model >() )
      : base_type( aName, aIHAQRSpace, aIdleToCostRatio ), m_model( aModel ), m_proxy_env_2D(), m_proxy_env_3D(){};

  /**
   * Default constructor.
   * \param aBaseSpace A MEAQR_topology object to copy.
   * \param aModel The direct-kinematics model for the system (this is for the application of
   *               state-space points to synchronize the collision-detection code which rely on KTE-based systems).
   */
  MEAQR_topology_with_CD( const base_type& aBaseSpace, const shared_ptr< kte::direct_kinematics_model >& aModel )
      : base_type( aBaseSpace ), m_model( aModel ), m_proxy_env_2D(), m_proxy_env_3D(){};

  virtual ~MEAQR_topology_with_CD(){};

  super_space_type& get_super_space() { return *this; };

  const super_space_type& get_super_space() const { return *this; };


  /*************************************************************************
   *                         for PointDistributionConcept
   * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance( const point_type& a, const point_type& b ) const {
    point_type result = this->move_position_toward_impl( a, 1.0, b, true );
    if( get( distance_metric, this->get_state_space() )( a.x, b.x, this->get_state_space() ) * 0.05
        > get( distance_metric, this->get_state_space() )( result.x, b.x, this->get_state_space() ) ) {
      return base_type::distance( a, b );
    } else {
      return std::numeric_limits< double >::infinity();
    };
  };

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const { return this->random_point_impl( true ); };


  bool is_free( const point_type& a ) const { return this->is_free_impl( a.x ); };

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward( const point_type& a, double fraction, const point_type& b ) const {
    return this->move_position_toward_impl( a, fraction, b, true );
  };

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  std::pair< point_type, steer_record_type > steer_position_toward( const point_type& a, double fraction,
                                                                    const point_type& b ) const {
    std::pair< point_type, steer_record_type > result(
      point_type(), steer_record_type(
                      shared_ptr< StateSpace >( this->m_IHAQR_space, &( this->m_IHAQR_space->get_state_space() ) ) ) );
    result.first = this->move_position_toward_impl( a, fraction, b, true, &result.second );
    return result;
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    base_type::save( A, base_type::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( m_model ) & RK_SERIAL_SAVE_WITH_NAME( m_proxy_env_2D )
      & RK_SERIAL_SAVE_WITH_NAME( m_proxy_env_3D );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    base_type::load( A, base_type::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( m_model ) & RK_SERIAL_LOAD_WITH_NAME( m_proxy_env_2D )
      & RK_SERIAL_LOAD_WITH_NAME( m_proxy_env_3D );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2400036, 1, "MEAQR_topology_with_CD", base_type )
};


template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_metric_space< MEAQR_topology_with_CD< StateSpace, StateSpaceSystem, StateSpaceSampler > >
  : boost::mpl::true_ {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_point_distribution< MEAQR_topology_with_CD< StateSpace, StateSpaceSystem, StateSpaceSampler > >
  : is_point_distribution< StateSpace > {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_metric_symmetric< MEAQR_topology_with_CD< StateSpace, StateSpaceSystem, StateSpaceSampler > >
  : boost::mpl::false_ {};

template < typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler >
struct is_steerable_space< MEAQR_topology_with_CD< StateSpace, StateSpaceSystem, StateSpaceSampler > >
  : boost::mpl::true_ {};
};
};


#endif
