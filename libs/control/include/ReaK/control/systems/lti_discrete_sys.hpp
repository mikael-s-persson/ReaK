/**
 * \file lti_discrete_sys.hpp
 *
 * This library provides a class template which can be used to create a simple discrete-time LTI
 * state-space system, as used in ReaK::ctrl.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_LTI_DISCRETE_SYS_HPP
#define REAK_LTI_DISCRETE_SYS_HPP

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

#include "discrete_linear_sss_concept.hpp"

namespace ReaK {

namespace ctrl {

/**
 * This class template can be used to create a simple discrete-time LTI
 * state-space system, as used in ReaK::ctrl. A discrete-time LTI state-space system is
 * basically described by four system matrices (A,B,C,D) which make the linear mapping
 * between the current state and input and the next state and current output.
 * \tparam T The value-type of the system matrices.
 */
template < typename T >
class lti_discrete_sys : public named_object {
private:
  T dt;
  mat< T, mat_structure::square > A;
  mat< T, mat_structure::rectangular > B;
  mat< T, mat_structure::rectangular > C;
  mat< T, mat_structure::rectangular > D;

public:
  typedef lti_discrete_sys< T > self;
  typedef T value_type;
  typedef std::size_t size_type;

  typedef vect_n< T > point_type;
  typedef vect_n< T > point_difference_type;
  typedef self topology;

  typedef T time_type;
  typedef T time_difference_type;

  typedef vect_n< T > input_type;
  typedef vect_n< T > output_type;

  typedef mat< T, mat_structure::square > matrixA_type;
  typedef mat< T, mat_structure::rectangular > matrixB_type;
  typedef mat< T, mat_structure::rectangular > matrixC_type;
  typedef mat< T, mat_structure::rectangular > matrixD_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, input_dimensions = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, output_dimensions = 0 );

  /**
   * Parametrized and default constructor.
   * \param aX_size The size of the state-vector.
   * \param aU_size The size of the input-vector.
   * \param aY_size The size of the output-vector.
   * \param aDt The time-step of the discrete-time system.
   */
  lti_discrete_sys( size_type aX_size = 0, size_type aU_size = 0, size_type aY_size = 0,
                    const time_difference_type& aDt = 1, const std::string& aName = "" )
      : A( aX_size ), B( aX_size, aU_size ), C( aY_size, aX_size ), D( aY_size, aU_size ), dt( aDt ) {
    setName( aName );
  };

  /**
   * Standard copy-constructor.
   */
  lti_discrete_sys( const self& rhs ) : A( rhs.A ), B( rhs.B ), C( rhs.C ), D( rhs.D ), dt( rhs.dt ) {
    setName( rhs.getName() );
  };

  /**
   * Parametrized constructor.
   * \param aA The discrete-time system matrix A.
   * \param aB The discrete-time system matrix B.
   * \param aC The discrete-time system matrix C.
   * \param aD The discrete-time system matrix D.
   * \param aDt The time-step of the discrete-time system.
   */
  template < typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD >
  lti_discrete_sys(
    const MatrixA& aA, const MatrixB& aB, const MatrixC& aC, const MatrixD& aD, const time_difference_type& aDt,
    const std::string& aName = "",
    typename boost::enable_if_c< is_readable_matrix< MatrixA >::value && is_readable_matrix< MatrixB >::value
                                 && is_readable_matrix< MatrixC >::value && is_readable_matrix< MatrixD >::value,
                                 void* >::type dummy = NULL )
      : A( aA ), B( aB ), C( aC ), D( aD ), dt( aDt ) {
    setName( aName );
  };


  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.A, rhs.A );
    swap( lhs.B, rhs.B );
    swap( lhs.C, rhs.C );
    swap( lhs.D, rhs.D );
    swap( lhs.dt, rhs.dt );
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) {
    swap( *this, rhs );
    return *this;
  };

  /**
   * Sets the discrete-time system matrix A of the discrete-time system.
   * \param aA The new discrete-time system matrix A for the system.
   */
  template < typename MatrixA >
  typename boost::enable_if_c< is_readable_matrix< MatrixA >::value, void >::type setA( const MatrixA& aA ) {
    A = aA;
  };

  /**
   * Sets the discrete-time system matrix B of the discrete-time system.
   * \param aB The new discrete-time system matrix B for the system.
   */
  template < typename MatrixB >
  typename boost::enable_if_c< is_readable_matrix< MatrixB >::value, void >::type setB( const MatrixB& aB ) {
    B = aB;
  };

  /**
   * Sets the discrete-time system matrix C of the discrete-time system.
   * \param aC The new discrete-time system matrix C for the system.
   */
  template < typename MatrixC >
  typename boost::enable_if_c< is_readable_matrix< MatrixC >::value, void >::type setC( const MatrixC& aC ) {
    C = aC;
  };

  /**
   * Sets the discrete-time system matrix D of the discrete-time system.
   * \param aD The new discrete-time system matrix D for the system.
   */
  template < typename MatrixD >
  typename boost::enable_if_c< is_readable_matrix< MatrixD >::value, void >::type setD( const MatrixD& aD ) {
    D = aD;
  };

  /**
   * Sets the time-step of the discrete-time system.
   * \param aDt The new time-step for the system.
   */
  void setDt( const time_difference_type& aDt ) { dt = aDt; };

  /**
   * Returns the discrete-time system matrix A of the system.
   * \return The discrete-time system matrix A of the system.
   */
  const matrixA_type& getA() const { return A; };
  /**
   * Returns the discrete-time system matrix B of the system.
   * \return The discrete-time system matrix B of the system.
   */
  const matrixB_type& getB() const { return B; };
  /**
   * Returns the discrete-time system matrix C of the system.
   * \return The discrete-time system matrix C of the system.
   */
  const matrixC_type& getC() const { return C; };
  /**
   * Returns the discrete-time system matrix D of the system.
   * \return The discrete-time system matrix D of the system.
   */
  const matrixD_type& getD() const { return D; };
  /**
   * Returns the time-step of the discrete-time system.
   */
  const time_difference_type& getDt() const { return dt; };

  /**
   * Returns the dimensions of the state vectors.
   * \return The dimensions of the state vectors.
   */
  size_type get_state_dimensions() const { return A.get_row_count(); };

  /**
   * Returns the dimensions of the input vectors.
   * \return The dimensions of the input vectors.
   */
  size_type get_input_dimensions() const { return B.get_col_count(); };

  /**
   * Returns the dimensions of the output vectors.
   * \return The dimensions of the output vectors.
   */
  size_type get_output_dimensions() const { return C.get_row_count(); };

  /**
   * Fills the given matrices with the discrete-time system's state transition matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   */
  template < typename MatrixA, typename MatrixB >
  typename boost::enable_if_c< is_writable_matrix< MatrixA >::value && is_writable_matrix< MatrixB >::value,
                               void >::type
    get_state_transition_blocks( MatrixA& aA, MatrixB& aB ) const {
    aA = A;
    aB = B;
  };

  /**
   * Fills the given matrices with the discrete-time system's state transition matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   */
  template < typename MatrixA, typename MatrixB, typename StateSpaceType >
  typename boost::enable_if_c< is_writable_matrix< MatrixA >::value && is_writable_matrix< MatrixB >::value,
                               void >::type
    get_state_transition_blocks( MatrixA& aA, MatrixB& aB, const StateSpaceType&, const time_type& = time_type(),
                                 const time_type& = time_type(), const point_type& = point_type(),
                                 const point_type& = point_type(), const input_type& = input_type(),
                                 const input_type& = input_type() ) const {
    aA = A;
    aB = B;
  };

  /**
   * Fills the given matrices with the discrete-time system's output function matrices.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template < typename MatrixC, typename MatrixD >
  typename boost::enable_if_c< is_writable_matrix< MatrixC >::value && is_writable_matrix< MatrixD >::value,
                               void >::type
    get_output_function_blocks( MatrixC& aC, MatrixD& aD ) const {
    aC = C;
    aD = D;
  };

  /**
   * Fills the given matrices with the discrete-time system's output function matrices.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template < typename MatrixC, typename MatrixD, typename StateSpaceType >
  typename boost::enable_if_c< is_writable_matrix< MatrixC >::value && is_writable_matrix< MatrixD >::value,
                               void >::type
    get_output_function_blocks( MatrixC& aC, MatrixD& aD, const StateSpaceType&, const time_type& = time_type(),
                                const point_type& = point_type(), const input_type& = input_type() ) const {
    aC = C;
    aD = D;
  };

  /**
   * Returns the time-step of the discrete-time system.
   */
  time_difference_type get_time_step() const { return dt; };

  /**
   * Returns next state given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The next state, at t + get_time_step().
   */
  template < typename StateSpaceType >
  point_type get_next_state( const StateSpaceType&, const point_type& p, const input_type& u,
                             const time_type& = 0 ) const {
    return A * p + B * u;
  };

  /**
   * Returns output of the system given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The current output.
   */
  template < typename StateSpaceType >
  output_type get_output( const StateSpaceType&, const point_type& p, const input_type& u,
                          const time_type& = 0 ) const {
    return C * p + D * u;
  };


  virtual void RK_CALL save( ReaK::serialization::oarchive& aA, unsigned int ) const {
    ReaK::named_object::save( aA, ReaK::named_object::getStaticObjectType()->TypeVersion() );
    aA& RK_SERIAL_SAVE_WITH_NAME( dt ) & RK_SERIAL_SAVE_WITH_NAME( A ) & RK_SERIAL_SAVE_WITH_NAME( B )
      & RK_SERIAL_SAVE_WITH_NAME( C ) & RK_SERIAL_SAVE_WITH_NAME( D );
  };
  virtual void RK_CALL load( ReaK::serialization::iarchive& aA, unsigned int ) {
    ReaK::named_object::load( aA, ReaK::named_object::getStaticObjectType()->TypeVersion() );
    aA& RK_SERIAL_LOAD_WITH_NAME( dt ) & RK_SERIAL_LOAD_WITH_NAME( A ) & RK_SERIAL_LOAD_WITH_NAME( B )
      & RK_SERIAL_LOAD_WITH_NAME( C ) & RK_SERIAL_LOAD_WITH_NAME( D );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2300004, 1, "lti_discrete_sys", named_object )
};
};
};

#endif
