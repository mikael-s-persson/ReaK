/**
 * \file rotations_3D.hpp
 *
 * This library declares all geometric 3D rotation classes for fixed (2,3) and variable dimensions.
 *
 * Note: All matrix memory is organized, by default, such that columns are concatenated. This
 *       was found to be a more efficient representation since columns often have
 *       more significances than rows (representing basis vectors for example).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_ROTATIONS_3D_HPP
#define REAK_ROTATIONS_3D_HPP

#include <ReaK/math/lin_alg/mat_concepts.hpp>
#include <ReaK/math/lin_alg/mat_alg_square.hpp>
#include <ReaK/math/lin_alg/mat_alg_symmetric.hpp>
#include <ReaK/math/lin_alg/mat_alg_skew_symmetric.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

#include <cassert>

namespace ReaK {


// Forward declaration
template < class T >
class quaternion;

template < class T >
class unit_quat;

template < class T >
class euler_angles_TB;

template < class T >
class axis_angle;

template < class T >
class trans_mat_3D;


/**
 * This class is a rotation matrix 3 by 3.
 * \test All tests for this class have been passed!
 */
template < typename T >
class rot_mat_3D {
public:
  typedef rot_mat_3D< T > self;
  typedef void allocator_type;

  typedef T value_type;
  typedef void container_type;

  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  typedef void col_iterator;
  typedef void const_col_iterator;
  typedef void row_iterator;
  typedef void const_row_iterator;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 3 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 3 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_alignment::column_major );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::orthogonal );

private:
  value_type q[9];

  /**
   * Constructor from the components of the rotation matrix.
   * \test PASSED
   */
  explicit rot_mat_3D( const_reference a11, const_reference a12, const_reference a13, const_reference a21,
                       const_reference a22, const_reference a23, const_reference a31, const_reference a32,
                       const_reference a33 ) BOOST_NOEXCEPT {
    q[0] = a11;
    q[1] = a21;
    q[2] = a31;
    q[3] = a12;
    q[4] = a22;
    q[5] = a32;
    q[6] = a13;
    q[7] = a23;
    q[8] = a33;
  };

public:
  friend class quaternion< value_type >;
  friend class euler_angles_TB< value_type >;
  friend class axis_angle< value_type >;
  friend class trans_mat_3D< value_type >;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default constructor, sets the matrix to identity (no rotation).
   * \test PASSED
   */
  rot_mat_3D() BOOST_NOEXCEPT {
    q[0] = 1.0;
    q[1] = 0.0;
    q[2] = 0.0;
    q[3] = 0.0;
    q[4] = 1.0;
    q[5] = 0.0;
    q[6] = 0.0;
    q[7] = 0.0;
    q[8] = 1.0;
  };

  /**
   * Constructor from an array of components.
   * \test PASSED
   */
  rot_mat_3D( const_pointer M ) BOOST_NOEXCEPT {
    vect< value_type, 3 > v1 = unit( vect< value_type, 3 >( M ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    vect< value_type, 3 > v2( &M[3] );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[3] = v2[0];
    q[4] = v2[1];
    q[5] = v2[2];
    v2 = v1 % v2;
    q[6] = v2[0];
    q[7] = v2[1];
    q[8] = v2[2];
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  rot_mat_3D( const self& R ) BOOST_NOEXCEPT { std::copy( R.q, R.q + 9, q ); };
#else
  rot_mat_3D( const self& ) BOOST_NOEXCEPT = default;
#endif

  template < typename Matrix >
  explicit rot_mat_3D(
    const Matrix& M,
    typename boost::enable_if_c< is_readable_matrix< Matrix >::value && !boost::is_same< Matrix, self >::value,
                                 void* >::type dummy = nullptr ) {
    if( ( M.get_col_count() != 3 ) || ( M.get_row_count() != 3 ) )
      throw std::range_error( "Right-hand-side of assignment to a 3D rotation matrix is not of dimension 3x3!" );
    vect< value_type, 3 > v1( M( 0, 0 ), M( 1, 0 ), M( 2, 0 ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    vect< value_type, 3 > v2( M( 0, 1 ), M( 1, 1 ), M( 2, 1 ) );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[3] = v2[0];
    q[4] = v2[1];
    q[5] = v2[2];
    v2 = v1 % v2;
    q[6] = v2[0];
    q[7] = v2[1];
    q[8] = v2[2];
  };

  // Copy constructor. Default is good. \test PASSED

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides a copy of the rotation matrix as an ordinary 3x3 matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::square > getMat() const {
    return mat< value_type, mat_structure::square >( q[0], q[3], q[6], q[1], q[4], q[7], q[2], q[5], q[8] );
  };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  value_type operator[]( size_type i ) const BOOST_NOEXCEPT {
    assert( i < 9 );
    return q[i];
  };

  /**
   * Array double-indexing operator, ith row and jth column, accessor for read only.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const BOOST_NOEXCEPT {
    assert( ( i < 3 ) && ( j < 3 ) );
    return q[j * 3 + i];
  };

  size_type get_row_count() const BOOST_NOEXCEPT { return 3; };
  size_type get_col_count() const BOOST_NOEXCEPT { return 3; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  self& operator=( const self& M ) BOOST_NOEXCEPT {
    std::copy( M.q, M.q + 9, q );
    return *this;
  };
#else
  self& operator=( const self& M ) BOOST_NOEXCEPT = default;
#endif

  template < typename Matrix >
  typename boost::enable_if_c< is_readable_matrix< Matrix >::value && !boost::is_same< Matrix, self >::value,
                               self& >::type
    operator=( const Matrix& M ) {
    if( ( M.get_col_count() != 3 ) || ( M.get_row_count() != 3 ) )
      throw std::range_error( "Right-hand-side of assignment to a 3D rotation matrix is not of dimension 3x3!" );
    vect< value_type, 3 > v1( M( 0, 0 ), M( 1, 0 ), M( 2, 0 ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    vect< value_type, 3 > v2( M( 0, 1 ), M( 1, 1 ), M( 2, 1 ) );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[3] = v2[0];
    q[4] = v2[1];
    q[5] = v2[2];
    v2 = v1 % v2;
    q[6] = v2[0];
    q[7] = v2[1];
    q[8] = v2[2];
  };

  /**
   * Assignment operator from a quaternion representation.
   */
  self& operator=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT { return ( *this = Q.getRotMat() ); };

  /**
   * Assignment operator from a euler angles TB representation.
   */
  self& operator=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { return ( *this = E.getRotMat() ); };

  /**
   * Assignment operator from an axis / angle representation.
   */
  self& operator=( const axis_angle< value_type >& A ) BOOST_NOEXCEPT { return ( *this = A.getRotMat() ); };

  /**
   * Multiplication by a rotation matrix and store.
   * \test PASSED
   */
  self& operator*=( const self& M ) BOOST_NOEXCEPT {
    *this = self( q[0] * M.q[0] + q[3] * M.q[1] + q[6] * M.q[2], q[0] * M.q[3] + q[3] * M.q[4] + q[6] * M.q[5],
                  q[0] * M.q[6] + q[3] * M.q[7] + q[6] * M.q[8], q[1] * M.q[0] + q[4] * M.q[1] + q[7] * M.q[2],
                  q[1] * M.q[3] + q[4] * M.q[4] + q[7] * M.q[5], q[1] * M.q[6] + q[4] * M.q[7] + q[7] * M.q[8],
                  q[2] * M.q[0] + q[5] * M.q[1] + q[8] * M.q[2], q[2] * M.q[3] + q[5] * M.q[4] + q[8] * M.q[5],
                  q[2] * M.q[6] + q[5] * M.q[7] + q[8] * M.q[8] );
    return *this;
  };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiplication by a rotation matrix.
   * \test PASSED
   */
  friend self operator*(const self& M1, const self& M2)BOOST_NOEXCEPT {
    return self( M1.q[0] * M2.q[0] + M1.q[3] * M2.q[1] + M1.q[6] * M2.q[2],
                 M1.q[0] * M2.q[3] + M1.q[3] * M2.q[4] + M1.q[6] * M2.q[5],
                 M1.q[0] * M2.q[6] + M1.q[3] * M2.q[7] + M1.q[6] * M2.q[8],
                 M1.q[1] * M2.q[0] + M1.q[4] * M2.q[1] + M1.q[7] * M2.q[2],
                 M1.q[1] * M2.q[3] + M1.q[4] * M2.q[4] + M1.q[7] * M2.q[5],
                 M1.q[1] * M2.q[6] + M1.q[4] * M2.q[7] + M1.q[7] * M2.q[8],
                 M1.q[2] * M2.q[0] + M1.q[5] * M2.q[1] + M1.q[8] * M2.q[2],
                 M1.q[2] * M2.q[3] + M1.q[5] * M2.q[4] + M1.q[8] * M2.q[5],
                 M1.q[2] * M2.q[6] + M1.q[5] * M2.q[7] + M1.q[8] * M2.q[8] );
  };

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::enable_if_c< is_fully_writable_matrix< Matrix >::value, Matrix >::type
    operator*( const self& M1, const Matrix& M2 ) {
    if( M2.get_row_count() != 3 )
      throw std::range_error( "Matrix M's row count is not 3, 3D rotation impossible!" );
    Matrix result( M2 );
    for( size_type i = 0; i < 3; ++i )
      for( size_type jj = 0; jj < result.get_col_count(); ++jj ) {
        result( i, jj ) = 0;
        for( size_type j = 0; j < 3; ++j )
          result( i, jj ) += M1.q[j * 3 + i] * M2( j, jj );
      };
    return result;
  };

  template < typename Matrix >
  friend typename boost::enable_if_c< is_fully_writable_matrix< Matrix >::value, Matrix >::type
    operator*( const Matrix& M1, const self& M2 ) {
    if( M1.get_col_count() != 3 )
      throw std::range_error( "Matrix M1's column count is not 3, 3D rotation impossible!" );
    Matrix result( M1 );
    for( size_type i = 0; i < result.get_row_count(); ++i )
      for( size_type jj = 0; jj < 3; ++jj ) {
        result( i, jj ) = 0;
        for( size_type j = 0; j < 3; ++j )
          result( i, jj ) += M1( i, j ) * M2.q[jj * 3 + j];
      };
    return result;
  };

  /**
   * Multiplication with a column vector.
   * \test PASSED
   */
  friend vect< value_type, 3 > operator*(const self& R, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( R.q[0] * V[0] + R.q[3] * V[1] + R.q[6] * V[2],
                                  R.q[1] * V[0] + R.q[4] * V[1] + R.q[7] * V[2],
                                  R.q[2] * V[0] + R.q[5] * V[1] + R.q[8] * V[2] );
  };

  friend vect< value_type, 3 > operator*(const vect< value_type, 3 >& V, const self& R)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( R.q[0] * V[0] + R.q[1] * V[1] + R.q[2] * V[2],
                                  R.q[3] * V[0] + R.q[4] * V[1] + R.q[5] * V[2],
                                  R.q[6] * V[0] + R.q[7] * V[1] + R.q[8] * V[2] );
  };


  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality operator for a rotation matrix.
   * \test PASSED
   */
  friend bool operator==( const self& M1, const self& M2 ) BOOST_NOEXCEPT {
    return ( ( M1.q[0] == M2.q[0] ) && ( M1.q[1] == M2.q[1] ) && ( M1.q[2] == M2.q[2] ) && ( M1.q[3] == M2.q[3] )
             && ( M1.q[4] == M2.q[4] ) && ( M1.q[5] == M2.q[5] ) );
  };

  /**
   * Inequality operator for a rotation matrix.
   * \test PASSED
   */
  friend bool operator!=( const self& M1, const self& M2 ) BOOST_NOEXCEPT {
    return ( ( M1.q[0] != M2.q[0] ) || ( M1.q[1] != M2.q[1] ) || ( M1.q[2] != M2.q[2] ) || ( M1.q[3] != M2.q[3] )
             || ( M1.q[4] != M2.q[4] ) || ( M1.q[5] != M2.q[5] ) );
  };


  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Produces a transpose matrix which is the inverse rotation.
   * \test PASSED
   */
  friend self transpose( const self& R ) BOOST_NOEXCEPT {
    return self( R.q[0], R.q[1], R.q[2], R.q[3], R.q[4], R.q[5], R.q[6], R.q[7], R.q[8] );
  };

  friend self transpose_move( const self& R ) BOOST_NOEXCEPT {
    return self( R.q[0], R.q[1], R.q[2], R.q[3], R.q[4], R.q[5], R.q[6], R.q[7], R.q[8] );
  };

  /**
   * Produces a cofactor matrix which is the same as the rotation matrix itself.
   * \test PASSED
   */
  friend self cofactor( const self& R ) BOOST_NOEXCEPT { return R; };

  /**
   * Invert the transformation.
   * \test PASSED
   */
  friend self invert( const self& R ) BOOST_NOEXCEPT {
    return self( R.q[0], R.q[1], R.q[2], R.q[3], R.q[4], R.q[5], R.q[6], R.q[7], R.q[8] );
  };

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace( const self& R ) BOOST_NOEXCEPT { return R.q[0] + R.q[4] + R.q[8]; };

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant( const self& ) BOOST_NOEXCEPT { return value_type( 1.0 ); };

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::symmetric > getSymPart() const {
    return mat< value_type, mat_structure::symmetric >( *this );
  };

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::skew_symmetric > getSkewSymPart() const {
    return mat< value_type, mat_structure::skew_symmetric >( *this );
  };


  /// Loading a rot_mat_3D value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, rot_mat_3D< T >& >& R ) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r11", R.second.q[0] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r21", R.second.q[1] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r31", R.second.q[2] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r12", R.second.q[3] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r22", R.second.q[4] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r32", R.second.q[5] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r13", R.second.q[6] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r23", R.second.q[7] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_r33", R.second.q[8] );
  };

  /// Loading a rot_mat_3D value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, rot_mat_3D< T >& R ) {
    return in & RK_SERIAL_LOAD_WITH_NAME( R );
  };

  /// Saving a rot_mat_3D value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const rot_mat_3D< T >& >& R ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r11", R.second.q[0] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r21", R.second.q[1] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r31", R.second.q[2] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r12", R.second.q[3] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r22", R.second.q[4] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r32", R.second.q[5] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r13", R.second.q[6] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r23", R.second.q[7] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_r33", R.second.q[8] );
  };

  /// Saving a rot_mat_3D value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const rot_mat_3D< T >& R ) {
    return out & RK_SERIAL_SAVE_WITH_NAME( R );
  };
};

namespace rtti {

template < typename T >
struct get_type_id< rot_mat_3D< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000018 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::rot_mat_3D" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::rot_mat_3D"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const rot_mat_3D< T >& save_type;
  typedef rot_mat_3D< T >& load_type;
};
};


/**
 * Prints a rotation matrix to a standard output stream (<<) as
 * "((a11; a12; a13); (a21; a22; a23); (a31; a32; a33))".
 * \test PASSED
 */
template < class T >
std::ostream& operator<<( std::ostream& out_stream, const rot_mat_3D< T >& R ) {
  return out_stream << "((" << R( 0, 0 ) << "; " << R( 0, 1 ) << "; " << R( 0, 2 ) << "); (" << R( 1, 0 ) << "; "
                    << R( 1, 1 ) << "; " << R( 1, 2 ) << "); (" << R( 2, 0 ) << "; " << R( 2, 1 ) << "; " << R( 2, 2 )
                    << "))";
};


template < typename T >
struct is_readable_matrix< rot_mat_3D< T > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< rot_mat_3D< T > > type;
};


/**
 * This class represents a rotation using quaternions (or Euler-Rodriguez parameters).
 * The convention used is with the leading scalar.
 * \test All tests for this class have been passed!
 */
template < typename T >
class quaternion {
public:
  typedef quaternion< T > self;
  typedef void allocator_type;

  typedef T value_type;
  typedef void container_type;

  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

private:
  value_type q[4];

  quaternion( const_reference q0, const_reference q1, const_reference q2, const_reference q3 ) BOOST_NOEXCEPT {
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
    return;
  };

public:
  friend class euler_angles_TB< value_type >;
  friend class axis_angle< value_type >;
  friend class trans_mat_3D< value_type >;
  friend class unit_quat< value_type >;

  class xrot {
  private:
    value_type q0;
    value_type qx;

    xrot( const_reference Q0, const_reference QX ) BOOST_NOEXCEPT : q0( Q0 ), qx( QX ){};

  public:
    friend class quaternion< value_type >; // befriend parent.

    xrot( const_reference ang = value_type( 0.0 ) ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qx = sin( ang * 0.5 );
    };

    value_type s() const BOOST_NOEXCEPT { return q0; };
    value_type v() const BOOST_NOEXCEPT { return qx; };

    value_type get_angle() const BOOST_NOEXCEPT {
      using std::atan2;
      return value_type( 2.0 ) * atan2( qx, q0 );
    };

    void set_angle( const_reference ang ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qx = sin( ang * 0.5 );
    };

    vect< value_type, 3 > get_axis() const BOOST_NOEXCEPT {
      return vect< value_type, 3 >( value_type( 1.0 ), value_type( 0.0 ), value_type( 0.0 ) );
    };

    quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT {
      return quaternion< value_type >( q0, qx, value_type( 0.0 ), value_type( 0.0 ) );
    };

    xrot& operator*=( const xrot& Q2 ) BOOST_NOEXCEPT {
      value_type tmp = Q2.q0 * q0 - Q2.qx * qx;
      qx = Q2.q0 * qx + Q2.qx * q0;
      q0 = tmp;
      return *this;
    };

    friend xrot operator*(const xrot& Q1, const xrot& Q2)BOOST_NOEXCEPT {
      return xrot( Q2.q0 * Q1.q0 - Q2.qx * Q1.qx, Q2.q0 * Q1.qx + Q2.qx * Q1.q0 );
    };

    friend vect< value_type, 3 > operator*(const xrot& Q, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
      value_type t0 = Q.q0 * Q.qx;
      value_type t3 = -Q.qx * Q.qx;
      return vect< value_type, 3 >( V[0], value_type( 2.0 ) * ( t3 * V[1] - t0 * V[2] ) + V[1],
                                    value_type( 2.0 ) * ( t0 * V[1] + t3 * V[2] ) + V[2] );
    };
  };

  class yrot {
  private:
    value_type q0;
    value_type qy;

    yrot( const_reference Q0, const_reference QY ) BOOST_NOEXCEPT : q0( Q0 ), qy( QY ){};

  public:
    friend class quaternion< value_type >; // befriend parent.

    yrot( const_reference ang = value_type( 0.0 ) ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qy = sin( ang * 0.5 );
    };

    value_type s() const BOOST_NOEXCEPT { return q0; };
    value_type v() const BOOST_NOEXCEPT { return qy; };

    value_type get_angle() const BOOST_NOEXCEPT {
      using std::atan2;
      return value_type( 2.0 ) * atan2( qy, q0 );
    };

    void set_angle( const_reference ang ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qy = sin( ang * 0.5 );
    };

    vect< value_type, 3 > get_axis() const BOOST_NOEXCEPT {
      return vect< value_type, 3 >( value_type( 0.0 ), value_type( 1.0 ), value_type( 0.0 ) );
    };

    quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT {
      return quaternion< value_type >( q0, value_type( 0.0 ), qy, value_type( 0.0 ) );
    };

    yrot& operator*=( const yrot& Q2 ) BOOST_NOEXCEPT {
      value_type tmp = Q2.q0 * q0 - Q2.qy * qy;
      qy = Q2.q0 * qy + Q2.qy * q0;
      q0 = tmp;
      return *this;
    };

    friend yrot operator*(const yrot& Q1, const yrot& Q2)BOOST_NOEXCEPT {
      return yrot( Q2.q0 * Q1.q0 - Q2.qy * Q1.qy, Q2.q0 * Q1.qy + Q2.qy * Q1.q0 );
    };

    friend vect< value_type, 3 > operator*(const yrot& Q, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
      value_type t1 = Q.q0 * Q.qy;
      value_type t6 = -Q.qy * Q.qy;
      return vect< value_type, 3 >( value_type( 2.0 ) * ( t6 * V[0] + t1 * V[2] ) + V[0], V[1],
                                    value_type( 2.0 ) * ( -t1 * V[0] + t6 * V[2] ) + V[2] );
    };
  };

  class zrot {
  private:
    value_type q0;
    value_type qz;

    zrot( const_reference Q0, const_reference QZ ) BOOST_NOEXCEPT : q0( Q0 ), qz( QZ ){};

  public:
    friend class quaternion< value_type >; // befriend parent.

    zrot( const_reference ang = value_type( 0.0 ) ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qz = sin( ang * 0.5 );
    };

    value_type s() const BOOST_NOEXCEPT { return q0; };
    value_type v() const BOOST_NOEXCEPT { return qz; };

    value_type get_angle() const BOOST_NOEXCEPT {
      using std::atan2;
      return value_type( 2.0 ) * atan2( qz, q0 );
    };

    void set_angle( const_reference ang ) BOOST_NOEXCEPT {
      using std::sin;
      using std::cos;
      q0 = cos( ang * 0.5 );
      qz = sin( ang * 0.5 );
    };

    vect< value_type, 3 > get_axis() const BOOST_NOEXCEPT {
      return vect< value_type, 3 >( value_type( 0.0 ), value_type( 0.0 ), value_type( 1.0 ) );
    };

    quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT {
      return quaternion< value_type >( q0, value_type( 0.0 ), value_type( 0.0 ), qz );
    };

    zrot& operator*=( const zrot& Q2 ) BOOST_NOEXCEPT {
      value_type tmp = Q2.q0 * q0 - Q2.qz * qz;
      qz = Q2.q0 * qz + Q2.qz * q0;
      q0 = tmp;
      return *this;
    };

    friend zrot operator*(const zrot& Q1, const zrot& Q2)BOOST_NOEXCEPT {
      return zrot( Q2.q0 * Q1.q0 - Q2.qz * Q1.qz, Q2.q0 * Q1.qz + Q2.qz * Q1.q0 );
    };

    friend vect< value_type, 3 > operator*(const zrot& Q, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
      value_type t2 = Q.q0 * Q.qz;
      value_type t8 = -Q.qz * Q.qz;
      return vect< value_type, 3 >( value_type( 2.0 ) * ( t8 * V[0] - t2 * V[1] ) + V[0],
                                    value_type( 2.0 ) * ( t2 * V[0] + t8 * V[1] ) + V[1], V[2] );
    };
  };

  friend self operator*(const xrot& q1, const yrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), q2.s() * q1.v(), q2.v() * q1.s(), q2.v() * q1.v() );
  };

  friend self operator*(const yrot& q1, const xrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), q2.v() * q1.s(), q2.s() * q1.v(), -q2.v() * q1.v() );
  };

  friend self operator*(const zrot& q1, const xrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), q2.v() * q1.s(), q2.v() * q1.v(), q2.s() * q1.v() );
  };

  friend self operator*(const xrot& q1, const zrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), q2.s() * q1.v(), -q2.v() * q1.v(), q2.v() * q1.s() );
  };

  friend self operator*(const yrot& q1, const zrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), q2.v() * q1.v(), q2.s() * q1.v(), q2.v() * q1.s() );
  };

  friend self operator*(const zrot& q1, const yrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.s(), -q2.v() * q1.v(), q2.v() * q1.s(), q2.s() * q1.v() );
  };

  friend self operator*(const self& q1, const xrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.q[0] - q2.v() * q1.q[1], q2.s() * q1.q[1] + q2.v() * q1.q[0],
                 q2.s() * q1.q[2] + q2.v() * q1.q[3], q2.s() * q1.q[3] - q2.v() * q1.q[2] );
  };

  friend self operator*(const self& q1, const yrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.q[0] - q2.v() * q1.q[2], q2.s() * q1.q[1] - q2.v() * q1.q[3],
                 q2.s() * q1.q[2] + q2.v() * q1.q[0], q2.s() * q1.q[3] + q2.v() * q1.q[1] );
  };

  friend self operator*(const self& q1, const zrot& q2)BOOST_NOEXCEPT {
    return self( q2.s() * q1.q[0] - q2.v() * q1.q[3], q2.s() * q1.q[1] + q2.v() * q1.q[2],
                 q2.s() * q1.q[2] - q2.v() * q1.q[1], q2.s() * q1.q[3] + q2.v() * q1.q[0] );
  };

  self& operator*=( const xrot& q2 ) BOOST_NOEXCEPT {
    value_type tmp = q2.s() * q[0] - q2.v() * q[1];
    q[1] = q2.s() * q[1] + q2.v() * q[0];
    q[0] = tmp;
    tmp = q2.s() * q[2] + q2.v() * q[3];
    q[3] = q2.s() * q[3] - q2.v() * q[2];
    q[2] = tmp;
    return *this;
  };

  self& operator*=( const yrot& q2 ) BOOST_NOEXCEPT {
    value_type tmp = q2.s() * q[0] - q2.v() * q[2];
    q[2] = q2.s() * q[2] + q2.v() * q[0];
    q[0] = tmp;
    tmp = q2.s() * q[1] - q2.v() * q[3];
    q[3] = q2.s() * q[3] + q2.v() * q[1];
    q[1] = tmp;
    return *this;
  };

  self& operator*=( const zrot& q2 ) BOOST_NOEXCEPT {
    value_type tmp = q2.s() * q[0] - q2.v() * q[3];
    q[3] = q2.s() * q[3] + q2.v() * q[0];
    q[0] = tmp;
    tmp = q2.s() * q[1] + q2.v() * q[2];
    q[2] = q2.s() * q[2] - q2.v() * q[1];
    q[1] = tmp;
    return *this;
  };

  friend self operator*(const xrot& q1, const self& q2)BOOST_NOEXCEPT {
    return self( q2.q[0] * q1.s() - q2.q[1] * q1.v(), q2.q[1] * q1.s() + q2.q[0] * q1.v(),
                 q2.q[2] * q1.s() - q2.q[3] * q1.v(), q2.q[3] * q1.s() + q2.q[2] * q1.v() );
  };

  friend self operator*(const yrot& q1, const self& q2)BOOST_NOEXCEPT {
    return self( q2.q[0] * q1.s() - q2.q[2] * q1.v(), q2.q[1] * q1.s() + q2.q[3] * q1.v(),
                 q2.q[2] * q1.s() + q2.q[0] * q1.v(), q2.q[3] * q1.s() - q2.q[1] * q1.v() );
  };

  friend self operator*(const zrot& q1, const self& q2)BOOST_NOEXCEPT {
    return self( q2.q[0] * q1.s() - q2.q[3] * q1.v(), q2.q[1] * q1.s() - q2.q[2] * q1.v(),
                 q2.q[2] * q1.s() + q2.q[1] * q1.v(), q2.q[3] * q1.s() + q2.q[0] * q1.v() );
  };


  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default Constructor.
   * \test PASSED
   */
  quaternion() BOOST_NOEXCEPT {
    q[0] = 1.0;
    q[1] = 0.0;
    q[2] = 0.0;
    q[3] = 0.0;
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Copy-constructor.
   * \test PASSED
   */
  quaternion( const self& Q ) BOOST_NOEXCEPT { std::copy( Q.q, Q.q + 4, q ); };
#else
  quaternion( const self& ) BOOST_NOEXCEPT = default;
#endif

  template < typename Vector >
  explicit quaternion( const Vector& aV,
                       typename boost::enable_if_c< is_readable_vector< Vector >::value, void* >::type dummy
                       = nullptr ) BOOST_NOEXCEPT {
    RK_UNUSED( dummy );
    vect< value_type, 4 > v = unit( vect< value_type, 4 >( aV[0], aV[1], aV[2], aV[3] ) );
    q[0] = v[0];
    q[1] = v[1];
    q[2] = v[2];
    q[3] = v[3];
  };

  /**
   * Constructor from a rotation matrix.
   * \test PASSED
   */
  explicit quaternion( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    using std::sqrt;
    value_type tra = R.q[0] + R.q[4] + R.q[8];
    if( tra > 0.01 ) {
      q[0] = value_type( 0.5 ) * sqrt( value_type( 1.0 ) + tra );
      q[1] = value_type( 0.25 ) * ( R.q[5] - R.q[7] ) / q[0];
      q[2] = value_type( 0.25 ) * ( R.q[6] - R.q[2] ) / q[0];
      q[3] = value_type( 0.25 ) * ( R.q[1] - R.q[3] ) / q[0];
    } else if( ( R.q[0] > R.q[4] ) && ( R.q[0] > R.q[8] ) ) {
      q[1] = value_type( 0.5 ) * sqrt( value_type( 1.0 ) + R.q[0] - R.q[4] - R.q[8] );
      q[0] = value_type( 0.25 ) * ( R.q[7] - R.q[5] ) / q[1];
      q[2] = value_type( 0.25 ) * ( R.q[3] + R.q[1] ) / q[1];
      q[3] = value_type( 0.25 ) * ( R.q[6] + R.q[2] ) / q[1];
    } else if( R.q[4] > R.q[8] ) {
      q[2] = value_type( 0.5 ) * sqrt( value_type( 1.0 ) + R.q[4] - R.q[0] - R.q[8] );
      q[0] = value_type( 0.25 ) * ( R.q[6] - R.q[2] ) / q[2];
      q[1] = value_type( 0.25 ) * ( R.q[3] + R.q[1] ) / q[2];
      q[3] = value_type( 0.25 ) * ( R.q[7] + R.q[5] ) / q[2];
    } else {
      q[3] = value_type( 0.5 ) * sqrt( value_type( 1.0 ) + R.q[8] - R.q[0] - R.q[4] );
      q[0] = value_type( 0.25 ) * ( R.q[3] - R.q[1] ) / q[3];
      q[1] = value_type( 0.25 ) * ( R.q[6] + R.q[2] ) / q[3];
      q[2] = value_type( 0.25 ) * ( R.q[7] + R.q[5] ) / q[3];
    };
    return;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides the rotation matrix as an ordinary 3x3 matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::square > getMat() const {
    value_type t01( value_type( 2.0 ) * q[0] * q[1] );
    value_type t02( value_type( 2.0 ) * q[0] * q[2] );
    value_type t03( value_type( 2.0 ) * q[0] * q[3] );
    value_type t11( value_type( 2.0 ) * q[1] * q[1] );
    value_type t12( value_type( 2.0 ) * q[1] * q[2] );
    value_type t13( value_type( 2.0 ) * q[1] * q[3] );
    value_type t22( value_type( 2.0 ) * q[2] * q[2] );
    value_type t23( value_type( 2.0 ) * q[2] * q[3] );
    value_type t33( value_type( 2.0 ) * q[3] * q[3] );
    return mat< value_type, mat_structure::square >( value_type( 1.0 ) - t22 - t33, t12 - t03, t02 + t13, t12 + t03,
                                                     value_type( 1.0 ) - t11 - t33, t23 - t01, t13 - t02, t01 + t23,
                                                     value_type( 1.0 ) - t11 - t22 );
  };

  /**
   * Provides the rotation matrix corresponding to the quaternion.
   * \test PASSED
   */
  rot_mat_3D< value_type > getRotMat() const BOOST_NOEXCEPT {
    value_type t01( value_type( 2.0 ) * q[0] * q[1] );
    value_type t02( value_type( 2.0 ) * q[0] * q[2] );
    value_type t03( value_type( 2.0 ) * q[0] * q[3] );
    value_type t11( value_type( 2.0 ) * q[1] * q[1] );
    value_type t12( value_type( 2.0 ) * q[1] * q[2] );
    value_type t13( value_type( 2.0 ) * q[1] * q[3] );
    value_type t22( value_type( 2.0 ) * q[2] * q[2] );
    value_type t23( value_type( 2.0 ) * q[2] * q[3] );
    value_type t33( value_type( 2.0 ) * q[3] * q[3] );
    return rot_mat_3D< value_type >( value_type( 1.0 ) - t22 - t33, t12 - t03, t02 + t13, t12 + t03,
                                     value_type( 1.0 ) - t11 - t33, t23 - t01, t13 - t02, t01 + t23,
                                     value_type( 1.0 ) - t11 - t22 );
  };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT {
    assert( i < 4 );
    return q[i];
  };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Assignment operator from a quaternion.
   * \test PASSED
   */
  self& operator=( const self& Q ) BOOST_NOEXCEPT {
    std::copy( Q.q, Q.q + 4, q );
    return *this;
  };
#else
  self& operator=( const self& ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Assignment operator from a rotation matrix.
   * \test PASSED
   */
  self& operator=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT { return ( *this = self( R ) ); };

  /**
   * Assignment operator from a euler angles TB representation.
   * \test PASSED
   */
  self& operator=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { return ( *this = E.getQuaternion() ); };

  /**
   * Assignment operator from an axis / angle representation.
   * \test PASSED
   */
  self& operator=( const axis_angle< value_type >& A ) BOOST_NOEXCEPT { return ( *this = A.getQuaternion() ); };

  /**
   * Multiply-and-store operator from a quaternion.
   * \test PASSED
   */
  self& operator*=( const self& Q ) BOOST_NOEXCEPT { return ( *this = *this * Q ); };

  /**
   * Multiply-and-store operator from a rotation matrix.
   * \test PASSED
   */
  self& operator*=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT { return ( *this *= self( R ) ); };

  /**
   * Multiply-and-store operator from a euler angles TB representation.
   * \test PASSED
   */
  self& operator*=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { return ( *this *= E.getQuaternion() ); };

  /**
   * Multiply-and-store operator from an axis / angle representation.
   * \test PASSED
   */
  self& operator*=( const axis_angle< value_type >& A ) BOOST_NOEXCEPT { return ( *this *= A.getQuaternion() ); };


  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiplication by a quaternion.
   * \test PASSED
   */
  friend self operator*(const self& Q1, const self& Q2)BOOST_NOEXCEPT {
    return self( Q2.q[0] * Q1.q[0] - Q2.q[1] * Q1.q[1] - Q2.q[2] * Q1.q[2] - Q2.q[3] * Q1.q[3],
                 Q2.q[0] * Q1.q[1] + Q2.q[3] * Q1.q[2] - Q2.q[2] * Q1.q[3] + Q2.q[1] * Q1.q[0],
                 Q2.q[0] * Q1.q[2] - Q2.q[3] * Q1.q[1] + Q2.q[1] * Q1.q[3] + Q2.q[2] * Q1.q[0],
                 Q2.q[0] * Q1.q[3] + Q2.q[2] * Q1.q[1] - Q2.q[1] * Q1.q[2] + Q2.q[3] * Q1.q[0] );
  };

  /**
   * Multiplication by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const self& Q, const Matrix& M ) {
    return Q.getRotMat() * M;
  };

  /**
   * Multiplication by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const Matrix& M, const self& Q ) {
    return M * Q.getRotMat();
  };

  /**
   * Multiplication by a column vector.
   * \test PASSED
   */
  friend vect< value_type, 3 > operator*(const self& Q, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
    value_type t[9];
    t[0] = Q.q[0] * Q.q[1];
    t[1] = Q.q[0] * Q.q[2];
    t[2] = Q.q[0] * Q.q[3];
    t[3] = -Q.q[1] * Q.q[1];
    t[4] = Q.q[1] * Q.q[2];
    t[5] = Q.q[1] * Q.q[3];
    t[6] = -Q.q[2] * Q.q[2];
    t[7] = Q.q[2] * Q.q[3];
    t[8] = -Q.q[3] * Q.q[3];
    return vect< T, 3 >(
      value_type( 2.0 ) * ( ( t[6] + t[8] ) * V[0] + ( t[4] - t[2] ) * V[1] + ( t[1] + t[5] ) * V[2] ) + V[0],
      value_type( 2.0 ) * ( ( t[2] + t[4] ) * V[0] + ( t[3] + t[8] ) * V[1] + ( t[7] - t[0] ) * V[2] ) + V[1],
      value_type( 2.0 ) * ( ( t[5] - t[1] ) * V[0] + ( t[0] + t[7] ) * V[1] + ( t[3] + t[6] ) * V[2] ) + V[2] );
  };

  friend vect< value_type, 3 > operator*(const self& Q, const vect_component< value_type, 0 >& x_value)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( x_value.q - 2.0 * ( Q.q[2] * Q.q[2] + Q.q[3] * Q.q[3] ) * x_value.q,
                                  2.0 * ( Q.q[0] * Q.q[3] + Q.q[1] * Q.q[2] ) * x_value.q,
                                  2.0 * ( Q.q[1] * Q.q[3] - Q.q[0] * Q.q[2] ) * x_value.q );
  };

  friend vect< value_type, 3 > operator*(const self& Q, const vect_component< value_type, 1 >& y_value)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( 2.0 * ( Q.q[1] * Q.q[2] - Q.q[0] * Q.q[3] ) * y_value.q,
                                  y_value.q - 2.0 * ( Q.q[1] * Q.q[1] + Q.q[3] * Q.q[3] ) * y_value.q,
                                  2.0 * ( Q.q[0] * Q.q[1] + Q.q[2] * Q.q[3] ) * y_value.q );
  };

  friend vect< value_type, 3 > operator*(const self& Q, const vect_component< value_type, 2 >& z_value)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( 2.0 * ( Q.q[0] * Q.q[2] + Q.q[1] * Q.q[3] ) * z_value.q,
                                  2.0 * ( Q.q[2] * Q.q[3] - Q.q[0] * Q.q[1] ) * z_value.q,
                                  z_value.q - 2.0 * ( Q.q[1] * Q.q[1] + Q.q[2] * Q.q[2] ) * z_value.q );
  };

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality operator with a quaternion representation.
   * \test PASSED
   */
  friend bool operator==( const self& Q1, const self& Q2 ) BOOST_NOEXCEPT {
    return ( ( Q1.q[0] == Q2.q[0] ) && ( Q1.q[1] == Q2.q[1] ) && ( Q1.q[2] == Q2.q[2] ) && ( Q1.q[3] == Q2.q[3] ) );
  };

  /**
   * Inequality operator with a quaternion representation.
   * \test PASSED
   */
  friend bool operator!=( const self& Q1, const self& Q2 ) BOOST_NOEXCEPT {
    return ( ( Q1.q[0] != Q2.q[0] ) || ( Q1.q[1] != Q2.q[1] ) || ( Q1.q[2] != Q2.q[2] ) || ( Q1.q[3] != Q2.q[3] ) );
  };


  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /**
   * Gets the time-derivative of the quaternion that corresponds to the angular velocity Omega.
   * \test PASSED
   */
  vect< value_type, 4 > getQuaternionDot( const vect< value_type, 3 >& Omega ) const BOOST_NOEXCEPT {
    return vect< value_type, 4 >( -value_type( 0.5 ) * ( q[1] * Omega.q[0] + q[2] * Omega.q[1] + q[3] * Omega.q[2] ),
                                  value_type( 0.5 ) * ( q[0] * Omega.q[0] - q[3] * Omega.q[1] + q[2] * Omega.q[2] ),
                                  value_type( 0.5 ) * ( q[0] * Omega.q[1] + q[3] * Omega.q[0] - q[1] * Omega.q[2] ),
                                  value_type( 0.5 ) * ( q[0] * Omega.q[2] - q[2] * Omega.q[0] + q[1] * Omega.q[1] ) );
  };

  /**
   * Gets the angular velocity that corresponds to the time-derivative of the quaternion.
   * \test PASSED
   */
  vect< value_type, 3 > getOmega( const vect< value_type, 4 >& QuaternionDot ) const BOOST_NOEXCEPT {
    return vect< value_type, 3 >( value_type( 2.0 ) * ( -q[1] * QuaternionDot.q[0] + q[0] * QuaternionDot.q[1]
                                                        + q[3] * QuaternionDot.q[2] - q[2] * QuaternionDot.q[3] ),
                                  value_type( 2.0 ) * ( -q[2] * QuaternionDot.q[0] - q[3] * QuaternionDot.q[1]
                                                        + q[0] * QuaternionDot.q[2] + q[1] * QuaternionDot.q[3] ),
                                  value_type( 2.0 ) * ( -q[3] * QuaternionDot.q[0] + q[2] * QuaternionDot.q[1]
                                                        - q[1] * QuaternionDot.q[2] + q[0] * QuaternionDot.q[3] ) );
  };

  /**
   * Gets the 2-time-derivative of the quaternion that corresponds to the angular velocity Omega.
   * \test PASSED
   */
  vect< value_type, 4 > getQuaternionDotDot( const vect< value_type, 4 >& QD, const vect< value_type, 3 >& W,
                                             const vect< value_type, 3 >& WD ) const BOOST_NOEXCEPT {
    return vect< value_type, 4 >( -value_type( 0.5 ) * ( q[1] * WD.q[0] + q[2] * WD.q[1] + q[3] * WD.q[2]
                                                         + QD.q[1] * W.q[0] + QD.q[2] * W.q[1] + QD.q[3] * W.q[2] ),
                                  value_type( 0.5 ) * ( q[0] * WD.q[0] - q[3] * WD.q[1] + q[2] * WD.q[2]
                                                        + QD.q[0] * W.q[0] - QD.q[3] * W.q[1] + QD.q[2] * W.q[2] ),
                                  value_type( 0.5 ) * ( q[3] * WD.q[0] + q[0] * WD.q[1] - q[1] * WD.q[2]
                                                        + QD.q[3] * W.q[0] + QD.q[0] * W.q[1] - QD.q[1] * W.q[2] ),
                                  value_type( 0.5 ) * ( -q[2] * WD.q[0] + q[1] * WD.q[1] + q[0] * WD.q[2]
                                                        - QD.q[2] * W.q[0] + QD.q[1] * W.q[1] + QD.q[0] * W.q[2] ) );
  };

  /**
   * Gets the angular acceleration that corresponds to the 2-time-derivative of the quaternion.
   * \test PASSED
   */
  vect< value_type, 3 > getOmegaDot( const vect< value_type, 4 >& QD,
                                     const vect< value_type, 4 >& QDD ) const BOOST_NOEXCEPT {
    return vect< value_type, 3 >(
      value_type( 2.0 ) * ( -q[1] * QDD.q[0] + q[0] * QDD.q[1] + q[3] * QDD.q[2] - q[2] * QDD.q[3] - QD.q[1] * QD.q[0]
                            + QD.q[0] * QD.q[1] + QD.q[3] * QD.q[2] - QD.q[2] * QD.q[3] ),
      value_type( 2.0 ) * ( -q[2] * QDD.q[0] - q[3] * QDD.q[1] + q[0] * QDD.q[2] + q[1] * QDD.q[3] - QD.q[2] * QD.q[0]
                            - QD.q[3] * QD.q[1] + QD.q[0] * QD.q[2] + QD.q[1] * QD.q[3] ),
      value_type( 2.0 ) * ( -q[3] * QDD.q[0] + q[2] * QDD.q[1] - q[1] * QDD.q[2] + q[0] * QDD.q[3] - QD.q[3] * QD.q[0]
                            + QD.q[2] * QD.q[1] - QD.q[1] * QD.q[2] + QD.q[0] * QD.q[3] ) );
  };

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Produces a transpose quaternion which is the inverse rotation.
   * \test PASSED
   */
  friend self transpose( const self& Q ) BOOST_NOEXCEPT { return self( Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3] ); };

  /**
   * Produces a transpose quaternion which is the inverse rotation.
   * \test PASSED
   */
  friend self transpose_move( const self& Q ) BOOST_NOEXCEPT { return self( Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3] ); };

  /**
   * Produces a cofactor matrix which is the same as the rotation matrix itself.
   * \test PASSED
   */
  friend self cofactor( const self& Q ) BOOST_NOEXCEPT { return Q; };

  /**
   * Invert the rotation.
   * \test PASSED
   */
  friend self invert( const self& Q ) BOOST_NOEXCEPT { return self( Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3] ); };

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace( const self& Q ) BOOST_NOEXCEPT {
    return value_type( 4.0 ) * Q.q[0] * Q.q[0] - value_type( 1.0 );
  };

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant( const self& Q ) BOOST_NOEXCEPT { return value_type( 1.0 ); };

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::symmetric > getSymPart() const {
    value_type t11( value_type( 2.0 ) * q[1] * q[1] );
    value_type t12( value_type( 2.0 ) * q[1] * q[2] );
    value_type t13( value_type( 2.0 ) * q[1] * q[3] );
    value_type t22( value_type( 2.0 ) * q[2] * q[2] );
    value_type t23( value_type( 2.0 ) * q[2] * q[3] );
    value_type t33( value_type( 2.0 ) * q[3] * q[3] );
    return mat< value_type, mat_structure::symmetric >(
      value_type( 1.0 ) - t22 - t33, t12, t13, value_type( 1.0 ) - t11 - t33, t23, value_type( 1.0 ) - t11 - t22 );
  };

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::skew_symmetric > getSkewSymPart() const {
    value_type t01( value_type( 2.0 ) * q[0] * q[1] );
    value_type t02( value_type( 2.0 ) * q[0] * q[2] );
    value_type t03( value_type( 2.0 ) * q[0] * q[3] );
    return mat< value_type, mat_structure::skew_symmetric >( -t03, t02, -t01 );
  };


  /// Loading a quaternion value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, quaternion< T >& >& R ) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_q0", R.second.q[0] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_q1", R.second.q[1] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_q2", R.second.q[2] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_q3", R.second.q[3] );
  };

  /// Loading a quaternion value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, quaternion< T >& R ) {
    return in & RK_SERIAL_LOAD_WITH_NAME( R );
  };

  /// Saving a quaternion value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const quaternion< T >& >& R ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_q0", R.second.q[0] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_q1", R.second.q[1] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_q2", R.second.q[2] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_q3", R.second.q[3] );
  };

  /// Saving a quaternion value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const quaternion< T >& R ) {
    return out & RK_SERIAL_SAVE_WITH_NAME( R );
  };
};

namespace rtti {

template < typename T >
struct get_type_id< quaternion< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x0000001A );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::quaternion" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::quaternion"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const quaternion< T >& save_type;
  typedef quaternion< T >& load_type;
};
};


/**
 * Prints a quaternion to a standard output stream (<<) as "(q0; q1; q2; q3)".
 * \test PASSED
 */
template < class T >
std::ostream& operator<<( std::ostream& out_stream, const quaternion< T >& Q ) {
  return ( out_stream << "(" << Q[0] << "; " << Q[1] << "; " << Q[2] << "; " << Q[3] << ")" );
};


/**
 * This class repressents a rotation using Euler angles (Tait-Bryan), 321-body-fixed, in body frame.
 */
template < class T >
class euler_angles_TB {
public:
  typedef euler_angles_TB< T > self;
  typedef void allocator_type;

  typedef T value_type;
  typedef void container_type;

  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

private:
  value_type q[3];

public:
  friend class axis_angle< value_type >;
  friend class trans_mat_3D< value_type >;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default Constructor.
   * \test PASSED
   */
  euler_angles_TB() BOOST_NOEXCEPT {
    q[0] = 0.0;
    q[1] = 0.0;
    q[2] = 0.0;
  };

  /**
   * Constructor from three euler angles.
   * \test PASSED
   */
  euler_angles_TB( const_reference Yaw_, const_reference Pitch_, const_reference Roll_ ) BOOST_NOEXCEPT {
    q[0] = Yaw_;
    q[1] = Pitch_;
    q[2] = Roll_;
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Copy-constructor.
   * \test PASSED
   */
  euler_angles_TB( const self& E ) BOOST_NOEXCEPT {
    q[0] = E.q[0];
    q[1] = E.q[1];
    q[2] = E.q[2];
  };
#else
  euler_angles_TB( const self& ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Constructor from a quaternion.
   * \test PASSED
   */
  explicit euler_angles_TB( const quaternion< value_type >& Q ) BOOST_NOEXCEPT {
    using std::asin;
    using std::atan2;
    using std::cos;
    q[1] = value_type( 2.0 ) * ( Q.q[0] * Q.q[2] - Q.q[1] * Q.q[3] );
    if( ( q[1] != value_type( 1.0 ) ) && ( q[1] != value_type( -1.0 ) ) ) {
      q[1] = asin( q[1] );
      value_type cp = value_type( 1.0 ) / cos( q[1] );
      q[2] = atan2( value_type( 2.0 ) * cp * ( Q.q[2] * Q.q[3] + Q.q[0] * Q.q[1] ),
                    cp * ( value_type( 1.0 ) - value_type( 2.0 ) * ( Q.q[1] * Q.q[1] + Q.q[2] * Q.q[2] ) ) );
      q[0] = atan2( value_type( 2.0 ) * cp * ( Q.q[1] * Q.q[2] + Q.q[0] * Q.q[3] ),
                    cp * ( value_type( 1.0 ) - value_type( 2.0 ) * ( Q.q[2] * Q.q[2] + Q.q[3] * Q.q[3] ) ) );
    } else {
      q[0] = value_type( 0.0 );
      q[2] = atan2( q[1] * value_type( 2.0 ) * ( Q.q[1] * Q.q[2] - Q.q[0] * Q.q[3] ),
                    q[1] * value_type( 2.0 ) * ( Q.q[1] * Q.q[3] + Q.q[0] * Q.q[2] ) );
      q[1] *= value_type( 1.57079632679489662 );
    };
    return;
  };

  /**
   * Constructor from a rotation matrix.
   * \test PASSED
   */
  explicit euler_angles_TB( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    using std::asin;
    using std::atan2;
    using std::cos;
    if( ( R.q[2] != value_type( 1.0 ) ) && ( R.q[2] != value_type( -1.0 ) ) ) {
      q[1] = asin( -R.q[2] );
      value_type cp = value_type( 1.0 ) / cos( q[1] );
      q[2] = atan2( cp * R.q[5], cp * R.q[8] );
      q[0] = atan2( cp * R.q[1], cp * R.q[0] );
    } else {
      q[0] = value_type( 0.0 );
      q[2] = atan2( -R.q[2] * R.q[3], -R.q[2] * R.q[6] );
      q[1] = -R.q[2] * value_type( 1.57079632679489662 );
    };
    return;
  };

  euler_angles_TB( const rot_mat_3D< value_type >& R, const euler_angles_TB< value_type >& Predicted ) BOOST_NOEXCEPT {
    using std::asin;
    using std::atan2;
    using std::cos;
    if( ( R.q[2] != value_type( 1.0 ) ) && ( R.q[2] != value_type( -1.0 ) ) ) {
      q[1] = asin( -R.q[2] );
      value_type cp = value_type( 1.0 ) / cos( q[1] );
      q[2] = atan2( cp * R.q[5], cp * R.q[8] );
      q[0] = atan2( cp * R.q[1], cp * R.q[0] );
    } else {
      q[0] = value_type( 0.0 );
      q[2] = atan2( -R.q[2] * R.q[3], -R.q[2] * R.q[6] );
      q[1] = -R.q[2] * value_type( 1.57079632679489662 );
    };
    return;

    // UNTESTED CODE: ...
    /*T s, c;
    T d1, d2, d3, d4, d5, m1;

    vect<T,3> Turns(floor((Predicted.q[0] + M_PI) / T(2.0*M_PI)), floor((Predicted.q[1] + M_PI) / T(2.0*M_PI)),
    floor((Predicted.q[2] + M_PI) / T(2.0*M_PI)));
    vect<T,3> Pred = vect<T,3>( Predicted.q[0] - Turns.q[0]*T(2.0*M_PI), Predicted.q[1] - Turns.q[1]*T(2.0*M_PI),
    Predicted.q[2] - Turns.q[2]*T(2.0*M_PI));

    vect<T,3> Result(0.0,asin(-m[6]),0.0);
    c = cos(Result.q[1]);

    if(c != 0.0) {
      Result.q[0] = atan2(m[3],m[0]);
      Result.q[2] = atan2(m[7],m[8]);
      d1 = fabs(Result.q[2] - Predicted.q[2]) + fabs(Result.q[0] - Predicted.q[0]);
      d2 = fabs(Result.q[2] + T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] + T(M_PI) - Predicted.q[0]);
      d3 = fabs(Result.q[2] - T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] - T(M_PI) - Predicted.q[0]);
      d4 = fabs(Result.q[2] - T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] + T(M_PI) - Predicted.q[0]);
      d5 = fabs(Result.q[2] + T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] - T(M_PI) - Predicted.q[0]);
      m1 = MIN(MIN(d1,d2),MIN(d3,MIN(d4,d5)));
      if(m1 == d5)
        Result = TVect3<T>(Result.q[0] - T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) +
    Turns.q[1]*T(2.0*M_PI), Result.q[2] + T(M_PI) + Turns.q[2]*T(2.0*M_PI));
      else if(m1 == d4)
        Result = TVect3<T>(Result.q[0] + T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) +
    Turns.q[1]*T(2.0*M_PI), Result.q[2] - T(M_PI) + Turns.q[2]*T(2.0*M_PI));
      else if(m1 == d3)
        Result = TVect3<T>(Result.q[0] - T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) +
    Turns.q[1]*T(2.0*M_PI), Result.q[2] - T(M_PI) + Turns.q[2]*T(2.0*M_PI));
      else if(m1 == d2)
        Result = TVect3<T>(Result.q[0] + T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) +
    Turns.q[1]*T(2.0*M_PI), Result.q[2] + T(M_PI) + Turns.q[2]*T(2.0*M_PI));
      else
        Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] +
    Turns.q[2]*T(2.0*M_PI));
    } else {
      //Singularity Arises here because Roll becomes equivalent to yaw.
      //Here, one of them (Yaw) is set to the prediction.
      Result.q[0] = Predicted.q[0];
      s = sin(Result.q[0]);
      c = cos(Result.q[0]);
      Result.q[2] = atan2(m[2]*s-m[5]*c,m[2]*c+m[5]*s);
      d1 = fabs(Result.q[2] - Predicted.q[2]);
      d2 = fabs(Result.q[2] + T(2.0*M_PI) - Predicted.q[2]);
      d3 = fabs(Result.q[2] - T(2.0*M_PI) - Predicted.q[2]);
      m1 = MIN(d1,MIN(d2,d3));
      if(m1 == d3)
        Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] -
    T(2.0*M_PI) + Turns.q[2]*T(2.0*M_PI));
      else if(m1 == d2)
        Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] +
    T(2.0*M_PI) + Turns.q[2]*T(2.0*M_PI));
      else
        Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] +
    Turns.q[2]*T(2.0*M_PI));
    };
    return Result;*/
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Get yaw, read-write.
   * \test PASSED
   */
  reference yaw() BOOST_NOEXCEPT { return q[0]; };
  /**
   * Get pitch, read-write.
   * \test PASSED
   */
  reference pitch() BOOST_NOEXCEPT { return q[1]; };
  /**
   * Get roll, read-write.
   * \test PASSED
   */
  reference roll() BOOST_NOEXCEPT { return q[2]; };

  /**
   * Get yaw, read-only.
   * \test PASSED
   */
  const_reference yaw() const BOOST_NOEXCEPT { return q[0]; };
  /**
   * Get pitch, read-only.
   * \test PASSED
   */
  const_reference pitch() const BOOST_NOEXCEPT { return q[1]; };
  /**
   * Get roll, read-only.
   * \test PASSED
   */
  const_reference roll() const BOOST_NOEXCEPT { return q[2]; };

  /**
   * Provides a quaternion corresponding to this rotation.
   * \test PASSED
   */
  quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT {
    using std::cos;
    using std::sin;
    value_type cpsi = cos( value_type( 0.5 ) * q[0] );
    value_type spsi = sin( value_type( 0.5 ) * q[0] );
    value_type ctheta = cos( value_type( 0.5 ) * q[1] );

    value_type stheta = sin( value_type( 0.5 ) * q[1] );
    value_type cphi = cos( value_type( 0.5 ) * q[2] );
    value_type sphi = sin( value_type( 0.5 ) * q[2] );

    return quaternion< value_type >(
      cphi * ctheta * cpsi + sphi * stheta * spsi, sphi * ctheta * cpsi - cphi * stheta * spsi,
      cphi * stheta * cpsi + sphi * ctheta * spsi, cphi * ctheta * spsi - sphi * stheta * cpsi );
  };

  /**
   * Provides a rotation matrix corresponding to this rotation.
   * \test PASSED
   */
  rot_mat_3D< value_type > getRotMat() const BOOST_NOEXCEPT {
    using std::sin;
    using std::cos;
    value_type s1( sin( q[0] ) );
    value_type c1( cos( q[0] ) );
    value_type s2( sin( q[1] ) );
    value_type c2( cos( q[1] ) );
    value_type s3( sin( q[2] ) );
    value_type c3( cos( q[2] ) );

    return rot_mat_3D< value_type >( c1 * c2, -( s1 * c3 ) + ( c1 * s2 * s3 ), ( s1 * s3 ) + ( c1 * s2 * c3 ), s1 * c2,
                                     ( c1 * c3 ) + ( s1 * s2 * s3 ), -( c1 * s3 ) + ( s1 * s2 * c3 ), -s2, c2 * s3,
                                     c2 * c3 );
  };

  /**
   * Provides a rotation matrix as a regular 3x3 matris corresponding to this rotation.
   * \test PASSED
   */
  mat< value_type, mat_structure::square > getMat() const {
    using std::sin;
    using std::cos;
    value_type s1( sin( q[0] ) );
    value_type c1( cos( q[0] ) );
    value_type s2( sin( q[1] ) );
    value_type c2( cos( q[1] ) );
    value_type s3( sin( q[2] ) );
    value_type c3( cos( q[2] ) );

    return mat< value_type, mat_structure::square >(
      c1 * c2, -( s1 * c3 ) + ( c1 * s2 * s3 ), ( s1 * s3 ) + ( c1 * s2 * c3 ), s1 * c2, ( c1 * c3 ) + ( s1 * s2 * s3 ),
      -( c1 * s3 ) + ( s1 * s2 * c3 ), -s2, c2 * s3, c2 * c3 );
  };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Standard assignment operator.
   * \test PASSED
   */
  self& operator=( const self& E ) BOOST_NOEXCEPT {
    q[0] = E.q[0];
    q[1] = E.q[1];
    q[2] = E.q[2];
    return *this;
  };
#else
  self& operator=( const self& ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Assignment from a quaternion.
   * \test PASSED
   */
  self& operator=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT { return *this = self( Q ); };

  /**
   * Assignment from a rotation matrix.
   * \test PASSED
   */
  self& operator=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT { return *this = self( R ); };

  /**
   * Assignment from an axis / angle representation.
   * \test PASSED
   */
  self& operator=( const axis_angle< T >& A ) BOOST_NOEXCEPT { return ( *this = A.getEulerAnglesTB() ); };

  /**
   * Multiply-and-store from a euler angles.
   * \test PASSED
   */
  self& operator*=( const self& E ) BOOST_NOEXCEPT { return ( *this = ( this->getRotMat() * E.getRotMat() ) ); };

  /**
   * Multiply-and-store from a rotation matrix.
   * \test PASSED
   */
  self& operator*=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    return ( *this = ( this->getRotMat() * R ) );
  };

  /**
   * Multiply-and-store from a quaternion.
   * \test PASSED
   */
  self& operator*=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT {
    return ( *this = ( this->getQuaternion() * Q ) );
  };

  /**
   * Assignment from an axis / angle representation.
   * \test PASSED
   */
  self& operator*=( const axis_angle< T >& A ) BOOST_NOEXCEPT { return ( *this = this->getRotMat() * A.getRotMat() ); };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiply by a euler angle representation.
   * \test PASSED
   */
  friend rot_mat_3D< value_type > operator*(const self& E1, const self& E2)BOOST_NOEXCEPT {
    return E1.getRotMat() * E2.getRotMat();
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const self& E, const Matrix& M ) {
    return E.getRotMat() * M;
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const Matrix& M, const self& E ) {
    return M * E.getRotMat();
  };

  /**
   * Multiply by a vector, rotating it.
   * \test PASSED
   */
  friend vect< value_type, 3 > operator*(const self& E, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
    return E.getRotMat() * V;
  };

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality comparison operator with a euler angle representation.
   * \test PASSED
   */
  friend bool operator==( const self& E1, const self& E2 ) BOOST_NOEXCEPT {
    return ( ( E1.q[0] == E2.q[0] ) && ( E1.q[1] == E2.q[1] ) && ( E1.q[2] == E2.q[2] ) );
  };

  /**
   * Inequality comparison operator with a euler angle representation.
   * \test PASSED
   */
  friend bool operator!=( const self& E1, const self& E2 ) BOOST_NOEXCEPT {
    return ( ( E1.q[0] != E2.q[0] ) || ( E1.q[1] != E2.q[1] ) || ( E1.q[2] != E2.q[2] ) );
  };


  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/
  /**
   * Produces a transpose quaternion which is the inverse rotation.
   * \test PASSED
   */
  friend self transpose( const self& E ) BOOST_NOEXCEPT {
    using std::sin;
    using std::cos;
    using std::asin;
    using std::atan2;

    value_type s1( sin( E.q[0] ) );
    value_type c1( cos( E.q[0] ) );
    value_type s2( sin( E.q[1] ) );
    value_type c2( cos( E.q[1] ) );
    value_type s3( sin( E.q[2] ) );
    value_type c3( cos( E.q[2] ) );

    value_type R2( ( s1 * s3 ) + ( c1 * s2 * c3 ) );

    self result;

    if( ( R2 != value_type( 1.0 ) ) && ( R2 != value_type( -1.0 ) ) ) {
      value_type R0( c1 * c2 );
      value_type R1( -( s1 * c3 ) + ( c1 * s2 * s3 ) );
      value_type R5( -( c1 * s3 ) + ( s1 * s2 * c3 ) );
      value_type R8( c2 * c3 );
      result.q[1] = asin( -R2 );
      value_type cp = value_type( 1.0 ) / cos( result.q[1] );
      result.q[2] = atan2( cp * R5, cp * R8 );
      result.q[0] = atan2( cp * R1, cp * R0 );
    } else {
      value_type R3( s1 * c2 );
      result.q[0] = value_type( 0.0 );
      result.q[2] = atan2( -R2 * R3, R2 * s2 );
      result.q[1] = -R2 * value_type( 1.57079632679489662 );
    };
    return result;
  };

  /**
   * Produces a cofactor matrix which is the same as the rotation matrix itself.
   * \test PASSED
   */
  friend self cofactor( const self& E ) BOOST_NOEXCEPT { return E; };

  /**
   * Invert the rotation.
   * \test PASSED
   */
  friend self invert( const self& E ) BOOST_NOEXCEPT { return transpose( E ); };

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace( const self& E ) BOOST_NOEXCEPT {
    using std::sin;
    using std::cos;
    value_type t
      = cos( value_type( 0.5 ) * E.q[2] ) * cos( value_type( 0.5 ) * E.q[1] ) * cos( value_type( 0.5 ) * E.q[0] )
        + sin( value_type( 0.5 ) * E.q[2] ) * sin( value_type( 0.5 ) * E.q[1] ) * sin( value_type( 0.5 ) * E.q[0] );
    return value_type( 4.0 ) * t * t - value_type( 1.0 );
  };

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant( const self& ) BOOST_NOEXCEPT { return value_type( 1.0 ); };

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::symmetric > getSymPart() const {
    return mat< value_type, mat_structure::symmetric >( this->getMat() );
  };

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::skew_symmetric > getSkewSymPart() const {
    return mat< value_type, mat_structure::skew_symmetric >( this->getMat() );
  };


  /// Loading a euler_angles_TB value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, euler_angles_TB< T >& >& R ) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_yaw", R.second.q[0] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_pitch", R.second.q[1] )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_roll", R.second.q[2] );
  };

  /// Loading a euler_angles_TB value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, euler_angles_TB< T >& R ) {
    return in & RK_SERIAL_LOAD_WITH_NAME( R );
  };

  /// Saving a euler_angles_TB value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const euler_angles_TB< T >& >& R ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_yaw", R.second.q[0] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_pitch", R.second.q[1] )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_roll", R.second.q[2] );
  };

  /// Saving a euler_angles_TB value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const euler_angles_TB< T >& R ) {
    return out & RK_SERIAL_SAVE_WITH_NAME( R );
  };
};

namespace rtti {

template < typename T >
struct get_type_id< euler_angles_TB< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x0000001B );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::euler_angles_TB" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::euler_angles_TB"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const euler_angles_TB< T >& save_type;
  typedef euler_angles_TB< T >& load_type;
};
};


/**
 * Prints a euler angles to a standard output stream (<<) as "(Yaw = value; Pitch = value; Roll = value)".
 * \test PASSED
 */
template < class T >
std::ostream& operator<<( std::ostream& out_stream, const euler_angles_TB< T >& E ) {
  return out_stream << "(Yaw = " << E.yaw() << "; Pitch = " << E.pitch() << "; Roll = " << E.roll() << ")";
};


/**
 * This class is a 3D rotation represented by an axis and angle.
 * \test All tests for this class have been passed!
 */
template < class T >
class axis_angle {
public:
  typedef axis_angle< T > self;
  typedef void allocator_type;

  typedef T value_type;
  typedef void container_type;

  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

private:
  value_type mAngle;
  vect< value_type, 3 > mAxis;

public:
  friend class trans_mat_3D< value_type >;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default Constructor.
   * \test PASSED
   */
  axis_angle() BOOST_NOEXCEPT : mAngle( 0.0 ), mAxis( 1.0, 0.0, 0.0 ){};

  /**
   * Constructor from angle and axis.
   * \test PASSED
   */
  axis_angle( const value_type& aAngle, const vect< value_type, 3 >& aAxis ) BOOST_NOEXCEPT : mAngle( aAngle ) {
    using std::sqrt;
    value_type tmp = norm_2( aAxis );
    if( tmp > value_type( 0.0000001 ) ) {
      mAxis.q[0] = aAxis[0] / tmp;
      mAxis.q[1] = aAxis[1] / tmp;
      mAxis.q[2] = aAxis[2] / tmp;
    } else {
      mAxis.q[0] = value_type( 1.0 );
      mAxis.q[1] = value_type( 0.0 );
      mAxis.q[2] = value_type( 0.0 );
    };
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Copy-constructor.
   * \test PASSED
   */
  axis_angle( const self& A ) BOOST_NOEXCEPT : mAngle( A.mAngle ), mAxis( A.mAxis ){};
#else
  axis_angle( const self& A ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Constructor from a quaternion.
   * \test PASSED
   */
  explicit axis_angle( const quaternion< value_type >& Q ) BOOST_NOEXCEPT {
    using std::sqrt;
    using std::acos;
    vect< value_type, 4 > v( Q.q[0], Q.q[1], Q.q[2], Q.q[3] );
    v = unit( v );
    value_type tmp( sqrt( v[1] * v[1] + v[2] * v[2] + v[3] * v[3] ) );
    if( tmp > value_type( 0.0000001 ) ) {
      mAxis.q[0] = v[1] / tmp;
      mAxis.q[1] = v[2] / tmp;
      mAxis.q[2] = v[3] / tmp;
      if( v[0] < value_type( 0.0 ) ) {
        mAngle = value_type( 2.0 ) * acos( -v[0] );
        mAxis = -mAxis;
      } else
        mAngle = value_type( 2.0 ) * acos( v[0] );
    } else {
      mAxis.q[0] = value_type( 1.0 );
      mAxis.q[1] = value_type( 0.0 );
      mAxis.q[2] = value_type( 0.0 );
      mAngle = value_type( 0.0 );
    };
  };

  /**
   * Constructor from a rotation matrix.
   * \test PASSED
   */
  explicit axis_angle( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    using std::sin;
    using std::acos;
    value_type tmp( value_type( 0.5 ) * ( trace( R ) - value_type( 1.0 ) ) );
    if( tmp > value_type( 0.0000001 ) ) {
      mAngle = acos( tmp );
      value_type cosec_a = value_type( 0.5 ) / sin( mAngle );
      mAxis.q[0] = ( R.q[5] - R.q[7] ) * cosec_a;
      mAxis.q[1] = ( R.q[6] - R.q[2] ) * cosec_a;
      mAxis.q[2] = ( R.q[1] - R.q[3] ) * cosec_a;
    } else {
      mAxis.q[0] = value_type( 1.0 );
      mAxis.q[1] = value_type( 0.0 );
      mAxis.q[2] = value_type( 0.0 );
      mAngle = value_type( 0.0 );
    };
  };

  /**
   * Constructor from euler angles.
   * \test PASSED
   */
  explicit axis_angle( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT {
    using std::sin;
    using std::cos;
    using std::sqrt;
    using std::acos;
    value_type cpsi = cos( value_type( 0.5 ) * E.q[0] );
    value_type spsi = sin( value_type( 0.5 ) * E.q[0] );
    value_type ctheta = cos( value_type( 0.5 ) * E.q[1] );
    value_type stheta = sin( value_type( 0.5 ) * E.q[1] );
    value_type cphi = cos( value_type( 0.5 ) * E.q[2] );
    value_type sphi = sin( value_type( 0.5 ) * E.q[2] );

    value_type q[4];
    q[0] = cphi * ctheta * cpsi + sphi * stheta * spsi;
    q[1] = sphi * ctheta * cpsi - cphi * stheta * spsi;
    q[2] = cphi * stheta * cpsi + sphi * ctheta * spsi;
    q[3] = cphi * ctheta * spsi - sphi * stheta * cpsi;

    value_type tmp( sqrt( value_type( 1.0 ) - q[0] * q[0] ) );
    if( tmp > value_type( 0.0000001 ) ) {
      mAxis.q[0] = q[1] / tmp;
      mAxis.q[1] = q[2] / tmp;
      mAxis.q[2] = q[3] / tmp;
      mAngle = value_type( 2.0 ) * acos( q[0] );
    } else {
      mAxis.q[0] = value_type( 1.0 );
      mAxis.q[1] = value_type( 0.0 );
      mAxis.q[2] = value_type( 0.0 );
      mAngle = value_type( 0.0 );
    };
  };


  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides the angle, read-write.
   * \test PASSED
   */
  reference angle() BOOST_NOEXCEPT { return mAngle; };

  /**
   * Provides the axis, read-write.
   * \test PASSED
   */
  vect< value_type, 3 >& axis() BOOST_NOEXCEPT { return mAxis; };

  /**
   * Provides the angle, read-only.
   * \test PASSED
   */
  const_reference angle() const BOOST_NOEXCEPT { return mAngle; };

  /**
   * Provides the axis, read-only.
   * \test PASSED
   */
  const vect< value_type, 3 >& axis() const BOOST_NOEXCEPT { return mAxis; };

  /**
   * Provides a quaternion representation of this rotation.
   * \test PASSED
   */
  quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT {
    using std::cos;
    using std::sin;
    value_type t = norm_2( mAxis );
    if( t == value_type( 0.0 ) )
      return quaternion< value_type >( value_type( 1.0 ), value_type( 0.0 ), value_type( 0.0 ), value_type( 0.0 ) );
    t = sin( value_type( 0.5 ) * mAngle );
    return quaternion< value_type >( cos( value_type( 0.5 ) * mAngle ), mAxis.q[0] * t, mAxis.q[1] * t,
                                     mAxis.q[2] * t );
  };

  /**
   * Provides a euler angles representation of this rotation.
   * \test PASSED
   */
  euler_angles_TB< value_type > getEulerAnglesTB() const BOOST_NOEXCEPT {
    using std::sin;
    using std::cos;
    using std::asin;
    using std::atan2;
    value_type quat[4];
    value_type t = norm_2( mAxis );
    if( t == value_type( 0.0 ) ) {
      quat[0] = value_type( 1.0 );
      quat[1] = value_type( 0.0 );
      quat[2] = value_type( 0.0 );
      quat[3] = value_type( 0.0 );
    } else {
      quat[0] = sin( mAngle / value_type( 2.0 ) );
      quat[1] = mAxis.q[0] * quat[0];
      quat[2] = mAxis.q[1] * quat[0];
      quat[3] = mAxis.q[2] * quat[0];
      quat[0] = cos( mAngle / value_type( 2.0 ) );
    };
    euler_angles_TB< value_type > result;
    result.q[1] = value_type( 2.0 ) * ( quat[0] * quat[2] - quat[1] * quat[3] );
    if( ( result.q[1] != value_type( 1.0 ) ) && ( result.q[1] != value_type( -1.0 ) ) ) {
      result.q[1] = asin( result.q[1] );
      value_type cp = value_type( 1.0 ) / cos( result.q[1] );
      result.q[2] = atan2( value_type( 2.0 ) * cp * ( quat[2] * quat[3] + quat[0] * quat[1] ),
                           cp * ( value_type( 1.0 ) - value_type( 2.0 ) * ( quat[1] * quat[1] + quat[2] * quat[2] ) ) );
      result.q[0] = atan2( value_type( 2.0 ) * cp * ( quat[1] * quat[2] + quat[0] * quat[3] ),
                           cp * ( value_type( 1.0 ) - value_type( 2.0 ) * ( quat[2] * quat[2] + quat[3] * quat[3] ) ) );
    } else {
      result.q[0] = value_type( 0.0 );
      result.q[2] = atan2( result.q[1] * value_type( 2.0 ) * ( quat[1] * quat[2] - quat[0] * quat[3] ),
                           result.q[1] * value_type( 2.0 ) * ( quat[1] * quat[3] + quat[0] * quat[2] ) );
      result.q[1] *= value_type( 1.57079632679489662 );
    };
    return result;
  };

  /**
   * Provides a rotation matrix representation of this rotation.
   * \test PASSED
   */
  rot_mat_3D< value_type > getRotMat() const BOOST_NOEXCEPT {
    using std::cos;
    using std::sin;
    value_type ca( cos( mAngle ) );
    value_type one_minus_ca( value_type( 1.0 ) - ca );
    value_type t11( ca + one_minus_ca * mAxis.q[0] * mAxis.q[0] );
    value_type t22( ca + one_minus_ca * mAxis.q[1] * mAxis.q[1] );
    value_type t33( ca + one_minus_ca * mAxis.q[2] * mAxis.q[2] );
    value_type t12( one_minus_ca * mAxis.q[0] * mAxis.q[1] );
    value_type t13( one_minus_ca * mAxis.q[0] * mAxis.q[2] );
    value_type t23( one_minus_ca * mAxis.q[1] * mAxis.q[2] );
    value_type sin_a( sin( mAngle ) );
    value_type t01( sin_a * mAxis.q[0] );
    value_type t02( sin_a * mAxis.q[1] );
    value_type t03( sin_a * mAxis.q[2] );

    return rot_mat_3D< value_type >( t11, t12 - t03, t13 + t02, t12 + t03, t22, t23 - t01, t13 - t02, t23 + t01, t33 );
  };

  /**
   * Provides a 3x3 matrix representation of this rotation.
   * \test PASSED
   */
  mat< value_type, mat_structure::square > getMat() const {
    using std::cos;
    using std::sin;
    value_type ca( cos( mAngle ) );
    value_type one_minus_ca( value_type( 1.0 ) - ca );
    value_type t11( ca + one_minus_ca * mAxis.q[0] * mAxis.q[0] );
    value_type t22( ca + one_minus_ca * mAxis.q[1] * mAxis.q[1] );
    value_type t33( ca + one_minus_ca * mAxis.q[2] * mAxis.q[2] );
    value_type t12( one_minus_ca * mAxis.q[0] * mAxis.q[1] );
    value_type t13( one_minus_ca * mAxis.q[0] * mAxis.q[2] );
    value_type t23( one_minus_ca * mAxis.q[1] * mAxis.q[2] );
    value_type sin_a( sin( mAngle ) );
    value_type t01( sin_a * mAxis.q[0] );
    value_type t02( sin_a * mAxis.q[1] );
    value_type t03( sin_a * mAxis.q[2] );

    return mat< value_type, mat_structure::square >( t11, t12 - t03, t13 + t02, t12 + t03, t22, t23 - t01, t13 - t02,
                                                     t23 + t01, t33 );
  };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Assignment from an axis / angle representation.
   * \test PASSED
   */
  self& operator=( const self& A ) BOOST_NOEXCEPT {
    mAngle = A.mAngle;
    mAxis = A.mAxis;
    return *this;
  };
#else
  self& operator=( const self& A ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  self& operator=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { return ( *this = self( E ) ); };

  /**
   * Assignment from a quaternion.
   * \test PASSED
   */
  self& operator=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT { return ( *this = self( Q ) ); };

  /**
   * Assignment from a rotation matrix.
   * \test PASSED
   */
  self& operator=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT { return ( *this = self( R ) ); };

  /**
   * Multiply-and-store from a axis / angle.
   * \test PASSED
   */
  self& operator*=( const self& A ) BOOST_NOEXCEPT { return ( *this = this->getQuaternion() * A.getQuaternion() ); };

  /**
   * Multiply-and-store from a euler angles.
   * \test PASSED
   */
  self& operator*=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT {
    return ( *this = ( this->getRotMat() * E.getRotMat() ) );
  };

  /**
   * Multiply-and-store from a rotation matrix.
   * \test PASSED
   */
  self& operator*=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    return ( *this = ( this->getRotMat() * R ) );
  };

  /**
   * Multiply-and-store from a quaternion.
   * \test PASSED
   */
  self& operator*=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT {
    return ( *this = ( this->getQuaternion() * Q ) );
  };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiplication with an axis / angle representation.
   * \test PASSED
   */
  friend quaternion< value_type > operator*(const self& A1, const self& A2)BOOST_NOEXCEPT {
    return A1.getQuaternion() * A2.getQuaternion();
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const self& A, const Matrix& M ) {
    return A.getRotMat() * M;
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const Matrix& M, const self& A ) {
    return M * A.getRotMat();
  };

  /**
   * Multiplication with a column vector.
   * \test PASSED
   */
  friend vect< value_type, 3 > operator*(const self& A, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
    return A.getRotMat() * V;
  };

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality comparison operator with a axis / angle representation.
   * \test PASSED
   */
  friend bool operator==( const self& A1, const self& A2 ) BOOST_NOEXCEPT {
    return ( ( A1.mAngle == A2.mAngle ) && ( A1.mAxis == A2.mAxis ) );
  };

  /**
   * Inequality comparison operator with a axis / angle representation.
   * \test PASSED
   */
  friend bool operator!=( const self& A1, const self& A2 ) BOOST_NOEXCEPT {
    return ( ( A1.mAngle != A2.mAngle ) || ( A1.mAxis != A2.mAxis ) );
  };

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Produces a transpose axis/angle which is the inverse rotation.
   * \test PASSED
   */
  friend self transpose( const self& A ) BOOST_NOEXCEPT { return self( -A.mAngle, A.mAxis ); };

  /**
   * Produces a cofactor matrix which is the same as the rotation itself.
   * \test PASSED
   */
  friend self cofactor( const self& A ) BOOST_NOEXCEPT { return A; };

  /**
   * Invert the rotation.
   * \test PASSED
   */
  friend self invert( const self& A ) BOOST_NOEXCEPT { return self( -A.mAngle, A.mAxis ); };

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace( const self& A ) BOOST_NOEXCEPT {
    using std::cos;
    return value_type( 2.0 ) * cos( A.mAngle ) + value_type( 1.0 );
  };

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant( const self& ) BOOST_NOEXCEPT { return value_type( 1.0 ); };

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::symmetric > getSymPart() const {
    using std::cos;
    value_type ca( cos( mAngle ) );
    value_type one_minus_ca( value_type( 1.0 ) - ca );
    value_type t11( ca + one_minus_ca * mAxis.q[0] * mAxis.q[0] );
    value_type t22( ca + one_minus_ca * mAxis.q[1] * mAxis.q[1] );
    value_type t33( ca + one_minus_ca * mAxis.q[2] * mAxis.q[2] );
    value_type t12( one_minus_ca * mAxis.q[0] * mAxis.q[1] );
    value_type t13( one_minus_ca * mAxis.q[0] * mAxis.q[2] );
    value_type t23( one_minus_ca * mAxis.q[1] * mAxis.q[2] );
    return mat< value_type, mat_structure::symmetric >( t11, t12, t13, t22, t23, t33 );
  };

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::skew_symmetric > getSkewSymPart() const {
    using std::sin;
    value_type sin_a( sin( mAngle ) );
    value_type t01( sin_a * mAxis.q[0] );
    value_type t02( sin_a * mAxis.q[1] );
    value_type t03( sin_a * mAxis.q[2] );
    return mat< value_type, mat_structure::skew_symmetric >( -t03, t02, -t01 );
  };


  /// Loading a axis_angle value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, axis_angle< T >& >& R ) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_angle", R.second.mAngle )
           & RK_SERIAL_LOAD_WITH_ALIAS( R.first + "_axis", R.second.mAxis );
  };

  /// Loading a axis_angle value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, axis_angle< T >& R ) {
    return in & RK_SERIAL_LOAD_WITH_NAME( R );
  };

  /// Saving a axis_angle value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const axis_angle< T >& >& R ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_angle", R.second.mAngle )
           & RK_SERIAL_SAVE_WITH_ALIAS( R.first + "_axis", R.second.mAxis );
  };

  /// Saving a axis_angle value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const axis_angle< T >& R ) {
    return out & RK_SERIAL_SAVE_WITH_NAME( R );
  };
};

namespace rtti {

template < typename T >
struct get_type_id< axis_angle< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x0000001C );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::axis_angle" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::axis_angle"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const axis_angle< T >& save_type;
  typedef axis_angle< T >& load_type;
};
};


/**
 * Prints a axis / angle to a standard output stream (<<) as "(Angle = value; Axis = (a1; a2; a3))".
 * \test PASSED
 */
template < class T >
std::ostream& operator<<( std::ostream& out_stream, const axis_angle< T >& A ) {
  return out_stream << "(Angle = " << A.angle() << "; Axis = " << A.axis() << ")";
};


/**
 * This class is a transformation matrix 4 by 4, i.e. to rotate and translate a 3D vector.
 * \test All tests for this class have been passed!
 */
template < class T >
class trans_mat_3D {
public:
  typedef trans_mat_3D< T > self;
  typedef void allocator_type;

  typedef T value_type;
  typedef void container_type;

  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;

  typedef void col_iterator;
  typedef void const_col_iterator;
  typedef void row_iterator;
  typedef void const_row_iterator;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 4 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 4 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_alignment::column_major );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::orthogonal );

  typedef vect< value_type, 3 > translation_type;

private:
  value_type q[16];

  trans_mat_3D( const_reference a11, const_reference a12, const_reference a13, const_reference a14, const_reference a21,
                const_reference a22, const_reference a23, const_reference a24, const_reference a31, const_reference a32,
                const_reference a33, const_reference a34 ) BOOST_NOEXCEPT {
    q[0] = a11;
    q[1] = a21;
    q[2] = a31;
    q[3] = value_type( 0.0 );
    q[4] = a12;
    q[5] = a22;
    q[6] = a32;
    q[7] = value_type( 0.0 );
    q[8] = a13;
    q[9] = a23;
    q[10] = a33;
    q[11] = value_type( 0.0 );
    q[12] = a14;
    q[13] = a24;
    q[14] = a34;
    q[15] = value_type( 1.0 );
  };

public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default Constructor.
   * \test PASSED
   */
  trans_mat_3D() BOOST_NOEXCEPT {
    std::fill( q, q + 16, value_type( 0.0 ) );
    q[0] = value_type( 1.0 );
    q[5] = value_type( 1.0 );
    q[10] = value_type( 1.0 );
    q[15] = value_type( 1.0 );
  };

  /**
   * Constructor from a 4x4 array (16 values).
   * \test PASSED
   */
  explicit trans_mat_3D( const_pointer M ) BOOST_NOEXCEPT {
    vect< value_type, 3 > v1 = unit( vect< value_type, 3 >( M ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    q[3] = value_type( 0.0 );
    vect< value_type, 3 > v2( M[4], M[5], M[6] );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[4] = v2[0];
    q[5] = v2[1];
    q[6] = v2[2];
    q[7] = value_type( 0.0 );
    v2 = v1 % v2;
    q[8] = v2[0];
    q[9] = v2[1];
    q[10] = v2[2];
    q[11] = value_type( 0.0 );
    q[12] = M[12];
    q[13] = M[13];
    q[14] = M[14];
    q[15] = value_type( 1.0 );
    return;
  };

  /**
   * Constructor from a regular matrix.
   * \test PASSED
   */
  template < typename Matrix >
  explicit trans_mat_3D(
    const Matrix& M,
    typename boost::enable_if_c< is_readable_matrix< Matrix >::value && !boost::is_same< Matrix, self >::value
                                 && !boost::is_same< Matrix, rot_mat_3D< value_type > >::value,
                                 void* >::type dummy = nullptr ) {
    if( ( M.get_row_count() != 4 ) || ( M.get_col_count() != 4 ) )
      throw std::range_error( "Matrix for creating the 3D transformation matrix is not of correct dimensions!" );
    vect< value_type, 3 > v1 = unit( vect< value_type, 3 >( M( 0, 0 ), M( 1, 0 ), M( 2, 0 ) ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    q[3] = value_type( 0.0 );
    vect< value_type, 3 > v2( M( 0, 1 ), M( 1, 1 ), M( 2, 1 ) );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[4] = v2[0];
    q[5] = v2[1];
    q[6] = v2[2];
    q[7] = value_type( 0.0 );
    v2 = v1 % v2;
    q[8] = v2[0];
    q[9] = v2[1];
    q[10] = v2[2];
    q[11] = value_type( 0.0 );
    q[12] = M( 0, 3 );
    q[13] = M( 1, 3 );
    q[14] = M( 2, 3 );
    q[15] = value_type( 1.0 );
    return;
  };

  /**
   * Constructor from a rotation matrix and an optional translation vector V.
   * \test PASSED
   */
  explicit trans_mat_3D( const rot_mat_3D< value_type >& R,
                         const translation_type& V = translation_type( value_type( 0.0 ), value_type( 0.0 ),
                                                                       value_type( 0.0 ) ) ) BOOST_NOEXCEPT {
    setRotMat( R );
    q[3] = value_type( 0.0 );
    q[7] = value_type( 0.0 );
    q[11] = value_type( 0.0 );
    setTranslation( V );
    q[15] = value_type( 1.0 );
  };

  /**
   * Constructor from a quaternion representation and an optional translation vector V.
   * \test PASSED
   */
  explicit trans_mat_3D( const quaternion< value_type >& Q,
                         const translation_type& V = translation_type( value_type( 0.0 ), value_type( 0.0 ),
                                                                       value_type( 0.0 ) ) ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( Q.getRotMat() );
    setRotMat( R );
    q[3] = value_type( 0.0 );
    q[7] = value_type( 0.0 );
    q[11] = value_type( 0.0 );
    setTranslation( V );
    q[15] = value_type( 1.0 );
  };

  /**
   * Constructor from a euler angles TB representation and an optional translation vector V.
   * \test PASSED
   */
  explicit trans_mat_3D( const euler_angles_TB< value_type >& E,
                         const translation_type& V = translation_type( value_type( 0.0 ), value_type( 0.0 ),
                                                                       value_type( 0.0 ) ) ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( E.getRotMat() );
    setRotMat( R );
    q[3] = value_type( 0.0 );
    q[7] = value_type( 0.0 );
    q[11] = value_type( 0.0 );
    setTranslation( V );
    q[15] = value_type( 1.0 );
  };

  /**
   * Constructor from an axis / angle representation and an optional translation vector V.
   * \test PASSED
   */
  explicit trans_mat_3D( const axis_angle< value_type >& A,
                         const translation_type& V = translation_type( value_type( 0.0 ), value_type( 0.0 ),
                                                                       value_type( 0.0 ) ) ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( A.getRotMat() );
    setRotMat( R );
    q[3] = value_type( 0.0 );
    q[7] = value_type( 0.0 );
    q[11] = value_type( 0.0 );
    setTranslation( V );
    q[15] = value_type( 1.0 );
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  // Copy-constructor. Default is good. \test PASSED
  trans_mat_3D( const self& T ) BOOST_NOEXCEPT { std::copy( T.q, T.q + 16, q ); };
#else
  trans_mat_3D( const self& ) BOOST_NOEXCEPT = default;
#endif

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides a copy of the transformation matrix as an ordinary 4x4 matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::square > getMat() const {
    return mat< value_type, mat_structure::square >( q[0], q[4], q[8], q[12], q[1], q[5], q[9], q[13], q[2], q[6],
                                                     q[10], q[14], value_type( 0.0 ), value_type( 0.0 ),
                                                     value_type( 0.0 ), value_type( 1.0 ) );
  };

  /**
   * Provides the rotation part of the transformation as a rotation matrix.
   * \test PASSED
   */
  rot_mat_3D< value_type > getRotMat() const BOOST_NOEXCEPT {
    return rot_mat_3D< value_type >( q[0], q[4], q[8], q[1], q[5], q[9], q[2], q[6], q[10] );
  };

  /**
   * Sets the rotation part of the transformation from a rotation matrix.
   * \test PASSED
   */
  void setRotMat( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
  };

  /**
   * Returns the quaternion of the rotation matrix.
   * \test PASSED
   */
  quaternion< value_type > getQuaternion() const BOOST_NOEXCEPT { return quaternion< value_type >( getRotMat() ); };

  /**
   * Sets the quaternion of the rotation matrix.
   * \test PASSED
   */
  void setQuaternion( const quaternion< value_type >& Q ) BOOST_NOEXCEPT { setRotMat( Q.getRotMat() ); };

  /**
   * Returns the euler angles TB of the rotation matrix.
   * \test PASSED
   */
  euler_angles_TB< value_type > getEulerAnglesTB() const BOOST_NOEXCEPT {
    return euler_angles_TB< value_type >( getRotMat() );
  };

  /**
   * Sets the euler angles TB of the rotation matrix.
   * \test PASSED
   */
  void setEulerAnglesTB( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { setRotMat( E.getRotMat() ); };

  /**
   * Returns the axis / angle of the rotation matrix.
   * \test PASSED
   */
  axis_angle< value_type > getAxisAngle() const BOOST_NOEXCEPT { return axis_angle< value_type >( getRotMat() ); };

  /**
   * Sets the axis / angle of the rotation matrix.
   * \test PASSED
   */
  void setAxisAngle( const axis_angle< value_type >& A ) BOOST_NOEXCEPT { setRotMat( A.getRotMat() ); };

  /**
   * Provides the translation part of the transformation matrix as a vector.
   * \test PASSED
   */
  translation_type getTranslation() const BOOST_NOEXCEPT { return translation_type( q[12], q[13], q[14] ); };

  /**
   * Sets the translation part of the transformation matrix to a vector.
   * \test PASSED
   */
  void setTranslation( const translation_type& Translation ) BOOST_NOEXCEPT {
    q[12] = Translation[0];
    q[13] = Translation[1];
    q[14] = Translation[2];
  };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT {
    assert( i < 16 );
    return q[i];
  };

  /**
   * Array double-indexing operator, ith row and jth column, accessor for read only.
   * \test PASSED
   */
  const_reference operator()( size_type i, size_type j ) const BOOST_NOEXCEPT {
    assert( ( i < 4 ) || ( j < 4 ) );
    return q[j * 4 + i];
  };

  size_type get_row_count() const BOOST_NOEXCEPT { return 4; };
  size_type get_col_count() const BOOST_NOEXCEPT { return 4; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Assignment operator with transformation matrix.
   * \test PASSED
   */
  self& operator=( const self& M ) BOOST_NOEXCEPT {
    std::copy( M.q, M.q + 16, q );
    return *this;
  };
#else
  self& operator=( const self& ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Assignment operator with regular matrix.
   * \test PASSED
   */
  template < typename Matrix >
  typename boost::enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                               boost::mpl::not_< boost::is_same< Matrix, self > >,
                                               boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
                             self& >::type
    operator=( const Matrix& M ) {
    if( ( M.get_row_count() != 4 ) || ( M.get_col_count() != 4 ) )
      throw std::range_error( "Matrix for creating the 3D transformation matrix is not of correct dimensions!" );
    vect< value_type, 3 > v1 = unit( vect< value_type, 3 >( M( 0, 0 ), M( 1, 0 ), M( 2, 0 ) ) );
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    q[3] = value_type( 0.0 );
    vect< value_type, 3 > v2( M( 0, 1 ), M( 1, 1 ), M( 2, 1 ) );
    v2 = unit( v2 - ( v2 * v1 ) * v1 );
    q[4] = v2[0];
    q[5] = v2[1];
    q[6] = v2[2];
    q[7] = value_type( 0.0 );
    v2 = v1 % v2;
    q[8] = v2[0];
    q[9] = v2[1];
    q[10] = v2[2];
    q[11] = value_type( 0.0 );
    q[12] = M( 0, 3 );
    q[13] = M( 1, 3 );
    q[14] = M( 2, 3 );
    q[15] = value_type( 1.0 );
    return *this;
  };

  /**
   * Assignment operator with rotation matrix.
   * \test PASSED
   */
  self& operator=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[3] = value_type( 0.0 );
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[7] = value_type( 0.0 );
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
    q[11] = value_type( 0.0 );
    q[12] = value_type( 0.0 );
    q[13] = value_type( 0.0 );
    q[14] = value_type( 0.0 );
    q[15] = value_type( 1.0 );
    return *this;
  };

  /**
   * Assignment operator with a quaternion representation.
   * \test PASSED
   */
  self& operator=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( Q.getRotMat() );
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[3] = value_type( 0.0 );
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[7] = value_type( 0.0 );
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
    q[11] = value_type( 0.0 );
    q[12] = value_type( 0.0 );
    q[13] = value_type( 0.0 );
    q[14] = value_type( 0.0 );
    q[15] = value_type( 1.0 );
    return *this;
  };

  /**
   * Assignment operator with euler angles TB representation.
   * \test PASSED
   */
  self& operator=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( E.getRotMat() );
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[3] = value_type( 0.0 );
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[7] = value_type( 0.0 );
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
    q[11] = value_type( 0.0 );
    q[12] = value_type( 0.0 );
    q[13] = value_type( 0.0 );
    q[14] = value_type( 0.0 );
    q[15] = value_type( 1.0 );
    return *this;
  };

  /**
   * Assignment operator with an axis / angle representation.
   * \test PASSED
   */
  self& operator=( const axis_angle< value_type >& A ) BOOST_NOEXCEPT {
    rot_mat_3D< value_type > R( A.getRotMat() );
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[3] = value_type( 0.0 );
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[7] = value_type( 0.0 );
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
    q[11] = value_type( 0.0 );
    q[12] = value_type( 0.0 );
    q[13] = value_type( 0.0 );
    q[14] = value_type( 0.0 );
    q[15] = value_type( 1.0 );
    return *this;
  };

  /**
   * Multiply-and-store with a transformation matrix.
   * \test PASSED
   */
  self& operator*=( const self& M ) BOOST_NOEXCEPT {
    ( *this ) = self(
      q[0] * M.q[0] + q[4] * M.q[1] + q[8] * M.q[2], q[0] * M.q[4] + q[4] * M.q[5] + q[8] * M.q[6],
      q[0] * M.q[8] + q[4] * M.q[9] + q[8] * M.q[10], q[0] * M.q[12] + q[4] * M.q[13] + q[8] * M.q[14] + q[12],
      q[1] * M.q[0] + q[5] * M.q[1] + q[9] * M.q[2], q[1] * M.q[4] + q[5] * M.q[5] + q[9] * M.q[6],
      q[1] * M.q[8] + q[5] * M.q[9] + q[9] * M.q[10], q[1] * M.q[12] + q[5] * M.q[13] + q[9] * M.q[14] + q[13],
      q[2] * M.q[0] + q[6] * M.q[1] + q[10] * M.q[2], q[2] * M.q[4] + q[6] * M.q[5] + q[10] * M.q[6],
      q[2] * M.q[8] + q[6] * M.q[9] + q[10] * M.q[10], q[2] * M.q[12] + q[6] * M.q[13] + q[10] * M.q[14] + q[14] );
    return *this;
  };

  /**
   * Multiply-and-store with a rotation matrix.
   * \test PASSED
   */
  self& operator*=( const rot_mat_3D< value_type >& R ) BOOST_NOEXCEPT {
    ( *this )
      = self( q[0] * R.q[0] + q[4] * R.q[1] + q[8] * R.q[2], q[0] * R.q[3] + q[4] * R.q[4] + q[8] * R.q[5],
              q[0] * R.q[6] + q[4] * R.q[7] + q[8] * R.q[8], q[12], q[1] * R.q[0] + q[5] * R.q[1] + q[9] * R.q[2],
              q[1] * R.q[3] + q[5] * R.q[4] + q[9] * R.q[5], q[1] * R.q[6] + q[5] * R.q[7] + q[9] * R.q[8], q[13],
              q[2] * R.q[0] + q[6] * R.q[1] + q[10] * R.q[2], q[2] * R.q[3] + q[6] * R.q[4] + q[10] * R.q[5],
              q[2] * R.q[6] + q[6] * R.q[7] + q[10] * R.q[8], q[14] );
    return *this;
  };

  /**
   * Multiply-and-store with a quaternion representation.
   * \test PASSED
   */
  self& operator*=( const quaternion< value_type >& Q ) BOOST_NOEXCEPT { return ( *this *= Q.getRotMat() ); };

  /**
   * Multiply-and-store with a euler angles TB representation.
   * \test PASSED
   */
  self& operator*=( const euler_angles_TB< value_type >& E ) BOOST_NOEXCEPT { return ( *this *= E.getRotMat() ); };

  /**
   * Multiply-and-store with an axis / angle representation.
   * \test PASSED
   */
  self& operator*=( const axis_angle< value_type >& A ) BOOST_NOEXCEPT { return ( *this *= A.getRotMat() ); };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiplication with a transformation matrix.
   * \test PASSED
   */
  friend self operator*(const self& M1, const self& M2)BOOST_NOEXCEPT {
    return self( M1.q[0] * M2.q[0] + M1.q[4] * M2.q[1] + M1.q[8] * M2.q[2],
                 M1.q[0] * M2.q[4] + M1.q[4] * M2.q[5] + M1.q[8] * M2.q[6],
                 M1.q[0] * M2.q[8] + M1.q[4] * M2.q[9] + M1.q[8] * M2.q[10],
                 M1.q[0] * M2.q[12] + M1.q[4] * M2.q[13] + M1.q[8] * M2.q[14] + M1.q[12],
                 M1.q[1] * M2.q[0] + M1.q[5] * M2.q[1] + M1.q[9] * M2.q[2],
                 M1.q[1] * M2.q[4] + M1.q[5] * M2.q[5] + M1.q[9] * M2.q[6],
                 M1.q[1] * M2.q[8] + M1.q[5] * M2.q[9] + M1.q[9] * M2.q[10],
                 M1.q[1] * M2.q[12] + M1.q[5] * M2.q[13] + M1.q[9] * M2.q[14] + M1.q[13],
                 M1.q[2] * M2.q[0] + M1.q[6] * M2.q[1] + M1.q[10] * M2.q[2],
                 M1.q[2] * M2.q[4] + M1.q[6] * M2.q[5] + M1.q[10] * M2.q[6],
                 M1.q[2] * M2.q[8] + M1.q[6] * M2.q[9] + M1.q[10] * M2.q[10],
                 M1.q[2] * M2.q[12] + M1.q[6] * M2.q[13] + M1.q[10] * M2.q[14] + M1.q[14] );
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const self& M1, const Matrix& M2 ) {
    return Matrix( M1.getMat() * M2 );
  };

  /**
   * Multiply by a matrix.
   * \test PASSED
   */
  template < typename Matrix >
  friend typename boost::
    enable_if< boost::mpl::and_< is_readable_matrix< Matrix >,
                                 boost::mpl::not_< boost::is_same< Matrix, trans_mat_3D< value_type > > >,
                                 boost::mpl::not_< boost::is_same< Matrix, rot_mat_3D< value_type > > > >,
               Matrix >::type
    operator*( const Matrix& M1, const self& M2 ) {
    return Matrix( M1 * M2.getMat() );
  };

  /**
   * Multiplication with a rotation matrix.
   * \test PASSED
   */
  friend self operator*(const self& M, const rot_mat_3D< T >& R)BOOST_NOEXCEPT {
    return self( M.q[0] * R[0] + M.q[4] * R[1] + M.q[8] * R[2], M.q[0] * R[3] + M.q[4] * R[4] + M.q[8] * R[5],
                 M.q[0] * R[6] + M.q[4] * R[7] + M.q[8] * R[8], M.q[12], M.q[1] * R[0] + M.q[5] * R[1] + M.q[9] * R[2],
                 M.q[1] * R[3] + M.q[5] * R[4] + M.q[9] * R[5], M.q[1] * R[6] + M.q[5] * R[7] + M.q[9] * R[8], M.q[13],
                 M.q[2] * R[0] + M.q[6] * R[1] + M.q[10] * R[2], M.q[2] * R[3] + M.q[6] * R[4] + M.q[10] * R[5],
                 M.q[2] * R[6] + M.q[6] * R[7] + M.q[10] * R[8], M.q[14] );
  };

  /**
   * Multiplication with a 3D column vector.
   * \test PASSED
   */
  friend vect< value_type, 3 > operator*(const self& M, const vect< value_type, 3 >& V)BOOST_NOEXCEPT {
    return vect< value_type, 3 >( M.q[0] * V[0] + M.q[4] * V[1] + M.q[8] * V[2] + M.q[12],
                                  M.q[1] * V[0] + M.q[5] * V[1] + M.q[9] * V[2] + M.q[13],
                                  M.q[2] * V[0] + M.q[6] * V[1] + M.q[10] * V[2] + M.q[14] );
  };

  /**
   * Multiplication with a 4D column vector.
   * \test PASSED
   */
  friend vect< value_type, 4 > operator*(const self& M, const vect< value_type, 4 >& V)BOOST_NOEXCEPT {
    return vect< value_type, 4 >( M.q[0] * V[0] + M.q[4] * V[1] + M.q[8] * V[2] + M.q[12] * V[3],
                                  M.q[1] * V[0] + M.q[5] * V[1] + M.q[9] * V[2] + M.q[13] * V[3],
                                  M.q[2] * V[0] + M.q[6] * V[1] + M.q[10] * V[2] + M.q[14] * V[3], V[3] );
  };

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality operator with a transformation matrix.
   * \test PASSED
   */
  friend bool operator==( const self& M1, const self& M2 ) BOOST_NOEXCEPT {
    return ( ( M1.q[0] == M2.q[0] ) && ( M1.q[1] == M2.q[1] ) && ( M1.q[2] == M2.q[2] ) && ( M1.q[4] == M2.q[4] )
             && ( M1.q[5] == M2.q[5] ) && ( M1.q[6] == M2.q[6] ) && ( M1.q[12] == M2.q[12] ) && ( M1.q[13] == M2.q[13] )
             && ( M1.q[14] == M2.q[14] ) );
  };

  /**
   * Inequality operator with a transformation matrix.
   * \test PASSED
   */
  friend bool operator!=( const self& M1, const self& M2 ) BOOST_NOEXCEPT {
    return ( ( M1.q[0] != M2.q[0] ) || ( M1.q[1] != M2.q[1] ) || ( M1.q[2] != M2.q[2] ) || ( M1.q[4] != M2.q[4] )
             || ( M1.q[5] != M2.q[5] ) || ( M1.q[6] != M2.q[6] ) || ( M1.q[12] != M2.q[12] ) || ( M1.q[13] != M2.q[13] )
             || ( M1.q[14] != M2.q[14] ) );
  };

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /**
   * Rotate-only a 3D column vector.
   * \test PASSED
   */
  vect< value_type, 3 > rotate( const vect< value_type, 3 >& V ) const BOOST_NOEXCEPT {
    return vect< value_type, 3 >( q[0] * V[0] + q[4] * V[1] + q[8] * V[2], q[1] * V[0] + q[5] * V[1] + q[9] * V[2],
                                  q[2] * V[0] + q[6] * V[1] + q[10] * V[2] );
  };

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Creates the transpose matrix.
   * \note the matrix is no longer a transformation matrix.
   * \test PASSED
   */
  friend mat< value_type, mat_structure::square > transpose( const self& M ) BOOST_NOEXCEPT {
    return mat< value_type, mat_structure::square >( M.q[0], M.q[1], M.q[2], 0.0, M.q[4], M.q[5], M.q[6], 0.0, M.q[8],
                                                     M.q[9], M.q[10], 0.0, M.q[12], M.q[13], M.q[14], 1.0 );
  };

  /**
   * Creates the transpose matrix.
   * \note the matrix is no longer a transformation matrix.
   * \test PASSED
   */
  friend mat< value_type, mat_structure::square > transpose_move( const self& M ) BOOST_NOEXCEPT {
    return mat< value_type, mat_structure::square >( M.q[0], M.q[1], M.q[2], 0.0, M.q[4], M.q[5], M.q[6], 0.0, M.q[8],
                                                     M.q[9], M.q[10], 0.0, M.q[12], M.q[13], M.q[14], 1.0 );
  };

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace( const self& M ) BOOST_NOEXCEPT { return M.q[0] + M.q[5] + M.q[10] + value_type( 1.0 ); };

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant( const self& ) BOOST_NOEXCEPT { return value_type( 1.0 ); };

  /**
   * Invert the transformation.
   * \test PASSED
   */
  friend self invert( const self& M ) BOOST_NOEXCEPT {
    return self( M.q[0], M.q[1], M.q[2], -M.q[0] * M.q[12] - M.q[1] * M.q[13] - M.q[2] * M.q[14], M.q[4], M.q[5],
                 M.q[6], -M.q[4] * M.q[12] - M.q[5] * M.q[13] - M.q[6] * M.q[14], M.q[8], M.q[9], M.q[10],
                 -M.q[8] * M.q[12] - M.q[9] * M.q[13] - M.q[10] * M.q[14] );
  };

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::symmetric > getSymPart() const {
    return mat< value_type, mat_structure::symmetric >( getMat() );
  };

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::skew_symmetric > getSkewSymPart() const {
    return mat< value_type, mat_structure::skew_symmetric >( getMat() );
  };


  /// Loading a trans_mat_3D value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, trans_mat_3D< T >& >& M ) {
    in& RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r11", M.second.q[0] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r21", M.second.q[1] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r31", M.second.q[2] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r12", M.second.q[4] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r22", M.second.q[5] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r32", M.second.q[6] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r13", M.second.q[8] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r23", M.second.q[9] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_r33", M.second.q[10] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_t_x", M.second.q[12] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_t_y", M.second.q[13] )
      & RK_SERIAL_LOAD_WITH_ALIAS( M.first + "_t_z", M.second.q[14] );
    M.second.q[3] = 0.0;
    M.second.q[7] = 0.0;
    M.second.q[11] = 0.0;
    M.second.q[15] = 1.0;
    return in;
  };

  /// Loading a trans_mat_3D value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, trans_mat_3D< T >& M ) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS( "T", M );
  };

  /// Saving a trans_mat_3D value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const trans_mat_3D< T >& >& M ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r11", M.second.q[0] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r21", M.second.q[1] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r31", M.second.q[2] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r12", M.second.q[4] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r22", M.second.q[5] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r32", M.second.q[6] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r13", M.second.q[8] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r23", M.second.q[9] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_r33", M.second.q[10] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_t_x", M.second.q[12] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_t_y", M.second.q[13] )
           & RK_SERIAL_SAVE_WITH_ALIAS( M.first + "_t_z", M.second.q[14] );
  };

  /// Saving a trans_mat_3D value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const trans_mat_3D< T >& M ) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS( "T", M );
  };
};

namespace rtti {

template < typename T >
struct get_type_id< trans_mat_3D< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000019 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::trans_mat_3D" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::trans_mat_3D"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const trans_mat_3D< T >& save_type;
  typedef trans_mat_3D< T >& load_type;
};
};

/**
 * Prints a 3D transformation matrix to a standard output stream (<<)
 * as "(quaternion = (q0; q1; q2; q3); translation = (tx; ty; tz))".
 * \test PASSED
 */
template < class T >
std::ostream& operator<<( std::ostream& out_stream, const trans_mat_3D< T >& M ) {
  out_stream << "(quaternion = " << M.getQuaternion() << "; translation = " << M.getTranslation() << ")";
  return out_stream;
};


template < typename T >
struct is_readable_matrix< trans_mat_3D< T > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< trans_mat_3D< T > > type;
};


/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const rot_mat_3D< T >& R, const quaternion< T >& Q)BOOST_NOEXCEPT {
  return R * Q.getRotMat();
};

/**
 * Multiplication by a rotation matrix.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const quaternion< T >& Q, const rot_mat_3D< T >& R)BOOST_NOEXCEPT {
  return Q.getRotMat() * R;
};

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const rot_mat_3D< T >& R, const euler_angles_TB< T >& E)BOOST_NOEXCEPT {
  return R * E.getRotMat();
};

/**
 * Multiply by a rotation matrix representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const euler_angles_TB< T >& E, const rot_mat_3D< T >& R)BOOST_NOEXCEPT {
  return E.getRotMat() * R;
};

/**
 * Multiplication by a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
quaternion< T > operator*(const quaternion< T >& Q, const euler_angles_TB< T >& E)BOOST_NOEXCEPT {
  return Q * E.getQuaternion();
};

/**
 * Multiply by a quaternion representation.
 * \test PASSED
 */
template < typename T >
quaternion< T > operator*(const euler_angles_TB< T >& E, const quaternion< T >& Q)BOOST_NOEXCEPT {
  return E.getQuaternion() * Q;
};

/**
 * Multiplication with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const rot_mat_3D< T >& R, const axis_angle< T >& A)BOOST_NOEXCEPT {
  return R * A.getRotMat();
};

/**
 * Multiplication with a rotation matrix.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const axis_angle< T >& A, const rot_mat_3D< T >& R)BOOST_NOEXCEPT {
  return A.getRotMat() * R;
};

/**
 * Multiplication by an axis / angle representation.
 * \test PASSED
 */
template < typename T >
quaternion< T > operator*(const quaternion< T >& Q, const axis_angle< T >& A)BOOST_NOEXCEPT {
  return Q * A.getQuaternion();
};

/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
template < typename T >
quaternion< T > operator*(const axis_angle< T >& A, const quaternion< T >& Q)BOOST_NOEXCEPT {
  return A.getQuaternion() * Q;
};

/**
 * Multiply by an axis / angle representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const euler_angles_TB< T >& E, const axis_angle< T >& A)BOOST_NOEXCEPT {
  return E.getRotMat() * A.getRotMat();
};

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
rot_mat_3D< T > operator*(const axis_angle< T >& A, const euler_angles_TB< T >& E)BOOST_NOEXCEPT {
  return A.getRotMat() * E.getRotMat();
};

/**
 * Multiplication with a transformation matrix.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const rot_mat_3D< T >& R, const trans_mat_3D< T >& M)BOOST_NOEXCEPT {
  return trans_mat_3D< T >( R ) * M;
};

/**
 * Multiplication by a transformation matrix.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const quaternion< T >& Q, const trans_mat_3D< T >& M)BOOST_NOEXCEPT {
  return trans_mat_3D< T >( Q.getRotMat() ) * M;
};

/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const trans_mat_3D< T >& M, const quaternion< T >& Q)BOOST_NOEXCEPT {
  return M * Q.getRotMat();
};

/**
 * Multiply by a transformation matrix.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const euler_angles_TB< T >& E, const trans_mat_3D< T >& M)BOOST_NOEXCEPT {
  return trans_mat_3D< T >( E ) * M;
};

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const trans_mat_3D< T >& M, const euler_angles_TB< T >& E)BOOST_NOEXCEPT {
  return M * E.getRotMat();
};

/**
 * Multiplication with a transformation matrix.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const axis_angle< T >& A, const trans_mat_3D< T >& M)BOOST_NOEXCEPT {
  return trans_mat_3D< T >( A ) * M;
};

/**
 * Multiplication with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
trans_mat_3D< T > operator*(const trans_mat_3D< T >& M, const axis_angle< T >& A)BOOST_NOEXCEPT {
  return M * A.getRotMat();
};


/**
 * Equality operator for a quaternion representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const rot_mat_3D< T >& R, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q.getRotMat() == R;
};

/**
 * Inequality operator for a quaternion representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const rot_mat_3D< T >& R, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q.getRotMat() != R;
};

/**
 * Equality operator with a rotation matrix.
 * \test PASSED
 */
template < typename T >
bool operator==( const quaternion< T >& Q, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return Q.getRotMat() == R;
};

/**
 * Inequality operator with a rotation matrix.
 * \test PASSED
 */
template < typename T >
bool operator!=( const quaternion< T >& Q, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return Q.getRotMat() != R;
};

/**
 * Equality operator with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const quaternion< T >& Q, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return Q == E.getQuaternion();
};

/**
 * Inequality operator with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const quaternion< T >& Q, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return Q != E.getQuaternion();
};

/**
 * Equality operator with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const euler_angles_TB< T >& E, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q == E.getQuaternion();
};

/**
 * Inequality operator with a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const euler_angles_TB< T >& E, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q != E.getQuaternion();
};

/**
 * Equality operator for a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const rot_mat_3D< T >& R, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return R == E.getRotMat();
};

/**
 * Inequality operator for a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const rot_mat_3D< T >& R, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return R != E.getRotMat();
};

/**
 * Equality operator for a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const euler_angles_TB< T >& E, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return R == E.getRotMat();
};

/**
 * Inequality operator for a euler angles TB representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const euler_angles_TB< T >& E, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return R != E.getRotMat();
};

/**
 * Equality operator for an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const rot_mat_3D< T >& R, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return R == A.getRotMat();
};

/**
 * Inequality operator for an axis /angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const rot_mat_3D< T >& R, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return R != A.getRotMat();
};

/**
 * Equality operator for an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const axis_angle< T >& A, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return R == A.getRotMat();
};

/**
 * Inequality operator for an axis /angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const axis_angle< T >& A, const rot_mat_3D< T >& R ) BOOST_NOEXCEPT {
  return R != A.getRotMat();
};

/**
 * Equality operator with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const quaternion< T >& Q, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return Q == A.getQuaternion();
};

/**
 * Inequality operator with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const quaternion< T >& Q, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return Q != A.getQuaternion();
};

/**
 * Equality operator with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const axis_angle< T >& A, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q == A.getQuaternion();
};

/**
 * Inequality operator with an axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const axis_angle< T >& A, const quaternion< T >& Q ) BOOST_NOEXCEPT {
  return Q != A.getQuaternion();
};

/**
 * Equality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const euler_angles_TB< T >& E, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return E == A.getEulerAnglesTB();
};

/**
 * Inequality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const euler_angles_TB< T >& E, const axis_angle< T >& A ) BOOST_NOEXCEPT {
  return E != A.getEulerAnglesTB();
};

/**
 * Equality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator==( const axis_angle< T >& A, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return E == A.getEulerAnglesTB();
};

/**
 * Inequality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template < typename T >
bool operator!=( const axis_angle< T >& A, const euler_angles_TB< T >& E ) BOOST_NOEXCEPT {
  return E != A.getEulerAnglesTB();
};


#if 0
// #ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template class rot_mat_3D<double>;
extern template class quaternion<double>;
extern template class euler_angles_TB<double>;
extern template class axis_angle<double>;
extern template class trans_mat_3D<double>;

extern template std::ostream& operator << <double>(std::ostream& out_stream, const rot_mat_3D<double>& R);
extern template std::ostream& operator << <double>(std::ostream& out_stream, const quaternion<double>& Q);
extern template std::ostream& operator << <double>(std::ostream& out_stream, const euler_angles_TB<double>& E);
extern template std::ostream& operator << <double>(std::ostream& out_stream, const axis_angle<double>& A);
extern template std::ostream& operator << <double>(std::ostream& out_stream, const trans_mat_3D<double>& M);


extern template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const quaternion<double>& Q);
extern template rot_mat_3D<double> operator *(const quaternion<double>& Q, const rot_mat_3D<double>& R);
extern template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
extern template rot_mat_3D<double> operator *(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
extern template quaternion<double> operator *(const quaternion<double>& Q, const euler_angles_TB<double>& E);
extern template quaternion<double> operator *(const euler_angles_TB<double>& E, const quaternion<double>& Q);
extern template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const axis_angle<double>& A);
extern template rot_mat_3D<double> operator *(const axis_angle<double>& A, const rot_mat_3D<double>& R);
extern template quaternion<double> operator *(const quaternion<double>& Q, const axis_angle<double>& A);
extern template quaternion<double> operator *(const axis_angle<double>& A, const quaternion<double>& Q);
extern template rot_mat_3D<double> operator *(const euler_angles_TB<double>& E, const axis_angle<double>& A);
extern template rot_mat_3D<double> operator *(const axis_angle<double>& A, const euler_angles_TB<double>& E);
extern template trans_mat_3D<double> operator *(const rot_mat_3D<double>& R, const trans_mat_3D<double>& M);
extern template trans_mat_3D<double> operator *(const quaternion<double>& Q, const trans_mat_3D<double>& M);
extern template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const quaternion<double>& Q);
extern template trans_mat_3D<double> operator *(const euler_angles_TB<double>& E, const trans_mat_3D<double>& M);
extern template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const euler_angles_TB<double>& E);
extern template trans_mat_3D<double> operator *(const axis_angle<double>& A, const trans_mat_3D<double>& M);
extern template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const axis_angle<double>& A);

extern template bool operator ==(const rot_mat_3D<double>& R, const quaternion<double>& Q);
extern template bool operator !=(const rot_mat_3D<double>& R, const quaternion<double>& Q);
extern template bool operator ==(const quaternion<double>& Q, const rot_mat_3D<double>& R);
extern template bool operator !=(const quaternion<double>& Q, const rot_mat_3D<double>& R);
extern template bool operator ==(const quaternion<double>& Q, const euler_angles_TB<double>& E);
extern template bool operator !=(const quaternion<double>& Q, const euler_angles_TB<double>& E);
extern template bool operator ==(const euler_angles_TB<double>& E, const quaternion<double>& Q);
extern template bool operator !=(const euler_angles_TB<double>& E, const quaternion<double>& Q);
extern template bool operator ==(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
extern template bool operator !=(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
extern template bool operator ==(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
extern template bool operator !=(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
extern template bool operator ==(const rot_mat_3D<double>& R, const axis_angle<double>& A);
extern template bool operator !=(const rot_mat_3D<double>& R, const axis_angle<double>& A);
extern template bool operator ==(const axis_angle<double>& A, const rot_mat_3D<double>& R);
extern template bool operator !=(const axis_angle<double>& A, const rot_mat_3D<double>& R);
extern template bool operator ==(const quaternion<double>& Q, const axis_angle<double>& A);
extern template bool operator !=(const quaternion<double>& Q, const axis_angle<double>& A);
extern template bool operator ==(const axis_angle<double>& A, const quaternion<double>& Q);
extern template bool operator !=(const axis_angle<double>& A, const quaternion<double>& Q);
extern template bool operator ==(const euler_angles_TB<double>& E, const axis_angle<double>& A);
extern template bool operator !=(const euler_angles_TB<double>& E, const axis_angle<double>& A);
extern template bool operator ==(const axis_angle<double>& A, const euler_angles_TB<double>& E);
extern template bool operator !=(const axis_angle<double>& A, const euler_angles_TB<double>& E);



extern template class rot_mat_3D<float>;
extern template class quaternion<float>;
extern template class euler_angles_TB<float>;
extern template class axis_angle<float>;
extern template class trans_mat_3D<float>;

extern template std::ostream& operator << <float>(std::ostream& out_stream, const rot_mat_3D<float>& R);
extern template std::ostream& operator << <float>(std::ostream& out_stream, const quaternion<float>& Q);
extern template std::ostream& operator << <float>(std::ostream& out_stream, const euler_angles_TB<float>& E);
extern template std::ostream& operator << <float>(std::ostream& out_stream, const axis_angle<float>& A);
extern template std::ostream& operator << <float>(std::ostream& out_stream, const trans_mat_3D<float>& M);


extern template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const quaternion<float>& Q);
extern template rot_mat_3D<float> operator *(const quaternion<float>& Q, const rot_mat_3D<float>& R);
extern template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
extern template rot_mat_3D<float> operator *(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
extern template quaternion<float> operator *(const quaternion<float>& Q, const euler_angles_TB<float>& E);
extern template quaternion<float> operator *(const euler_angles_TB<float>& E, const quaternion<float>& Q);
extern template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const axis_angle<float>& A);
extern template rot_mat_3D<float> operator *(const axis_angle<float>& A, const rot_mat_3D<float>& R);
extern template quaternion<float> operator *(const quaternion<float>& Q, const axis_angle<float>& A);
extern template quaternion<float> operator *(const axis_angle<float>& A, const quaternion<float>& Q);
extern template rot_mat_3D<float> operator *(const euler_angles_TB<float>& E, const axis_angle<float>& A);
extern template rot_mat_3D<float> operator *(const axis_angle<float>& A, const euler_angles_TB<float>& E);
extern template trans_mat_3D<float> operator *(const rot_mat_3D<float>& R, const trans_mat_3D<float>& M);
extern template trans_mat_3D<float> operator *(const quaternion<float>& Q, const trans_mat_3D<float>& M);
extern template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const quaternion<float>& Q);
extern template trans_mat_3D<float> operator *(const euler_angles_TB<float>& E, const trans_mat_3D<float>& M);
extern template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const euler_angles_TB<float>& E);
extern template trans_mat_3D<float> operator *(const axis_angle<float>& A, const trans_mat_3D<float>& M);
extern template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const axis_angle<float>& A);

extern template bool operator ==(const rot_mat_3D<float>& R, const quaternion<float>& Q);
extern template bool operator !=(const rot_mat_3D<float>& R, const quaternion<float>& Q);
extern template bool operator ==(const quaternion<float>& Q, const rot_mat_3D<float>& R);
extern template bool operator !=(const quaternion<float>& Q, const rot_mat_3D<float>& R);
extern template bool operator ==(const quaternion<float>& Q, const euler_angles_TB<float>& E);
extern template bool operator !=(const quaternion<float>& Q, const euler_angles_TB<float>& E);
extern template bool operator ==(const euler_angles_TB<float>& E, const quaternion<float>& Q);
extern template bool operator !=(const euler_angles_TB<float>& E, const quaternion<float>& Q);
extern template bool operator ==(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
extern template bool operator !=(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
extern template bool operator ==(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
extern template bool operator !=(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
extern template bool operator ==(const rot_mat_3D<float>& R, const axis_angle<float>& A);
extern template bool operator !=(const rot_mat_3D<float>& R, const axis_angle<float>& A);
extern template bool operator ==(const axis_angle<float>& A, const rot_mat_3D<float>& R);
extern template bool operator !=(const axis_angle<float>& A, const rot_mat_3D<float>& R);
extern template bool operator ==(const quaternion<float>& Q, const axis_angle<float>& A);
extern template bool operator !=(const quaternion<float>& Q, const axis_angle<float>& A);
extern template bool operator ==(const axis_angle<float>& A, const quaternion<float>& Q);
extern template bool operator !=(const axis_angle<float>& A, const quaternion<float>& Q);
extern template bool operator ==(const euler_angles_TB<float>& E, const axis_angle<float>& A);
extern template bool operator !=(const euler_angles_TB<float>& E, const axis_angle<float>& A);
extern template bool operator ==(const axis_angle<float>& A, const euler_angles_TB<float>& E);
extern template bool operator !=(const axis_angle<float>& A, const euler_angles_TB<float>& E);


#endif
};


#endif
