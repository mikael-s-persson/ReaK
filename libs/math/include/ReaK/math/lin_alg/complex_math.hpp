/**
 * \file complex_math.hpp
 *
 * This library implements complex-valued mathematics. Not the most useful thing since C++ standard libraries
 * already include a complex type.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date january 2010
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

#ifndef REAK_COMPLEX_MATH_HPP
#define REAK_COMPLEX_MATH_HPP

#include <ReaK/core/serialization/archiver.hpp>

#include <cmath>


namespace ReaK {

/**
 * This template class defines a complex-valued variable in cartesian native representation.
 */
template < class T >
class complex {
public:
  typedef complex< T > self;

  T Re; ///< Holds the real part of the complex number.
  T Im; ///< Holds the imaginary part of the complex number.

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Explicit cast from a real number to a complex class.
   */
  explicit complex( T aRe ) BOOST_NOEXCEPT : Re( aRe ), Im( 0.0 ){};
  /**
   * Constructor from cartesian complex values, real and imaginary parts.
   */
  complex( T aRe, T aIm ) BOOST_NOEXCEPT : Re( aRe ), Im( aIm ){};
  /**
   * Default constructor, zero-valued.
   */
  complex() BOOST_NOEXCEPT : Re( 0.0 ), Im( 0.0 ){};

  // Copy-constructor is default.
  // Assignment operator is default.

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Assignment operator.
   */
  self& operator=( const T& R ) BOOST_NOEXCEPT {
    Re = R;
    Im = 0.0;
    return *this;
  };

  /** Addition-assignment operator. */
  self& operator+=( const self& C ) BOOST_NOEXCEPT {
    Re += C.Re;
    Im += C.Im;
    return *this;
  };

  /** Addition-assignment operator with a real value. */
  self& operator+=( const T& R ) BOOST_NOEXCEPT {
    Re += R;
    return *this;
  };

  /** Substraction-assignment operator. */
  self& operator-=( const self& C ) BOOST_NOEXCEPT {
    Re -= C.Re;
    Im -= C.Im;
    return *this;
  };

  /** Substraction-assignment operator with a real value. */
  self& operator-=( const T& R ) BOOST_NOEXCEPT {
    Re -= R;
    return *this;
  };

  /** Multiplication-assignment operator. */
  self& operator*=( const self& C ) BOOST_NOEXCEPT {
    T Re_tmp = Re * C.Re - Im * C.Im;
    Im = Re * C.Im + Im * C.Re;
    Re = Re_tmp;
    return *this;
  };

  /** Multiplication-assignment operator with a real value. */
  self& operator*=( const T& R ) BOOST_NOEXCEPT {
    Re *= R;
    Im *= R;
    return *this;
  };

  /** Division-assignment operator. */
  self& operator/=( const self& C ) BOOST_NOEXCEPT {
    T sqr_mag_inv = 1.0 / ( C.Re * C.Re + C.Im * C.Im );
    T Re_tmp = Re * C.Re + Im * C.Im;
    Im = sqr_mag_inv * ( Im * C.Re - Re * C.Im );
    Re = sqr_mag_inv * Re_tmp;
    return *this;
  };

  /** Division-assignment operator with a real value. */
  self& operator/=( const T& R ) BOOST_NOEXCEPT {
    Re /= R;
    Im /= R;
    return *this;
  };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /** Addition operator. */
  friend self operator+( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return self( C1.Re + C2.Re, C1.Im + C2.Im );
  };

  /** Addition operator with a real value. */
  friend self operator+( const self& C, const T& R ) BOOST_NOEXCEPT { return self( C.Re + R, C.Im ); };

  /** Addition operator with a real value. */
  friend self operator+( const T& R, const self& C ) BOOST_NOEXCEPT { return self( C.Re + R, C.Im ); };

  /** Negation operator. */
  self operator-() const BOOST_NOEXCEPT { return self( -Re, -Im ); };

  /** Substraction operator. */
  friend self operator-( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return self( C1.Re - C2.Re, C1.Im - C2.Im );
  };

  /** Substraction operator with a real value. */
  friend self operator-( const self& C, const T& R ) BOOST_NOEXCEPT { return self( C.Re - R, C.Im ); };

  /** Substraction operator with a real value. */
  friend self operator-( const T& R, const self& C ) BOOST_NOEXCEPT { return self( R - C.Re, -C.Im ); };

  /** Multiplication operator. */
  friend self operator*(const self& C1, const self& C2)BOOST_NOEXCEPT {
    return self( C1.Re * C2.Re - C1.Im * C2.Im, C1.Re * C2.Im + C1.Im * C2.Re );
  };

  /** Multiplication operator with a real value. */
  friend self operator*(const self& C, const T& R)BOOST_NOEXCEPT { return self( C.Re * R, C.Im * R ); };

  /** Division operator. */
  friend self operator/( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    T sqr_mag_inv = 1.0 / ( C2.Re * C2.Re + C2.Im * C2.Im );
    return self( sqr_mag_inv * ( C1.Re * C2.Re + C1.Im * C2.Im ), sqr_mag_inv * ( C1.Im * C2.Re - C1.Re * C2.Im ) );
  };

  /** Division operator with a real value. */
  friend self operator/( const self& C, const T& R ) BOOST_NOEXCEPT { return self( C.Re / R, C.Im / R ); };

  /** Division operator for a real and complex. */
  friend self operator/( const double& R, const self& C ) BOOST_NOEXCEPT {
    T sqr_mag_inv = R / ( C.Re * C.Re + C.Im * C.Im );
    return self( sqr_mag_inv * C.Re, -sqr_mag_inv * C.Im );
  };

  /*******************************************************************************
                             Comparison Operators
  *******************************************************************************/
  /** Equality-comparison operator. */
  friend bool operator==( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re == C2.Re ) && ( C1.Im == C2.Im ) );
  };

  /** Equality-comparison operator with a real value. */
  friend bool operator==( const self& C, const T& R ) BOOST_NOEXCEPT { return ( ( C.Re == R ) && ( C.Im == 0.0 ) ); };

  /** Equality-comparison operator with a real value. */
  friend bool operator==( const T& R, const self& C ) BOOST_NOEXCEPT { return ( ( C.Re == R ) && ( C.Im == 0.0 ) ); };

  /** Inequality-comparison operator. */
  friend bool operator!=( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re != C2.Re ) || ( C1.Im != C2.Im ) );
  };

  /** Inequality-comparison operator with a real value. */
  friend bool operator!=( const self& C, const T& R ) BOOST_NOEXCEPT { return ( ( C.Re != R ) || ( C.Im != 0.0 ) ); };

  /** Inequality-comparison operator with a real value. */
  friend bool operator!=( const T& R, const self& C ) BOOST_NOEXCEPT { return ( ( C.Re != R ) || ( C.Im != 0.0 ) ); };

  /** Less-than-comparison operator, magnitude-wise. */
  friend bool operator<( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re * C1.Re + C1.Im * C1.Im ) < ( C2.Re * C2.Re + C2.Im * C2.Im ) );
  };

  /** Less-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator<( const self& C, const T& R ) BOOST_NOEXCEPT {
    return ( ( C.Re * C.Re + C.Im * C.Im ) < R * R );
  };

  /** Less-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator<( const T& R, const self& C ) BOOST_NOEXCEPT {
    return ( R * R < ( C.Re * C.Re + C.Im * C.Im ) );
  };

  /** Less-or-equal-comparison operator, magnitude-wise. */
  friend bool operator<=( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re * C1.Re + C1.Im * C1.Im ) <= ( C2.Re * C2.Re + C2.Im * C2.Im ) );
  };

  /** Less-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator<=( const self& C, const T& R ) BOOST_NOEXCEPT {
    return ( ( C.Re * C.Re + C.Im * C.Im ) <= R * R );
  };

  /** Less-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator<=( const T& R, const self& C ) BOOST_NOEXCEPT {
    return ( R * R <= ( C.Re * C.Re + C.Im * C.Im ) );
  };

  /** Greater-than-comparison operator, magnitude-wise. */
  friend bool operator>( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re * C1.Re + C1.Im * C1.Im ) > ( C2.Re * C2.Re + C2.Im * C2.Im ) );
  };

  /** Greater-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator>( const self& C, const T& R ) BOOST_NOEXCEPT {
    return ( ( C.Re * C.Re + C.Im * C.Im ) > R * R );
  };

  /** Greater-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator>( const T& R, const self& C ) BOOST_NOEXCEPT {
    return ( R * R > ( C.Re * C.Re + C.Im * C.Im ) );
  };

  /** Greater-or-equal-comparison operator, magnitude-wise. */
  friend bool operator>=( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return ( ( C1.Re * C1.Re + C1.Im * C1.Im ) >= ( C2.Re * C2.Re + C2.Im * C2.Im ) );
  };

  /** Greater-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator>=( const self& C, const T& R ) BOOST_NOEXCEPT {
    return ( ( C.Re * C.Re + C.Im * C.Im ) >= R * R );
  };

  /** Greater-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator>=( const T& R, const self& C ) BOOST_NOEXCEPT {
    return ( R * R >= ( C.Re * C.Re + C.Im * C.Im ) );
  };


  /** Complex conjugate for a complex value. */
  friend self conj( const self& x ) BOOST_NOEXCEPT { return self( x.Re, -x.Im ); };

  /** Complex square magnitude for a complex value. */
  friend T norm_2( const self& x ) BOOST_NOEXCEPT {
    using std::sqrt;
    return sqrt( x.Re * x.Re + x.Im * x.Im );
  };

  /** Complex square magnitude for a complex value. */
  friend T norm_2_sqr( const self& x ) BOOST_NOEXCEPT { return T( x.Re * x.Re + x.Im * x.Im ); };

  /** Complex square magnitude for a complex value. */
  friend T sqr_mag( const self& x ) BOOST_NOEXCEPT { return T( x.Re * x.Re + x.Im * x.Im ); };


  // Exponential and logarithmic functions:

  /** Compute exponential function (function), for a complex value. */
  friend self exp( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::sin;
    if( x.Re == 0.0 )
      return self( cos( x.Im ), sin( x.Im ) );
    else
      return exp( x.Re ) * self( cos( x.Im ), sin( x.Im ) );
  };

  /** Compute natural logarithm (function), for a complex value. */
  friend self log( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    using std::sqrt;
    using std::atan2;
    return self( log( sqrt( x.Re * x.Re + x.Im * x.Im ) ), atan2( x.Im, x.Re ) );
  };

  /** Compute common logarithm (function), for a complex value. */
  friend self log10( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    return log( x ) / log( T( 10.0 ) );
  };

  // Trigonometric functions:

  /** Compute cosine (function), for a complex value.*/
  friend self cos( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    return self( cos( x.Re ) * cosh( x.Im ), -sin( x.Re ) * sinh( x.Im ) );
  };

  /** Compute sine (function), for a complex value.*/
  friend self sin( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    return self( sin( x.Re ) * cosh( x.Im ), cos( x.Re ) * sinh( x.Im ) );
  };

  /** Compute tangent (function), for a complex value.*/
  friend self tan( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    T tmp = 1.0 / ( cos( 2.0 * x.Re ) + cosh( 2.0 * x.Im ) );
    return self( sin( 2.0 * x.Re ) * tmp, sinh( 2.0 * x.Im ) * tmp );
  };

  /** Compute arc cosine (function), for a complex value.*/
  friend self acos( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    using std::sqrt;
    return self( 0.0, -1.0 ) * log( x + sqrt( x * x - 1.0 ) );
  };

  /** Compute arc sine (function), for a complex value.*/
  friend self asin( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    using std::sqrt;
    return self( 0.0, -1.0 ) * log( self( -x.Im, x.Re ) + sqrt( 1.0 - x * x ) );
  };

  /** Compute arc tangent (function), for a complex value.*/
  friend self atan( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    return self( 0.0, 0.5 ) * log( self( x.Im + 1.0, -x.Re ) / self( 1.0 - x.Im, x.Re ) );
  };

  /** Compute arc tangent with two parameters (function), for a complex value.*/
  friend self atan2( const self& y, const self& x ) BOOST_NOEXCEPT { return atan( y / x ); };

  // Hyperbolic functions:

  /** Compute hyperbolic cosine (function), for a complex value.*/
  friend self cosh( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    return self( cosh( x.Re ) * cos( x.Im ), sinh( x.Re ) * sin( x.Im ) );
  };

  /** Compute hyperbolic sine (function), for a complex value.*/
  friend self sinh( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    return self( cos( x.Im ) * sinh( x.Re ), sin( x.Im ) * cosh( x.Re ) );
  };

  /** Compute hyperbolic tangent (function), for a complex value.*/
  friend self tanh( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    T tmp = 1.0 / ( cosh( 2.0 * x.Re ) + cos( 2.0 * x.Im ) );
    return self( sinh( 2.0 * x.Re ) * tmp, sin( 2.0 * x.Im ) * tmp );
  };


  // Power functions

  /** Raise to power (function), for a complex value.*/
  friend self pow( const self& base, const self& exponent ) BOOST_NOEXCEPT { return exp( exponent * log( base ) ); };

  /** Compute square root (function), for a complex value.*/
  friend self sqrt( const self& x ) BOOST_NOEXCEPT {
    using std::cos;
    using std::sin;
    using std::pow;
    using std::atan2;
    T angle = atan2( x.Im, x.Re );
    return pow( x.Re * x.Re + x.Im * x.Im, 0.25 ) * self( cos( 0.5 * angle ), sin( 0.5 * angle ) );
  };


  // Rounding, absolute value and remainder functions:

  /** Round up value (function), for a complex value. */
  friend self ceil( const self& x ) BOOST_NOEXCEPT {
    using std::ceil;
    return self( ceil( x.Re ), ceil( x.Im ) );
  };

  /** Compute absolute value (function), for a complex value. */
  friend T fabs( const self& x ) BOOST_NOEXCEPT {
    using std::sqrt;
    return sqrt( x.Re * x.Re + x.Im * x.Im );
  };

  /** Round down value (function), for a complex value.*/
  friend self floor( const self& x ) BOOST_NOEXCEPT {
    using std::floor;
    return self( floor( x.Re ), floor( x.Im ) );
  };


  /// Loading a complex value.
  friend serialization::iarchive& RK_CALL operator>>( serialization::iarchive& in, self& C ) {
    return in >> C.Re >> C.Im;
  };

  /// Loading a complex value with a name.
  friend serialization::iarchive& RK_CALL
    operator&( serialization::iarchive& in, const std::pair< std::string, self& >& f ) {
    return in & std::pair< std::string, T& >( f.first + "_Re", f.second.Re )
           & std::pair< std::string, T& >( f.first + "_Im", f.second.Im );
  };

  /// Saving a complex value.
  friend serialization::oarchive& RK_CALL operator<<( serialization::oarchive& out, const self& C ) {
    return out << C.Re << C.Im;
  };

  /// Saving a complex value with a name.
  friend serialization::oarchive& RK_CALL
    operator&( serialization::oarchive& out, const std::pair< std::string, const self& >& C ) {
    return out & std::pair< std::string, T >( C.first + "_Re", C.second.Re )
           & std::pair< std::string, T >( C.first + "_Im", C.second.Im );
  };
};


namespace rtti {

template < typename T >
struct get_type_id< complex< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000007 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "ReaK::complex" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "ReaK::complex"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const complex< T >& save_type;
  typedef complex< T >& load_type;
};
};


/**
 * This template class defines a complex-valued variable in cartesian native representation.
 */
template < class T >
class unit_complex : public complex< T > {
public:
  typedef complex< T > base_type;
  typedef unit_complex< T > self;

  typedef T scalar_type;
  typedef T vector_type;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Constructor from cartesian complex values, real and imaginary parts.
   */
  unit_complex( scalar_type aRe, scalar_type aIm ) BOOST_NOEXCEPT : base_type( aRe, aIm ) {
    using std::sqrt;
    scalar_type mag = sqrt( this->Re * this->Re + this->Im * this->Im );
    this->Re /= mag;
    this->Im /= mag;
  };
  /**
   * Default constructor, zero-valued.
   */
  unit_complex() BOOST_NOEXCEPT : Re( 1.0 ), Im( 0.0 ){};

  // Copy-constructor is default.
  // Assignment operator is default.

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /** Multiplication-assignment operator. */
  self& operator*=( const self& C ) BOOST_NOEXCEPT {
    scalar_type Re_tmp = this->Re * C.Re - this->Im * C.Im;
    this->Im = this->Re * C.Im + this->Im * C.Re;
    this->Re = Re_tmp;
    return *this;
  };

  /** Division-assignment operator. */
  self& operator/=( const self& C ) BOOST_NOEXCEPT {
    scalar_type Re_tmp = this->Re * C.Re + this->Im * C.Im;
    this->Im = this->Im * C.Re - this->Re * C.Im;
    this->Re = Re_tmp;
    return *this;
  };

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /** Negation operator. */
  self operator-() const BOOST_NOEXCEPT { return self( -this->Re, -this->Im ); };

  /** Multiplication operator. */
  friend self operator*(const self& C1, const self& C2)BOOST_NOEXCEPT {
    return self( C1.Re * C2.Re - C1.Im * C2.Im, C1.Re * C2.Im + C1.Im * C2.Re );
  };

  /** Division operator. */
  friend self operator/( const self& C1, const self& C2 ) BOOST_NOEXCEPT {
    return self( C1.Re * C2.Re + C1.Im * C2.Im, C1.Im * C2.Re - C1.Re * C2.Im );
  };

  /** Complex conjugate for a complex value. */
  friend self conj( const self& x ) BOOST_NOEXCEPT { return self( x.Re, -x.Im ); };

  /** Complex square magnitude for a complex value. */
  friend scalar_type norm_2( const self& x ) BOOST_NOEXCEPT { return scalar_type( 1.0 ); };

  /** Complex square magnitude for a complex value. */
  friend scalar_type norm_2_sqr( const self& x ) BOOST_NOEXCEPT { return scalar_type( 1.0 ); };

  /** Complex square magnitude for a complex value. */
  friend scalar_type sqr_mag( const self& x ) BOOST_NOEXCEPT { return scalar_type( 1.0 ); };

  // Exponential and logarithmic functions:

  /** Compute exponential function (function), for a complex value. */
  //   friend self exp( const vector_type& x ) BOOST_NOEXCEPT {
  //     using std::cos; using std::sin;
  //     return self( cos( x ), sin( x ) );
  //   };

  /** Compute natural logarithm (function), for a complex value. */
  friend vector_type log( const self& x ) BOOST_NOEXCEPT {
    using std::atan2;
    return atan2( x.Im, x.Re );
  };

  /** Compute common logarithm (function), for a complex value. */
  friend vector_type log10( const self& x ) BOOST_NOEXCEPT {
    using std::log;
    return log( x ) / log( scalar_type( 10.0 ) );
  };

  // Power functions

  /** Raise to power (function), for a complex value.*/
  friend self pow( const self& base, const self& exponent ) BOOST_NOEXCEPT {
    scalar_type angle = atan2( x.Im, x.Re );
    return self( cos( exponent * angle ), sin( exponent * angle ) );
    //     return exp( exponent * log( base ) );
  };

  /** Compute square root (function), for a complex value.*/
  friend self sqrt( const self& x ) BOOST_NOEXCEPT {
    scalar_type angle = atan2( x.Im, x.Re );
    return self( cos( 0.5 * angle ), sin( 0.5 * angle ) );
    //     return exp( scalar_type(0.5) * log( x ) );
  };

  // Rounding, absolute value and remainder functions:

  /** Compute absolute value (function), for a complex value. */
  friend scalar_type fabs( const self& x ) BOOST_NOEXCEPT { return scalar_type( 1.0 ); };
};


/** Complex conjugate for a real value. */
inline float conj( float x ) BOOST_NOEXCEPT { return x; };
/** Complex conjugate for a real value. */
inline double conj( double x ) BOOST_NOEXCEPT { return x; };

/** Complex square magnitude for a real value. */
inline float norm_2( float x ) BOOST_NOEXCEPT {
  using std::fabs;
  return fabs( x );
};
/** Complex square magnitude for a real value. */
inline double norm_2( double x ) BOOST_NOEXCEPT {
  using std::fabs;
  return fabs( x );
};

/** Complex square magnitude for a real value. */
inline float norm_2_sqr( float x ) BOOST_NOEXCEPT { return x * x; };
/** Complex square magnitude for a real value. */
inline double norm_2_sqr( double x ) BOOST_NOEXCEPT { return x * x; };

/** Complex square magnitude for a real value. */
inline float sqr_mag( float x ) BOOST_NOEXCEPT { return x * x; };
/** Complex square magnitude for a real value. */
inline double sqr_mag( double x ) BOOST_NOEXCEPT { return x * x; };
};


#endif // RK_COMPLEX_MATH_HPP
