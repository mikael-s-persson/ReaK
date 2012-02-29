/**
 * \file quat_alg.hpp
 *
 * This library implements quaternionic algebra.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
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

#ifndef REAK_QUAT_ALG_HPP
#define REAK_QUAT_ALG_HPP

#include "lin_alg/vect_alg.hpp"
#include "serialization/archiver.hpp"

#include "rotations_3D.hpp"

#include <cmath>

namespace ReaK {


/**
 * This template class defines a quaternion-valued variable (not a unit-quaternion for representing rotations).
 */
template <class T>
class quat {
  public:
    typedef quat<T> self;
    
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef void allocator_type;
    
    typedef pointer iterator;
    typedef const_pointer const_iterator;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    typedef T scalar_type;
    typedef vect<T,3> vector_type;
          
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 4);
    
    
    value_type q[4]; ///< Holds the four values of the quaternion (q[0] is the real part).
  
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
  
    /**
     * Explicit cast from a real number to a quat.
     */
    explicit quat(const scalar_type& aRe = scalar_type()) { q[0] = aRe; q[3] = q[2] = q[1] = value_type(); };
    
    /**
     * Constructor from cartesian complex values, real and imaginary parts.
     */
    quat(const scalar_type& aRe, const vector_type& aIm) { q[0] = aRe; q[1] = aIm[0]; q[2] = aIm[1]; q[3] = aIm[2]; };
  
    /**
     * Converts a vector into a pure quaternion.
     */
    explicit quat(const vector_type& aIm) { q[0] = scalar_type(); q[1] = aIm[0]; q[2] = aIm[1]; q[3] = aIm[2]; };
    
    /**
     * Converts a 4D vector into a quaternion.
     */
    explicit quat(const vect<value_type,4>& V) { q[0] = V[0]; q[1] = V[1]; q[2] = V[2]; q[3] = V[3]; };
    
    /**
     * Convstructs a quaternion from 4 components.
     */
    quat(const_reference q0, const_reference q1, const_reference q2, const_reference q3) { q[0] = q0; q[1] = q1; q[2] = q2; q[3] = q3; };
    
    //Copy-constructor is default.
    //Assignment operator is default.
    
    
    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      if(i >= 4)
	throw std::range_error("Quaternion index out of range.");
      return q[i];
    };
    
    /**
     * Array indexing operator, accessor for lvalue.
     * \test PASSED
     */
    reference operator [](size_type i) {
      if(i >= 4)
	throw std::range_error("Quaternion index out of range.");
      return q[i];
    };
    
    /**
     * Returns the size of the quaternion (viewed as a 4D vector).
     */
    size_type size() const { return 4; };
    
    /**
     * Returns an iterator to the first element of the quaternion (viewed as a 4D vector).
     */
    iterator begin() { return q; };
    /**
     * Returns a const-iterator to the first element of the quaternion (viewed as a 4D vector).
     */
    const_iterator begin() const { return q; };
    /**
     * Returns an iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
     */
    iterator end() { return q + 4; };
    /**
     * Returns a const-iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
     */
    const_iterator end() const { return q + 4; };
    
      
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Assignment operator.
     */
    self& operator =(const scalar_type& R) {
      q[0] = R;
      q[3] = q[2] = q[1] = value_type();
      return *this; 
    };
    
    /**
     * Assignment operator.
     */
    self& operator =(const vector_type& I) {
      q[0] = value_type();
      q[1] = I[0];
      q[2] = I[1];
      q[3] = I[2];
      return *this; 
    };
  
    /** Addition-assignment operator. */
    self& operator +=(const self& C) {
      q[0] += C.q[0];
      q[1] += C.q[1];
      q[1] += C.q[2];
      q[1] += C.q[3];
      return *this;
    };
  
    /** Addition-assignment operator with a real value. */
    self& operator +=(const scalar_type& R) {
      q[0] += R;
      return *this;
    };

    /** Addition-assignment operator with a real value. */
    self& operator +=(const vector_type& V) {
      q[1] += V[0];
      q[2] += V[1];
      q[3] += V[2];
      return *this;
    };

    /** Substraction-assignment operator. */
    self& operator -=(const self& C) {
      q[0] -= C.q[0];
      q[1] -= C.q[1];
      q[1] -= C.q[2];
      q[1] -= C.q[3];
      return *this;
    };
  
    /** Substraction-assignment operator with a real value. */
    self& operator -=(const scalar_type& R) {
      q[0] -= R;
      return *this;
    };

    /** Substraction-assignment operator with a real value. */
    self& operator -=(const vector_type& V) {
      q[1] -= V[0];
      q[2] -= V[1];
      q[3] -= V[2];
      return *this;
    };

    /** Multiplication-assignment operator. */
    self& operator *=(const self& C) {
      return (*this = ((*this) * C));
    };

    /** Multiplication-assignment operator with a real value. */
    self& operator *=(const scalar_type& R) {
      q[0] *= R;
      q[1] *= R;
      q[2] *= R;
      q[3] *= R;
      return *this;
    };
    
  
/*******************************************************************************
                         Basic Operators
*******************************************************************************/
  
    /** Addition operator. */
    friend self operator +(const self& C1, const self& C2) {
      return self(C1.q[0] + C2.q[0], C1.q[1] + C2.q[1], C1.q[2] + C2.q[2], C1.q[3] + C2.q[3]);
    };

    /** Addition operator. */
    friend self operator +(const self& C1, const scalar_type& C2) {
      return self(C1.q[0] + C2, C1.q[1], C1.q[2], C1.q[3]);
    };

    /** Addition operator. */
    friend self operator +(const scalar_type& C1, const self& C2) {
      return self(C2.q[0] + C1, C2.q[1], C2.q[2], C2.q[3]);
    };
    
    /** Addition operator. */
    friend self operator +(const self& C1, const vector_type& C2) {
      return self(C1.q[0], C1.q[1] + C2[0], C1.q[2] + C2[1], C1.q[3] + C2[2]);
    };
    
    /** Addition operator. */
    friend self operator +(const vector_type& C1, const self& C2) {
      return self(C2.q[0], C1[0] + C2.q[1], C1[1] + C2.q[2], C1[2] + C2.q[3]);
    };
  
    /** Negation operator. */
    friend self operator -(const self& Q) {
      return self(-Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3]);
    };
  
    /** Substraction operator. */
    friend self operator -(const self& C1, const self& C2) {
      return self(C1.q[0] - C2.q[0], C1.q[1] - C2.q[1], C1.q[2] - C2.q[2], C1.q[3] - C2.q[3]);
    };

    /** Substraction operator. */
    friend self operator -(const self& C1, const scalar_type& C2) {
      return self(C1.q[0] - C2, C1.q[1], C1.q[2], C1.q[3]);
    };

    /** Substraction operator. */
    friend self operator -(const scalar_type& C1, const self& C2) {
      return self(C1 - C2.q[0], -C2.q[1], -C2.q[2], -C2.q[3]);
    };
  
    /** Substraction operator. */
    friend self operator -(const self& C1, const vector_type& C2) {
      return self(C1.q[0], C1.q[1] - C2[0], C1.q[2] - C2[1], C1.q[3] - C2[2]);
    };

    /** Substraction operator. */
    friend self operator -(const vector_type& C1, const self& C2) {
      return self(-C2.q[0], C1[0] - C2.q[1], C1[1] - C2.q[2], C1[2] - C2.q[3]);
    };

    /**
     * Multiplication by a quaternion.
     * \test PASSED
     */
    friend
    self operator *(const self& Q1, const self& Q2) {
      return self(Q2.q[0] * Q1.q[0] - Q2.q[1] * Q1.q[1] - Q2.q[2] * Q1.q[2] - Q2.q[3] * Q1.q[3],
                  Q2.q[0] * Q1.q[1] + Q2.q[3] * Q1.q[2] - Q2.q[2] * Q1.q[3] + Q2.q[1] * Q1.q[0],
                  Q2.q[0] * Q1.q[2] - Q2.q[3] * Q1.q[1] + Q2.q[1] * Q1.q[3] + Q2.q[2] * Q1.q[0],
                  Q2.q[0] * Q1.q[3] + Q2.q[2] * Q1.q[1] - Q2.q[1] * Q1.q[2] + Q2.q[3] * Q1.q[0]);
    };
  
    /**
     * Multiplication by a scalar.
     * \test PASSED
     */
    friend
    self operator *(const self& Q1, const scalar_type& Q2) {
      return self(Q2 * Q1.q[0], Q2 * Q1.q[1], Q2 * Q1.q[2], Q2 * Q1.q[3]);
    };
  
    /**
     * Multiplication by a scalar.
     * \test PASSED
     */
    friend
    self operator *(const scalar_type& Q1, const self& Q2) {
      return self(Q1 * Q2.q[0], Q1 * Q2.q[1], Q1 * Q2.q[2], Q1 * Q2.q[3]);
    };
  
    /**
     * Multiplication by a quaternion.
     * \test PASSED
     */
    friend
    self operator *(const self& Q1, const vector_type& Q2) {
      return self(-Q2[0] * Q1.q[1] - Q2[1] * Q1.q[2] - Q2[2] * Q1.q[3],
                   Q2[2] * Q1.q[2] - Q2[1] * Q1.q[3] + Q2[0] * Q1.q[0],
                  -Q2[2] * Q1.q[1] + Q2[0] * Q1.q[3] + Q2[1] * Q1.q[0],
                   Q2[1] * Q1.q[1] - Q2[0] * Q1.q[2] + Q2[2] * Q1.q[0]);
    };
    
    /**
     * Multiplication by a quaternion.
     * \test PASSED
     */
    friend
    self operator *(const vector_type& Q1, const self& Q2) {
      return self(-Q2.q[1] * Q1[0] - Q2.q[2] * Q1[1] - Q2.q[3] * Q1[2],
                   Q2.q[0] * Q1[0] + Q2.q[3] * Q1[1] - Q2.q[2] * Q1[2],
                   Q2.q[0] * Q1[1] - Q2.q[3] * Q1[0] + Q2.q[1] * Q1[2],
                   Q2.q[0] * Q1[2] + Q2.q[2] * Q1[0] - Q2.q[1] * Q1[1]);
    };
  
/*******************************************************************************
                           Comparison Operators
*******************************************************************************/
 
    /** Equality-comparison operator. */
    friend bool operator ==(const self& C1, const self& C2) {
      return ((C1.q[0] == C2.q[0]) && (C1.q[1] == C2.q[1]) && (C1.q[2] == C2.q[2]) && (C1.q[3] == C2.q[3]));
    };
  
    /** Equality-comparison operator with a real value. */
    friend bool operator ==(const self& C, const scalar_type& R) {
      return ((C.q[0] == R) && (C.q[1] == value_type()) && (C.q[2] == value_type()) && (C.q[3] == value_type()));
    };
  
    /** Equality-comparison operator with a real value. */
    friend bool operator ==(const scalar_type& R, const self& C) {
      return ((C.q[0] == R) && (C.q[1] == value_type()) && (C.q[2] == value_type()) && (C.q[3] == value_type()));
    };
  
    /** Equality-comparison operator with a vector value. */
    friend bool operator ==(const self& C, const vector_type& V) {
      return ((C.q[0] == scalar_type()) && (C.q[1] == V[0]) && (C.q[2] == V[1]) && (C.q[3] == V[2]));
    };
  
    /** Equality-comparison operator with a vector value. */
    friend bool operator ==(const vector_type& V, const self& C) {
      return ((C.q[0] == scalar_type()) && (C.q[1] == V[0]) && (C.q[2] == V[1]) && (C.q[3] == V[2]));
    };

    /** Inequality-comparison operator. */
    friend bool operator !=(const self& C1, const self& C2) {
      return ((C1.q[0] != C2.q[0]) || (C1.q[1] != C2.q[1]) || (C1.q[2] != C2.q[2]) || (C1.q[3] != C2.q[3]));
    };
  
    /** Inequality-comparison operator with a real value. */
    friend bool operator !=(const self& C, const scalar_type& R) {
      return ((C.q[0] != R) || (C.q[1] != value_type()) || (C.q[2] != value_type()) || (C.q[3] != value_type()));
    };
  
    /** Inequality-comparison operator with a real value. */
    friend bool operator !=(const scalar_type& R, const self& C) {
      return ((C.q[0] != R) || (C.q[1] != value_type()) || (C.q[2] != value_type()) || (C.q[3] != value_type()));
    };
  
    /** Inequality-comparison operator with a vector value. */
    friend bool operator !=(const self& C, const vector_type& V) {
      return ((C.q[0] == scalar_type()) || (C.q[1] != V[0]) || (C.q[2] != V[1]) || (C.q[3] != V[2]));
    };
  
    /** Inequality-comparison operator with a vector value. */
    friend bool operator !=(const vector_type& V, const self& C) {
      return ((C.q[0] != scalar_type()) || (C.q[1] != V[0]) || (C.q[2] != V[1]) || (C.q[3] != V[2]));
    };    
    
  
    /** Quaternionic conjugate for a quaternion value. */
    friend self conj(const self& x) {
      return self(x.q[0],-x.q[1],-x.q[2],-x.q[3]);
    };
    
    /**
     * Square magnitude of the quaternion.
     * \test PASSED
     */
    friend value_type norm_2_sqr(const self& v) {
      return v.q[0] * v.q[0] + v.q[1] * v.q[1] + v.q[2] * v.q[2] + v.q[3] * v.q[3];
    };

    /**
     * Magnitude of the quaternion.
     * \test PASSED
     */
    friend value_type norm_2(const self& v) {
      using std::sqrt;
      return sqrt( norm_2_sqr(v) );
    };


    /**
     * Unit quaternion in the same direction.
     * \test PASSED
     */
    friend self unit(const self& v) {
      return v * (1.0 / norm_2(v));
    };

    /**
     * Checks if two quaternions are colinear.
     * \test PASSED
     */
    friend bool colinear(const self& v1, const self& v2) {
      using std::fabs;
      T tmp_mag1 = norm_2(v1);
      T tmp_mag2 = norm_2(v2);
      T tmp_comb = norm_2(v1 + v2);
      return (((tmp_mag1 + tmp_mag2) * (T(1.0) - std::numeric_limits<T>::epsilon()) <= tmp_comb) || 
              (fabs(tmp_mag1 - tmp_mag2) * (T(1.0) + std::numeric_limits<T>::epsilon()) >= tmp_comb));
    };
  
  

    //Exponential and logarithmic functions:

    /** Compute exponential function (function), for a quaternion value. */
    friend self exp(const self& x) {
      using std::sin; using std::cos;
      using std::exp;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(exp(x.q[0]));
      value_type fact = sin(theta) / theta;
      return exp(x.q[0]) * self(cos(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute natural logarithm (function), for a quaternion value. */
    friend self log(const self& x) {
      using std::atan2;
      using std::log;
      using std::fabs;
      using std::sqrt;
      value_type st = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(st < std::numeric_limits<value_type>::epsilon())
	return self(log(fabs(x.q[0])));
      value_type fact = atan2(st,x.q[0]) / st;
      return self(log(sqrt(x.q[0] * x.q[0] + st * st)), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };
    
    
    //Power functions

    /** Raise to power (function), for a quaternion value.*/
    friend self pow(const self& base, const self& exponent) {
      return exp(exponent * log(base));
    };  

    /** Compute square root (function), for a quaternion value.*/
    friend self sqrt(const self& x) {
      return exp( 0.5 * log(x) );
      //T angle = atan2(x.Im,x.Re);
      //return pow(x.Re * x.Re + x.Im * x.Im,0.25) * complex<T>(cos(0.5 * angle),sin(0.5 * angle));
    };  
    
    /**
     * Inverts the quaternion.
     */
    friend self invert(const self& x) {
      value_type tmp = 1.0 / norm_2_sqr(x);
      return self(x.q[0] * tmp, -x.q[1] * tmp, -x.q[2] * tmp, -x.q[3] * tmp);
    };


    //Rounding, absolute value and remainder functions:

    /** Round up value (function), for a quaternion value. */
    friend self ceil(const self& x) {
      return self(ceil(x.q[0]),ceil(x.q[1]),ceil(x.q[2]),ceil(x.q[3]));
    };  

    /** Compute absolute value (function), for a quaternion value. */
    friend value_type fabs(const self& x) {
      return norm_2(x);
    };  

    /** Round down value (function), for a quaternion value.*/
    friend self floor(const self& x) {
      return self(floor(x.q[0]),floor(x.q[1]),floor(x.q[2]),floor(x.q[3]));
    };
  

    //Trigonometric functions:

    /** Compute cosine (function), for a quaternion value.*/
    friend self cos(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(cos(x.q[0]));
      value_type fact = - sin(x.q[0]) * sinh(theta) / theta;
      return self(cos(x.q[0]) * cosh(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    }; 

    /** Compute sine (function), for a quaternion value.*/
    friend self sin(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(sin(x.q[0]));
      value_type fact = cos(x.q[0]) * sinh(theta) / theta;
      return self(sin(x.q[0]) * cosh(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute tangent (function), for a quaternion value.*/
    friend self tan(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt; using std::tan;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(tan(x.q[0]));
      value_type tmp = 1.0 / (cos(2.0 * x.q[0]) + cosh(2.0 * theta));
      value_type fact = sinh(2.0*theta) * tmp / theta;
      return self(sin(2.0*x.q[0]) * tmp, fact * x.q[1], fact * x.q[2], fact * x.q[3]); 
    };  

    /** Compute arc cosine (function), for a quaternion value.*/
    friend self acos(const self& x) {
      using std::sqrt;
      using std::pow; using std::log;
      using std::atan2; using std::cos; using std::sin; using std::acos;
      value_type ss_sht = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(ss_sht < std::numeric_limits<value_type>::epsilon())
	return self(acos(x.q[0]));
      
      //complex acos...
      value_type Re1 = x.q[0] * x.q[0] - ss_sht * ss_sht - 1.0; // x * x - 1.0
      value_type Im1 = 2.0 * x.q[0] * ss_sht;                   // ...........
      value_type angle = atan2(Im1,Re1);       //sqrt( x * x - 1.0 ) + x
      Im1 = pow(Re1 * Re1 + Im1 * Im1, 0.25);  //.......................
      Re1 = cos(0.5 * angle) * Im1 + x.q[0];   //.......................
      Im1 = Im1 * sin(0.5 * angle) + ss_sht;   //.......................
      
      angle = atan2(Im1,Re1);                  //(-1.0i) * log()
      Im1 = -log(sqrt(Re1 * Re1 + Im1 * Im1)); //...............
      Re1 = angle;                             //...............
      //end complex acos.
      
      value_type fact = Im1 / ss_sht;
      return  self(Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };

    /** Compute arc sine (function), for a quaternion value.*/
    friend self asin(const self& x) {
      using std::sqrt;
      using std::pow; using std::log;
      using std::atan2; using std::cos; using std::sin; using std::asin;
      value_type cs_sht = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(cs_sht < std::numeric_limits<value_type>::epsilon())
	return self(asin(x.q[0]));
      
      //complex asin...
      value_type Re1 = 1.0 - x.q[0] * x.q[0] + cs_sht * cs_sht; // 1.0 - x * x
      value_type Im1 = -2.0 * x.q[0] * cs_sht;                   // ...........
      value_type angle = atan2(Im1,Re1);       //sqrt( 1.0 - x * x ) + x.Re i - x.Im
      Im1 = pow(Re1 * Re1 + Im1 * Im1, 0.25);  //.......................
      Re1 = cos(0.5 * angle) * Im1 - cs_sht;   //.......................
      Im1 = Im1 * sin(0.5 * angle) + x.q[0];   //.......................
      
      angle = atan2(Im1,Re1);                  //(-1.0i) * log()
      Im1 = -log(sqrt(Re1 * Re1 + Im1 * Im1)); //...............
      Re1 = angle;                             //...............
      //end complex asin.
      
      value_type fact = Im1 / cs_sht;
      return  self(Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute arc tangent (function), for a quaternion value.*/
    friend self atan(const self& x) {
      using std::sqrt;
      using std::pow; using std::log;
      using std::atan2; using std::cos; using std::sin; using std::atan;
      value_type tmp = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(tmp < std::numeric_limits<value_type>::epsilon())
	return self(atan(x.q[0]));
            
      //complex atan...
      value_type sqr_mag_inv = 1.0 / ((1.0 - tmp) * (1.0 - tmp) + x.q[0] * x.q[0]); //complex division
      value_type Re1 = sqr_mag_inv * (1.0 - tmp * tmp - x.q[0] * x.q[0]);           //................
      value_type Im1 = -2.0 * sqr_mag_inv * x.q[0];                                 //................
      value_type angle = atan2(Im1,Re1);                  //(0.5i) * log()
      Im1 = 0.5 * log(sqrt(Re1 * Re1 + Im1 * Im1));       //..............
      Re1 = -0.5 * angle;                                 //..............
      //end complex atan.
      
      value_type fact = Im1 / tmp;
      return  self(Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute arc tangent with two parameters (function), for a quaternion value.*/
    friend self atan2(const self& y, const self& x) {
      return atan(y / x);
    };  

//Hyperbolic functions:

    /** Compute hyperbolic cosine (function), for a quaternion value.*/
    friend self cosh(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(cosh(x.q[0]));
      value_type fact = sinh(x.q[0]) * sin(theta) / theta;
      return self(cosh(x.q[0]) * cos(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute hyperbolic sine (function), for a quaternion value.*/
    friend self sinh(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(sinh(x.q[0]));
      value_type fact = cosh(x.q[0]) * sin(theta) / theta;
      return self(sinh(x.q[0]) * cos(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };  

    /** Compute hyperbolic tangent (function), for a quaternion value.*/
    friend self tanh(const self& x) {
      using std::cos; using std::sin; 
      using std::cosh; using std::sinh;
      using std::sqrt;
      value_type theta = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(theta < std::numeric_limits<value_type>::epsilon())
	return self(tanh(x.q[0]));
      value_type tmp = 1.0 / (cosh(2.0 * x.q[0]) + cos(2.0 * theta));
      value_type fact = sin(2.0*theta) * tmp / theta;
      return self(sinh(2.0*x.q[0]) * tmp, fact * x.q[1], fact * x.q[2], fact * x.q[3]); 
    };  
  
    
};










/**
 * This template class defines a quaternion-valued variable (not a unit-quaternion for representing rotations).
 */
template <class T>
class unit_quat : public quat<T> {
  public:
    typedef unit_quat<T> self;
    
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef void allocator_type;
    
    typedef pointer iterator;
    typedef const_pointer const_iterator;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    typedef T scalar_type;
    typedef vect<T,3> vector_type;
          
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 4);
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
  
    /**
     * Default constructor, always yields (1.0, 0.0, 0.0, 0.0).
     */
    unit_quat() : quat<T>(scalar_type(1.0)) { };
    
    /**
     * Constructor from quaternion value.
     */
    explicit unit_quat(const quat<T>& aQ) : quat<T>(scalar_type(1.0)) { 
      using std::sqrt;
      scalar_type factor = sqrt(aQ.q[0] * aQ.q[0] + aQ.q[1] * aQ.q[1] + aQ.q[2] * aQ.q[2] + aQ.q[3] * aQ.q[3]); 
      if( factor > std::numeric_limits<scalar_type>::epsilon() ) {
	factor = 1.0 / factor;
        this->q[0] = aQ.q[0] * factor; 
        this->q[1] = aQ.q[1] * factor; 
        this->q[2] = aQ.q[2] * factor; 
        this->q[3] = aQ.q[3] * factor;
      };
    };
    
    /**
     * Converts a 4D vector into a quaternion.
     */
    explicit unit_quat(const vect<value_type,4>& V): quat<T>(scalar_type(1.0)) { 
      using std::sqrt;
      scalar_type factor = norm_2(V); 
      if( factor > std::numeric_limits<scalar_type>::epsilon() ) {
	factor = 1.0 / factor;
        this->q[0] = V[0] * factor; 
        this->q[1] = V[1] * factor; 
        this->q[2] = V[2] * factor; 
        this->q[3] = V[3] * factor;
      };
    };
    
    /**
     * Convstructs a quaternion from 4 components.
     */
    unit_quat(const_reference q0, const_reference q1, const_reference q2, const_reference q3) : quat<T>(scalar_type(1.0)) { 
      using std::sqrt;
      scalar_type factor = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3); 
      if( factor > std::numeric_limits<scalar_type>::epsilon() ) {
	factor = 1.0 / factor;
        this->q[0] = q0 * factor; 
        this->q[1] = q1 * factor; 
        this->q[2] = q2 * factor; 
        this->q[3] = q3 * factor;
      };
    };
    
    //Copy-constructor is default.
    //Assignment operator is default.
    
    // NOTE hiding the non-const overloads in the base class is intentional here:
    
    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      if(i >= 4)
	throw std::range_error("Quaternion index out of range.");
      return this->q[i];
    };
    
    /**
     * Returns a const-iterator to the first element of the quaternion (viewed as a 4D vector).
     */
    const_iterator begin() const { return this->q; };
    /**
     * Returns a const-iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
     */
    const_iterator end() const { return this->q + 4; };
    
    /**
     * Returns the same unit-quaternion, but now as a quaternion class, which is used to represent 3D rotations.
     */
    quaternion< value_type > as_rotation() const {
      return quaternion< value_type >(this->q[0],this->q[1],this->q[2],this->q[3]);
    };
      
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Assignment operator.
     */
    self& operator =(const self& Q) {
      this->q[0] = Q.q[0];
      this->q[1] = Q.q[1];
      this->q[2] = Q.q[2];
      this->q[3] = Q.q[3];
      return *this; 
    };
    
    /** Multiplication-assignment operator. */
    self& operator *=(const self& C) {
      return (*this = ((*this) * C));
    };
  
/*******************************************************************************
                         Basic Operators
*******************************************************************************/
  
    /** Negation operator. */
    friend self operator -(self Q) {
      Q.q[0] = -Q.q[0];
      Q.q[1] = -Q.q[1];
      Q.q[2] = -Q.q[2];
      Q.q[3] = -Q.q[3];
      return Q;
    };
  
    /**
     * Multiplication by a quaternion.
     * \test PASSED
     */
    friend
    self operator *(const self& Q1, const self& Q2) {
      return self(Q2.q[0] * Q1.q[0] - Q2.q[1] * Q1.q[1] - Q2.q[2] * Q1.q[2] - Q2.q[3] * Q1.q[3],
                  Q2.q[0] * Q1.q[1] + Q2.q[3] * Q1.q[2] - Q2.q[2] * Q1.q[3] + Q2.q[1] * Q1.q[0],
                  Q2.q[0] * Q1.q[2] - Q2.q[3] * Q1.q[1] + Q2.q[1] * Q1.q[3] + Q2.q[2] * Q1.q[0],
                  Q2.q[0] * Q1.q[3] + Q2.q[2] * Q1.q[1] - Q2.q[1] * Q1.q[2] + Q2.q[3] * Q1.q[0]);
    };
  
    
  
    /** Quaternionic conjugate for a quaternion value. */
    friend self conj(self x) {
      x.q[1] = -x.q[1];
      x.q[2] = -x.q[2];
      x.q[3] = -x.q[3];
      return x;
    };
    
    /**
     * Square magnitude of the quaternion.
     * \test PASSED
     */
    friend value_type norm_2_sqr(const self& v) {
      return value_type(1.0);
    };

    /**
     * Magnitude of the quaternion.
     * \test PASSED
     */
    friend value_type norm_2(const self& v) {
      return value_type(1.0);
    };


    /**
     * Unit quaternion in the same direction.
     * \test PASSED
     */
    friend self unit(const self& v) {
      return v;
    };
  

    //Exponential and logarithmic functions:

    
    /** Compute natural logarithm (function), for a quaternion value. */
    friend vector_type log(const self& x) {
      using std::atan2;
      using std::sqrt;
      value_type st = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
      if(st < std::numeric_limits<value_type>::epsilon())
	return vector_type(value_type(0.0),value_type(0.0),value_type(0.0));
      value_type fact = atan2(st,x.q[0]) / st;
      return vector_type(fact * x.q[1], fact * x.q[2], fact * x.q[3]);
    };
    
    
    //Power functions

    /** Raise to power (function), for a quaternion value.*/
    friend self pow(const self& base, const scalar_type& exponent) {
      return exp(exponent * log(base));
    };  

    /** Compute square root (function), for a quaternion value.*/
    friend self sqrt(const self& x) {
      return exp( 0.5 * log(x) );
    };  
    
    /**
     * Inverts the quaternion.
     */
    friend self invert(const self& x) {
      return conj(x);
    };


    //Rounding, absolute value and remainder functions:

    /** Compute absolute value (function), for a quaternion value. */
    friend value_type fabs(const self& x) {
      return 1.0;
    };  
  
    
};




/** Compute exponential function (function), for a quaternion value. */
template <typename T>
unit_quat<T> exp(const vect<T,3>& x) {
  using std::sin; using std::cos;
  using std::exp;
  using std::sqrt;
  T theta = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  if(theta < std::numeric_limits<T>::epsilon())
    return unit_quat<T>();
  T fact = sin(theta) / theta;
  return unit_quat<T>(cos(theta), fact * x[0], fact * x[1], fact * x[2]);
};  






namespace serialization {
  
  /// Loading a quaternion value.
  template <typename T>
  iarchive& RK_CALL operator >>(iarchive& in, quat<T>& C) {
    return in >> C.q[0] >> C.q[1] >> C.q[2] >> C.q[3];
  };

  /// Loading a quaternion value with a name.
  template <typename T>
  iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, quat<T>& >& f) {
    return in & std::pair<std::string, T& >(f.first + "_q0",f.second.q[0])
              & std::pair<std::string, T& >(f.first + "_q1",f.second.q[1])
              & std::pair<std::string, T& >(f.first + "_q2",f.second.q[2])
              & std::pair<std::string, T& >(f.first + "_q3",f.second.q[3]);
  };
  
  /// Saving a quaternion value.
  template <typename T>
  oarchive& RK_CALL operator <<(oarchive& out, const quat<T>& C) {
    return out << C.q[0] << C.q[1] << C.q[2] << C.q[3];
  };
  
  /// Saving a quaternion value with a name.
  template <typename T>
  oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, const quat<T>& >& C) {
    return out & std::pair<std::string, T >(C.first + "_q0",C.second.q[0])
               & std::pair<std::string, T >(C.first + "_q1",C.second.q[1])
               & std::pair<std::string, T >(C.first + "_q2",C.second.q[2])
               & std::pair<std::string, T >(C.first + "_q3",C.second.q[3]);
  };
    
};
    
namespace rtti {

template <typename T>
struct get_type_id< quat<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000002B);
  static std::string type_name() { return "ReaK::quat"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const quat<T>& save_type;
  typedef quat<T>& load_type;
};

template <typename T>
struct get_type_id< unit_quat<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000002D);
  static std::string type_name() { return "ReaK::unit_quat"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const quat<T>& save_type;
  typedef quat<T>& load_type;
};

};



template <typename T>
struct is_readable_vector< quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< quat<T> > type;
};

template <typename T>
struct is_writable_vector< quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< quat<T> > type;
};

template <typename T>
struct is_resizable_vector< quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< quat<T> > type;
};


template <typename T>
struct has_allocator_vector< quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector< quat<T> > type;
};


template <typename T>
struct is_readable_vector< unit_quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< unit_quat<T> > type;
};

template <typename T>
struct is_writable_vector< unit_quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector< unit_quat<T> > type;
};

template <typename T>
struct is_resizable_vector< unit_quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< unit_quat<T> > type;
};


template <typename T>
struct has_allocator_vector< unit_quat<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector< unit_quat<T> > type;
};




};

#endif








