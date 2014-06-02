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
template <class T>
class complex {
  public:
  T Re; ///< Holds the real part of the complex number.
  T Im; ///< Holds the imaginary part of the complex number.
  
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
  
  /**
   * Explicit cast from a real number to a complex class.
   */
  explicit complex(T aRe) : Re(aRe), Im(0.0) { };
  /**
   * Constructor from cartesian complex values, real and imaginary parts.
   */
  complex(T aRe,T aIm) : Re(aRe), Im(aIm) { };
  /**
   * Default constructor, zero-valued.
   */
  complex() : Re(0.0), Im(0.0) { };
  
  //Copy-constructor is default.
  //Assignment operator is default.
  
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  /**
   * Assignment operator.
   */
  complex<T>& operator =(const T& R) {
    Re = R;
    Im = 0.0;
    return *this; 
  };
  
  /** Addition-assignment operator. */
  complex<T>& operator +=(const complex<T>& C) {
    Re += C.Re;
    Im += C.Im;
    return *this;
  };
  
  /** Addition-assignment operator with a real value. */
  complex<T>& operator +=(const T& R) {
    Re += R;
    return *this;
  };

  /** Substraction-assignment operator. */
  complex<T>& operator -=(const complex<T>& C) {
    Re -= C.Re;
    Im -= C.Im;
    return *this;
  };
  
  /** Substraction-assignment operator with a real value. */
  complex<T>& operator -=(const T& R) {
    Re -= R;
    return *this;
  };

  /** Multiplication-assignment operator. */
  complex<T>& operator *=(const complex<T>& C) {
    T Re_tmp = Re * C.Re - Im * C.Im;
    Im = Re * C.Im + Im * C.Re;
    Re = Re_tmp;
    return *this;
  };

  /** Multiplication-assignment operator with a real value. */
  complex<T>& operator *=(const T& R) {
    Re *= R;
    Im *= R;
    return *this;
  };
  
  /** Division-assignment operator. */
  complex<T>& operator /=(const complex<T>& C) {
    T sqr_mag_inv = 1.0 / (C.Re * C.Re + C.Im * C.Im);
    T Re_tmp = Re * C.Re + Im * C.Im;
    Im = sqr_mag_inv * (Im * C.Re - Re * C.Im);
    Re = sqr_mag_inv * Re_tmp;
    return *this;
  };
  
  /** Division-assignment operator with a real value. */
  complex<T>& operator /=(const T& R) {
    Re /= R;
    Im /= R;
    return *this;
  };
  
/*******************************************************************************
                         Basic Operators
*******************************************************************************/
  
  /** Addition operator. */
  friend complex<T> operator +(const complex<T>& C1, const complex<T>& C2) {
    return complex<T>(C1.Re + C2.Re, C1.Im + C2.Im);
  };

  /** Addition operator with a real value. */
  friend complex<T> operator +(const complex<T>& C, const T& R) {
    return complex<T>(C.Re + R, C.Im);
  };
  
  /** Addition operator with a real value. */
  friend complex<T> operator +(const T& R, const complex<T>& C) {
    return complex<T>(C.Re + R, C.Im);
  };
  
  /** Negation operator. */
  complex<T> operator -() const {
    return complex<T>(-Re,-Im);
  };
  
  /** Substraction operator. */
  friend complex<T> operator -(const complex<T>& C1, const complex<T>& C2) {
    return complex<T>(C1.Re - C2.Re, C1.Im - C2.Im);
  };
  
  /** Substraction operator with a real value. */
  friend complex<T> operator -(const complex<T>& C, const T& R) {
    return complex<T>(C.Re - R,C.Im);
  };
  
  /** Substraction operator with a real value. */
  friend complex<T> operator -(const T& R, const complex<T>& C) {
    return complex<T>(R - C.Re,-C.Im);
  };
  
  /** Multiplication operator. */
  friend complex<T> operator *(const complex<T>& C1, const complex<T>& C2) {
    return complex<T>(C1.Re * C2.Re - C1.Im * C2.Im, C1.Re * C2.Im + C1.Im * C2.Re);
  };
  
  /** Multiplication operator with a real value. */
  friend complex<T> operator *(const complex<T>& C, const T& R) {
    return complex<T>(C.Re * R,C.Im * R);
  };
  
  /** Division operator. */
  friend complex<T> operator /(const complex<T>& C1, const complex<T>& C2) {
    T sqr_mag_inv = 1.0 / (C2.Re * C2.Re + C2.Im * C2.Im);
    return complex<T>(sqr_mag_inv * (C1.Re * C2.Re + C1.Im * C2.Im), sqr_mag_inv * (C1.Im * C2.Re - C1.Re * C2.Im));
  };
  
  /** Division operator with a real value. */
  friend complex<T> operator /(const complex<T>& C, const T& R) {
    return complex<T>(C.Re / R, C.Im / R);
  };
  
  /** Division operator for a real and complex. */
  friend complex<T> operator /(const double& R, const complex<T>& C) {
    T sqr_mag_inv = R / (C.Re * C.Re + C.Im * C.Im);
    return complex<T>(sqr_mag_inv * C.Re, -sqr_mag_inv * C.Im);
  };
  
/*******************************************************************************
                           Comparison Operators
*******************************************************************************/
  /** Equality-comparison operator. */
  friend bool operator ==(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re == C2.Re) && (C1.Im == C2.Im));
  };
  
  /** Equality-comparison operator with a real value. */
  friend bool operator ==(const complex<T>& C, const T& R) {
    return ((C.Re == R) && (C.Im == 0.0));
  };
  
  /** Equality-comparison operator with a real value. */
  friend bool operator ==(const T& R, const complex<T>& C) {
    return ((C.Re == R) && (C.Im == 0.0));
  };
  
  /** Inequality-comparison operator. */
  friend bool operator !=(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re != C2.Re) || (C1.Im != C2.Im));
  };
  
  /** Inequality-comparison operator with a real value. */
  friend bool operator !=(const complex<T>& C, const T& R) {
    return ((C.Re != R) || (C.Im != 0.0));
  };
    
  /** Inequality-comparison operator with a real value. */
  friend bool operator !=(const T& R, const complex<T>& C) {
    return ((C.Re != R) || (C.Im != 0.0));
  };
  
  /** Less-than-comparison operator, magnitude-wise. */
  friend bool operator <(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re * C1.Re + C1.Im * C1.Im) < (C2.Re * C2.Re + C2.Im * C2.Im));
  };
  
  /** Less-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator <(const complex<T>& C, const T& R) {
    return ((C.Re * C.Re + C.Im * C.Im) < R*R);
  };
  
  /** Less-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator <(const T& R, const complex<T>& C) {
    return (R*R < (C.Re * C.Re + C.Im * C.Im));
  };
  
  /** Less-or-equal-comparison operator, magnitude-wise. */
  friend bool operator <=(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re * C1.Re + C1.Im * C1.Im) <= (C2.Re * C2.Re + C2.Im * C2.Im));
  };
  
  /** Less-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator <=(const complex<T>& C, const T& R) {
    return ((C.Re * C.Re + C.Im * C.Im) <= R*R);
  };
  
  /** Less-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator <=(const T& R, const complex<T>& C) {
    return (R*R <= (C.Re * C.Re + C.Im * C.Im));
  };
  
  /** Greater-than-comparison operator, magnitude-wise. */
  friend bool operator >(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re * C1.Re + C1.Im * C1.Im) > (C2.Re * C2.Re + C2.Im * C2.Im));
  };
  
  /** Greater-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator >(const complex<T>& C, const T& R) {
    return ((C.Re * C.Re + C.Im * C.Im) > R*R);
  };
  
  /** Greater-than-comparison operator, magnitude-wise, with a real value. */
  friend bool operator >(const T& R, const complex<T>& C) {
    return (R*R > (C.Re * C.Re + C.Im * C.Im));
  };
  
  /** Greater-or-equal-comparison operator, magnitude-wise. */
  friend bool operator >=(const complex<T>& C1, const complex<T>& C2) {
    return ((C1.Re * C1.Re + C1.Im * C1.Im) >= (C2.Re * C2.Re + C2.Im * C2.Im));
  };
  
  /** Greater-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator >=(const complex<T>& C, const T& R) {
    return ((C.Re * C.Re + C.Im * C.Im) >= R*R);
  };
  
  /** Greater-or-equal-comparison operator, magnitude-wise, with a real value. */
  friend bool operator >=(const T& R, const complex<T>& C) {
    return (R*R >= (C.Re * C.Re + C.Im * C.Im));
  };
  
  
  
  
  /** Complex conjugate for a complex value. */
  friend complex<T> conj(const complex<T>& x) {
    return complex<T>(x.Re,-x.Im);
  };
  
  /** Complex square magnitude for a complex value. */
  friend T sqr_mag(const complex<T>& x) {
    return x.Re*x.Re + x.Im*x.Im;
  };
  
  
  
  


//Exponential and logarithmic functions:

  /** Compute exponential function (function), for a complex value. */
  friend complex<T> exp(const complex<T>& x) {
    if(x.Re == 0.0)
      return complex<T>(cos(x.Im),sin(x.Im));
    else 
      return exp(x.Re) * complex<T>(cos(x.Im),sin(x.Im));
  };  

  /** Compute natural logarithm (function), for a complex value. */
  friend complex<T> log(const complex<T>& x) {
    return complex<T>(log(sqrt(x.Re * x.Re + x.Im * x.Im)),atan2(x.Im,x.Re));
  };  

  /** Compute common logarithm (function), for a complex value. */
  friend complex<T> log10(const complex<T>& x) {
    return log(x) / log(T(10.0));
  };

//Trigonometric functions:

  /** Compute cosine (function), for a complex value.*/
  friend complex<T> cos(const complex<T>& x) {
    return complex<T>(cos(x.Re) * cosh(x.Im),-sin(x.Re) * sinh(x.Im)); 
  }; 

  /** Compute sine (function), for a complex value.*/
  friend complex<T> sin(const complex<T>& x) {
    return complex<T>(sin(x.Re) * cosh(x.Im),cos(x.Re) * sinh(x.Im));
  };  

  /** Compute tangent (function), for a complex value.*/
  friend complex<T> tan(const complex<T>& x) {
  T tmp = 1.0/(cos(2.0 * x.Re) + cosh(2.0 * x.Im));
    return complex<T>(sin(2.0*x.Re) * tmp, sinh(2.0*x.Im) * tmp); 
  };  

  /** Compute arc cosine (function), for a complex value.*/
  friend complex<T> acos(const complex<T>& x) {
    return  complex<T>(0.0,-1.0) * log(x + sqrt(x * x - 1.0));
  };  

  /** Compute arc sine (function), for a complex value.*/
  friend complex<T> asin(const complex<T>& x) {
    return complex<T>(0.0,-1.0) * log(complex<T>(-x.Im,x.Re) + sqrt(1.0 - x * x));
  };  

  /** Compute arc tangent (function), for a complex value.*/
  friend complex<T> atan(const complex<T>& x) {
    return complex<T>(0.0,0.5) * log(complex<T>(x.Im + 1.0,-x.Re)/complex<T>(1.0-x.Im,x.Re));
  };  

  /** Compute arc tangent with two parameters (function), for a complex value.*/
  friend complex<T> atan2(const complex<T>& y, const complex<T>& x) {
    return atan(y / x);
  };  

//Hyperbolic functions:

  /** Compute hyperbolic cosine (function), for a complex value.*/
  friend complex<T> cosh(const complex<T>& x) {
    return complex<T>(cosh(x.Re) * cos(x.Im),sinh(x.Re) * sin(x.Im));
  };  

  /** Compute hyperbolic sine (function), for a complex value.*/
  friend complex<T> sinh(const complex<T>& x) {
    return complex<T>(cos(x.Im) * sinh(x.Re),sin(x.Im) * cosh(x.Re)); 
  };  

  /** Compute hyperbolic tangent (function), for a complex value.*/
  friend complex<T> tanh(const complex<T>& x) {
    T tmp = 1.0/(cosh(2.0 * x.Re) + cos(2.0 * x.Im));
    return complex<T>(sinh(2.0*x.Re) * tmp, sin(2.0*x.Im) * tmp);
  };  


//Power functions

  /** Raise to power (function), for a complex value.*/
  friend complex<T> pow(const complex<T>& base, const complex<T>& exponent) {
    return exp(exponent * log(base));
  };  

  /** Compute square root (function), for a complex value.*/
  friend complex<T> sqrt(const complex<T>& x) {
    T angle = atan2(x.Im,x.Re);
    return pow(x.Re * x.Re + x.Im * x.Im,0.25) * complex<T>(cos(0.5 * angle),sin(0.5 * angle));
  };  


//Rounding, absolute value and remainder functions:

  /** Round up value (function), for a complex value. */
  friend complex<T> ceil(const complex<T>& x) {
    return complex<T>(ceil(x.Re),ceil(x.Im));
  };  

  /** Compute absolute value (function), for a complex value. */
  friend T fabs(const complex<T>& x) {
    return sqrt(x.Re * x.Re + x.Im * x.Im);
  };  

  /** Round down value (function), for a complex value.*/
  friend complex<T> floor(const complex<T>& x) {
    return complex<T>(floor(x.Re),floor(x.Im));
  };  
  
  
  
  
  
  
  
  /// Loading a complex value.
  friend serialization::iarchive& RK_CALL operator >>(serialization::iarchive& in, complex<T>& C) {
    return in >> C.Re >> C.Im;
  };

  /// Loading a complex value with a name.
  friend serialization::iarchive& RK_CALL operator &(serialization::iarchive& in, const std::pair<std::string, complex<T>& >& f) {
    return in & std::pair<std::string, T& >(f.first + "_Re",f.second.Re)
              & std::pair<std::string, T& >(f.first + "_Im",f.second.Im);
  };
  
  /// Saving a complex value.
  friend serialization::oarchive& RK_CALL operator <<(serialization::oarchive& out, const complex<T>& C) {
    return out << C.Re << C.Im;
  };
  
  /// Saving a complex value with a name.
  friend serialization::oarchive& RK_CALL operator &(serialization::oarchive& out, const std::pair<std::string, const complex<T>& >& C) {
    return out & std::pair<std::string, T >(C.first + "_Re",C.second.Re)
               & std::pair<std::string, T >(C.first + "_Im",C.second.Im);
  };
  
};


namespace rtti {

template <typename T>
struct get_type_id< complex<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000007);
  static std::string type_name() { return "ReaK::complex"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const complex<T>& save_type;
  typedef complex<T>& load_type;
};

};


/** Complex conjugate for a real value. */
inline float conj(float x) { return x; };
/** Complex conjugate for a real value. */
inline double conj(double x) { return x; };
  
/** Complex square magnitude for a real value. */
inline float sqr_mag(float x) { return x*x; };
/** Complex square magnitude for a real value. */
inline double sqr_mag(double x) { return x*x; };

  
};




#endif //RK_COMPLEX_MATH_HPP









