/**
 * \file mat_exp_methods.hpp
 * 
 * This library provides a function to approximate a matrix exponential using Pade approximant (classic Pade 
 * Square-and-Sum method). This function uses the help of a linear equation solver, which is a functor that 
 * wraps any method to solve a linear system of equations involving square matrices (this is a customization 
 * point for the user, where any of the many linear equation solvers can be used, whichever is deemed 
 * appropriate).
 * 
 * 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_MAT_EXP_METHODS_HPP
#define REAK_MAT_EXP_METHODS_HPP

#include "mat_traits.hpp"
#include "mat_concepts.hpp"

#include "mat_norms.hpp"
#include "mat_alg.hpp"


namespace ReaK {

/**
 * This function approximates a matrix exponential using Pade approximant (classic Pade 
 * Square-and-Sum method). This function uses the help of a linear equation solver, which is a functor that 
 * wraps any method to solve a linear system of equations involving square matrices (this is a customization 
 * point for the user, where any of the many linear equation solvers can be used, whichever is deemed 
 * appropriate).
 * \tparam Matrix1 Should be a square, readable matrix type.
 * \tparam Matrix2 Should be a fully-writable matrix type (which can be square).
 * \tparam LinearEqSolver Should be a functor that can be called to solve a linear system of equation (see examples linked to existing linear equation solvers).
 * \param A The matrix that is the operand to the exponential.
 * \param X The matrix that will store, as output, the result.
 * \param linsolve The functor that can solve a linear system of equations.
 * \param NumTol The tolerance at which a value is considered zero, used for singularity detection.
 * 
 * \throw singularity_error if any of the linear systems involved in the computation turns out to involve a singular matrix, 
 *                          meaning that the problem is ill-conditioned.
 * 
 */
template <typename Matrix1, typename Matrix2, typename LinearEqSolver>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type exp_PadeSAS(const Matrix1& A, Matrix2& X, LinearEqSolver linsolve, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::exp;
  
  SizeType N = A.get_row_count();
  ValueType mu = trace(A) / N;
  mat<ValueType,mat_structure::identity> eye_N = mat<ValueType,mat_structure::identity>(N);
  
  mat<ValueType,mat_structure::square> A_tmp = A - (mu * eye_N);
  ValueType n1 = norm_1(A_tmp);
    
  if(n1 <= 1.495585217958292e-2) { //m = 3
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> U = A_tmp * (A2
                                                    + ValueType(60) * eye_N);
    mat<ValueType,mat_structure::square> V = ValueType(12) * A2
                                           + ValueType(120) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 2.539398330063230e-1) { //m = 5
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> U = A_tmp * (A4
                                                    + ValueType(420) * A2
                                                    + ValueType(15120) * eye_N);
    mat<ValueType,mat_structure::square> V = ValueType(30) * A4
                                           + ValueType(3360) * A2
                                           + ValueType(30240) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 9.504178996162932e-1) { //m = 7
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * (A6
                                                    + ValueType(1512) * A4
                                                    + ValueType(277200) * A2
                                                    + ValueType(8648640) * eye_N);
    mat<ValueType,mat_structure::square> V = ValueType(56) * A6
                                           + ValueType(25200) * A4
                                           + ValueType(1995840) * A2
                                           + ValueType(17297280) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 2.097847961257068e0) { //m = 9
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> A8 = A4 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * (A8
                                                    + ValueType(3960) * A6
                                                    + ValueType(2162160) * A4
                                                    + ValueType(302702400) * A2
                                                    + ValueType(8821612800) * eye_N);
    mat<ValueType,mat_structure::square> V = ValueType(90) * A8
                                           + ValueType(110880) * A6 
                                           + ValueType(30270240) * A4
                                           + ValueType(2075673600) * A2
                                           + ValueType(17643225600) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 5.371920351148152e0) { //m = 13
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * ( A6 * (A6
                                                           + ValueType(16380) * A4 
                                                           + ValueType(40840800) * A2) 
                                                     + ValueType(33522128640) * A6
                                                     + ValueType(10559470521600) * A4 
                                                     + ValueType(1187353796428800) * A2 
                                                     + ValueType(32382376266240000) * eye_N);
    mat<ValueType,mat_structure::square> V = A6 * (ValueType(182) * A6 
                                                 + ValueType(960960) * A4 
                                                 + ValueType(1323241920) * A2) 
                                           + ValueType(670442572800) * A6 
                                           + ValueType(129060195264000) * A4 
                                           + ValueType(7771770303897600) * A2 
                                           + ValueType(64764752532480000) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
  } else {
    SizeType s = 1;
    while(5.371920351148152e0 * (1 << s) < n1) ++s;
    A_tmp *= 1.0 / ValueType(1 << s);
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * ( A6 * (A6 
                                                           + ValueType(16380) * A4 
                                                           + ValueType(40840800) * A2) 
                                                     + ValueType(33522128640) * A6
                                                     + ValueType(10559470521600) * A4 
                                                     + ValueType(1187353796428800) * A2 
                                                     + ValueType(32382376266240000) * eye_N);
    mat<ValueType,mat_structure::square> V = A6 * (ValueType(182) * A6 
                                                 + ValueType(960960) * A4 
                                                 + ValueType(1323241920) * A2) 
                                           + ValueType(670442572800) * A6 
                                           + ValueType(129060195264000) * A4 
                                           + ValueType(7771770303897600) * A2 
                                           + ValueType(64764752532480000) * eye_N;
    linsolve(V - U, X, U + V, NumTol);
    for(;s != 0; --s)
      X *= X;
  };
  X *= exp(mu);
  
//   for m = [3 5 7 9 13]
//     if A 1 ≤ θm
//       X = rm (A) % rm (A) = [m/m] Pad´ approximant to A.
//       X = eμ X % Undo preprocessing.
//     end
//   end
//   A ← A/2^s with s a minimal integer such that A/2s 1 ≤ θ13 (i.e., s = log2 ( A 1 /θ13 ) ).
//   
//   % Form [13/13] Pad´ approximant to eA .
//   A2 = A^2 , A4 = A2^2, A6 = A2 * A4
//   U = A A6 (b13 A6 + b11 A4 + b9 A2 ) + b7 A6 + b5 A4 + b3 A2 + b1 I)
//   V = A6 (b12 A6 + b10 A4 + b8 A2 ) + b6 A6 + b4 A4 + b2 A2 + b0 I
//   Solve (−U + V )r13 = U + V for r13 .
//   X = (r13^2)^2 .. s-times  by repeated squaring.
//   X = eμ X % Undo preprocessing.

};

// m	Theta_m 
// 3	1.495585217958292e-2
// 5	2.539398330063230e-1
// 7	9.504178996162932e-1
// 9	2.097847961257068e0
// 13	5.371920351148152e0

// b(0: 13) = [64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1]

// b_j = (2*m - j)! / (m - j)! / j!



#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template void exp_PadeSAS(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& X, QR_linlsqsolver linsolve, double NumTol);
extern template void exp_PadeSAS(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& X, RRQR_linlsqsolver linsolve, double NumTol);


extern template void exp_PadeSAS(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& X, QR_linlsqsolver linsolve, float NumTol);
extern template void exp_PadeSAS(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& X, RRQR_linlsqsolver linsolve, float NumTol);

#endif




};


#endif








