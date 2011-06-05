
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

#ifndef MAT_EXP_METHODS_HPP
#define MAT_EXP_METHODS_HPP

#include "mat_traits.hpp"
#include "mat_concepts.hpp"

#include "mat_norms.hpp"
#include "mat_alg.hpp"


namespace ReaK {


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
  
  mat<ValueType,mat_structure::square> A_tmp = A - (mu * mat<ValueType,mat_structure::identity>(N));
  ValueType n1 = norm_1(A_tmp);
    
  if(n1 <= 1.495585217958292e-2) { //m = 3
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> U = A_tmp * (A2
                                                    + 60 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = 12 * A2
                                           + 120 * mat<ValueType,mat_structure::identity>(N);
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 2.539398330063230e-1) { //m = 5
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> U = A_tmp * (A4
                                                    + 420 * A2
                                                    + 15120 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = 30 * A4
                                           + 3360 * A2
                                           + 30240 * mat<ValueType,mat_structure::identity>(N);
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 9.504178996162932e-1) { //m = 7
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * (A6
                                                    + 1512 * A4
                                                    + 277200 * A2
                                                    + 8648640 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = 56 * A6
                                           + 25200 * A4
                                           + 1995840 * A2
                                           + 17297280 * mat<ValueType,mat_structure::identity>(N);
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 2.097847961257068e0) { //m = 9
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> A8 = A4 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * (A8
                                                    + 3960 * A6
                                                    + 2162160 * A4
                                                    + 302702400 * A2
                                                    + 8821612800 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = 90 * A8
                                           + 110880 * A6 
                                           + 30270240 * A4
                                           + 2075673600 * A2
                                           + 17643225600 * mat<ValueType,mat_structure::identity>(N);
    linsolve(V - U, X, U + V, NumTol);
  } else if(n1 <= 5.371920351148152e0) { //m = 13
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * ( A6 * (A6
                                                           + 16380 * A4 
                                                           + 40840800 * A2) 
                                                     + 33522128640 * A6
                                                     + 10559470521600 * A4 
                                                     + 1187353796428800 * A2 
                                                     + 32382376266240000 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = A6 * (182 * A6 
                                                 + 960960 * A4 
                                                 + 1323241920 * A2) 
                                           + 670442572800 * A6 
                                           + 129060195264000 * A4 
                                           + 7771770303897600 * A2 
                                           + 64764752532480000 * mat<ValueType,mat_structure::identity>(N);
    linsolve(V - U, X, U + V, NumTol);
  } else {
    SizeType s = 1;
    while(5.371920351148152e0 * (1 << s) < n1) ++s;
    A_tmp *= 1.0 / ValueType(1 << s);
    mat<ValueType,mat_structure::square> A2 = A_tmp * A_tmp;
    mat<ValueType,mat_structure::square> A4 = A2 * A2;
    mat<ValueType,mat_structure::square> A6 = A2 * A4;
    mat<ValueType,mat_structure::square> U = A_tmp * ( A6 * (A6 
                                                           + 16380 * A4 
                                                           + 40840800 * A2) 
                                                     + 33522128640 * A6
                                                     + 10559470521600 * A4 
                                                     + 1187353796428800 * A2 
                                                     + 32382376266240000 * mat<ValueType,mat_structure::identity>(N));
    mat<ValueType,mat_structure::square> V = A6 * (182 * A6 
                                                 + 960960 * A4 
                                                 + 1323241920 * A2) 
                                           + 670442572800 * A6 
                                           + 129060195264000 * A4 
                                           + 7771770303897600 * A2 
                                           + 64764752532480000 * mat<ValueType,mat_structure::identity>(N);
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





};


#endif








