
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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_are_solver.hpp>
#include <ReaK/math/lin_alg/mat_norms.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/*
 * This program is just a little test program that is used temporarily to test little
 * bits of code here and there related to the lin-alg library. In other words, if I
 * don't know if a particular expression is going to compile and compute correctly,
 * I can just test it in this "sand-box". Useful for user-space testing as well as
 * debugging of new additions to the library.
 */

int main() {

  using namespace ReaK;

  {
    /*
    mat<double,mat_structure::rectangular> A(4,4);
    mat<double,mat_structure::rectangular> B(4,2);
    A(0,0) = 0.0; A(0,1) = 0.0; A(0,2) = 1.0; A(0,3) = 0.0;
    A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.0; A(1,3) = 1.0;
    A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
    A(3,0) = 0.0; A(3,1) = 0.0; A(3,2) = 0.0; A(3,3) = 0.0;

    B(0,0) = 0.0; B(0,1) = 0.0;
    B(1,0) = 0.0; B(1,1) = 0.0;
    B(2,0) = 0.01; B(2,1) = 0.0;
    B(3,0) = 0.0; B(3,1) = 1.0;

    mat<double,mat_structure::rectangular> R =
    mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(2));
    R(0,0) = 0.01;
    mat<double,mat_structure::rectangular> Q =
    mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(4));
    Q(0,0) = 0.01;
    Q(2,2) = 0.01;
    */

    mat<double, mat_structure::rectangular> A(4, 4);
    mat<double, mat_structure::rectangular> B(4, 1);
    A(0, 0) = 0.0;
    A(0, 1) = 0.0;
    A(0, 2) = 1.0;
    A(0, 3) = 0.0;
    A(1, 0) = 0.0;
    A(1, 1) = 0.0;
    A(1, 2) = 0.0;
    A(1, 3) = 1.0;
    A(2, 0) = 0.0;
    A(2, 1) = 0.0;
    A(2, 2) = 0.0;
    A(2, 3) = 0.0;
    A(3, 0) = 0.0;
    A(3, 1) = 0.0;
    A(3, 2) = 0.0;
    A(3, 3) = 0.0;

    B(0, 0) = 0.0;
    B(1, 0) = 0.0;
    B(2, 0) = 0.0;
    B(3, 0) = 1.0;

    mat<double, mat_structure::rectangular> R =
        mat<double, mat_structure::rectangular>(
            mat<double, mat_structure::identity>(1));
    mat<double, mat_structure::rectangular> Q =
        mat<double, mat_structure::rectangular>(
            mat<double, mat_structure::identity>(4));
    Q(0, 0) = 0.0;
    Q(2, 2) = 0.0;

    std::cout << " A = " << A << std::endl;
    std::cout << " B = " << B << std::endl;
    mat<double, mat_structure::rectangular> Br(B);
    mat<double, mat_structure::square> Tr(
        mat<double, mat_structure::identity>(B.get_row_count()));
    mat<double, mat_structure::permutation> PCr(B.get_row_count());
    std::size_t r = detail::decompose_RRQR_impl(Br, &Tr, PCr, 1e-4);
    mat<double, mat_structure::rectangular> Ar = transpose_view(Tr) * A * Tr;
    std::cout << " Ar = " << Ar << std::endl;
    std::cout << " Br = " << Br << std::endl;
    std::cout << " Tr = " << Tr << std::endl;
    std::cout << " PCr = " << PCr << std::endl;
    // look for the biggest zero-block at the lower left corner of Ar, cannot be bigger than N-r.
    std::size_t N = Ar.get_row_count();
    for (std::size_t i = N; i > r;) {
      --i;
      for (std::size_t j = 0; j < N - r; ++j)
        if (std::abs(Ar(i, j)) > 1e-6)
          r = N - j;
    };
    // now, r contains the number of controllable states.
    std::cout << " r = " << r << std::endl;
    mat<double, mat_structure::square> Rr(transpose(PCr) * R * PCr);
    mat<double, mat_structure::rectangular> Trr(
        sub(Tr)(range(0, Tr.get_row_count()), range(0, r)));
    mat<double, mat_structure::square> Qr(transpose_view(Trr) * Q * Trr);
    std::cout << " Trr = " << Trr << std::endl;
    std::cout << " Rr = " << Rr << std::endl;
    std::cout << " Qr = " << Qr << std::endl;

    mat<double, mat_structure::rectangular> Pr(r, r);
    mat<double, mat_structure::rectangular> Kr(B.get_col_count(), r);
    solve_IHCT_LQR(sub(Ar)(range(0, r), range(0, r)),
                   sub(Br)(range(0, r), range(0, Br.get_col_count())), Qr, Rr,
                   Kr, Pr, 1e-6, false);

    mat<double, mat_structure::rectangular> P = Trr * Pr * transpose_view(Trr);
    mat<double, mat_structure::rectangular> K = PCr * Kr * transpose_view(Trr);
    std::cout << " Pr = " << Pr << std::endl;
    std::cout << " Kr = " << Kr << std::endl;
    std::cout << " P = " << P << std::endl;
    std::cout << " K = " << K << std::endl;

#if 0
  mat<double,mat_structure::rectangular> P(A.get_row_count(),A.get_col_count());
  mat<double,mat_structure::rectangular> K(B.get_col_count(),A.get_col_count());
  solve_IHCT_LQR(A, B, Q, R, K, P, 1e-4, false);
#endif

    mat<double, mat_structure::rectangular> M_tmp = R;
    mat<double, mat_structure::rectangular> M2_tmp(B.get_col_count(),
                                                   A.get_col_count());
    M2_tmp = transpose_view(B) * P;
    mat<double, mat_structure::rectangular> Msol_tmp;
    linlsq_QR(M_tmp, Msol_tmp, M2_tmp);
    mat<double, mat_structure::rectangular> X =
        (transpose_view(A) * P + P * A - P * B * Msol_tmp + Q);
    std::cout << "Q + A^T P + P A - P B R^{-1} B^T P = " << X << std::endl;
    double err_norm = norm_1(X);
    double P_norm = norm_1(P);
    std::cout << "\t|| Err || = " << std::setw(14) << err_norm
              << " relative to || P || = " << std::setw(14) << P_norm
              << " for a relative error of " << std::setw(14)
              << (err_norm / P_norm) << std::endl;
  };

  return 0;
};

// These were just some tests for debugging the Schur decomposition code.
#if 0
  {
  mat<double,mat_structure::rectangular> A(6,6);
  mat<double,mat_structure::rectangular> B(6,6);
  A(0,0) = 50.0; A(0,1) = -60.0; A(0,2) = 50.0; A(0,3) = -27.0; A(0,4) = 6.0; A(0,5) = 6.0;
  A(1,0) = 38.0; A(1,1) = -28.0; A(1,2) = 27.0; A(1,3) = -17.0; A(1,4) = 5.0; A(1,5) = 5.0;
  A(2,0) = 27.0; A(2,1) = -17.0; A(2,2) = 27.0; A(2,3) = -17.0; A(2,4) = 5.0; A(2,5) = 5.0;
  A(3,0) = 27.0; A(3,1) = -28.0; A(3,2) = 38.0; A(3,3) = -17.0; A(3,4) = 5.0; A(3,5) = 5.0;
  A(4,0) = 27.0; A(4,1) = -28.0; A(4,2) = 27.0; A(4,3) = -17.0; A(4,4) = 16.0; A(4,5) = 5.0;
  A(5,0) = 27.0; A(5,1) = -28.0; A(5,2) = 27.0; A(5,3) = -17.0; A(5,4) = 5.0; A(5,5) = 16.0;
  
  B(0,0) = 16.0; B(0,1) = 5.0; B(0,2) = 5.0; B(0,3) = 5.0; B(0,4) = -6.0; B(0,5) = 5.0;
  B(1,0) = 5.0; B(1,1) = 16.0; B(1,2) = 5.0; B(1,3) = 5.0; B(1,4) = -6.0; B(1,5) = 5.0;
  B(2,0) = 5.0; B(2,1) = 5.0; B(2,2) = 16.0; B(2,3) = 5.0; B(2,4) = -6.0; B(2,5) = 5.0;
  B(3,0) = 5.0; B(3,1) = 5.0; B(3,2) = 5.0; B(3,3) = 16.0; B(3,4) = -6.0; B(3,5) = 5.0;
  B(4,0) = 5.0; B(4,1) = 5.0; B(4,2) = 5.0; B(4,3) = 5.0; B(4,4) = -6.0; B(4,5) = 16.0;
  B(5,0) = 6.0; B(5,1) = 6.0; B(5,2) = 6.0; B(5,3) = 6.0; B(5,4) = -5.0; B(5,5) = 6.0;
  
  mat<double,mat_structure::rectangular> T(6,6);
  mat<double,mat_structure::rectangular> R(6,6);
  mat<double,mat_structure::rectangular> Q = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(6));
  mat<double,mat_structure::rectangular> Z = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(6));
  decompose_GenRealSchur(A,B,Q,Z,T,R);
  std::cout << "T = " << T << std::endl;
  std::cout << "R = " << R << std::endl;
  std::cout << "Q = " << Q << std::endl;
  std::cout << "Z = " << Z << std::endl;
  std::cout << "QTZ' = " << (Q * T * transpose_view(Z)) << std::endl;
  std::cout << "QRZ' = " << (Q * R * transpose_view(Z)) << std::endl;
  
  };
  
  {
  mat<double,mat_structure::rectangular> A(3,3);
  mat<double,mat_structure::rectangular> B(3,3);
  A(0,0) = 2.0; A(0,1) = 3.0; A(0,2) = 5.0; 
  A(1,0) = 8.0; A(1,1) = -2.0; A(1,2) = -6.0; 
  A(2,0) = 7.0; A(2,1) = -10.0; A(2,2) = 3.0; 
  
  B(0,0) = 6.0; B(0,1) = 3.0; B(0,2) = 2.0; 
  B(1,0) = 5.0; B(1,1) = -3.0; B(1,2) = 2.0; 
  B(2,0) = 2.0; B(2,1) = 1.0; B(2,2) = 1.0; 
  
  mat<double,mat_structure::rectangular> T(3,3);
  mat<double,mat_structure::rectangular> R(3,3);
  mat<double,mat_structure::rectangular> Q = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(3));
  mat<double,mat_structure::rectangular> Z = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(3));
  decompose_GenRealSchur(A,B,Q,Z,T,R);
  std::cout << "T = " << T << std::endl;
  std::cout << "R = " << R << std::endl;
  std::cout << "Q = " << Q << std::endl;
  std::cout << "Z = " << Z << std::endl;
  std::cout << "QTZ' = " << (Q * T * transpose_view(Z)) << std::endl;
  std::cout << "QRZ' = " << (Q * R * transpose_view(Z)) << std::endl;
  
  };
  
  {
  mat<double,mat_structure::rectangular> A(4,4);
  mat<double,mat_structure::rectangular> B(4,4);
  A(0,0) = 2.0; A(0,1) = 3.0; A(0,2) = 5.0; A(0,3) = 2.0; 
  A(1,0) = 8.0; A(1,1) = -2.0; A(1,2) = -6.0; A(1,3) = 3.0; 
  A(2,0) = 7.0; A(2,1) = -10.0; A(2,2) = 3.0; A(2,3) = 5.0; 
  A(3,0) = 7.0; A(3,1) = -10.0; A(3,2) = 3.0; A(3,3) = 3.0; 
  
  B(0,0) = 6.0; B(0,1) = 3.0; B(0,2) = 2.0; B(0,3) = 3.0;
  B(1,0) = 5.0; B(1,1) = -3.0; B(1,2) = 2.0; B(1,3) = 6.0;
  B(2,0) = 2.0; B(2,1) = 1.0; B(2,2) = 1.0; B(2,3) = 1.0;
  B(3,0) = 2.0; B(3,1) = 1.0; B(3,2) = 1.0; B(3,3) = 2.0;
  
  mat<double,mat_structure::rectangular> T(4,4);
  mat<double,mat_structure::rectangular> R(4,4);
  mat<double,mat_structure::rectangular> Q = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(4));
  mat<double,mat_structure::rectangular> Z = mat<double,mat_structure::rectangular>(mat<double,mat_structure::identity>(4));
  decompose_GenRealSchur(A,B,Q,Z,T,R);
  std::cout << "T = " << T << std::endl;
  std::cout << "R = " << R << std::endl;
  std::cout << "Q = " << Q << std::endl;
  std::cout << "Z = " << Z << std::endl;
  std::cout << "QTZ' = " << (Q * T * transpose_view(Z)) << std::endl;
  std::cout << "QRZ' = " << (Q * R * transpose_view(Z)) << std::endl;
  
  };
  
  
  
  {
  std::cout << "****** Problem 1 ****** " << std::endl;
  mat<double,mat_structure::rectangular> F(2,2);
  mat<double,mat_structure::rectangular> G(2,1);
  F(0,0) = 2.0; F(0,1) = -1.0; 
  F(1,0) = 1.0; F(1,1) = 0.0; 
  
  G(0,0) = 1.0; 
  G(1,0) = 0.0; 
  
  mat<double,mat_structure::rectangular> R(1,1);
  R(0,0) = 0.0;
  mat<double,mat_structure::rectangular> Q(2,2);
  Q(0,0) = 0.0; Q(0,1) = 0.0;
  Q(1,0) = 0.0; Q(1,1) = 1.0;
  
  mat<double,mat_structure::rectangular> P(2,2);
  solve_dare_problem(F,G,Q,R,P);
  std::cout << "P = " << P << std::endl;
  
  
  mat<double,mat_structure::rectangular> M_tmp(1,1);
  M_tmp = ( R + transpose_view(G) * P * G );
  mat<double,mat_structure::rectangular> M2_tmp(1,2);
  M2_tmp = transpose_view(G) * P * F;
  mat<double,mat_structure::rectangular> Msol_tmp(1,2);
  linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
  std::cout << "F' P F - F' P G ( R + G' P G )^-1 G' P F + Q = " 
            << (transpose_view(F) * P * F - transpose_view(F) * P * G * Msol_tmp + Q) << std::endl;
  
  };
  
  {
  
  std::cout << "****** Problem 2 ****** " << std::endl;
  mat<double,mat_structure::rectangular> F(2,2);
  mat<double,mat_structure::rectangular> G(2,1);
  F(0,0) = 1.0; F(0,1) = 0.2; 
  F(1,0) = 0.0; F(1,1) = 1.0; 
  
  G(0,0) = 0.0; 
  G(1,0) = 1.0; 
  
  mat<double,mat_structure::rectangular> R(1,1);
  R(0,0) = 1.0;
  mat<double,mat_structure::rectangular> Q(2,2);
  Q(0,0) = 1.0; Q(0,1) = 0.0;
  Q(1,0) = 0.0; Q(1,1) = 1.0;
  
  mat<double,mat_structure::rectangular> P(2,2);
  solve_dare_problem(F,G,Q,R,P,1e-4);
  std::cout << "P = " << P << std::endl;
  
  
  mat<double,mat_structure::rectangular> M_tmp(1,1);
  M_tmp = ( R + transpose_view(G) * P * G );
  mat<double,mat_structure::rectangular> M2_tmp(1,2);
  M2_tmp = transpose_view(G) * P * F;
  mat<double,mat_structure::rectangular> Msol_tmp(1,2);
  linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
  std::cout << "F' P F - F' P G ( R + G' P G )^-1 G' P F + Q = " 
            << (transpose_view(F) * P * F - transpose_view(F) * P * G * Msol_tmp + Q) << std::endl;
  };

#endif
