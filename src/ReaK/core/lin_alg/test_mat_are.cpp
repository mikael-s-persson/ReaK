
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

#include "base/defs.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

#include "mat_alg.hpp"

#include "mat_are_solver.hpp"

#include "mat_norms.hpp"

int main() {

  using namespace ReaK;
  
  
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
  
  
  
  std::size_t passed = 0;
  std::size_t possible_passes = 17;
  
  {
  std::vector< mat<double, mat_structure::rectangular> > F_list(9);
  std::vector< mat<double, mat_structure::rectangular> > G_list(9);
  std::vector< mat<double, mat_structure::rectangular> > Q_list(9);
  std::vector< mat<double, mat_structure::rectangular> > R_list(9);
  
  F_list[0] = mat<double,mat_structure::rectangular>(2,2);
  F_list[0](0,0) =  4.0; F_list[0](0,1) = 3.0; 
  F_list[0](1,0) = -4.5; F_list[0](1,1) = -3.5; 
  
  G_list[0] = mat<double,mat_structure::rectangular>(2,1);
  G_list[0](0,0) =  1.0; 
  G_list[0](1,0) = -1.0; 
  
  R_list[0] = mat<double,mat_structure::rectangular>(1,1);
  R_list[0](0,0) = 1.0; 
  
  Q_list[0] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[0](0,0) = 9.0; Q_list[0](0,1) = 6.0; 
  Q_list[0](1,0) = 6.0; Q_list[0](1,1) = 4.0; 
  
  
  F_list[1] = mat<double,mat_structure::rectangular>(2,2);
  F_list[1](0,0) = 0.9512; F_list[1](0,1) = 0.0; 
  F_list[1](1,0) = 0.0; F_list[1](1,1) = 0.9048; 
  
  G_list[1] = mat<double,mat_structure::rectangular>(2,2);
  G_list[1](0,0) = 4.877; G_list[1](0,1) = 4.877; 
  G_list[1](1,0) = -1.1895; G_list[1](1,1) = 3.569; 
  
  R_list[1] = mat<double,mat_structure::rectangular>(2,2);
  R_list[1](0,0) = 1.0 / 3.0; R_list[1](0,1) = 0.0; 
  R_list[1](1,0) = 0.0; R_list[1](1,1) = 3.0; 
  
  Q_list[1] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[1](0,0) = 0.005; Q_list[1](0,1) = 0.0; 
  Q_list[1](1,0) = 0.0; Q_list[1](1,1) = 0.02; 
  
  
  F_list[2] = mat<double,mat_structure::rectangular>(2,2);
  F_list[2](0,0) = 2.0; F_list[2](0,1) = -1.0; 
  F_list[2](1,0) = 1.0; F_list[2](1,1) = 0.0; 
  
  G_list[2] = mat<double,mat_structure::rectangular>(2,1);
  G_list[2](0,0) = 1.0; 
  G_list[2](1,0) = 0.0; 
  
  R_list[2] = mat<double,mat_structure::rectangular>(1,1);
  R_list[2](0,0) = 0.0; 
  
  Q_list[2] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[2](0,0) = 0.0; Q_list[2](0,1) = 0.0; 
  Q_list[2](1,0) = 0.0; Q_list[2](1,1) = 1.0; 
  
  
  F_list[3] = mat<double,mat_structure::rectangular>(2,2);
  F_list[3](0,0) = 0.0; F_list[3](0,1) =  1.0; 
  F_list[3](1,0) = 0.0; F_list[3](1,1) = -1.0; 
  
  G_list[3] = mat<double,mat_structure::rectangular>(2,2);
  G_list[3](0,0) = 1.0; G_list[3](0,1) = 0.0; 
  G_list[3](1,0) = 2.0; G_list[3](1,1) = 1.0; 
  
  R_list[3] = mat<double,mat_structure::rectangular>(2,2);
  R_list[3](0,0) = 9.0; R_list[3](0,1) = 3.0; 
  R_list[3](1,0) = 3.0; R_list[3](1,1) = 1.0; 
  
  Q_list[3] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[3](0,0) = -4.0/11.0; Q_list[3](0,1) = -4.0/11.0; 
  Q_list[3](1,0) = -4.0/11.0; Q_list[3](1,1) =  7.0/11.0; 
  
  
  F_list[4] = mat<double,mat_structure::rectangular>(2,2);
  F_list[4](0,0) = 1.0; F_list[4](0,1) = 0.2; 
  F_list[4](1,0) = 0.0; F_list[4](1,1) = 1.0; 
  
  G_list[4] = mat<double,mat_structure::rectangular>(2,1);
  G_list[4](0,0) = 0.0; 
  G_list[4](1,0) = 1.0; 
  
  R_list[4] = mat<double,mat_structure::rectangular>(1,1);
  R_list[4](0,0) = 1.0; 
  
  Q_list[4] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[4](0,0) = 1.0; Q_list[4](0,1) = 0.0; 
  Q_list[4](1,0) = 0.0; Q_list[4](1,1) = 1.0; 
  
  
  F_list[5] = mat<double,mat_structure::rectangular>(2,2);
  F_list[5](0,0) = 0.0; F_list[5](0,1) = 1.0; 
  F_list[5](1,0) = 0.0; F_list[5](1,1) = 0.0; 
  
  G_list[5] = mat<double,mat_structure::rectangular>(2,1);
  G_list[5](0,0) = 0.0; 
  G_list[5](1,0) = 1.0; 
  
  R_list[5] = mat<double,mat_structure::rectangular>(1,1);
  R_list[5](0,0) = 1.0; 
  
  Q_list[5] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[5](0,0) = 1.0; Q_list[5](0,1) = 2.0; 
  Q_list[5](1,0) = 2.0; Q_list[5](1,1) = 4.0; 
  
  
  F_list[6] = mat<double,mat_structure::rectangular>(4,4);
  F_list[6](0,0) =  0.998; F_list[6](0,1) = 0.067; F_list[6](0,2) = 0.0; F_list[6](0,3) = 0.0; 
  F_list[6](1,0) = -0.067; F_list[6](1,1) = 0.998; F_list[6](1,2) = 0.0; F_list[6](1,3) = 0.0;
  F_list[6](2,0) = 0.0; F_list[6](2,1) = 0.0; F_list[6](2,2) =  0.998; F_list[6](2,3) = 0.153; 
  F_list[6](3,0) = 0.0; F_list[6](3,1) = 0.0; F_list[6](3,2) = -0.153; F_list[6](3,3) = 0.998; 
  
  G_list[6] = mat<double,mat_structure::rectangular>(4,2);
  G_list[6](0,0) = 0.0033; G_list[6](0,1) = 0.02; 
  G_list[6](1,0) = 0.1; G_list[6](1,1) = -0.0007; 
  G_list[6](2,0) = 0.04; G_list[6](2,1) = 0.0073; 
  G_list[6](3,0) = -0.0028; G_list[6](3,1) = 0.1; 
  
  R_list[6] = mat<double,mat_structure::rectangular>(2,2);
  R_list[6](0,0) = 1.0; R_list[6](0,1) = 0.0; 
  R_list[6](1,0) = 0.0; R_list[6](1,1) = 1.0; 
  
  Q_list[6] = mat<double,mat_structure::rectangular>(4,4);
  Q_list[6](0,0) = 1.87; Q_list[6](0,1) = 0.0; Q_list[6](0,2) = 0.0; Q_list[6](0,3) = -0.244; 
  Q_list[6](1,0) = 0.0; Q_list[6](1,1) = 0.744; Q_list[6](1,2) = 0.205; Q_list[6](1,3) = 0.0; 
  Q_list[6](2,0) = 0.0; Q_list[6](2,1) = 0.205; Q_list[6](2,2) = 0.589; Q_list[6](2,3) = 0.0; 
  Q_list[6](3,0) = -0.244; Q_list[6](3,1) = 0.0; Q_list[6](3,2) = 0.0; Q_list[6](3,3) = 1.048; 
  
  
  F_list[7] = mat<double,mat_structure::rectangular>(4,4);
  F_list[7](0,0) =  0.98475; F_list[7](0,1) = -0.079903; F_list[7](0,2) = 0.0009054; F_list[7](0,3) = -0.0010765; 
  F_list[7](1,0) = 0.041588; F_list[7](1,1) = 0.99899; F_list[7](1,2) = -0.035855; F_list[7](1,3) = 0.012684;
  F_list[7](2,0) = -0.54662; F_list[7](2,1) = 0.044916; F_list[7](2,2) =  -0.32991; F_list[7](2,3) = 0.19318; 
  F_list[7](3,0) = 2.6624; F_list[7](3,1) = -0.10045; F_list[7](3,2) = -0.92455; F_list[7](3,3) = -0.26325; 
  
  G_list[7] = mat<double,mat_structure::rectangular>(4,2);
  G_list[7](0,0) = 0.0037112; G_list[7](0,1) = 0.0007361; 
  G_list[7](1,0) = -0.087051; G_list[7](1,1) = 9.3411e-6; 
  G_list[7](2,0) = -1.19844; G_list[7](2,1) = -4.1378e-4; 
  G_list[7](3,0) = -3.1927; G_list[7](3,1) = 9.2535e-4; 
  
  R_list[7] = mat<double,mat_structure::rectangular>(2,2);
  R_list[7](0,0) = 1.0; R_list[7](0,1) = 0.0; 
  R_list[7](1,0) = 0.0; R_list[7](1,1) = 1.0; 
  
  Q_list[7] = mat<double,mat_structure::rectangular>(4,4);
  Q_list[7](0,0) = 0.01; Q_list[7](0,1) = 0.0; Q_list[7](0,2) = 0.0; Q_list[7](0,3) = 0.0; 
  Q_list[7](1,0) = 0.0; Q_list[7](1,1) = 0.01; Q_list[7](1,2) = 0.0; Q_list[7](1,3) = 0.0; 
  Q_list[7](2,0) = 0.0; Q_list[7](2,1) = 0.0; Q_list[7](2,2) = 0.01; Q_list[7](2,3) = 0.0; 
  Q_list[7](3,0) = 0.0; Q_list[7](3,1) = 0.0; Q_list[7](3,2) = 0.0; Q_list[7](3,3) = 0.01; 
  
  
  F_list[8] = mat<double,mat_structure::rectangular>(5,5);
  F_list[8](0,0) =  95.407; F_list[8](0,1) = 1.9643; F_list[8](0,2) = 0.3597; F_list[8](0,3) = 0.0673; F_list[8](0,4) = 0.019;
  F_list[8](1,0) = 40.849; F_list[8](1,1) = 41.317; F_list[8](1,2) = 16.084; F_list[8](1,3) = 4.4679; F_list[8](1,4) = 1.1971;
  F_list[8](2,0) = 12.217; F_list[8](2,1) = 26.326; F_list[8](2,2) =  36.149; F_list[8](2,3) = 15.93; F_list[8](2,4) = 12.383;
  F_list[8](3,0) = 4.1118; F_list[8](3,1) = 12.858; F_list[8](3,2) = 27.209; F_list[8](3,3) = 21.442; F_list[8](3,4) = 40.976;
  F_list[8](4,0) = 0.1305; F_list[8](4,1) = 0.5808; F_list[8](4,2) = 1.875; F_list[8](4,3) = 3.6162; F_list[8](4,4) = 94.28;
  F_list[8] *= 0.01;
  
  G_list[8] = mat<double,mat_structure::rectangular>(5,2);
  G_list[8](0,0) = 0.0434; G_list[8](0,1) = -0.0122; 
  G_list[8](1,0) = 2.6606; G_list[8](1,1) = -1.0453; 
  G_list[8](2,0) = 3.753; G_list[8](2,1) = -5.51; 
  G_list[8](3,0) = 3.6076; G_list[8](3,1) = -6.6; 
  G_list[8](4,0) = 0.4617; G_list[8](4,1) = -0.9148; 
  G_list[8] *= 0.01;
  
  R_list[8] = mat<double,mat_structure::rectangular>(2,2);
  R_list[8](0,0) = 1.0; R_list[8](0,1) = 0.0; 
  R_list[8](1,0) = 0.0; R_list[8](1,1) = 1.0; 
  
  Q_list[8] = mat<double,mat_structure::rectangular>(5,5);
  Q_list[8](0,0) = 1.0; Q_list[8](0,1) = 0.0; Q_list[8](0,2) = 0.0; Q_list[8](0,3) = 0.0; Q_list[8](0,4) = 0.0;
  Q_list[8](1,0) = 0.0; Q_list[8](1,1) = 1.0; Q_list[8](1,2) = 0.0; Q_list[8](1,3) = 0.0; Q_list[8](1,4) = 0.0;
  Q_list[8](2,0) = 0.0; Q_list[8](2,1) = 0.0; Q_list[8](2,2) = 1.0; Q_list[8](2,3) = 0.0; Q_list[8](2,4) = 0.0;
  Q_list[8](3,0) = 0.0; Q_list[8](3,1) = 0.0; Q_list[8](3,2) = 0.0; Q_list[8](3,3) = 1.0; Q_list[8](3,4) = 0.0; 
  Q_list[8](4,0) = 0.0; Q_list[8](4,1) = 0.0; Q_list[8](4,2) = 0.0; Q_list[8](4,3) = 0.0; Q_list[8](4,4) = 1.0; 
  
  
  for(std::size_t i = 0; i < F_list.size(); ++i) {
    
    std::cout << "****** Problem " << i << " ****** " << std::endl;
    
    mat<double,mat_structure::rectangular> P(F_list[i].get_row_count(),F_list[i].get_col_count());
    solve_dare_problem(F_list[i], G_list[i], Q_list[i], R_list[i], P, 1e-6);
    std::cout << "P = " << P << std::endl;
    
    mat<double,mat_structure::rectangular> M_tmp = R_list[i];
    M_tmp += transpose_view(G_list[i]) * P * G_list[i];
    mat<double,mat_structure::rectangular> M2_tmp(G_list[i].get_col_count(),F_list[i].get_col_count());
    M2_tmp = transpose_view(G_list[i]) * P * F_list[i];
    mat<double,mat_structure::rectangular> Msol_tmp;
    linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
    mat<double,mat_structure::rectangular> X = 
      (transpose_view(F_list[i]) * P * F_list[i] 
     - transpose_view(F_list[i]) * P * G_list[i] * Msol_tmp + Q_list[i]);
    std::cout << "F' P F - F' P G ( R + G' P G )^-1 G' P F + Q = " 
              << X << std::endl;
    std::cout << "P_err = " 
              << (X - P) << std::endl;
    
    if(norm_1( X - P ) < 1e-5 * norm_1(P)) 
      ++passed;
    else
      RK_WARNING("DARE problem #" << i << " suffered a loss of precision!");
    if(norm_1( X - P ) < 1e-3 * norm_1(P)) 
      ++passed;
    else
      RK_ERROR("DARE problem #" << i << " did not pass!");
  };
  };
  
#if 0
  // These are the benchmark tests for the Continuous-time Algebraic Riccatic Equations (CARE)
  // NOTE So, far it seems that these are not working. Needs investigation.
  //      It is a bit weird that it wouldn't work since DARE is exactly the same algorithm
  {
  std::vector< mat<double, mat_structure::rectangular> > A_list(3);
  std::vector< mat<double, mat_structure::rectangular> > B_list(3);
  std::vector< mat<double, mat_structure::rectangular> > Q_list(3);
  std::vector< mat<double, mat_structure::rectangular> > R_list(3);
  
  A_list[0] = mat<double,mat_structure::rectangular>(2,2);
  A_list[0](0,0) = 0.0; A_list[0](0,1) = 1.0; 
  A_list[0](1,0) = 0.0; A_list[0](1,1) = 1.0; 
  
  B_list[0] = mat<double,mat_structure::rectangular>(2,1);
  B_list[0](0,0) = 0.0; 
  B_list[0](1,0) = 1.0; 
  
  R_list[0] = mat<double,mat_structure::rectangular>(1,1);
  R_list[0](0,0) = 1.0; 
  
  Q_list[0] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[0](0,0) = 1.0; Q_list[0](0,1) = 0.0; 
  Q_list[0](1,0) = 0.0; Q_list[0](1,1) = 2.0; 
  
  
  
  A_list[1] = mat<double,mat_structure::rectangular>(2,2);
  A_list[1](0,0) = 4.0; A_list[1](0,1) = 3.0; 
  A_list[1](1,0) = -4.5; A_list[1](1,1) = -3.5; 
  
  B_list[1] = mat<double,mat_structure::rectangular>(2,1);
  B_list[1](0,0) = 1.0;  
  B_list[1](1,0) = -1.0; 
  
  R_list[1] = mat<double,mat_structure::rectangular>(1,1);
  R_list[1](0,0) = 1.0;
  
  Q_list[1] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[1](0,0) = 9.0; Q_list[1](0,1) = 6.0; 
  Q_list[1](1,0) = 6.0; Q_list[1](1,1) = 4.0; 
  
  A = [4 3; -4.5 -3.5];
  G = [1 -1; -1 1];
  Q = [9 6; 6 4];
  X = (1+sqrt(2))*[9 6; 6 4];
  parout = [2, 1, 2];
  if nargout > 5
    B = [1; -1];  R = 1;  C = eye(2);  Q0 = Q;
  
  
  A_list[2] = mat<double,mat_structure::rectangular>(4,4);
  A_list[2](0,0) = 0.0; A_list[2](0,1) = 1.0; A_list[2](0,2) = 0.0; A_list[2](0,3) = 0.0; 
  A_list[2](1,0) = 0.0; A_list[2](1,1) = -1.89; A_list[2](1,2) = 0.39; A_list[2](1,3) = -5.53; 
  A_list[2](2,0) = 0.0; A_list[2](2,1) = -0.034; A_list[2](2,2) = -2.98; A_list[2](2,3) = 2.43; 
  A_list[2](3,0) = 0.034; A_list[2](3,1) = -0.0011; A_list[2](3,2) = -0.99; A_list[2](3,3) = -0.21; 
  
  B_list[2] = mat<double,mat_structure::rectangular>(4,2);
  B_list[2](0,0) = 0.0; B_list[2](0,1) = 0.0; 
  B_list[2](1,0) = 0.36; B_list[2](1,1) = -1.6; 
  B_list[2](2,0) = -0.95; B_list[2](2,1) = -0.032; 
  B_list[2](3,0) = 0.03; B_list[2](3,1) = 0.0; 
  
  R_list[2] = mat<double,mat_structure::rectangular>(2,2);
  R_list[2](0,0) = 1.0; R_list[2](0,1) = 0.0; 
  R_list[2](1,0) = 0.0; R_list[2](1,1) = 1.0; 
  
  Q_list[2] = mat<double,mat_structure::rectangular>(4,4);
  Q_list[2](0,0) = 2.313; Q_list[2](0,1) = 2.727; Q_list[2](0,2) = 0.688; Q_list[2](0,3) = 0.023; 
  Q_list[2](1,0) = 2.727; Q_list[2](1,1) = 4.271; Q_list[2](1,2) = 1.148; Q_list[2](1,3) = 0.323; 
  Q_list[2](2,0) = 0.688; Q_list[2](2,1) = 1.148; Q_list[2](2,2) = 0.313; Q_list[2](2,3) = 0.102; 
  Q_list[2](3,0) = 0.023; Q_list[2](3,1) = 0.323; Q_list[2](3,2) = 0.102; Q_list[2](3,3) = 0.083; 
  
  /*
  A_list[3] = mat<double,mat_structure::rectangular>(2,2);
  A_list[3](0,0) = 0.0; A_list[3](0,1) =  1.0; 
  A_list[3](1,0) = 0.0; A_list[3](1,1) = -1.0; 
  
  B_list[3] = mat<double,mat_structure::rectangular>(2,2);
  B_list[3](0,0) = 1.0; B_list[3](0,1) = 0.0; 
  B_list[3](1,0) = 2.0; B_list[3](1,1) = 1.0; 
  
  R_list[3] = mat<double,mat_structure::rectangular>(2,2);
  R_list[3](0,0) = 9.0; R_list[3](0,1) = 3.0; 
  R_list[3](1,0) = 3.0; R_list[3](1,1) = 1.0; 
  
  Q_list[3] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[3](0,0) = -4.0/11.0; Q_list[3](0,1) = -4.0/11.0; 
  Q_list[3](1,0) = -4.0/11.0; Q_list[3](1,1) =  7.0/11.0; 
  
  
  A_list[4] = mat<double,mat_structure::rectangular>(2,2);
  A_list[4](0,0) = 1.0; A_list[4](0,1) = 0.2; 
  A_list[4](1,0) = 0.0; A_list[4](1,1) = 1.0; 
  
  B_list[4] = mat<double,mat_structure::rectangular>(2,1);
  B_list[4](0,0) = 0.0; 
  B_list[4](1,0) = 1.0; 
  
  R_list[4] = mat<double,mat_structure::rectangular>(1,1);
  R_list[4](0,0) = 1.0; 
  
  Q_list[4] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[4](0,0) = 1.0; Q_list[4](0,1) = 0.0; 
  Q_list[4](1,0) = 0.0; Q_list[4](1,1) = 1.0; 
  
  
  A_list[5] = mat<double,mat_structure::rectangular>(2,2);
  A_list[5](0,0) = 0.0; A_list[5](0,1) = 1.0; 
  A_list[5](1,0) = 0.0; A_list[5](1,1) = 0.0; 
  
  B_list[5] = mat<double,mat_structure::rectangular>(2,1);
  B_list[5](0,0) = 0.0; 
  B_list[5](1,0) = 1.0; 
  
  R_list[5] = mat<double,mat_structure::rectangular>(1,1);
  R_list[5](0,0) = 1.0; 
  
  Q_list[5] = mat<double,mat_structure::rectangular>(2,2);
  Q_list[5](0,0) = 1.0; Q_list[5](0,1) = 2.0; 
  Q_list[5](1,0) = 2.0; Q_list[5](1,1) = 4.0; 
  */
  
  for(std::size_t i = 0; i < A_list.size(); ++i) {
    
    std::cout << "****** Problem " << i << " ****** " << std::endl;
    
    mat<double,mat_structure::rectangular> P(A_list[i].get_row_count(),A_list[i].get_col_count());
    solve_care_problem(A_list[i], B_list[i], Q_list[i], R_list[i], P, 1e-6);
    std::cout << "P = " << P << std::endl;
    
    mat<double,mat_structure::rectangular> M_tmp = R_list[i];
    mat<double,mat_structure::rectangular> M2_tmp(B_list[i].get_col_count(),A_list[i].get_col_count());
    M2_tmp = transpose_view(B_list[i]) * P;
    mat<double,mat_structure::rectangular> Msol_tmp;
    linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
    mat<double,mat_structure::rectangular> X = 
      (transpose_view(A_list[i]) * P + P * A_list[i] 
     - P * B_list[i] * Msol_tmp + Q_list[i]);
    std::cout << "Q + A^T P + P A - P B R^{-1} B^T P = " 
              << X << std::endl;
    
    if(norm_1( X ) < 1e-6 * norm_1(P)) 
      ++passed;
    else
      RK_WARNING("CARE problem #" << i << " suffered a loss of precision!");
    if(norm_1( X ) < 1e-3 * norm_1(P)) 
      ++passed;
    else
      RK_ERROR("CARE problem #" << i << " did not pass!");
  };
  };
#endif
  
  
  
  std::cout << "Algebraic Riccati Equation Tests results: " << passed << " out of " << possible_passes << std::endl;
  
  if(passed == possible_passes) {
    return 0;
  } else
    return 1;
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
#endif







