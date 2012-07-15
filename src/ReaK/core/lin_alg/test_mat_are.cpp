
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

#include "mat_alg.hpp"

#include "mat_are_solver.hpp"

#include "mat_norms.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

int main() {

  using namespace ReaK;
  
  
  
  std::size_t passed = 0;
  std::size_t expected_passes = 0;
  
  {
  
  std::vector< mat<double, mat_structure::rectangular> > F_list;
  std::vector< mat<double, mat_structure::rectangular> > G_list;
  std::vector< mat<double, mat_structure::rectangular> > R_list;
  std::vector< mat<double, mat_structure::rectangular> > Q_list;
  
  std::ifstream infile("darex_data.txt");
  while(infile) {
    std::string str_tmp;
    std::stringstream ss;
    std::getline(infile,str_tmp);
    if(!infile)
      break;
    std::size_t N, M;
    ss.str(str_tmp);
    ss >> N >> M;
    
    mat<double, mat_structure::rectangular> F_tmp(N,N);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < N; ++j)
        ss >> F_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> G_tmp(N,M);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < M; ++j)
        ss >> G_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> R_tmp(M,M);
    for(std::size_t i = 0; i < M; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < M; ++j)
        ss >> R_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> Q_tmp(N,N);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < N; ++j)
        ss >> Q_tmp(i,j);
    };
    
    F_list.push_back(F_tmp);
    G_list.push_back(G_tmp);
    R_list.push_back(R_tmp);
    Q_list.push_back(Q_tmp);
  };
  
  expected_passes += F_list.size();
  
  std::cout << "****** Discrete-time Algebraic Riccati Equations Tests with Tolerance of 1e-6 ****** " << std::endl;
  
  for(std::size_t i = 0; i < F_list.size(); ++i) {
    
    std::cout << "Problem " << i << ": " << std::endl;
    
    mat<double,mat_structure::rectangular> P(F_list[i].get_row_count(),F_list[i].get_col_count());
    solve_dare_problem(F_list[i], G_list[i], Q_list[i], R_list[i], P, 1e-6);
//     std::cout << "P = " << P << std::endl;
    
    mat<double,mat_structure::rectangular> M_tmp = R_list[i];
    M_tmp += transpose_view(G_list[i]) * P * G_list[i];
    mat<double,mat_structure::rectangular> M2_tmp(G_list[i].get_col_count(),F_list[i].get_col_count());
    M2_tmp = transpose_view(G_list[i]) * P * F_list[i];
    mat<double,mat_structure::rectangular> Msol_tmp;
    linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
    mat<double,mat_structure::rectangular> X = 
      (transpose_view(F_list[i]) * P * F_list[i] 
     - transpose_view(F_list[i]) * P * G_list[i] * Msol_tmp + Q_list[i]);
//     std::cout << "F' P F - F' P G ( R + G' P G )^-1 G' P F + Q = " 
//               << X << std::endl;
//     std::cout << "P_err = " 
//               << (X - P) << std::endl;
    double err_norm = norm_1( X - P );
    double P_norm = norm_1(P);
    std::cout << "\t|| Err || = " << std::setw(14) << err_norm 
              << " relative to || P || = " << std::setw(14) << P_norm 
              << " for a relative error of " << std::setw(14) << (err_norm / P_norm) << std::endl;
    
    if(err_norm > 2e-5 * P_norm) 
      std::cout << "\tWARNING: Significant loss of precision!" << std::endl;
    if(err_norm < 1e-3 * P_norm) 
      ++passed;
    else
      std::cout << "\tERROR: Such a loss of precision indicates a failure to solve the problem!" << std::endl;
  };
  };
  
#if 1
  // These are the benchmark tests for the Continuous-time Algebraic Riccatic Equations (CARE)
  // NOTE So, far it seems that these are not working. Needs investigation.
  //      It is a bit weird that it wouldn't work since DARE is exactly the same algorithm
  {
  
  std::vector< mat<double, mat_structure::rectangular> > A_list;
  std::vector< mat<double, mat_structure::rectangular> > B_list;
  std::vector< mat<double, mat_structure::rectangular> > R_list;
  std::vector< mat<double, mat_structure::rectangular> > Q_list;
  
  std::ifstream infile("carex_data.txt");
  while(infile) {
    std::string str_tmp;
    std::stringstream ss;
    std::getline(infile,str_tmp);
    if(!infile)
      break;
    std::size_t N, M;
    ss.str(str_tmp);
    ss >> N >> M;
    
    mat<double, mat_structure::rectangular> A_tmp(N,N);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < N; ++j)
        ss >> A_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> B_tmp(N,M);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < M; ++j)
        ss >> B_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> R_tmp(M,M);
    for(std::size_t i = 0; i < M; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < M; ++j)
        ss >> R_tmp(i,j);
    };
    mat<double, mat_structure::rectangular> Q_tmp(N,N);
    for(std::size_t i = 0; i < N; ++i) {
      std::getline(infile,str_tmp);
      ss.clear();
      ss.str(str_tmp);
      for(std::size_t j = 0; j < N; ++j)
        ss >> Q_tmp(i,j);
    };
    
    A_list.push_back(A_tmp);
    B_list.push_back(B_tmp);
    R_list.push_back(R_tmp);
    Q_list.push_back(Q_tmp);
  };
  
  expected_passes += A_list.size();
  
  std::cout << "****** Continuous-time Algebraic Riccati Equations Tests with Tolerance of 1e-6 ****** " << std::endl;
  
  for(std::size_t i = 0; i < A_list.size(); ++i) {
    
    std::cout << "Problem " << i << ": " << std::endl;
    
    mat<double,mat_structure::rectangular> P(A_list[i].get_row_count(),A_list[i].get_col_count());
    solve_care_problem(A_list[i], B_list[i], Q_list[i], R_list[i], P, 1e-6);
    //std::cout << "P = " << P << std::endl;
    
    mat<double,mat_structure::rectangular> M_tmp = R_list[i];
    mat<double,mat_structure::rectangular> M2_tmp(B_list[i].get_col_count(),A_list[i].get_col_count());
    M2_tmp = transpose_view(B_list[i]) * P;
    mat<double,mat_structure::rectangular> Msol_tmp;
    linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
    mat<double,mat_structure::rectangular> X = 
      (transpose_view(A_list[i]) * P + P * A_list[i] 
     - P * B_list[i] * Msol_tmp + Q_list[i]);
    //std::cout << "Q + A^T P + P A - P B R^{-1} B^T P = " 
    //          << X << std::endl;
    double err_norm = norm_1( X );
    double P_norm = norm_1(P);
    std::cout << "\t|| Err || = " << std::setw(14) << err_norm 
              << " relative to || P || = " << std::setw(14) << P_norm 
              << " for a relative error of " << std::setw(14) << (err_norm / P_norm) << std::endl;
              
    if(err_norm > 2e-5 * P_norm) 
      std::cout << "\tWARNING: Significant loss of precision!" << std::endl;
    if(err_norm < 1e-3 * P_norm) 
      ++passed;
    else
      std::cout << "\tERROR: Such a loss of precision indicates a failure to solve the problem!" << std::endl;
  };
  };
#endif
  
  
  std::cout << "*******************************************************************************************************" << std::endl;
  std::cout << "Algebraic Riccati Equation Tests results: " << passed << " out of " << expected_passes << " expected successful tests." << std::endl;
  
  
  if(passed == expected_passes) { 
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







