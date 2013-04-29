
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

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE mat_alg_riccati_eqs
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( mat_darex_tests )
{

  using namespace ReaK;
  
  std::vector< mat<double, mat_structure::rectangular> > F_list;
  std::vector< mat<double, mat_structure::rectangular> > G_list;
  std::vector< mat<double, mat_structure::rectangular> > R_list;
  std::vector< mat<double, mat_structure::rectangular> > Q_list;
  
  std::ifstream infile("are_data/darex_data.txt");
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
  
  for(std::size_t i = 0; i < F_list.size(); ++i) {
    
    try {
      mat<double,mat_structure::rectangular> P(F_list[i].get_row_count(),F_list[i].get_col_count());
      solve_dare_problem(F_list[i], G_list[i], Q_list[i], R_list[i], P, 1e-6);
      
      mat<double,mat_structure::rectangular> M_tmp = R_list[i];
      M_tmp += transpose_view(G_list[i]) * P * G_list[i];
      mat<double,mat_structure::rectangular> M2_tmp(G_list[i].get_col_count(),F_list[i].get_col_count());
      M2_tmp = transpose_view(G_list[i]) * P * F_list[i];
      mat<double,mat_structure::rectangular> Msol_tmp;
      linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
      mat<double,mat_structure::rectangular> X = 
        (transpose_view(F_list[i]) * P * F_list[i] 
       - transpose_view(F_list[i]) * P * G_list[i] * Msol_tmp + Q_list[i]);
      double err_norm = norm_1( X - P );
      double P_norm = norm_1(P);
      BOOST_WARN_MESSAGE(  (err_norm < 2e-5 * P_norm), "Significant loss of precision on DAREX problem " << i << "!" );
      BOOST_CHECK_MESSAGE( (err_norm < 1e-3 * P_norm), "Failed to solve DAREX problem " << i << "! Loss of precision unacceptable!" );
    } catch(...) {
      BOOST_ERROR( "DAREX problem " << i << " caused an exception to be thrown!" );
    };
    
  };
  
};

BOOST_AUTO_TEST_CASE( mat_carex_tests )
{

  using namespace ReaK;
  
  // These are the benchmark tests for the Continuous-time Algebraic Riccatic Equations (CARE)
  // NOTE So, far it seems that these are not working. Needs investigation.
  //      It is a bit weird that it wouldn't work since DARE is exactly the same algorithm
  
  std::vector< mat<double, mat_structure::rectangular> > A_list;
  std::vector< mat<double, mat_structure::rectangular> > B_list;
  std::vector< mat<double, mat_structure::rectangular> > R_list;
  std::vector< mat<double, mat_structure::rectangular> > Q_list;
  
  std::ifstream infile("are_data/carex_data.txt");
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
  
  for(std::size_t i = 0; i < A_list.size(); ++i) {
    
    try {
      mat<double,mat_structure::rectangular> P(A_list[i].get_row_count(),A_list[i].get_col_count());
      solve_care_problem(A_list[i], B_list[i], Q_list[i], R_list[i], P, 1e-6);
      
      mat<double,mat_structure::rectangular> M_tmp = R_list[i];
      mat<double,mat_structure::rectangular> M2_tmp(B_list[i].get_col_count(),A_list[i].get_col_count());
      M2_tmp = transpose_view(B_list[i]) * P;
      mat<double,mat_structure::rectangular> Msol_tmp;
      linlsq_QR(M_tmp,Msol_tmp,M2_tmp);
      mat<double,mat_structure::rectangular> X = 
        (transpose_view(A_list[i]) * P + P * A_list[i] 
      - P * B_list[i] * Msol_tmp + Q_list[i]);
      double err_norm = norm_1( X );
      double P_norm = norm_1(P);
      BOOST_WARN_MESSAGE(  (err_norm < 2e-5 * P_norm), "Significant loss of precision on CAREX problem " << i << "!" );
      BOOST_CHECK_MESSAGE( (err_norm < 1e-3 * P_norm), "Failed to solve CAREX problem " << i << "! Loss of precision unacceptable!" );
    } catch(...) {
      BOOST_ERROR( "CAREX problem " << i << " caused an exception to be thrown!" );
    };
  };
  
};







