
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

int main() {

  using namespace ReaK;

  unsigned int passed = 0;

#if 0
  try {
 
    if(true){
      RK_NOTICE(2,"/*********************************************/");
      RK_NOTICE(2,"/********* MATRIX NUM-METHODS TESTS **********/");
      RK_NOTICE(2,"/*********************************************/");

      mat<double,mat_structure::symmetric> m_gauss(2.0,-1.0,0.0,2.0,-1.0,2.0);
      RK_NOTICE(2,"Testing the following matrix: " << m_gauss);
    
      mat<double,mat_structure::square> m_gauss_inv(3);
      invert_gaussian(m_gauss,m_gauss_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Gaussian: " << m_gauss_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_gauss_inv));

      mat<double,mat_structure::square> m_plu_inv(mat<double,mat_structure::identity>(3));
      invert_PLU(m_gauss,m_plu_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with PLU: " << m_plu_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_plu_inv));
    
      mat<double,mat_structure::symmetric> m_gauss2_inv(3);
      invert_gaussian(m_gauss,m_gauss2_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Gaussian Symmetric: " << m_gauss2_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_gauss2_inv));

      mat<double,mat_structure::symmetric> m_plu2_inv(mat<double,mat_structure::identity>(3));
      invert_PLU(m_gauss,m_plu2_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with PLU Symmetric: " << m_plu2_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_plu2_inv));

      mat<double,mat_structure::symmetric> m_cholesky_inv(mat<double,mat_structure::identity>(3));
      invert_Cholesky(m_gauss,m_cholesky_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Cholesky: " << m_cholesky_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_cholesky_inv));
    
      mat<double,mat_structure::square> m_cholesky2_inv(mat<double,mat_structure::identity>(3));
      invert_Cholesky(m_gauss,m_cholesky2_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Cholesky Square: " << m_cholesky2_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_cholesky2_inv));
     
      mat<double,mat_structure::diagonal> m_gauss_diag(vect<double,3>(2,1,0.5));
      RK_NOTICE(2,"Testing the following matrix: " << m_gauss_diag);
    
      mat<double,mat_structure::square> m_cholesky3_inv(mat<double,mat_structure::identity>(3));
      invert_Cholesky(m_gauss_diag,m_cholesky3_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Diag Cholesky Square: " << m_cholesky3_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss_diag * m_cholesky3_inv));
    
      mat<double,mat_structure::diagonal> m_cholesky4_inv(mat<double,mat_structure::identity>(3));
      invert_Cholesky(m_gauss_diag,m_cholesky4_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Diag Cholesky Square: " << m_cholesky4_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss_diag * m_cholesky4_inv));
     
      mat<double,mat_structure::symmetric> m_jacobi_pinv(3);
      pseudoinvert_Jacobi(m_gauss,m_jacobi_pinv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with Jacobi: " << m_jacobi_pinv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_jacobi_pinv));
   
      mat<double,mat_structure::square> m_qr_inv(3);
      pseudoinvert_QR(m_gauss,m_qr_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with QR: " << m_qr_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_qr_inv));

      mat<double,mat_structure::square> m_svd_inv(3);
      pseudoinvert_SVD(m_gauss,m_svd_inv,double(1E-15));
      RK_NOTICE(2,"The inverse was found with SVD: " << m_svd_inv);
      RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_svd_inv));

      mat<double,mat_structure::diagonal> m_jacobi_E(3);
      mat<double,mat_structure::square> m_jacobi_Q(3);
      eigensolve_Jacobi(m_gauss,m_jacobi_E,m_jacobi_Q,double(1E-15));
      RK_NOTICE(2,"The eigenvalues were found with Jacobi.");
      RK_NOTICE(2,"The eigenvalues are = " << m_jacobi_E);
      RK_NOTICE(2,"The eigenvectors are = " << m_jacobi_Q);
      RK_NOTICE(2,"Q * E * Qt = " << m_jacobi_Q * m_jacobi_E * transpose(m_jacobi_Q));

      mat<double,mat_structure::diagonal> m_qr_E(3,true);
      mat<double,mat_structure::square> m_qr_Q(3);
      //eigensolve_QR(m_gauss,m_qr_E,m_qr_Q,50,double(1E-6));
      RK_NOTICE(2,"The eigenvalues were found with QR.");
      RK_NOTICE(2,"The eigenvalues are = " << m_qr_E);
      RK_NOTICE(2,"The eigenvectors are = " << m_qr_Q);
      RK_NOTICE(2,"Q * E * Qt = " << m_qr_Q * m_qr_E * transpose(m_qr_Q));

      mat<double,mat_structure::diagonal> m_svd_E(3);
      mat<double,mat_structure::square> m_svd_U(3);
      mat<double,mat_structure::square> m_svd_V(3);
      decompose_SVD(m_gauss,m_svd_U,m_svd_E,m_svd_V,double(1E-15));
      RK_NOTICE(2,"The eigenvalues were found with SVD.");
      RK_NOTICE(2,"The eigenvalues are = " << m_svd_E);
      RK_NOTICE(2,"V = " << m_svd_V);
      RK_NOTICE(2,"U = " << m_svd_U);
      RK_NOTICE(2,"Q * E * Qt = " << m_svd_U * m_svd_E * transpose(m_svd_V));
    
    };
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
  };
  
  RK_NOTICE(2,"There were " << passed << " successful tests passed on the math_gen library, out of 45 possible successes.");
  
#endif
  
  
  mat<double,mat_structure::rectangular> m_test(2,3);
  m_test(0,0) = 1.0; m_test(0,1) = 2.0; m_test(0,2) = 3.0; 
  m_test(1,0) = 4.0; m_test(1,1) = 5.0; m_test(1,2) = 6.0; 
  
  std::cout << m_test << std::endl;
  mat<double,mat_structure::rectangular> m_test_R(2,3);
  mat<double,mat_structure::square> m_test_Q(mat<double,mat_structure::identity>(2));
  detail::decompose_QR_impl< mat<double,mat_structure::rectangular>, mat<double,mat_structure::square> >(m_test,&m_test_Q,1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test) << std::endl;
  
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
  
  return 0;
};



