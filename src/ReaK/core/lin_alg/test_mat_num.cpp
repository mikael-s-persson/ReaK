
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

#include "mat_gaussian_elim.hpp"
#include "mat_cholesky.hpp"
#include "mat_jacobi_method.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_svd_method.hpp"

int main() {

  using namespace ReaK;

  unsigned int passed = 0;

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
  
  
  
  return 0;
};



