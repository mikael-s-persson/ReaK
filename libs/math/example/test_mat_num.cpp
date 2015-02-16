
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
#include <ReaK/math/lin_alg/mat_gaussian_elim.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>
#include <ReaK/math/lin_alg/mat_jacobi_method.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/math/lin_alg/mat_svd_method.hpp>
#include <ReaK/math/lin_alg/mat_schur_decomp.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>


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
  };
  
  {
  mat<double,mat_structure::rectangular> m_test(3,3);
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 0.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 0.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  
  std::cout << "Testing symmetric QR steps.." << std::endl;
  std::cout << m_test << std::endl;
  mat<double,mat_structure::square> m_test_Q(mat<double,mat_structure::identity>(3));
  detail::symmetric_QR_step(m_test, &m_test_Q, 1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  detail::symmetric_QR_step(m_test, &m_test_Q, 1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  detail::symmetric_QR_step(m_test, &m_test_Q, 1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 0.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 0.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0;
  m_test_Q = mat<double,mat_structure::identity>(3);
  std::cout << "Testing symmetric QR algorithm.." << std::endl;
  std::cout << m_test << std::endl;
  detail::symmetric_QRalg_impl(m_test, &m_test_Q, 1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 2.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 2.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  m_test_Q = mat<double,mat_structure::identity>(3);
  std::cout << "Testing tri-diagonal decomposition.." << std::endl;
  std::cout << m_test << std::endl;
  detail::decompose_TriDiag_impl(m_test,&m_test_Q,1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 2.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 2.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  m_test_Q = mat<double,mat_structure::identity>(3);
  std::cout << "Testing symmetric QR algorithm.." << std::endl;
  std::cout << m_test << std::endl;
  detail::symmetric_QRalg_impl(m_test,&m_test_Q,1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test << std::endl;
  std::cout << (m_test_Q * m_test * transpose_view(m_test_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  mat<double,mat_structure::rectangular> m_test4(4,4);
  mat<double,mat_structure::square> m_test4_Q(mat<double,mat_structure::identity>(4));
  m_test4(0,0) = 1.0; m_test4(0,1) = 3.0; m_test4(0,2) = 2.0; m_test4(0,3) =-1.0; 
  m_test4(1,0) = 3.0; m_test4(1,1) = 5.0; m_test4(1,2) = 2.0; m_test4(1,3) = 5.0; 
  m_test4(2,0) = 2.0; m_test4(2,1) = 2.0; m_test4(2,2) = 4.0; m_test4(2,3) = 3.0; 
  m_test4(3,0) =-1.0; m_test4(3,1) = 5.0; m_test4(3,2) = 3.0; m_test4(3,3) = -2.0;
  m_test4_Q = mat<double,mat_structure::identity>(4);
  std::cout << "Testing symmetric QR algorithm.." << std::endl;
  std::cout << m_test4 << std::endl;
  detail::symmetric_QRalg_impl(m_test4,&m_test4_Q,1e-6);
  std::cout << m_test4_Q << std::endl;
  std::cout << m_test4 << std::endl;
  std::cout << (m_test4_Q * m_test4 * transpose_view(m_test4_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  mat<double,mat_structure::square> m_test_L(3,0.0);
  mat<double,mat_structure::square> m_test_D(3,0.0);
  mat<double,mat_structure::square> m_test_sqr(3,0.0);
  mat<double,mat_structure::square> m_test_inv(3,0.0);
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing LDL decomposition.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  decompose_LDL(m_test_sqr,m_test_L,m_test_D,1e-6);
  std::cout << m_test_L << std::endl;
  std::cout << m_test_D << std::endl;
  std::cout << (m_test_L * m_test_D * transpose_view(m_test_L)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing LDL inversion.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  invert_LDL(m_test_sqr,m_test_inv,1e-6);
  std::cout << m_test_inv << std::endl;
  std::cout << (m_test_inv * m_test_sqr) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 0.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 0.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing Band-Cholesky decomposition.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  decompose_BandCholesky(m_test_sqr,m_test_L,1,1e-6);
  std::cout << m_test_L << std::endl;
  std::cout << (m_test_L * transpose_view(m_test_L)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 0.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 0.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing Band-Cholesky inversion.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  invert_BandCholesky(m_test_sqr,m_test_inv,1,1e-6);
  std::cout << m_test_inv << std::endl;
  std::cout << (m_test_inv * m_test_sqr) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_D = mat<double,mat_structure::identity>(3);
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 0.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 0.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing TriDiagLDL decomposition.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  decompose_TriDiagLDL(m_test_sqr,m_test_L,m_test_D,1e-6);
  std::cout << m_test_L << std::endl;
  std::cout << m_test_D << std::endl;
  std::cout << (m_test_L * m_test_D * transpose_view(m_test_L)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 1.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  m_test_Q = mat<double,mat_structure::identity>(3);
  m_test_D = mat<double,mat_structure::identity>(3);
  std::cout << "Testing SymQR eigensolve.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  eigensolve_SymQR(m_test_sqr,m_test_Q,m_test_D,1e-6);
  std::cout << m_test_Q << std::endl;
  std::cout << m_test_D << std::endl;
  std::cout << (m_test_Q * m_test_D * transpose_view(m_test_Q)) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 1.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing SymQR inversion.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  pseudoinvert_SymQR(m_test_sqr,m_test_inv,1e-6);
  std::cout << m_test_inv << std::endl;
  std::cout << (m_test_inv * m_test_sqr) << std::endl;
  std::cout << "Done!" << std::endl;
  
  m_test_sqr(0,0) = 1.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  std::cout << "Testing QR inversion.." << std::endl;
  std::cout << m_test_sqr << std::endl;
  pseudoinvert_QR(m_test_sqr,m_test_inv,1e-6);
  std::cout << m_test_inv << std::endl;
  std::cout << (m_test_inv * m_test_sqr) << std::endl;
  std::cout << "Done!" << std::endl;
  
  };
  
  
  return 0;
};



