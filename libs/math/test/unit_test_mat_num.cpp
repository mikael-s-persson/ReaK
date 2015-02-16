
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
#include <ReaK/math/lin_alg/mat_ctrl_decomp.hpp>
#include <ReaK/math/lin_alg/mat_balance.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE mat_num
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( mat_fullrank_inversion_tests )
{
  
  using namespace ReaK;
  
  
  mat<double,mat_structure::symmetric> m_gauss(2.0, -1.0,  0.0,  
                                                     2.0, -1.0, 
                                                           2.0);
  mat<double,mat_structure::symmetric> m_gauss_trueinv(0.75, 0.5,  0.25,  
                                                             1.0,  0.5, 
                                                                   0.75);
  
  mat<double,mat_structure::square> m_gauss_inv(3);
  BOOST_CHECK_NO_THROW( invert_gaussian(m_gauss,m_gauss_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_gauss_inv - m_gauss_trueinv,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_gauss_inv),std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_plu_inv(mat<double,mat_structure::identity>(3));
  BOOST_CHECK_NO_THROW( invert_PLU(m_gauss,m_plu_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_plu_inv - m_gauss_trueinv,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_plu_inv),std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::symmetric> m_gauss2_inv(3);
  BOOST_CHECK_NO_THROW( invert_gaussian(m_gauss,m_gauss2_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_gauss2_inv - m_gauss_trueinv,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_gauss2_inv),std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::symmetric> m_plu2_inv(mat<double,mat_structure::identity>(3));
  BOOST_CHECK_NO_THROW( invert_PLU(m_gauss,m_plu2_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_plu2_inv - m_gauss_trueinv,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_plu2_inv),std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::symmetric> m_cholesky_inv(mat<double,mat_structure::identity>(3));
  BOOST_CHECK_NO_THROW( invert_Cholesky(m_gauss,m_cholesky_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_cholesky_inv - m_gauss_trueinv, 2.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_cholesky_inv), 4.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_cholesky2_inv(mat<double,mat_structure::identity>(3));
  invert_Cholesky(m_gauss,m_cholesky2_inv,double(1E-15));
  BOOST_CHECK_NO_THROW( invert_Cholesky(m_gauss,m_cholesky2_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_cholesky2_inv - m_gauss_trueinv, 2.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_cholesky2_inv), 4.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::diagonal> m_gauss_diag(vect<double,3>(2,1,0.5));
  mat<double,mat_structure::diagonal> m_gauss_diag_trueinv(vect<double,3>(0.5,1,2.0));
  
  mat<double,mat_structure::square> m_cholesky3_inv(mat<double,mat_structure::identity>(3));
  BOOST_CHECK_NO_THROW( invert_Cholesky(m_gauss_diag,m_cholesky3_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_cholesky3_inv - m_gauss_diag_trueinv, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss_diag * m_cholesky3_inv), std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::diagonal> m_cholesky4_inv(mat<double,mat_structure::identity>(3));
  BOOST_CHECK_NO_THROW( invert_Cholesky(m_gauss_diag,m_cholesky4_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_cholesky4_inv - m_gauss_diag_trueinv, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss_diag * m_cholesky4_inv), std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::symmetric> m_jacobi_pinv(3);
  BOOST_CHECK_NO_THROW( pseudoinvert_Jacobi(m_gauss,m_jacobi_pinv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_jacobi_pinv - m_gauss_trueinv, 32.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_jacobi_pinv), 32.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_qr_inv(3);
  BOOST_CHECK_NO_THROW( pseudoinvert_QR(m_gauss,m_qr_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_qr_inv - m_gauss_trueinv, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_qr_inv), std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_svd_inv(3);
  BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_gauss,m_svd_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_svd_inv - m_gauss_trueinv, 4.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_svd_inv), 4.0 * std::numeric_limits<double>::epsilon()) ) );
  
  
  mat<double,mat_structure::square> m_ldl_inv(3);
  BOOST_CHECK_NO_THROW( invert_LDL(m_gauss,m_ldl_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_ldl_inv - m_gauss_trueinv, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_ldl_inv), std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_bandchol_inv(3);
  BOOST_CHECK_NO_THROW( invert_BandCholesky(m_gauss, m_bandchol_inv, 1, double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_bandchol_inv - m_gauss_trueinv, 4.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_bandchol_inv), 4.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_symqr_pinv(3);
  BOOST_CHECK_NO_THROW( pseudoinvert_SymQR(m_gauss,m_symqr_pinv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_symqr_pinv - m_gauss_trueinv, 4.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_symqr_pinv), 4.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::square> m_tridiagldl_inv(3);
  BOOST_CHECK_NO_THROW( invert_TriDiagLDL(m_gauss,m_tridiagldl_inv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_tridiagldl_inv - m_gauss_trueinv, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_tridiagldl_inv), std::numeric_limits<double>::epsilon()) ) );
  
};


BOOST_AUTO_TEST_CASE( mat_left_pseudoinversion_tests )
{
  using namespace ReaK;
  
  mat<double,mat_structure::rectangular> m_test(3,2);
  m_test(0,0) = 1.0; m_test(0,1) = 3.0;
  m_test(1,0) = 4.0; m_test(1,1) = 5.0;
  m_test(2,0) = 7.0; m_test(2,1) = 8.0;
  
  mat<double,mat_structure::rectangular> m_qr_inv(2,3);
  BOOST_CHECK_NO_THROW( pseudoinvert_QR(m_test,m_qr_inv, 1e-15) );
  BOOST_CHECK( ( is_identity_mat((m_qr_inv * m_test), 1e-14) ) );
  
  mat<double,mat_structure::rectangular> m_svd_inv(2,3);
  BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_test,m_svd_inv, 1e-15) );
  BOOST_CHECK( ( is_identity_mat((m_svd_inv * m_test), 1e-14) ) );
  
  mat<double,mat_structure::rectangular> m_rrqr_inv(2,3);
  try {
    linlsq_RRQR(m_test, m_rrqr_inv, mat<double,mat_structure::identity>(3), 1e-15);
  } catch(std::exception& e) {
    std::cout << "Exception thrown: " << e.what() << std::endl;
  };
//   BOOST_CHECK_NO_THROW( linlsq_RRQR(m_test, m_rrqr_inv, mat<double,mat_structure::identity>(3), 1e-15) );
  BOOST_CHECK( ( is_identity_mat((m_rrqr_inv * m_test), 1e-14) ) );
  
  // last column is equal to first column + 0.5 * second column (i.e., it is not full column-rank)
  mat<double,mat_structure::rectangular> m_rankdef_test(4,3);
  m_rankdef_test(0,0) = 1.0; m_rankdef_test(0,1) = 3.0; m_rankdef_test(0,2) = 2.5;
  m_rankdef_test(1,0) = 4.0; m_rankdef_test(1,1) = 5.0; m_rankdef_test(1,2) = 6.5;
  m_rankdef_test(2,0) = 7.0; m_rankdef_test(2,1) = 8.0; m_rankdef_test(2,2) = 11.0;
  m_rankdef_test(3,0) = 2.0; m_rankdef_test(3,1) = 4.0; m_rankdef_test(3,2) = 4.0;
  
  mat<double,mat_structure::rectangular> m_svd_rankdef_inv(3,4);
  BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_rankdef_test,m_svd_rankdef_inv, 1e-10) );
//   std::cout << "SVD rankdef pinv = " << m_svd_rankdef_inv << std::endl;
//   std::cout << "SVD rankdef pinv * A = " << (m_svd_rankdef_inv * m_rankdef_test) << std::endl;
//   BOOST_CHECK( ( is_identity_mat((m_svd_rankdef_inv * m_rankdef_test), 1e-14) ) );
  
  // NOTE: I don't know what else to do to check that the solution is correct, besides comparing to SVD solution.
  
  mat<double,mat_structure::rectangular> m_rrqr_rankdef_inv(3,4);
  BOOST_CHECK_NO_THROW( linlsq_RRQR(m_rankdef_test, m_rrqr_rankdef_inv, mat<double,mat_structure::identity>(4), 1e-10) );
//   std::cout << "RRQR rankdef pinv = " << m_rrqr_rankdef_inv << std::endl;
//   std::cout << "RRQR rankdef pinv * A = " << (m_rrqr_rankdef_inv * m_rankdef_test) << std::endl;
//   BOOST_CHECK( ( is_identity_mat((m_rrqr_rankdef_inv * m_rankdef_test), 1e-14) ) );
  
  BOOST_CHECK( ( is_null_mat((m_svd_rankdef_inv - m_rrqr_rankdef_inv), 1e-14) ) );
  
  
};


BOOST_AUTO_TEST_CASE( mat_right_pseudoinversion_tests )
{
  using namespace ReaK;
  
  mat<double,mat_structure::rectangular> m_test(2,3);
  m_test(0,0) = 1.0; m_test(0,1) = 2.0; m_test(0,2) = 3.0; 
  m_test(1,0) = 4.0; m_test(1,1) = 5.0; m_test(1,2) = 6.0; 
  
  mat<double,mat_structure::rectangular> m_test_R = m_test;
  mat<double,mat_structure::square> m_test_Q(mat<double,mat_structure::identity>(2));
  BOOST_CHECK_NO_THROW( (detail::decompose_QR_impl< mat<double,mat_structure::rectangular>, mat<double,mat_structure::square> >(m_test_R,&m_test_Q,1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_upper_triangular(m_test_R, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R - m_test), 4.0 * std::numeric_limits<double>::epsilon()) );
    
  mat<double,mat_structure::rectangular> m_test_minnorm_x(3,2);
  mat<double,mat_structure::rectangular> m_test_minnorm_b(2,2);
  m_test_minnorm_b(0,0) = 4.0; m_test_minnorm_b(0,1) = 5.0; 
  m_test_minnorm_b(1,0) = 7.0; m_test_minnorm_b(1,1) = 8.0; 
  BOOST_CHECK_NO_THROW( (minnorm_QR(m_test, m_test_minnorm_x, m_test_minnorm_b, 1e-6)) );
  BOOST_CHECK( is_null_mat( m_test * m_test_minnorm_x - m_test_minnorm_b, 1e-6) );
  
  mat<double,mat_structure::rectangular> m_test_minnormRRQR_x(3,2);
  BOOST_CHECK_NO_THROW( (minnorm_RRQR(m_test, m_test_minnormRRQR_x, m_test_minnorm_b, 1e-6)) );
  BOOST_CHECK( is_null_mat( m_test * m_test_minnormRRQR_x - m_test_minnorm_b, 1e-6) );
  
  
  /* SVD pseudo-invert for a min-norm problem. */
  mat<double,mat_structure::rectangular> m_test_svd_pinv(3,2);
  BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_test, m_test_svd_pinv, 1e-6) );
  mat<double,mat_structure::rectangular> m_test_svd_x = m_test_svd_pinv * m_test_minnorm_b;
  BOOST_CHECK( is_null_mat( m_test * m_test_svd_x - m_test_minnorm_b, 1e-6) );
  
  
  // last row is equal to first column + 0.5 * second column (i.e., it is not full column-rank)
  mat<double,mat_structure::rectangular> m_rankdef_test(3,4);
  m_rankdef_test(0,0) = 1.0; m_rankdef_test(0,1) = 4.0; m_rankdef_test(0,2) = 7.0;  m_rankdef_test(0,3) = 2.0; 
  m_rankdef_test(1,0) = 3.0; m_rankdef_test(1,1) = 5.0; m_rankdef_test(1,2) = 8.0;  m_rankdef_test(1,3) = 4.0; 
  m_rankdef_test(2,0) = 2.5; m_rankdef_test(2,1) = 6.5; m_rankdef_test(2,2) = 11.0; m_rankdef_test(2,3) = 4.0;
  
  mat<double,mat_structure::rectangular> m_rankdef_minnorm_b(3,2);
  m_rankdef_minnorm_b(0,0) = 4.0; m_rankdef_minnorm_b(0,1) = 5.0; 
  m_rankdef_minnorm_b(1,0) = 7.0; m_rankdef_minnorm_b(1,1) = 8.0; 
  m_rankdef_minnorm_b(2,0) = 4.0; m_rankdef_minnorm_b(2,1) = 2.0; 
  
  mat<double,mat_structure::rectangular> m_rankdef_minnormRRQR_x(4,2);
  BOOST_CHECK_NO_THROW( (minnorm_RRQR(m_rankdef_test, m_rankdef_minnormRRQR_x, m_rankdef_minnorm_b, 1e-6)) );
//   std::cout << "RRQR rank-def x = " << m_rankdef_minnormRRQR_x << std::endl;
//   std::cout << "RRQR rank-def M * x = " << (m_rankdef_test * m_rankdef_minnormRRQR_x) << std::endl;
//   BOOST_CHECK( is_null_mat( m_rankdef_test * m_rankdef_minnormRRQR_x - m_rankdef_minnorm_b, 1e-6) );
  
  // NOTE: I don't know what else to do to check that the solution is correct, besides comparing to SVD solution.
  
  mat<double,mat_structure::rectangular> m_rankdef_svd_pinv(4,3);
  BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_rankdef_test, m_rankdef_svd_pinv, 1e-6) );
  mat<double,mat_structure::rectangular> m_rankdef_svd_x = m_rankdef_svd_pinv * m_rankdef_minnorm_b;
//   std::cout << "SVD rank-def x = " << m_rankdef_svd_x << std::endl;
//   std::cout << "SVD rank-def M * x = " << (m_rankdef_test * m_rankdef_svd_x) << std::endl;
//   BOOST_CHECK( is_null_mat( m_rankdef_test * m_rankdef_svd_x - m_rankdef_minnorm_b, 1e-6) );
  
  BOOST_CHECK( ( is_null_mat((m_rankdef_svd_x - m_rankdef_minnormRRQR_x), 1e-14) ) );
  
};


BOOST_AUTO_TEST_CASE( mat_eigenvalues_tests )
{
  
  using namespace ReaK;
  
  
  mat<double,mat_structure::symmetric> m_gauss(2.0, -1.0,  0.0,  
                                                     2.0, -1.0, 
                                                           2.0);
  mat<double,mat_structure::symmetric> m_gauss_trueinv(0.75, 0.5,  0.25,  
                                                             1.0,  0.5, 
                                                                   0.75);
  
  vect_n<double> m_gauss_trueeig( 1.0 / (1.0 + std::sqrt(0.5)), 2.0, 2.0 + std::sqrt(2.0));
  
  mat<double,mat_structure::diagonal> m_jacobi_E(3);
  mat<double,mat_structure::square> m_jacobi_Q(3);
  BOOST_CHECK_NO_THROW( eigensolve_Jacobi(m_gauss,m_jacobi_E,m_jacobi_Q,double(1E-15)) );
  vect_n<double> m_jacobi_Es(m_jacobi_E(0,0), m_jacobi_E(1,1), m_jacobi_E(2,2));
  std::sort(m_jacobi_Es.begin(), m_jacobi_Es.end());
  BOOST_CHECK( ( norm_inf(m_jacobi_Es - m_gauss_trueeig) < 1e-14 ) );
  BOOST_CHECK( ( is_null_mat((m_jacobi_Q * m_jacobi_E * transpose(m_jacobi_Q) - m_gauss), 1e-14) ) );
  
  mat<double,mat_structure::diagonal> m_qr_E(3,true);
  mat<double,mat_structure::square> m_qr_Q(3);
  BOOST_CHECK_NO_THROW( eigensolve_SymQR(m_gauss,m_qr_Q,m_qr_E,double(1E-15)) );
  vect_n<double> m_qr_Es(m_qr_E(0,0), m_qr_E(1,1), m_qr_E(2,2));
  std::sort(m_qr_Es.begin(), m_qr_Es.end());
  BOOST_CHECK( ( norm_inf(m_qr_Es - m_gauss_trueeig) < 1e-14 ) );
  BOOST_CHECK( ( is_null_mat((m_qr_Q * m_qr_E * transpose(m_qr_Q) - m_gauss), 1e-14) ) );
  
  mat<double,mat_structure::diagonal> m_svd_E(3);
  mat<double,mat_structure::square> m_svd_U(3);
  mat<double,mat_structure::square> m_svd_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(m_gauss,m_svd_U,m_svd_E,m_svd_V,double(1E-15)) );
  vect_n<double> m_svd_Es(m_svd_E(0,0), m_svd_E(1,1), m_svd_E(2,2));
  std::sort(m_svd_Es.begin(), m_svd_Es.end());
  BOOST_CHECK( ( norm_inf(m_svd_Es - m_gauss_trueeig) < 10.0 * std::numeric_limits<double>::epsilon() ) );
  BOOST_CHECK( ( is_null_mat((m_svd_U * m_svd_E * transpose(m_svd_V) - m_gauss), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  
};



BOOST_AUTO_TEST_CASE( mat_decompositions_tests )
{
  
  using namespace ReaK;
  
  mat<double,mat_structure::rectangular> m_test(3,3);
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 0.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 0.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  
  mat<double,mat_structure::rectangular> m_test_Gram = m_test;
  BOOST_CHECK_NO_THROW( (orthogonalize_StableGramSchmidt(m_test_Gram, true, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Gram * transpose(m_test_Gram), 1e-6) );
  
  mat<double,mat_structure::rectangular> m_test_R = m_test;
  mat<double,mat_structure::square> m_test_Q(mat<double,mat_structure::identity>(3));
  
  BOOST_CHECK_NO_THROW( (detail::symmetric_QR_step(m_test_R, &m_test_Q, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 10.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R * transpose_view(m_test_Q) - m_test), 50.0 * std::numeric_limits<double>::epsilon()) );
  
  BOOST_CHECK_NO_THROW( (detail::symmetric_QR_step(m_test_R, &m_test_Q, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 10.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R * transpose_view(m_test_Q) - m_test), 50.0 * std::numeric_limits<double>::epsilon()) );
  
  BOOST_CHECK_NO_THROW( (detail::symmetric_QR_step(m_test_R, &m_test_Q, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 10.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R * transpose_view(m_test_Q) - m_test), 50.0 * std::numeric_limits<double>::epsilon()) );
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 0.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 0.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0;
  m_test_R = m_test;
  m_test_Q = mat<double,mat_structure::identity>(3);
  BOOST_CHECK_NO_THROW( (detail::symmetric_QRalg_impl(m_test_R, &m_test_Q, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_upper_triangular(m_test_R, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R * transpose_view(m_test_Q) - m_test), 1e-6) );
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 2.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 2.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  m_test_R = m_test;
  m_test_Q = mat<double,mat_structure::identity>(3);
  BOOST_CHECK_NO_THROW( (detail::decompose_TriDiag_impl(m_test_R,&m_test_Q,1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_tri_diagonal(m_test_R, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat(((m_test_Q * m_test_R * transpose_view(m_test_Q)) - m_test), 1e-6) );
  
  m_test(0,0) = 1.0; m_test(0,1) = 3.0; m_test(0,2) = 2.0; 
  m_test(1,0) = 3.0; m_test(1,1) = 5.0; m_test(1,2) = 2.0; 
  m_test(2,0) = 2.0; m_test(2,1) = 2.0; m_test(2,2) = 4.0; 
  m_test_R = m_test;
  m_test_Q = mat<double,mat_structure::identity>(3);
  BOOST_CHECK_NO_THROW( (detail::symmetric_QRalg_impl(m_test_R, &m_test_Q, 1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test_Q * transpose(m_test_Q), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_upper_triangular(m_test_R, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test_Q * m_test_R * transpose_view(m_test_Q) - m_test), 1e-6) );
  
  mat<double,mat_structure::rectangular> m_test4(4,4);
  mat<double,mat_structure::square> m_test4_Q(mat<double,mat_structure::identity>(4));
  m_test4(0,0) = 1.0; m_test4(0,1) = 3.0; m_test4(0,2) = 2.0; m_test4(0,3) =-1.0; 
  m_test4(1,0) = 3.0; m_test4(1,1) = 5.0; m_test4(1,2) = 2.0; m_test4(1,3) = 5.0; 
  m_test4(2,0) = 2.0; m_test4(2,1) = 2.0; m_test4(2,2) = 4.0; m_test4(2,3) = 3.0; 
  m_test4(3,0) =-1.0; m_test4(3,1) = 5.0; m_test4(3,2) = 3.0; m_test4(3,3) = -2.0;
  mat<double,mat_structure::rectangular> m_test4_R = m_test4;
  m_test4_Q = mat<double,mat_structure::identity>(4);
  BOOST_CHECK_NO_THROW( (detail::symmetric_QRalg_impl(m_test4_R,&m_test4_Q,1e-6)) );
  BOOST_CHECK( is_identity_mat(m_test4_Q * transpose(m_test4_Q), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_upper_triangular(m_test4_R, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat((m_test4_Q * m_test4_R * transpose_view(m_test4_Q) - m_test4), 1e-6) );
  
  
  mat<double,mat_structure::square> m_hess_Q( mat<double,mat_structure::identity>(4) );
  mat<double,mat_structure::square> m_hess_H( m_test4 );
  BOOST_CHECK_NO_THROW( (decompose_Hess(m_test4, m_hess_Q, m_hess_H, 1E-6)) );
  BOOST_CHECK( is_identity_mat(m_hess_Q * transpose(m_hess_Q), 1e-8) );
  BOOST_CHECK( is_upper_hessenberg(m_hess_H, 1e-8) );
  // A = Q T Q^T
  BOOST_CHECK( is_null_mat((m_hess_Q * m_hess_H * transpose_view(m_hess_Q) - m_test4), 1e-6) );
  
  
  mat<double,mat_structure::square> m_realshur_Q( mat<double,mat_structure::identity>(4) );
  mat<double,mat_structure::square> m_realshur_T( m_test4 );
  BOOST_CHECK_NO_THROW( (decompose_RealSchur(m_test4, m_realshur_Q, m_realshur_T, 1E-6)) );
  BOOST_CHECK( is_identity_mat(m_realshur_Q * transpose(m_realshur_Q), 1e-8) );
  BOOST_CHECK( is_upper_hessenberg(m_realshur_T, 1e-8) );
  // A = Q T Q^T
  BOOST_CHECK( is_null_mat((m_realshur_Q * m_realshur_T * transpose(m_realshur_Q) - m_test4), 1e-6) );
  
  
  mat<double,mat_structure::square> m_test_L(3,0.0);
  mat<double,mat_structure::square> m_test_D(3,0.0);
  mat<double,mat_structure::square> m_test_sqr(3,0.0);
  mat<double,mat_structure::square> m_test_inv(3,0.0);
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  BOOST_CHECK_NO_THROW( (decompose_LDL(m_test_sqr,m_test_L,m_test_D,1e-6)) );
  BOOST_CHECK( is_lower_triangular(m_test_L, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_diagonal(m_test_D, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat(((m_test_L * m_test_D * transpose_view(m_test_L)) - m_test_sqr), 1e-6) );
  
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 0.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 0.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  BOOST_CHECK_NO_THROW( (decompose_BandCholesky(m_test_sqr,m_test_L,1,1e-6)) );
  BOOST_CHECK( is_lower_triangular(m_test_L, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat(((m_test_L * transpose_view(m_test_L)) - m_test_sqr), 1e-6) );
  
  m_test_D = mat<double,mat_structure::identity>(3);
  m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 0.0; 
  m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
  m_test_sqr(2,0) = 0.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
  BOOST_CHECK_NO_THROW( (decompose_TriDiagLDL(m_test_sqr,m_test_L,m_test_D,1e-6)) );
  BOOST_CHECK( is_lower_triangular(m_test_L, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_diagonal(m_test_D, std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_null_mat(((m_test_L * m_test_D * transpose_view(m_test_L)) - m_test_sqr), 1e-6) );
  
};



BOOST_AUTO_TEST_CASE( mat_ctrl_reduction_tests )
{
  
  using namespace ReaK;
  
  std::size_t N = 6;
  std::size_t M = 2;
  
  mat<double,mat_structure::rectangular> A(N,N);
  mat<double,mat_structure::rectangular> B(N,M);
  A(0,0) = 0.0; A(0,1) = 0.0; A(0,2) = 0.0; A(0,3) = 1.0; A(0,4) = 0.0; A(0,5) = 0.0;
  A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.0; A(1,3) = 0.0; A(1,4) = 1.0; A(1,5) = 0.0; 
  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0; A(2,4) = 0.0; A(2,5) = 1.0; 
  A(3,0) = 0.0; A(3,1) = 0.0; A(3,2) = 0.0; A(3,3) = 0.0; A(3,4) = 0.0; A(3,5) = 0.0; 
  A(4,0) = 0.0; A(4,1) = 0.0; A(4,2) = 0.0; A(4,3) = 0.0; A(4,4) = 0.0; A(4,5) = 0.0; 
  A(5,0) = 0.0; A(5,1) = 0.0; A(5,2) = 0.0; A(5,3) = 0.0; A(5,4) = 0.0; A(5,5) = 0.0; 
  
  B(0,0) = 0.0; B(0,1) = 0.0;
  B(1,0) = 0.0; B(1,1) = 0.0;
  B(2,0) = 0.0; B(2,1) = 0.0;
  B(3,0) = 0.0; B(3,1) = 0.0;
  B(4,0) = 1.0; B(4,1) = 0.0;
  B(5,0) = 0.0; B(5,1) = 1.0;
  
  mat<double,mat_structure::rectangular> Ar = A;
  mat<double,mat_structure::rectangular> Br = B;
  
  mat<double,mat_structure::square> Q( (mat<double,mat_structure::identity>(N)) );
  mat<double,mat_structure::square> Z( (mat<double,mat_structure::identity>(M)) );
  std::size_t num_rank = N;
  BOOST_CHECK_NO_THROW( num_rank = ctrl_reduction(Ar, Br, Q, Z, double(1E-15)) );
  
  // check basic structures, Q and Z being orthogonal, and equality is preserved (A = Q Ar Qt and B = Q Br Zt).
  BOOST_CHECK( is_identity_mat(Q * transpose(Q), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( is_identity_mat(Z * transpose(Z), 4.0 * std::numeric_limits<double>::epsilon()) );
  BOOST_CHECK( ( is_null_mat((Q * Ar * transpose(Q) - A), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat((Q * Br * transpose(Z) - B), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  
  // Check the reduced controllability matrix:
  mat<double,mat_structure::rectangular> Cr(N,M*N);
  sub(Cr)(range(0,N),range(0,2))   = Br;
  sub(Cr)(range(0,N),range(2,4))   = Ar * Br;
  sub(Cr)(range(0,N),range(4,6))   = Ar * Ar * Br;
  sub(Cr)(range(0,N),range(6,8))   = Ar * Ar * Ar * Br;
  sub(Cr)(range(0,N),range(8,10))   = Ar * Ar * Ar * Ar * Br;
  sub(Cr)(range(0,N),range(10,12)) = Ar * Ar * Ar * Ar * Ar * Br;
  BOOST_CHECK_EQUAL( num_rank, 4 );
  BOOST_CHECK( ( is_null_mat((sub(Cr)(range(num_rank,6),range(0,12))), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  
  /* Use row-rank SVD to check that the controllable parts are full row-rank. */
  mat<double,mat_structure::diagonal> Er(num_rank);
  mat<double,mat_structure::rectangular> Ur(num_rank, num_rank);
  mat<double,mat_structure::rectangular> Vr(M*N, num_rank);
  BOOST_CHECK_NO_THROW( decompose_SVD(sub(Cr)(range(0,num_rank),range(0,12)),Ur,Er,Vr,double(1E-15)) );
  BOOST_CHECK_EQUAL( ( numrank_SVD(Er, double(1E-15)) ), num_rank );
  
};




BOOST_AUTO_TEST_CASE( mat_balancing_tests )
{
  using namespace ReaK;
  
  mat<double,mat_structure::square> A(3);
  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 1e-4;
  A(1,0) = 1.0; A(1,1) = 1.0; A(1,2) = 1e-2;
  A(2,0) = 1e4; A(2,1) = 1e2; A(2,2) = 1.0;
  
  mat<double,mat_structure::square> A_bal = A;
  vect_n<int> D_bal(3);
  BOOST_CHECK_NO_THROW( balance(A_bal, D_bal) );

//   std::cout << " A = " << A << std::endl;
  mat<double,mat_structure::diagonal> A_E(3);
  mat<double,mat_structure::square> A_U(3);
  mat<double,mat_structure::square> A_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(A, A_U, A_E, A_V, double(1E-15)) );
  double cond_before = condition_number_SVD(A_E);
//   std::cout << " cond(A) = " << cond_before << std::endl;
  
//   std::cout << " A_bal = " << A_bal << std::endl;
  mat<double,mat_structure::diagonal> A_bal_E(3);
  mat<double,mat_structure::square> A_bal_U(3);
  mat<double,mat_structure::square> A_bal_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(A_bal, A_bal_U, A_bal_E, A_bal_V, double(1E-15)) );
  double cond_after = condition_number_SVD(A_bal_E);
//   std::cout << " cond(A_bal) = " << cond_after << std::endl;
  
  // The condition number must have been reduced:
  BOOST_CHECK( cond_after < cond_before );
  
  // Condition number must be with the N-th power of 2 (at most the difference between 
  // any two eigen-value is a factor of 2, so, at most, the difference between the largest
  // and smallest eigen-value is 2^N. Here, N = 3.
  BOOST_CHECK( cond_after < std::ldexp(1,3) );
  
  apply_left_bal_exp(D_bal, A_bal);
  apply_right_bal_inv_exp(A_bal, D_bal);
//   std::cout << " A_restored = " << A_bal << std::endl;
  
  // must be exact because exact floating-point arithmetic is used in balancing:
  BOOST_CHECK( ( is_null_mat((A - A_bal), 1e-200) ) );  
  
  
  
  mat<double,mat_structure::square> B(3);
  B(0,0) = 1.0; B(0,1) = 1e-2; B(0,2) = 1e-4;
  B(1,0) = 1e3; B(1,1) = 1.0; B(1,2) = 1.0;
  B(2,0) = 1e4; B(2,1) = 0.0; B(2,2) = 2.0;
  
  
  mat<double,mat_structure::square> A_penbal = A;
  mat<double,mat_structure::square> B_penbal = B;
  vect_n<int> Dl_penbal(3);
  vect_n<int> Dr_penbal(3);
  BOOST_CHECK_NO_THROW( balance_pencil(A_penbal, B_penbal, Dl_penbal, Dr_penbal) );
  
//   std::cout << " A = " << A << std::endl;
//   std::cout << " cond(A) = " << cond_before << std::endl;
//   std::cout << " B = " << B << std::endl;
  mat<double,mat_structure::diagonal> B_E(3);
  mat<double,mat_structure::square> B_U(3);
  mat<double,mat_structure::square> B_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(B, B_U, B_E, B_V, double(1E-15)) );
  double cond_B_before = condition_number_SVD(B_E);
//   std::cout << " cond(B) = " << cond_B_before << std::endl;
  
  
  
  mat<double, mat_structure::square> M_pen(3);
  for(std::size_t i = 0; i < 3; ++i)
    for(std::size_t j = 0; j < 3; ++j)
      M_pen(i,j) = A(i,j) * A(i,j) + B(i,j) * B(i,j);
//   std::cout << " M_pen = " << M_pen << std::endl;
  mat<double,mat_structure::diagonal> M_pen_E(3);
  mat<double,mat_structure::square> M_pen_U(3);
  mat<double,mat_structure::square> M_pen_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(M_pen, M_pen_U, M_pen_E, M_pen_V, double(1E-15)) );
  double cond_M_before = condition_number_SVD(M_pen_E);
//   std::cout << " cond(M_pen) = " << cond_M_before << std::endl;
  
  
//   std::cout << " A_penbal = " << A_penbal << std::endl;
  mat<double,mat_structure::diagonal> A_penbal_E(3);
  mat<double,mat_structure::square> A_penbal_U(3);
  mat<double,mat_structure::square> A_penbal_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(A_penbal, A_penbal_U, A_penbal_E, A_penbal_V, double(1E-15)) );
  double cond_A_penbal = condition_number_SVD(A_penbal_E);
//   std::cout << " cond(A_penbal) = " << cond_A_penbal << std::endl;
  
//   std::cout << " B_penbal = " << B_penbal << std::endl;
  mat<double,mat_structure::diagonal> B_penbal_E(3);
  mat<double,mat_structure::square> B_penbal_U(3);
  mat<double,mat_structure::square> B_penbal_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(B_penbal, B_penbal_U, B_penbal_E, B_penbal_V, double(1E-15)) );
  double cond_B_penbal = condition_number_SVD(B_penbal_E);
//   std::cout << " cond(B_penbal) = " << cond_B_penbal << std::endl;
  
  
  mat<double, mat_structure::square> M_penbal(3);
  for(std::size_t i = 0; i < 3; ++i)
    for(std::size_t j = 0; j < 3; ++j)
      M_penbal(i,j) = A_penbal(i,j) * A_penbal(i,j) + B_penbal(i,j) * B_penbal(i,j);
//   std::cout << " M_penbal = " << M_penbal << std::endl;
  mat<double,mat_structure::diagonal> M_penbal_E(3);
  mat<double,mat_structure::square> M_penbal_U(3);
  mat<double,mat_structure::square> M_penbal_V(3);
  BOOST_CHECK_NO_THROW( decompose_SVD(M_penbal, M_penbal_U, M_penbal_E, M_penbal_V, double(1E-15)) );
  double cond_M_after = condition_number_SVD(M_penbal_E);
//   std::cout << " cond(M_penbal) = " << cond_M_after << std::endl;
  
  
  // (A,B)_balanced = Dl^-1 (A,B) Dr
  
  apply_left_bal_exp(Dl_penbal, A_penbal);
  apply_left_bal_exp(Dl_penbal, B_penbal);
  apply_right_bal_inv_exp(A_penbal, Dr_penbal);
  apply_right_bal_inv_exp(B_penbal, Dr_penbal);
//   std::cout << " A_penbal_restored = " << A_penbal << std::endl;
//   std::cout << " B_penbal_restored = " << B_penbal << std::endl;
  
  
  // The condition numbers must have been reduced:
  BOOST_CHECK( cond_A_penbal < cond_before );
  BOOST_CHECK( cond_B_penbal < cond_B_before );
  BOOST_CHECK( cond_M_after < cond_M_before );
  
  BOOST_CHECK( cond_A_penbal < 150.0 );
  BOOST_CHECK( cond_B_penbal < 10.0 );
  BOOST_CHECK( cond_M_after < 600.0 );
  
  // must be exact because exact floating-point arithmetic is used in balancing:
  BOOST_CHECK( ( is_null_mat((A - A_penbal), 1e-200) ) );  
  BOOST_CHECK( ( is_null_mat((B - B_penbal), 1e-200) ) );  
  
  
};








