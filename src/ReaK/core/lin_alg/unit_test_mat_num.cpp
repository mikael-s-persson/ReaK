
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

#include "mat_schur_decomp.hpp"

#include "mat_ctrl_decomp.hpp"

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE mat_num
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( mat_inversion_tests )
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
  
#if 0
  mat<double,mat_structure::symmetric> m_jacobi_pinv(3);
  BOOST_CHECK_NO_THROW( pseudoinvert_Jacobi(m_gauss,m_jacobi_pinv,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_jacobi_pinv - m_gauss_trueinv, 32.0 * std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_identity_mat((m_gauss * m_jacobi_pinv), 32.0 * std::numeric_limits<double>::epsilon()) ) );
#endif
  
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
  
#if 0
  mat<double,mat_structure::diagonal> m_jacobi_E(3);
  mat<double,mat_structure::square> m_jacobi_Q(3);
  BOOST_CHECK_NO_THROW( eigensolve_Jacobi(m_gauss,m_jacobi_E,m_jacobi_Q,double(1E-15)) );
  BOOST_CHECK( ( is_null_mat(m_jacobi_E - m_gauss_trueeig, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat((m_jacobi_Q * m_jacobi_E * transpose(m_jacobi_Q) - m_gauss), std::numeric_limits<double>::epsilon()) ) );
#endif
  
  mat<double,mat_structure::diagonal> m_qr_E(3,true);
  mat<double,mat_structure::square> m_qr_Q(3);
  BOOST_CHECK_NO_THROW( eigensolve_SymQR(m_gauss,m_qr_Q,m_qr_E,double(1E-15)) );
  vect_n<double> m_qr_Es(m_qr_E(0,0), m_qr_E(1,1), m_qr_E(2,2));
  std::sort(m_qr_Es.begin(), m_qr_Es.end());
  BOOST_CHECK( ( norm_inf(m_qr_Es - m_gauss_trueeig) < 10.0 * std::numeric_limits<double>::epsilon() ) );
  BOOST_CHECK( ( is_null_mat((m_qr_Q * m_qr_E * transpose(m_qr_Q) - m_gauss), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  
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
  
  {
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
    
    
#if 0
    // it seems that this SVD implementation cannot be used for min-norm pseudo-invert
    mat<double,mat_structure::rectangular> m_test_svd_pinv(3,2);
    
    int nu = (m_test.get_row_count() > m_test.get_col_count() ? m_test.get_col_count() : m_test.get_row_count());
    mat<double,mat_structure::square> U(m_test.get_row_count());
    mat<double,mat_structure::diagonal> E(nu);
    mat<double,mat_structure::square> V(m_test.get_col_count());
    
    decompose_SVD(m_test,U,E,V,1e-6);
    m_test_svd_pinv = V * invert(E) * transpose(U);
    
//     BOOST_CHECK_NO_THROW( pseudoinvert_SVD(m_test, m_test_svd_pinv, 1e-6) );
    mat<double,mat_structure::rectangular> m_test_svd_x = m_test_svd_pinv * m_test_minnorm_b;
    std::cout << "A x = " << (m_test * m_test_svd_x) << std::endl;
    BOOST_CHECK( is_null_mat( m_test * m_test_svd_x - m_test_minnorm_b, 1e-6) );
#endif
    
  };
  
  {
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
    
    mat<double,mat_structure::square> m_test_P(mat<double,mat_structure::identity>(3));
    m_test_sqr(0,0) = 6.0; m_test_sqr(0,1) = 3.0; m_test_sqr(0,2) = 2.0; 
    m_test_sqr(1,0) = 3.0; m_test_sqr(1,1) = 5.0; m_test_sqr(1,2) = 2.0; 
    m_test_sqr(2,0) = 2.0; m_test_sqr(2,1) = 2.0; m_test_sqr(2,2) = 4.0; 
    BOOST_CHECK_NO_THROW( (decompose_PLTLP(m_test_sqr,m_test_L,m_test_D,m_test_P,1e-6)) );
    BOOST_CHECK( is_lower_triangular(m_test_L, std::numeric_limits<double>::epsilon()) );
    // THIS TEST FAILS:
    //BOOST_CHECK( is_diagonal(m_test_D, 4.0 * std::numeric_limits<double>::epsilon()) );
    BOOST_CHECK( is_identity_mat(m_test_P * transpose(m_test_P), 4.0 * std::numeric_limits<double>::epsilon()) );
    // THIS TEST FAILS:
    //BOOST_CHECK( is_null_mat(((transpose_view(m_test_P) * m_test_L * m_test_D * transpose_view(m_test_L) * m_test_P) - m_test_sqr), 1e-6) );
    
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
    // THIS TEST FAILS:
    //BOOST_CHECK( is_null_mat(((m_test_L * m_test_D * transpose_view(m_test_L)) - m_test_sqr), 1e-6) );
    
    
  };
  
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
  sub(Cr)(range(0,N-1),range(0,1))   = Br;
  sub(Cr)(range(0,N-1),range(2,3))   = Ar * Br;
  sub(Cr)(range(0,N-1),range(4,5))   = Ar * Ar * Br;
  sub(Cr)(range(0,N-1),range(6,7))   = Ar * Ar * Ar * Br;
  sub(Cr)(range(0,N-1),range(8,9))   = Ar * Ar * Ar * Ar * Br;
  sub(Cr)(range(0,N-1),range(10,11)) = Ar * Ar * Ar * Ar * Ar * Br;
  BOOST_CHECK_EQUAL( num_rank, 4 );
  BOOST_CHECK( ( is_null_mat((sub(Cr)(range(num_rank,5),range(0,11))), 10.0 * std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::diagonal> Er(num_rank);
  mat<double,mat_structure::square> Ur(num_rank);
  mat<double,mat_structure::square> Vr(M*N);
  BOOST_CHECK_NO_THROW( decompose_SVD(sub(Cr)(range(0,num_rank-1),range(0,11)),Ur,Er,Vr,double(1E-15)) );
  BOOST_CHECK_EQUAL( ( numrank_SVD(Er, double(1E-15)) ), num_rank );
  
};









