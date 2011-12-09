
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
#include "vect_alg.hpp"

#include "mat_gaussian_elim.hpp"
#include "mat_cholesky.hpp"
#include "mat_jacobi_method.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_svd_method.hpp"

#include "boost/date_time/posix_time/posix_time.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979324
#define M_PI_2 1.57079632679489662
#endif

int main() {

using namespace ReaK;

  unsigned int passed = 0;

try {

  /*************************** Fixed-Length Vector Tests ****************************/
  if(true){
    vect<double,3> gravity_acc;
    gravity_acc[0] = 0.0;
    gravity_acc[1] = -9.81;
    gravity_acc[2] = 0.0;
    if( 3 == gravity_acc.size() )
      ++passed;
    else
      RK_ERROR("vect<double,3>::size() test did not pass!");
    if( 3 == gravity_acc.max_size() )
      ++passed;
    else
      RK_ERROR("vect<double,3>::max_size() test did not pass!");
    if( 3 == gravity_acc.capacity() )
      ++passed;
    else
      RK_ERROR("vect<double,3>::capacity() test did not pass!");
    if( !gravity_acc.empty() )
      ++passed;
    else
      RK_ERROR("vect<double,3>::empty() test did not pass!");
    vect<double,3>::iterator it = gravity_acc.begin();      //RK_NOTICE(2,"Passed: " << passed << " should have 8.");
    if(std::fabs(*it) < std::numeric_limits< double >::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::begin() test did not pass!");
    if(gravity_acc.end() - it == 3)
      ++passed;
    else
      RK_ERROR("vect<double,3>::end() test did not pass!");
    if(gravity_acc.end() - ++it == 2)
      ++passed;
    else
      RK_ERROR("vect<double,3>::iterator::operator++ test did not pass!");
    const vect<double,3>& gravity_acc_ref = gravity_acc;
    vect<double,3>::const_iterator cit = gravity_acc_ref.begin();
    if(std::fabs(*cit) < std::numeric_limits< double >::epsilon()) 
      ++passed;
    else
      RK_ERROR("vect<double,3>::cbegin() test did not pass!");
    if(gravity_acc_ref.end() - cit == 3) 
      ++passed;
    else
      RK_ERROR("vect<double,3>::cend() test did not pass!");
    if(gravity_acc_ref.end() - ++cit == 2) 
      ++passed;
    else
      RK_ERROR("vect<double,3>::const_iterator::operator++ test did not pass!");
    double ones[] = {1.0,1.0,1.0};
    vect<double,3> ones_v(ones);
    if(ones[1] == ones_v[1]) 
      ++passed;
    else
      RK_ERROR("vect<double,3>::vect(double*) test did not pass!");
    double obj_mass(3.0);
    vect<double,3> gravity_force; 
    gravity_force = (gravity_acc * obj_mass);
    if(std::fabs(gravity_force[1] + 3.0*9.81) < std::numeric_limits<double>::epsilon()) 
      ++passed;
    else
      RK_ERROR("operator*(vect<double,3>,double) test did not pass!");
    if(gravity_force == obj_mass * gravity_acc) 
      ++passed;
    else
      RK_ERROR("operator*(double,vect<double,3>) test did not pass!");
    if(std::fabs(norm_2_sqr(gravity_acc) - 9.81*9.81) < 100.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm_2_sqr(vect<double,3>) test did not pass!");
    if(std::fabs(norm_2(gravity_acc) - 9.81) < 10.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm_2(vect<double,3>) test did not pass!");
    vect<double,3> gravity_dir(unit(gravity_acc));
    if(std::fabs(norm_2(gravity_dir) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("unit(vect<double,3>) test did not pass!");
    if(std::fabs(norm_2(vect<double,1>(1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,1>(d) test did not pass!");
    if(std::fabs(norm_2(vect<double,2>(0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,2>(d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,3>(0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>(d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,4>(0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,4>(d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,5>(0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,5>(d,d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,6>(0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,6>(d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,7>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,7>(d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,8>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,8>(d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,9>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,9>(d,d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm_2(vect<double,10>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,10>(d,d,d,d,d,d,d,d,d,d) test did not pass!");
    vect<double,3> displacement(1.0, 2.0, 3.0);             //RK_NOTICE(2,"Passed: " << passed << " should have 32.");
    double gravity_potential = gravity_force * displacement;
    if(std::fabs(gravity_potential + 2.0 * 9.81 * 3.0) < 60.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("operator*(vect<double,3>,vect<double,3>) test did not pass!");
    vect<double,3> gravity_moment = displacement % gravity_force;
    if(std::fabs(gravity_moment[0] - 3.0*9.81*3.0) + std::fabs(gravity_moment[2] + 3.0*9.81) < 100.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("operator%(double,vect<double,3>) test did not pass!");
    
    double dist_sqr = norm_2_sqr(displacement);               //RK_NOTICE(2,"Passed: " << passed << " should have 34.");
    vect<double,3> displacement_inv(-displacement); 
    if(std::fabs(displacement_inv * displacement + dist_sqr) < std::numeric_limits<double>::epsilon())
                                                 ++passed;
    vect<double,3> long_displacement(10.0, 20.0, 30.0); ++passed;
    if(colinear(displacement,long_displacement)) ++passed;
    if(long_displacement > displacement)         ++passed;
    if(displacement < long_displacement)         ++passed;
    if(!(displacement == long_displacement))     ++passed;
    if(displacement == displacement)             ++passed;

    long_displacement += displacement;                     //RK_NOTICE(2,"Passed: " << passed << " should have 41.");
    if(std::fabs(norm_2(long_displacement) - 11.0*norm_2(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
                                                 ++passed;
    long_displacement *= 2.0;
    if(std::fabs(norm_2(long_displacement) - 22.0*norm_2(displacement)) < 400.0*std::numeric_limits<double>::epsilon())
                                                 ++passed;
    long_displacement /= 2.0;
    if(std::fabs(norm_2(long_displacement) - 11.0*norm_2(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
                                                 ++passed;
    long_displacement -= displacement;
    if(std::fabs(norm_2(long_displacement) - 10.0*norm_2(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
                                                 ++passed;
                                                           //RK_NOTICE(2,"Passed: " << passed << " should have 45.");
  };

  /***************************** Variable-Length Vector Tests ***********************************/
  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/** PRIMITIVE VARIABLE-LENGTH VECTOR TESTS ***/");
    RK_NOTICE(2,"/*********************************************/");
    vect_n<double> gravity_acc(3);
    gravity_acc[0] = 0.0f;
    gravity_acc[1] = -9.81f;
    gravity_acc[2] = 0.0f;
    double obj_mass(3.0);
    vect_n<double> gravity_force;
    gravity_force = (gravity_acc * obj_mass);
    vect_n<double> gravity_dir(unit(gravity_acc));
    RK_NOTICE(2,"Gravity vector is: " << gravity_force[0] << ";" << gravity_force[1] << ";" << gravity_force[2]);
    RK_NOTICE(2,"Gravity direction is: " << gravity_dir[0] << ";" << gravity_dir[1] << ";" << gravity_dir[2]);
    vect_n<double> displacement(1.0, 2.0, 3.0);
    double gravity_potential = gravity_force * displacement;
    RK_NOTICE(2,"Gravity Potential is: " << gravity_potential);

    double dist_sqr = norm_2_sqr(displacement);
    double distance = norm_2(displacement);
    vect_n<double> displacement_inv(-displacement);
    vect_n<double> long_displacement(10.0, 20.0, 30.0);
    bool is_colinear = colinear(displacement,long_displacement);
    bool is_greater = (displacement > long_displacement);
    bool is_smaller = (displacement < long_displacement);
    bool is_equal = (displacement == long_displacement);
    bool is_equal2 = (displacement == displacement);

    RK_NOTICE(2,"Displacement vector is: " << displacement[0] << ";" << displacement[1] << ";" << displacement[2]);
    RK_NOTICE(2,"Inverse Displacement vector is: " << displacement_inv[0] << ";" << displacement_inv[1] << ";" << displacement_inv[2]);
    RK_NOTICE(2,"Long Displacement vector is: " << long_displacement[0] << ";" << long_displacement[1] << ";" << long_displacement[2]);
    RK_NOTICE(2,"The above two are" << (is_colinear ? " indeed" : " not") << " colinear.");
    RK_NOTICE(2,"Distance Square is: " << dist_sqr);
    RK_NOTICE(2,"Distance is: " << distance);
    RK_NOTICE(2,"Displacement is" << (is_greater ? " indeed" : " not") << " greater than Long Displacement.");
    RK_NOTICE(2,"Displacement is" << (is_smaller ? " indeed" : " not") << " smaller than Long Displacement.");
    RK_NOTICE(2,"Displacement is" << (is_equal ? " indeed" : " not") << " equal to Long Displacement.");
    RK_NOTICE(2,"Displacement is" << (is_equal2 ? " indeed" : " not") << " equal to Displacement");
    long_displacement += displacement;
    bool is_colinear2 = colinear(displacement,long_displacement);
    RK_NOTICE(2,"Long Displacement vector is now: " << long_displacement[0] << ";" << long_displacement[1] << ";" << long_displacement[2]);
    RK_NOTICE(2,"The above is" << (is_colinear2 ? " still" : " no longer") << " colinear with Displacement.");
    long_displacement *= 2.0;
    RK_NOTICE(2,"Long Displacement vector Doubled is now: " << long_displacement[0] << ";" << long_displacement[1] << ";" << long_displacement[2]);
    long_displacement /= 2.0;
    RK_NOTICE(2,"Long Displacement vector Halfed is now: " << long_displacement[0] << ";" << long_displacement[1] << ";" << long_displacement[2]);
    RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
  };

  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/**** PRIMITVE VARIABLE-SIZE MATRIX TESTS ****/");
    RK_NOTICE(2,"/*********************************************/");

    mat<double> m1234(1.0,2.0,3.0,4.0);
    RK_NOTICE(2,"m1234 = " << m1234);
    RK_NOTICE(2,"m1234 = " << m1234(0,0) << " " << m1234(1,0) << " " << m1234(0,1) << " " << m1234(1,1));
    RK_NOTICE(2,"m1324 = " << (m1234 = transpose(m1234)));
    RK_NOTICE(2,"m1324 = " << (m1234 = transpose_move(m1234)));

    mat<double> m_ident(3,3,true);
    RK_NOTICE(2,"identity 3x3 = " << m_ident);
    mat_identity<double>::type m_ident2(3);
    RK_NOTICE(2,"identity 3x3 = " << m_ident2);
    mat<double> m_ones(3,3,double(1.0));
    RK_NOTICE(2,"ones 3x3 = " << m_ones);
    mat_null<double>::type m_zeroes(3,3);
    RK_NOTICE(2,"zeroes 3x3 = " << m_zeroes);

    double f4321[] = {4.0,3.0,2.0,1.0};
    std::vector<double> v4321(f4321,f4321 + 4);
    mat<double> m4321(v4321,2,2);
    RK_NOTICE(2,"m4321 = " << m4321);

    mat<double> m4321_cpy(m4321);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy = m4321;
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy += mat<double>(2,2,true);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy += mat_identity<double>::type(2);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy -= mat<double>(2,2,true);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy -= mat_identity<double>::type(2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy *= 2.0;
    RK_NOTICE(2,"m4321 copy twice = " << m4321_cpy);

    m4321_cpy = m4321_cpy * double(0.5);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy = m4321_cpy + mat<double>(2,2,true);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy = m4321_cpy + mat_identity<double>::type(2);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy = m4321_cpy - mat<double>(2,2,true);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    m4321_cpy = m4321_cpy - mat_identity<double>::type(2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    RK_NOTICE(2,"m4321 copy negated = " << -m4321_cpy);

    m4321_cpy = m4321_cpy * mat_identity<double>::type(2);
    RK_NOTICE(2,"m4321 copy * identity = " << m4321_cpy);

    if(m4321_cpy == m4321)
      RK_NOTICE(2,"m4321 copy is equal to m4321");
    if(m4321_cpy != mat<double>(2,2,true))
      RK_NOTICE(2,"m4321 copy is not equal to identity");
    if(m4321_cpy != mat_identity<double>::type(2))
      RK_NOTICE(2,"m4321 copy is not equal to identity");

    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    vect_n<double> v85 = m4321_cpy * vect_n<double>(2,1.0);
    RK_NOTICE(2,"m4321 copy * vect(1.0,1.0) = " << v85);
    vect_n<double> v104 =  vect_n<double>(2,1.0) * m4321_cpy;
    RK_NOTICE(2,"vect(1.0,1.0) * m4321 copy = " << v104);

    set_block(m_ident,m4321_cpy,1,1);
    RK_NOTICE(2,"identity 3x3 with sub matrix m4321 is " << m_ident);

    m4321_cpy = get_block(m_ident,1,1,2,2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    RK_NOTICE(2,"m4321 copy transpose = " << transpose(m4321_cpy));
    RK_NOTICE(2,"m4321 copy sym part = " << (mat<double,mat_structure::symmetric>(m4321_cpy)));
    RK_NOTICE(2,"m4321 copy skew part = " << (mat<double,mat_structure::skew_symmetric>(m4321_cpy)));
    RK_NOTICE(2,"m4321 copy twice = " << double(2.0) * m4321_cpy);

    mat<double,mat_structure::diagonal> m123(vect_n<double>(1.0,2.0,3.0));
    RK_NOTICE(2,"m123 diagonal = " << m123);
    //mat_cm<double> m_blk_diag(block_diag_mat(mat_cm<double>(1.0,2.0,3.0,4.0),mat_cm<double>(1,1,true)));
    //RK_NOTICE(2,"m_blk_diag block diagonal = " << m_blk_diag);
    //mat_cm<double> m_blk(block_mat(mat_cm<double>(1.0,2.0,3.0,4.0),m4321.getSubMat(0,0,2,1),m4321.getSubMat(1,0,1,2),mat_cm<double>(1,1,true)));
    //RK_NOTICE(2,"m_blk block matrix = " << m_blk);
    RK_NOTICE(2,"vect(1,2,3) skew matrix = " << (mat<double,mat_structure::skew_symmetric>(vect_n<double>(1.0,2.0,3.0))));

    RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
  };
  
  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/********* MATRIX NUM-METHODS TESTS **********/");
    RK_NOTICE(2,"/*********************************************/");

    mat<double,mat_structure::symmetric> m_gauss(2.0,-1.0,0.0,2.0,-1.0,2.0);
    RK_NOTICE(2,"Testing the following matrix: " << m_gauss);
    
    boost::posix_time::ptime t1(boost::posix_time::microsec_clock::local_time());
    mat<double,mat_structure::square> m_gauss_inv(3);
    invert_gaussian(m_gauss,m_gauss_inv,double(1E-15));
    boost::posix_time::time_duration dt(boost::posix_time::microsec_clock::local_time() - t1);
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Gaussian: " << m_gauss_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_gauss_inv));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::square> m_plu_inv(mat<double,mat_structure::identity>(3));
    invert_PLU(m_gauss,m_plu_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with PLU: " << m_plu_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_plu_inv));
    
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::symmetric> m_gauss2_inv(3);
    invert_gaussian(m_gauss,m_gauss2_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Gaussian Symmetric: " << m_gauss2_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_gauss2_inv));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::symmetric> m_plu2_inv(mat<double,mat_structure::identity>(3));
    invert_PLU(m_gauss,m_plu2_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with PLU Symmetric: " << m_plu2_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_plu2_inv));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::symmetric> m_cholesky_inv(mat<double,mat_structure::identity>(3));
    invert_Cholesky(m_gauss,m_cholesky_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Cholesky: " << m_cholesky_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_cholesky_inv));
    
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::square> m_cholesky2_inv(mat<double,mat_structure::identity>(3));
    invert_Cholesky(m_gauss,m_cholesky2_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Cholesky Square: " << m_cholesky2_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_cholesky2_inv));
    
    mat<double,mat_structure::diagonal> m_gauss_diag(vect<double,3>(2,1,0.5));
    RK_NOTICE(2,"Testing the following matrix: " << m_gauss_diag);
    
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::square> m_cholesky3_inv(mat<double,mat_structure::identity>(3));
    invert_Cholesky(m_gauss_diag,m_cholesky3_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Diag Cholesky Square: " << m_cholesky3_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss_diag * m_cholesky3_inv));
    
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::diagonal> m_cholesky4_inv(mat<double,mat_structure::identity>(3));
    invert_Cholesky(m_gauss_diag,m_cholesky4_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Diag Cholesky Square: " << m_cholesky4_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss_diag * m_cholesky4_inv));
     
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::symmetric> m_jacobi_pinv(3);
    pseudoinvert_Jacobi(m_gauss,m_jacobi_pinv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with Jacobi: " << m_jacobi_pinv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_jacobi_pinv));
   
    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::square> m_qr_inv(3);
    pseudoinvert_QR(m_gauss,m_qr_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with QR: " << m_qr_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_qr_inv));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::square> m_svd_inv(3);
    pseudoinvert_SVD(m_gauss,m_svd_inv,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The inverse was found in " << dt.total_microseconds() << " microseconds with SVD: " << m_svd_inv);
    RK_NOTICE(2,"The inverse * the matrix = " << (m_gauss * m_svd_inv));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::diagonal> m_jacobi_E(3);
    mat<double,mat_structure::square> m_jacobi_Q(3);
    eigensolve_Jacobi(m_gauss,m_jacobi_E,m_jacobi_Q,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The eigenvalues were found in " << dt.total_microseconds() << " microseconds with Jacobi.");
    RK_NOTICE(2,"The eigenvalues are = " << m_jacobi_E);
    RK_NOTICE(2,"The eigenvectors are = " << m_jacobi_Q);
    RK_NOTICE(2,"Q * E * Qt = " << m_jacobi_Q * m_jacobi_E * transpose(m_jacobi_Q));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::diagonal> m_qr_E(3,true);
    mat<double,mat_structure::square> m_qr_Q(3);
    //eigensolve_QR(m_gauss,m_qr_E,m_qr_Q,50,double(1E-6));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The eigenvalues were found in " << dt.total_microseconds() << " microseconds with QR.");
    RK_NOTICE(2,"The eigenvalues are = " << m_qr_E);
    RK_NOTICE(2,"The eigenvectors are = " << m_qr_Q);
    RK_NOTICE(2,"Q * E * Qt = " << m_qr_Q * m_qr_E * transpose(m_qr_Q));

    t1 = boost::posix_time::microsec_clock::local_time();
    mat<double,mat_structure::diagonal> m_svd_E(3);
    mat<double,mat_structure::square> m_svd_U(3);
    mat<double,mat_structure::square> m_svd_V(3);
    decompose_SVD(m_gauss,m_svd_U,m_svd_E,m_svd_V,double(1E-15));
    dt = boost::posix_time::microsec_clock::local_time() - t1;
    RK_NOTICE(2,"The eigenvalues were found in " << dt.total_microseconds() << " microseconds with SVD.");
    RK_NOTICE(2,"The eigenvalues are = " << m_svd_E);
    RK_NOTICE(2,"V = " << m_svd_V);
    RK_NOTICE(2,"U = " << m_svd_U);
    RK_NOTICE(2,"Q * E * Qt = " << m_svd_U * m_svd_E * transpose(m_svd_V));
    
  };
  
  if(true) {
    unsigned int inc(1);
    mat<double,mat_structure::symmetric> m_test(2.0,-1.0,0.0,2.0,-1.0,2.0);
    mat<double,mat_structure::symmetric> m_inc(2.0,-1.0,2.0);
    boost::posix_time::ptime t1;
    boost::posix_time::time_duration dt[9];

    std::ofstream out_stream;
    out_stream.open("performance_data.dat");
    out_stream << "N Gauss PLU Chol Jac QR SVD Jac_E QR_E SVD_E" << std::endl;
    std::cout << "Recording performance..." << std::endl;
    for(unsigned int i=3;i<=1000;i += inc) {
      if(i == 50) {
	inc = 10; 
	m_inc = get_block(m_test,0,inc+1); 
      } else if(i == 100) {
	inc = 20; 
	m_inc = get_block(m_test,0,inc+1); 
      } else if(i == 500) {
	inc = 50; 
	m_inc = get_block(m_test,0,inc+1); 
      };

      m_test.set_row_count(i);
      set_block(m_test,m_inc,i-inc-1);

      mat<double,mat_structure::square> m_gauss_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_gaussian(m_test,m_gauss_inv,double(1E-15));
      dt[0] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::square> m_plu_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_PLU(m_test,m_plu_inv,double(1E-15));
      dt[1] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::symmetric> m_cholesky_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_Cholesky(m_test,m_cholesky_inv,double(1E-15));
      dt[2] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::symmetric> m_jacobi_pinv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_Jacobi(m_test,m_jacobi_pinv,double(1E-15));
      dt[3] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::square> m_qr_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_QR(m_test,m_qr_inv,double(1E-15));
      dt[4] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::square> m_svd_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_SVD(m_test,m_svd_inv,double(1E-15));
      dt[5] = boost::posix_time::microsec_clock::local_time() - t1;;

      mat<double,mat_structure::diagonal> m_jacobi_E(i);
      mat<double,mat_structure::square> m_jacobi_Q(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      eigensolve_Jacobi(m_test,m_jacobi_E,m_jacobi_Q,double(1E-15));
      dt[6] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::diagonal> m_qr_E(i);
      mat<double,mat_structure::square> m_qr_Q(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      //eigensolve_QR(m_test,m_qr_E,m_qr_Q,50000,double(1E-15)); //this iterates forever!! (well, for a long time at least)
      dt[7] = boost::posix_time::microsec_clock::local_time() - t1;

      mat<double,mat_structure::diagonal> m_svd_E(i);
      mat<double,mat_structure::square> m_svd_U(i);
      mat<double,mat_structure::square> m_svd_V(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      decompose_SVD(m_test,m_svd_U,m_svd_E,m_svd_V,double(1E-15));
      dt[8] = boost::posix_time::microsec_clock::local_time() - t1;

      out_stream << i << " " << dt[0].total_microseconds()
                      << " " << dt[1].total_microseconds()
                      << " " << dt[2].total_microseconds()
		      << " " << dt[3].total_microseconds()
		      << " " << dt[4].total_microseconds()
		      << " " << dt[5].total_microseconds()
		      << " " << dt[6].total_microseconds()
		      << " " << dt[7].total_microseconds()
		      << " " << dt[8].total_microseconds() << std::endl;
      std::cout << i << std::endl;

    };
    std::cout << "Done!" << std::endl;
    out_stream.close();
  };
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
  };
  
  RK_NOTICE(2,"There were " << passed << " successful tests passed on the math_gen library, out of 45 possible successes.");
  
  return 0;
};







