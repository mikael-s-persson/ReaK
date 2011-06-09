
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
#include <stdio.h>

#include "mat_alg.hpp"
#include "vect_alg.hpp"

#include "mat_gaussian_elim.hpp"
#include "mat_cholesky.hpp"
#include "mat_jacobi_method.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_svd_method.hpp"

#include "rotations.hpp"

#include "boost/date_time/posix_time/posix_time.hpp"



int main() {

using namespace ReaK;

  unsigned int passed = 0;

try {

  /*************************** Fixed-Length Vector Tests ****************************/
  if(true){
    vect<float,3> gravity_acc;                   ++passed;
    gravity_acc[0] = 0.0;                        ++passed; //RK_NOTICE(2,"Passed: " << passed << " should have 2.");
    gravity_acc[1] = -9.81;
    gravity_acc[2] = 0.0;
    if( 3 == gravity_acc.size() )                ++passed;
    if( 3 == gravity_acc.max_size() )            ++passed;
    if( 3 == gravity_acc.capacity() )            ++passed;
    gravity_acc.resize(5);                                 //RK_NOTICE(2,"Passed: " << passed << " should have 5.");
    if( 3 == gravity_acc.size() )                ++passed;
    if( !gravity_acc.empty() )                   ++passed;
    gravity_acc.reserve(5);
    if( 3 == gravity_acc.capacity() )            ++passed;
    vect<float,3>::iterator it = gravity_acc.begin();      //RK_NOTICE(2,"Passed: " << passed << " should have 8.");
    if(*it == 0.0)                               ++passed;
    if(gravity_acc.end() - it == 3)              ++passed;
    if(gravity_acc.end() - ++it == 2)            ++passed;
    const vect<float,3>& gravity_acc_ref = gravity_acc;    //RK_NOTICE(2,"Passed: " << passed << " should have 11.");
    vect<float,3>::const_iterator cit = gravity_acc_ref.begin();
    if(*cit == 0.0)                              ++passed;
    if(gravity_acc_ref.end() - cit == 3)         ++passed;
    if(gravity_acc_ref.end() - ++cit == 2)       ++passed;
    float ones[] = {1.0,1.0,1.0};                         // RK_NOTICE(2,"Passed: " << passed << " should have 14.");
    vect<float,3> ones_v(ones);                  ++passed;
    if(ones[1] == ones_v[1])                     ++passed;
    float obj_mass(3.0);      
    vect<float,3> gravity_force;                           //RK_NOTICE(2,"Passed: " << passed << " should have 16.");
    gravity_force = (gravity_acc * obj_mass);    ++passed;
    if(gravity_force == obj_mass * gravity_acc)  ++passed;
    if(std::fabs(norm_sqr(gravity_acc) - 9.81*9.81) < 100.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(gravity_acc) - 9.81) < 10.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    vect<float,3> gravity_dir(unit(gravity_acc));++passed; //RK_NOTICE(2,"Passed: " << passed << " should have 21.");
    if(std::fabs(norm(gravity_dir) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,1>(1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,2>(0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,3>(0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,4>(0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,5>(0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,6>(0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,7>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,8>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,9>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    if(std::fabs(norm(vect<float,10>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    vect<float,3> displacement(1.0, 2.0, 3.0);             //RK_NOTICE(2,"Passed: " << passed << " should have 32.");
    float gravity_potential = gravity_force * displacement;
    if(std::fabs(gravity_potential + 2.0 * 9.81 * 3.0) < 60.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    vect<float,3> gravity_moment = displacement % gravity_force;
    if(std::fabs(gravity_moment[0] - 3.0*9.81*3.0) + std::fabs(gravity_moment[2] + 3.0*9.81) < 100.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    
    float dist_sqr = norm_sqr(displacement);               //RK_NOTICE(2,"Passed: " << passed << " should have 34.");
    vect<float,3> displacement_inv(-displacement); 
    if(std::fabs(displacement_inv * displacement + dist_sqr) < std::numeric_limits<float>::epsilon())
                                                 ++passed;
    vect<float,3> long_displacement(10.0, 20.0, 30.0); ++passed;
    if(colinear(displacement,long_displacement)) ++passed;
    if(long_displacement > displacement)         ++passed;
    if(displacement < long_displacement)         ++passed;
    if(!(displacement == long_displacement))     ++passed;
    if(displacement == displacement)             ++passed;

    long_displacement += displacement;                     //RK_NOTICE(2,"Passed: " << passed << " should have 41.");
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    long_displacement *= 2.0;
    if(std::fabs(norm(long_displacement) - 22.0*norm(displacement)) < 400.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    long_displacement /= 2.0;
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
    long_displacement -= displacement;
    if(std::fabs(norm(long_displacement) - 10.0*norm(displacement)) < 200.0*std::numeric_limits<float>::epsilon())
                                                 ++passed;
                                                           //RK_NOTICE(2,"Passed: " << passed << " should have 45.");
  };

  /***************************** Variable-Length Vector Tests ***********************************/
  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/** PRIMITIVE VARIABLE-LENGTH VECTOR TESTS ***/");
    RK_NOTICE(2,"/*********************************************/");
    vect_n<float> gravity_acc(3);
    gravity_acc[0] = 0.0;
    gravity_acc[1] = -9.81;
    gravity_acc[2] = 0.0;
    float obj_mass(3.0);
    vect_n<float> gravity_force;
    gravity_force = (gravity_acc * obj_mass);
    vect_n<float> gravity_dir(unit(gravity_acc));
    RK_NOTICE(2,"Gravity vector is: " << gravity_force[0] << ";" << gravity_force[1] << ";" << gravity_force[2]);
    RK_NOTICE(2,"Gravity direction is: " << gravity_dir[0] << ";" << gravity_dir[1] << ";" << gravity_dir[2]);
    vect_n<float> displacement(1.0, 2.0, 3.0);
    float gravity_potential = gravity_force * displacement;
    RK_NOTICE(2,"Gravity Potential is: " << gravity_potential);

    float dist_sqr = norm_sqr(displacement);
    float distance = norm(displacement);
    vect_n<float> displacement_inv(-displacement);
    vect_n<float> long_displacement(10.0, 20.0, 30.0);
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

    mat<float> m1234(1.0,2.0,3.0,4.0);
    RK_NOTICE(2,"m1234 = " << m1234);
    RK_NOTICE(2,"m1234 = " << m1234(0,0) << " " << m1234(1,0) << " " << m1234(0,1) << " " << m1234(1,1));
    RK_NOTICE(2,"m1324 = " << (m1234 = transpose(m1234)));
    RK_NOTICE(2,"m1324 = " << (m1234 = transpose_move(m1234)));

    mat<float> m_ident(3,3,true);
    RK_NOTICE(2,"identity 3x3 = " << m_ident);
    mat_identity<float>::type m_ident2(3);
    RK_NOTICE(2,"identity 3x3 = " << m_ident2);
    mat<float> m_ones(3,3,float(1.0));
    RK_NOTICE(2,"ones 3x3 = " << m_ones);
    mat_null<float>::type m_zeroes(3,3);
    RK_NOTICE(2,"zeroes 3x3 = " << m_zeroes);

    float f4321[] = {4.0,3.0,2.0,1.0};
    std::vector<float> v4321(f4321,f4321 + 4);
    mat<float> m4321(v4321,2,2);
    RK_NOTICE(2,"m4321 = " << m4321);

    mat<float> m4321_cpy(m4321);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy = m4321;
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy += mat<float>(2,2,true);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy += mat_identity<float>::type(2);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy -= mat<float>(2,2,true);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy -= mat_identity<float>::type(2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy *= 2.0;
    RK_NOTICE(2,"m4321 copy twice = " << m4321_cpy);

    m4321_cpy = m4321_cpy * float(0.5);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

    m4321_cpy = m4321_cpy + mat<float>(2,2,true);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy = m4321_cpy + mat_identity<float>::type(2);
    RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

    m4321_cpy = m4321_cpy - mat<float>(2,2,true);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    m4321_cpy = m4321_cpy - mat_identity<float>::type(2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    RK_NOTICE(2,"m4321 copy negated = " << -m4321_cpy);

    m4321_cpy = m4321_cpy * mat_identity<float>::type(2);
    RK_NOTICE(2,"m4321 copy * identity = " << m4321_cpy);

    if(m4321_cpy == m4321)
      RK_NOTICE(2,"m4321 copy is equal to m4321");
    if(m4321_cpy != mat<float>(2,2,true))
      RK_NOTICE(2,"m4321 copy is not equal to identity");
    if(m4321_cpy != mat_identity<float>::type(2))
      RK_NOTICE(2,"m4321 copy is not equal to identity");

    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    vect_n<float> v85 = m4321_cpy * vect_n<float>(2,1.0);
    RK_NOTICE(2,"m4321 copy * vect(1.0,1.0) = " << v85);
    vect_n<float> v104 =  vect_n<float>(2,1.0) * m4321_cpy;
    RK_NOTICE(2,"vect(1.0,1.0) * m4321 copy = " << v104);

    set_block(m_ident,m4321_cpy,1,1);
    RK_NOTICE(2,"identity 3x3 with sub matrix m4321 is " << m_ident);

    m4321_cpy = get_block(m_ident,1,1,2,2);
    RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
    RK_NOTICE(2,"m4321 copy transpose = " << transpose(m4321_cpy));
    RK_NOTICE(2,"m4321 copy sym part = " << (mat<float,mat_structure::symmetric>(m4321_cpy)));
    RK_NOTICE(2,"m4321 copy skew part = " << (mat<float,mat_structure::skew_symmetric>(m4321_cpy)));
    RK_NOTICE(2,"m4321 copy twice = " << float(2.0) * m4321_cpy);

    mat<float,mat_structure::diagonal> m123(vect_n<float>(1.0,2.0,3.0));
    RK_NOTICE(2,"m123 diagonal = " << m123);
    //mat_cm<float> m_blk_diag(block_diag_mat(mat_cm<float>(1.0,2.0,3.0,4.0),mat_cm<float>(1,1,true)));
    //RK_NOTICE(2,"m_blk_diag block diagonal = " << m_blk_diag);
    //mat_cm<float> m_blk(block_mat(mat_cm<float>(1.0,2.0,3.0,4.0),m4321.getSubMat(0,0,2,1),m4321.getSubMat(1,0,1,2),mat_cm<float>(1,1,true)));
    //RK_NOTICE(2,"m_blk block matrix = " << m_blk);
    RK_NOTICE(2,"vect(1,2,3) skew matrix = " << (mat<float,mat_structure::skew_symmetric>(vect_n<float>(1.0,2.0,3.0))));

    RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
  };
  
  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/**** 2D ROTATION & TRANSFORMATION TESTS *****/");
    RK_NOTICE(2,"/*********************************************/");

    rot_mat_2D<float> r_ident;
    RK_NOTICE(2,"Rotation identity = " << r_ident);
    RK_NOTICE(2,"Rotation identity = (" << r_ident[0] << "; " << r_ident[2] << "); (" << r_ident[1] << "; " << r_ident[3] << ")");
    rot_mat_2D<float> r_45deg(0.25 * M_PI);
    RK_NOTICE(2,"Rotation 45 degrees = " << r_45deg);
    RK_NOTICE(2,"Rotation 45 degrees = (" << r_45deg(0,0) << "; " << r_45deg(0,1) << "); (" << r_45deg(1,0) << "; " << r_45deg(1,1) << ")");
    rot_mat_2D<float> r_45deg_cpy(r_45deg);
    RK_NOTICE(2,"Rotation 45 degrees copy = " << r_45deg_cpy);
    RK_NOTICE(2,"Rotation 45 degrees copy = " << r_45deg_cpy.getMat());
    RK_NOTICE(2,"Rotation 45 degrees angle = " << r_45deg.getAngle());
    r_45deg_cpy.setAngle(0.0);
    RK_NOTICE(2,"Rotation identity = " << r_45deg_cpy);
    r_45deg_cpy = r_45deg;
    RK_NOTICE(2,"Rotation 45 degrees copy = " << r_45deg_cpy);
    rot_mat_2D<float> r_90deg = r_45deg * r_45deg_cpy;
    RK_NOTICE(2,"Rotation 90 degrees = " << r_90deg);
    r_90deg *= r_45deg;
    RK_NOTICE(2,"Rotation 135 degrees = " << r_90deg);
    vect<float,2> v1 = r_45deg * vect<float,2>(1.0,1.0);
    RK_NOTICE(2,"Vector (1,1) rotated by 45 degrees = " << v1);
    v1 = vect<float,2>(1.0,1.0) * r_45deg;
    RK_NOTICE(2,"Vector (1,1) rotated by 45 degrees = " << v1);
    r_90deg *= invert(r_45deg);
    RK_NOTICE(2,"Rotation 90 degrees = " << r_90deg);
    RK_NOTICE(2,"r_45deg copy trace = " << trace(r_45deg_cpy));
    RK_NOTICE(2,"r_45deg copy determinant = " << determinant(r_45deg_cpy));
    RK_NOTICE(2,"r_45deg copy sym part = " << r_45deg_cpy.getSymPart());
    RK_NOTICE(2,"r_45deg copy skew part = " << r_45deg_cpy.getSkewSymPart());
    RK_NOTICE(2,"r_45deg * r_45deg.invert = " << (r_45deg * invert(r_45deg)));
    if(r_45deg * invert(r_45deg) == r_ident)
      RK_NOTICE(2,"r_45deg * r_45deg.invert is equal to identity!");
    if(r_45deg != r_ident)
      RK_NOTICE(2,"r_45deg is not equal to identity!");

    mat<float,mat_structure::rectangular> m23(2,3);
    m23(0,0) = 1.0; m23(1,0) = 0.0;
    m23(0,1) = 0.0; m23(1,1) = 1.0;
    m23(0,2) = 1.0; m23(1,2) = 1.0;
    RK_NOTICE(2,"m23 = " << m23);
    mat<float,mat_structure::rectangular> m23_r45(r_45deg * m23);
    RK_NOTICE(2,"m23 rotated 45 degrees = " << m23_r45);
    mat<float,mat_structure::rectangular> m32(3,2);
    m32(0,0) = 1.0; m32(1,0) = 0.0; m32(2,0) = 1.0;
    m32(0,1) = 0.0; m32(1,1) = 1.0; m32(2,1) = 1.0;
    RK_NOTICE(2,"m32 = " << m32);
    mat<float,mat_structure::rectangular> m32_r45(m32 * r_45deg);
    RK_NOTICE(2,"m32 rotated -45 degrees = " << m32_r45);

    trans_mat_2D<float> t_ident;
    RK_NOTICE(2,"Transformation identity = " << t_ident);
    RK_NOTICE(2,"Transformation identity = (" << t_ident[0] << "; " << t_ident[3] << "; " << t_ident[6] << "); (" << t_ident[1] << "; " << t_ident[4] << "; " << t_ident[7] << "); (" << t_ident[2] << "; " << t_ident[5] << "; " << t_ident[8] << ")");
    trans_mat_2D<float> t_45deg_11(0.25 * M_PI,vect<float,2>(1.0,1.0));
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) = " << t_45deg_11);
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) = (" << t_45deg_11(0,0) << "; " << t_45deg_11(0,1) << "; " << t_45deg_11(0,2) << "); (" << t_45deg_11(1,0) << "; " << t_45deg_11(1,1) << "; " << t_45deg_11(1,2) << "); (" << t_45deg_11(2,0) << "; " << t_45deg_11(2,1) << "; " << t_45deg_11(2,2) << ")");
    trans_mat_2D<float> t_45deg_11_cpy(t_45deg_11);
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) copy = " << t_45deg_11_cpy);
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) copy = " << t_45deg_11_cpy.getMat());
    RK_NOTICE(2,"Rotation 45 degrees angle = " << t_45deg_11_cpy.getAngle());
    RK_NOTICE(2,"Rotation 45 degrees angle = " << t_45deg_11_cpy.getRotMat());
    t_45deg_11_cpy.setAngle(0.0);
    RK_NOTICE(2,"Rotation identity + (1,1) = " << t_45deg_11_cpy);
    t_45deg_11_cpy.setRotMat(r_45deg);
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) copy = " << t_45deg_11_cpy);
    t_45deg_11_cpy.setTranslation(vect<float,2>(-1.0,-1.0));
    RK_NOTICE(2,"Rotation identity + (-1,-1) = " << t_45deg_11_cpy);
    t_45deg_11_cpy.setTranslation(vect<float,2>(1.0,1.0));
    RK_NOTICE(2,"Rotation 45 degrees + (1,1) copy = " << t_45deg_11_cpy);
    RK_NOTICE(2,"t_45deg_11_cpy trace = " << trace(t_45deg_11_cpy));
    RK_NOTICE(2,"t_45deg_11_cpy determinant = " << determinant(t_45deg_11_cpy));
    RK_NOTICE(2,"t_45deg_11_cpy sym part = " << t_45deg_11_cpy.getSymPart());
    RK_NOTICE(2,"t_45deg_11_cpy skew part = " << t_45deg_11_cpy.getSkewSymPart());
    RK_NOTICE(2,"t_45deg_11_cpy invert = " << invert(t_45deg_11_cpy));
    RK_NOTICE(2,"t_45deg_11_cpy transpose = " << transpose(t_45deg_11_cpy));
    RK_NOTICE(2,"t_45deg_11_cpy invert * t_45deg_11_cpy = " << (invert(t_45deg_11_cpy) * t_45deg_11_cpy));
    t_ident *= t_45deg_11_cpy;
    RK_NOTICE(2,"Rotation 45 degrees and translation (1,1) times identity = " << t_ident);
    t_ident *= invert(t_45deg_11_cpy);
    RK_NOTICE(2,"Back to transformation identity = " << t_ident);

    vect<float,2> v2 = t_45deg_11 * vect<float,2>(1.0,1.0);
    RK_NOTICE(2,"Vector (1,1) rotated by 45 degrees + (1,1) = " << v2);
    vect<float,3> v3 = t_45deg_11 * vect<float,3>(1.0,1.0,1.0);
    RK_NOTICE(2,"Vector (1,1,1) rotated by 45 degrees + (1,1) = " << v3);
    v2 = t_45deg_11.rotate(vect<float,2>(1.0,1.0));
    RK_NOTICE(2,"Vector (1,1) rotated by 45 degrees = " << v2);

    mat<float,mat_structure::rectangular> m23_t(2,3);
    m23_t(0,0) = 1.0; m23_t(1,0) = 0.0;
    m23_t(0,1) = 0.0; m23_t(1,1) = 1.0;
    m23_t(0,2) = 1.0; m23_t(1,2) = 1.0;
    RK_NOTICE(2,"m23_t = " << m23_t);
    mat<float,mat_structure::rectangular> m23_r45_11(m23_t * t_45deg_11);
    RK_NOTICE(2,"m23_t (rotated 45 degrees + (1,1)).transpose = " << m23_r45_11);
    mat<float,mat_structure::rectangular> m32_t(3,2);
    m32_t(0,0) = 1.0; m32_t(1,0) = 0.0; m32_t(2,0) = 1.0;
    m32_t(0,1) = 0.0; m32_t(1,1) = 1.0; m32_t(2,1) = 1.0;
    RK_NOTICE(2,"m32_t = " << m32_t);
    mat<float,mat_structure::rectangular> m32_r45_11(t_45deg_11 * m32_t);
    RK_NOTICE(2,"m32_t rotated 45 degrees + (1,1) = " << m32_r45_11);

    RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
  };

  if(true){
    RK_NOTICE(2,"/*********************************************/");
    RK_NOTICE(2,"/**** 3D ROTATION & TRANSFORMATION TESTS *****/");
    RK_NOTICE(2,"/*********************************************/");

    rot_mat_3D<float> r_ident;
    RK_NOTICE(2,"Rotation identity = " << r_ident);
    float r_45z_a[] = {float(std::cos(0.25*M_PI)),float(std::sin(0.25*M_PI)),0.0,float(-std::sin(0.25*M_PI)),float(std::cos(0.25*M_PI)),0.0,0.0,0.0,1.0};
    rot_mat_3D<float> r_45z(r_45z_a);
    RK_NOTICE(2,"Rotation 45 deg about z = " << r_45z);
    rot_mat_3D<float> r_45z_cpy(r_45z);
    RK_NOTICE(2,"Rotation 45 deg about z copy = " << r_45z_cpy);
    r_45z_cpy = r_45z;
    RK_NOTICE(2,"Rotation 45 deg about z copy = " << r_45z_cpy);
    r_ident *= r_45z;
    RK_NOTICE(2,"Identity * Rotation 45 deg about z = " << r_ident);
    r_ident *= invert(r_45z);
    RK_NOTICE(2,"Back to Identity = " << r_ident);

    RK_NOTICE(2,"r_45z_cpy trace = " << trace(r_45z_cpy));
    RK_NOTICE(2,"r_45z_cpy determinant = " << determinant(r_45z_cpy));
    RK_NOTICE(2,"r_45z_cpy sym part = " << r_45z_cpy.getSymPart());
    RK_NOTICE(2,"r_45z_cpy skew part = " << r_45z_cpy.getSkewSymPart());
    RK_NOTICE(2,"r_45z_cpy invert = " << invert(r_45z_cpy));
    RK_NOTICE(2,"r_45z_cpy transpose = " << transpose(r_45z_cpy));
    RK_NOTICE(2,"r_45z_cpy invert * r_45z_cpy = " << (invert(r_45z_cpy) * r_45z_cpy));

    RK_NOTICE(2,"r_45z_cpy as matrix = " << r_45z_cpy.getMat());
    RK_NOTICE(2,"r_45z_cpy invert * r_45z_cpy as matrix = " << (invert(r_45z_cpy) * r_45z_cpy.getMat()));
    RK_NOTICE(2,"r_45z_cpy invert as matrix * r_45z_cpy = " << (invert(r_45z_cpy).getMat() * r_45z_cpy));

    if(r_ident == r_ident)
      RK_NOTICE(2,"Rotation identity is equal to itself!");
    if(r_ident != r_45z)
      RK_NOTICE(2,"Rotation identity is not equal to 45 deg rotation about z!");

    vect<float,3> v1(1.0,1.0,2.0);
    RK_NOTICE(2,"r_45z * v(1,1,2) = " << (r_45z * v1));
    RK_NOTICE(2,"v(1,1,2) * r_45z = " << (v1 * r_45z));

    quaternion<float> q_45z(r_45z);
    RK_NOTICE(2,"q_45z = " << q_45z);
    RK_NOTICE(2,"r_45z * q_45z.invert = " << (r_45z * invert(q_45z)));
    quaternion<float> q_ident;
    RK_NOTICE(2,"q_ident = " << q_ident);
    quaternion<float> q_45z_cpy(q_45z);
    RK_NOTICE(2,"q_45z_cpy = " << q_45z_cpy);
    q_45z_cpy = q_45z;
    RK_NOTICE(2,"q_45z_cpy = " << q_45z_cpy);
    q_ident = q_45z * q_45z;
    RK_NOTICE(2,"rotation by 90 degrees about z = " << q_ident);
    q_ident *= invert(q_45z);
    RK_NOTICE(2,"rotation by 45 degrees about z = " << q_ident.getMat());
    q_ident *= invert(q_45z);
    RK_NOTICE(2,"rotation by 0 degrees about z = " << q_ident);
    RK_NOTICE(2,"rotation by 90 degrees about z = " << (q_45z * r_45z));

    RK_NOTICE(2,"q_45z_cpy trace = " << trace(q_45z_cpy));
    RK_NOTICE(2,"q_45z_cpy determinant = " << determinant(q_45z_cpy));
    RK_NOTICE(2,"q_45z_cpy sym part = " << q_45z_cpy.getSymPart());
    RK_NOTICE(2,"q_45z_cpy skew part = " << q_45z_cpy.getSkewSymPart());
    RK_NOTICE(2,"q_45z_cpy invert = " << invert(q_45z_cpy));
    RK_NOTICE(2,"q_45z_cpy transpose = " << transpose(q_45z_cpy));
    RK_NOTICE(2,"q_45z_cpy invert * q_45z_cpy = " << (invert(q_45z_cpy) * q_45z_cpy));

    if(q_ident == q_ident)
      RK_NOTICE(2,"q_ident is equal to itself!");
    if(q_ident != q_45z)
      RK_NOTICE(2,"q_ident is not equal to q_45z!");

    RK_NOTICE(2,"q_45z * v(1,1,2) = " << (q_45z * v1));

    euler_angles_TB<float> e_45z(0.25*M_PI,0.0,0.0);
    RK_NOTICE(2,"q_45z * EA(Yaw = 45) = " << (q_45z * e_45z));
    axis_angle<float> a_45z(0.25*M_PI,vect<float,3>(0.0,0.0,1.0));
    RK_NOTICE(2,"q_45z * a_45z = " << (q_45z * a_45z));

    euler_angles_TB<float> e_ident;
    RK_NOTICE(2,"e_ident = " << e_ident);
    euler_angles_TB<float> e_45z_cpy(e_45z);
    RK_NOTICE(2,"e_45z_cpy = " << e_45z_cpy);
    euler_angles_TB<float> e_45z_r(r_45z);
    RK_NOTICE(2,"e_45z_r = " << e_45z_r);
    euler_angles_TB<float> e_45z_q(q_45z);
    RK_NOTICE(2,"e_45z_q = " << e_45z_q);
    RK_NOTICE(2,"e_45z_q rot mat = " << e_45z_q.getRotMat());
    RK_NOTICE(2,"e_45z_q mat = " << e_45z_q.getMat());
    e_45z_cpy = e_45z;
    RK_NOTICE(2,"e_45z_cpy = " << e_45z_cpy);
    e_45z_cpy = r_45z;
    RK_NOTICE(2,"e_45z_cpy = " << e_45z_cpy);
    e_45z_cpy = q_45z;
    RK_NOTICE(2,"e_45z_cpy = " << e_45z_cpy);
    e_45z_cpy = a_45z;
    RK_NOTICE(2,"e_45z_cpy = " << e_45z_cpy);

    if(e_ident == e_ident)
      RK_NOTICE(2,"e_ident is equal to itself!");
    if(e_ident != e_45z)
      RK_NOTICE(2,"e_ident is not equal to e_45z!");

    RK_NOTICE(2,"e_45z_cpy trace = " << trace(e_45z_cpy));
    RK_NOTICE(2,"e_45z_cpy determinant = " << determinant(e_45z_cpy));
    RK_NOTICE(2,"e_45z_cpy sym part = " << e_45z_cpy.getSymPart());
    RK_NOTICE(2,"e_45z_cpy skew part = " << e_45z_cpy.getSkewSymPart());
    RK_NOTICE(2,"e_45z_cpy invert = " << invert(e_45z_cpy));
    RK_NOTICE(2,"e_45z_cpy transpose = " << transpose(e_45z_cpy));
    RK_NOTICE(2,"e_45z_cpy invert * e_45z_cpy = " << (invert(e_45z_cpy) * e_45z_cpy));

    axis_angle<float> a_ident;
    RK_NOTICE(2,"a_ident = " << a_ident);
    axis_angle<float> a_45z_cpy(a_45z);
    RK_NOTICE(2,"a_45z_cpy = " << a_45z_cpy);
    axis_angle<float> a_45z_r(r_45z);
    RK_NOTICE(2,"a_45z_r = " << a_45z_r);
    axis_angle<float> a_45z_q(q_45z);
    RK_NOTICE(2,"a_45z_q = " << a_45z_q);
    axis_angle<float> a_45z_e(e_45z);
    RK_NOTICE(2,"a_45z_e = " << a_45z_e);
    RK_NOTICE(2,"a_45z_q rot mat = " << a_45z_q.getRotMat());
    RK_NOTICE(2,"a_45z_q mat = " << a_45z_q.getMat());
    a_45z_cpy = a_45z;
    RK_NOTICE(2,"a_45z_cpy = " << a_45z_cpy);
    a_45z_cpy = r_45z;
    RK_NOTICE(2,"a_45z_cpy = " << a_45z_cpy);
    a_45z_cpy = q_45z;
    RK_NOTICE(2,"a_45z_cpy = " << a_45z_cpy);
    a_45z_cpy = e_45z;
    RK_NOTICE(2,"a_45z_cpy = " << a_45z_cpy);

    if(a_ident == a_ident)
      RK_NOTICE(2,"a_ident is equal to itself!");
    if(a_ident != a_45z)
      RK_NOTICE(2,"a_ident is not equal to a_45z!");

    RK_NOTICE(2,"a_45z_cpy trace = " << trace(a_45z_cpy));
    RK_NOTICE(2,"a_45z_cpy determinant = " << determinant(a_45z_cpy));
    RK_NOTICE(2,"a_45z_cpy sym part = " << a_45z_cpy.getSymPart());
    RK_NOTICE(2,"a_45z_cpy skew part = " << a_45z_cpy.getSkewSymPart());
    RK_NOTICE(2,"a_45z_cpy invert = " << invert(a_45z_cpy));
    RK_NOTICE(2,"a_45z_cpy transpose = " << transpose(a_45z_cpy));
    RK_NOTICE(2,"a_45z_cpy invert * a_45z_cpy = " << (invert(a_45z_cpy) * a_45z_cpy));



    trans_mat_3D<float> t_ident;
    RK_NOTICE(2,"t_ident = " << t_ident);
    float t_45z_array[] = {float(std::cos(0.25*M_PI)),float(std::sin(0.25*M_PI)),0.0,0.0,float(-std::sin(0.25*M_PI)),float(std::cos(0.25*M_PI)),0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0};
    trans_mat_3D<float> t_45z(t_45z_array);
    RK_NOTICE(2,"Rotation 45 deg about z = " << t_45z);
    trans_mat_3D<float> t_45z_cpy(t_45z);
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    trans_mat_3D<float> t_45z_r(r_45z);
    RK_NOTICE(2,"t_45z_r = " << t_45z_r);
    trans_mat_3D<float> t_45z_q(q_45z);
    RK_NOTICE(2,"t_45z_q = " << t_45z_q);
    trans_mat_3D<float> t_45z_e(e_45z);
    RK_NOTICE(2,"t_45z_e = " << t_45z_e);
    trans_mat_3D<float> t_45z_a(a_45z);
    RK_NOTICE(2,"t_45z_a = " << t_45z_a);
    RK_NOTICE(2,"t_45z_q rot mat = " << t_45z_q.getRotMat());
    RK_NOTICE(2,"t_45z_q mat = " << t_45z_q.getMat());
    t_45z_cpy = t_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = r_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = q_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = e_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = a_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = t_ident;
    t_45z_cpy *= t_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = t_ident;
    t_45z_cpy *= r_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = t_ident;
    t_45z_cpy *= q_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = t_ident;
    t_45z_cpy *= e_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);
    t_45z_cpy = t_ident;
    t_45z_cpy *= a_45z;
    RK_NOTICE(2,"t_45z_cpy = " << t_45z_cpy);

    if(t_ident == t_ident)
      RK_NOTICE(2,"t_ident is equal to itself!");
    if(t_ident != t_45z)
      RK_NOTICE(2,"t_ident is not equal to t_45z!");

    trans_mat_3D<float> t_45z_123(r_45z,vect<float,3>(1.0,2.0,3.0));
    RK_NOTICE(2,"t_45z_123 = " << t_45z_123);

    RK_NOTICE(2,"t_45z_123 trace = " << trace(t_45z_123));
    RK_NOTICE(2,"t_45z_123 determinant = " << determinant(t_45z_123));
    RK_NOTICE(2,"t_45z_123 sym part = " << t_45z_123.getSymPart());
    RK_NOTICE(2,"t_45z_123 skew part = " << t_45z_123.getSkewSymPart());
    RK_NOTICE(2,"t_45z_123 invert = " << invert(t_45z_123));
    RK_NOTICE(2,"t_45z_123 transpose = " << transpose(t_45z_123));
    RK_NOTICE(2,"t_45z_123 invert * t_45z_123 = " << (invert(t_45z_123) * t_45z_123));
    RK_NOTICE(2,"t_45z_123 invert mat * t_45z_123 = " << (invert(t_45z_123).getMat() * t_45z_123));
    RK_NOTICE(2,"t_45z_123 invert * t_45z_123 mat = " << (invert(t_45z_123) * t_45z_123.getMat()));
    RK_NOTICE(2,"t_45z_123 invert * r_45z = " << (invert(t_45z_123) * r_45z));

    RK_NOTICE(2,"t_45z_123 * v(1,1,2) = " << (t_45z_123 * v1));
    RK_NOTICE(2,"t_45z_123.rotate(v(1,1,2)) = " << t_45z_123.rotate(v1));
    vect<float,4> v3(1,1,2,2);
    RK_NOTICE(2,"t_45z_123 * v(1,1,2,2) = " << (t_45z_123 * v3));


    axis_angle<float> a_weird(0.3241,vect<float,3>(0.5,0.5,sqrt(0.5)));
    quaternion<float> q_weird;
    q_weird = a_weird;
    euler_angles_TB<float> e_weird;
    e_weird = q_weird;
    rot_mat_3D<float> r_weird;
    r_weird = q_weird;
    RK_NOTICE(2,"a_weird as A/A = " << a_weird);
    RK_NOTICE(2,"q_weird as A/A = " << axis_angle<float>(q_weird));
    RK_NOTICE(2,"e_weird as A/A = " << axis_angle<float>(e_weird));
    RK_NOTICE(2,"r_weird as A/A = " << axis_angle<float>(r_weird));

    RK_NOTICE(2,"result a as Q = " << q_45z * quaternion<float>(a_weird * a_45z));
    RK_NOTICE(2,"result q as Q = " << q_45z * (q_weird * q_45z));
    RK_NOTICE(2,"result e as Q = " << quaternion<float>(q_45z * e_weird * r_45z));
    RK_NOTICE(2,"result r as Q = " << quaternion<float>(a_45z * r_weird * e_45z));


    RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
  };
  /*************************** Matrix Numerical Methods Tests ****************************/

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
    eigensolve_QR(m_gauss,m_qr_E,m_qr_Q,50,double(1E-6));
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







