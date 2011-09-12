
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

#include "vect_alg.hpp"

int main() {

  using namespace ReaK;

  unsigned int passed = 0;

  try {

    {
  /*************************** Fixed-Length Vector Tests ****************************/
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
    if(std::fabs(norm_sqr(gravity_acc) - 9.81*9.81) < 100.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm_sqr(vect<double,3>) test did not pass!");
    if(std::fabs(norm(gravity_acc) - 9.81) < 10.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm(vect<double,3>) test did not pass!");
    vect<double,3> gravity_dir(unit(gravity_acc));
    if(std::fabs(norm(gravity_dir) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("unit(vect<double,3>) test did not pass!");
    if(std::fabs(norm(vect<double,1>(1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,1>(d) test did not pass!");
    if(std::fabs(norm(vect<double,2>(0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,2>(d,d) test did not pass!");
    if(std::fabs(norm(vect<double,3>(0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>(d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,4>(0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,4>(d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,5>(0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,5>(d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,6>(0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,6>(d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,7>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,7>(d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,8>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,8>(d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,9>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,9>(d,d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect<double,10>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
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
    
    double dist_sqr = norm_sqr(displacement);               //RK_NOTICE(2,"Passed: " << passed << " should have 34.");
    vect<double,3> displacement_inv(-displacement); 
    if(std::fabs(displacement_inv * displacement + dist_sqr) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::operator-() test did not pass!");
    vect<double,3> long_displacement(10.0, 20.0, 30.0);
    if(colinear(displacement,long_displacement))   
      ++passed;
    else
      RK_ERROR("colinear(vect<double,3>,vect<double,3>) test did not pass!");
    if(long_displacement > displacement)
      ++passed;
    else
      RK_ERROR("operator>(vect<double,3>,vect<double,3>) test did not pass!");
    if(displacement < long_displacement)
      ++passed;
    else
      RK_ERROR("operator<(vect<double,3>,vect<double,3>) test did not pass!");
    if(displacement != long_displacement)  
      ++passed;
    else
      RK_ERROR("operator!=(vect<double,3>,vect<double,3>) test did not pass!");
    if(displacement == displacement) 
      ++passed;
    else
      RK_ERROR("operator==(vect<double,3>,vect<double,3>) test did not pass!");
    
    long_displacement += displacement;
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::operator+=(vect<double,3>) test did not pass!");
    long_displacement *= 2.0;
    if(std::fabs(norm(long_displacement) - 22.0*norm(displacement)) < 400.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::operator*=(double) test did not pass!");
    long_displacement /= 2.0;
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::operator/=(double) test did not pass!");
    long_displacement -= displacement;
    if(std::fabs(norm(long_displacement) - 10.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect<double,3>::operator-=(vect<double,3>) test did not pass!");
    };

    {
  /*************************** Fixed-Length Vector Tests ****************************/
    vect_n<double> gravity_acc(3);
    if( 3 == gravity_acc.size() )
      ++passed;
    else
      RK_ERROR("vect_n<double>::size() test did not pass!");
    gravity_acc.resize(4);
    if( 4 == gravity_acc.size() )
      ++passed;
    else
      RK_ERROR("vect_n<double>::resize() test did not pass!");
    gravity_acc.resize(3);
    gravity_acc[0] = 0.0;
    gravity_acc[1] = -9.81;
    gravity_acc[2] = 0.0;
    if( !gravity_acc.empty() )
      ++passed;
    else
      RK_ERROR("vect_n<double>::empty() test did not pass!");
    vect_n<double>::iterator it = gravity_acc.begin();      //RK_NOTICE(2,"Passed: " << passed << " should have 8.");
    if(std::fabs(*it) < std::numeric_limits< double >::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::begin() test did not pass!");
    if(gravity_acc.end() - it == 3)
      ++passed;
    else
      RK_ERROR("vect_n<double>::end() test did not pass!");
    if(gravity_acc.end() - ++it == 2)
      ++passed;
    else
      RK_ERROR("vect_n<double>::iterator::operator++ test did not pass!");
    const vect_n<double>& gravity_acc_ref = gravity_acc;
    vect_n<double>::const_iterator cit = gravity_acc_ref.begin();
    if(std::fabs(*cit) < std::numeric_limits< double >::epsilon()) 
      ++passed;
    else
      RK_ERROR("vect_n<double>::cbegin() test did not pass!");
    if(gravity_acc_ref.end() - cit == 3) 
      ++passed;
    else
      RK_ERROR("vect_n<double>::cend() test did not pass!");
    if(gravity_acc_ref.end() - ++cit == 2) 
      ++passed;
    else
      RK_ERROR("vect_n<double>::const_iterator::operator++ test did not pass!");
    double ones[] = {1.0,1.0,1.0};
    vect_n<double> ones_v(ones, ones + 3);
    if(ones[1] == ones_v[1]) 
      ++passed;
    else
      RK_ERROR("vect_n<double>::vect(double*,double*) test did not pass!");
    double obj_mass(3.0);
    vect_n<double> gravity_force; 
    gravity_force = (gravity_acc * obj_mass);
    if(std::fabs(gravity_force[1] + 3.0*9.81) < std::numeric_limits<double>::epsilon()) 
      ++passed;
    else
      RK_ERROR("operator*(vect_n<double>,double) test did not pass!");
    if(gravity_force == obj_mass * gravity_acc) 
      ++passed;
    else
      RK_ERROR("operator*(double,vect_n<double>) test did not pass!");
    if(std::fabs(norm_sqr(gravity_acc) - 9.81*9.81) < 100.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm_sqr(vect_n<double>) test did not pass!");
    if(std::fabs(norm(gravity_acc) - 9.81) < 10.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("norm(vect_n<double>) test did not pass!");
    vect_n<double> gravity_dir(unit(gravity_acc));
    if(std::fabs(norm(gravity_dir) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("unit(vect_n<double>) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d,d,d,d,d) test did not pass!");
    if(std::fabs(norm(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)) - 1.0) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>(d,d,d,d,d,d,d,d,d,d) test did not pass!");
    vect_n<double> displacement(1.0, 2.0, 3.0);             
    double gravity_potential = gravity_force * displacement;
    if(std::fabs(gravity_potential + 2.0 * 9.81 * 3.0) < 60.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("operator*(vect_n<double>,vect_n<double>) test did not pass!");
    
    double dist_sqr = norm_sqr(displacement);               
    vect_n<double> displacement_inv(-displacement); 
    if(std::fabs(displacement_inv * displacement + dist_sqr) < std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::operator-() test did not pass!");
    vect_n<double> long_displacement(10.0, 20.0, 30.0);
    if(colinear(displacement,long_displacement))   
      ++passed;
    else
      RK_ERROR("colinear(vect_n<double>,vect_n<double>) test did not pass!");
    if(long_displacement > displacement)
      ++passed;
    else
      RK_ERROR("operator>(vect_n<double>,vect_n<double>) test did not pass!");
    if(displacement < long_displacement)
      ++passed;
    else
      RK_ERROR("operator<(vect_n<double>,vect_n<double>) test did not pass!");
    if(displacement != long_displacement)  
      ++passed;
    else
      RK_ERROR("operator!=(vect_n<double>,vect_n<double>) test did not pass!");
    if(displacement == displacement) 
      ++passed;
    else
      RK_ERROR("operator==(vect_n<double>,vect_n<double>) test did not pass!");
    
    long_displacement += displacement;
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::operator+=(vect_n<double>) test did not pass!");
    long_displacement *= 2.0;
    if(std::fabs(norm(long_displacement) - 22.0*norm(displacement)) < 400.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::operator*=(double) test did not pass!");
    long_displacement /= 2.0;
    if(std::fabs(norm(long_displacement) - 11.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::operator/=(double) test did not pass!");
    long_displacement -= displacement;
    if(std::fabs(norm(long_displacement) - 10.0*norm(displacement)) < 200.0*std::numeric_limits<double>::epsilon())
      ++passed;
    else
      RK_ERROR("vect_n<double>::operator-=(vect_n<double>) test did not pass!");
    };
    
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
    return 1;
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
    return 1;
  };
  
  RK_NOTICE(2,"There were " << passed << " / 72 successful tests passed on the vect_alg library.");
  
  if( passed >= 72 )
    return 0;
  else
    return 1;
};







