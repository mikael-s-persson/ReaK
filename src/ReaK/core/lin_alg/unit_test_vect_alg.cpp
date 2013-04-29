
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


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE vect_alg
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>



typedef boost::mpl::list< 
  ReaK::vect<float,3>,
  ReaK::vect<double,3> > vect3d_tests_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( vect3d_tests, Vector, vect3d_tests_types )
{
  typedef typename ReaK::vect_traits<Vector>::value_type ValueType;
  using std::fabs;
  
  Vector gravity_acc(0.0, -9.81, 0.0);
  BOOST_CHECK_EQUAL( gravity_acc.size(), 3 );
  BOOST_CHECK( gravity_acc.max_size() );
  BOOST_CHECK( gravity_acc.capacity() );
  BOOST_CHECK( !gravity_acc.empty() );
  
  typename Vector::iterator it = gravity_acc.begin();      //RK_NOTICE(2,"Passed: " << passed << " should have 8.");
  BOOST_CHECK_SMALL( *it, std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_EQUAL( gravity_acc.end() - it, 3 );
  ++it;
  BOOST_CHECK_EQUAL( gravity_acc.end() - it, 2 );
  
  const Vector& gravity_acc_ref = gravity_acc;
  typename Vector::const_iterator cit = gravity_acc_ref.begin();
  BOOST_CHECK_SMALL( *cit, std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_EQUAL( gravity_acc.end() - cit, 3 );
  ++cit;
  BOOST_CHECK_EQUAL( gravity_acc.end() - cit, 2 );
  
  ValueType ones[] = {1.0,1.0,1.0};
  Vector ones_v(ones);
  BOOST_CHECK_CLOSE( ones[1], ones_v[1], std::numeric_limits< ValueType >::epsilon() );
  
  ValueType obj_mass(3.0);
  Vector gravity_force; 
  gravity_force = (gravity_acc * obj_mass);
  BOOST_CHECK_CLOSE( gravity_force[1], ValueType(-3.0*9.81), std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK( (gravity_force == obj_mass * gravity_acc) );
  BOOST_CHECK( (gravity_force == obj_mass * gravity_acc) );
  BOOST_CHECK_CLOSE( norm_2_sqr(gravity_acc), ValueType(9.81*9.81), 100.0*std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_CLOSE( norm_2(gravity_acc), ValueType(9.81), 10.0*std::numeric_limits< ValueType >::epsilon() );
  
  Vector gravity_dir(unit(gravity_acc));
  BOOST_CHECK_CLOSE( norm_2(gravity_dir), ValueType(1.0), std::numeric_limits< ValueType >::epsilon() );
  
  
  Vector displacement(1.0, 2.0, 3.0);
  ValueType gravity_potential = gravity_force * displacement;
  BOOST_CHECK_CLOSE( gravity_potential, ValueType(-2.0 * 9.81 * 3.0), 100.0*std::numeric_limits< ValueType >::epsilon() );
  
  Vector gravity_moment = displacement % gravity_force;
  BOOST_CHECK_SMALL( (fabs(gravity_moment[0] - 3.0*9.81*3.0) + fabs(gravity_moment[2] + 3.0*9.81)), 100.0*std::numeric_limits< ValueType >::epsilon() );
  
  ValueType dist_sqr = norm_2_sqr(displacement);
  Vector displacement_inv(-displacement); 
  BOOST_CHECK_CLOSE( displacement_inv * displacement, -dist_sqr, std::numeric_limits< ValueType >::epsilon() );
  
  Vector long_displacement(10.0, 20.0, 30.0);
  BOOST_CHECK( colinear(displacement,long_displacement) );
  BOOST_CHECK( (long_displacement > displacement) );
  BOOST_CHECK( (displacement < long_displacement) );
  BOOST_CHECK( (displacement != long_displacement) );
  BOOST_CHECK( (displacement == displacement) );
  
  long_displacement += displacement;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 11.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement *= 2.0;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 22.0*norm_2(displacement), 400.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement /= 2.0;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 11.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement -= displacement;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 10.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  
};


typedef boost::mpl::list< 
  ReaK::vect_n<float>,
  ReaK::vect_n<double> > vectn_tests_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( vectn_tests, Vector, vectn_tests_types )
{
  typedef typename ReaK::vect_traits<Vector>::value_type ValueType;
  using std::fabs;
  
  Vector gravity_acc(0.0, -9.81, 0.0);
  BOOST_CHECK_EQUAL( gravity_acc.size(), 3 );
  BOOST_CHECK( gravity_acc.max_size() );
  BOOST_CHECK( gravity_acc.capacity() );
  BOOST_CHECK( !gravity_acc.empty() );
  
  typename Vector::iterator it = gravity_acc.begin();      //RK_NOTICE(2,"Passed: " << passed << " should have 8.");
  BOOST_CHECK_SMALL( *it, std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_EQUAL( gravity_acc.end() - it, 3 );
  ++it;
  BOOST_CHECK_EQUAL( gravity_acc.end() - it, 2 );
  
  const Vector& gravity_acc_ref = gravity_acc;
  typename Vector::const_iterator cit = gravity_acc_ref.begin();
  BOOST_CHECK_SMALL( *cit, std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_EQUAL( gravity_acc.end() - cit, 3 );
  ++cit;
  BOOST_CHECK_EQUAL( gravity_acc.end() - cit, 2 );
  
  ValueType obj_mass(3.0);
  Vector gravity_force; 
  gravity_force = (gravity_acc * obj_mass);
  BOOST_CHECK_CLOSE( gravity_force[1], ValueType(-3.0*9.81), std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK( (gravity_force == obj_mass * gravity_acc) );
  BOOST_CHECK( (gravity_force == obj_mass * gravity_acc) );
  BOOST_CHECK_CLOSE( norm_2_sqr(gravity_acc), ValueType(9.81*9.81), 100.0*std::numeric_limits< ValueType >::epsilon() );
  BOOST_CHECK_CLOSE( norm_2(gravity_acc), ValueType(9.81), 10.0*std::numeric_limits< ValueType >::epsilon() );
  
  Vector gravity_dir(unit(gravity_acc));
  BOOST_CHECK_CLOSE( norm_2(gravity_dir), ValueType(1.0), std::numeric_limits< ValueType >::epsilon() );
  
  Vector displacement(1.0, 2.0, 3.0);
  ValueType gravity_potential = gravity_force * displacement;
  BOOST_CHECK_CLOSE( gravity_potential, ValueType(-2.0 * 9.81 * 3.0), 100.0*std::numeric_limits< ValueType >::epsilon() );
  
  ValueType dist_sqr = norm_2_sqr(displacement);
  Vector displacement_inv(-displacement); 
  BOOST_CHECK_CLOSE( displacement_inv * displacement, -dist_sqr, std::numeric_limits< ValueType >::epsilon() );
  
  Vector long_displacement(10.0, 20.0, 30.0);
  BOOST_CHECK( colinear(displacement,long_displacement) );
  BOOST_CHECK( (long_displacement > displacement) );
  BOOST_CHECK( (displacement < long_displacement) );
  BOOST_CHECK( (displacement != long_displacement) );
  BOOST_CHECK( (displacement == displacement) );
  
  long_displacement += displacement;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 11.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement *= 2.0;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 22.0*norm_2(displacement), 400.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement /= 2.0;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 11.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  long_displacement -= displacement;
  BOOST_CHECK_CLOSE( norm_2(long_displacement), 10.0*norm_2(displacement), 200.0*std::numeric_limits< ValueType >::epsilon() );
  
};


BOOST_AUTO_TEST_CASE( vect_construction_tests )
{
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,1>(1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,2>(0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,3>(0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,4>(0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,5>(0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,7>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,8>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,9>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect<double,10>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  BOOST_CHECK_CLOSE( norm_2(ReaK::vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0)), 1.0, std::numeric_limits<double>::epsilon() ); 
  
};







