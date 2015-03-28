
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


#include <iostream>
#include <iomanip>

#include <set>
#include <map>

#include <ReaK/planning/graph_alg/avl_tree.hpp>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE assoc_containers
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


typedef boost::mpl::list< std::map< int, int >, ReaK::graph::avlbfl_map< int, int >
                          //   , ReaK::graph::avlvebl_map<int, int>
                          > intint_maptest_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( intint_map_test, Map, intint_maptest_types ) {
  typedef typename Map::key_type KeyType;
  typedef typename Map::mapped_type MappedType;
  typedef std::pair< KeyType, MappedType > ValueType;

  Map m;
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  BOOST_CHECK_THROW( m.at( KeyType( 5 ) ), std::out_of_range );
#endif

  BOOST_CHECK( ( m[KeyType( 5 )] = MappedType( 10 ) ) ); // check insertion through [] operator.
  BOOST_CHECK(
    ( m[KeyType( 5 )] == MappedType( 10 ) ) ); // validate the insertion (still has the expected mapped-value).

  BOOST_CHECK( m.insert( ValueType( KeyType( 10 ), MappedType( 20 ) ) ).second );

  std::vector< ValueType > tmp_v;
  tmp_v.push_back( ValueType( KeyType( 7 ), MappedType( 14 ) ) );
  tmp_v.push_back( ValueType( KeyType( 12 ), MappedType( 24 ) ) );
  tmp_v.push_back( ValueType( KeyType( 15 ), MappedType( 30 ) ) );
  tmp_v.push_back( ValueType( KeyType( 17 ), MappedType( 34 ) ) );
  tmp_v.push_back( ValueType( KeyType( 22 ), MappedType( 44 ) ) );
  m.insert( tmp_v.begin(), tmp_v.end() );
  BOOST_CHECK_MESSAGE( ( ( m[7] == 14 ) && ( m[12] == 24 ) && ( m[15] == 30 ) && ( m[17] == 34 ) && ( m[22] == 44 ) ),
                       "insert iterator range" );

#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  m.insert( {ValueType( KeyType( 8 ), MappedType( 16 ) ), ValueType( KeyType( 13 ), MappedType( 26 ) ),
             ValueType( KeyType( 16 ), MappedType( 32 ) )} );
  BOOST_CHECK_MESSAGE( ( ( m[8] == 16 ) && ( m[13] == 26 ) && ( m[16] == 32 ) ), "insert std::initializer_list" );
#endif
};

typedef boost::mpl::list< std::map< int, int >, std::multimap< int, int >, ReaK::graph::avlbfl_map< int, int >,
                          ReaK::graph::avlbfl_multimap< int, int >
                          //   , ReaK::graph::avlvebl_map<int, int>
                          //   , ReaK::graph::avlvebl_multimap<int, int>
                          > intint_multimaptest_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( intint_multimap_test, Map, intint_multimaptest_types ) {
  using std::swap;
  typedef typename Map::key_type KeyType;
  typedef typename Map::mapped_type MappedType;
  //   typedef typename Map::value_type ValueType;
  typedef std::pair< KeyType, MappedType > ValueType;
  typedef typename Map::iterator Iter;
  typedef typename Map::const_iterator ConstIter;
  typedef typename Map::const_reverse_iterator ConstRevIter;

  Map m;
  BOOST_CHECK_EQUAL( m.size(), 0 );
  BOOST_CHECK( m.empty() );
  BOOST_CHECK( m.begin() == m.end() );
  BOOST_CHECK( m.rbegin() == m.rend() );
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  BOOST_CHECK( m.cbegin() == m.cend() );
  BOOST_CHECK( m.crbegin() == m.crend() );
#endif

  m.insert( ValueType( KeyType( 5 ), MappedType( 10 ) ) );
  m.insert( ValueType( KeyType( 10 ), MappedType( 20 ) ) );
  BOOST_CHECK_MESSAGE( m.count( 10 ), "insert single element" );
  BOOST_CHECK( m.insert( m.end(), ValueType( KeyType( 20 ), MappedType( 40 ) ) ) != m.end() );
  BOOST_CHECK( m.insert( m.begin(), ValueType( KeyType( 4 ), MappedType( 8 ) ) ) != m.end() );

  std::vector< ValueType > tmp_v;
  tmp_v.push_back( ValueType( KeyType( 12 ), MappedType( 24 ) ) );
  tmp_v.push_back( ValueType( KeyType( 22 ), MappedType( 44 ) ) );
  tmp_v.push_back( ValueType( KeyType( 17 ), MappedType( 34 ) ) );
  tmp_v.push_back( ValueType( KeyType( 7 ), MappedType( 14 ) ) );
  tmp_v.push_back( ValueType( KeyType( 15 ), MappedType( 30 ) ) );

  m.insert( tmp_v.begin(), tmp_v.end() );
  BOOST_CHECK_MESSAGE( ( ( m.count( 7 ) == 1 ) && ( m.count( 12 ) == 1 ) && ( m.count( 15 ) == 1 )
                         && ( m.count( 17 ) == 1 ) && ( m.count( 22 ) == 1 ) ),
                       "insert iterator range" );

#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  m.insert( {ValueType( KeyType( 16 ), MappedType( 32 ) ), ValueType( KeyType( 8 ), MappedType( 16 ) ),
             ValueType( KeyType( 13 ), MappedType( 26 ) )} );
  BOOST_CHECK_MESSAGE( ( ( m.count( 8 ) == 1 ) && ( m.count( 13 ) == 1 ) && ( m.count( 16 ) == 1 ) ),
                       "insert std::initializer_list" );
#endif

  ConstIter it8_lo = m.lower_bound( 8 );
  ConstIter it8_hi = m.upper_bound( 8 );
  BOOST_CHECK( ( ( it8_lo != it8_hi ) && ( it8_lo != m.end() ) && ( std::distance( it8_lo, it8_hi ) == 1 ) ) );
  BOOST_CHECK( ( ( it8_lo->first == 8 ) && ( it8_lo->second == 16 ) ) );

  std::pair< ConstIter, ConstIter > it17_eq_range = m.equal_range( 17 );
  BOOST_CHECK( ( ( it17_eq_range.first != it17_eq_range.second ) && ( it17_eq_range.first != m.end() )
                 && ( std::distance( it17_eq_range.first, it17_eq_range.second ) == 1 ) ) );
  BOOST_CHECK( ( ( it17_eq_range.first->first == 17 ) && ( it17_eq_range.first->second == 34 ) ) );

  BOOST_CHECK_EQUAL( m.count( 15 ), 1 );

  ConstIter it15 = m.find( 15 );
  BOOST_CHECK( it15 != m.end() );
  BOOST_CHECK( ( ( it15->first == 15 ) && ( it15->second == 30 ) ) );

#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  Iter it16 = m.erase( it15 );
  BOOST_CHECK( m.count( 15 ) == 0 );
  BOOST_CHECK( it16 != m.end() );
  BOOST_CHECK( ( ( it16->first == 16 ) && ( it16->second == 32 ) ) );
#else
  m.erase( it15 );
  BOOST_CHECK( m.count( 15 ) == 0 );
#endif
  m.erase( m.find( 7 ), m.find( 13 ) );
  BOOST_CHECK_MESSAGE( ( ( m.count( 7 ) == 0 ) && ( m.count( 8 ) == 0 ) && ( m.count( 9 ) == 0 )
                         && ( m.count( 10 ) == 0 ) && ( m.count( 11 ) == 0 ) && ( m.count( 12 ) == 0 ) ),
                       "erase iterator range" );
  BOOST_CHECK( m.erase( 16 ) == 1 );

  // at this point the map contains: 4 5 13 17 20 22
  ConstIter it = m.begin();
  bool in_order = ( it != m.end() ) && ( it->first == 4 );
  ++it;
  in_order = in_order && ( it != m.end() ) && ( it->first == 5 );
  ++it;
  in_order = in_order && ( it != m.end() ) && ( it->first == 13 );
  ++it;
  in_order = in_order && ( it != m.end() ) && ( it->first == 17 );
  ++it;
  in_order = in_order && ( it != m.end() ) && ( it->first == 20 );
  ++it;
  in_order = in_order && ( it != m.end() ) && ( it->first == 22 );
  ++it;
  in_order = in_order && ( it == m.end() );
  BOOST_CHECK_MESSAGE( in_order, "in-order element traversal" );
  ConstRevIter rit = m.rbegin();
  in_order = ( rit != m.rend() ) && ( rit->first == 22 );
  ++rit;
  in_order = in_order && ( rit != m.rend() ) && ( rit->first == 20 );
  ++rit;
  in_order = in_order && ( rit != m.rend() ) && ( rit->first == 17 );
  ++rit;
  in_order = in_order && ( rit != m.rend() ) && ( rit->first == 13 );
  ++rit;
  in_order = in_order && ( rit != m.rend() ) && ( rit->first == 5 );
  ++rit;
  in_order = in_order && ( rit != m.rend() ) && ( rit->first == 4 );
  ++rit;
  in_order = in_order && ( rit == m.rend() );
  BOOST_CHECK_MESSAGE( in_order, "reverse in-order element traversal" );

  Map m2( tmp_v.begin(), tmp_v.end() );
  BOOST_CHECK_MESSAGE( ( ( m2.count( 7 ) == 1 ) && ( m2.count( 12 ) == 1 ) && ( m2.count( 15 ) == 1 )
                         && ( m2.count( 17 ) == 1 ) && ( m2.count( 22 ) == 1 ) ),
                       "constructor iterator range" );
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  Map m3{ValueType( KeyType( 8 ), MappedType( 16 ) ), ValueType( KeyType( 13 ), MappedType( 26 ) ),
         ValueType( KeyType( 16 ), MappedType( 32 ) )};
  BOOST_CHECK_MESSAGE( ( ( m3.count( 8 ) == 1 ) && ( m3.count( 13 ) == 1 ) && ( m3.count( 16 ) == 1 ) ),
                       "constructor std::initializer_list" );

  swap( m2, m3 );
  BOOST_CHECK_MESSAGE( ( ( m3.count( 7 ) == 1 ) && ( m3.count( 12 ) == 1 ) && ( m3.count( 15 ) == 1 )
                         && ( m3.count( 17 ) == 1 ) && ( m3.count( 22 ) == 1 ) && ( m2.count( 8 ) == 1 )
                         && ( m2.count( 13 ) == 1 ) && ( m2.count( 16 ) == 1 ) ),
                       "swap free function" );
#endif

  m.swap( m2 );
  BOOST_CHECK_MESSAGE( ( ( m.count( 8 ) == 1 ) && ( m.count( 13 ) == 1 ) && ( m.count( 16 ) == 1 ) ),
                       "swap member function" );

  m.clear();
  BOOST_CHECK_EQUAL( m.size(), 0 );
};


typedef boost::mpl::list< std::set< int >, std::multiset< int >, ReaK::graph::avlbfl_set< int >,
                          ReaK::graph::avlbfl_multiset< int >
                          //   , ReaK::graph::avlvebl_set<int>
                          //   , ReaK::graph::avlvebl_multiset<int>
                          > int_multisettest_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( int_multiset_test, Set, int_multisettest_types ) {
  using std::swap;
  typedef typename Set::value_type ValueType;
  typedef typename Set::iterator Iter;
  typedef typename Set::const_iterator ConstIter;
  typedef typename Set::const_reverse_iterator ConstRevIter;

  Set s;
  BOOST_CHECK_EQUAL( s.size(), 0 );
  BOOST_CHECK( s.empty() );
  BOOST_CHECK( s.begin() == s.end() );
  BOOST_CHECK( s.rbegin() == s.rend() );
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  BOOST_CHECK( s.cbegin() == s.cend() );
  BOOST_CHECK( s.crbegin() == s.crend() );
#endif

  s.insert( ValueType( 5 ) );
  BOOST_CHECK_MESSAGE( s.count( 5 ), "insert single element" );
  s.insert( ValueType( 10 ) );
  BOOST_CHECK_MESSAGE( s.count( 10 ), "insert single element" );
  BOOST_CHECK( s.insert( s.end(), ValueType( 20 ) ) != s.end() );
  BOOST_CHECK( s.insert( s.begin(), ValueType( 4 ) ) != s.end() );

  std::vector< ValueType > tmp_v;
  tmp_v.push_back( ValueType( 12 ) );
  tmp_v.push_back( ValueType( 22 ) );
  tmp_v.push_back( ValueType( 17 ) );
  tmp_v.push_back( ValueType( 7 ) );
  tmp_v.push_back( ValueType( 15 ) );

  s.insert( tmp_v.begin(), tmp_v.end() );
  BOOST_CHECK_MESSAGE( ( ( s.count( 7 ) == 1 ) && ( s.count( 12 ) == 1 ) && ( s.count( 15 ) == 1 )
                         && ( s.count( 17 ) == 1 ) && ( s.count( 22 ) == 1 ) ),
                       "insert iterator range" );

#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  s.insert( {ValueType( 16 ), ValueType( 8 ), ValueType( 13 )} );
  BOOST_CHECK_MESSAGE( ( ( s.count( 8 ) == 1 ) && ( s.count( 13 ) == 1 ) && ( s.count( 16 ) == 1 ) ),
                       "insert std::initializer_list" );
#endif

  ConstIter it8_lo = s.lower_bound( 8 );
  ConstIter it8_hi = s.upper_bound( 8 );
  BOOST_CHECK( ( ( it8_lo != it8_hi ) && ( it8_lo != s.end() ) && ( std::distance( it8_lo, it8_hi ) == 1 ) ) );
  BOOST_CHECK( ( *it8_lo == 8 ) );

  std::pair< ConstIter, ConstIter > it17_eq_range = s.equal_range( 17 );
  BOOST_CHECK( ( ( it17_eq_range.first != it17_eq_range.second ) && ( it17_eq_range.first != s.end() )
                 && ( std::distance( it17_eq_range.first, it17_eq_range.second ) == 1 ) ) );
  BOOST_CHECK( ( *( it17_eq_range.first ) == 17 ) );

  BOOST_CHECK_EQUAL( s.count( 15 ), 1 );

  ConstIter it15 = s.find( 15 );
  BOOST_CHECK( it15 != s.end() );
  BOOST_CHECK( ( *it15 == 15 ) );

#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  Iter it16 = s.erase( it15 );
  BOOST_CHECK( s.count( 15 ) == 0 );
  BOOST_CHECK( it16 != s.end() );
  BOOST_CHECK( ( *it16 == 16 ) );
#else
  s.erase( it15 );
  BOOST_CHECK( s.count( 15 ) == 0 );
#endif
  s.erase( s.find( 7 ), s.find( 13 ) );
  BOOST_CHECK_MESSAGE( ( ( s.count( 7 ) == 0 ) && ( s.count( 8 ) == 0 ) && ( s.count( 9 ) == 0 )
                         && ( s.count( 10 ) == 0 ) && ( s.count( 11 ) == 0 ) && ( s.count( 12 ) == 0 ) ),
                       "erase iterator range" );
  BOOST_CHECK( s.erase( 16 ) == 1 );

  // at this point the map contains: 4 5 13 17 20 22
  ConstIter it = s.begin();
  bool in_order = ( it != s.end() ) && ( *it == 4 );
  ++it;
  in_order = in_order && ( it != s.end() ) && ( *it == 5 );
  ++it;
  in_order = in_order && ( it != s.end() ) && ( *it == 13 );
  ++it;
  in_order = in_order && ( it != s.end() ) && ( *it == 17 );
  ++it;
  in_order = in_order && ( it != s.end() ) && ( *it == 20 );
  ++it;
  in_order = in_order && ( it != s.end() ) && ( *it == 22 );
  ++it;
  in_order = in_order && ( it == s.end() );
  BOOST_CHECK_MESSAGE( in_order, "in-order element traversal" );
  ConstRevIter rit = s.rbegin();
  in_order = ( rit != s.rend() ) && ( *rit == 22 );
  ++rit;
  in_order = in_order && ( rit != s.rend() ) && ( *rit == 20 );
  ++rit;
  in_order = in_order && ( rit != s.rend() ) && ( *rit == 17 );
  ++rit;
  in_order = in_order && ( rit != s.rend() ) && ( *rit == 13 );
  ++rit;
  in_order = in_order && ( rit != s.rend() ) && ( *rit == 5 );
  ++rit;
  in_order = in_order && ( rit != s.rend() ) && ( *rit == 4 );
  ++rit;
  in_order = in_order && ( rit == s.rend() );
  BOOST_CHECK_MESSAGE( in_order, "reverse in-order element traversal" );

  Set s2( tmp_v.begin(), tmp_v.end() );
  BOOST_CHECK_MESSAGE( ( ( s2.count( 7 ) == 1 ) && ( s2.count( 12 ) == 1 ) && ( s2.count( 15 ) == 1 )
                         && ( s2.count( 17 ) == 1 ) && ( s2.count( 22 ) == 1 ) ),
                       "constructor iterator range" );
#ifndef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
  Set s3{ValueType( 8 ), ValueType( 13 ), ValueType( 16 )};
  BOOST_CHECK_MESSAGE( ( ( s3.count( 8 ) == 1 ) && ( s3.count( 13 ) == 1 ) && ( s3.count( 16 ) == 1 ) ),
                       "constructor std::initializer_list" );

  swap( s2, s3 );
  BOOST_CHECK_MESSAGE( ( ( s3.count( 7 ) == 1 ) && ( s3.count( 12 ) == 1 ) && ( s3.count( 15 ) == 1 )
                         && ( s3.count( 17 ) == 1 ) && ( s3.count( 22 ) == 1 ) && ( s2.count( 8 ) == 1 )
                         && ( s2.count( 13 ) == 1 ) && ( s2.count( 16 ) == 1 ) ),
                       "swap free function" );
#endif

  s.swap( s2 );
  BOOST_CHECK_MESSAGE( ( ( s.count( 8 ) == 1 ) && ( s.count( 13 ) == 1 ) && ( s.count( 16 ) == 1 ) ),
                       "swap member function" );

  s.clear();
  BOOST_CHECK_EQUAL( s.size(), 0 );
};
