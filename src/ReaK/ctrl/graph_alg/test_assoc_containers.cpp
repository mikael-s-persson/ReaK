
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

#include "graph_alg/avl_tree.hpp"

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE assoc_containers
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>



typedef boost::mpl::list< std::map<int, int> > intint_maptest_types;


BOOST_AUTO_TEST_CASE_TEMPLATE( intint_map_test, Map, intint_maptest_types )
{
  typedef typename Map::key_type KeyType;
  typedef typename Map::mapped_type MappedType;
  typedef typename Map::value_type ValueType;
  typedef typename Map::iterator Iter;
  typedef typename Map::const_iterator ConstIter;
  typedef typename Map::reverse_iterator RevIter;
  typedef typename Map::const_reverse_iterator ConstRevIter;
  
  Map m;
#ifdef RK_ENABLE_CXX11_FEATURES
  BOOST_CHECK_THROW( m.at(KeyType(5)), std::out_of_range);
#endif
  
  BOOST_CHECK( (m[KeyType(5)] =  MappedType(10)) );  // check insertion through [] operator.
  BOOST_CHECK( (m[KeyType(5)] == MappedType(10)) );  // validate the insertion (still has the expected mapped-value).
  
  BOOST_CHECK( m.insert( ValueType(KeyType(10),MappedType(20))).second );
  
  std::vector< ValueType > tmp_v;
  tmp_v.push_back(ValueType(KeyType(7), MappedType(14)));
  tmp_v.push_back(ValueType(KeyType(12),MappedType(24)));
  tmp_v.push_back(ValueType(KeyType(15),MappedType(30)));
  tmp_v.push_back(ValueType(KeyType(17),MappedType(34)));
  tmp_v.push_back(ValueType(KeyType(22),MappedType(44)));
  m.insert(tmp_v.begin(), tmp_v.end());
  BOOST_CHECK_MESSAGE( ((m[7] == 14) && (m[12] == 24) && (m[15] == 30) && (m[17] == 34) && (m[22] == 44)), "insert iterator range" );
  
#ifdef RK_ENABLE_CXX11_FEATURES
  m.insert({ValueType(KeyType(8), MappedType(16)), ValueType(KeyType(13), MappedType(26)), ValueType(KeyType(16), MappedType(32))});
  BOOST_CHECK_MESSAGE( ((m[8] == 16) && (m[13] == 26) && (m[16] == 32)), "insert std::initializer_list" );
#endif
  
};

// typedef boost::mpl::list< std::map<int, int>, std::multimap<int, int>, ReaK::graph::avlbf_map<int, int> > intint_multimaptest_types;
typedef boost::mpl::list< ReaK::graph::avlbf_map<int, int> > intint_multimaptest_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( intint_multimap_test, Map, intint_multimaptest_types )
{
  using std::swap;
  typedef typename Map::key_type KeyType;
  typedef typename Map::mapped_type MappedType;
  typedef typename Map::value_type ValueType;
  typedef typename Map::iterator Iter;
  typedef typename Map::const_iterator ConstIter;
  typedef typename Map::reverse_iterator RevIter;
  typedef typename Map::const_reverse_iterator ConstRevIter;
  
  Map m;
  BOOST_CHECK_EQUAL( m.size(), 0);
  BOOST_CHECK( m.empty() ); 
  BOOST_CHECK( m.begin() == m.end() ); 
  BOOST_CHECK( m.rbegin() == m.rend() ); 
#ifdef RK_ENABLE_CXX11_FEATURES
  BOOST_CHECK( m.cbegin() == m.cend() ); 
  BOOST_CHECK( m.crbegin() == m.crend() ); 
#endif
  
  m.insert( ValueType(KeyType(5),MappedType(10))); 
  m.insert( ValueType(KeyType(10),MappedType(20))); 
  BOOST_CHECK_MESSAGE( m.count(10), "insert single element" ); 
  BOOST_CHECK( m.insert( m.end(), ValueType(KeyType(20),MappedType(40))) != m.end() ); 
  BOOST_CHECK( m.insert( m.begin(), ValueType(KeyType(4),MappedType(8))) != m.end() ); 
  
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
  
  std::vector< ValueType > tmp_v; 
  tmp_v.push_back(ValueType(KeyType(12),MappedType(24)));
  tmp_v.push_back(ValueType(KeyType(22),MappedType(44)));
  tmp_v.push_back(ValueType(KeyType(17),MappedType(34)));
  tmp_v.push_back(ValueType(KeyType(7), MappedType(14)));
  tmp_v.push_back(ValueType(KeyType(15),MappedType(30)));
  
  m.insert(tmp_v.begin(), tmp_v.end()); 
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
  
  BOOST_CHECK_MESSAGE( ((m.count(7) == 1) && (m.count(12) == 1) && (m.count(15) == 1) && (m.count(17) == 1) && (m.count(22) == 1)), "insert iterator range" );
   
#ifdef RK_ENABLE_CXX11_FEATURES
  m.insert({ValueType(KeyType(16), MappedType(32)), ValueType(KeyType(8), MappedType(16)), ValueType(KeyType(13), MappedType(26))});
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
  BOOST_CHECK_MESSAGE( ((m.count(8) == 1) && (m.count(13) == 1) && (m.count(16) == 1)), "insert std::initializer_list" );
#endif
   
  ConstIter it8_lo = m.lower_bound(8);
  ConstIter it8_hi = m.upper_bound(8);
  BOOST_CHECK( ((it8_lo != it8_hi) && (it8_lo != m.end()) && (std::distance(it8_lo, it8_hi) == 1)) );
  BOOST_CHECK( ((it8_lo->first == 8) && (it8_lo->second == 16)) );
   
  std::pair< ConstIter, ConstIter > it17_eq_range = m.equal_range(17);
  BOOST_CHECK( ((it17_eq_range.first != it17_eq_range.second) && (it17_eq_range.first != m.end()) && (std::distance(it17_eq_range.first, it17_eq_range.second) == 1)) );
  BOOST_CHECK( ((it17_eq_range.first->first == 17) && (it17_eq_range.first->second == 34)) );
   
  BOOST_CHECK_EQUAL( m.count(15), 1 );
   
  ConstIter it15 = m.find(15);
  BOOST_CHECK( it15 != m.end() );
  BOOST_CHECK( ((it15->first == 15) && (it15->second == 30)) );
   
#ifdef RK_ENABLE_CXX11_FEATURES
  Iter it16 = m.erase(it15);
  BOOST_CHECK( m.count(15) == 0 );
  BOOST_CHECK( it16 != m.end() );
  BOOST_CHECK( ((it16->first == 16) && (it16->second == 32)) ); 
#else
  m.erase(it15);
  BOOST_CHECK( m.count(15) == 0 ); 
#endif
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
  m.erase( m.find(7), m.find(13) );
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
  BOOST_CHECK_MESSAGE( ((m.count(7) == 0) && (m.count(8) == 0) && (m.count(9) == 0) && (m.count(10) == 0) && (m.count(11) == 0) && (m.count(12) == 0)), "erase iterator range" );
  BOOST_CHECK( m.erase(16) == 1 );
  for(ConstIter it = m.begin(); it != m.end(); ++it) {
    std::cout << std::setw(4) << it->first << std::setw(4) << it->second << std::endl;
  };
   
  // at this point the map contains: 4 5 13 17 20 22
  ConstIter it = m.begin();
  bool in_order  = (it != m.end()) && (it->first == 4); 
  ++it; in_order = in_order && (it != m.end()) && (it->first == 5); 
  ++it; in_order = in_order && (it != m.end()) && (it->first == 13); 
  ++it; in_order = in_order && (it != m.end()) && (it->first == 17); 
  ++it; in_order = in_order && (it != m.end()) && (it->first == 20); 
  ++it; in_order = in_order && (it != m.end()) && (it->first == 22); 
  ++it; in_order = in_order && (it == m.end()); 
  BOOST_CHECK_MESSAGE( in_order, "in-order element traversal"); 
  for(ConstRevIter rit = m.rbegin(); rit != m.rend(); ++rit) {
    std::cout << std::setw(4) << rit->first << std::setw(4) << rit->second << std::endl;
  };
  ConstRevIter rit = m.rbegin();
  in_order  = (rit != m.rend()) && (rit->first == 22); 
  ++rit; in_order = in_order && (rit != m.rend()) && (rit->first == 20); 
  ++rit; in_order = in_order && (rit != m.rend()) && (rit->first == 17); 
  ++rit; in_order = in_order && (rit != m.rend()) && (rit->first == 13); 
  ++rit; in_order = in_order && (rit != m.rend()) && (rit->first == 5); 
  ++rit; in_order = in_order && (rit != m.rend()) && (rit->first == 4); 
  ++rit; in_order = in_order && (rit == m.rend()); 
  BOOST_CHECK_MESSAGE( in_order, "reverse in-order element traversal"); 
  
  Map m2(tmp_v.begin(), tmp_v.end()); 
  BOOST_CHECK_MESSAGE( ((m2.count(7) == 1) && (m2.count(12) == 1) && (m2.count(15) == 1) && (m2.count(17) == 1) && (m2.count(22) == 1)), "constructor iterator range" );
#ifdef RK_ENABLE_CXX11_FEATURES
  Map m3{ValueType(KeyType(8), MappedType(16)), ValueType(KeyType(13), MappedType(26)), ValueType(KeyType(16), MappedType(32))};
  BOOST_CHECK_MESSAGE( ((m3.count(8) == 1) && (m3.count(13) == 1) && (m3.count(16) == 1)), "constructor std::initializer_list" ); 
#endif
  
  swap(m2, m3); 
  BOOST_CHECK_MESSAGE( ((m3.count(7) == 1) && (m3.count(12) == 1) && (m3.count(15) == 1) && (m3.count(17) == 1) && (m3.count(22) == 1) && (m2.count(8) == 1) && (m2.count(13) == 1) && (m2.count(16) == 1)), "swap free function" );
  
  m.swap(m2); 
  BOOST_CHECK_MESSAGE( ((m.count(8) == 1) && (m.count(13) == 1) && (m.count(16) == 1)), "swap member function" );
  
  m.clear(); 
  BOOST_CHECK_EQUAL( m.size(), 0);
  
};


// typedef boost::mpl::list< std::map<int, int> > intint_test_types;



