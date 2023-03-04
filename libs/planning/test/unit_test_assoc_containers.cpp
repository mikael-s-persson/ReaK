
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

#include <iomanip>
#include <iostream>

#include <map>
#include <set>

#include <ReaK/planning/graph_alg/avl_tree.hpp>

#include "gtest/gtest.h"

namespace ReaK::graph {
namespace {

template <typename T>
class IntIntMapTest : public ::testing::Test {};

using IntIntMapTestTypes =
    ::testing::Types<std::map<int, int>, avlbfl_map<int, int>>;
TYPED_TEST_SUITE(IntIntMapTest, IntIntMapTestTypes);

TYPED_TEST(IntIntMapTest, IntIntMapOperations) {
  using Map = TypeParam;
  using KeyType = typename Map::key_type;
  using MappedType = typename Map::mapped_type;
  using ValueType = std::pair<KeyType, MappedType>;

  Map m;
  EXPECT_THROW(m.at(KeyType(5)), std::out_of_range);

  EXPECT_TRUE((m[KeyType(5)] =
                   MappedType(10)));  // check insertion through [] operator.
  EXPECT_TRUE((
      m[KeyType(5)] ==
      MappedType(
          10)));  // validate the insertion (still has the expected mapped-value).

  EXPECT_TRUE(m.insert(ValueType(KeyType(10), MappedType(20))).second);

  std::vector<ValueType> tmp_v;
  tmp_v.push_back(ValueType(KeyType(7), MappedType(14)));
  tmp_v.push_back(ValueType(KeyType(12), MappedType(24)));
  tmp_v.push_back(ValueType(KeyType(15), MappedType(30)));
  tmp_v.push_back(ValueType(KeyType(17), MappedType(34)));
  tmp_v.push_back(ValueType(KeyType(22), MappedType(44)));
  m.insert(tmp_v.begin(), tmp_v.end());
  EXPECT_TRUE(((m[7] == 14) && (m[12] == 24) && (m[15] == 30) &&
               (m[17] == 34) && (m[22] == 44)))
      << "insert iterator range";

  m.insert({ValueType(KeyType(8), MappedType(16)),
            ValueType(KeyType(13), MappedType(26)),
            ValueType(KeyType(16), MappedType(32))});
  EXPECT_TRUE(((m[8] == 16) && (m[13] == 26) && (m[16] == 32)))
      << "insert std::initializer_list";
}

template <typename T>
class IntIntMultiMapTest : public ::testing::Test {};

using IntIntMultiMapTestTypes =
    ::testing::Types<std::map<int, int>, std::multimap<int, int>,
                     avlbfl_map<int, int>, avlbfl_multimap<int, int>>;
TYPED_TEST_SUITE(IntIntMultiMapTest, IntIntMultiMapTestTypes);

TYPED_TEST(IntIntMultiMapTest, IntIntMultiMapOperations) {
  using Map = TypeParam;
  using std::swap;
  using KeyType = typename Map::key_type;
  using MappedType = typename Map::mapped_type;
  //   typedef typename Map::value_type ValueType;
  using ValueType = std::pair<KeyType, MappedType>;
  using Iter = typename Map::iterator;
  using ConstIter = typename Map::const_iterator;
  using ConstRevIter = typename Map::const_reverse_iterator;

  Map m;
  EXPECT_EQ(m.size(), 0);
  EXPECT_TRUE(m.empty());
  EXPECT_TRUE(m.begin() == m.end());
  EXPECT_TRUE(m.rbegin() == m.rend());
  EXPECT_TRUE(m.cbegin() == m.cend());
  EXPECT_TRUE(m.crbegin() == m.crend());

  m.insert(ValueType(KeyType(5), MappedType(10)));
  m.insert(ValueType(KeyType(10), MappedType(20)));
  EXPECT_TRUE(m.count(10)) << "insert single element";
  EXPECT_TRUE(m.insert(m.end(), ValueType(KeyType(20), MappedType(40))) !=
              m.end());
  EXPECT_TRUE(m.insert(m.begin(), ValueType(KeyType(4), MappedType(8))) !=
              m.end());

  std::vector<ValueType> tmp_v;
  tmp_v.push_back(ValueType(KeyType(12), MappedType(24)));
  tmp_v.push_back(ValueType(KeyType(22), MappedType(44)));
  tmp_v.push_back(ValueType(KeyType(17), MappedType(34)));
  tmp_v.push_back(ValueType(KeyType(7), MappedType(14)));
  tmp_v.push_back(ValueType(KeyType(15), MappedType(30)));

  m.insert(tmp_v.begin(), tmp_v.end());
  EXPECT_TRUE((m.count(7) == 1) && (m.count(12) == 1) && (m.count(15) == 1) &&
              (m.count(17) == 1) && (m.count(22) == 1))
      << "insert iterator range";

  m.insert({ValueType(KeyType(16), MappedType(32)),
            ValueType(KeyType(8), MappedType(16)),
            ValueType(KeyType(13), MappedType(26))});
  EXPECT_TRUE((m.count(8) == 1) && (m.count(13) == 1) && (m.count(16) == 1))
      << "insert std::initializer_list";

  auto it8_lo = m.lower_bound(8);
  auto it8_hi = m.upper_bound(8);
  EXPECT_TRUE(((it8_lo != it8_hi) && (it8_lo != m.end()) &&
               (std::distance(it8_lo, it8_hi) == 1)));
  EXPECT_TRUE(((it8_lo->first == 8) && (it8_lo->second == 16)));

  std::pair<ConstIter, ConstIter> it17_eq_range = m.equal_range(17);
  EXPECT_TRUE(
      ((it17_eq_range.first != it17_eq_range.second) &&
       (it17_eq_range.first != m.end()) &&
       (std::distance(it17_eq_range.first, it17_eq_range.second) == 1)));
  EXPECT_TRUE(((it17_eq_range.first->first == 17) &&
               (it17_eq_range.first->second == 34)));

  EXPECT_EQ(m.count(15), 1);

  auto it15 = m.find(15);
  EXPECT_TRUE(it15 != m.end());
  EXPECT_TRUE(((it15->first == 15) && (it15->second == 30)));

  auto it16 = m.erase(it15);
  EXPECT_TRUE(m.count(15) == 0);
  EXPECT_TRUE(it16 != m.end());
  EXPECT_TRUE(((it16->first == 16) && (it16->second == 32)));
  m.erase(m.find(7), m.find(13));
  EXPECT_TRUE((m.count(7) == 0) && (m.count(8) == 0) && (m.count(9) == 0) &&
              (m.count(10) == 0) && (m.count(11) == 0) && (m.count(12) == 0))
      << "erase iterator range";
  EXPECT_TRUE(m.erase(16) == 1);

  // at this point the map contains: 4 5 13 17 20 22
  auto it = m.begin();
  bool in_order = (it != m.end()) && (it->first == 4);
  ++it;
  in_order = in_order && (it != m.end()) && (it->first == 5);
  ++it;
  in_order = in_order && (it != m.end()) && (it->first == 13);
  ++it;
  in_order = in_order && (it != m.end()) && (it->first == 17);
  ++it;
  in_order = in_order && (it != m.end()) && (it->first == 20);
  ++it;
  in_order = in_order && (it != m.end()) && (it->first == 22);
  ++it;
  in_order = in_order && (it == m.end());
  EXPECT_TRUE(in_order) << "in-order element traversal";
  auto rit = m.rbegin();
  in_order = (rit != m.rend()) && (rit->first == 22);
  ++rit;
  in_order = in_order && (rit != m.rend()) && (rit->first == 20);
  ++rit;
  in_order = in_order && (rit != m.rend()) && (rit->first == 17);
  ++rit;
  in_order = in_order && (rit != m.rend()) && (rit->first == 13);
  ++rit;
  in_order = in_order && (rit != m.rend()) && (rit->first == 5);
  ++rit;
  in_order = in_order && (rit != m.rend()) && (rit->first == 4);
  ++rit;
  in_order = in_order && (rit == m.rend());
  EXPECT_TRUE(in_order) << "reverse in-order element traversal";

  Map m2(tmp_v.begin(), tmp_v.end());
  EXPECT_TRUE((m2.count(7) == 1) && (m2.count(12) == 1) &&
              (m2.count(15) == 1) && (m2.count(17) == 1) && (m2.count(22) == 1))
      << "constructor iterator range";
  Map m3{ValueType(KeyType(8), MappedType(16)),
         ValueType(KeyType(13), MappedType(26)),
         ValueType(KeyType(16), MappedType(32))};
  EXPECT_TRUE((m3.count(8) == 1) && (m3.count(13) == 1) && (m3.count(16) == 1))
      << "constructor std::initializer_list";

  swap(m2, m3);
  EXPECT_TRUE((m3.count(7) == 1) && (m3.count(12) == 1) &&
              (m3.count(15) == 1) && (m3.count(17) == 1) &&
              (m3.count(22) == 1) && (m2.count(8) == 1) &&
              (m2.count(13) == 1) && (m2.count(16) == 1))
      << "swap free function";

  m.swap(m2);
  EXPECT_TRUE((m.count(8) == 1) && (m.count(13) == 1) && (m.count(16) == 1))
      << "swap member function";

  m.clear();
  EXPECT_EQ(m.size(), 0);
}

template <typename T>
class IntIntMultiSetTest : public ::testing::Test {};

using IntIntMultiSetTestTypes =
    ::testing::Types<std::set<int>, std::multiset<int>, avlbfl_set<int>,
                     avlbfl_multiset<int>>;
TYPED_TEST_SUITE(IntIntMultiSetTest, IntIntMultiSetTestTypes);

TYPED_TEST(IntIntMultiSetTest, IntIntMultiSetOperations) {
  using Set = TypeParam;
  using std::swap;
  using ValueType = typename Set::value_type;
  using Iter = typename Set::iterator;
  using ConstIter = typename Set::const_iterator;
  using ConstRevIter = typename Set::const_reverse_iterator;

  Set s;
  EXPECT_EQ(s.size(), 0);
  EXPECT_TRUE(s.empty());
  EXPECT_TRUE(s.begin() == s.end());
  EXPECT_TRUE(s.rbegin() == s.rend());
  EXPECT_TRUE(s.cbegin() == s.cend());
  EXPECT_TRUE(s.crbegin() == s.crend());

  s.insert(ValueType(5));
  EXPECT_TRUE(s.count(5)) << "insert single element";
  s.insert(ValueType(10));
  EXPECT_TRUE(s.count(10)) << "insert single element";
  EXPECT_TRUE(s.insert(s.end(), ValueType(20)) != s.end());
  EXPECT_TRUE(s.insert(s.begin(), ValueType(4)) != s.end());

  std::vector<ValueType> tmp_v;
  tmp_v.push_back(ValueType(12));
  tmp_v.push_back(ValueType(22));
  tmp_v.push_back(ValueType(17));
  tmp_v.push_back(ValueType(7));
  tmp_v.push_back(ValueType(15));

  s.insert(tmp_v.begin(), tmp_v.end());
  EXPECT_TRUE((s.count(7) == 1) && (s.count(12) == 1) && (s.count(15) == 1) &&
              (s.count(17) == 1) && (s.count(22) == 1))
      << "insert iterator range";

  s.insert({ValueType(16), ValueType(8), ValueType(13)});
  EXPECT_TRUE((s.count(8) == 1) && (s.count(13) == 1) && (s.count(16) == 1))
      << "insert std::initializer_list";

  auto it8_lo = s.lower_bound(8);
  auto it8_hi = s.upper_bound(8);
  EXPECT_TRUE(((it8_lo != it8_hi) && (it8_lo != s.end()) &&
               (std::distance(it8_lo, it8_hi) == 1)));
  EXPECT_TRUE((*it8_lo == 8));

  std::pair<ConstIter, ConstIter> it17_eq_range = s.equal_range(17);
  EXPECT_TRUE(
      ((it17_eq_range.first != it17_eq_range.second) &&
       (it17_eq_range.first != s.end()) &&
       (std::distance(it17_eq_range.first, it17_eq_range.second) == 1)));
  EXPECT_TRUE((*(it17_eq_range.first) == 17));

  EXPECT_EQ(s.count(15), 1);

  auto it15 = s.find(15);
  EXPECT_TRUE(it15 != s.end());
  EXPECT_TRUE((*it15 == 15));

  auto it16 = s.erase(it15);
  EXPECT_TRUE(s.count(15) == 0);
  EXPECT_TRUE(it16 != s.end());
  EXPECT_TRUE((*it16 == 16));
  s.erase(s.find(7), s.find(13));
  EXPECT_TRUE((s.count(7) == 0) && (s.count(8) == 0) && (s.count(9) == 0) &&
              (s.count(10) == 0) && (s.count(11) == 0) && (s.count(12) == 0))
      << "erase iterator range";
  EXPECT_TRUE(s.erase(16) == 1);

  // at this point the map contains: 4 5 13 17 20 22
  auto it = s.begin();
  bool in_order = (it != s.end()) && (*it == 4);
  ++it;
  in_order = in_order && (it != s.end()) && (*it == 5);
  ++it;
  in_order = in_order && (it != s.end()) && (*it == 13);
  ++it;
  in_order = in_order && (it != s.end()) && (*it == 17);
  ++it;
  in_order = in_order && (it != s.end()) && (*it == 20);
  ++it;
  in_order = in_order && (it != s.end()) && (*it == 22);
  ++it;
  in_order = in_order && (it == s.end());
  EXPECT_TRUE(in_order) << "in-order element traversal";
  auto rit = s.rbegin();
  in_order = (rit != s.rend()) && (*rit == 22);
  ++rit;
  in_order = in_order && (rit != s.rend()) && (*rit == 20);
  ++rit;
  in_order = in_order && (rit != s.rend()) && (*rit == 17);
  ++rit;
  in_order = in_order && (rit != s.rend()) && (*rit == 13);
  ++rit;
  in_order = in_order && (rit != s.rend()) && (*rit == 5);
  ++rit;
  in_order = in_order && (rit != s.rend()) && (*rit == 4);
  ++rit;
  in_order = in_order && (rit == s.rend());
  EXPECT_TRUE(in_order) << "reverse in-order element traversal";

  Set s2(tmp_v.begin(), tmp_v.end());
  EXPECT_TRUE((s2.count(7) == 1) && (s2.count(12) == 1) &&
              (s2.count(15) == 1) && (s2.count(17) == 1) && (s2.count(22) == 1))
      << "constructor iterator range";
  Set s3{ValueType(8), ValueType(13), ValueType(16)};
  EXPECT_TRUE((s3.count(8) == 1) && (s3.count(13) == 1) && (s3.count(16) == 1))
      << "constructor std::initializer_list";

  swap(s2, s3);
  EXPECT_TRUE((s3.count(7) == 1) && (s3.count(12) == 1) &&
              (s3.count(15) == 1) && (s3.count(17) == 1) &&
              (s3.count(22) == 1) && (s2.count(8) == 1) &&
              (s2.count(13) == 1) && (s2.count(16) == 1))
      << "swap free function";

  s.swap(s2);
  EXPECT_TRUE((s.count(8) == 1) && (s.count(13) == 1) && (s.count(16) == 1))
      << "swap member function";

  s.clear();
  EXPECT_EQ(s.size(), 0);
}

}  // namespace
}  // namespace ReaK::graph
