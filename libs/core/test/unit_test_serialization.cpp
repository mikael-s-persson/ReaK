
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

#include "ReaK/core/base/named_object.hpp"

#include "ReaK/core/serialization/bin_archiver.hpp"
#include "ReaK/core/serialization/objtree_archiver.hpp"
#include "ReaK/core/serialization/protobuf_archiver.hpp"
#include "ReaK/core/serialization/xml_archiver.hpp"

#include <sstream>

#include "gtest/gtest.h"

namespace ReaK {

class obj_with_named_members : public ReaK::named_object {
 public:
  unsigned int m_uint;
  int m_int;
  float m_float;
  double m_double;
  char m_char;
  bool m_bool;
  std::string m_str;
  std::vector<int> m_vect;
  std::list<int> m_list;
  std::set<int> m_set;
  std::map<int, std::string> m_map;

  obj_with_named_members() {
    setName("object_with_named_data_members");

    m_uint = 42;
    m_int = 69;
    m_float = 1.5;
    m_double = 9.5;
    m_char = 'f';
    m_bool = true;
    m_str = "dude";

    m_vect.push_back(89);
    m_vect.push_back(78);
    m_vect.push_back(67);
    m_vect.push_back(56);

    m_list.push_back(45);
    m_list.push_back(34);
    m_list.push_back(23);
    m_list.push_back(12);

    m_set.insert(5);
    m_set.insert(10);
    m_set.insert(15);
    m_set.insert(20);

    m_map.insert(std::pair<int, std::string>(5, "five"));
    m_map.insert(std::pair<int, std::string>(10, "ten"));
    m_map.insert(std::pair<int, std::string>(15, "fifteen"));
    m_map.insert(std::pair<int, std::string>(20, "twenty"));
  }

  bool check_uint() const { return m_uint == 42; }
  bool check_int() const { return m_int == 69; }
  bool check_float() const { return std::abs(m_float - 1.5) < 1e-6; }
  bool check_double() const { return std::abs(m_double - 9.5) < 1e-6; }
  bool check_char() const { return m_char == 'f'; }
  bool check_bool() const { return m_bool; }
  bool check_str() const { return m_str == "dude"; }
  bool check_vect() const {
    if (m_vect.size() != 4) {
      return false;
    }
    return ((m_vect[0] == 89) && (m_vect[1] == 78) && (m_vect[2] == 67) &&
            (m_vect[3] == 56));
  }
  bool check_list() const {
    if (m_list.size() != 4) {
      return false;
    }
    auto it = m_list.begin();
    if (*it != 45) {
      return false;
    }
    ++it;
    if (*it != 34) {
      return false;
    }
    ++it;
    if (*it != 23) {
      return false;
    }
    ++it;
    return *it == 12;
  }
  bool check_set() const {
    if (m_set.size() != 4) {
      return false;
    }
    auto it = m_set.begin();
    if (*it != 5) {
      return false;
    }
    ++it;
    if (*it != 10) {
      return false;
    }
    ++it;
    if (*it != 15) {
      return false;
    }
    ++it;
    return *it == 20;
  }
  bool check_map() const {
    if (m_map.size() != 4) {
      return false;
    }
    auto it = m_map.begin();
    if ((it->first != 5) || (it->second != "five")) {
      return false;
    }
    ++it;
    if ((it->first != 10) || (it->second != "ten")) {
      return false;
    }
    ++it;
    if ((it->first != 15) || (it->second != "fifteen")) {
      return false;
    }
    ++it;
    return !((it->first != 20) || (it->second != "twenty"));
  }

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, ReaK::named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_uint) & RK_SERIAL_SAVE_WITH_NAME(m_int) &
        RK_SERIAL_SAVE_WITH_NAME(m_float) & RK_SERIAL_SAVE_WITH_NAME(m_double) &
        RK_SERIAL_SAVE_WITH_NAME(m_char) & RK_SERIAL_SAVE_WITH_NAME(m_bool) &
        RK_SERIAL_SAVE_WITH_NAME(m_str) & RK_SERIAL_SAVE_WITH_NAME(m_vect) &
        RK_SERIAL_SAVE_WITH_NAME(m_list) & RK_SERIAL_SAVE_WITH_NAME(m_set) &
        RK_SERIAL_SAVE_WITH_NAME(m_map);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, ReaK::named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_uint) & RK_SERIAL_LOAD_WITH_NAME(m_int) &
        RK_SERIAL_LOAD_WITH_NAME(m_float) & RK_SERIAL_LOAD_WITH_NAME(m_double) &
        RK_SERIAL_LOAD_WITH_NAME(m_char) & RK_SERIAL_LOAD_WITH_NAME(m_bool) &
        RK_SERIAL_LOAD_WITH_NAME(m_str) & RK_SERIAL_LOAD_WITH_NAME(m_vect) &
        RK_SERIAL_LOAD_WITH_NAME(m_list) & RK_SERIAL_LOAD_WITH_NAME(m_set) &
        RK_SERIAL_LOAD_WITH_NAME(m_map);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(obj_with_named_members, 0xFFFFFFFE, 1,
                              "obj_with_named_members", ReaK::named_object)
};

class obj_with_unnamed_members : public ReaK::named_object {
 public:
  unsigned int m_uint;
  int m_int;
  float m_float;
  double m_double;
  char m_char;
  bool m_bool;
  std::string m_str;
  std::vector<int> m_vect;
  std::list<int> m_list;
  std::set<int> m_set;
  std::map<int, std::string> m_map;

  obj_with_unnamed_members() {
    setName("object_with_unnamed_data_members");

    m_uint = 42;
    m_int = 69;
    m_float = 1.5;
    m_double = 9.5;
    m_char = 'f';
    m_bool = true;
    m_str = "dude";

    m_vect.push_back(89);
    m_vect.push_back(78);
    m_vect.push_back(67);
    m_vect.push_back(56);

    m_list.push_back(45);
    m_list.push_back(34);
    m_list.push_back(23);
    m_list.push_back(12);

    m_set.insert(5);
    m_set.insert(10);
    m_set.insert(15);
    m_set.insert(20);

    m_map.insert(std::pair<int, std::string>(5, "five"));
    m_map.insert(std::pair<int, std::string>(10, "ten"));
    m_map.insert(std::pair<int, std::string>(15, "fifteen"));
    m_map.insert(std::pair<int, std::string>(20, "twenty"));
  }

  bool check_uint() const { return m_uint == 42; }
  bool check_int() const { return m_int == 69; }
  bool check_float() const { return std::abs(m_float - 1.5) < 1e-6; }
  bool check_double() const { return std::abs(m_double - 9.5) < 1e-6; }
  bool check_char() const { return m_char == 'f'; }
  bool check_bool() const { return m_bool; }
  bool check_str() const { return m_str == "dude"; }
  bool check_vect() const {
    if (m_vect.size() != 4) {
      return false;
    }
    return ((m_vect[0] == 89) && (m_vect[1] == 78) && (m_vect[2] == 67) &&
            (m_vect[3] == 56));
  }
  bool check_list() const {
    if (m_list.size() != 4) {
      return false;
    }
    auto it = m_list.begin();
    if (*it != 45) {
      return false;
    }
    ++it;
    if (*it != 34) {
      return false;
    }
    ++it;
    if (*it != 23) {
      return false;
    }
    ++it;
    return *it == 12;
  }
  bool check_set() const {
    if (m_set.size() != 4) {
      return false;
    }
    auto it = m_set.begin();
    if (*it != 5) {
      return false;
    }
    ++it;
    if (*it != 10) {
      return false;
    }
    ++it;
    if (*it != 15) {
      return false;
    }
    ++it;
    return *it == 20;
  }
  bool check_map() const {
    if (m_map.size() != 4) {
      return false;
    }
    auto it = m_map.begin();
    if ((it->first != 5) || (it->second != "five")) {
      return false;
    }
    ++it;
    if ((it->first != 10) || (it->second != "ten")) {
      return false;
    }
    ++it;
    if ((it->first != 15) || (it->second != "fifteen")) {
      return false;
    }
    ++it;
    return !((it->first != 20) || (it->second != "twenty"));
  }

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, ReaK::named_object::getStaticObjectType()->TypeVersion());
    A << m_uint << m_int << m_float << m_double << m_char << m_bool << m_str
      << m_vect << m_list << m_set << m_map;
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, ReaK::named_object::getStaticObjectType()->TypeVersion());
    A >> m_uint >> m_int >> m_float >> m_double >> m_char >> m_bool >> m_str >>
        m_vect >> m_list >> m_set >> m_map;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(obj_with_unnamed_members, 0xFFFFFFFF, 1,
                              "obj_with_unnamed_members", ReaK::named_object)
};

}  // namespace ReaK

namespace ReaK::serialization {
namespace {

TEST(SerializationTests, BinSerializers) {
  {
    std::stringstream ss;
    {
      bin_oarchive output_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names(
          new obj_with_named_members());
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names(
          new obj_with_unnamed_members());

      EXPECT_NO_THROW(output_arc << obj_with_names);
      EXPECT_NO_THROW(output_arc << ptr_with_names);
      EXPECT_NO_THROW(output_arc << obj_with_no_names);
      EXPECT_NO_THROW(output_arc << ptr_with_no_names);
    }

    {
      bin_iarchive input_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names;
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names;

      EXPECT_NO_THROW(input_arc >> obj_with_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_names);
      EXPECT_NO_THROW(input_arc >> obj_with_no_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_no_names);

      EXPECT_TRUE(obj_with_names.check_uint());
      EXPECT_TRUE(obj_with_names.check_int());
      EXPECT_TRUE(obj_with_names.check_float());
      EXPECT_TRUE(obj_with_names.check_double());
      EXPECT_TRUE(obj_with_names.check_char());
      EXPECT_TRUE(obj_with_names.check_bool());
      EXPECT_TRUE(obj_with_names.check_str());
      EXPECT_TRUE(obj_with_names.check_vect());
      EXPECT_TRUE(obj_with_names.check_list());
      EXPECT_TRUE(obj_with_names.check_set());
      EXPECT_TRUE(obj_with_names.check_map());

      EXPECT_TRUE(ptr_with_names);
      EXPECT_TRUE(ptr_with_names->check_uint());
      EXPECT_TRUE(ptr_with_names->check_int());
      EXPECT_TRUE(ptr_with_names->check_float());
      EXPECT_TRUE(ptr_with_names->check_double());
      EXPECT_TRUE(ptr_with_names->check_char());
      EXPECT_TRUE(ptr_with_names->check_bool());
      EXPECT_TRUE(ptr_with_names->check_str());
      EXPECT_TRUE(ptr_with_names->check_vect());
      EXPECT_TRUE(ptr_with_names->check_list());
      EXPECT_TRUE(ptr_with_names->check_set());
      EXPECT_TRUE(ptr_with_names->check_map());

      EXPECT_TRUE(obj_with_no_names.check_uint());
      EXPECT_TRUE(obj_with_no_names.check_int());
      EXPECT_TRUE(obj_with_no_names.check_float());
      EXPECT_TRUE(obj_with_no_names.check_double());
      EXPECT_TRUE(obj_with_no_names.check_char());
      EXPECT_TRUE(obj_with_no_names.check_bool());
      EXPECT_TRUE(obj_with_no_names.check_str());
      EXPECT_TRUE(obj_with_no_names.check_vect());
      EXPECT_TRUE(obj_with_no_names.check_list());
      EXPECT_TRUE(obj_with_no_names.check_set());
      EXPECT_TRUE(obj_with_no_names.check_map());

      EXPECT_TRUE(ptr_with_no_names);
      EXPECT_TRUE(ptr_with_no_names->check_uint());
      EXPECT_TRUE(ptr_with_no_names->check_int());
      EXPECT_TRUE(ptr_with_no_names->check_float());
      EXPECT_TRUE(ptr_with_no_names->check_double());
      EXPECT_TRUE(ptr_with_no_names->check_char());
      EXPECT_TRUE(ptr_with_no_names->check_bool());
      EXPECT_TRUE(ptr_with_no_names->check_str());
      EXPECT_TRUE(ptr_with_no_names->check_vect());
      EXPECT_TRUE(ptr_with_no_names->check_list());
      EXPECT_TRUE(ptr_with_no_names->check_set());
      EXPECT_TRUE(ptr_with_no_names->check_map());
    }
  }
}

TEST(SerializationTests, XmlSerializers) {
  {
    std::stringstream ss;
    {
      xml_oarchive output_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names(
          new obj_with_named_members());
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names(
          new obj_with_unnamed_members());

      EXPECT_NO_THROW(output_arc << obj_with_names);
      EXPECT_NO_THROW(output_arc << ptr_with_names);
      EXPECT_NO_THROW(output_arc << obj_with_no_names);
      EXPECT_NO_THROW(output_arc << ptr_with_no_names);
    }

    {
      xml_iarchive input_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names;
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names;

      EXPECT_NO_THROW(input_arc >> obj_with_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_names);
      EXPECT_NO_THROW(input_arc >> obj_with_no_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_no_names);

      EXPECT_TRUE(obj_with_names.check_uint());
      EXPECT_TRUE(obj_with_names.check_int());
      EXPECT_TRUE(obj_with_names.check_float());
      EXPECT_TRUE(obj_with_names.check_double());
      EXPECT_TRUE(obj_with_names.check_char());
      EXPECT_TRUE(obj_with_names.check_bool());
      EXPECT_TRUE(obj_with_names.check_str());
      EXPECT_TRUE(obj_with_names.check_vect());
      EXPECT_TRUE(obj_with_names.check_list());
      EXPECT_TRUE(obj_with_names.check_set());
      EXPECT_TRUE(obj_with_names.check_map());

      EXPECT_TRUE(ptr_with_names);
      EXPECT_TRUE(ptr_with_names->check_uint());
      EXPECT_TRUE(ptr_with_names->check_int());
      EXPECT_TRUE(ptr_with_names->check_float());
      EXPECT_TRUE(ptr_with_names->check_double());
      EXPECT_TRUE(ptr_with_names->check_char());
      EXPECT_TRUE(ptr_with_names->check_bool());
      EXPECT_TRUE(ptr_with_names->check_str());
      EXPECT_TRUE(ptr_with_names->check_vect());
      EXPECT_TRUE(ptr_with_names->check_list());
      EXPECT_TRUE(ptr_with_names->check_set());
      EXPECT_TRUE(ptr_with_names->check_map());

      EXPECT_TRUE(obj_with_no_names.check_uint());
      EXPECT_TRUE(obj_with_no_names.check_int());
      EXPECT_TRUE(obj_with_no_names.check_float());
      EXPECT_TRUE(obj_with_no_names.check_double());
      EXPECT_TRUE(obj_with_no_names.check_char());
      EXPECT_TRUE(obj_with_no_names.check_bool());
      EXPECT_TRUE(obj_with_no_names.check_str());
      EXPECT_TRUE(obj_with_no_names.check_vect());
      EXPECT_TRUE(obj_with_no_names.check_list());
      EXPECT_TRUE(obj_with_no_names.check_set());
      EXPECT_TRUE(obj_with_no_names.check_map());

      EXPECT_TRUE(ptr_with_no_names);
      EXPECT_TRUE(ptr_with_no_names->check_uint());
      EXPECT_TRUE(ptr_with_no_names->check_int());
      EXPECT_TRUE(ptr_with_no_names->check_float());
      EXPECT_TRUE(ptr_with_no_names->check_double());
      EXPECT_TRUE(ptr_with_no_names->check_char());
      EXPECT_TRUE(ptr_with_no_names->check_bool());
      EXPECT_TRUE(ptr_with_no_names->check_str());
      EXPECT_TRUE(ptr_with_no_names->check_vect());
      EXPECT_TRUE(ptr_with_no_names->check_list());
      EXPECT_TRUE(ptr_with_no_names->check_set());
      EXPECT_TRUE(ptr_with_no_names->check_map());
    }
  }
}

TEST(SerializationTests, ProtobufSerializers) {
  {
    std::stringstream ss;
    {
      protobuf_oarchive output_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names(
          new obj_with_named_members());
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names(
          new obj_with_unnamed_members());

      EXPECT_NO_THROW(output_arc << obj_with_names);
      EXPECT_NO_THROW(output_arc << ptr_with_names);
      EXPECT_NO_THROW(output_arc << obj_with_no_names);
      EXPECT_NO_THROW(output_arc << ptr_with_no_names);
    }

    {
      protobuf_iarchive input_arc(ss);

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names;
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names;

      EXPECT_NO_THROW(input_arc >> obj_with_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_names);
      EXPECT_NO_THROW(input_arc >> obj_with_no_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_no_names);

      EXPECT_TRUE(obj_with_names.check_uint());
      EXPECT_TRUE(obj_with_names.check_int());
      EXPECT_TRUE(obj_with_names.check_float());
      EXPECT_TRUE(obj_with_names.check_double());
      EXPECT_TRUE(obj_with_names.check_char());
      EXPECT_TRUE(obj_with_names.check_bool());
      EXPECT_TRUE(obj_with_names.check_str());
      EXPECT_TRUE(obj_with_names.check_vect());
      EXPECT_TRUE(obj_with_names.check_list());
      EXPECT_TRUE(obj_with_names.check_set());
      EXPECT_TRUE(obj_with_names.check_map());

      EXPECT_TRUE(ptr_with_names);
      EXPECT_TRUE(ptr_with_names->check_uint());
      EXPECT_TRUE(ptr_with_names->check_int());
      EXPECT_TRUE(ptr_with_names->check_float());
      EXPECT_TRUE(ptr_with_names->check_double());
      EXPECT_TRUE(ptr_with_names->check_char());
      EXPECT_TRUE(ptr_with_names->check_bool());
      EXPECT_TRUE(ptr_with_names->check_str());
      EXPECT_TRUE(ptr_with_names->check_vect());
      EXPECT_TRUE(ptr_with_names->check_list());
      EXPECT_TRUE(ptr_with_names->check_set());
      EXPECT_TRUE(ptr_with_names->check_map());

      EXPECT_TRUE(obj_with_no_names.check_uint());
      EXPECT_TRUE(obj_with_no_names.check_int());
      EXPECT_TRUE(obj_with_no_names.check_float());
      EXPECT_TRUE(obj_with_no_names.check_double());
      EXPECT_TRUE(obj_with_no_names.check_char());
      EXPECT_TRUE(obj_with_no_names.check_bool());
      EXPECT_TRUE(obj_with_no_names.check_str());
      EXPECT_TRUE(obj_with_no_names.check_vect());
      EXPECT_TRUE(obj_with_no_names.check_list());
      EXPECT_TRUE(obj_with_no_names.check_set());
      EXPECT_TRUE(obj_with_no_names.check_map());

      EXPECT_TRUE(ptr_with_no_names);
      EXPECT_TRUE(ptr_with_no_names->check_uint());
      EXPECT_TRUE(ptr_with_no_names->check_int());
      EXPECT_TRUE(ptr_with_no_names->check_float());
      EXPECT_TRUE(ptr_with_no_names->check_double());
      EXPECT_TRUE(ptr_with_no_names->check_char());
      EXPECT_TRUE(ptr_with_no_names->check_bool());
      EXPECT_TRUE(ptr_with_no_names->check_str());
      EXPECT_TRUE(ptr_with_no_names->check_vect());
      EXPECT_TRUE(ptr_with_no_names->check_list());
      EXPECT_TRUE(ptr_with_no_names->check_set());
      EXPECT_TRUE(ptr_with_no_names->check_map());
    }
  }
}

TEST(SerializationTests, ObjTreeSerializers) {
  {
    objtree_editor otree;
    {
      objtree_oarchive output_arc(otree.get_object_graph(),
                                  otree.get_root_node());

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names(
          new obj_with_named_members());
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names(
          new obj_with_unnamed_members());

      EXPECT_NO_THROW(output_arc << obj_with_names);
      EXPECT_NO_THROW(output_arc << ptr_with_names);
      EXPECT_NO_THROW(output_arc << obj_with_no_names);
      EXPECT_NO_THROW(output_arc << ptr_with_no_names);
    }

    {
      objtree_iarchive input_arc(otree.get_object_graph(),
                                 otree.get_root_node());

      obj_with_named_members obj_with_names;
      std::shared_ptr<obj_with_named_members> ptr_with_names;
      obj_with_unnamed_members obj_with_no_names;
      std::shared_ptr<obj_with_unnamed_members> ptr_with_no_names;

      EXPECT_NO_THROW(input_arc >> obj_with_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_names);
      EXPECT_NO_THROW(input_arc >> obj_with_no_names);
      EXPECT_NO_THROW(input_arc >> ptr_with_no_names);

      EXPECT_TRUE(obj_with_names.check_uint());
      EXPECT_TRUE(obj_with_names.check_int());
      EXPECT_TRUE(obj_with_names.check_float());
      EXPECT_TRUE(obj_with_names.check_double());
      EXPECT_TRUE(obj_with_names.check_char());
      EXPECT_TRUE(obj_with_names.check_bool());
      EXPECT_TRUE(obj_with_names.check_str());
      EXPECT_TRUE(obj_with_names.check_vect());
      EXPECT_TRUE(obj_with_names.check_list());
      EXPECT_TRUE(obj_with_names.check_set());
      EXPECT_TRUE(obj_with_names.check_map());

      EXPECT_TRUE(ptr_with_names);
      EXPECT_TRUE(ptr_with_names->check_uint());
      EXPECT_TRUE(ptr_with_names->check_int());
      EXPECT_TRUE(ptr_with_names->check_float());
      EXPECT_TRUE(ptr_with_names->check_double());
      EXPECT_TRUE(ptr_with_names->check_char());
      EXPECT_TRUE(ptr_with_names->check_bool());
      EXPECT_TRUE(ptr_with_names->check_str());
      EXPECT_TRUE(ptr_with_names->check_vect());
      EXPECT_TRUE(ptr_with_names->check_list());
      EXPECT_TRUE(ptr_with_names->check_set());
      EXPECT_TRUE(ptr_with_names->check_map());

      EXPECT_TRUE(obj_with_no_names.check_uint());
      EXPECT_TRUE(obj_with_no_names.check_int());
      EXPECT_TRUE(obj_with_no_names.check_float());
      EXPECT_TRUE(obj_with_no_names.check_double());
      EXPECT_TRUE(obj_with_no_names.check_char());
      EXPECT_TRUE(obj_with_no_names.check_bool());
      EXPECT_TRUE(obj_with_no_names.check_str());
      EXPECT_TRUE(obj_with_no_names.check_vect());
      EXPECT_TRUE(obj_with_no_names.check_list());
      EXPECT_TRUE(obj_with_no_names.check_set());
      EXPECT_TRUE(obj_with_no_names.check_map());

      EXPECT_TRUE(ptr_with_no_names);
      EXPECT_TRUE(ptr_with_no_names->check_uint());
      EXPECT_TRUE(ptr_with_no_names->check_int());
      EXPECT_TRUE(ptr_with_no_names->check_float());
      EXPECT_TRUE(ptr_with_no_names->check_double());
      EXPECT_TRUE(ptr_with_no_names->check_char());
      EXPECT_TRUE(ptr_with_no_names->check_bool());
      EXPECT_TRUE(ptr_with_no_names->check_str());
      EXPECT_TRUE(ptr_with_no_names->check_vect());
      EXPECT_TRUE(ptr_with_no_names->check_list());
      EXPECT_TRUE(ptr_with_no_names->check_set());
      EXPECT_TRUE(ptr_with_no_names->check_map());
    }
  }
}

}  // namespace
}  // namespace ReaK::serialization
