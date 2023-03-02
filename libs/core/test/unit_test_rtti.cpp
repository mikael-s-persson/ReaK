
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include <ReaK/core/base/shared_object.hpp>
#include <ReaK/core/rtti/rtti.hpp>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE rtti
#include <boost/test/unit_test.hpp>

class dummy_rtti_test : public ReaK::shared_object {
 public:
  int value;

  dummy_rtti_test() : ReaK::shared_object(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    ReaK::shared_object::save(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    ReaK::shared_object::load(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(dummy_rtti_test, 0xFFFFFFF0, 1, "dummy_rtti_test",
                              ReaK::shared_object)
};

template <int N, typename T>
class dummy_rtti_template : public ReaK::shared_object {
 public:
  typedef dummy_rtti_template<N, T> self;

  dummy_rtti_template() : ReaK::shared_object(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    ReaK::shared_object::save(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    ReaK::shared_object::load(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xFFFFFFF1, 1, "dummy_rtti_template",
                              ReaK::shared_object)
};

// NOTE: Because the template has a non-type template parameter, we need to specialize the get_type_info template:
namespace ReaK {
namespace rtti {

template <int N, typename T, typename Tail>
struct get_type_info<dummy_rtti_template<N, T>, Tail> {
  typedef type_id<dummy_rtti_template<N, T>,
                  typename get_type_info<std::integral_constant<int, N>,
                                         get_type_info<T, Tail>>::type>
      type;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<dummy_rtti_template<N, T>>::type_name,
                  lsl_left_bracket,
                  get_type_id<std::integral_constant<int, N>>::type_name,
                  lsl_comma, get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
};  // namespace rtti
};  // namespace ReaK

template <typename N, typename T>
class dummy_rtti_template2 : public ReaK::shared_object {
 public:
  typedef dummy_rtti_template2<N, T> self;

  N n_value;
  T t_value;

  dummy_rtti_template2() : ReaK::shared_object(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    ReaK::shared_object::save(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    ReaK::shared_object::load(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xFFFFFFF2, 1, "dummy_rtti_template2",
                              ReaK::shared_object)
};

class dummy_rtti_base1 : public virtual ReaK::shared_object {
 public:
  int b1_value;

  dummy_rtti_base1(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    ReaK::shared_object::save(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(b1_value);
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    ReaK::shared_object::load(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(b1_value);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(dummy_rtti_base1, 0xFFFFFFF3, 1,
                              "dummy_rtti_base1", ReaK::shared_object)
};

class dummy_rtti_base2 : public virtual ReaK::shared_object {
 public:
  int b2_value;

  dummy_rtti_base2(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    ReaK::shared_object::save(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(b2_value);
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    ReaK::shared_object::load(
        A, ReaK::shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(b2_value);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(dummy_rtti_base2, 0xFFFFFFF4, 1,
                              "dummy_rtti_base2", ReaK::shared_object)
};

class dummy_rtti_multi_inherit : public dummy_rtti_base1,
                                 public dummy_rtti_base2 {
 public:
  dummy_rtti_multi_inherit() : dummy_rtti_base1(), dummy_rtti_base2(){};

  void save(ReaK::serialization::oarchive& A, unsigned int) const override {
    dummy_rtti_base1::save(
        A, dummy_rtti_base1::getStaticObjectType()->TypeVersion());
    dummy_rtti_base2::save(
        A, dummy_rtti_base2::getStaticObjectType()->TypeVersion());
  };
  void load(ReaK::serialization::iarchive& A, unsigned int) override {
    dummy_rtti_base1::load(
        A, dummy_rtti_base1::getStaticObjectType()->TypeVersion());
    dummy_rtti_base2::load(
        A, dummy_rtti_base2::getStaticObjectType()->TypeVersion());
  };

  RK_RTTI_MAKE_CONCRETE_2BASE(dummy_rtti_multi_inherit, 0xFFFFFFF5, 1,
                              "dummy_rtti_multi_inherit", dummy_rtti_base1,
                              dummy_rtti_base2)
};

BOOST_AUTO_TEST_CASE(rtti_test) {
  using namespace ReaK;

  {
    rtti::so_type* p_tp = dummy_rtti_test::getStaticObjectType();
    BOOST_REQUIRE(p_tp != nullptr);
    BOOST_CHECK_EQUAL(p_tp->TypeName(), "dummy_rtti_test");
    BOOST_CHECK_EQUAL(p_tp->TypeVersion(), 1);
    BOOST_CHECK(p_tp->isConcrete());
    const unsigned int* t_id = p_tp->TypeID_begin();
    BOOST_CHECK_EQUAL(*t_id, 0xFFFFFFF0);
    BOOST_CHECK_EQUAL(*(t_id + 1), 0);

    rtti::so_type* p_so_tp =
        p_tp->findAncestor(shared_object::getStaticObjectType());
    BOOST_REQUIRE(p_so_tp);

    rtti::so_type* p_so_tp_2 = p_tp->findAncestor(
        shared_object::getStaticObjectType()->TypeID_begin());
    BOOST_CHECK_EQUAL(p_so_tp_2, p_so_tp);

    rtti::so_type* p_tp_2 = p_so_tp->findDescendant(p_tp);
    BOOST_CHECK_EQUAL(p_tp_2, p_tp);

    rtti::so_type* p_tp_3 = p_so_tp->findDescendant(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_3, p_tp);

    rtti::so_type* p_tp_4 =
        rtti::so_type_repo::getInstance().findType(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_4, p_tp);

    rtti::so_type* p_tp_5 = rtti::so_type_repo::getInstance().findType(p_tp);
    BOOST_CHECK_EQUAL(p_tp_5, p_tp);

    unsigned int desc_count = p_so_tp->getDirectDescendantCount();
    bool found_descendant_in_so_type = false;
    for (unsigned int i = 0; i < desc_count; ++i) {
      rtti::so_type* tmp_tp = p_so_tp->getDirectDescendant(i);
      if (tmp_tp == p_tp)
        found_descendant_in_so_type = true;
    };
    BOOST_CHECK(found_descendant_in_so_type);

    std::shared_ptr<shared_object> p_obj = p_tp->CreateObject();
    BOOST_REQUIRE(p_obj);

    void* p_obj_raw = p_obj->castTo(p_tp);
    BOOST_CHECK(p_obj_raw);

    std::shared_ptr<dummy_rtti_test> p_obj_statcast =
        rtti::rk_static_ptr_cast<dummy_rtti_test>(p_obj);
    BOOST_CHECK(p_obj_statcast);
    BOOST_CHECK_EQUAL(p_obj_statcast.get(), p_obj_raw);

    std::shared_ptr<dummy_rtti_test> p_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_test>(p_obj);
    BOOST_CHECK(p_obj_dyncast);
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_statcast.get());
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_raw);
  };

  {
    rtti::so_type* p_tp =
        dummy_rtti_template<10, dummy_rtti_test>::getStaticObjectType();
    BOOST_REQUIRE(p_tp != nullptr);
    BOOST_CHECK_EQUAL(p_tp->TypeName(),
                      "dummy_rtti_template<10,dummy_rtti_test>");
    BOOST_CHECK_EQUAL(p_tp->TypeVersion(), 1);
    BOOST_CHECK(p_tp->isConcrete());
    const unsigned int* t_id = p_tp->TypeID_begin();
    BOOST_CHECK_EQUAL(*t_id, 0xFFFFFFF1);
    BOOST_CHECK_EQUAL(*(t_id + 1), 10);
    BOOST_CHECK_EQUAL(*(t_id + 2), 0xFFFFFFF0);
    BOOST_CHECK_EQUAL(*(t_id + 3), 0);

    rtti::so_type* p_so_tp =
        p_tp->findAncestor(shared_object::getStaticObjectType());
    BOOST_CHECK(p_so_tp);

    rtti::so_type* p_so_tp_2 = p_tp->findAncestor(
        shared_object::getStaticObjectType()->TypeID_begin());
    BOOST_CHECK_EQUAL(p_so_tp_2, p_so_tp);

    rtti::so_type* p_tp_2 = p_so_tp->findDescendant(p_tp);
    BOOST_CHECK_EQUAL(p_tp_2, p_tp);

    rtti::so_type* p_tp_3 = p_so_tp->findDescendant(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_3, p_tp);

    rtti::so_type* p_tp_4 =
        rtti::so_type_repo::getInstance().findType(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_4, p_tp);

    rtti::so_type* p_tp_5 = rtti::so_type_repo::getInstance().findType(p_tp);
    BOOST_CHECK_EQUAL(p_tp_5, p_tp);

    unsigned int desc_count = p_so_tp->getDirectDescendantCount();
    bool found_descendant_in_so_type = false;
    for (unsigned int i = 0; i < desc_count; ++i) {
      rtti::so_type* tmp_tp = p_so_tp->getDirectDescendant(i);
      if (tmp_tp == p_tp)
        found_descendant_in_so_type = true;
    };
    BOOST_CHECK(found_descendant_in_so_type);

    std::shared_ptr<shared_object> p_obj = p_tp->CreateObject();
    BOOST_REQUIRE(p_obj);

    void* p_obj_raw = p_obj->castTo(p_tp);
    BOOST_CHECK(p_obj_raw);

    std::shared_ptr<dummy_rtti_template<10, dummy_rtti_test>> p_obj_statcast =
        rtti::rk_static_ptr_cast<dummy_rtti_template<10, dummy_rtti_test>>(
            p_obj);
    BOOST_CHECK(p_obj_statcast);
    BOOST_CHECK_EQUAL(p_obj_statcast.get(), p_obj_raw);

    std::shared_ptr<dummy_rtti_template<10, dummy_rtti_test>> p_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_template<10, dummy_rtti_test>>(
            p_obj);
    BOOST_CHECK(p_obj_dyncast);
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_statcast.get());
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_raw);
  };

  {
    rtti::so_type* p_tp =
        dummy_rtti_template2<int, dummy_rtti_test>::getStaticObjectType();
    BOOST_REQUIRE(p_tp != nullptr);
    BOOST_CHECK_EQUAL(p_tp->TypeName(),
                      "dummy_rtti_template2<int,dummy_rtti_test>");
    BOOST_CHECK_EQUAL(p_tp->TypeVersion(), 1);
    BOOST_CHECK(p_tp->isConcrete());
    const unsigned int* t_id = p_tp->TypeID_begin();
    BOOST_CHECK_EQUAL(*t_id, 0xFFFFFFF2);
    BOOST_CHECK_EQUAL(*(t_id + 1), 0x00000001);
    BOOST_CHECK_EQUAL(*(t_id + 2), 0xFFFFFFF0);
    BOOST_CHECK_EQUAL(*(t_id + 3), 0);

    rtti::so_type* p_so_tp =
        p_tp->findAncestor(shared_object::getStaticObjectType());
    BOOST_CHECK(p_so_tp);

    rtti::so_type* p_so_tp_2 = p_tp->findAncestor(
        shared_object::getStaticObjectType()->TypeID_begin());
    BOOST_CHECK_EQUAL(p_so_tp_2, p_so_tp);

    rtti::so_type* p_tp_2 = p_so_tp->findDescendant(p_tp);
    BOOST_CHECK_EQUAL(p_tp_2, p_tp);

    rtti::so_type* p_tp_3 = p_so_tp->findDescendant(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_3, p_tp);

    rtti::so_type* p_tp_4 =
        rtti::so_type_repo::getInstance().findType(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_4, p_tp);

    rtti::so_type* p_tp_5 = rtti::so_type_repo::getInstance().findType(p_tp);
    BOOST_CHECK_EQUAL(p_tp_5, p_tp);

    unsigned int desc_count = p_so_tp->getDirectDescendantCount();
    bool found_descendant_in_so_type = false;
    for (unsigned int i = 0; i < desc_count; ++i) {
      rtti::so_type* tmp_tp = p_so_tp->getDirectDescendant(i);
      if (tmp_tp == p_tp)
        found_descendant_in_so_type = true;
    };
    BOOST_CHECK(found_descendant_in_so_type);

    std::shared_ptr<shared_object> p_obj = p_tp->CreateObject();
    BOOST_REQUIRE(p_obj);

    void* p_obj_raw = p_obj->castTo(p_tp);
    BOOST_CHECK(p_obj_raw);

    std::shared_ptr<dummy_rtti_template2<int, dummy_rtti_test>> p_obj_statcast =
        rtti::rk_static_ptr_cast<dummy_rtti_template2<int, dummy_rtti_test>>(
            p_obj);
    BOOST_CHECK(p_obj_statcast);
    BOOST_CHECK_EQUAL(p_obj_statcast.get(), p_obj_raw);

    std::shared_ptr<dummy_rtti_template2<int, dummy_rtti_test>> p_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_template2<int, dummy_rtti_test>>(
            p_obj);
    BOOST_CHECK(p_obj_dyncast);
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_statcast.get());
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_raw);
  };

  {
    rtti::so_type* p_tp = dummy_rtti_multi_inherit::getStaticObjectType();
    BOOST_REQUIRE(p_tp != nullptr);
    BOOST_CHECK_EQUAL(p_tp->TypeName(), "dummy_rtti_multi_inherit");
    BOOST_CHECK_EQUAL(p_tp->TypeVersion(), 1);
    BOOST_CHECK(p_tp->isConcrete());
    const unsigned int* t_id = p_tp->TypeID_begin();
    BOOST_CHECK_EQUAL(*t_id, 0xFFFFFFF5);
    BOOST_CHECK_EQUAL(*(t_id + 1), 0);

    rtti::so_type* p_so_tp =
        p_tp->findAncestor(shared_object::getStaticObjectType());
    BOOST_CHECK(p_so_tp);

    rtti::so_type* p_b1_tp =
        p_tp->findAncestor(dummy_rtti_base1::getStaticObjectType());
    BOOST_CHECK(p_b1_tp);

    rtti::so_type* p_b2_tp =
        p_tp->findAncestor(dummy_rtti_base2::getStaticObjectType());
    BOOST_CHECK(p_b2_tp);

    rtti::so_type* p_so_tp_2 = p_tp->findAncestor(
        shared_object::getStaticObjectType()->TypeID_begin());
    BOOST_CHECK_EQUAL(p_so_tp_2, p_so_tp);

    rtti::so_type* p_tp_2 = p_so_tp->findDescendant(p_tp);
    BOOST_CHECK_EQUAL(p_tp_2, p_tp);

    rtti::so_type* p_tp_3 = p_so_tp->findDescendant(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_3, p_tp);

    rtti::so_type* p_tp_4 =
        rtti::so_type_repo::getInstance().findType(p_tp->TypeID_begin());
    BOOST_CHECK_EQUAL(p_tp_4, p_tp);

    rtti::so_type* p_tp_5 = rtti::so_type_repo::getInstance().findType(p_tp);
    BOOST_CHECK_EQUAL(p_tp_5, p_tp);

    unsigned int b1_desc_count = p_b1_tp->getDirectDescendantCount();
    bool found_descendant_in_b1_so_type = false;
    for (unsigned int i = 0; i < b1_desc_count; ++i) {
      rtti::so_type* tmp_tp = p_b1_tp->getDirectDescendant(i);
      if (tmp_tp == p_tp)
        found_descendant_in_b1_so_type = true;
    };
    BOOST_CHECK(found_descendant_in_b1_so_type);

    unsigned int b2_desc_count = p_b2_tp->getDirectDescendantCount();
    bool found_descendant_in_b2_so_type = false;
    for (unsigned int i = 0; i < b2_desc_count; ++i) {
      rtti::so_type* tmp_tp = p_b2_tp->getDirectDescendant(i);
      if (tmp_tp == p_tp)
        found_descendant_in_b2_so_type = true;
    };
    BOOST_CHECK(found_descendant_in_b2_so_type);

    std::shared_ptr<shared_object> p_obj = p_tp->CreateObject();
    BOOST_REQUIRE(p_obj);

    void* p_obj_raw = p_obj->castTo(p_tp);
    BOOST_CHECK(p_obj_raw);

    std::shared_ptr<dummy_rtti_multi_inherit> p_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_multi_inherit>(p_obj);
    BOOST_CHECK(p_obj_dyncast);
    BOOST_CHECK_EQUAL(p_obj_dyncast.get(), p_obj_raw);

    p_obj_dyncast->b1_value = 42;
    p_obj_dyncast->b2_value = 69;

    void* p_b1_obj_raw = p_obj->castTo(p_b1_tp);
    BOOST_CHECK(p_b1_obj_raw);

    std::shared_ptr<dummy_rtti_base1> p_b1_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_base1>(p_obj);
    BOOST_CHECK(p_b1_obj_dyncast);
    BOOST_CHECK_EQUAL(p_b1_obj_dyncast.get(), p_b1_obj_raw);
    BOOST_CHECK_EQUAL(p_obj_dyncast->b1_value, p_b1_obj_dyncast->b1_value);

    void* p_b2_obj_raw = p_obj->castTo(p_b2_tp);
    BOOST_CHECK(p_b2_obj_raw);

    std::shared_ptr<dummy_rtti_base2> p_b2_obj_dyncast =
        rtti::rk_dynamic_ptr_cast<dummy_rtti_base2>(p_obj);
    BOOST_CHECK(p_b2_obj_dyncast);
    BOOST_CHECK_EQUAL(p_b2_obj_dyncast.get(), p_b2_obj_raw);
    BOOST_CHECK_EQUAL(p_obj_dyncast->b2_value, p_b2_obj_dyncast->b2_value);
  };
};
