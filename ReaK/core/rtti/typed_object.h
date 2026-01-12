/**
 * \file typed_object.h
 *
 * This library declares the basic type for all objects in ReaK which are registered to the
 * ReaK::rtti (note that special template techniques can be used to register other types too,
 * without the need to derive from typed_object or contain virtual member functions, see "so_type.hpp").
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date april 2011 (orginally january 2010)
 */

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

#ifndef REAK_CORE_RTTI_TYPED_OBJECT_H_
#define REAK_CORE_RTTI_TYPED_OBJECT_H_

#include "ReaK/core/rtti/rtti.h"

#include <memory>

namespace ReaK {

/**
 * The basic class "ReaK::typed_object" allows all descendants
 * to be associated with at shared object type (so_type) structure. All descendant classes can
 * be registered in the ReaK::rtti system and thus, enjoy the dynamic casting
 * features of the ReaK platform, across executable modules. This feature is also required for
 * serializability of the classes as the ReaK::rtti system is used to identity the types of the
 * serialized objects.
 */
class typed_object {
 public:
  /**
   * This method is used to perform up- and down- casting of object pointers via a virtual call.
   */
  virtual void* cast_to(rtti::so_type* tp) {
    if (*(tp->id_begin()) == 0) {
      return static_cast<void*>(this);
    }
    return nullptr;
  }

  /**
   * This method is used to perform up- and down- casting of const-object pointers via a virtual call.
   */
  virtual const void* cast_to(rtti::so_type* tp) const {
    if (*(tp->id_begin()) == 0) {
      return static_cast<const void*>(this);
    }
    return nullptr;
  }

  virtual ~typed_object() = default;

  /** This method fetches the object type structure from the ReaK::rtti system or creates it if it has not been
   * registered yet. */
  [[nodiscard]] virtual rtti::so_type* get_object_type() const {
    return nullptr;
  }
  /** This method fetches the object type structure from the ReaK::rtti system or creates it if it has not been
   * registered yet. */
  [[nodiscard]] static rtti::so_type* get_static_object_type() {
    return nullptr;
  }
};

namespace rtti {

template <typename T>
bool rk_is_of_type(const typed_object& obj) {
  return (obj.cast_to(T::get_static_object_type()) != nullptr);
}

template <typename T>
bool rk_is_of_type(const typed_object* obj) {
  return (obj->cast_to(T::get_static_object_type()) != nullptr);
}

template <typename T>
bool rk_is_of_type(const std::shared_ptr<const typed_object>& obj) {
  return (obj->cast_to(T::get_static_object_type()) != nullptr);
}

/**
 * This function replaces the standard C++ static cast for pointers (i.e. static_cast<>()) and furthermore, also
 * replaces the std::shared_ptr static cast (i.e. std::static_pointer_cast<>()). This new function
 * for dynamic casting is required in order to overload between C++0x smart-pointers and Boost smart-pointers.
 */
template <typename Y, typename U>
std::shared_ptr<Y> rk_static_ptr_cast(const std::shared_ptr<U>& p) {
  return std::static_pointer_cast<Y>(p);
}

/**
 * This function replaces the standard C++ dynamic cast (i.e. dynamic_cast<>()) and furthermore, also
 * replaces the std::shared_ptr dynamic cast (i.e. std::dynamic_pointer_cast<>()). This new function
 * for dynamic casting is required in order for ReaK::rtti system to take precedence over the C++ RTTI
 * because, unlike the C++ standard version of RTTI, this implementation will work across executable modules,
 * and thus, allow objects to be shared between modules with full dynamic up- and down- casting capabilities.
 */
template <typename Y, typename U>
std::shared_ptr<Y> rk_dynamic_ptr_cast(const std::shared_ptr<U>& p) {
  if (!p) {
    return std::shared_ptr<Y>();
  }
  return std::shared_ptr<Y>(
      p, static_cast<Y*>(p->cast_to(Y::get_static_object_type())));
}

template <typename Y, typename U, typename Deleter>
std::unique_ptr<Y, Deleter> rk_dynamic_ptr_cast(
    std::unique_ptr<U, Deleter>&& p) {
  if (!p) {
    return std::unique_ptr<Y, Deleter>();
  }
  void* tmp = p->cast_to(Y::get_static_object_type());
  if (tmp) {
    std::unique_ptr<Y, Deleter> r(tmp, p.get_deleter());
    p.release();
    return std::move(r);
  }
  return std::unique_ptr<Y, Deleter>();
}

template <typename T, typename Deleter>
bool rk_is_of_type(const std::unique_ptr<const typed_object, Deleter>& obj) {
  return (obj->cast_to(T::get_static_object_type()) != nullptr);
}

/// This MACRO creates a static (no-parameter) factory function for the current class CLASS_NAME.
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)            \
  static std::shared_ptr<::ReaK::shared_object> Create() {  \
    return std::shared_ptr<CLASS_NAME>(new CLASS_NAME());   \
  }                                                         \
  static ::ReaK::rtti::construct_ptr rk_rtti_create_ptr() { \
    return &Create;                                         \
  }

/// This MACRO registers a custom (no-parameter) factory function pointer for the current class.
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CUSTOM_FACTORY(CLASS_FACTORY)      \
  static ::ReaK::rtti::construct_ptr rk_rtti_create_ptr() { \
    return CLASS_FACTORY;                                   \
  }

/// This MACRO creates the static elements for the current class that registers the CLASS_ID and CLASS_NAME.
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_ID(CLASS_NAME, CLASS_ID)     \
  static constexpr std::uint32_t rk_rtti_type_id = CLASS_ID; \
  static constexpr auto rk_rtti_type_name = std::string_view{CLASS_NAME};

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is
/// guaranteed to be added if the class is instantiated).
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME, CLASS_VERSION)            \
  static ::ReaK::rtti::so_type* get_static_object_type() {                 \
    static ::ReaK::rtti::so_type_ptr ptr(                                  \
        ::ReaK::rtti::so_type_details::create_type_descriptor<CLASS_NAME>( \
            CLASS_VERSION));                                               \
    return ptr.ptr;                                                        \
  }                                                                        \
  ::ReaK::rtti::so_type* get_object_type() const override {                \
    (void)::ReaK::rtti::register_type<                                     \
        CLASS_NAME>::impl; /* force instantiate */                         \
    return CLASS_NAME::get_static_object_type();                           \
  }                                                                        \
  void* cast_to(::ReaK::rtti::so_type* tp) override {                      \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                  \
      return static_cast<void*>(this);                                     \
    }                                                                      \
    return nullptr;                                                        \
  }                                                                        \
  const void* cast_to(::ReaK::rtti::so_type* tp) const override {          \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                  \
      return static_cast<const void*>(this);                               \
    }                                                                      \
    return nullptr;                                                        \
  }

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is
/// guaranteed to be added if the class is instantiated).
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME) \
  static ::ReaK::rtti::so_type* get_static_object_type() {                 \
    static ::ReaK::rtti::so_type_ptr ptr(                                  \
        ::ReaK::rtti::so_type_details::create_type_descriptor<CLASS_NAME>( \
            CLASS_VERSION));                                               \
    static bool first_pass = true;                                         \
    if (first_pass) {                                                      \
      ptr.ptr->add_ancestor(BASE_NAME::get_static_object_type());          \
    }                                                                      \
    return ptr.ptr;                                                        \
  }                                                                        \
  ::ReaK::rtti::so_type* get_object_type() const override {                \
    (void)::ReaK::rtti::register_type<                                     \
        CLASS_NAME>::impl; /* force instantiate */                         \
    return CLASS_NAME::get_static_object_type();                           \
  }                                                                        \
  void* cast_to(::ReaK::rtti::so_type* tp) override {                      \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                  \
      return static_cast<void*>(this);                                     \
    }                                                                      \
    return BASE_NAME::cast_to(tp);                                         \
  }                                                                        \
  const void* cast_to(::ReaK::rtti::so_type* tp) const override {          \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                  \
      return static_cast<const void*>(this);                               \
    }                                                                      \
    return BASE_NAME::cast_to(tp);                                         \
  }

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is
/// guaranteed to be added if the class is instantiated).
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1, \
                                     BASE_NAME2)                            \
  static ::ReaK::rtti::so_type* get_static_object_type() {                  \
    static ::ReaK::rtti::so_type_ptr ptr(                                   \
        ::ReaK::rtti::so_type_details::create_type_descriptor<CLASS_NAME>(  \
            CLASS_VERSION));                                                \
    static bool first_pass = true;                                          \
    if (first_pass) {                                                       \
      ptr.ptr->add_ancestor(BASE_NAME1::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME2::get_static_object_type());          \
    }                                                                       \
    return ptr.ptr;                                                         \
  }                                                                         \
  ::ReaK::rtti::so_type* get_object_type() const override {                 \
    (void)::ReaK::rtti::register_type<                                      \
        CLASS_NAME>::impl; /* force instantiate */                          \
    return CLASS_NAME::get_static_object_type();                            \
  }                                                                         \
  void* cast_to(::ReaK::rtti::so_type* tp) override {                       \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return static_cast<void*>(this);                                      \
    }                                                                       \
    void* result = BASE_NAME1::cast_to(tp);                                 \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME2::cast_to(tp);                                         \
  }                                                                         \
  const void* cast_to(::ReaK::rtti::so_type* tp) const override {           \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return static_cast<const void*>(this);                                \
    }                                                                       \
    const void* result = BASE_NAME1::cast_to(tp);                           \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME2::cast_to(tp);                                         \
  }

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is
/// guaranteed to be added if the class is instantiated).
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1, \
                                     BASE_NAME2, BASE_NAME3)                \
  static ::ReaK::rtti::so_type* get_static_object_type() {                  \
    static ::ReaK::rtti::so_type_ptr ptr(                                   \
        ::ReaK::rtti::so_type_details::create_type_descriptor<CLASS_NAME>(  \
            CLASS_VERSION));                                                \
    static bool first_pass = true;                                          \
    if (first_pass) {                                                       \
      ptr.ptr->add_ancestor(BASE_NAME1::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME2::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME3::get_static_object_type());          \
    }                                                                       \
    return ptr.ptr;                                                         \
  }                                                                         \
  ::ReaK::rtti::so_type* get_object_type() const override {                 \
    (void)::ReaK::rtti::register_type<                                      \
        CLASS_NAME>::impl; /* force instantiate */                          \
    return CLASS_NAME::get_static_object_type();                            \
  }                                                                         \
  void* cast_to(::ReaK::rtti::so_type* tp) override {                       \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return static_cast<void*>(this);                                      \
    }                                                                       \
    void* result = BASE_NAME1::cast_to(tp);                                 \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME2::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME3::cast_to(tp);                                         \
  }                                                                         \
  const void* cast_to(::ReaK::rtti::so_type* tp) const override {           \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return static_cast<const void*>(this);                                \
    }                                                                       \
    const void* result = BASE_NAME1::cast_to(tp);                           \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME2::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME3::cast_to(tp);                                         \
  }

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is
/// guaranteed to be added if the class is instantiated).
// NOLINTNEXTLINE
#define RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1, \
                                     BASE_NAME2, BASE_NAME3, BASE_NAME4)    \
  static ::ReaK::rtti::so_type* get_static_object_type() {                  \
    static ::ReaK::rtti::so_type_ptr ptr(                                   \
        ::ReaK::rtti::so_type_details::create_type_descriptor<CLASS_NAME>(  \
            CLASS_VERSION));                                                \
    static bool first_pass = true;                                          \
    if (first_pass) {                                                       \
      ptr.ptr->add_ancestor(BASE_NAME1::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME2::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME3::get_static_object_type());          \
      ptr.ptr->add_ancestor(BASE_NAME4::get_static_object_type());          \
    }                                                                       \
    return ptr.ptr;                                                         \
  }                                                                         \
  ::ReaK::rtti::so_type* get_object_type() const override {                 \
    (void)::ReaK::rtti::register_type<                                      \
        CLASS_NAME>::impl; /* force instantiate */                          \
    return CLASS_NAME::get_static_object_type();                            \
  }                                                                         \
  void* cast_to(::ReaK::rtti::so_type* tp) override {                       \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return reinterpret_cast<void*>(this);                                 \
    }                                                                       \
    void* result = BASE_NAME1::cast_to(tp);                                 \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME2::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME3::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME4::cast_to(tp);                                         \
  }                                                                         \
  const void* cast_to(::ReaK::rtti::so_type* tp) const override {           \
    if (*tp == *(CLASS_NAME::get_static_object_type())) {                   \
      return reinterpret_cast<const void*>(this);                           \
    }                                                                       \
    const void* result = BASE_NAME1::cast_to(tp);                           \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME2::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    result = BASE_NAME3::cast_to(tp);                                       \
    if (result) {                                                           \
      return result;                                                        \
    }                                                                       \
    return BASE_NAME4::cast_to(tp);                                         \
  }

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class
 * declaration section (public), and all classes derived from typed_object are required to have this macro (or another
 * version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_CONCRETE_0BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION, \
                                    CLASS_STR_NAME)                      \
  RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)                               \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                    \
  RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME, CLASS_VERSION)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the
 * class declaration section (public), and all classes derived from typed_object are required to have this macro (or
 * another version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_ABSTRACT_0BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION, \
                                    CLASS_STR_NAME)                      \
  RK_RTTI_REGISTER_CUSTOM_FACTORY(0)                                     \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                    \
  RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME, CLASS_VERSION)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class
 * declaration section (public), and all classes derived from typed_object are required to have this macro (or another
 * version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_CONCRETE_1BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION, \
                                    CLASS_STR_NAME, BASE_NAME)           \
  RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)                               \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                    \
  RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the
 * class declaration section (public), and all classes derived from typed_object are required to have this macro (or
 * another version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_ABSTRACT_1BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION, \
                                    CLASS_STR_NAME, BASE_NAME)           \
  RK_RTTI_REGISTER_CUSTOM_FACTORY(0)                                     \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                    \
  RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class
 * declaration section (public), and all classes derived from typed_object are required to have this macro (or another
 * version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_CONCRETE_2BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2) \
  RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)                                  \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the
 * class declaration section (public), and all classes derived from typed_object are required to have this macro (or
 * another version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_ABSTRACT_2BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2) \
  RK_RTTI_REGISTER_CUSTOM_FACTORY(0)                                        \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class
 * declaration section (public), and all classes derived from typed_object are required to have this macro (or another
 * version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_CONCRETE_3BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2, \
                                    BASE_NAME3)                             \
  RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)                                  \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2, BASE_NAME3)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the
 * class declaration section (public), and all classes derived from typed_object are required to have this macro (or
 * another version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_ABSTRACT_3BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2, \
                                    BASE_NAME3)                             \
  RK_RTTI_REGISTER_CUSTOM_FACTORY(0)                                        \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2, BASE_NAME3)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class
 * declaration section (public), and all classes derived from typed_object are required to have this macro (or another
 * version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_CONCRETE_4BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2, \
                                    BASE_NAME3, BASE_NAME4)                 \
  RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME)                                  \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2, BASE_NAME3, BASE_NAME4)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the
 * class declaration section (public), and all classes derived from typed_object are required to have this macro (or
 * another version of it).*/
// NOLINTNEXTLINE
#define RK_RTTI_MAKE_ABSTRACT_4BASE(CLASS_NAME, CLASS_ID, CLASS_VERSION,    \
                                    CLASS_STR_NAME, BASE_NAME1, BASE_NAME2, \
                                    BASE_NAME3, BASE_NAME4)                 \
  RK_RTTI_REGISTER_CUSTOM_FACTORY(0)                                        \
  RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME, CLASS_ID)                       \
  RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME, CLASS_VERSION, BASE_NAME1,       \
                               BASE_NAME2, BASE_NAME3, BASE_NAME4)

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_CORE_RTTI_TYPED_OBJECT_H_
