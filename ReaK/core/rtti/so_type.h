/**
 * \file so_type.h
 *
 * This library defines a set of template classes used to create the type descriptors for
 * the types registered to the ReaK::rtti system.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
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

#ifndef REAK_CORE_RTTI_SO_TYPE_H_
#define REAK_CORE_RTTI_SO_TYPE_H_

#include "ReaK/core/base/cnst_string.h"

#include <memory>
#include <string>
#include <type_traits>

/** Main namespace for ReaK */
namespace ReaK {

class shared_object;  // forward-declaration.

/** Main namespace for ReaK's Serialization */
class serializable;  // forward-declaration

/** Main namespace for ReaK's Run-time Type Identification (RTTI) */
namespace rtti {

using construct_ptr = std::shared_ptr<shared_object> (*)();

// this is really the only thing that the class needs to define
template <typename T>
struct get_type_id {
  static constexpr unsigned int ID = T::rk_rtti_ID;
  static constexpr auto type_name = T::rk_rtti_TypeName;

  static construct_ptr CreatePtr() noexcept { return T::rk_rtti_CreatePtr(); }

  using save_type = const serializable&;
  using load_type = serializable&;
};

static constexpr auto lsl_left_bracket = std::string_view{"<"};
static constexpr auto lsl_right_bracket = std::string_view{">"};
static constexpr auto lsl_comma = std::string_view{","};

template <unsigned int U>
struct get_type_id<std::integral_constant<unsigned int, U>> {
  static constexpr unsigned int ID = U;
  static constexpr auto type_name = ct_itoa_v<ID>;

  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<std::true_type> {
  static constexpr unsigned int ID = 2;
  static constexpr auto type_name = std::string_view{"true"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<std::false_type> {
  static constexpr unsigned int ID = 1;
  static constexpr auto type_name = std::string_view{"false"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <int I>
struct get_type_id<std::integral_constant<int, I>> {
  static constexpr unsigned int ID = I;
  static constexpr auto type_name = ct_itoa_v<ID>;
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

struct null_type_id {
  static constexpr unsigned int ID = 0;
};

struct null_type_info {
  using type = null_type_id;
  static constexpr std::string_view type_name = {};
};

template <>
struct get_type_id<null_type_info> {
  static constexpr unsigned int ID = 0;
  static constexpr std::string_view type_name = {};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

namespace {

template <typename T, typename Tail = null_type_id>
struct type_id {
  using tail = Tail;
  static constexpr unsigned int ID = ::ReaK::rtti::get_type_id<T>::ID;
};

}  // namespace

template <typename Tail>
struct get_type_name_tail {
  static constexpr auto value = ct_concat_v<lsl_comma, Tail::type_name>;
};

template <>
struct get_type_name_tail<null_type_info> {
  static constexpr std::string_view value = {};
};

template <typename T, typename Tail = null_type_info>
struct get_type_info {
  using type = type_id<T, typename Tail::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<T>::type_name, get_type_name_tail<Tail>::value>;
};

namespace {

template <typename... T>
struct get_type_info_seq;

template <typename T1>
struct get_type_info_seq<T1> {
  template <typename Tail = null_type_info>
  struct with_tail {
    using type = get_type_info<T1, Tail>;
  };
  static constexpr auto type_name = get_type_id<T1>::type_name;
};

template <typename T1, typename... T>
struct get_type_info_seq<T1, T...> {
  template <typename Tail = null_type_info>
  struct with_tail {
    using type = get_type_info<
        T1, typename get_type_info_seq<T...>::template with_tail<Tail>::type>;
  };

  static constexpr auto type_name =
      ct_concat_v<get_type_id<T1>::type_name, lsl_comma,
                  get_type_info_seq<T...>::type_name>;
};

}  // namespace

template <template <typename...> class U, typename Tail, typename... T>
struct get_type_info<U<T...>, Tail> {
  using type = type_id<
      U<T...>,
      typename get_type_info_seq<T...>::template with_tail<Tail>::type::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<U<T...>>::type_name, lsl_left_bracket,
                  get_type_info_seq<T...>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};

template <unsigned int U, typename Tail>
struct get_type_info<std::integral_constant<unsigned int, U>, Tail> {
  using type =
      type_id<std::integral_constant<unsigned int, U>, typename Tail::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<std::integral_constant<unsigned int, U>>::type_name,
      get_type_name_tail<Tail>::value>;
};

template <typename Tail>
struct get_type_info<std::true_type, Tail> {
  using type = type_id<std::true_type, typename Tail::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::true_type>::type_name,
                  get_type_name_tail<Tail>::value>;
};

template <typename Tail>
struct get_type_info<std::false_type, Tail> {
  using type = type_id<std::false_type, typename Tail::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::false_type>::type_name,
                  get_type_name_tail<Tail>::value>;
};

template <int I, typename Tail>
struct get_type_info<std::integral_constant<int, I>, Tail> {
  using type = type_id<std::integral_constant<int, I>, typename Tail::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::integral_constant<int, I>>::type_name,
                  get_type_name_tail<Tail>::value>;
};

/**
 * This class is used to identify and register in the platform a new object
 * type. When a plugin is loaded, at initialization of the plugin, it is responsible
 * for registering all the new types it contains.
 */
class so_type {
 public:
  so_type(const so_type&) = delete;
  so_type& operator=(const so_type&) = delete;

  static so_type* createTypeInfo(unsigned int aTypeVersion,
                                 unsigned int* aTypeID,
                                 std::string_view aTypeName,
                                 construct_ptr aConstruct);

  static unsigned int* createTypeID(unsigned int aTypeIDSize);

 protected:
  so_type() = default;

  static bool compare_equal(const unsigned int* pid1, const unsigned int* pid2);

 public:
  /// This function adds a Descendant of this.
  so_type* addDescendant(so_type* aObj);

  so_type* addAncestor(so_type* aObj);

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type* findDescendant(const unsigned int* aTypeID);

  /// This function gets the number of direct descendants of this.
  unsigned int getDirectDescendantCount();

  /// This function gets a type record by index in the direct descendants of this.
  so_type* getDirectDescendant(unsigned int aIndex);

  /// This function checks if a typeID is parent to this.
  so_type* findAncestor(const unsigned int* aTypeID);

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type* findDescendant(so_type* aTypeID);

  /// This function checks if a typeID is parent to this.
  so_type* findAncestor(so_type* aTypeID);

  /// This function inserts this into a global repo.
  void insertToRepo(so_type* aRepo);

  [[nodiscard]] const unsigned int* TypeID_begin() const;

  [[nodiscard]] unsigned int TypeVersion() const;

  [[nodiscard]] const std::string& TypeName() const;

  [[nodiscard]] std::shared_ptr<shared_object> CreateObject() const;

  [[nodiscard]] bool isConcrete() const;

  friend bool operator==(const so_type& t1, const so_type& t2) {
    return compare_equal(t1.TypeID_begin(), t2.TypeID_begin());
  }

  friend bool operator!=(const so_type& t1, const so_type& t2) {
    return !compare_equal(t1.TypeID_begin(), t2.TypeID_begin());
  }
};

struct so_type_ptr {
  so_type* ptr;
  explicit so_type_ptr(so_type* aPtr) : ptr(aPtr){};
  ~so_type_ptr();
};

namespace {

template <typename T>
struct get_type_id_prop;

template <>
struct get_type_id_prop<null_type_id> {
  static constexpr unsigned int count = 1;
  static unsigned int at(unsigned int /*unused*/) { return 0; }
};

template <typename T>
struct get_type_id_prop {
  static constexpr unsigned int count =
      get_type_id_prop<typename T::tail>::count + 1;
  static unsigned int at(unsigned int i) {
    if (i == 0) {
      return T::ID;
    }
    return get_type_id_prop<typename T::tail>::at(--i);
  }
};

template <typename T>
so_type_ptr create_type_descriptor(unsigned int aVersion = 1) {
  using SizeType = unsigned int;
  const SizeType TypeIDLength =
      get_type_id_prop<typename get_type_info<T>::type>::count;
  SizeType* typeID = so_type::createTypeID(TypeIDLength);
  for (SizeType i = 0; i < TypeIDLength; ++i) {
    typeID[i] = get_type_id_prop<typename get_type_info<T>::type>::at(i);
  }
  constexpr auto tname = get_type_info<T>::type_name;
  return so_type_ptr(so_type::createTypeInfo(aVersion, typeID, tname,
                                             get_type_id<T>::CreatePtr()));
}

}  // namespace

so_type_ptr create_dummy_so_type(const unsigned int* aTypeID);

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_CORE_RTTI_SO_TYPE_H_
