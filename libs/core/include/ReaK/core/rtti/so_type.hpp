/**
 * \file so_type.hpp
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

#ifndef REAK_SO_TYPE_HPP
#define REAK_SO_TYPE_HPP

#include <ReaK/core/base/defs.hpp>

// This would be nice but the compilation times / mem are insane when type names are computed as constexpr strings.
// #define RK_RTTI_USE_CONSTEXPR_STRINGS

#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
// #include <ReaK/core/base/literal_string_list.hpp>
#include <ReaK/core/base/cnst_string.hpp>
#endif

#include <string>
#include <cstdio>

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>


/** Main namespace for ReaK */
namespace ReaK {


class shared_object; //forward-declaration.

/** Main namespace for ReaK's Serialization */
class serializable; //forward-declaration

/** Main namespace for ReaK's Run-time Type Identification (RTTI) */
namespace rtti {

typedef ReaK::shared_ptr<shared_object> (RK_CALL *construct_ptr)();

//this is really the only thing that the class needs to define
template <typename T>
struct get_type_id {
  BOOST_STATIC_CONSTANT(unsigned int, ID = T::rk_rtti_ID);
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = T::rk_rtti_TypeName;
#else
  static const char* type_name() BOOST_NOEXCEPT { return T::rk_rtti_TypeName(); };
#endif
  
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return T::rk_rtti_CreatePtr(); };
  
  typedef const serializable& save_type;
  typedef serializable& load_type;
};


#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS

BOOST_CONSTEXPR_OR_CONST auto lsl_left_bracket = RK_LSA("<");
BOOST_CONSTEXPR_OR_CONST auto lsl_right_bracket = RK_LSA(">");
BOOST_CONSTEXPR_OR_CONST auto lsl_comma = RK_LSA(",");

#endif


template <unsigned int U>
struct get_type_id< boost::mpl::integral_c<unsigned int,U> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = U);
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = ct_itoa<ID>::text;
#else
  static const char* type_name() BOOST_NOEXCEPT { static char type_name_buf[32]; std::sprintf(type_name_buf,"%d",ID); return type_name_buf; };
#endif
  
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::true_ > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("true");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "true"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::false_ > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("false");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "false"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <int I>
struct get_type_id< boost::mpl::int_<I> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = I);
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = ct_itoa<ID>::text;
#else
  static const char* type_name() BOOST_NOEXCEPT { static char type_name_buf[32]; std::sprintf(type_name_buf,"%d",ID); return type_name_buf; };
#endif
  
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};


struct null_type_id { 
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0);
};

struct null_type_info {
  typedef null_type_id type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("");
#else
  static const char* type_name() BOOST_NOEXCEPT { return ""; };
#endif
};

template <>
struct get_type_id< null_type_info > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("");
#else
  static const char* type_name() BOOST_NOEXCEPT { return ""; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};


namespace {
  
  template <typename T, typename Tail = null_type_id>
  struct type_id {
    typedef Tail tail;
    BOOST_STATIC_CONSTANT(unsigned int, ID = ::ReaK::rtti::get_type_id<T>::ID);
  };
  
};

template <typename Tail>
struct get_type_name_tail {
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto value = lsl_comma + Tail::type_name;
#else
  static std::string value() { 
    std::string result = ",";
    result += Tail::type_name();
    return result; // NRVO
  };
#endif
};

template <>
struct get_type_name_tail<null_type_info> {
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto value = RK_LSA("");
#else
  static const char* value() BOOST_NOEXCEPT { return ""; };
#endif
};

template <typename T, typename Tail = null_type_info>
struct get_type_info {
  typedef type_id<T, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<T>::type_name + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id<T>::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO 
  };
#endif
};


namespace {

#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES

template <typename... T> struct get_type_info_seq;

template <typename T1>
struct get_type_info_seq<T1> {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef get_type_info<T1, Tail> type;
  };
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<T1>::type_name;
#else
  static std::string type_name() { return get_type_id<T1>::type_name(); };
#endif
};


template <typename T1, typename... T>
struct get_type_info_seq<T1, T...> {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef get_type_info<T1, 
      typename get_type_info_seq< T... >::template with_tail<Tail>::type > type;
  };
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<T1>::type_name + lsl_comma + get_type_info_seq< T... >::type_name;
#else
  static std::string type_name() { 
    std::string result = get_type_id<T1>::type_name();
    result += ",";
    result += get_type_info_seq< T... >::type_name();
    return result; //NRVO
  };
#endif
};

#else

template <typename T1, typename T2 = void, typename T3 = void, typename T4 = void, typename T5 = void, 
          typename T6 = void, typename T7 = void, typename T8 = void, typename T9 = void, typename T10 = void>
struct get_type_info_seq;

template <typename T1>
struct get_type_info_seq<T1,void,void,void,void,void,void,void,void,void> {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef get_type_info<T1, Tail> type;
  };
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<T1>::type_name;
#else
  static std::string type_name() { return get_type_id<T1>::type_name(); };
#endif
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10>
struct get_type_info_seq {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef get_type_info<T1, 
      typename get_type_info_seq<T2,T3,T4,T5,T6,T7,T8,T9,T10>::template with_tail<Tail>::type > type;
  };
  
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<T1>::type_name + lsl_comma + get_type_info_seq<T2,T3,T4,T5,T6,T7,T8,T9,T10>::type_name;
#else
  static std::string type_name() { 
    std::string result = get_type_id<T1>::type_name();
    result += ",";
    result += get_type_info_seq<T2,T3,T4,T5,T6,T7,T8,T9,T10>::type_name();
    return result; //NRVO
  };
#endif
};

#endif

};



#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES

template <template <typename...> class U, typename Tail, typename... T>
struct get_type_info< U<T...>, Tail > {
  typedef type_id< U<T...>, 
    typename get_type_info_seq<T...>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T...> >::type_name + lsl_left_bracket + get_type_info_seq<T...>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T...> >::type_name();
    result += "<";
    result += get_type_info_seq<T...>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

#else

template <template <typename> class U, typename T, typename Tail>
struct get_type_info< U<T>, Tail > {
  typedef type_id< U<T>, 
    typename get_type_info_seq<T>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T> >::type_name + lsl_left_bracket + get_type_info_seq<T>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T> >::type_name();
    result += "<";
    result += get_type_info_seq<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename> class U, 
          typename T1, typename T2, typename Tail>
struct get_type_info< U<T1,T2>, Tail > {
  typedef type_id< U<T1,T2>, 
    typename get_type_info_seq<T1,T2>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename Tail>
struct get_type_info< U<T1,T2,T3>, Tail > {
  typedef type_id< U<T1,T2,T3>, 
    typename get_type_info_seq<T1,T2,T3>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename Tail>
struct get_type_info< U<T1,T2,T3,T4>, Tail > {
  typedef type_id< U<T1,T2,T3,T4>, 
    typename get_type_info_seq<T1,T2,T3,T4>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename T5, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename,typename> class U, 
                    typename T1, typename T2, typename T3, typename T4, typename T5, 
                    typename T6, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5,T6>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5,T6>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5,T6>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5,T6> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5,T6>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5,T6> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5,T6>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5,T6,T7>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5,T6,T7>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5,T6,T7>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5,T6,T7> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5,T6,T7>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5,T6,T7> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5,T6,T7>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5,T6,T7,T8>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5,T6,T7,T8>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename,typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5,T6,T7,T8,T9>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

template <template <typename,typename,typename,typename,typename,typename,typename,typename,typename,typename> class U, 
          typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10, typename Tail>
struct get_type_info< U<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>, Tail > {
  typedef type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>, 
    typename get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> >::type_name + lsl_left_bracket + get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< U<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> >::type_name();
    result += "<";
    result += get_type_info_seq<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

#endif



template <unsigned int U, typename Tail>
struct get_type_info< boost::mpl::integral_c<unsigned int,U>, Tail > {
  typedef type_id<boost::mpl::integral_c<unsigned int,U>, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< boost::mpl::integral_c<unsigned int,U> >::type_name + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< boost::mpl::integral_c<unsigned int,U> >::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

template <typename Tail>
struct get_type_info< boost::mpl::true_, Tail > {
  typedef type_id<boost::mpl::true_, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<boost::mpl::true_>::type_name + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id<boost::mpl::true_>::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

template <typename Tail>
struct get_type_info< boost::mpl::false_, Tail > {
  typedef type_id<boost::mpl::false_, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id<boost::mpl::false_>::type_name + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id<boost::mpl::false_>::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

template <int I, typename Tail>
struct get_type_info< boost::mpl::int_<I>, Tail > {
  typedef type_id<boost::mpl::int_<I>, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< boost::mpl::int_<I> >::type_name + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< boost::mpl::int_<I> >::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};



/**
 * This class is used to identify and register in the platform a new object
 * type. When a plugin is loaded, at initialization of the plugin, it is responsible
 * for registering all the new types it contains.
 */
class so_type {
private:
  so_type(const so_type&);
  so_type& operator =(const so_type&);
  
public:
  static so_type* createTypeInfo(unsigned int aTypeVersion, unsigned int* aTypeID, 
                                 const std::string& aTypeName, construct_ptr aConstruct);
  
  static unsigned int* createTypeID(unsigned int aTypeIDSize);
  
protected:
  
  so_type() { };
  
  static bool compare_equal(const unsigned int* pid1, const unsigned int* pid2);
  
public:
  
  ///This function adds a Descendant of this.
  so_type* addDescendant(so_type* aObj);
  
  so_type* addAncestor(so_type* aObj);
  
  ///This function finds a TypeID in the descendants (recusively) of this.
  so_type* findDescendant(const unsigned int* aTypeID );
  
  ///This function gets the number of direct descendants of this.
  unsigned int getDirectDescendantCount();
  
  ///This function gets a type record by index in the direct descendants of this.
  so_type* getDirectDescendant(unsigned int aIndex);
  
  ///This function checks if a typeID is parent to this.
  so_type* findAncestor(const unsigned int* aTypeID );
  
  ///This function finds a TypeID in the descendants (recusively) of this.
  so_type* findDescendant(so_type* aTypeID );
  
  ///This function checks if a typeID is parent to this.
  so_type* findAncestor(so_type* aTypeID );
  
  ///This function inserts this into a global repo.
  void insertToRepo(so_type* aRepo);
  
  const unsigned int* TypeID_begin() const;
  
  unsigned int TypeVersion() const;
  
  const std::string& TypeName() const;
  
  ReaK::shared_ptr<shared_object> CreateObject() const;
  
  bool isConcrete() const;
  
  friend bool operator ==(const so_type& t1, const so_type& t2) {
    return compare_equal(t1.TypeID_begin(),t2.TypeID_begin());
  };
  
  friend bool operator !=(const so_type& t1, const so_type& t2) {
    return !compare_equal(t1.TypeID_begin(),t2.TypeID_begin());
  };
  
};


struct so_type_ptr {
  so_type* ptr;
  so_type_ptr(so_type* aPtr) : ptr(aPtr) { };
  ~so_type_ptr();
};


namespace {
  
  template <typename T> struct get_type_id_prop;
  
  template <>
  struct get_type_id_prop<null_type_id> {
    BOOST_STATIC_CONSTANT(unsigned int, count = 1);
    static unsigned int at(unsigned int) { return 0; };
  };
  
  template <typename T>
  struct get_type_id_prop {
    BOOST_STATIC_CONSTANT(unsigned int, count = get_type_id_prop< typename T::tail >::count + 1);
    static unsigned int at(unsigned int i) {
      if(i == 0)
        return T::ID;
      else
        return get_type_id_prop< typename T::tail >::at(--i);
    };
  };
  
  
  template <typename T>
  so_type_ptr create_type_descriptor(unsigned int aVersion = 1) {
    typedef unsigned int SizeType;
    const SizeType TypeIDLength = get_type_id_prop< typename get_type_info<T>::type >::count;
    SizeType* typeID = so_type::createTypeID(TypeIDLength);
    for(SizeType i = 0; i < TypeIDLength; ++i)
      typeID[i] = get_type_id_prop< typename get_type_info<T>::type >::at(i);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
    constexpr auto tname = get_type_info<T>::type_name;
    return so_type_ptr(so_type::createTypeInfo(aVersion, typeID, tname.to_string(), get_type_id<T>::CreatePtr()));
#else
    return so_type_ptr(so_type::createTypeInfo(aVersion, typeID, get_type_info<T>::type_name(), get_type_id<T>::CreatePtr()));
#endif
  };
  
};

so_type_ptr create_dummy_so_type(const unsigned int* aTypeID);



};

};

#endif


#include "typed_primitives.hpp"
#include "typed_containers.hpp"










