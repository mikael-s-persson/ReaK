/**
 * \file typed_containers.hpp
 * 
 * This library associates type information to STL containers types of C++ standard libraries.
 * This allows STL containers to be integrated to the ReaK::rtti system.
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

#ifndef REAK_TYPED_CONTAINERS_HPP
#define REAK_TYPED_CONTAINERS_HPP

#include "so_type.hpp"


#include <vector>

namespace ReaK { namespace rtti {

template <typename T>
struct get_type_id< std::vector<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000008);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::vector");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::vector"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::vector<T>& save_type;
  typedef std::vector<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::vector<T>, Tail > {
  typedef type_id<  std::vector<T> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::vector<T> >::type_name + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::vector<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };


#include <list>

namespace ReaK { namespace rtti {

template <typename T>
struct get_type_id< std::list<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000009);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::list");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::list"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::list<T>& save_type;
  typedef std::list<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::list<T>, Tail > {
  typedef type_id<  std::list<T> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::list<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::list<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

}; };


#include <map>

namespace ReaK { namespace rtti {

template <typename Key, typename T>
struct get_type_id< std::map<Key,T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000A);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::map");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::map"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::map<Key,T>& save_type;
  typedef std::map<Key,T>& load_type;
};

template <typename Key, typename T, typename Tail>
struct get_type_info< std::map<Key,T>, Tail > {
  typedef type_id< std::map<Key,T> , typename get_type_info<Key, get_type_info<T, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::map<Key,T> >::type_name
    + lsl_left_bracket + get_type_id<Key>::type_name + lsl_comma + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::map<Key,T> >::type_name();
    result += "<";
    result += get_type_id<Key>::type_name();
    result += ",";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


template <typename Key, typename T>
struct get_type_id< std::multimap<Key,T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000D);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::multimap");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::multimap"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::multimap<Key,T>& save_type;
  typedef std::multimap<Key,T>& load_type;
};

template <typename Key, typename T, typename Tail>
struct get_type_info< std::multimap<Key,T>, Tail > {
  typedef type_id< std::multimap<Key,T> , typename get_type_info<Key, get_type_info<T, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::multimap<Key,T> >::type_name
    + lsl_left_bracket + get_type_id<Key>::type_name + lsl_comma + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::multimap<Key,T> >::type_name();
    result += "<";
    result += get_type_id<Key>::type_name();
    result += ",";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };


#include <set>

namespace ReaK { namespace rtti {

template <typename T>
struct get_type_id< std::set<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000B);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::set");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::set"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::set<T>& save_type;
  typedef std::set<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::set<T>, Tail > {
  typedef type_id< std::set<T> , typename get_type_info<T, Tail>::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::set<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::set<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


template <typename T>
struct get_type_id< std::multiset<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000E);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::multiset");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::multiset"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::multiset<T>& save_type;
  typedef std::multiset<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::multiset<T>, Tail > {
  typedef type_id< std::multiset<T> , typename get_type_info<T, Tail>::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::multiset<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::multiset<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };


#include <utility>

namespace ReaK { namespace rtti {

template <typename T1, typename T2>
struct get_type_id< std::pair<T1,T2> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000C);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::pair");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::pair"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::pair<T1,T2>& save_type;
  typedef std::pair<T1,T2>& load_type;
};

template <typename T1, typename T2, typename Tail>
struct get_type_info< std::pair<T1,T2>, Tail > {
  typedef type_id< std::pair<T1,T2> , typename get_type_info<T1, get_type_info<T2, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::pair<T1,T2> >::type_name
    + lsl_left_bracket + get_type_id<T1>::type_name + lsl_comma + get_type_id<T2>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::pair<T1,T2> >::type_name();
    result += "<";
    result += get_type_id<T1>::type_name();
    result += ",";
    result += get_type_id<T2>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };



#ifndef BOOST_NO_CXX11_HDR_FORWARD_LIST

#include <forward_list>

namespace ReaK { namespace rtti {

template <typename T>
struct get_type_id< std::forward_list<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000040);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::forward_list");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::forward_list"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::forward_list<T>& save_type;
  typedef std::forward_list<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::forward_list<T>, Tail > {
  typedef type_id< std::forward_list<T> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::forward_list<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::forward_list<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

}; };

#endif

#ifndef BOOST_NO_CXX11_HDR_ARRAY

#include <array>

namespace ReaK { namespace rtti {

template <typename T, std::size_t N>
struct get_type_id< std::array<T,N> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000041);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::array");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::array"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::array<T,N>& save_type;
  typedef std::array<T,N>& load_type;
};

template <typename T, std::size_t N, typename Tail>
struct get_type_info< std::array<T,N>, Tail > {
  typedef type_id< std::array<T,N> , typename get_type_info<T, 
                                                   get_type_info<boost::mpl::integral_c<unsigned int,N>, Tail> >::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::array<T,N> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_comma 
    + get_type_id< boost::mpl::integral_c<unsigned int,N> >::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::array<T,N> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ",";
    result += get_type_id< boost::mpl::integral_c<unsigned int,N> >::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

}; };

#endif

#if !defined(BOOST_NO_CXX11_HDR_TUPLE) && !defined(BOOST_NO_CXX11_VARIADIC_TEMPLATES)

#include <tuple>

namespace ReaK { namespace rtti {

template <typename... T>
struct get_type_id< std::tuple< T... > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000042);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::tuple");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::tuple"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::tuple< T... >& save_type;
  typedef std::tuple< T... >& load_type;
};

template <typename Tail, typename... T>
struct get_type_info< std::tuple< T... >, Tail > {
  typedef type_id< std::tuple<T...>, 
    typename get_type_info_seq<T...>::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::tuple< T... > >::type_name
    + lsl_left_bracket + get_type_info_seq< T... >::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::tuple< T... > >::type_name();
    result += "<";
    result += get_type_info_seq< T... >::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; // NRVO
  };
#endif
};

}; };


#endif

#ifndef BOOST_NO_CXX11_HDR_UNORDERED_MAP

#include <unordered_map>

namespace ReaK { namespace rtti {

template <typename Key, typename T>
struct get_type_id< std::unordered_map<Key,T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000044);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::unordered_map");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::unordered_map"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::unordered_map<Key,T>& save_type;
  typedef std::unordered_map<Key,T>& load_type;
};

template <typename Key, typename T, typename Tail>
struct get_type_info< std::unordered_map<Key,T>, Tail > {
  typedef type_id< std::unordered_map<Key,T> , typename get_type_info<Key, get_type_info<T, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::unordered_map<Key,T> >::type_name
    + lsl_left_bracket + get_type_id<Key>::type_name + lsl_comma + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::unordered_map<Key,T> >::type_name();
    result += "<";
    result += get_type_id<Key>::type_name();
    result += ",";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


template <typename Key, typename T>
struct get_type_id< std::unordered_multimap<Key,T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000046);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::unordered_multimap");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::unordered_multimap"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::unordered_multimap<Key,T>& save_type;
  typedef std::unordered_multimap<Key,T>& load_type;
};

template <typename Key, typename T, typename Tail>
struct get_type_info< std::unordered_multimap<Key,T>, Tail > {
  typedef type_id< std::unordered_multimap<Key,T> , typename get_type_info<Key, get_type_info<T, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::unordered_multimap<Key,T> >::type_name
    + lsl_left_bracket + get_type_id<Key>::type_name + lsl_comma + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::unordered_multimap<Key,T> >::type_name();
    result += "<";
    result += get_type_id<Key>::type_name();
    result += ",";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };

#endif

#ifndef BOOST_NO_CXX11_HDR_UNORDERED_SET

#include <unordered_set>

namespace ReaK { namespace rtti {

template <typename T>
struct get_type_id< std::unordered_set<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000045);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::unordered_set");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::unordered_set"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::unordered_set<T>& save_type;
  typedef std::unordered_set<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::unordered_set<T>, Tail > {
  typedef type_id< std::unordered_set<T> , typename get_type_info<T, Tail>::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::unordered_set<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::unordered_set<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


template <typename T>
struct get_type_id< std::unordered_multiset<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000047);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::unordered_multiset");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::unordered_multiset"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::unordered_multiset<T>& save_type;
  typedef std::unordered_multiset<T>& load_type;
};

template <typename T, typename Tail>
struct get_type_info< std::unordered_multiset<T>, Tail > {
  typedef type_id< std::unordered_multiset<T> , typename get_type_info<T, Tail>::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::unordered_multiset<T> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::unordered_multiset<T> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

}; };

#endif


#endif
