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
#include <list>
#include <set>
#include <map>


namespace ReaK {

namespace rtti {


template <typename T, typename Allocator>
struct get_type_id< std::vector<T,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000008);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::vector");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::vector"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::vector<T,Allocator>& save_type;
  typedef std::vector<T,Allocator>& load_type;
};

template <typename T, typename Allocator, typename Tail>
struct get_type_info< std::vector<T,Allocator>, Tail > {
  typedef type_id<  std::vector<T,Allocator> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::vector<T,Allocator> >::type_name + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::vector<T,Allocator> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


template <typename T, typename Allocator>
struct get_type_id< std::list<T,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000009);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::list");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::list"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::list<T,Allocator>& save_type;
  typedef std::list<T,Allocator>& load_type;
};

template <typename T, typename Allocator, typename Tail>
struct get_type_info< std::list<T,Allocator>, Tail > {
  typedef type_id<  std::list<T,Allocator> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::list<T,Allocator> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::list<T,Allocator> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};


template <typename Key, typename T, typename Compare, typename Allocator>
struct get_type_id< std::map<Key,T,Compare,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000A);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::map");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::map"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::map<Key,T,Compare,Allocator>& save_type;
  typedef std::map<Key,T,Compare,Allocator>& load_type;
};

template <typename Key, typename T, typename Compare, typename Allocator, typename Tail>
struct get_type_info< std::map<Key,T,Compare,Allocator>, Tail > {
  typedef type_id< std::map<Key,T,Compare,Allocator> , typename get_type_info<Key, get_type_info<T, Tail> >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::map<Key,T,Compare,Allocator> >::type_name
    + lsl_left_bracket + get_type_id<Key>::type_name + lsl_comma + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::map<Key,T,Compare,Allocator> >::type_name();
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


template <typename T, typename Compare, typename Allocator>
struct get_type_id< std::set<T,Compare,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000B);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("std::set");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "std::set"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const std::set<T,Compare,Allocator>& save_type;
  typedef std::set<T,Compare,Allocator>& load_type;
};

template <typename T, typename Compare, typename Allocator, typename Tail>
struct get_type_info< std::set<T,Compare,Allocator>, Tail > {
  typedef type_id< std::set<T,Compare,Allocator> , typename get_type_info<T, Tail>::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::set<T,Compare,Allocator> >::type_name
    + lsl_left_bracket + get_type_id<T>::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< std::set<T,Compare,Allocator> >::type_name();
    result += "<";
    result += get_type_id<T>::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};



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






};

};


#endif
