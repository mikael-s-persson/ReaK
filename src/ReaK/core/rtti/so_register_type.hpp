/**
 * \file so_register_type.hpp
 * 
 * This library declares the object which, when instantiated, registers a type to the ReaK::rtti system.
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

#ifndef REAK_SO_REGISTER_TYPE_HPP
#define REAK_SO_REGISTER_TYPE_HPP

#include <ReaK/core/base/defs.hpp>

#include "so_type.hpp"
#include "so_type_repo.hpp"

namespace ReaK {

namespace rtti {

namespace detail {
  
  struct null_base_type { };
  
  template <typename Base, typename Tail = null_base_type>
  struct base_type_list {
    typedef Tail tail;
    typedef Base type; 
  };
  
//   template <typename Base>
//   struct base_type_list<Base,null_base_type> {
//     typedef Base type;
//   };
  
  template <typename T>
  struct base_type_count {
    BOOST_STATIC_CONSTANT(unsigned int, value = base_type_count< typename T::tail >::value + 1);
  };
  
  template <>
  struct base_type_count<null_base_type> {
    BOOST_STATIC_CONSTANT(unsigned int, value = 1);
  };
  
  template <typename T>
  struct add_base_type {
    static void to(so_type::shared_pointer& aObj) {
      aObj->addAncestor(aObj,T::type::getStaticObjectType());
      add_base_type<typename T::tail>::to(aObj);
    };
  };
  
  template <>
  struct add_base_type<null_base_type> {
    static void to(so_type::shared_pointer&) { };
  };
  
};

template <typename T, unsigned int Version, typename BaseList = detail::null_base_type>
struct register_type {
public:
  
  struct register_type_impl {
    so_type::shared_pointer ptr;
    register_type_impl() : ptr(new so_type_descriptor<T,Version>(),scoped_deleter()) {
      detail::add_base_type<BaseList>::to(ptr);
      so_type_repo::getInstance().addType(ptr);
    };
  };
  static const register_type_impl impl;
  
  register_type() {  impl; /* force the instantiations! */ };
  
  
};

template <typename T, unsigned int Version, typename BaseList>
const typename register_type<T,Version,BaseList>::register_type_impl register_type<T,Version,BaseList>::impl;



};

};

#endif










