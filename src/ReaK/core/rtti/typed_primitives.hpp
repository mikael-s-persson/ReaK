/**
 * \file typed_primitives.hpp
 * 
 * This library associates type information to primitive "built-in" types of C++.
 * This allows built-in types to be integrated to the ReaK::rtti system.
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

#ifndef REAK_TYPED_PRIMITIVES_HPP
#define REAK_TYPED_PRIMITIVES_HPP

#include "so_type.hpp"

#ifndef RK_ENABLE_CXX0X_FEATURES

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#else

#include <memory>

#endif

namespace ReaK {

namespace rtti {


template <>
struct get_type_id<int> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000001);
  static std::string type_name() { return "int"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef int save_type;
  typedef int& load_type;
};

template <>
struct get_type_id<unsigned int> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000002);
  static std::string type_name() { return "unsigned int"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef unsigned int save_type;
  typedef unsigned int& load_type;
};

template <>
struct get_type_id<long unsigned int> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000000F);
  static std::string type_name() { return "long unsigned int"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const long unsigned int& save_type;
  typedef long unsigned int& load_type;
};

template <>
struct get_type_id<float> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000003);
  static std::string type_name() { return "float"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef float save_type;
  typedef float& load_type;
};
  
template <>
struct get_type_id<double> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000004);
  static std::string type_name() { return "double"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef double save_type;
  typedef double& load_type;
};

template <>
struct get_type_id<bool> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000005);
  static std::string type_name() { return "bool"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef bool save_type;
  typedef bool& load_type;
};

template <>
struct get_type_id<std::string> {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000006);
  static std::string type_name() { return "std::string"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const std::string& save_type;
  typedef std::string& load_type;
};

#ifndef RK_ENABLE_CXX0X_FEATURES

template <typename T>
struct get_type_id< boost::shared_ptr<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
  static std::string type_name() { return "boost::shared_ptr"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const boost::shared_ptr<T>& save_type;
  typedef boost::shared_ptr<T>& load_type;
};

template <typename T>
struct get_type_id< boost::weak_ptr<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
  static std::string type_name() { return "boost::weak_ptr"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const boost::weak_ptr<T>& save_type;
  typedef boost::weak_ptr<T>& load_type;
};

#else

template <typename T>
struct get_type_id< std::shared_ptr<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
  static std::string type_name() { return "std::shared_ptr"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const std::shared_ptr<T>& save_type;
  typedef std::shared_ptr<T>& load_type;
};

template <typename T>
struct get_type_id< std::weak_ptr<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
  static std::string type_name() { return "std::weak_ptr"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const std::weak_ptr<T>& save_type;
  typedef std::weak_ptr<T>& load_type;
};

template <typename T>
struct get_type_id< std::unique_ptr<T> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
  static std::string type_name() { return "std::unique_ptr"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const std::unique_ptr<T>& save_type;
  typedef std::unique_ptr<T>& load_type;
};

#endif


};

};


#endif






