/**
 * \file so_register_type.h
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

#ifndef REAK_CORE_RTTI_SO_REGISTER_TYPE_H_
#define REAK_CORE_RTTI_SO_REGISTER_TYPE_H_

#include "ReaK/core/rtti/so_type_repo.h"

namespace ReaK::rtti {

namespace {

template <typename T>
struct register_type {

  struct register_type_impl {
    explicit register_type_impl(int /*unused*/) {
      so_type_repo::getInstance().addType(T::getStaticObjectType());
    }
  };
  static const register_type_impl impl;

  /*
   * Explanation of registration scheme:
   * The register_type template is instantiated for every class that has a RK_RTTI_REGISTER
   * MACRO in its declaration. That MACRO inserts an invocation of the "impl" static member
   * within a virtual function (the getObjectType() function), which has the effect of forcing
   * this static member to be considered as "used somewhere", which means the compiler must
   * generate its static initialization code, i.e., the constructor of register_type_impl
   * must get executed before entering main, and in that constructor, the type is registered
   * to the global repository of types. This sequence is the key to making this work because
   * the presence of "impl" invocation somewhere within the body of a virtual function makes it
   * so that even when the class being registered is a template, it must be created, because
   * only virtual functions are guaranteed to be instantiated (even if never used) in a class
   * template instantiation.
   *
   * The use of get_ptr() is there to avoid a static initialization order fiasco between a class
   * and its base-classes, which might not have been initialized yet.
   */
};

template <typename T>
const typename register_type<T>::register_type_impl register_type<T>::impl(0);
}  // namespace

}  // namespace ReaK::rtti

#endif  // REAK_CORE_RTTI_SO_REGISTER_TYPE_H_
