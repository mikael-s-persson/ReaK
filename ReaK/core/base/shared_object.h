/**
 *\file shared_object.h
 *
 * This library declares the basic class "ReaK::shared_object" which allows all descendants
 * to be safely shared across executable modules (given the binary compatibility of the boost
 * smart pointers is kept, this is out of my control). All classes in the ReaK platform
 * should at least be a descendant of this class (except some low-level classes, but no end-user
 * defined classes). Usually, however, end-users will prefer the ReaK::named_object class which
 * adds a name to the objects (and it is the next descendant after ReaK::shared_object).
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date february 2010
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

#ifndef REAK_CORE_BASE_SHARED_OBJECT_H_
#define REAK_CORE_BASE_SHARED_OBJECT_H_

#include <memory>

#include "ReaK/core/serialization/serializable.h"

/** Main namespace for ReaK */
namespace ReaK {

/**
 * This basic class allows all descendants
 * to be safely shared across executable modules (given the binary compatibility of the boost
 * smart pointers is kept, but this is out of the my control). All classes in the ReaK platform
 * should at least be a descendant of this class (except some low-level classes, but no end-user
 * defined classes). Usually, however, end-users will prefer the ReaK::named_object class which
 * adds a name to the objects (and it is the next descendant after ReaK::shared_object).
 */
class shared_object : public serializable {
 public:
  ~shared_object() override = default;

  RK_RTTI_MAKE_ABSTRACT_1BASE(shared_object, 0x80000001, 1, "shared_object",
                              serializable)
};

class empty_base_object {};

/**
 * This structure is a simple callable structure that does nothing. It acts as a place-holder for
 * a special deleter for a std::shared_ptr, in this case this "null deleter" simply does not delete anything.
 * This is useful for an object for which you want a std::shared_ptr to, but that is not going to be
 * deleted by this std::shared_ptr branch (i.e. will be deleted by another std::shared_ptr branch). This can
 * be used to break cycles in the object hierarchy.
 */
struct null_deleter {
  void operator()(void const* /*unused*/) const {}
};

}  // namespace ReaK

#endif  // REAK_CORE_BASE_SHARED_OBJECT_H_
