/**
 *\file shared_object.hpp
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

#ifndef REAK_SHARED_OBJECT_HPP
#define REAK_SHARED_OBJECT_HPP

#include "defs.hpp"

#include "serializable.hpp"
#include "shared_object_base.hpp"

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
class shared_object : public shared_object_base, public serialization::serializable {
  public:
    virtual void RK_CALL destroy() { delete this; };
    
    virtual ~shared_object() { };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(shared_object,0x80000001,1,"shared_object",serialization::serializable)
    
};

class empty_base_object { };

#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES

template <typename T>
ReaK::shared_ptr<T> rk_create() {
  return ReaK::shared_ptr<T>(new T(), ReaK::scoped_deleter());
};

#else

template <typename T, typename... Args>
ReaK::shared_ptr<T> rk_create(Args&&... args) {
  return ReaK::shared_ptr<T>(new T(std::forward<Args>(args)...), ReaK::scoped_deleter());
};

#endif

template <typename T>
ReaK::shared_ptr<T> rk_share(T& t) {
  return ReaK::shared_ptr<T>(&t, ReaK::null_deleter());
};


namespace rtti {

#ifdef BOOST_NO_CXX11_SMART_PTR

template <typename Y>
ReaK::shared_ptr<Y> rk_static_ptr_cast(const ReaK::shared_ptr<ReaK::shared_object_base>& p) {
  return boost::static_pointer_cast<Y>(p);
};

#else 

template <typename Y>
ReaK::shared_ptr<Y> rk_static_ptr_cast(const ReaK::shared_ptr<ReaK::shared_object_base>& p) {
  return std::static_pointer_cast<Y>(p);
};

template <typename Y, typename Deleter>
ReaK::unique_ptr<Y,Deleter> rk_dynamic_ptr_cast(ReaK::unique_ptr<ReaK::shared_object_base,Deleter>&& p) {
  void* tmp = static_cast<ReaK::shared_object*>(p.get())->castTo(Y::getStaticObjectType());
  if(tmp) {
    ReaK::unique_ptr<Y,Deleter> r(tmp, p.get_deleter());
    p.release();
    return std::move(r);
  } else
    return ReaK::unique_ptr<Y,Deleter>();
};

#endif

/**
 * This function replaces the standard C++ dynamic cast (i.e. dynamic_cast<>()) and furthermore, also
 * replaces the boost::shared_ptr dynamic cast (i.e. boost::dynamic_pointer_cast<>()). This new function
 * for dynamic casting is required in order for ReaK::rtti system to take precedence over the C++ RTTI
 * because, unlike the C++ standard version of RTTI, this implementation will work across executable modules,
 * and thus, allow objects to be shared between modules with full dynamic up- and down- casting capabilities.
 * \note this function is a special overload for the shared_object_base class pointers, the general version is
 *       found in the "typed_object.hpp" library, as part of the ReaK::rtti system.
 */
template <typename Y>
shared_ptr<Y> rk_dynamic_ptr_cast(const shared_ptr<shared_object_base>& p) {
  return shared_ptr<Y>(p,reinterpret_cast<Y*>(rk_static_ptr_cast<shared_object>(p)->castTo(Y::getStaticObjectType())));
};


};


};

#endif
















