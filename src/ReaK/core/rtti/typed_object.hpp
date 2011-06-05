/**
 * \file typed_object.hpp
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

#ifndef TYPED_OBJECT_HPP
#define TYPED_OBJECT_HPP

#include "base/defs.hpp"

#include "so_type_repo.hpp"
#include "so_register_type.hpp"

#include "typed_primitives.hpp"
#include "typed_containers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/stringize.hpp>


namespace ReaK {

namespace rtti {


/**
 * The basic class "ReaK::rtti::typed_object" allows all descendants
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
    virtual void* RK_CALL castTo(const boost::shared_ptr<so_type>& aTypeID) { 
      if(*(aTypeID->TypeID_begin()) == 0) 
	return reinterpret_cast<void*>(this); 
      else 
	return NULL;
    };
    
    /**
     * This method is used to perform up- and down- casting of const-object pointers via a virtual call.
     */
    virtual const void* RK_CALL castTo(const boost::shared_ptr<so_type>& aTypeID) const { 
      if(*(aTypeID->TypeID_begin()) == 0) 
	return reinterpret_cast<const void*>(this); 
      else 
	return NULL;
    };

    virtual ~typed_object() { RK_NOTICE(8,"Typed object destructor reached!"); };

       
    /** This method fetches the object type structure from the ReaK::rtti system or creates it if it has not been registered yet. */
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const {
      return boost::shared_ptr<so_type>();
    };
    /** This method fetches the object type structure from the ReaK::rtti system or creates it if it has not been registered yet. */
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() {
      return boost::shared_ptr<so_type>();
    };

};

/**
 * This function replaces the standard C++ dynamic cast (i.e. dynamic_cast<>()) and furthermore, also
 * replaces the boost::shared_ptr dynamic cast (i.e. boost::dynamic_pointer_cast<>()). This new function
 * for dynamic casting is required in order for ReaK::rtti system to take precedence over the C++ RTTI
 * because, unlike the C++ standard version of RTTI, this implementation will work across executable modules,
 * and thus, allow objects to be shared between modules with full dynamic up- and down- casting capabilities.
 */
template <class Y,class U>
boost::shared_ptr<Y> rk_dynamic_ptr_cast(boost::shared_ptr<U> p) {
  return boost::shared_ptr<Y>(p,reinterpret_cast<Y*>(p->castTo(Y::getStaticObjectType())));
};

/// This MACRO creates a static (no-parameter) factory function for the current class CLASS_NAME.
#define RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    static boost::shared_ptr<ReaK::shared_object> RK_CALL Create() { \
      return boost::shared_ptr< CLASS_NAME >(new CLASS_NAME(), ReaK::scoped_deleter()); \
    }; \
    static ReaK::rtti::construct_ptr rk_rtti_CreatePtr() { return &Create; };
    
/// This MACRO registers a custom (no-parameter) factory function pointer for the current class.
#define RK_RTTI_REGISTER_CUSTOM_FACTORY(CLASS_FACTORY) \
    static ReaK::rtti::construct_ptr rk_rtti_CreatePtr() { return CLASS_FACTORY; };

/// This MACRO creates the static elements for the current class that registers the CLASS_ID and CLASS_NAME.
#define RK_RTTI_REGISTER_CLASS_ID(CLASS_NAME, CLASS_ID) \
    static const unsigned int rk_rtti_ID = CLASS_ID; \
    static std::string rk_rtti_TypeName() { return CLASS_NAME; };
    
/// This MACRO creates the static elements for the current class to be added to the global type registry (it is guaranteed to be added if the class is instantiated).
#define RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME,CLASS_VERSION) \
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() { \
      return ReaK::rtti::register_type< CLASS_NAME , CLASS_VERSION >::impl.ptr; \
    }; \
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const { \
      return CLASS_NAME::getStaticObjectType(); \
    };\
    virtual void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<void*>(this); \
      else \
	return NULL; \
    };\
    virtual const void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) const { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<const void*>(this); \
      else \
	return NULL; \
    };

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is guaranteed to be added if the class is instantiated).
#define RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME) \
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() { \
      return ReaK::rtti::register_type< CLASS_NAME , CLASS_VERSION , ReaK::rtti::detail::base_type_list< BASE_NAME > >::impl.ptr; \
    }; \
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const { \
      return CLASS_NAME::getStaticObjectType(); \
    };\
    virtual void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<void*>(this); \
      else \
	return BASE_NAME::castTo(aTypeID); \
    };\
    virtual const void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) const { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<const void*>(this); \
      else \
	return BASE_NAME::castTo(aTypeID); \
    };

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is guaranteed to be added if the class is instantiated).
#define RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2) \
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() { \
      return ReaK::rtti::register_type< CLASS_NAME , CLASS_VERSION , ReaK::rtti::detail::base_type_list< BASE_NAME1, ReaK::rtti::detail::base_type_list< BASE_NAME2 > > >::impl.ptr; \
    }; \
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const { \
      return CLASS_NAME::getStaticObjectType(); \
    };\
    virtual void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<void*>(this); \
      else { \
	void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	return BASE_NAME2::castTo(aTypeID); \
      }; \
    };\
    virtual const void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) const { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<const void*>(this); \
      else { \
	const void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	return BASE_NAME2::castTo(aTypeID); \
      }; \
    };
    
/// This MACRO creates the static elements for the current class to be added to the global type registry (it is guaranteed to be added if the class is instantiated).
#define RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3) \
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() { \
      return ReaK::rtti::register_type< CLASS_NAME , CLASS_VERSION , ReaK::rtti::detail::base_type_list< BASE_NAME1, ReaK::rtti::detail::base_type_list< BASE_NAME2 , ReaK::rtti::detail::base_type_list< BASE_NAME3 > > > >::impl.ptr; \
    }; \
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const { \
      return CLASS_NAME::getStaticObjectType(); \
    };\
    virtual void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<void*>(this); \
      else { \
	void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME2::castTo(aTypeID); \
	if(result) return result; \
        return BASE_NAME3::castTo(aTypeID); \
      }; \
    };\
    virtual const void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) const { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<const void*>(this); \
      else { \
	const void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME2::castTo(aTypeID); \
	if(result) return result; \
        return BASE_NAME3::castTo(aTypeID); \
      }; \
    };

/// This MACRO creates the static elements for the current class to be added to the global type registry (it is guaranteed to be added if the class is instantiated).
#define RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3,BASE_NAME4) \
    static boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getStaticObjectType() { \
      return ReaK::rtti::register_type< CLASS_NAME , CLASS_VERSION , ReaK::rtti::detail::base_type_list< BASE_NAME1, ReaK::rtti::detail::base_type_list< BASE_NAME2 , ReaK::rtti::detail::base_type_list< BASE_NAME3 , ReaK::rtti::detail::base_type_list< BASE_NAME4 > > > > >::impl.ptr; \
    }; \
    virtual boost::shared_ptr<ReaK::rtti::so_type> RK_CALL getObjectType() const { \
      return CLASS_NAME::getStaticObjectType(); \
    };\
    virtual void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<void*>(this); \
      else { \
	void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME2::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME3::castTo(aTypeID); \
	if(result) return result; \
        return BASE_NAME4::castTo(aTypeID); \
      }; \
    };\
    virtual const void* RK_CALL castTo(const boost::shared_ptr<ReaK::rtti::so_type>& aTypeID) const { \
      if(*aTypeID == *(CLASS_NAME::getStaticObjectType())) \
	return reinterpret_cast<const void*>(this); \
      else { \
	const void* result = BASE_NAME1::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME2::castTo(aTypeID); \
	if(result) return result; \
	result = BASE_NAME3::castTo(aTypeID); \
	if(result) return result; \
        return BASE_NAME4::castTo(aTypeID); \
      }; \
    };
    
/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_CONCRETE_0BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME) \
    RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME,CLASS_VERSION)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_ABSTRACT_0BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME)  \
    RK_RTTI_REGISTER_CUSTOM_FACTORY(0) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_0BASE(CLASS_NAME,CLASS_VERSION)


/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_CONCRETE_1BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME) \
    RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_ABSTRACT_1BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME)  \
    RK_RTTI_REGISTER_CUSTOM_FACTORY(0) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_1BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME)


/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_CONCRETE_2BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2) \
    RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_ABSTRACT_2BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2)  \
    RK_RTTI_REGISTER_CUSTOM_FACTORY(0) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_2BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2)
    
    
/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_CONCRETE_3BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2,BASE_NAME3) \
    RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_ABSTRACT_3BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2,BASE_NAME3)  \
    RK_RTTI_REGISTER_CUSTOM_FACTORY(0) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_3BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3)
    
    
/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_CONCRETE_4BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2,BASE_NAME3,BASE_NAME4) \
    RK_RTTI_MAKE_DEFAULT_FACTORY(CLASS_NAME) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3,BASE_NAME4)

/** This is a macro to setup the casting and ReaK::rtti registry related methods in a descendant abstract class. For the class declaration section (public), and all classes derived from typed_object are required to have this macro (or another version of it).*/
#define RK_RTTI_MAKE_ABSTRACT_4BASE(CLASS_NAME,CLASS_ID,CLASS_VERSION,CLASS_STR_NAME,BASE_NAME1,BASE_NAME2,BASE_NAME3,BASE_NAME4)  \
    RK_RTTI_REGISTER_CUSTOM_FACTORY(0) \
    RK_RTTI_REGISTER_CLASS_ID(CLASS_STR_NAME,CLASS_ID) \
    RK_RTTI_REGISTER_CLASS_4BASE(CLASS_NAME,CLASS_VERSION,BASE_NAME1,BASE_NAME2,BASE_NAME3,BASE_NAME4)








};

};

#endif



