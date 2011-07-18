/**
 *\file shared_object_base.hpp
 *
 * This library declares an underlying class that is the base for the ReaK::shared_object class.
 * This header and its contents are not meant to be included are directly used by end users 
 * except for the null_deleter and scoped_deleter for use in setting up shared pointers with
 * either no deletion on destruction of the last shared pointer or a deletion that is in the 
 * same module scope as the code that created the shared pointer (this is the most usual,
 * although for a single module project it has no real effect).
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

#ifndef SHARED_OBJECT_BASE_HPP
#define SHARED_OBJECT_BASE_HPP

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "defs.hpp"


namespace ReaK {

/** 
 * This structure is a simple callable structure that does nothing. It acts as a place-holder for
 * a special deleter for a shared_ptr, in this case this "null deleter" simply does not delete anything.
 * This is useful for an object for which you want a shared_ptr to, but that is not going to be
 * deleted by this shared_ptr branch (i.e. will be deleted by another shared_ptr branch). This can
 * be used to break cycles in the object hierarchy.
 */
struct null_deleter
{
    void operator()(void const *) const
    {
    }
};

/**
 * This is a base class for the ReaK::shared_object. This holds a "null deleted" shared_ptr to
 * itself as well as a basic virtual destroy function that is used to ensure proper scoping of the
 * object deletion. Those two features allow the descendant-class objects to deliver a weak_ptr
 * to themselves (that will expire when the object is destroyed) and to be shared across modules 
 * without deleter module scope issues (i.e. it will be deleted from the same compiled library 
 * from which is got created with its own vtable).
 */
class shared_object_base {
  protected:
    boost::shared_ptr<shared_object_base> mThis;
  
  public:
    virtual void RK_CALL destroy() = 0;
    
    /**
     * This method returns a weak_ptr to this object. The weak pointer will expire as the object gets
     * deleted.
     * \note This pointer is weaker than a regular weak pointer because locking this pointer does not
     *       guarantee that it remains for as long as the locked pointer exists so this is more intended
     *       for use as a test pointer to see if the object still exists. If a real shared ownership 
     *       scheme is desired, the shared pointer should be obtained from an owner of this object.
     * \return weak pointer to this object. A lock on this pointer will not give shared ownership and 
     *         should be regarded as a momentary access to the object with no guarantee that it will 
     *         not get deleted in the meantime. A real shared ownership can only be obtained from the 
     *         actual owner of this object.
     */
    boost::weak_ptr<shared_object_base> RK_CALL getWeakPtr() const { return mThis; };
    
    shared_object_base() { mThis = boost::shared_ptr<shared_object_base>(this,null_deleter()); };

    virtual ~shared_object_base() { RK_NOTICE(8,"Shared object base destructor reached!");};
};

/** 
 * This structure is a simple callable structure that deletes an object by calling the virtual "destroy" 
 * method. It acts as a special deleter for a shared_ptr, in this case this "scoped deleter" simply makes
 * sure the object is deleted via its vtable and thus will go to the code in the same executable scope
 * from which it was created and thus making the shared_ptr movable between executable modules.
 */
struct scoped_deleter {
  void operator()(shared_object_base * p) const {
    RK_NOTICE(8,"Shared object base scoped deleter reached!");
    p->destroy();
  };
};

};

#endif
