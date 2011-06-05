/**
 * \file named_object.hpp
 *
 * This library defines the ReaK::named_object class and its interface. This serves as the
 * basic class for all named objects in the ReaK platform. It can also serve as an example
 * of the simplest full-fledged class that implements the RTTI and serialization interfaces.
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

#ifndef NAMED_OBJECT_HPP
#define NAMED_OBJECT_HPP

#include <string.h>
#include "shared_object.hpp"
#include "rtti/so_type_repo.hpp"

namespace ReaK {

/**
 * This class declares an interface for an object that has a name.
 */
class named_interface {
  public:
    virtual ~named_interface() { };

    /**
     * This method returns the name of the object.
     * \pre The name of the object is known to the object (privately).
     * \post The name of the object is given to the caller as constant.
     * \return The name of the object.
     */
    virtual const std::string& getName() const = 0;
    /**
     * This method sets the name of the object.
     * \pre Any state.
     * \post The name of the object is known to the object, or by default, this method has no effect, i.e., the name is read-only.
     * \param aName The new name of the object.
     */
    virtual void setName(const std::string& aName) = 0;
};

/**
 * This class declares a base class for an object that stores a name string (read- and writable).
 */
class named_object : public virtual shared_object , public named_interface {
  protected:
    std::string mName;
  public:
    /**
     * Default Constructor. The name must be set by the derived class using the setName function.
     */
    named_object() { };

    virtual ~named_object() {
      RK_NOTICE(4,"Object " << mName << " deleted.");
    };

    /**
     * This method returns the name of the object.
     * \pre The name of the object is known to the object (privately stored in mName).
     * \post The name of the object is given to the caller as constant.
     * \return The name of the object.
     */
    virtual const std::string& getName() const {
      return mName;
    };
    /**
     * This method sets the name of the object.
     * \pre Any state, mName is empty or not.
     * \post The name of the object is written to the data member mName.
     * \param aName The new name of the object.
     */
    virtual void setName(const std::string& aName) {
      mName = aName;
    };

    /**
     * This method saves the content of the object to a serial archive of any type.
     *
     * \pre any valid object state.
     * \post the object state has been saved to the archive such that it can be completely reconstructed at load-time.
     * \param A any type of output archive.
     * \param Version the version of this object that is to be saved (always the latest version).
     */
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & std::pair<std::string, const std::string&>("name",mName);
    };

    /**
     * This method loads the content of the object from a serial archive of any type.
     *
     * \pre any object state.
     * \post the object state has been loaded from the archive such that it is completely reconstructed and valid.
     * \param A any type of input archive.
     * \param Version the version of this object that was saved (it is the user's responsability to maintain backward compatibility as much as desired).
     */
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & std::pair<std::string, std::string&>("name",mName);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(named_object,0xC0000000,1,"named_object",shared_object)

};




};

#endif


