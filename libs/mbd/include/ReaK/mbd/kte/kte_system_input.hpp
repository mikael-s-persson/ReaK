/**
 * \file kte_system_input.hpp
 * 
 * This library defines the base class for system inputs to a KTE model. A system input 
 * is simply a vector of values which serve as an input to a KTE model. This model is useful
 * when using a KTE model into a state-space system definition.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_KTE_SYSTEM_INPUT_HPP
#define REAK_KTE_SYSTEM_INPUT_HPP

#include <ReaK/core/base/named_object.hpp>

#include <vector>


namespace ReaK {

namespace kte {
  
/**
 * This class is a base class for system inputs to a KTE model. A system input 
 * is simply a vector of values which serve as an input to a KTE model. This model is useful
 * when using a KTE model into a state-space system definition.
 */
class system_input : public virtual named_object {
  public:
    
    /**
     * Constructs a system input class with the given name.
     */
    system_input(const std::string& aName = "") {
      this->setName(aName);
    };
    
    /**
     * Destructor.
     */
    virtual ~system_input() { };
    
    /**
     * Returns the number of input variables provided by this system input.
     * \return the number of input variables provided by this system input.
     */
    virtual unsigned int getInputCount() const = 0;
    
    /**
     * Returns the input variable at index i, with read-write access.
     * \param i The index of the input variable.
     * \return The variable at index i.
     */
    virtual void setInput(unsigned int i, double val) = 0;
    /**
     * Returns the input variable at index i, with read-only access.
     * \param i The index of the input variable.
     * \return The variable at index i.
     */
    virtual double getInput(unsigned int i) const = 0;
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(system_input,0xC2100033,1,"system_input",named_object)
    
};


};

};

#endif













