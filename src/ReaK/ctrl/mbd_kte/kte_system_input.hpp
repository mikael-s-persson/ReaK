
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

#ifndef KTE_SYSTEM_INPUT_HPP
#define KTE_SYSTEM_INPUT_HPP

#include "base/named_object.hpp"

#include <vector>


namespace ReaK {

namespace kte {


class system_input : public virtual named_object {
  public:
    
    system_input(const std::string& aName = "") {
      this->setName(aName);
    };
    
    virtual ~system_input() { };
    
    virtual unsigned int getInputCount() const = 0;
    
    virtual double& getInput(unsigned int i) = 0;
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













