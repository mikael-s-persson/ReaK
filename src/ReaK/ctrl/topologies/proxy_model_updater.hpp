/**
 * \file proxy_model_updater.hpp
 * 
 * This library defines an abstract base class for to "update" a proximity-query model for a given time (for 
 * dynamic models). A typical derived class implementation would be to query a spatial-trajectory (see SpatialTrajectoryConcept)
 * for the state at the given time and then apply that state to the geometry used by the proximity-query method.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PROXY_MODEL_UPDATER_HPP
#define REAK_PROXY_MODEL_UPDATER_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"


namespace ReaK {

namespace pp {


/**
 * This base-class is used to update dynamic proximity-query models to a given time value.
 */
class proxy_model_updater : public shared_object {
  public:
    typedef proxy_model_updater self;
    
    proxy_model_updater() { };
    
    virtual void synchronize_proxy_model(double t) const = 0;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400029,1,"proxy_model_updater",shared_object)
    
    
};



};

};

#endif

