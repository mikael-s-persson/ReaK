/**
 * \file kte_chain_geometry.hpp
 *
 * This library declares a class 
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_KTE_CHAIN_GEOMETRY_HPP
#define REAK_KTE_CHAIN_GEOMETRY_HPP

#include "colored_model.hpp"
#include "proximity/proxy_query_model.hpp"
#include "mbd_kte/kte_map_chain.hpp"

#include <string>
#include <map>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class defines a colored model for 2D geometries. */
class kte_chain_geometry_2D : public named_object {
  public:
    std::map< std::string, colored_geometry_2D > mGeomList;
    
    /**
     * Default constructor.
     */
    kte_chain_geometry_2D(const std::string& aName = "") : named_object(), mGeomList() { this->setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~kte_chain_geometry_2D() { };
    
    kte_chain_geometry_2D& addElement(const std::string& aKTEObjName, const color& aColor, const shared_ptr< geometry_2D >& aGeom) {
      mGeomList[aKTEObjName] = colored_geometry_2D(aColor, aGeom);
      return *this;
    };
    
    shared_ptr< colored_model_2D > attachGeomToKTEChain(const kte::kte_map_chain& aKTEChain) const;
    
    shared_ptr< proxy_query_model_2D > attachProxyToKTEChain(const kte::kte_map_chain& aKTEChain) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mGeomList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mGeomList);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(kte_chain_geometry_2D,0xC3100024,1,"kte_chain_geometry_2D",named_object)
};



/** This class defines a colored model for 3D geometries. */
class kte_chain_geometry_3D : public named_object {
  public:
    std::map< std::string, colored_geometry_3D > mGeomList;
    
    /**
     * Default constructor.
     */
    kte_chain_geometry_3D(const std::string& aName = "") : named_object(), mGeomList() { this->setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~kte_chain_geometry_3D() { };
    
    kte_chain_geometry_3D& addElement(const std::string& aKTEObjName, const color& aColor, const shared_ptr< geometry_3D >& aGeom) {
      mGeomList[aKTEObjName] = colored_geometry_3D(aColor, aGeom);
      return *this;
    };
    
    shared_ptr< colored_model_3D > attachGeomToKTEChain(const kte::kte_map_chain& aKTEChain) const;
    
    shared_ptr< proxy_query_model_3D > attachProxyToKTEChain(const kte::kte_map_chain& aKTEChain) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mGeomList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mGeomList);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(kte_chain_geometry_3D,0xC3100025,1,"kte_chain_geometry_3D",named_object)
};

};

};

#endif










