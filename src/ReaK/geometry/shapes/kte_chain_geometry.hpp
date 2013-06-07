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

#include <vector>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class defines a colored model for 2D geometries. */
class colored_model_2D : public named_object {
  public:
    struct element {
      color mColor;
      shared_ptr< geometry_2D > mGeom;
      element(const color& aColor = color(), 
              const shared_ptr< geometry_2D >& aGeom = shared_ptr< geometry_2D >()) : 
              mColor(aColor), mGeom(aGeom) { };
    };
    
    std::vector< shared_ptr< pose_2D<double> > > mAnchorList;
    std::vector< element > mGeomList;
    
    
    /**
     * Default constructor.
     */
    colored_model_2D(const std::string& aName = "") : named_object(), mAnchorList(), mGeomList() { this->setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~colored_model_2D() { };
    
    
    colored_model_2D& addAnchor(const shared_ptr< pose_2D<double> >& aAnchor) {
      mAnchorList.push_back(aAnchor);
      return *this;
    };
    
    colored_model_2D& addElement(const color& aColor, const shared_ptr< geometry_2D >& aGeom) {
      mGeomList.push_back(element(aColor, aGeom));
      return *this;
    };
    
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchorList);
      A & RK_SERIAL_SAVE_WITH_ALIAS("GeomCount",mGeomList.size());
      for(std::size_t i = 0; i < mGeomList.size(); ++i) {
        { std::stringstream s_stream;
        s_stream << "color[" << i << "]";
        A & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), mGeomList[i].mColor);
        };
        
        { std::stringstream s_stream;
        s_stream << "geom[" << i << "]";
        A & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), mGeomList[i].mGeom);
        };
      };
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchorList);
      std::size_t geom_count = 0;
      A & RK_SERIAL_LOAD_WITH_ALIAS("GeomCount",geom_count);
      mGeomList.resize(geom_count);
      for(std::size_t i = 0; i < mGeomList.size(); ++i) {
        { std::stringstream s_stream;
        s_stream << "color[" << i << "]";
        A & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), mGeomList[i].mColor);
        };
        
        { std::stringstream s_stream;
        s_stream << "geom[" << i << "]";
        A & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), mGeomList[i].mGeom);
        };
      };
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(colored_model_2D,0xC3100020,1,"colored_model_2D",named_object)
    
    
};



/** This class defines a colored model for 3D geometries. */
class colored_model_3D : public named_object {
  public:
    struct element {
      color mColor;
      shared_ptr< geometry_3D > mGeom;
      element(const color& aColor = color(), 
              const shared_ptr< geometry_3D >& aGeom = shared_ptr< geometry_3D >()) : 
              mColor(aColor), mGeom(aGeom) { };
    };
    
    std::vector< shared_ptr< pose_3D<double> > > mAnchorList;
    std::vector< element > mGeomList;
    
    
    /**
     * Default constructor.
     */
    colored_model_3D(const std::string& aName = "") : named_object(), mAnchorList(), mGeomList() { setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~colored_model_3D() { };
    
    
    colored_model_3D& addAnchor(const shared_ptr< pose_3D<double> >& aAnchor) {
      mAnchorList.push_back(aAnchor);
      return *this;
    };
    
    colored_model_3D& addElement(const color& aColor, const shared_ptr< geometry_3D >& aGeom) {
      mGeomList.push_back(element(aColor, aGeom));
      return *this;
    };
    
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchorList);
      A & RK_SERIAL_SAVE_WITH_ALIAS("GeomCount",mGeomList.size());
      for(std::size_t i = 0; i < mGeomList.size(); ++i) {
        { std::stringstream s_stream;
        s_stream << "color[" << i << "]";
        A & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), mGeomList[i].mColor);
        };
        
        { std::stringstream s_stream;
        s_stream << "geom[" << i << "]";
        A & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), mGeomList[i].mGeom);
        };
      };
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchorList);
      std::size_t geom_count = 0;
      A & RK_SERIAL_LOAD_WITH_ALIAS("GeomCount",geom_count);
      mGeomList.resize(geom_count);
      for(std::size_t i = 0; i < mGeomList.size(); ++i) {
        { std::stringstream s_stream;
        s_stream << "color[" << i << "]";
        A & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), mGeomList[i].mColor);
        };
        
        { std::stringstream s_stream;
        s_stream << "geom[" << i << "]";
        A & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), mGeomList[i].mGeom);
        };
      };
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(colored_model_3D,0xC3100021,1,"colored_model_3D",named_object)
    
    
};

};

};

#endif










