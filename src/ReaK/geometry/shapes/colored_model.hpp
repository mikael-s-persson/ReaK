/**
 * \file colored_model.hpp
 *
 * This library declares a class that can act as a collection of geometries with associated colors and "external" anchors.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date September 2012
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

#ifndef REAK_COLORED_MODEL_HPP
#define REAK_COLORED_MODEL_HPP

#include "geometry_2D.hpp"
#include "geometry_3D.hpp"
#include "color.hpp"

#include <vector>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {

  
class colored_geometry_2D : public shared_object {
  public:
    color mColor;
    shared_ptr< geometry_2D > mGeom;
    
    colored_geometry_2D(const color& aColor = color(),
                        const shared_ptr< geometry_2D >& aGeom = shared_ptr< geometry_2D >()) :
                        mColor(aColor), mGeom(aGeom) { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mColor)
        & RK_SERIAL_SAVE_WITH_NAME(mGeom);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mColor)
        & RK_SERIAL_LOAD_WITH_NAME(mGeom);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(colored_geometry_2D,0xC3100022,1,"colored_geometry_2D",shared_object)
    
};
  
  
/** This class defines a colored model for 2D geometries. */
class colored_model_2D : public named_object {
  public:
    std::vector< shared_ptr< pose_2D<double> > > mAnchorList;
    std::vector< colored_geometry_2D > mGeomList;
    
    
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
      mGeomList.push_back(colored_geometry_2D(aColor, aGeom));
      return *this;
    };
    
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchorList)
        & RK_SERIAL_SAVE_WITH_NAME(mGeomList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchorList)
        & RK_SERIAL_LOAD_WITH_NAME(mGeomList);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(colored_model_2D,0xC3100020,1,"colored_model_2D",named_object)
    
    
};


class colored_geometry_3D : public shared_object {
  public:
    color mColor;
    shared_ptr< geometry_3D > mGeom;
    
    colored_geometry_3D(const color& aColor = color(),
                        const shared_ptr< geometry_3D >& aGeom = shared_ptr< geometry_3D >()) :
                        mColor(aColor), mGeom(aGeom) { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mColor)
        & RK_SERIAL_SAVE_WITH_NAME(mGeom);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mColor)
        & RK_SERIAL_LOAD_WITH_NAME(mGeom);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(colored_geometry_3D,0xC3100023,1,"colored_geometry_3D",shared_object)
    
};


/** This class defines a colored model for 3D geometries. */
class colored_model_3D : public named_object {
  public:
    std::vector< shared_ptr< pose_3D<double> > > mAnchorList;
    std::vector< colored_geometry_3D > mGeomList;
    
    
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
      mGeomList.push_back(colored_geometry_3D(aColor, aGeom));
      return *this;
    };
    
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchorList)
        & RK_SERIAL_SAVE_WITH_NAME(mGeomList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchorList)
        & RK_SERIAL_LOAD_WITH_NAME(mGeomList);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(colored_model_3D,0xC3100021,1,"colored_model_3D",named_object)
    
    
};

};

};

#endif










