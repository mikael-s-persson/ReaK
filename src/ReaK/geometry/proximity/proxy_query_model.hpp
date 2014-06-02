/**
 * \file proxy_query_model.hpp
 *
 * This library declares a class that can act as a collection of geometries used to represent a model used for 
 * proximity queries with other models.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
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

#ifndef REAK_PROXY_QUERY_MODEL_HPP
#define REAK_PROXY_QUERY_MODEL_HPP

#include <ReaK/geometry/shapes/shape_2D.hpp>
#include <ReaK/geometry/shapes/shape_3D.hpp>
#include "proximity_finder_2D.hpp"
#include "proximity_finder_3D.hpp"

#include <vector>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class defines a proximity-query model for 2D shapes. */
class proxy_query_model_2D : public named_object {
  public:
    
    std::vector< shared_ptr< shape_2D > > mShapeList;
    
    /**
     * Default constructor.
     */
    proxy_query_model_2D(const std::string& aName = "") : named_object(), mShapeList() { this->setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~proxy_query_model_2D() { };
    
    proxy_query_model_2D& addShape(const shared_ptr< shape_2D >& aShape) {
      mShapeList.push_back(aShape);
      return *this;
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mShapeList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mShapeList);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(proxy_query_model_2D,0xC320001A,1,"proxy_query_model_2D",named_object)
    
};



/** This class defines a proximity-query pair for 2D models. */
class proxy_query_pair_2D : public named_object {
  protected:
    
    shared_ptr< proxy_query_model_2D > mModel1;
    shared_ptr< proxy_query_model_2D > mModel2;
    
    std::vector< shared_ptr< proximity_finder_2D > > mProxFinders;
    
    void createProxFinderList();
    
  public:
    
    void setModelPair(const shared_ptr< proxy_query_model_2D >& aModel1, const shared_ptr< proxy_query_model_2D >& aModel2) {
      mModel1 = aModel1;
      mModel2 = aModel2;
      createProxFinderList();
    };
    
    /**
     * Default constructor.
     */
    proxy_query_pair_2D(const std::string& aName = "",
                        const shared_ptr< proxy_query_model_2D >& aModel1 = shared_ptr< proxy_query_model_2D >(), 
                        const shared_ptr< proxy_query_model_2D >& aModel2 = shared_ptr< proxy_query_model_2D >()) : 
                        named_object(), mModel1(aModel1), mModel2(aModel2), mProxFinders() { 
      this->setName(aName); 
      createProxFinderList();
    };
    
    /**
     * Default destructor.
     */
    virtual ~proxy_query_pair_2D() { };
    
    virtual shared_ptr< proximity_finder_2D > findMinimumDistance() const;
    
    virtual bool gatherCollisionPoints(std::vector< proximity_record_2D >& aOutput) const;
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mModel1)
        & RK_SERIAL_SAVE_WITH_NAME(mModel2)
        & RK_SERIAL_SAVE_WITH_NAME(mProxFinders);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mModel1)
        & RK_SERIAL_LOAD_WITH_NAME(mModel2)
        & RK_SERIAL_LOAD_WITH_NAME(mProxFinders);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(proxy_query_pair_2D,0xC320001C,1,"proxy_query_pair_2D",named_object)
    
};


/** This class defines a colored model for 3D geometries. */
class proxy_query_model_3D : public named_object {
  public:
    
    std::vector< shared_ptr< shape_3D > > mShapeList;
    
    /**
     * Default constructor.
     */
    proxy_query_model_3D(const std::string& aName = "") : named_object(), mShapeList() { setName(aName); };
    
    /**
     * Default destructor.
     */
    virtual ~proxy_query_model_3D() { };
    
    proxy_query_model_3D& addShape(const shared_ptr< shape_3D >& aShape) {
      mShapeList.push_back(aShape);
      return *this;
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mShapeList);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mShapeList);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(proxy_query_model_3D,0xC320001B,1,"proxy_query_model_3D",named_object)
    
};



/** This class defines a proximity-query pair for 2D models. */
class proxy_query_pair_3D : public named_object {
  protected:
    
    shared_ptr< proxy_query_model_3D > mModel1;
    shared_ptr< proxy_query_model_3D > mModel2;
    
    std::vector< shared_ptr< proximity_finder_3D > > mProxFinders;
    
    void createProxFinderList();
    
  public:
    
    void setModelPair(const shared_ptr< proxy_query_model_3D >& aModel1, const shared_ptr< proxy_query_model_3D >& aModel2) {
      mModel1 = aModel1;
      mModel2 = aModel2;
      createProxFinderList();
    };
    
    /**
     * Default constructor.
     */
    proxy_query_pair_3D(const std::string& aName = "",
                        const shared_ptr< proxy_query_model_3D >& aModel1 = shared_ptr< proxy_query_model_3D >(), 
                        const shared_ptr< proxy_query_model_3D >& aModel2 = shared_ptr< proxy_query_model_3D >()) : 
                        named_object(), mModel1(aModel1), mModel2(aModel2), mProxFinders() { 
      this->setName(aName); 
      createProxFinderList();
    };
    
    /**
     * Default destructor.
     */
    virtual ~proxy_query_pair_3D() { };
    
    virtual shared_ptr< proximity_finder_3D > findMinimumDistance() const;
    
    virtual bool gatherCollisionPoints(std::vector< proximity_record_3D >& aOutput) const;
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mModel1)
        & RK_SERIAL_SAVE_WITH_NAME(mModel2)
        & RK_SERIAL_SAVE_WITH_NAME(mProxFinders);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mModel1)
        & RK_SERIAL_LOAD_WITH_NAME(mModel2)
        & RK_SERIAL_LOAD_WITH_NAME(mProxFinders);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(proxy_query_pair_3D,0xC320001D,1,"proxy_query_pair_3D",named_object)
    
};


};

};

#endif










