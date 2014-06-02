/**
 * \file X8_quadrotor_geom.hpp
 *
 * This library has constructs of the geometric models (visualization and proxy) for the X8 quadrotor.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
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


#ifndef RK_X8_QUADROTOR_GEOM_HPP
#define RK_X8_QUADROTOR_GEOM_HPP

#include <ReaK/ctrl/kte_models/uav_kinematics.hpp>

#include <string>

namespace ReaK {

namespace geom {

class colored_model_3D;
class proxy_query_model_3D;


/**
 * This class serves to store (load / save) the data that models the CRS A465 robot and its geometry.
 */
class X8_quadrotor_geom : public named_object {
  public:
    
    shared_ptr< colored_model_3D > geom_model;
    shared_ptr< proxy_query_model_3D > proxy_model;
    
    
    /**
     * Default constructor.
     */
    X8_quadrotor_geom() : named_object() { setName("X8_quadrotor_geom"); };
    
    virtual ~X8_quadrotor_geom() { };
    
    /**
     * This function will create all the model data from a preset definition of the model.
     * This function may be useful to either create the preset X8 quadrotor model or to create a 
     * complete model that can be saved to an xml file and modified subsequently.
     */
    void create_geom_from_preset(const kte::UAV_kinematics& aModel);
    
    /**
     * This function returns the geometric model that represents the CRS A465 robot in a rendered scene.
     * \return The geometric model that represents the CRS A465 robot in a rendered scene.
     */
    shared_ptr< colored_model_3D > get_geometric_model() const { return geom_model; };
    
    /**
     * This function returns the proximity model that represents the CRS A465 robot geometry for proximity 
     * query purposes (collision checking).
     * \return The proximity model that represents the CRS A465 robot geometry for proximity query purposes (collision checking).
     */
    shared_ptr< proxy_query_model_3D > get_proximity_model() const { return proxy_model; };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(X8_quadrotor_geom,0xC3300001,1,"X8_quadrotor_geom",named_object)
    
};


};
  
  
};


#endif






