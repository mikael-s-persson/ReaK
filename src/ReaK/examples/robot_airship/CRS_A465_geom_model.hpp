/**
 * \file CRS_A465_geom_model.hpp
 *
 * This library has constructs of the KTE-based kinematics and geometric models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2012
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


#ifndef RK_CRS_A465_GEOM_MODEL_HPP
#define RK_CRS_A465_GEOM_MODEL_HPP

#include "CRS_A465_models.hpp"

#include <string>

namespace ReaK {

namespace geom { 
  class box;
  class capped_cylinder;
  class colored_model_3D;
  class coord_arrows_3D;
  class proxy_query_model_3D;
  class sphere;
}

namespace robot_airship {


/**
 * This class serves to store (load / save) the data that models the CRS A465 robot and its geometry.
 */
class CRS_A465_geom_builder : public CRS_A465_model_builder {
  protected:
    
    void load_geom_from_archive(serialization::iarchive& aInput);
    
    void save_geom_to_archive(serialization::oarchive& aOutput) const;
    
  public:
    
    shared_ptr< geom::colored_model_3D > geom_model;
    shared_ptr< geom::proxy_query_model_3D > proxy_model;
    
    shared_ptr< geom::coord_arrows_3D > global_frame_arrows;
    shared_ptr< geom::coord_arrows_3D > robot_base_arrows;
    shared_ptr< geom::coord_arrows_3D > track_joint_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_1_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_2_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_3_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_4_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_5_arrows;
    shared_ptr< geom::coord_arrows_3D > arm_joint_6_arrows;
    
    shared_ptr< geom::capped_cylinder > link1_cyl;
    shared_ptr< geom::capped_cylinder > joint2_cyl;
    shared_ptr< geom::capped_cylinder > link2_cyl;
    shared_ptr< geom::capped_cylinder > link3_cyl;
    shared_ptr< geom::capped_cylinder > link5_cyl;
    shared_ptr< geom::sphere > EE_sphere;
    
    shared_ptr< geom::box > EE_bumblebee;
    shared_ptr< geom::box > EE_bumblebee_support;
    shared_ptr< geom::box > EE_gripper_box;
    shared_ptr< geom::box > EE_gripper_fingers;
    
    
    
    /**
     * Default constructor.
     */
    CRS_A465_geom_builder() : CRS_A465_model_builder() { };
    
    /**
     * This function will load all the KTE model data and geometric data from the given file.
     * \param aFileName The absolute or relative path to the model file to load (xml format only).
     */
    void load_kte_and_geom(const std::string& aFileName);
    
    /**
     * This function will create all the model data from a preset definition of the model.
     * This function may be useful to either create the preset CRS A465 model or to create a 
     * complete model that can be saved to an xml file and modified subsequently.
     */
    void create_geom_from_preset();
    
    /**
     * This function will save all the KTE model data to the given file.
     * \param aFileName The absolute or relative path to the model file to which to save (xml format only).
     */
    void save_kte_and_geom(const std::string& aFileName) const;
    
    /**
     * This function returns the geometric model that represents the CRS A465 robot in a rendered scene.
     * \return The geometric model that represents the CRS A465 robot in a rendered scene.
     */
    shared_ptr< geom::colored_model_3D > get_geometric_model() const;
    
    /**
     * This function returns the proximity model that represents the CRS A465 robot geometry for proximity 
     * query purposes (collision checking).
     * \return The proximity model that represents the CRS A465 robot geometry for proximity query purposes (collision checking).
     */
    shared_ptr< geom::proxy_query_model_3D > get_proximity_model() const;
    
};


};
  
  
};


#endif






