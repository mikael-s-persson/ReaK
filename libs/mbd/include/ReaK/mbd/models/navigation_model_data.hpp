/**
 * \file navigation_model_data.hpp
 * 
 * This library defines a class that is a meant to hold a collection of modeling data related to 
 * navigation scenario (in 3D) where we have a robot model which is meant to navigate in a cluttered environment.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_UAV_NAVIGATION_MODEL_DATA_HPP
#define REAK_UAV_NAVIGATION_MODEL_DATA_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "joint_space_limits.hpp"

#include "direct_kinematics_model.hpp"

namespace ReaK {

/* forward-declarations */
namespace geom {
  class proxy_query_model_3D;
  class colored_model_3D;
  class proxy_query_pair_3D;
};

  
namespace kte {


/**
 * This class is a meant to hold a collection of modeling data related to a
 * navigation scenario (in 3D) where we have a UAV model that is meant to navigate 
 * in a cluttered environment.
 * This class exposes all the data members publicly because it is not a "functional" class
 * that exposes some sort of abstract functional interface, it is a data-collection class 
 * that just acts as a convenient repository for the data of this scenario, and also for 
 * the convenience of loading / saving the data.
 */
class navigation_scenario : public named_object {
  public:
    /** Holds a pointer to the base-frame (reference frame) of the robot's coordinate system. */
    shared_ptr< frame_3D<double> >                 robot_base_frame;
    /** Holds the robot's direct-kinematics model (polymorphic). */
    shared_ptr< direct_kinematics_model >          robot_kin_model;
    /** Holds the robot's joint limits, for velocity, acceleration and jerk. */
    shared_ptr< joint_limits_collection<double> >  robot_jt_limits;
    /** Holds the robot's proximity-query model, i.e., for checking collisions involving the robot. */
    shared_ptr< geom::proxy_query_model_3D >       robot_proxy;
    /** Holds the robot's geometric model, i.e., for displaying the robot's 3D representation. */
    shared_ptr< geom::colored_model_3D >           robot_geom_model;
    
    /** Holds the target's frame, i.e., the goal of the navigation is that the robot's end-effector frame meets this target frame. */
    shared_ptr< frame_3D<double> >           target_frame;
    
    /** Holds the environment's geometric models, i.e., for displaying, in 3D, objects of the environment (static in the robot's base-frame). */
    std::vector< shared_ptr< geom::colored_model_3D > > env_geom_models;
    /** Holds the environment's proximity-query models, i.e., for checking collisions involving objects of the environment (static in the robot's base-frame). */
    std::vector< shared_ptr< geom::proxy_query_model_3D > > env_proxy_models;
    
    /** Holds the robot vs. environment proximity-query pairs, i.e., for checking collisions between the robot and the environment. */
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > robot_env_proxies;
    
    /**
     * Default constructor. Leaves all data members empty.
     */
    navigation_scenario();
    
    /**
     * This function loads the robot models from a given filename. The file is expected to have 
     * all the models laid out in the following order: robot_base_frame, robot_kin_model, 
     * robot_jt_limits, robot_geom_model, and robot_proxy.
     * \param fileName The filename of the file from which to load the robot models.
     */
    void load_robot(const std::string& fileName);
    
    /**
     * This function saves the robot models to a given filename. The file will have 
     * all the models laid out in the following order: robot_base_frame, robot_kin_model, 
     * robot_jt_limits, robot_geom_model, and robot_proxy.
     * \param fileName The filename of the file to which to save the robot models.
     */
    void save_robot(const std::string& fileName) const;
    
    /**
     * This function loads (and adds) an environment object from a given filename. The file is 
     * expected to have the models laid out in the following order: env_geom_model and env_proxy.
     * \param fileName The filename of the file from which to load the environment object.
     */
    void load_environment(const std::string& fileName);
    
    /**
     * This function saves an environment object to a given filename. The file will 
     * have the models laid out in the following order: env_geom_model and env_proxy.
     * \param id The index of the environment object in the array of environment geometries.
     * \param fileName The filename of the file to which to save the environment object.
     */
    void save_environment(std::size_t id, const std::string& fileName) const;
    
    /**
     * This function clears the list of environment objects (and the related proximity models). 
     */
    void clear_environment();
    
  private:
    void create_robot_env_proxies();
    
  public:  
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const;
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(navigation_scenario,0xC210005D,1,"navigation_scenario",named_object)
    
    
};



};

};

#endif

