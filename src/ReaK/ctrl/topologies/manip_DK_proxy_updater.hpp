/**
 * \file manip_DK_proxy_updater.hpp
 * 
 * This library defines a class for path-planning for a manipulator moving inside an environment 
 * with obstacles and physical limits (joint limits). This class is essentially just an assembly 
 * of many of the building blocks in the ReaK path-planning library. This class can also draw the 
 * elements of the motion graph (edges) as end-effector trajectories in a Coin3D scene-graph.
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

#ifndef REAK_MANIP_FREE_WORKSPACE_HPP
#define REAK_MANIP_FREE_WORKSPACE_HPP

#include "base/defs.hpp"
#include <boost/config.hpp>

#include "proxy_model_updater.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/spatial_trajectory_concept.hpp"  // for SpatialTrajectoryConcept

#include "kte_models/direct_kinematics_model.hpp"
#include "direct_kinematics_topomap_detail.hpp"       // for write_joint_coordinates_impl


namespace ReaK {

namespace pp {




/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates (both 
 * generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename JointTrajectory>
class manip_DK_proxy_updater : public proxy_model_updater {
  public:
    
    typedef manip_DK_proxy_updater< JointTrajectory > self;
    typedef typename spatial_trajectory_traits< JointTrajectory >::topology temporal_space_type;
    typedef typename topology_traits< temporal_space_type >::point_type point_type;
    typedef typename spatial_trajectory_traits< JointTrajectory >::const_waypoint_descriptor wp_desc_type;
    
    BOOST_CONCEPT_ASSERT((SpatialTrajectoryConcept<JointTrajectory, temporal_space_type>));
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::direct_kinematics_model > model; 
  private:
    shared_ptr< JointTrajectory > traj;
    mutable std::pair<wp_desc_type, point_type> last_wp;
  public:
    
    void set_trajectory(const shared_ptr< JointTrajectory >& aTraj) {
      traj = aTraj;
      if(traj)
        last_wp = traj->get_waypoint_at_time(traj->get_start_time());
    };
    
    manip_DK_proxy_updater(const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >(),
                           const shared_ptr< JointTrajectory >& aTraj = shared_ptr< JointTrajectory >()) :
                           model(aModel), traj(aTraj) { };
    
    
    virtual void synchronize_proxy_model(double t) const {
      if(!traj)
        return;
      
      last_wp = traj->move_time_diff_from(last_wp, t - last_wp.second.time);
      
      detail::write_joint_coordinates_impl(last_wp, traj->get_temporal_space().get_space_topology(), model);
    
      model->doDirectMotion();  //NOTE: It is assumed that the motion of the proxy-model is linked to the manip-kin-model used here.
      
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(traj);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(traj);
      if(traj)
        last_wp = traj->get_waypoint_at_time(traj->get_start_time());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240002A,1,"manip_DK_proxy_updater",shared_object)
    
    
};




};

};

#endif

