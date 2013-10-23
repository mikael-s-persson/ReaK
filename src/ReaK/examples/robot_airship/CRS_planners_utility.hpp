/**
 * \file CRS_planners_utility.hpp
 * 
 * This library defines a class
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_CRS_PLANNERS_UTILITY_HPP
#define REAK_CRS_PLANNERS_UTILITY_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "shapes/colored_model.hpp"
#include "proximity/proxy_query_model.hpp"
#include "kte_models/manip_P3R3R_arm.hpp"
#include "topologies/joint_space_limits.hpp"

namespace ReaK {
  
namespace robot_airship {


class scenario_data : public named_object {
  public:
    shared_ptr< frame_3D<double> >                    chaser_base_frame;
    shared_ptr< kte::manip_P3R3R_kinematics >         chaser_kin_model;
    shared_ptr< pp::joint_limits_collection<double> > chaser_jt_limits;
    shared_ptr< geom::proxy_query_model_3D >          chaser_proxy;
    shared_ptr< geom::colored_model_3D >              chaser_geom_model;
    vect_n<double>                                    chaser_jt_positions;
    
    shared_ptr< kte::kte_map_chain >     target_kin_chain;
    shared_ptr< frame_3D<double> >       target_state;
    shared_ptr< frame_3D<double> >       target_frame;
    shared_ptr< geom::proxy_query_model_3D > target_proxy;
    shared_ptr< geom::colored_model_3D > target_geom_model;
    
    shared_ptr< geom::proxy_query_pair_3D > chaser_target_proxy;
    
    std::vector< shared_ptr< geom::colored_model_3D > > env_geom_models;
    
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > chaser_env_proxies;
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > target_env_proxies;
    
    scenario_data(const std::string& aChaserFileName);
    
    
    void load_target(const std::string& fileName);
    void load_positions(const std::string& fileName);
    
    void load_environment(const std::string& fileName);
    void clear_environment();
    
  private:
    
    vect_n<double> get_chaser_goal_config();
    
  public:
    
};





/*

typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o0_jt_space_type CRS3D_jspace_rl_o0_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o1_jt_space_type CRS3D_jspace_rl_o1_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o2_jt_space_type CRS3D_jspace_rl_o2_type;

typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o0_jt_space_type CRS3D_jspace_o0_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o1_jt_space_type CRS3D_jspace_o1_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o2_jt_space_type CRS3D_jspace_o2_type;


typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o0_workspace_type CRS3D_workspace_o0_i1_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_i1_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i1_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::cubic_hermite_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_i3_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::cubic_hermite_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i3_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::quintic_hermite_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i5_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::svp_Ndof_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_svp_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::svp_Ndof_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_svp_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::sap_Ndof_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_sap_type;


typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o0_map_type CRS3D_rlDK_o0_type;
typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o1_map_type CRS3D_rlDK_o1_type;
typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o2_map_type CRS3D_rlDK_o2_type;

*/



};

};

#endif

