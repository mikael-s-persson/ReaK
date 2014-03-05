/**
 * \file X8_run_planner.cpp
 *
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */


#include "kte_models/uav_kinematics.hpp"
#include "X8_quadrotor_geom.hpp"
#include "IHAQR_topology.hpp"
#include "MEAQR_topology.hpp"
#include "ss_systems/quadrotor_system.hpp"
#include "topologies/se3_random_samplers.hpp"

#include "MEAQR_rrtstar_planner.hpp"
#include "MEAQR_sbastar_planner.hpp"

#include "proximity/proxy_query_model.hpp"

#include "path_planning/sbmp_point_recorder.hpp"
#include "path_planning/basic_sbmp_reporters.hpp"

#include "shapes/colored_model.hpp"
#include "shapes/coord_arrows_3D.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"
#include "optimization/optim_exceptions.hpp"

using namespace ReaK;
using namespace pp;
using namespace ctrl;


typedef IHAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_IHAQR_space_type;
typedef MEAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_MEAQR_space_type;
typedef IHAQR_topology_with_CD< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_IHAQRCD_space_type;
typedef MEAQR_topology_with_CD< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_MEAQRCD_space_type;


int main(int argc, char ** argv) {
  using namespace ReaK;
  
  if(argc < 5) {
    std::cout << "Usage: ./X8_test_scene [world_filename.xml] [quadrotor_spaces.xml] [planner_parameters.xml]" << std::endl
              << "\t world_filename.xml: \tThe XML filename of the world geometry (including start and end position of planning)." << std::endl
              << "\t quadrotor_spaces.xml: \tThe XML filename of the definitions of the quad-rotor IHAQR and MEAQR topologies." << std::endl
              << "\t planner_params.xml: \tThe XML filename of the set of parameters for the planner." << std::endl
              << "\t path/to/output_records: \tThe file-path to the results of the planner." << std::endl;
    return 1;
  };
  
  shared_ptr<kte::UAV_kinematics> builder   = shared_ptr< kte::UAV_kinematics >(new kte::UAV_kinematics());
  shared_ptr< kte::kte_map_chain >         kin_chain = builder->getKTEChain();
  
  geom::X8_quadrotor_geom geom_builder;
  geom_builder.create_geom_from_preset(*builder);
  
  shared_ptr< frame_3D<double> > target_frame = builder->getFrame3D(0);
  
  shared_ptr< geom::colored_model_3D >     world_geom_model;
  shared_ptr< geom::proxy_query_model_3D > world_geom_proxy;
  vect<double,3>  start_position;
  vect<double,3>  end_position;
  {
    std::string world_filename(argv[1]);
    serialization::xml_iarchive in(world_filename);
    in >> world_geom_model >> world_geom_proxy >> start_position >> end_position;
  };
  
  shared_ptr< quadrotor_system >      X8_sys;
  shared_ptr< X8_IHAQRCD_space_type > X8_IHAQR_space;
  shared_ptr< X8_MEAQRCD_space_type > X8_MEAQR_space; 
  {
    std::string X8spaces_filename(argv[2]);
    serialization::xml_iarchive in(X8spaces_filename);
    
    shared_ptr< X8_IHAQR_space_type > tmp_IHAQR_space;
    shared_ptr< X8_MEAQR_space_type > tmp_MEAQR_space;
    
    in >> X8_sys >> tmp_IHAQR_space >> tmp_MEAQR_space;
    
    X8_IHAQR_space = shared_ptr< X8_IHAQRCD_space_type >(new X8_IHAQRCD_space_type(
      *tmp_IHAQR_space, 
      builder));
    X8_MEAQR_space = shared_ptr< X8_MEAQRCD_space_type >(new X8_MEAQRCD_space_type(
      *tmp_MEAQR_space, 
      builder));
  };
  
  shared_ptr< geom::proxy_query_pair_3D > X8_to_world_proxy(new geom::proxy_query_pair_3D(
    "X8_to_world_proxy", 
    geom_builder.get_proximity_model(),
    world_geom_proxy));
  
  X8_IHAQR_space->m_proxy_env_3D.push_back(X8_to_world_proxy);
  X8_MEAQR_space->m_proxy_env_3D.push_back(X8_to_world_proxy);
  
  std::string planner_kind  = "rrt_star";
  std::size_t max_vertices  = 1000;
  std::size_t prog_interval = 100;
  std::size_t max_results   = 50;
  std::size_t knn_method    = ReaK::pp::LINEAR_SEARCH_KNN;  // (0) or ReaK::pp::DVP_BF2_TREE_KNN (3)
  {
    std::string planner_param_filename(argv[3]);
    serialization::xml_iarchive in(planner_param_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("planner_kind", planner_kind)
       & RK_SERIAL_LOAD_WITH_ALIAS("max_vertices", max_vertices) 
       & RK_SERIAL_LOAD_WITH_ALIAS("prog_interval", prog_interval) 
       & RK_SERIAL_LOAD_WITH_ALIAS("max_results", max_results) 
       & RK_SERIAL_LOAD_WITH_ALIAS("knn_method", knn_method);
  };
  
  std::string output_path(argv[4]);
  
  
  
  try {
    quadrotor_system::state_space_type::point_type start_state, goal_state;
    set_position(start_state, start_position);
    set_position(goal_state, end_position);
    
    X8_MEAQRCD_space_type::point_type start_point( start_state );
    X8_MEAQRCD_space_type::point_type goal_point( goal_state );
    
//     typedef sbmp_point_recorder<
//       quadrotor_system::state_space_type,
//       MEAQR_to_state_mapper,
//       print_sbmp_progress<> > point_recorder_type;
    
    std::ofstream timer_output(output_path + "/X8_pp_times.txt");
    
    any_sbmp_reporter_chain< X8_MEAQRCD_space_type > report_chain;
    report_chain.add_reporter( timing_sbmp_report<>(timer_output) );
    report_chain.add_reporter( sbmp_point_recorder< quadrotor_system::state_space_type, MEAQR_to_state_mapper >(
      shared_ptr<quadrotor_system::state_space_type>(X8_MEAQR_space,&(X8_MEAQR_space->get_state_space())),
      MEAQR_to_state_mapper(),
      output_path + "/X8_pp_") );
    report_chain.add_reporter( print_sbmp_progress<>() );
    
    path_planning_p2p_query< X8_MEAQRCD_space_type > pp_query(
      "pp_query",
      X8_MEAQR_space,
      start_point,
      goal_point,
      max_results);
    
    typedef MEAQR_rrtstar_planner< 
      quadrotor_system::state_space_type, 
      quadrotor_system, 
      position_only_sampler > X8_rrtstar_planner_type;
    X8_rrtstar_planner_type X8_planner(
      X8_MEAQR_space, max_vertices, prog_interval, 
      ADJ_LIST_MOTION_GRAPH | knn_method,
      UNIDIRECTIONAL_PLANNING,
      0.1, 0.05, 3, report_chain);
    
//     typedef MEAQR_sbastar_planner< 
//       quadrotor_system::state_space_type, 
//       quadrotor_system, 
//       position_only_sampler > X8_sbastar_planner_type;
//     X8_sbastar_planner_type X8_planner(
//       X8_MEAQR_space, max_vertices, prog_interval,
//       ADJ_LIST_MOTION_GRAPH | knn_method,
//       LAZY_COLLISION_CHECKING | PLAN_WITH_VORONOI_PULL,
//       0.1, 0.05, 
//       0.25 * X8_MEAQR_space->get_max_time_horizon() * X8_MEAQR_space->get_idle_power_cost(start_point), 
//       3, report_chain);
//     X8_planner.set_initial_density_threshold(0.0);
//     X8_planner.set_initial_relaxation(10.0);
//     X8_planner.set_initial_SA_temperature(5.0);
    
    
    pp_query.reset_solution_records();
    X8_planner.solve_planning_query(pp_query);
    
  } catch( ReaK::optim::infeasible_problem& e) {
    std::cout << "ERROR: Planning was unsuccessful due to some infeasibility problem: " << e.what() << std::endl;
    return 1;
  };
  
  return 0;
};







