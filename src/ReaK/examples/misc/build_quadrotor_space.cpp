/**
 * \file quadrotor_scene.cpp
 *
 * This application 
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */


#include "IHAQR_topology.hpp"
#include "MEAQR_topology.hpp"
#include "ss_systems/quadrotor_system.hpp"

#include "topologies/se3_random_samplers.hpp"

#include "path_planning/path_planner_options.hpp"

#include "serialization/xml_archiver.hpp"

int main(int argc, char ** argv) {
  using namespace ReaK;
  using namespace serialization;
  using namespace pp;
  using namespace ctrl;
  
  shared_ptr< quadrotor_system > quad_sys( new quadrotor_system(
    "Quadrotor_system", 
    2.025, // aMass, 
    mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(vect<double,3>(0.0613, 0.0612, 0.1115))), // aInertiaMoment
    mat<double,mat_structure::diagonal>(vect<double,3>(0.1, 0.1, 0.1)),
    mat<double,mat_structure::diagonal>(vect<double,3>(0.1, 0.1, 0.1))));
  
  
  vect<double,3> min_corner(0.0, 0.0, -5.0);  // min corner
  vect<double,3> max_corner(5.0, 5.0, 0.0);  // max corner
  double v_max = 6.0;
  double w_max = M_PI;
  vect<double,4> u_max(35.0, 5.0,  5.0,  3.0);
  mat<double,mat_structure::diagonal> weightR_mat(vect<double,4>(25, 50, 50, 50));
//   mat<double,mat_structure::diagonal> weightR_mat(vect<double,4>(0.25, 0.5, 0.5, 0.5));
  mat<double,mat_structure::diagonal> Rscale(vect<double,4>(v_max / u_max[0], w_max / u_max[1], w_max / u_max[2], w_max / u_max[3]));
  
  double char_length_sqr = norm_2_sqr(max_corner - min_corner);
  mat<double,mat_structure::diagonal> weightQ_mat(vect_n<double>(12, 1.0));
  mat<double,mat_structure::diagonal> Qscale(vect_n<double>(
    u_max[0] * v_max / char_length_sqr,
    u_max[0] * v_max / char_length_sqr,
    u_max[0] * v_max / char_length_sqr,
    u_max[0] / v_max, 
    u_max[0] / v_max, 
    u_max[0] / v_max, 
    u_max[0] * v_max, 
    u_max[0] * v_max, 
    u_max[0] * v_max, 
    u_max[0] * char_length_sqr / v_max,
    u_max[0] * char_length_sqr / v_max,
    u_max[0] * char_length_sqr / v_max
  ));
  
  typedef IHAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > IHAQR_space_type;
  
  shared_ptr< IHAQR_space_type > quad_space( new IHAQR_space_type(
      "Quadrotor_IHAQR_topology",
      quad_sys,
      make_se3_space(
        "Quadrotor_state_space",
        min_corner,  // min corner
        max_corner,  // mix corner
        v_max,    // aMaxSpeed
        w_max),  // aMaxAngularSpeed
//       vect<double,4>(0.0, -5.0, -5.0, -3.0),  // aMinInput
      vect<double,4>(0.0, -5.0, -5.0, -5.0),  // aMinInput
//       u_max,  // aMaxInput
      vect<double,4>(35.0, 5.0,  5.0,  5.0),
      vect<double,4>(100.0, 25.0, 25.0, 25.0),  // aInputBandwidth
      mat<double,mat_structure::diagonal>(weightR_mat * Rscale),
      mat<double,mat_structure::diagonal>(weightQ_mat * Qscale),
      0.01, // aTimeStep = 0.1,
      20.0, // aMaxTimeHorizon = 10.0,
      0.1)); //aGoalProximityThreshold = 1.0)
  
  typedef MEAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > MEAQR_space_type;
  
  shared_ptr< MEAQR_space_type > quad_MEAQR_space(new MEAQR_space_type(
    "QuadRotor_MEAQR_topology",
    quad_space,
    0.01));  // IdleToCostRatio 
  
  xml_oarchive file_out("models/quadrotor_spaces.xml");
  
  file_out << quad_sys << quad_space << quad_MEAQR_space;
  
  {
    serialization::xml_oarchive out("models/X8_planner_params.xml");
    out & RK_SERIAL_SAVE_WITH_ALIAS("planner_kind", (std::string("rrt_star")))
        & RK_SERIAL_SAVE_WITH_ALIAS("max_vertices", (std::size_t(1000))) 
        & RK_SERIAL_SAVE_WITH_ALIAS("prog_interval", (std::size_t(100))) 
        & RK_SERIAL_SAVE_WITH_ALIAS("max_results", (std::size_t(50))) 
        & RK_SERIAL_SAVE_WITH_ALIAS("knn_method", (std::size_t(pp::DVP_BF2_TREE_KNN)));
  };
  
  return 0;
};







