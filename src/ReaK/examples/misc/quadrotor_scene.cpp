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
#include "quadrotor_system.hpp"

#include "MEAQR_rrtstar_planner.hpp"
#include "MEAQR_sbastar_planner.hpp"

int main(int argc, char ** argv) {
  using namespace ReaK;
  using namespace pp;
  using namespace ctrl;
  
  shared_ptr< quadrotor_system > quad_sys( new quadrotor_system(
    "Quadrotor_system", 
    2.025, // aMass, 
    mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(vect<double,3>(0.0613, 0.0612, 0.1115))), // aInertiaMoment
    mat<double,mat_structure::diagonal>(vect<double,3>(0.1, 0.1, 0.1)),
    mat<double,mat_structure::diagonal>(vect<double,3>(0.1, 0.1, 0.1))));
  
  
  vect<double,3> min_corner(0.0, 0.0, 0.0);  // min corner
  vect<double,3> max_corner(5.0, 5.0, 5.0);  // mix corner
  double v_max = 6.0;
  double w_max = M_PI;
  vect<double,4> u_max(35.0, 5.0,  5.0,  3.0);
  mat<double,mat_structure::diagonal> weightR_mat(vect<double,4>(0.25, 0.5, 0.5, 0.5));
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
  
  typedef IHAQR_topology< quadrotor_system::state_space_type, quadrotor_system > IHAQR_space_type;
  
  IHAQR_space_type quad_space(
      "Quadrotor_IHAQR_topology",
      quad_sys,
      make_se3_space(
        "Quadrotor_state_space",
        min_corner,  // min corner
        max_corner,  // mix corner
        v_max,    // aMaxSpeed
        w_max),  // aMaxAngularSpeed
      vect<double,4>(0.0, -5.0, -5.0, -3.0),  // aMinInput
      u_max,  // aMaxInput
      vect<double,4>(10.0, 10.0, 10.0, 10.0),  // aInputBandwidth
      mat<double,mat_structure::diagonal>(weightR_mat * Rscale),
      mat<double,mat_structure::diagonal>(weightQ_mat * Qscale),
      0.01, // aTimeStep = 0.1,
      10.0, // aMaxTimeHorizon = 10.0,
      0.1); //aGoalProximityThreshold = 1.0)
  {
    IHAQR_space_type::point_type p1 = quad_space.random_point();
    IHAQR_space_type::point_type p2 = quad_space.random_point();
    
    IHAQR_space_type::point_type p_inter = quad_space.move_position_toward(p1, 0.5, p2);
    
    double dist = quad_space.distance(p1, p2); RK_UNUSED(dist);
  }; 
  
  typedef MEAQR_topology< quadrotor_system::state_space_type, quadrotor_system > MEAQR_space_type;
  
  MEAQR_space_type quad_MEAQR_space(
    "QuadRotor_MEAQR_topology",
    shared_ptr< IHAQR_space_type >(&quad_space, null_deleter()),
    0.02,
    1.0);
  
  {
    MEAQR_space_type::point_type p1 = quad_MEAQR_space.random_point();
    MEAQR_space_type::point_type p2 = quad_MEAQR_space.random_point();
    
    std::cout << " p1 = " << p1.x << std::endl;
    std::cout << " p2 = " << p2.x << std::endl << std::endl;
    
    for(double t = 0.0; t < 1.0; t += 0.1) { 
      MEAQR_space_type::point_type p_inter = quad_MEAQR_space.move_position_toward(p1, t, p2);
      std::cout << t << " " << p_inter.x << std::endl;
    };
    
    double dist = quad_MEAQR_space.distance(p1, p2); RK_UNUSED(dist);
  };
  
  MEAQR_rrtstar_planner< quadrotor_system::state_space_type, quadrotor_system > planner1;
  MEAQR_sbastar_planner< quadrotor_system::state_space_type, quadrotor_system > planner2;
  
  return 0;
};







