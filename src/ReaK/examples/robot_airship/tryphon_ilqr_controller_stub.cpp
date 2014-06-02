// Standard C/C++ libraries
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <lin_alg/vect_alg.hpp>
#include <lin_alg/mat_alg.hpp>
#include <ss_systems/airship_assembled_models.hpp>
#include <serialization/archiver_factory.hpp>

#include <lin_alg/mat_are_solver.hpp>

ReaK::vect<double,3> current_pos = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::vect<double,3> current_vel = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::unit_quat<double> current_quat;
ReaK::vect<double,3> current_ang_vel = ReaK::vect<double,3>(0.0,0.0,0.0);

ReaK::vect<double,3> desired_pos = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::vect<double,3> desired_vel = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::unit_quat<double> desired_quat;
ReaK::vect<double,3> desired_ang_vel = ReaK::vect<double,3>(0.0,0.0,0.0);

void subState()
{
  // State estimation //
  current_pos[0] = 0.0;
  current_pos[1] = 0.0;
  current_pos[2] = 0.0;
  current_vel[0] = 0.0;
  current_vel[1] = 0.0;
  current_vel[2] = 0.0;
  
  current_quat = ReaK::unit_quat<double>(ReaK::vect<double,4>(0.0, 0.0, 0.0, 0.0));
  current_ang_vel[0] = 0.0;
  current_ang_vel[1] = 0.0;
  current_ang_vel[2] = 0.0;
}

void subDesiredState()
{
  // State estimation //
  desired_pos[0] = 0.0;
  desired_pos[1] = 0.0;
  desired_pos[2] = 0.0;
  desired_vel[0] = 0.0;
  desired_vel[1] = 0.0;
  desired_vel[2] = 0.0;
  
  desired_quat = ReaK::unit_quat<double>(ReaK::vect<double,4>(0.0, 0.0, 0.0, 0.0));
  desired_ang_vel[0] = 0.0;
  desired_ang_vel[1] = 0.0;
  desired_ang_vel[2] = 0.0;
}

int main(int argc, char **argv)
{
  
  // ----------- General airship configs -----------
  
  ReaK::ctrl::airship_parameter_pack a_params;
  
  a_params.mass = 1.0;
  a_params.J = ReaK::mat<double,ReaK::mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0);
  a_params.added_mass_factor = 0.5;
  
  // Load from a file:
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_sys_params.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(a_params);
  } catch(std::exception& e) {
    std::cerr << "Could not load the Tryphon's system parameters, with exception: " << e.what() << std::endl;
    return 2;
  };
  
  
  
  // ----------- Offline computation of control gains -----------
  
  ReaK::mat<double, ReaK::mat_structure::diagonal> LQR_Q(12, 1.0);
  LQR_Q(0,0) = 2.0;
  LQR_Q(1,1) = 2.0;
  LQR_Q(2,2) = 2.0;
  ReaK::mat<double, ReaK::mat_structure::diagonal> LQR_R(6, 2.0);
  
  // Load from a file:
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_lqr_costs.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(LQR_Q)
      & RK_SERIAL_LOAD_WITH_NAME(LQR_R);
  } catch(std::exception& e) {
    std::cerr << "Could not load the Tryphon's LQR cost (penalty) matrices, with exception: " << e.what() << std::endl;
    return 3;
  };
  
  
  ReaK::mat<double, ReaK::mat_structure::rectangular> LQR_K(6, 12); // <-- control gain
  
  {
  ReaK::mat<double, ReaK::mat_structure::square> LQR_A_ct(12, 0.0);
  ReaK::set_block(LQR_A_ct, ReaK::mat_ident<double>(3), 0, 3);
  ReaK::set_block(LQR_A_ct, ReaK::mat_ident<double>(3), 6, 9);
  
  ReaK::mat<double, ReaK::mat_structure::rectangular> LQR_B_ct(12, 6, 0.0);
  ReaK::set_block(LQR_B_ct, (1.0 / ((1.0 + a_params.added_mass_factor) * a_params.mass)) * ReaK::mat_ident<double>(3), 3, 0);
  try {
    ReaK::mat<double, ReaK::mat_structure::symmetric> J_inv;
    ReaK::invert_Cholesky(a_params.J, J_inv, 1e-4);
    ReaK::set_block(LQR_B_ct, J_inv, 9, 3);
  } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
    std::cerr << "Warning: Inertia tensor of the system is singular! .. continuing with identity matrix instead" << std::endl;
    ReaK::set_block(LQR_B_ct, ReaK::mat_ident<double>(3), 9, 3);
  };
  
  try {
    ReaK::solve_IHCT_LQR(LQR_A_ct, LQR_B_ct, LQR_Q, LQR_R, LQR_K, 1e-4);
  } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
    std::cerr << "Error: Could not solve the infinite-horizon continuous-time LQR control gain with the given system matrices!" << std::endl;
    return 2;
  };
  };
  
  {
  double time_step = 0.05;
  
  ReaK::mat<double, ReaK::mat_structure::square> LQR_A_dt(12, true);
  ReaK::set_block(LQR_A_dt, time_step * ReaK::mat_ident<double>(3), 0, 3);
  ReaK::set_block(LQR_A_dt, time_step * ReaK::mat_ident<double>(3), 6, 9);
  
  ReaK::mat<double, ReaK::mat_structure::rectangular> LQR_B_dt(12, 6, 0.0);
  ReaK::set_block(LQR_B_dt, (0.5 * time_step * time_step / ((1.0 + a_params.added_mass_factor) * a_params.mass)) * ReaK::mat_ident<double>(3), 0, 0);
  ReaK::set_block(LQR_B_dt, (time_step / ((1.0 + a_params.added_mass_factor) * a_params.mass)) * ReaK::mat_ident<double>(3), 3, 0);
  try {
    ReaK::mat<double, ReaK::mat_structure::symmetric> J_inv;
    ReaK::invert_Cholesky(a_params.J, J_inv, 1e-4);
    ReaK::set_block(LQR_B_dt, (0.5 * time_step * time_step) * J_inv, 6, 3);
    ReaK::set_block(LQR_B_dt, time_step * J_inv, 9, 3);
  } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
    std::cerr << "Warning: Inertia tensor of the system is singular! .. continuing with identity matrix instead" << std::endl;
    ReaK::set_block(LQR_B_dt, (0.5 * time_step * time_step) * ReaK::mat_ident<double>(3), 6, 3);
    ReaK::set_block(LQR_B_dt, time_step * ReaK::mat_ident<double>(3), 9, 3);
  };
  
  try {
    ReaK::solve_IHDT_LQR(LQR_A_dt, LQR_B_dt, LQR_Q, LQR_R, LQR_K, 1e-4, false);
  } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
    std::cerr << "Error: Could not solve the infinite-horizon discrete-time LQR control gain with the given system matrices!" << std::endl;
    return 2;
  };
  };
  
  
  using ReaK::log;
  ReaK::vect_n<double> err(12, 0.0);
  ReaK::quaternion<double> cur_inv_quat;
  ReaK::unit_quat<double> quat_diff;
  
  while (true)
  {
    ////////////////////////////////////
    ////       Controller           ////
    ////////////////////////////////////
    
    // First, compute the invariant error vector:
    cur_inv_quat = invert(current_quat).as_rotation();
    err[ReaK::range(0,3)]  = cur_inv_quat * (desired_pos - current_pos);
    err[ReaK::range(3,6)]  = cur_inv_quat * (desired_vel - current_vel);
    quat_diff = invert(current_quat) * desired_quat;
    err[ReaK::range(6,9)]  = 2.0 * log(quat_diff);
    err[ReaK::range(9,12)] = quat_diff.as_rotation() * desired_ang_vel - current_ang_vel;
    
    ReaK::vect_n<double> ft = LQR_K * err;
    
    double dummy;
    dummy = ft[0];
    dummy = ft[1];
    dummy = ft[2];
    dummy = ft[3];
    dummy = ft[4];
    dummy = ft[5];
    
    RK_UNUSED(dummy);
    
  }
  return 0;
}




