// Standard C/C++ libraries
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <lin_alg/vect_alg.hpp>
#include <ss_systems/airship_assembled_models.hpp>
#include <serialization/archiver_factory.hpp>

ReaK::vect<double,3> current_pos = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::vect<double,3> current_vel = ReaK::vect<double,3>(0.0,0.0,0.0);
ReaK::unit_quat<double> current_quat;
ReaK::vect<double,3> current_ang_vel = ReaK::vect<double,3>(0.0,0.0,0.0);

double current_dm = 0.0;
ReaK::vect<double,3> current_r = ReaK::vect<double,3>(0.0,0.0,0.0);

double current_df = 0.0;
double current_dt = 0.0;

void subState()
{
  // State estimation //
  current_pos[0] = 0.0;
  current_pos[1] = 0.0;
  current_pos[2] = 0.0;
  current_vel[0] = 0.0;
  current_vel[1] = 0.0;
  current_vel[2] = 0.0;
  
  current_quat = ReaK::unit_quat<double>(ReaK::vect<double,4>(1.0, 0.0, 0.0, 0.0));
  current_ang_vel[0] = 0.0;
  current_ang_vel[1] = 0.0;
  current_ang_vel[2] = 0.0;
}

void subDMEcc()
{
  // dm and ecc
  current_dm = 0.0;
  current_r[0] = 0.0;
  current_r[1] = 0.0;
  current_r[2] = 0.0;
}

void subDragCoefs()
{
  // linear and rotational drag coefs
  current_df = 0.0;
  current_dt = 0.0;
}


int main(int argc, char **argv)
{
  
  
  // ----------- General airship configs -----------
  
  ReaK::ctrl::airship_parameter_pack a_params;
  
  a_params.mass = 1.0;
  a_params.J = ReaK::mat<double,ReaK::mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0);
  a_params.added_mass_factor = 0.5;
  
  a_params.gravity_acc_vect = ReaK::vect<double,3>(0.0,0.0,-9.81);
  
  // Load from a file:
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_sys_params.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(a_params);
  } catch(std::exception& e) {
    std::cerr << "Could not load the Tryphon's system parameters, with exception: " << e.what() << std::endl;
    return 2;
  };
  
  
  ReaK::quaternion<double> cur_inv_quat;
  
  while (true)
  {
    ////////////////////////////////////
    ////       Controller           ////
    ////////////////////////////////////
    
    cur_inv_quat = invert(current_quat).as_rotation();
    
    ReaK::vect<double,3> comp_force;
    ReaK::vect<double,3> comp_torque;
    
    ReaK::vect<double,3> g_local = cur_inv_quat * a_params.gravity_acc_vect;
    comp_force  += current_dm * g_local;
    comp_torque += current_dm * (current_r % g_local);
    
    comp_force  += cur_inv_quat * ((-current_df * ReaK::norm_2(current_vel)) * current_vel);
    comp_torque += (-current_dt * ReaK::norm_2(current_ang_vel)) * current_ang_vel;
    
    double dummy;
    dummy = -comp_force[0];
    dummy = -comp_force[1];
    dummy = -comp_force[2];
    dummy = -comp_torque[0];
    dummy = -comp_torque[1];
    dummy = -comp_torque[2];
    
    RK_UNUSED(dummy);
    
  }
  return 0;
}




