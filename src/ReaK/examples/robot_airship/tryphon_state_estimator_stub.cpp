// Standard C/C++ libraries
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <lin_alg/vect_alg.hpp>
#include <ss_systems/airship_assembled_models.hpp>
#include <ctrl_sys/tsos_aug_inv_kalman_filter.hpp>
#include <serialization/archiver_factory.hpp>

ReaK::vect_n<double> sonar_meas = ReaK::vect_n<double>(4, 0.0);

ReaK::vect<double,3> acc_meas;
ReaK::vect<double,3> mag_meas;
ReaK::vect<double,3> gyro_meas;

ReaK::vect<double,3> last_force;
ReaK::vect<double,3> last_torque;

double compass_angle = 0.0;

void subSonar()
{
  for(std::size_t i = 0; i < sonar_meas.size(); ++i)
    sonar_meas[i] = 0.0;
  sonar_meas[0] += 1.0;
}

void subImu()
{
  acc_meas[0] = 0.0;
  acc_meas[1] = 0.0;
  acc_meas[2] = 0.0;
  gyro_meas[0] = 0.0;
  gyro_meas[1] = 0.0;
  gyro_meas[2] = 0.0;
  mag_meas[0] = 0.0;
  mag_meas[1] = 0.0;
  mag_meas[2] = 0.0;
}

void subForces()
{
  last_force[0] = 0.0;
  last_force[1] = 0.0;
  last_force[2] = 0.0;  // TODO verify that this is correct (I assume body-fixed force / torque).
  last_torque[0] = 0.0;
  last_torque[1] = 0.0;
  last_torque[2] = 0.0;
}


int main(int argc, char **argv)
{
  
  // ----------- General airship configs -----------
  
  ReaK::ctrl::airship_parameter_pack a_params;
  a_params.mass = 1.0;
  a_params.J = ReaK::mat<double,ReaK::mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0);
  a_params.added_mass_factor = 0.5;
  
  a_params.use_hot_del_q_terms = false;
  a_params.use_momentum_transfer_terms = true;
  
  a_params.gravity_acc_vect    = ReaK::vect<double,3>(0.0,0.0,-9.81);
  a_params.magnetic_field_vect = ReaK::vect<double,3>(0.0,0.0,0.0);
  
  a_params.IMU_position = ReaK::vect<double,3>(0.0,0.0,0.0);            // relative to body-fixed frame.
  a_params.IMU_orientation = ReaK::quaternion<double>(ReaK::vect<double,4>(1.0,0.0,0.0,0.0)); // relative to body-fixed frame.
  
  // Load from a file:
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_sys_params.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(a_params);
  } catch(std::exception& e) {
    std::cerr << "Could not load the Tryphon's system parameters, with exception: " << e.what() << std::endl;
    return 2;
  };
  
  
  // ----------- Configs of sonars -----------
  
  ReaK::ctrl::sonars_in_room_output_model sonars(6);
  sonars.lower_corner = ReaK::vect<double,3>(0.0,0.0,0.0); // lower-corner of room, in meters.
  sonars.upper_corner = ReaK::vect<double,3>(0.0,0.0,0.0); // upper-corner of room, in meters.
  
  sonars.sonar_pos[0] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[0] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[1] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[1] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[2] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[2] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[3] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[3] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[4] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[4] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[5] = ReaK::vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[5] = ReaK::vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  // Load from a file:
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_sonar_parameters.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(sonars);
    // TODO Add the loading of sonar IDs.
  } catch(std::exception& e) {
    std::cerr << "Could not load the Tryphon's sonar parameters, with exception: " << e.what() << std::endl;
    return 5;
  };
  
  typedef ReaK::ctrl::tryphon_sgam_me SysType;
  
  SysType sys(
    "tryphon_sgam_me", 
    a_params,
    ReaK::arithmetic_tuple<
      ReaK::ctrl::satellite_state_model, 
      ReaK::ctrl::near_buoyancy_state_model, 
      ReaK::ctrl::eccentricity_state_model>(),
    ReaK::arithmetic_tuple<ReaK::ctrl::airship3D_6dof_thrusters>(),
    ReaK::arithmetic_tuple< 
      ReaK::ctrl::sonars_in_room_output_model,
      ReaK::ctrl::sat_gyros_output_model,
      ReaK::ctrl::sat_accelerometer_output_model,
      ReaK::ctrl::sat_magnetometer_output_model
    >(sonars, ReaK::ctrl::sat_gyros_output_model(), ReaK::ctrl::sat_accelerometer_output_model(), ReaK::ctrl::sat_magnetometer_output_model()),
    0.01);
  
#if 0
  typedef ReaK::ctrl::tryphon_sgam_megam SysType;
  
  SysType sys(
    "tryphon_sgam_megam", 
    a_params,
    ReaK::arithmetic_tuple<
      ReaK::ctrl::satellite_state_model, 
      ReaK::ctrl::near_buoyancy_state_model, 
      ReaK::ctrl::eccentricity_state_model, 
      ReaK::ctrl::gyros_bias_state_model,
      ReaK::ctrl::accelerometer_bias_state_model,
      ReaK::ctrl::magnetometer_bias_state_model>(),
    ReaK::arithmetic_tuple<ReaK::ctrl::airship3D_6dof_thrusters>(),
    ReaK::arithmetic_tuple< 
      ReaK::ctrl::sonars_in_room_output_model,
      ReaK::ctrl::sat_gyros_output_model,
      ReaK::ctrl::sat_accelerometer_output_model,
      ReaK::ctrl::sat_magnetometer_output_model
    >(sonars, ReaK::ctrl::sat_gyros_output_model(), ReaK::ctrl::sat_accelerometer_output_model(), ReaK::ctrl::sat_magnetometer_output_model()),
    0.01);
#endif
  
  typedef typename SysType::state_space_type SpaceType;
  typedef typename SysType::covar_type CovarType;
  typedef typename CovarType::matrix_type CovarMatType;
  typedef typename SysType::state_belief_type StateBeliefType;
  typedef typename SysType::input_belief_type InputBeliefType;
  typedef typename SysType::output_belief_type OutputBeliefType;
  
  ReaK::shared_ptr< SpaceType > sat_space = sys.get_state_space();
  
  StateBeliefType b = sys.get_zero_state_belief(10.0);
  
  ReaK::mat<double,ReaK::mat_structure::diagonal> input_disturbance;
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_disturbances.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(input_disturbance);
  } catch(std::exception& e) {
    std::cerr << "Could not load the input-disturbance matrix, with exception: " << e.what() << std::endl;
  };
  
  InputBeliefType b_u = sys.get_zero_input_belief();
  b_u.set_covariance(CovarType(CovarMatType(input_disturbance)));
  
  ReaK::mat<double,ReaK::mat_structure::diagonal> measurement_noise;
  try {
    // TODO change the file name.
    *(ReaK::serialization::open_iarchive("tryphon_noise_4s.rkx"))
      & RK_SERIAL_LOAD_WITH_NAME(measurement_noise);
  } catch(std::exception& e) {
    std::cerr << "Could not load the measurement-noise matrix, with exception: " << e.what() << std::endl;
  };
  
  OutputBeliefType b_z = sys.get_zero_output_belief();
  b_z.set_covariance(CovarType(CovarMatType(measurement_noise)));
  
  ReaK::vect_n<double> z(sys.get_output_dimensions(), 0.0);
  ReaK::vect_n<double> u(sys.get_input_dimensions(), 0.0);
  
  using ReaK::get;
  
  while (true)
  {
    ////////////////////////////////////
    ////       State estimator      ////
    ////////////////////////////////////
    
    std::size_t i = 0;
    for(; i < sonar_meas.size(); ++i)
      z[i] = sonar_meas[i];
    z[i++] = gyro_meas[0];
    z[i++] = gyro_meas[1];
    z[i++] = gyro_meas[2];
    z[i++] = acc_meas[0];
    z[i++] = acc_meas[1];
    z[i++] = acc_meas[2];
    z[i++] = mag_meas[0];
    z[i++] = mag_meas[1];
    z[i++] = mag_meas[2];
    
    b_z.set_mean_state(z);
    
    
    u[0] = last_force[0];
    u[1] = last_force[1];
    u[2] = last_force[2];
    u[3] = last_torque[0];
    u[4] = last_torque[1];
    u[5] = last_torque[2];
    
    b_u.set_mean_state(u);
    
    
    ReaK::ctrl::tsos_aug_inv_kalman_filter_step(sys, *sat_space, b, b_u, b_z, 0.0);
    
    
    const ReaK::vect<double,3>& pos = get<0>(get<0>(get<0>(b.get_mean_state())));
    const ReaK::vect<double,3>& vel = get<1>(get<0>(get<0>(b.get_mean_state())));
    const ReaK::unit_quat<double>& quat = get<0>(get<1>(get<0>(b.get_mean_state())));
    const ReaK::vect<double,3>& omega = get<1>(get<1>(get<0>(b.get_mean_state())));
    
    double dummy;
    dummy = pos[0];
    dummy = pos[1];
    dummy = pos[2];
    dummy = quat[0];
    dummy = quat[1];
    dummy = quat[2];
    dummy = quat[3];
    dummy = vel[0];
    dummy = vel[1];
    dummy = vel[2];
    dummy = omega[0];
    dummy = omega[1];
    dummy = omega[2];
    
    dummy = get<1>(b.get_mean_state());
    dummy = get<2>(b.get_mean_state())[0];
    dummy = get<2>(b.get_mean_state())[1];
    dummy = get<2>(b.get_mean_state())[2];
    dummy = 0.0;
    dummy = 0.0;
    
    RK_UNUSED(dummy);
    
  }
  return 0;
}
