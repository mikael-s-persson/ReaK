
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


#include "airship3D_lin_model.hpp"

#include "serialization/xml_archiver.hpp"
#include "recorders/ssv_recorder.hpp"

#include "ctrl_sys/gaussian_belief_state.hpp"
#include "ctrl_sys/covariance_matrix.hpp"

#include "ctrl_sys/kalman_filter.hpp"
#include "ctrl_sys/kalman_bucy_filter.hpp"
#include "ctrl_sys/invariant_kalman_filter.hpp"
#include "ctrl_sys/invariant_kalman_bucy_filter.hpp"
#include "ctrl_sys/unscented_kalman_filter.hpp"


#include "integrators/fixed_step_integrators.hpp"


#include "mbd_kte/inertia.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/state_measures.hpp"
#include "mbd_kte/free_joints.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "ctrl_sys/kte_nl_system.hpp"
#include "ctrl_sys/num_int_dtnl_system.hpp"
#include "integrators/variable_step_integrators.hpp"


#include "boost/date_time/posix_time/posix_time.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

void print_usage() {
  std::cout << "Usage:\n"
            << "\t./estimate_airship3D [model_filename.xml] [init_conditions.xml] [inertia_data.xml] [result_filename] [time_step] [end_time] [Qu.xml] [R.xml] [skips_min] [skips_max]\n"
            << "\t\t model_filename.xml:\t The filename for the airship model to use.\n"
            << "\t\t init_conditions.xml:\t The filename for the airship's initial conditions.\n"
            << "\t\t inertia_data.xml:\t The filename for the airship's inertial data.\n"
            << "\t\t result_filename:\t The filename where to record the results as a space-separated values file.\n"
            << "\t\t time_step:\t\t The time-step of the data points in the model integration.\n"
            << "\t\t end_time:\t\t The end-time of the data points in the model integration.\n"
            << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
            << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix.\n"
            << "\t\t skips_min:\t\t Minimum number of data rows to skip from measurements.\n"
            << "\t\t skips_max:\t\t Maximum number of data rows to skip from measurements.\n"
            << "\t\t mc_count:\t\t Number of Monte-Carlo runs to perform." << std::endl;
};


bool load_parameters(int argc, char** argv, std::string& results_filename, 
                     double& time_step, double& end_time,
                     ReaK::ctrl::kte_nl_system& airship3D_system, 
                     ReaK::vect_n<double>& x_0, double& mass, 
                     ReaK::mat<double,ReaK::mat_structure::symmetric>& inertia_tensor,
                     ReaK::mat<double,ReaK::mat_structure::diagonal>& Qu, 
                     ReaK::mat<double,ReaK::mat_structure::diagonal>& R,
                     unsigned int& skips_min, unsigned int& skips_max, unsigned int& mc_count) {
  
  if(argc < 12) {
    print_usage();
    return false;
  };
  
  
  ReaK::shared_ptr< ReaK::frame_3D<double> >        airship3D_frame;
  ReaK::shared_ptr< ReaK::kte::position_measure_3D > airship3D_position;
  ReaK::shared_ptr< ReaK::kte::rotation_measure_3D > airship3D_rotation;
  
  ReaK::shared_ptr< ReaK::kte::driving_actuator_3D > airship3D_actuator;
  ReaK::shared_ptr< ReaK::kte::inertia_3D >          airship3D_inertia;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain >       airship3D_model;
  ReaK::shared_ptr< ReaK::kte::mass_matrix_calc >    airship3D_mass_calc;
  
  std::string model_filename(argv[1]);
  try {
    ReaK::serialization::xml_iarchive in(model_filename);
    in >> airship3D_frame
       >> airship3D_position
       >> airship3D_rotation
       >> airship3D_actuator
       >> airship3D_inertia
       >> airship3D_model
       >> airship3D_mass_calc;
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the airship's model!");
    return false;
  };
  
  /* initial states */
  std::string init_filename(argv[2]);
  try {
    ReaK::serialization::xml_iarchive in(init_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("initial_motion",(*airship3D_frame));
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return false;
  };
  x_0.resize(13);
  x_0[0] = airship3D_frame->Position[0];
  x_0[1] = airship3D_frame->Position[1];
  x_0[2] = airship3D_frame->Position[2];
  x_0[3] = airship3D_frame->Quat[0];
  x_0[4] = airship3D_frame->Quat[1];
  x_0[5] = airship3D_frame->Quat[2];
  x_0[6] = airship3D_frame->Quat[3];
  x_0[7] = airship3D_frame->Velocity[0];
  x_0[8] = airship3D_frame->Velocity[1];
  x_0[9] = airship3D_frame->Velocity[2];
  x_0[10] = airship3D_frame->AngVelocity[0];
  x_0[11] = airship3D_frame->AngVelocity[1];
  x_0[12] = airship3D_frame->AngVelocity[2];
  
  /* inertial data */
  std::string inertia_filename(argv[3]);
  try {
    ReaK::serialization::xml_iarchive in(inertia_filename);
    double tmp_mass;
    ReaK::mat<double,ReaK::mat_structure::symmetric> tmp_tensor;
    in & RK_SERIAL_LOAD_WITH_ALIAS("mass",tmp_mass)
       & RK_SERIAL_LOAD_WITH_ALIAS("inertia_tensor",tmp_tensor);
    airship3D_inertia->setMass(tmp_mass);
    airship3D_inertia->setInertiaTensor(tmp_tensor);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return false;
  };
  mass = airship3D_inertia->Mass();
  inertia_tensor = airship3D_inertia->InertiaTensor();
  
  airship3D_system.dofs_3D.push_back(airship3D_frame);
  airship3D_system.inputs.push_back(airship3D_actuator);
  airship3D_system.outputs.push_back(airship3D_position);
  airship3D_system.outputs.push_back(airship3D_rotation);
  airship3D_system.chain = airship3D_model;
  airship3D_system.mass_calc = airship3D_mass_calc;
  
  
  results_filename = std::string(argv[4]);
  
  time_step = 0.001;
  std::stringstream(argv[5]) >> time_step;
  
  end_time = 0.0;
  std::stringstream(argv[6]) >> end_time;
  
  
  /* input disturbance */
  std::string Qu_filename(argv[7]);
  try {
    ReaK::serialization::xml_iarchive in(Qu_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("input_disturbance",Qu);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return false;
  };
  
  /* measurement noise */
  std::string R_filename(argv[8]);
  try {
    ReaK::serialization::xml_iarchive in(R_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("measurement_noise",R);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return false;
  };
  
  skips_min = 0;
  std::stringstream(argv[9]) >> skips_min;
  ++skips_min;
  
  skips_max = 0;
  std::stringstream(argv[10]) >> skips_max;
  ++skips_max;
  
  mc_count = 100;
  std::stringstream(argv[11]) >> mc_count;
  
  return true;
};




int main(int argc, char** argv) {
  using namespace ReaK;
  
  std::string result_filename;
  double time_step, end_time;
  shared_ptr<ctrl::kte_nl_system> airship3D_system(new ctrl::kte_nl_system("airship3D_system"));
  vect_n<double> x_0;
  double mass;
  mat<double,mat_structure::symmetric> inertia_tensor;
  mat<double,mat_structure::diagonal> Qu;
  mat<double,mat_structure::diagonal> R;
  unsigned int skips_min, skips_max, mc_count;
  
  if( !load_parameters(argc, argv, result_filename, time_step, end_time, *airship3D_system, x_0, mass, inertia_tensor, Qu, R, skips_min, skips_max, mc_count) )
    return 1;
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  typedef ctrl::num_int_dtnl_sys< ctrl::kte_nl_system, dormand_prince45_integrator<double> > sys_type;
  sys_type 
    airship3D_dt_sys( airship3D_system, 
                      dormand_prince45_integrator<double>("ode45_integrator",
                                                          ReaK::vect_n<double>(),
                                                          0.0,
                                                          time_step * 0.01,
                                                          weak_ptr< state_rate_function<double> >(),
                                                          time_step,
                                                          time_step * 0.00001,
                                                          1e-2),
                      time_step,
                      "airship3D_dt_sys");
  
  ctrl::airship3D_lin_dt_system airship3D_mdl_dt("airship3D_linear_discrete",mass,inertia_tensor,time_step);
  
  ctrl::airship3D_lin_system airship3D_mdl("airship3D_linear",mass,inertia_tensor);
  
  typedef ctrl::num_int_dtnl_sys< ctrl::airship3D_lin_system, dormand_prince45_integrator<double> > sys2_type;
  sys2_type 
    airship3D_dt_sys2(shared_ptr<ctrl::airship3D_lin_system>(&airship3D_mdl, null_deleter()), 
                      dormand_prince45_integrator<double>("ode45_integrator",
                                                          ReaK::vect_n<double>(),
                                                          0.0,
                                                          time_step * 0.01,
                                                          weak_ptr< state_rate_function<double> >(),
                                                          time_step,
                                                          time_step * 0.00001,
                                                          1e-2),
                      time_step,
                      "airship3D_dt_sys2");
  
  std::vector<double> std_devs(12 * (1 + skips_max - skips_min));
  
  pp::vector_topology< vect_n<double> > mdl_state_space;
  
  for(unsigned int i = 0; i < mc_count; ++i) {
    
    
    std::cout << std::endl << "Starting Kalman filtering..." << std::endl;
    
    vect_n<double> x_init(13);
    x_init[0] = 0.0; x_init[1] = 0.0; x_init[2] = 0.0; 
    x_init[3] = 1.0; x_init[4] = 0.0; x_init[5] = 0.0; x_init[6] = 0.0;
    x_init[7] = 0.0; x_init[8] = 0.0; x_init[9] = 0.0; 
    x_init[10] = 0.0; x_init[11] = 0.0; x_init[12] = 0.0;
    ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
      b_init(x_init,
             ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(13,1000.0))));
  
    ctrl::covariance_matrix< vect_n<double> > Rcov = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R));
    
    mat<double,mat_structure::diagonal> Qu_avg = Qu;
  
    for(unsigned int j = skips_min; j <= skips_max; ++j) {
      
      // Qu_avg = Qu; 
      Qu_avg = (1.0 / double(j)) * Qu;
      // Qu_avg = (double(j)) * Qu;
      
      
      std::list< std::pair< double, vect_n<double> > > measurements;
      std::list< std::pair< double, vect_n<double> > > measurements_noisy;
      typedef std::list< std::pair< double, vect_n<double> > >::const_iterator MeasIter;
      
      std::cout << "Starting simulation run... " << i << std::endl;
      //airship3D_dt_sys.set_time_step(j * time_step);
      airship3D_dt_sys2.set_time_step(j * time_step);
      //airship3D_mdl_dt.set_time_step(j * time_step);
      try {
        sys_type::point_type x = x_0;
        for(double t = 0.0; t < end_time; t += j * time_step) {
          
          sys_type::input_type u = sys_type::input_type(6);
          u[0] = var_rnd() * sqrt(Qu_avg(0,0));
          u[1] = var_rnd() * sqrt(Qu_avg(1,1));
          u[2] = var_rnd() * sqrt(Qu_avg(2,2));
          u[3] = var_rnd() * sqrt(Qu_avg(3,3));
          u[4] = var_rnd() * sqrt(Qu_avg(4,4));
          u[5] = var_rnd() * sqrt(Qu_avg(5,5));
          
          //x = airship3D_dt_sys.get_next_state(mdl_state_space,x,u,t);
          //sys_type::output_type y = airship3D_dt_sys.get_output(mdl_state_space,x,u,t);
          
          x = airship3D_dt_sys2.get_next_state(mdl_state_space,x,u,t);
          sys_type::output_type y = airship3D_dt_sys2.get_output(mdl_state_space,x,u,t);
          
          //x = airship3D_mdl_dt.get_next_state(mdl_state_space, x, u, t);
          //sys_type::output_type y = airship3D_mdl_dt.get_output(mdl_state_space,x,u,t);
          
          std::cout << "\r" << std::setw(20) << t; std::cout.flush();
          
          measurements.push_back( std::make_pair(t, y) );
          sys_type::output_type y_noisy = y;
          y_noisy[0] += var_rnd() * sqrt(R(0,0));
          y_noisy[1] += var_rnd() * sqrt(R(1,1));
          y_noisy[2] += var_rnd() * sqrt(R(2,2));
          y_noisy[3] += var_rnd() * sqrt(R(3,3));
          y_noisy[4] += var_rnd() * sqrt(R(4,4));
          y_noisy[5] += var_rnd() * sqrt(R(5,5));
          y_noisy[6] += var_rnd() * sqrt(R(6,6));
          measurements_noisy.push_back( std::make_pair(t, y_noisy) );
        };
         
      } catch(impossible_integration& e) {
        std::cout << "Integration was deemed impossible, with message: '" << e.what() << "'" << std::endl;
        return 3;
      } catch(untolerable_integration& e) {
        std::cout << "Integration was deemed untolerable with message: '" << e.what() << "'" << std::endl;
        return 4;
      };
      std::cout << std::endl << "Done." << std::endl;
      
      //compute std-dev of the measurements:
      { 
        MeasIter it = measurements.begin();
        MeasIter it_noisy = measurements_noisy.begin();
        double std_dev_tmp[2] = {0.0,0.0}; unsigned int k = 0;
        while(it != measurements.end()) {
          std_dev_tmp[0] = ( k * std_dev_tmp[0] + ( (it_noisy->second[0] - it->second[0]) * (it_noisy->second[0] - it->second[0])
                                                  + (it_noisy->second[1] - it->second[1]) * (it_noisy->second[1] - it->second[1])
                                                  + (it_noisy->second[2] - it->second[2]) * (it_noisy->second[2] - it->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(it_noisy->second[3],it_noisy->second[4],it_noisy->second[5],it_noisy->second[6]));
          quaternion<double> q_true(vect<double,4>(it->second[3],it->second[4],it->second[5],it->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev_tmp[1] = ( k * std_dev_tmp[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++it; ++it_noisy; ++k;
        };
        std_devs[12 * (j - skips_min)]     = (i * std_devs[12 * (j - skips_min)]     + std::sqrt(std_dev_tmp[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 1] = (i * std_devs[12 * (j - skips_min) + 1] + std::sqrt(std_dev_tmp[1])) / (double(i+1));
      };
      
      
      
      
      std::cout << "\r" << std::setw(20) << j * time_step; std::cout.flush();

      ctrl::airship3D_lin2_dt_system mdl_lin2_dt("airship3D_linear2_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_dt_system mdl_inv_dt("airship3D_invariant_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv2_dt_system mdl_inv2_dt("airship3D_invariant2_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_mom_dt_system mdl_inv_mom_dt("airship3D_invariant_momentum_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_mid_dt_system mdl_inv_mid_dt("airship3D_invariant_midpoint_discrete",mass,inertia_tensor,time_step * j);
      
      
      std::cout << "Running Extended Kalman Filter version 2..." << std::endl;
      try {
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b = b_init;
        ctrl::covariance_matrix< vect_n<double> > Qcov;
        Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Rcov);
        
        double std_dev[2] = {0.0,0.0}; int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
          b_z.set_mean_state(it->second);
          ctrl::kalman_filter_step(mdl_lin2_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
          
          const vect_n<double>& x_mean = b.get_mean_state();
          
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[0]) * (x_mean[0] - it_orig->second[0])
                                          + (x_mean[1] - it_orig->second[1]) * (x_mean[1] - it_orig->second[1])
                                          + (x_mean[2] - it_orig->second[2]) * (x_mean[2] - it_orig->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(x_mean[3],x_mean[4],x_mean[5],x_mean[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[3],it_orig->second[4],it_orig->second[5],it_orig->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++k;
          
          //for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          ++it; ++it_orig;
        };
        std_devs[12 * (j - skips_min) + 2] = (i * std_devs[12 * (j - skips_min) + 2] + std::sqrt(std_dev[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 3] = (i * std_devs[12 * (j - skips_min) + 3] + std::sqrt(std_dev[1])) / (double(i+1));
      } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
        std::cout << "The Extended Kalman Filter (v2) has thrown a singularity error!" << std::endl;
        std_devs[12 * (j - skips_min) + 2] = std::numeric_limits<double>::infinity();
        std_devs[12 * (j - skips_min) + 3] = std::numeric_limits<double>::infinity();
      };
      
      std::cout << "Running Invariant Extended Kalman Filter..." << std::endl;
      try {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,1000.0))));
  
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
        ctrl::covariance_matrix< vect_n<double> > Qcov;
        Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Rcovinv);
        
        double std_dev[2] = {0.0,0.0}; int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
    
          const vect_n<double>& x_mean = b.get_mean_state();
          
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[0]) * (x_mean[0] - it_orig->second[0])
                                          + (x_mean[1] - it_orig->second[1]) * (x_mean[1] - it_orig->second[1])
                                          + (x_mean[2] - it_orig->second[2]) * (x_mean[2] - it_orig->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(x_mean[3],x_mean[4],x_mean[5],x_mean[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[3],it_orig->second[4],it_orig->second[5],it_orig->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++k;
          
          //for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          ++it; ++it_orig;
        };
        std_devs[12 * (j - skips_min) + 4] = (i * std_devs[12 * (j - skips_min) + 4] + std::sqrt(std_dev[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 5] = (i * std_devs[12 * (j - skips_min) + 5] + std::sqrt(std_dev[1])) / (double(i+1));
      } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
        std::cout << "The Invariant Extended Kalman Filter has thrown a singularity error!" << std::endl;
        std_devs[12 * (j - skips_min) + 4] = std::numeric_limits<double>::infinity();
        std_devs[12 * (j - skips_min) + 5] = std::numeric_limits<double>::infinity();
      };
      
      
      std::cout << "Running Invariant Extended Kalman Filter (version 2)..." << std::endl;
      try {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,1000.0))));
  
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
        ctrl::covariance_matrix< vect_n<double> > Qcov;
        Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Rcovinv);
        
        double std_dev[2] = {0.0,0.0}; int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv2_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
    
          const vect_n<double>& x_mean = b.get_mean_state();
          
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[0]) * (x_mean[0] - it_orig->second[0])
                                          + (x_mean[1] - it_orig->second[1]) * (x_mean[1] - it_orig->second[1])
                                          + (x_mean[2] - it_orig->second[2]) * (x_mean[2] - it_orig->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(x_mean[3],x_mean[4],x_mean[5],x_mean[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[3],it_orig->second[4],it_orig->second[5],it_orig->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++k;
          
          //for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          ++it; ++it_orig;
        };
        std_devs[12 * (j - skips_min) + 6] = (i * std_devs[12 * (j - skips_min) + 6] + std::sqrt(std_dev[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 7] = (i * std_devs[12 * (j - skips_min) + 7] + std::sqrt(std_dev[1])) / (double(i+1));
      } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
        std::cout << "The Invariant Extended Kalman Filter (version 2) has thrown a singularity error!" << std::endl;
        std_devs[12 * (j - skips_min) + 6] = std::numeric_limits<double>::infinity();
        std_devs[12 * (j - skips_min) + 7] = std::numeric_limits<double>::infinity();
      };
      
      
      
      std::cout << "Running Invariant-Momentum Kalman Filter..." << std::endl;
      try {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,1000.0))));
  
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
        ctrl::covariance_matrix< vect_n<double> > Qcov;
        Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Rcovinv);
        
        double std_dev[2] = {0.0,0.0}; int k = 0;
        MeasIter it_orig = measurements.begin();
        for(MeasIter it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_mom_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
          
          const vect_n<double>& x_mean = b.get_mean_state();
          
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[0]) * (x_mean[0] - it_orig->second[0])
                                          + (x_mean[1] - it_orig->second[1]) * (x_mean[1] - it_orig->second[1])
                                          + (x_mean[2] - it_orig->second[2]) * (x_mean[2] - it_orig->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(x_mean[3],x_mean[4],x_mean[5],x_mean[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[3],it_orig->second[4],it_orig->second[5],it_orig->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++k;
          
          //for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          ++it; ++it_orig;
        };
        std_devs[12 * (j - skips_min) + 8] = (i * std_devs[12 * (j - skips_min) + 8] + std::sqrt(std_dev[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 9] = (i * std_devs[12 * (j - skips_min) + 9] + std::sqrt(std_dev[1])) / (double(i+1));
      } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
        std::cout << "The Invariant-Momentum Kalman Filter (v1) has thrown a singularity error!" << std::endl;
        std_devs[12 * (j - skips_min) + 8] = std::numeric_limits<double>::infinity();
        std_devs[12 * (j - skips_min) + 9] = std::numeric_limits<double>::infinity();
      };
  
      std::cout << "Running Invariant-Midpoint Kalman Filter..." << std::endl;
      try {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,1000.0))));
  
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
        ctrl::covariance_matrix< vect_n<double> > Qcov;
        Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                     Rcovinv);
        
        double std_dev[2] = {0.0,0.0}; int k = 0;
        MeasIter it_orig = measurements.begin();
        for(MeasIter it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_mid_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
          
          const vect_n<double>& x_mean = b.get_mean_state();
          
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[0]) * (x_mean[0] - it_orig->second[0])
                                          + (x_mean[1] - it_orig->second[1]) * (x_mean[1] - it_orig->second[1])
                                          + (x_mean[2] - it_orig->second[2]) * (x_mean[2] - it_orig->second[2]) ) ) / (k + 1);
          
          quaternion<double> q_noisy(vect<double,4>(x_mean[3],x_mean[4],x_mean[5],x_mean[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[3],it_orig->second[4],it_orig->second[5],it_orig->second[6]));
          axis_angle<double> aa_diff( invert(q_true) * q_noisy );
          
          std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
          
          ++k;
          
          //for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          ++it; ++it_orig;
        };
        std_devs[12 * (j - skips_min) + 10] = (i * std_devs[12 * (j - skips_min) + 10] + std::sqrt(std_dev[0])) / (double(i+1));
        std_devs[12 * (j - skips_min) + 11] = (i * std_devs[12 * (j - skips_min) + 11] + std::sqrt(std_dev[1])) / (double(i+1));
      } catch(ReaK::singularity_error& e) { RK_UNUSED(e);
        std::cout << "The Invariant-Momentum Kalman Filter (v2) has thrown a singularity error!" << std::endl;
        std_devs[12 * (j - skips_min) + 10] = std::numeric_limits<double>::infinity();
        std_devs[12 * (j - skips_min) + 11] = std::numeric_limits<double>::infinity();
      };
    };
    std::cout << std::endl << "Done." << std::endl;
    
    
    
    
    
  };
  
  
  recorder::ssv_recorder results(result_filename + "_stddevs.ssv");
  results << "time_step" << "meas_p" << "meas_a" 
                         << "ekf2_p" << "ekf2_a" 
                         << "iekf_p" << "iekf_a" 
                         << "iekf2_p" << "iekf2_a" 
                         << "imkf_q" << "imkf_a" 
                         << "imkf2_p" << "imkf2_a" 
                         << recorder::data_recorder::end_name_row;
  
  for(unsigned int j = skips_min; j <= skips_max; ++j) {
    
    results << (time_step * j);
    
    results << std_devs[12 * (j - skips_min)]
            << std_devs[12 * (j - skips_min) + 1]
            << std_devs[12 * (j - skips_min) + 2]
            << std_devs[12 * (j - skips_min) + 3]
            << std_devs[12 * (j - skips_min) + 4]
            << std_devs[12 * (j - skips_min) + 5]
            << std_devs[12 * (j - skips_min) + 6]
            << std_devs[12 * (j - skips_min) + 7]
            << std_devs[12 * (j - skips_min) + 8]
            << std_devs[12 * (j - skips_min) + 9]
            << std_devs[12 * (j - skips_min) + 10]
            << std_devs[12 * (j - skips_min) + 11];
    
    results << recorder::data_recorder::end_value_row;
  };

  results << recorder::data_recorder::flush;
  
  
};









