
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


#include "ctrl_sys/kalman_filter.hpp"
#include "ctrl_sys/kalman_bucy_filter.hpp"
#include "ctrl_sys/invariant_kalman_filter.hpp"
#include "ctrl_sys/invariant_kalman_bucy_filter.hpp"
#include "ctrl_sys/unscented_kalman_filter.hpp"

#include "ctrl_sys/gaussian_belief_state.hpp"
#include "ctrl_sys/covariance_matrix.hpp"

#include "integrators/fixed_step_integrators.hpp"

#include "boost/date_time/posix_time/posix_time.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <iomanip>

int main(int argc, char** argv) {
  using namespace ReaK;
  
  if(argc < 8) {
    std::cout << "Usage:\n"
	      << "\t./estimate_airship3D [inertial_data.xml] [meas_filename] [result_filename] [Qu.xml] [R.xml] [added_R.xml] [test_count]\n"
	      << "\t\t inertial_data.xml:\t The filename of the inertial data of the airship3D.\n"
	      << "\t\t meas_filename:\t The filename prefix of space-sep. values files with the recorded states and measurements.\n"
	      << "\t\t result_filename:\t The filename prefix where to record the results as a space-separated values file.\n"
	      << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
	      << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix."
	      << "\t\t added_R.xml:\t\t Measurement noise covariance matrix to be artificially added to the data points."
	      << "\t\t test_count:\t\t Count of the number of test files to process." << std::endl;
    return 0;
  };
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  std::string inertia_filename(argv[1]);
  std::string meas_filename(argv[2]);
  std::string result_filename(argv[3]);
  
  double time_step = 0.01;
  
  RK_NOTICE(1," reached!");
  
  std::string Qu_filename(argv[4]);
  std::string R_filename(argv[5]);
  std::string R_added_filename(argv[6]);
  
  std::size_t test_count = 0;
  std::stringstream(argv[7]) >> test_count;
  
  RK_NOTICE(1," reached!");
  
  double mass;
  mat<double,mat_structure::symmetric> inertia_tensor;
  try {
    serialization::xml_iarchive in(inertia_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("mass",mass)
       & RK_SERIAL_LOAD_WITH_ALIAS("inertia_tensor",inertia_tensor);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the inertial data!");
    return 6;
  };
  
  RK_NOTICE(1," reached!");
  
  /* input disturbance */
  mat<double,mat_structure::diagonal> Qu;
  try {
    serialization::xml_iarchive in(Qu_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("input_disturbance",Qu);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return 2;
  };
  
  RK_NOTICE(1," reached!");
  
  /* measurement noise */
  mat<double,mat_structure::diagonal> R;
  try {
    serialization::xml_iarchive in(R_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("measurement_noise",R);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return 3;
  };
  
  RK_NOTICE(1," reached!");
  
  mat<double,mat_structure::diagonal> R_added;
  try {
    serialization::xml_iarchive in(R_added_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("artificial_noise",R_added);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the artificial measurement noise covariance matrix!");
    return 7;
  };
  
  RK_NOTICE(1," reached!");
  
  R += R_added;
  
  mat<double,mat_structure::diagonal> Qu_avg = Qu;
  
  RK_NOTICE(1," reached!");
  
  vect_n<double> x_init(13);
  x_init[0] = 0.0; x_init[1] = 0.0; x_init[2] = 0.0; 
  x_init[3] = 1.0; x_init[4] = 0.0; x_init[5] = 0.0; x_init[6] = 0.0;
  x_init[7] = 0.0; x_init[8] = 0.0; x_init[9] = 0.0; 
  x_init[10] = 0.0; x_init[11] = 0.0; x_init[12] = 0.0;
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
    b_init(x_init,
           ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(13,10.0))));
  RK_NOTICE(1," reached!");
  ctrl::covariance_matrix< vect_n<double> > Rcov = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R));
  
  RK_NOTICE(1," reached!");
  recorder::ssv_recorder results(result_filename + "_stddevs.ssv");
  results << "time_step" << "w_avg" << "std_p" << "std_q"
          << "EKF_p" << "EKF_q" << "EKF_w"
          << "IEKF_p" << "IEKF_q" << "IEKF_w"
          << "IMKF_p" << "IMKF_q" << "IMKF_w"
          << "IMKFv2_p" << "IMKFv2_q" << "IMKFv2_w"
	  << recorder::data_recorder::end_name_row;
  
  std::cout << "Starting estimation..." << std::endl;
  
  for(std::size_t test_num = 0; test_num < test_count; ++test_num) {
    
    std::cout << "\nTest " << std::setw(5) << test_num; std::cout.flush();
    
    std::list< std::pair< double, vect_n<double> > > measurements;
    std::list< std::pair< double, vect_n<double> > > measurements_noisy;
    double w_avg = 0.0;
    {
      std::stringstream ss_tmp;
      ss_tmp << test_num << ".ssv";
      recorder::ssv_extractor meas_file(meas_filename + ss_tmp.str());
      //std::cout << " loading file: '" << ss_tmp.str() << "'" << std::endl;
      try {
        while(true) {
	  double t;
	  meas_file >> t;
	  //std::cout << " time: " << t << std::endl;
          std::vector<double> meas;
	  try {
	    while(true) {
	      double dummy;
	      meas_file >> dummy;
	      //std::cout << " value: " << dummy << std::endl;
	      meas.push_back(dummy);
	    };
	  } catch(recorder::out_of_bounds&) { };
	  if(meas.size() < 26) {
	    RK_ERROR("The measurement file does not appear to have the required number of columns!");
            return 4;
          };
	  meas_file >> recorder::data_extractor::end_value_row;
	  w_avg *= measurements.size();
	  measurements.push_back(std::make_pair(t,vect_n<double>(meas.begin(),meas.end())));
	  w_avg += std::sqrt( meas[7] * meas[7] + meas[8] * meas[8] + meas[9] * meas[9] );
	  w_avg /= measurements.size();
	  //std::cout << " size: " << measurements.size() << std::endl;
	  vect_n<double> meas_v(7);
	  meas_v[0] = meas[4] + var_rnd() * std::sqrt(R_added(0,0));
	  meas_v[1] = meas[5] + var_rnd() * std::sqrt(R_added(1,1));
	  meas_v[2] = meas[6] + var_rnd() * std::sqrt(R_added(2,2));
	  meas_v[3] = meas[0] + var_rnd() * std::sqrt(R_added(3,3));
	  meas_v[4] = meas[1] + var_rnd() * std::sqrt(R_added(4,4));
	  meas_v[5] = meas[2] + var_rnd() * std::sqrt(R_added(5,5));
	  meas_v[6] = meas[3] + var_rnd() * std::sqrt(R_added(6,6));
	  measurements_noisy.push_back(std::make_pair(t,meas_v));
        };
      } catch(recorder::out_of_bounds&) {
        RK_ERROR("The measurement file does not appear to have the required number of columns!");
        return 4;
      } catch(recorder::end_of_record&) { };
    };
    
    std::cout << " data acquired " << std::endl;
    
    for(unsigned int j = 1; j <= 10; ++j) {
      
      results << (time_step * j) << w_avg;
      
      {
	vect<double,2> std_dev = vect<double,2>(0.0,0.0); int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
	  quaternion<double> q_mean(vect<double,4>(it->second[3],it->second[4],it->second[5],it->second[6]));
          quaternion<double> q_true(vect<double,4>(it_orig->second[0],it_orig->second[1],it_orig->second[2],it_orig->second[3]));
	  axis_angle<double> aa_diff( invert(q_true) * q_mean );
     
          std_dev[0] = ( k * std_dev[0] + ( (it->second[0] - it_orig->second[4]) * (it->second[0] - it_orig->second[4])
                                          + (it->second[1] - it_orig->second[5]) * (it->second[1] - it_orig->second[5])
                                          + (it->second[2] - it_orig->second[6]) * (it->second[2] - it_orig->second[6]) ) ) / (k + 1);
	  std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
	  ++k;

          for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
      
        };
        std_dev[0] = std::sqrt(std_dev[0]);
        std_dev[1] = std::sqrt(std_dev[1]);
        results << std_dev[0] << std_dev[1];
      };
      
      std::cout << "\r" << std::setw(3) << j; std::cout.flush();
      
      Qu_avg = (1.0 / double(j)) * Qu;

      ctrl::airship3D_lin_dt_system mdl_lin_dt("airship3D_linear_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_dt_system mdl_inv_dt("airship3D_invariant_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_mom_dt_system mdl_inv_mom_dt("airship3D_invariant_momentum_discrete",mass,inertia_tensor,time_step * j);
      ctrl::airship3D_inv_mid_dt_system mdl_inv_mid_dt("airship3D_invariant_midpoint_discrete",mass,inertia_tensor,time_step * j);
      pp::vector_topology< vect_n<double> > mdl_state_space;
      
      //std::cout << "Running Extended Kalman Filter..." << std::endl;
      {
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b = b_init;
        ctrl::covariance_matrix< vect_n<double> > Qcov;
	Qcov = Qu_avg;
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											           Qcov);
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											           Rcov);
        vect<double,3> std_dev = vect<double,3>(0.0,0.0,0.0); int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
          b_z.set_mean_state(it->second);
          ctrl::kalman_filter_step(mdl_lin_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
    
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
      
          const vect_n<double>& x_mean = b.get_mean_state();
	  
	  quaternion<double> q_true(vect<double,4>(it_orig->second[0],it_orig->second[1],it_orig->second[2],it_orig->second[3]));
	  axis_angle<double> aa_diff( invert(q_true) * q_mean );
     
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[4]) * (x_mean[0] - it_orig->second[4])
                                          + (x_mean[1] - it_orig->second[5]) * (x_mean[1] - it_orig->second[5])
                                          + (x_mean[2] - it_orig->second[6]) * (x_mean[2] - it_orig->second[6]) ) ) / (k + 1);
	  std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
	  vect<double,3> w_true(it_orig->second[7],it_orig->second[8],it_orig->second[9]);
	  w_true = invert(q_true) * w_true;
          std_dev[2] = ( k * std_dev[2] + ( (x_mean[10] - w_true[0]) * (x_mean[10] - w_true[0])
                                          + (x_mean[11] - w_true[1]) * (x_mean[11] - w_true[1])
                                          + (x_mean[12] - w_true[2]) * (x_mean[12] - w_true[2]) ) ) / (k + 1);
          ++k;

          for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
      
        };
        std_dev[0] = std::sqrt(std_dev[0]);
        std_dev[1] = std::sqrt(std_dev[1]);
        std_dev[2] = std::sqrt(std_dev[2]);
        results << std_dev[0] << std_dev[1] << std_dev[2];
      };
      std::cout << " EKF "; std::cout.flush();
  
      //std::cout << "Running Invariant Extended Kalman Filter..." << std::endl;
      {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
  
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
        
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											           ctrl::covariance_matrix< vect_n<double> >(Qu_avg));
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											           Rcovinv);
      
        vect<double,3> std_dev = vect<double,3>(0.0,0.0,0.0); int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
    
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
    
          const vect_n<double>& x_mean = b.get_mean_state();
	  
	  quaternion<double> q_true(vect<double,4>(it_orig->second[0],it_orig->second[1],it_orig->second[2],it_orig->second[3]));
	  axis_angle<double> aa_diff( invert(q_true) * q_mean );
     
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[4]) * (x_mean[0] - it_orig->second[4])
                                          + (x_mean[1] - it_orig->second[5]) * (x_mean[1] - it_orig->second[5])
                                          + (x_mean[2] - it_orig->second[6]) * (x_mean[2] - it_orig->second[6]) ) ) / (k + 1);
	  std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
	  vect<double,3> w_true(it_orig->second[7],it_orig->second[8],it_orig->second[9]);
	  w_true = invert(q_true) * w_true;
          std_dev[2] = ( k * std_dev[2] + ( (x_mean[10] - w_true[0]) * (x_mean[10] - w_true[0])
                                          + (x_mean[11] - w_true[1]) * (x_mean[11] - w_true[1])
                                          + (x_mean[12] - w_true[2]) * (x_mean[12] - w_true[2]) ) ) / (k + 1);
          ++k;

          for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
      
        };
        std_dev[0] = std::sqrt(std_dev[0]);
        std_dev[1] = std::sqrt(std_dev[1]);
        std_dev[2] = std::sqrt(std_dev[2]);
        results << std_dev[0] << std_dev[1] << std_dev[2];
      };
      std::cout << " IEKF "; std::cout.flush();
  
    //std::cout << "Running Invariant-Momentum Kalman Filter..." << std::endl;
      {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
    
        mat<double,mat_structure::diagonal> R_inv(6);
        R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
        R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
        ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
      
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											           ctrl::covariance_matrix< vect_n<double> >(Qu_avg));
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											           Rcovinv);
      
        vect<double,3> std_dev = vect<double,3>(0.0,0.0,0.0); int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
          b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_mom_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
    
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
    
          const vect_n<double>& x_mean = b.get_mean_state();
	  
	  quaternion<double> q_true(vect<double,4>(it_orig->second[0],it_orig->second[1],it_orig->second[2],it_orig->second[3]));
	  axis_angle<double> aa_diff( invert(q_true) * q_mean );
     
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[4]) * (x_mean[0] - it_orig->second[4])
                                          + (x_mean[1] - it_orig->second[5]) * (x_mean[1] - it_orig->second[5])
                                          + (x_mean[2] - it_orig->second[6]) * (x_mean[2] - it_orig->second[6]) ) ) / (k + 1);
	  std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
	  vect<double,3> w_true(it_orig->second[7],it_orig->second[8],it_orig->second[9]);
	  w_true = invert(q_true) * w_true;
          std_dev[2] = ( k * std_dev[2] + ( (x_mean[10] - w_true[0]) * (x_mean[10] - w_true[0])
                                          + (x_mean[11] - w_true[1]) * (x_mean[11] - w_true[1])
                                          + (x_mean[12] - w_true[2]) * (x_mean[12] - w_true[2]) ) ) / (k + 1);
          ++k;

          for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
    
        };
        std_dev[0] = std::sqrt(std_dev[0]);
        std_dev[1] = std::sqrt(std_dev[1]);
        std_dev[2] = std::sqrt(std_dev[2]);
        results << std_dev[0] << std_dev[1] << std_dev[2];
      };
      std::cout << " IMKF "; std::cout.flush();
  
      //std::cout << "Running Invariant-Midpoint Kalman Filter..." << std::endl;
      {
    
        ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
          b(b_init.get_mean_state(),
            ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
  
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
        
        vect<double,3> std_dev = vect<double,3>(0.0,0.0,0.0); int k = 0;
        std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
        for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
          
	  b_z.set_mean_state(it->second);
          ctrl::invariant_kalman_filter_step(mdl_inv_mid_dt,mdl_state_space,b,b_u,b_z,it->first - time_step * j);
          
          vect_n<double> b_mean = b.get_mean_state();
          quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
          b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
          b.set_mean_state(b_mean);
          
          const vect_n<double>& x_mean = b.get_mean_state();
	  
	  quaternion<double> q_true(vect<double,4>(it_orig->second[0],it_orig->second[1],it_orig->second[2],it_orig->second[3]));
	  axis_angle<double> aa_diff( invert(q_true) * q_mean );
     
          std_dev[0] = ( k * std_dev[0] + ( (x_mean[0] - it_orig->second[4]) * (x_mean[0] - it_orig->second[4])
                                          + (x_mean[1] - it_orig->second[5]) * (x_mean[1] - it_orig->second[5])
                                          + (x_mean[2] - it_orig->second[6]) * (x_mean[2] - it_orig->second[6]) ) ) / (k + 1);
	  std_dev[1] = ( k * std_dev[1] + aa_diff.angle() * aa_diff.angle() ) / (k + 1);
	  vect<double,3> w_true(it_orig->second[7],it_orig->second[8],it_orig->second[9]);
	  w_true = invert(q_true) * w_true;
          std_dev[2] = ( k * std_dev[2] + ( (x_mean[10] - w_true[0]) * (x_mean[10] - w_true[0])
                                          + (x_mean[11] - w_true[1]) * (x_mean[11] - w_true[1])
                                          + (x_mean[12] - w_true[2]) * (x_mean[12] - w_true[2]) ) ) / (k + 1);
          ++k;
          
          for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
          
        };
        std_dev[0] = std::sqrt(std_dev[0]);
        std_dev[1] = std::sqrt(std_dev[1]);
        std_dev[2] = std::sqrt(std_dev[2]);
        results << std_dev[0] << std_dev[1] << std_dev[2];
      };
      std::cout << " IMKFv2 "; std::cout.flush();
      
      results << recorder::data_recorder::end_value_row;
    };
  };

  results << recorder::data_recorder::flush;

};









