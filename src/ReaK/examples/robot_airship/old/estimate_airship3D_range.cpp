
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

int main(int argc, char** argv) {
  using namespace ReaK;
  
  if(argc < 10) {
    std::cout << "Usage:\n"
              << "\t./estimate_airship3D [inertial_data.xml] [meas_filename.ssv] [result_filename] [time_step] [Qu.xml] [R.xml] [skips_min] [skips_max] [added_R.xml]\n"
              << "\t\t inertial_data.xml:\t The filename of the inertial data of the airship3D.\n"
              << "\t\t meas_filename.ssv:\t The filename of a space-sep. values file with the recorded states and measurements.\n"
              << "\t\t result_filename:\t The filename prefix where to record the results as a space-separated values file.\n"
              << "\t\t time_step:\t\t The time-step of the data points in the measurement file.\n"
              << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
              << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix."
              << "\t\t skips_min:\t\t Minimum number of data rows to skip from meas_filename.ssv."
              << "\t\t skips_max:\t\t Maximum number of data rows to skip from meas_filename.ssv."
              << "\t\t added_R.xml:\t\t Measurement noise covariance matrix to be artificially added to the data points." << std::endl;
    return 0;
  };
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  std::string inertia_filename(argv[1]);
  std::string meas_filename(argv[2]);
  std::string result_filename(argv[3]);
  
  double time_step = 0.001;
  std::stringstream(argv[4]) >> time_step;
  
  std::string Qu_filename(argv[5]);
  std::string R_filename(argv[6]);
  
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
  
  /* input disturbance */
  mat<double,mat_structure::diagonal> Qu;
  try {
    serialization::xml_iarchive in(Qu_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("input_disturbance",Qu);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return 2;
  };
  
  /* measurement noise */
  mat<double,mat_structure::diagonal> R;
  try {
    serialization::xml_iarchive in(R_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("measurement_noise",R);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return 3;
  };
  
  unsigned int skips_min = 0;
  std::stringstream(argv[7]) >> skips_min;
  ++skips_min;
  
  unsigned int skips_max = 0;
  std::stringstream(argv[8]) >> skips_max;
  ++skips_max;

  mat<double,mat_structure::diagonal> R_added(7,0.0);
  try {
    std::string R_added_filename(argv[9]);
    serialization::xml_iarchive in(R_added_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("artificial_noise",R_added);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the artificial measurement noise covariance matrix!");
    return 3;
  };
  
  
  std::list< std::pair< double, vect_n<double> > > measurements;
  std::list< std::pair< double, vect_n<double> > > measurements_noisy;
  {
    recorder::ssv_extractor meas_file(meas_filename);
    try {
      while(true) {
        double t;
        meas_file >> t;
        std::vector<double> meas;
        try {
          while(true) {
            double dummy;
            meas_file >> dummy;
            meas.push_back(dummy);
          };
        } catch(recorder::out_of_bounds&) { };
        if(meas.size() < 7) {
          RK_ERROR("The measurement file does not appear to have the required number of columns!");
          return 4;
        };
        meas_file >> recorder::data_extractor::end_value_row;
        vect_n<double> meas_v(meas.end() - 7,meas.end());
        measurements.push_back(std::make_pair(t,meas_v));
        meas_v[0] += var_rnd() * sqrt(R_added(0,0));
        meas_v[1] += var_rnd() * sqrt(R_added(1,1));
        meas_v[2] += var_rnd() * sqrt(R_added(2,2));
        meas_v[3] += var_rnd() * sqrt(R_added(3,3));
        meas_v[4] += var_rnd() * sqrt(R_added(4,4));
        meas_v[5] += var_rnd() * sqrt(R_added(5,5));
        meas_v[6] += var_rnd() * sqrt(R_added(6,6));
        measurements_noisy.push_back(std::make_pair(t,meas_v));
      };
    } catch(recorder::out_of_bounds&) {
      RK_ERROR("The measurement file does not appear to have the required number of columns!");
      return 4;
    } catch(recorder::end_of_record&) { }
  };
  
  R += R_added;
  
  vect_n<double> x_init(13);
  x_init[0] = 0.0; x_init[1] = 0.0; x_init[2] = 0.0; 
  x_init[3] = 1.0; x_init[4] = 0.0; x_init[5] = 0.0; x_init[6] = 0.0;
  x_init[7] = 0.0; x_init[8] = 0.0; x_init[9] = 0.0; 
  x_init[10] = 0.0; x_init[11] = 0.0; x_init[12] = 0.0;
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
    b_init(x_init,
           ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(13,10.0))));
  
  ctrl::covariance_matrix< vect_n<double> > Rcov = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R));
  
  recorder::ssv_recorder results(result_filename + "_stddevs_quat.ssv");
  results << "time_step" << "ekf" << "iekf" << "imkf" << recorder::data_recorder::end_name_row;
  
  for(unsigned int j = skips_min+1; j<=skips_max+1; ++j) {

    ctrl::airship3D_lin_dt_system mdl_lin_dt("airship3D_linear_discrete",mass,inertia_tensor,time_step * j);
    ctrl::airship3D_inv_dt_system mdl_inv_dt("airship3D_invariant_discrete",mass,inertia_tensor,time_step * j);
    ctrl::airship3D_inv_mom_dt_system mdl_inv_mom_dt("airship3D_invariant_momentum_discrete",mass,inertia_tensor,time_step * j);
    pp::vector_topology< vect_n<double> > mdl_state_space;
    
    results << (time_step * j);

    //std::cout << "Running Extended Kalman Filter..." << std::endl;
    {
      ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b = b_init;
      ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                   ctrl::covariance_matrix< vect_n<double> >(Qu));
      ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                   Rcov);
      double std_dev = 0.0; int k = 0;
      std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
      for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
        b_z.set_mean_state(it->second);
        ctrl::kalman_filter_step(mdl_lin_dt,mdl_state_space,b,b_u,b_z,it->first);
    
        vect_n<double> b_mean = b.get_mean_state();
        quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
        b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
        b.set_mean_state(b_mean);
    
        const vect_n<double>& x_mean = b.get_mean_state();
     
        std_dev = ( k * std_dev + ( (x_mean[3] - it_orig->second[3]) * (x_mean[3] - it_orig->second[3])
                                  + (x_mean[4] - it_orig->second[4]) * (x_mean[4] - it_orig->second[4])
                                  + (x_mean[5] - it_orig->second[5]) * (x_mean[5] - it_orig->second[5])
                                  + (x_mean[6] - it_orig->second[6]) * (x_mean[6] - it_orig->second[6]) ) ) / (k + 1);
        ++k;

        for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
      
      };
      std_dev = std::sqrt(std_dev);
      results << std_dev;
    };
  
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
                                                                                                   ctrl::covariance_matrix< vect_n<double> >(Qu));
      ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                   Rcovinv);
      
      double std_dev = 0.0; int k = 0;
      std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
      for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
        b_z.set_mean_state(it->second);
        ctrl::invariant_kalman_filter_step(mdl_inv_dt,mdl_state_space,b,b_u,b_z,it->first);
    
        vect_n<double> b_mean = b.get_mean_state();
        quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
        b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
        b.set_mean_state(b_mean);
    
        const vect_n<double>& x_mean = b.get_mean_state();
      
        std_dev = ( k * std_dev + ( (x_mean[3] - it_orig->second[3]) * (x_mean[3] - it_orig->second[3])
                                  + (x_mean[4] - it_orig->second[4]) * (x_mean[4] - it_orig->second[4])
                                  + (x_mean[5] - it_orig->second[5]) * (x_mean[5] - it_orig->second[5])
                                  + (x_mean[6] - it_orig->second[6]) * (x_mean[6] - it_orig->second[6]) ) ) / (k + 1);
        ++k;

        for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
      
      };
      std_dev = std::sqrt(std_dev);
      results << std_dev;
    };
  
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
                                                                                                   ctrl::covariance_matrix< vect_n<double> >(Qu));
      ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
                                                                                                   Rcovinv);
      
      double std_dev = 0.0; int k = 0;
      std::list< std::pair< double, vect_n<double> > >::iterator it_orig = measurements.begin();
      for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements_noisy.begin(); it != measurements_noisy.end();) {
        
        b_z.set_mean_state(it->second);
        ctrl::invariant_kalman_filter_step(mdl_inv_mom_dt,mdl_state_space,b,b_u,b_z,it->first);
    
        vect_n<double> b_mean = b.get_mean_state();
        quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
        b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
        b.set_mean_state(b_mean);
    
        const vect_n<double>& x_mean = b.get_mean_state();
      
        std_dev = ( k * std_dev + ( (x_mean[3] - it_orig->second[3]) * (x_mean[3] - it_orig->second[3])
                                  + (x_mean[4] - it_orig->second[4]) * (x_mean[4] - it_orig->second[4])
                                  + (x_mean[5] - it_orig->second[5]) * (x_mean[5] - it_orig->second[5])
                                  + (x_mean[6] - it_orig->second[6]) * (x_mean[6] - it_orig->second[6]) ) ) / (k + 1);
        ++k;

        for(unsigned int i = 0; ((i < j) && (it != measurements_noisy.end())); ++i) { ++it; ++it_orig; };
    
      };
      std_dev = std::sqrt(std_dev);
      results << std_dev;
    };
    results << recorder::data_recorder::end_value_row;
  };

  results << recorder::data_recorder::flush;

};









