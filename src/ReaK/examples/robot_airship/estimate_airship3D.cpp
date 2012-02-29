
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
  
  if(argc < 7) {
    std::cout << "Usage:\n"
	      << "\t./estimate_airship3D [inertial_data.xml] [meas_filename.ssv] [result_filename] [time_step] [Qu.xml] [R.xml] [skips] [added_R.xml]\n"
	      << "\t\t inertial_data.xml:\t The filename of the inertial data of the airship3D.\n"
	      << "\t\t meas_filename.ssv:\t The filename of a space-sep. values file with the recorded states and measurements.\n"
	      << "\t\t result_filename:\t The filename prefix where to record the results as a space-separated values file.\n"
	      << "\t\t time_step:\t\t The time-step of the data points in the measurement file.\n"
	      << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
	      << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix."
	      << "\tOptional:"
	      << "\t\t skips:\t\t Number of data rows to skip from meas_filename.ssv."
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
  
  unsigned int skips = 0;
  if(argc >= 8) {
    std::stringstream(argv[7]) >> skips;
  };
  ++skips;
  
  mat<double,mat_structure::diagonal> R_added(7,0.0);
  try {
    if(argc >= 9) {
      std::cout << "address 1" << std::endl;
      std::string R_added_filename(argv[8]);
      std::cout << "address 2" << std::endl;
      serialization::xml_iarchive in(R_added_filename);
      std::cout << "address 3" << std::endl;
      in & RK_SERIAL_LOAD_WITH_ALIAS("artificial_noise",R_added);
      std::cout << "address 4" << std::endl;
    };
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the artificial measurement noise covariance matrix!");
    return 3;
  };
  
  
  std::list< std::pair< double, vect_n<double> > > measurements;
  {
    recorder::ssv_extractor meas_file(meas_filename);
    try {
      unsigned int j = 0;
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
	if(j == 0) {
	  vect_n<double> meas_v(meas.end() - 7,meas.end());
	  meas_v[0] += var_rnd() * sqrt(R_added(0,0));
	  meas_v[1] += var_rnd() * sqrt(R_added(1,1));
	  meas_v[2] += var_rnd() * sqrt(R_added(2,2));
	  meas_v[3] += var_rnd() * sqrt(R_added(3,3));
	  meas_v[4] += var_rnd() * sqrt(R_added(4,4));
	  meas_v[5] += var_rnd() * sqrt(R_added(5,5));
	  meas_v[6] += var_rnd() * sqrt(R_added(6,6));
	  measurements.push_back(std::make_pair(t,meas_v));
	};
	j = (j+1) % skips;
      };
    } catch(recorder::out_of_bounds&) {
      RK_ERROR("The measurement file does not appear to have the required number of columns!");
      return 4;
    } catch(recorder::end_of_record&) { }
  };
  
  R += R_added;
  
  ctrl::airship3D_lin_system mdl_lin("airship3D_linear",mass,inertia_tensor);
  ctrl::airship3D_inv_system mdl_inv("airship3D_invariant",mass,inertia_tensor);
  ctrl::airship3D_lin_dt_system mdl_lin_dt("airship3D_linear_discrete",mass,inertia_tensor,time_step);
  ctrl::airship3D_inv_dt_system mdl_inv_dt("airship3D_invariant_discrete",mass,inertia_tensor,time_step);
  ctrl::airship3D_inv_mom_dt_system mdl_inv_mom_dt("airship3D_invariant_momentum_discrete",mass,inertia_tensor,time_step);
  
  pp::vector_topology< vect_n<double> > mdl_state_space;
  
  boost::posix_time::ptime t1;
  boost::posix_time::time_duration dt[6];
  
  vect_n<double> x_init(13);
  x_init[0] = 0.0; x_init[1] = 0.0; x_init[2] = 0.0; 
  x_init[3] = 1.0; x_init[4] = 0.0; x_init[5] = 0.0; x_init[6] = 0.0;
  x_init[7] = 0.0; x_init[8] = 0.0; x_init[9] = 0.0; 
  x_init[10] = 0.0; x_init[11] = 0.0; x_init[12] = 0.0;
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
    b_init(x_init,
           ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(13,10.0))));
  
  ctrl::covariance_matrix< vect_n<double> > Rcov = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R));
    
  euler_integrator<double> integ;
  integ.setStepSize(0.001 * time_step);

#if 0
  std::cout << "Running Kalman-Bucy Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  recorder::ssv_recorder results(result_filename + "_kbf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship3D_lin_system::matrixA_type A;
    ctrl::airship3D_lin_system::matrixB_type B;
    ctrl::airship3D_lin_system::matrixC_type C;
    ctrl::airship3D_lin_system::matrixD_type D;
    mdl_lin.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::kalman_bucy_filter_step(mdl_lin,integ,b,vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0),it->second,Qcov,Rcov,time_step,it->first);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[0] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 0
  integ.setStepSize(0.0001 * time_step);
  std::cout << "Running Invariant Kalman-Bucy Filter..." << std::endl;
  {
  
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > 
    b(b_init.get_mean_state(),
      ctrl::covariance_matrix<double>(ctrl::covariance_matrix<double>::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
    
  ctrl::covariance_matrix<double> RcovInvar = ctrl::covariance_matrix<double>(ctrl::covariance_matrix<double>::matrix_type( mat_const_sub_sym_block< mat<double,mat_structure::diagonal> >(R,6,0) ));
  recorder::ssv_recorder results(result_filename + "_ikbf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship3D_inv_system::matrixA_type A;
    ctrl::airship3D_inv_system::matrixB_type B;
    ctrl::airship3D_inv_system::matrixC_type C;
    ctrl::airship3D_inv_system::matrixD_type D;
    mdl_inv.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::invariant_kalman_bucy_filter_step(mdl_inv,integ,b,vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0),it->second,Qcov,RcovInvar,time_step,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
    b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
    b.set_mean_state(b_mean);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;

  };
  results << recorder::data_recorder::flush;
  dt[1] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 1
  std::cout << "Running Extended Kalman Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b = b_init;
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											       ctrl::covariance_matrix< vect_n<double> >(Qu));
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											       Rcov);
  recorder::ssv_recorder results(result_filename + "_ekf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
     
    b_z.set_mean_state(it->second);
    ctrl::kalman_filter_step(mdl_lin_dt,mdl_state_space,b,b_u,b_z,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
    b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
    b.set_mean_state(b_mean);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[2] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 0
  std::cout << "Running Unscented Kalman Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  recorder::ssv_recorder results(result_filename + "_ukf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship3D_lin_dt_system::matrixA_type A;
    ctrl::airship3D_lin_dt_system::matrixB_type B;
    ctrl::airship3D_lin_dt_system::matrixC_type C;
    ctrl::airship3D_lin_dt_system::matrixD_type D;
    mdl_lin_dt.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    try {
      ctrl::unscented_kalman_filter_step(mdl_lin_dt,b,vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0),it->second,Qcov,Rcov,it->first);
    } catch(singularity_error& e) {
      RK_ERROR("The Unscented Kalman filtering was interupted by a singularity, at time " << it->first << " with message: " << e.what());
      break;
    } catch(std::exception& e) {
      RK_ERROR("The Unscented Kalman filtering was interupted by an exception, at time " << it->first << " with message: " << e.what());   
      break;
    };
    
    vect_n<double> b_mean = b.get_mean_state();
    quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
    b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
    b.set_mean_state(b_mean);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[3] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 1
  std::cout << "Running Invariant Extended Kalman Filter..." << std::endl;
  {
    
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
    b(b_init.get_mean_state(),
      ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
  
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											       ctrl::covariance_matrix< vect_n<double> >(Qu));
  
  mat<double,mat_structure::diagonal> R_inv(6);
  R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
  R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
  ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
  
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											       Rcovinv);
  
  recorder::ssv_recorder results(result_filename + "_iekf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    
    b_z.set_mean_state(it->second);
    ctrl::invariant_kalman_filter_step(mdl_inv_dt,mdl_state_space,b,b_u,b_z,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
    b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
    b.set_mean_state(b_mean);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[4] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 1
  std::cout << "Running Invariant-Momentum Kalman Filter..." << std::endl;
  {
    
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > 
    b(b_init.get_mean_state(),
      ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
  
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_u(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0), 
											       ctrl::covariance_matrix< vect_n<double> >(Qu));
  
  mat<double,mat_structure::diagonal> R_inv(6);
  R_inv(0,0) = R(0,0); R_inv(1,1) = R(1,1); R_inv(2,2) = R(2,2);
  R_inv(3,3) = 4*R(4,4); R_inv(4,4) = 4*R(5,5); R_inv(5,5) = 4*R(6,6);
  ctrl::covariance_matrix< vect_n<double> > Rcovinv = ctrl::covariance_matrix< vect_n<double> >(ctrl::covariance_matrix< vect_n<double> >::matrix_type(R_inv));
  
  ctrl::gaussian_belief_state< vect_n<double>, ctrl::covariance_matrix< vect_n<double> > > b_z(vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0), 
											       Rcovinv);
  
  recorder::ssv_recorder results(result_filename + "_imkf.ssv");
  results << "time" << "pos_x" << "pos_y" << "pos_z" << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    
    b_z.set_mean_state(it->second);
    ctrl::invariant_kalman_filter_step(mdl_inv_mom_dt,mdl_state_space,b,b_u,b_z,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    quaternion<double> q_mean(vect<double,4>(b_mean[3],b_mean[4],b_mean[5],b_mean[6]));
    b_mean[3] = q_mean[0]; b_mean[4] = q_mean[1]; b_mean[5] = q_mean[2]; b_mean[6] = q_mean[3];
    b.set_mean_state(b_mean);
    
    const vect_n<double>& x_mean = b.get_mean_state();
    results << it->first << x_mean[0] << x_mean[1] << x_mean[2] 
                         << x_mean[3] << x_mean[4] << x_mean[5] << x_mean[6] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[5] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
  {
  recorder::ssv_recorder results(result_filename + "_times.ssv");
  results << "step_count" << "kbf(ms)" << "ikbf(ms)" << "ekf(ms)" << "ukf(ms)" << "iekf(ms)" << "imkf(ms)" << recorder::data_recorder::end_name_row;
  results << double(measurements.size()) << double(dt[0].total_milliseconds())
                                         << double(dt[1].total_milliseconds())
					 << double(dt[2].total_milliseconds())
					 << double(dt[3].total_milliseconds())
					 << double(dt[4].total_milliseconds())
					 << double(dt[5].total_milliseconds()) 
					 << recorder::data_recorder::end_value_row << recorder::data_recorder::flush;
  };
  
  
};






