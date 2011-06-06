

#include "airship2D_lin_model.hpp"

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

int main(int argc, char** argv) {
  using namespace ReaK;
  
  if(argc < 6) {
    std::cout << "Usage:\n"
	      << "\t./estimate_airship2D [meas_filename.ssv] [result_filename] [time_step] [Qu.xml] [R.xml]\n"
	      << "\t\t meas_filename.ssv:\t The filename of a space-sep. values file with the recorded states and measurements.\n"
	      << "\t\t result_filename:\t The filename prefix where to record the results as a space-separated values file.\n"
	      << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
	      << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix." << std::endl;
    return 0;
  };
  
  std::string meas_filename(argv[1]);
  std::string result_filename(argv[2]);
  
  double time_step = 0.001;
  std::stringstream(argv[3]) >> time_step;
  
  std::string Qu_filename(argv[4]);
  std::string R_filename(argv[5]);
  
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
  
  
  std::list< std::pair< double, vect_n<double> > > measurements;
  {
    recorder::ssv_extractor meas_file(meas_filename);
    try {
      while(true) {
	double t;
	meas_file >> t;
	for(unsigned int i = 0; i < 3; ++i) {
	  double dummy;
	  meas_file >> dummy;
	};
	vect_n<double> meas(4);
	for(unsigned int i = 0; i < 4; ++i)
	  meas_file >> meas[i];
	meas_file >> recorder::data_extractor::end_value_row;
	measurements.push_back(std::make_pair(t,meas));
      };
    } catch(recorder::out_of_bounds& e) {
      RK_ERROR("The measurement file does not appear to have the required number of columns!");
      return 4;
    } catch(recorder::end_of_record&) { }
  };
  
  
  ctrl::airship2D_lin_system mdl_lin("airship2D_linear",1.0,1.0);
  ctrl::airship2D_inv_system mdl_inv("airship2D_invariant",1.0,1.0);
  ctrl::airship2D_lin_dt_system mdl_lin_dt("airship2D_linear_discrete",1.0,1.0,time_step);
  ctrl::airship2D_inv_dt_system mdl_inv_dt("airship2D_invariant_discrete",1.0,1.0,time_step);
  
  
  boost::posix_time::ptime t1;
  boost::posix_time::time_duration dt[5];
  
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > 
    b_init(vect_n<double>(0.0,0.0,1.0,0.0,0.0,0.0,0.0),
           ctrl::covariance_matrix<double>(ctrl::covariance_matrix<double>::matrix_type(mat<double,mat_structure::diagonal>(7,10.0))));
  
  ctrl::covariance_matrix<double> Rcov = ctrl::covariance_matrix<double>(ctrl::covariance_matrix<double>::matrix_type(R));
    
  euler_integrator<double> integ;
  integ.setStepSize(0.001 * time_step);

#if 0
  std::cout << "Running Kalman-Bucy Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  recorder::ssv_recorder results(result_filename + "_kbf.ssv");
  results << "time" << "pos_x" << "pos_y" << "cos(a)" << "sin(a)" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship2D_lin_system::matrixA_type A;
    ctrl::airship2D_lin_system::matrixB_type B;
    ctrl::airship2D_lin_system::matrixC_type C;
    ctrl::airship2D_lin_system::matrixD_type D;
    mdl_lin.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::kalman_bucy_filter_step(mdl_lin,integ,b,vect_n<double>(0.0,0.0,0.0),it->second,Qcov,Rcov,time_step,it->first);
    
    results << it->first << b.get_mean_state()[0] << b.get_mean_state()[1] << b.get_mean_state()[2] << b.get_mean_state()[3] << recorder::data_recorder::end_value_row;
    
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
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  ctrl::covariance_matrix<double> RcovInvar = ctrl::covariance_matrix<double>(ctrl::covariance_matrix<double>::matrix_type( mat_const_sub_sym_block< mat<double,mat_structure::diagonal> >(R,3,0) ));
  recorder::ssv_recorder results(result_filename + "_ikbf.ssv");
  results << "time" << "pos_x" << "pos_y" << "cos(a)" << "sin(a)" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship2D_inv_system::matrixA_type A;
    ctrl::airship2D_inv_system::matrixB_type B;
    ctrl::airship2D_inv_system::matrixC_type C;
    ctrl::airship2D_inv_system::matrixD_type D;
    mdl_inv.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::invariant_kalman_bucy_filter_step(mdl_inv,integ,b,vect_n<double>(0.0,0.0,0.0),it->second,Qcov,RcovInvar,time_step,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    vect<double,2> tmp = unit(vect<double,2>(b_mean[2],b_mean[3]));
    b_mean[2] = tmp[0]; b_mean[3] = tmp[1];
    b.set_mean_state(b_mean);
    
    results << it->first << b.get_mean_state()[0] << b.get_mean_state()[1] << b.get_mean_state()[2] << b.get_mean_state()[3] << recorder::data_recorder::end_value_row;
    
    std::cout << "\r" << std::setw(20) << it->first 
                      << std::setw(20) << b.get_mean_state()[0] 
                      << std::setw(20) << b.get_mean_state()[1] 
                      << std::setw(20) << b.get_mean_state()[2] 
                      << std::setw(20) << b.get_mean_state()[3];
  };
  results << recorder::data_recorder::flush;
  dt[1] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 1
  std::cout << "Running Extended Kalman Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  recorder::ssv_recorder results(result_filename + "_ekf.ssv");
  results << "time" << "pos_x" << "pos_y" << "cos(a)" << "sin(a)" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship2D_lin_dt_system::matrixA_type A;
    ctrl::airship2D_lin_dt_system::matrixB_type B;
    ctrl::airship2D_lin_dt_system::matrixC_type C;
    ctrl::airship2D_lin_dt_system::matrixD_type D;
    mdl_lin_dt.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::kalman_filter_step(mdl_lin_dt,b,vect_n<double>(0.0,0.0,0.0),it->second,Qcov,Rcov,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    vect<double,2> tmp = unit(vect<double,2>(b_mean[2],b_mean[3]));
    b_mean[2] = tmp[0]; b_mean[3] = tmp[1];
    b.set_mean_state(b_mean);
    
    results << it->first << b.get_mean_state()[0] << b.get_mean_state()[1] << b.get_mean_state()[2] << b.get_mean_state()[3] << recorder::data_recorder::end_value_row;
    
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
  results << "time" << "pos_x" << "pos_y" << "cos(a)" << "sin(a)" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship2D_lin_dt_system::matrixA_type A;
    ctrl::airship2D_lin_dt_system::matrixB_type B;
    ctrl::airship2D_lin_dt_system::matrixC_type C;
    ctrl::airship2D_lin_dt_system::matrixD_type D;
    mdl_lin_dt.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::unscented_kalman_filter_step(mdl_lin_dt,b,vect_n<double>(0.0,0.0,0.0),it->second,Qcov,Rcov,it->first);
    
    results << it->first << b.get_mean_state()[0] << b.get_mean_state()[1] << b.get_mean_state()[2] << b.get_mean_state()[3] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[3] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
#if 1
  std::cout << "Running Invariant Extended Kalman Filter..." << std::endl;
  {
  ctrl::gaussian_belief_state< ctrl::covariance_matrix<double> > b = b_init;
  recorder::ssv_recorder results(result_filename + "_iekf.ssv");
  results << "time" << "pos_x" << "pos_y" << "cos(a)" << "sin(a)" << recorder::data_recorder::end_name_row;
  t1 = boost::posix_time::microsec_clock::local_time();
  for(std::list< std::pair< double, vect_n<double> > >::iterator it = measurements.begin(); it != measurements.end(); ++it) {
    ctrl::airship2D_inv_dt_system::matrixA_type A;
    ctrl::airship2D_inv_dt_system::matrixB_type B;
    ctrl::airship2D_inv_dt_system::matrixC_type C;
    ctrl::airship2D_inv_dt_system::matrixD_type D;
    mdl_inv_dt.get_linear_blocks(A,B,C,D,it->first,b.get_mean_state(),vect_n<double>(0.0,0.0,0.0));
    ctrl::covariance_matrix<double> Qcov(ctrl::covariance_matrix<double>::matrix_type( B * Qu * transpose(B) ));
    
    ctrl::invariant_kalman_filter_step(mdl_inv_dt,b,vect_n<double>(0.0,0.0,0.0),it->second,Qcov,Rcov,it->first);
    
    vect_n<double> b_mean = b.get_mean_state();
    vect<double,2> tmp = unit(vect<double,2>(b_mean[2],b_mean[3]));
    b_mean[2] = tmp[0]; b_mean[3] = tmp[1];
    b.set_mean_state(b_mean);
    
    results << it->first << b.get_mean_state()[0] << b.get_mean_state()[1] << b.get_mean_state()[2] << b.get_mean_state()[3] << recorder::data_recorder::end_value_row;
    
  };
  results << recorder::data_recorder::flush;
  dt[4] = boost::posix_time::microsec_clock::local_time() - t1;
  };
  std::cout << "Done." << std::endl;
#endif
  
  {
  recorder::ssv_recorder results(result_filename + "_times.ssv");
  results << "step_count" << "kbf(ms)" << "ikbf(ms)" << "ekf(ms)" << "ukf(ms)" << "iekf(ms)" << recorder::data_recorder::end_name_row;
  results << double(measurements.size()) << dt[0].total_milliseconds() 
                                         << dt[1].total_milliseconds() 
					 << dt[2].total_milliseconds() 
					 << dt[3].total_milliseconds() 
					 << dt[4].total_milliseconds() 
					 << recorder::data_recorder::end_value_row << recorder::data_recorder::flush;
  };
  
  
};






