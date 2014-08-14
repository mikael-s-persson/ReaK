
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


#include <ReaK/ctrl/ss_systems/satellite_invar_models.hpp>
#include <ReaK/ctrl/ss_systems/satellite_modeling_po.hpp>

#include <ReaK/ctrl/ctrl_sys/kalman_filter.hpp>
#include <ReaK/ctrl/ctrl_sys/invariant_kalman_filter.hpp>
#include <ReaK/ctrl/ctrl_sys/tsos_aug_kalman_filter.hpp>
#include <ReaK/ctrl/ctrl_sys/tsos_aug_inv_kalman_filter.hpp>

#include <ReaK/ctrl/ctrl_sys/gaussian_belief_state.hpp>
#include <ReaK/ctrl/ctrl_sys/covariance_matrix.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>
#include <ReaK/core/recorders/data_record_po.hpp>

#include <ReaK/ctrl/topologies/temporal_space.hpp>
#include <ReaK/ctrl/interpolation/trajectory_base.hpp>
#include <ReaK/ctrl/interpolation/discrete_point_trajectory.hpp>

#include <ReaK/ctrl/ctrl_sys/augmented_to_state_mapping.hpp>

#include <boost/random/normal_distribution.hpp>
#include <ReaK/core/base/global_rng.hpp>
#include <ReaK/core/base/atomic_incl.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;



  
boost::variate_generator< ReaK::global_rng_type&, boost::normal_distribution<double> > var_rnd = 
  boost::variate_generator< ReaK::global_rng_type&, boost::normal_distribution<double> >(ReaK::get_global_rng(), boost::normal_distribution<double>());


struct sat3D_measurement_point {
  ReaK::vect_n<double> pose;
  ReaK::vect_n<double> gyro;
  ReaK::vect_n<double> IMU_a_m;
  ReaK::vect_n<double> u;
};


typedef ReaK::ctrl::satellite_model_options::state_space_type      sat3D_state_space_type;
typedef ReaK::ctrl::satellite_model_options::temp_state_space_type sat3D_temp_space_type;
typedef ReaK::ctrl::satellite_model_options::state_type            sat3D_state_type;
typedef ReaK::ctrl::satellite_model_options::input_type            sat3D_input_type;
typedef ReaK::ctrl::satellite_model_options::output_type           sat3D_output_type;
typedef ReaK::pp::topology_traits< sat3D_temp_space_type >::point_type sat3D_temp_point_type;

typedef ReaK::ctrl::satellite_model_options::covar_type cov_type;
typedef cov_type::matrix_type cov_matrix_type;
typedef ReaK::ctrl::satellite_model_options::state_belief_type  sat3D_state_belief_type;
typedef ReaK::ctrl::satellite_model_options::input_belief_type  sat3D_input_belief_type;
typedef ReaK::ctrl::satellite_model_options::output_belief_type sat3D_output_belief_type;



struct sat3D_meas_true_from_vectors {
  const std::vector< std::pair< double, sat3D_measurement_point > >* measurements;
  const std::vector< std::pair< double, sat3D_state_type > >* ground_truth;
  std::size_t skips;
  
  std::vector< std::pair< double, sat3D_measurement_point > >::const_iterator cur_meas;
  std::vector< std::pair< double, sat3D_state_type > >::const_iterator cur_true;
  
  double get_current_time() const { return cur_meas->first; };
  const sat3D_measurement_point& get_current_measurement() const { return cur_meas->second; };
  const sat3D_state_type* get_current_gnd_truth_ptr() const { 
    if(ground_truth && ground_truth->size())
      return &(cur_true->second); 
    else
      return NULL;
  };
  
  bool step_once() {
    for(std::size_t i = 0; i < skips; ++i) {
      ++cur_meas;
      if( cur_meas == measurements->end() )
        return false;
      if( ground_truth && cur_true != ground_truth->end() )
        ++cur_true;
    };
    return true;
  };
  
  sat3D_meas_true_from_vectors(const std::vector< std::pair< double, sat3D_measurement_point > >* aMeasurements,
                               const std::vector< std::pair< double, sat3D_state_type > >* aGroundTruth = NULL,
                               std::size_t aSkips = 1) :
                               measurements(aMeasurements), ground_truth(aGroundTruth), skips(aSkips) { 
    cur_meas = measurements->begin();
    if(ground_truth)
      cur_true = ground_truth->begin();
  };
  
  
};



struct sat3D_meas_true_from_extractor {
  ReaK::shared_ptr< ReaK::recorder::data_extractor > data_in;
  ReaK::ctrl::satellite_model_options sat_options;
  std::size_t skips;
  
  ReaK::recorder::named_value_row nvr_in;
  double time_val;
  sat3D_measurement_point meas_pt;
  sat3D_state_type gnd_pt;
  bool has_ground_truth;
  
  double get_current_time() const { return time_val; };
  const sat3D_measurement_point& get_current_measurement() const { return meas_pt; };
  const sat3D_state_type* get_current_gnd_truth_ptr() const { 
    if(has_ground_truth)
      return &gnd_pt; 
    else
      return NULL;
  };
  
  bool step_once() {
    using namespace ReaK;
    try {
      (*data_in) >> nvr_in;
      
      time_val = nvr_in["time"];
      
      meas_pt.pose.resize(7);
      meas_pt.pose[0] = nvr_in["p_x"];
      meas_pt.pose[1] = nvr_in["p_y"];
      meas_pt.pose[2] = nvr_in["p_z"];
      meas_pt.pose[3] = nvr_in["q_0"];
      meas_pt.pose[4] = nvr_in["q_1"];
      meas_pt.pose[5] = nvr_in["q_2"];
      meas_pt.pose[6] = nvr_in["q_3"];
      
      meas_pt.u.resize(6);
      meas_pt.u[0] = nvr_in["f_x"];
      meas_pt.u[1] = nvr_in["f_y"];
      meas_pt.u[2] = nvr_in["f_z"];
      meas_pt.u[3] = nvr_in["t_x"];
      meas_pt.u[4] = nvr_in["t_y"];
      meas_pt.u[5] = nvr_in["t_z"];
      
      try {
        meas_pt.gyro.resize(3);
        meas_pt.gyro[0] = nvr_in["w_x"];
        meas_pt.gyro[1] = nvr_in["w_y"];
        meas_pt.gyro[2] = nvr_in["w_z"];
      } catch(ReaK::recorder::out_of_bounds& e) { RK_UNUSED(e); 
        meas_pt.gyro.resize(0);
      };
      
      try {
        meas_pt.IMU_a_m.resize(6);
        meas_pt.IMU_a_m[0] = nvr_in["acc_x"];
        meas_pt.IMU_a_m[1] = nvr_in["acc_y"];
        meas_pt.IMU_a_m[2] = nvr_in["acc_z"];
        meas_pt.IMU_a_m[3] = nvr_in["mag_x"];
        meas_pt.IMU_a_m[4] = nvr_in["mag_y"];
        meas_pt.IMU_a_m[5] = nvr_in["mag_z"];
      } catch(ReaK::recorder::out_of_bounds& e) { RK_UNUSED(e); 
        meas_pt.IMU_a_m.resize(0);
      };
      
      try {
        set_position(gnd_pt, vect<double,3>(nvr_in["p_x_true"], nvr_in["p_y_true"], nvr_in["p_z_true"]));
        set_quaternion(gnd_pt, unit_quat<double>(nvr_in["q_0_true"], nvr_in["q_1_true"], nvr_in["q_2_true"], nvr_in["q_3_true"]));
        set_velocity(gnd_pt, vect<double,3>(nvr_in["v_x_true"], nvr_in["v_y_true"], nvr_in["v_z_true"]));
        set_ang_velocity(gnd_pt, vect<double,3>(nvr_in["w_x_true"], nvr_in["w_y_true"], nvr_in["w_z_true"]));
      } catch(ReaK::recorder::out_of_bounds& e) { RK_UNUSED(e); 
        set_position(gnd_pt, vect<double,3>(meas_pt.pose[0],meas_pt.pose[1],meas_pt.pose[2]));
        set_quaternion(gnd_pt, unit_quat<double>(meas_pt.pose[3],meas_pt.pose[4],meas_pt.pose[5],meas_pt.pose[6]));
        set_velocity(gnd_pt, vect<double,3>(0.0, 0.0, 0.0));
        set_ang_velocity(gnd_pt, vect<double,3>(0.0, 0.0, 0.0));
      };
      
      vect_n<double> added_noise = ctrl::sample_gaussian_point(vect_n<double>(sat_options.artificial_noise.get_row_count(), 0.0), sat_options.artificial_noise);
      if(sat_options.artificial_noise.get_row_count() < 6) 
        added_noise.resize(6, 0.0);
      
      meas_pt.pose[range(0,3)] = meas_pt.pose[range(0,3)] + added_noise[range(0,3)];
      vect<double,3> aa_noise(added_noise[3], added_noise[4], added_noise[5]);
      quaternion<double> y_quat(meas_pt.pose[range(3,6)]);
      y_quat *= axis_angle<double>(norm_2(aa_noise), aa_noise).getQuaternion();
      meas_pt.pose[range(3,7)] = vect<double,4>(y_quat[0], y_quat[1], y_quat[2], y_quat[3]);
      
      if( ( meas_pt.gyro.size() ) && ( added_noise.size() >= 9 ) )
        meas_pt.gyro += added_noise[range(6,9)];
      if( ( meas_pt.IMU_a_m.size() ) && ( added_noise.size() >= 15 ) )
        meas_pt.IMU_a_m += added_noise[range(9,15)];
      
    } catch(recorder::end_of_record&) { 
      return false;
    };
    
    return true;
  };
  
  sat3D_meas_true_from_extractor(const ReaK::shared_ptr< ReaK::recorder::data_extractor >& aDataIn,
                                 const ReaK::ctrl::satellite_model_options& aSatOptions,
                                 std::size_t aSkips = 1) :
                                 data_in(aDataIn), sat_options(aSatOptions), 
                                 skips(aSkips), nvr_in(data_in->getFreshNamedValueRow()), 
                                 time_val(0.0), has_ground_truth(false) { 
    step_once();
  };
  
  
};


namespace {
  
  static const sat3D_state_type& get_sat3D_state(const sat3D_state_type& x) { return x; };
  
  template <typename StateTuple>
  static const sat3D_state_type& get_sat3D_state(const StateTuple& x) { using ReaK::get; return get<0>(x); };
  
  
  static void set_sat3D_state(sat3D_state_type& x, const sat3D_state_type& val) { x = val; };
  
  template <typename StateTuple>
  static void set_sat3D_state(StateTuple& x, const sat3D_state_type& val) { using ReaK::get; get<0>(x) = val; };
  
};


struct sat3D_estimate_result_to_recorder {
  ReaK::shared_ptr< ReaK::recorder::data_recorder > rec;
  
  sat3D_estimate_result_to_recorder(const ReaK::shared_ptr< ReaK::recorder::data_recorder >& aRec) : rec(aRec) { };
  
  void initialize() { };
  void finalize() {
    (*rec) << ReaK::recorder::data_recorder::flush;
  };
  
  template <typename BeliefStateType, typename InputBeliefType, typename OutputBeliefType>
  void add_record(const BeliefStateType& b,
                  const InputBeliefType& b_u, 
                  const OutputBeliefType& b_z,
                  double time,
                  const sat3D_state_type* true_state = NULL) {
    using namespace ReaK;
    using ReaK::to_vect;
    
    const sat3D_state_type& x_mean = get_sat3D_state(b.get_mean_state());
    (*rec) << time 
           << get_position(x_mean) << get_quaternion(x_mean)
           << get_velocity(x_mean) << get_ang_velocity(x_mean);
    
    vect_n<double> all_x = to_vect<double>(b.get_mean_state());
    for(std::size_t l = 13; l < all_x.size(); ++l)
      (*rec) << all_x[l];
    
    if( true_state ) {
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) * get_quaternion(*true_state).as_rotation());
      (*rec) << (get_position(x_mean) - get_position(*true_state)) 
             << (aa_diff.angle() * aa_diff.axis())
             << (get_velocity(x_mean) - get_velocity(*true_state))
             << (get_ang_velocity(x_mean) - get_ang_velocity(*true_state));
    } else {
      const vect_n<double>& z = b_z.get_mean_state();
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) * quaternion<double>(z[range(3,7)]));
      (*rec) << (get_position(x_mean) - z[range(0,3)]) 
             << (aa_diff.angle() * aa_diff.axis())
             << vect<double,3>(0.0,0.0,0.0);
      if( z.size() >= 10 )
        (*rec) << (get_ang_velocity(x_mean) - z[range(7,10)]);
      else
        (*rec) << vect<double,3>(0.0,0.0,0.0);
    };
    
    const cov_matrix_type& P_xx = b.get_covariance().get_matrix();
    for(std::size_t l = 0; l < P_xx.get_row_count(); ++l)
      (*rec) << P_xx(l,l);
    
    (*rec) << recorder::data_recorder::end_value_row;
  };
  
};



struct sat3D_collect_stddevs {
  ReaK::vect_n< double > stddevs;
  std::size_t counter;
  ReaK::shared_ptr< ReaK::recorder::data_recorder > rec;
  
  sat3D_collect_stddevs(const ReaK::shared_ptr< ReaK::recorder::data_recorder >& aRec) : stddevs(aRec->getColCount(), 0.0), counter(0), rec(aRec) { };
  
  void initialize() { 
    stddevs = ReaK::vect_n< double >(rec->getColCount(), 0.0);
    counter = 0;
  };
  
  void finalize() {
    for(std::size_t j = 0; j < stddevs.size(); ++j)
      (*rec) << std::sqrt(stddevs[j]); // turn variances into std-devs.
    (*rec) << ReaK::recorder::data_recorder::end_value_row << ReaK::recorder::data_recorder::flush;
  };
  
  template <typename BeliefStateType, typename InputBeliefType, typename OutputBeliefType>
  void add_record(const BeliefStateType& b,
                  const InputBeliefType& b_u, 
                  const OutputBeliefType& b_z,
                  double time,
                  const sat3D_state_type* true_state = NULL) {
    using namespace ReaK;
    
    const sat3D_state_type& x_mean = get_sat3D_state(b.get_mean_state());
    vect<double,3> pos_err, aa_err, vel_err, ang_vel_err;
    if( true_state ) {
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) * get_quaternion(*true_state).as_rotation());
      pos_err     = get_position(x_mean) - get_position(*true_state);
      aa_err      = aa_diff.angle() * aa_diff.axis();
      vel_err     = get_velocity(x_mean) - get_velocity(*true_state);
      ang_vel_err = get_ang_velocity(x_mean) - get_ang_velocity(*true_state);
    } else {
      const vect_n<double>& z = b_z.get_mean_state();
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) * quaternion<double>(z[range(3,7)]));
      pos_err     = get_position(x_mean) - z[range(0,3)];
      aa_err      = aa_diff.angle() * aa_diff.axis();
      vel_err     = vect<double,3>(0.0,0.0,0.0);
      if( z.size() >= 10 )
        ang_vel_err = get_ang_velocity(x_mean) - z[range(7,10)];
      else
        ang_vel_err = vect<double,3>(0.0,0.0,0.0);
    };
    stddevs[range(0,3)]  = ( counter * stddevs[range(0,3)] + elem_product(pos_err, pos_err) ) / (counter + 1);
    stddevs[range(3,6)]  = ( counter * stddevs[range(3,6)] + elem_product(aa_err, aa_err) ) / (counter + 1);
    stddevs[range(6,9)]  = ( counter * stddevs[range(6,9)] + elem_product(vel_err, vel_err) ) / (counter + 1);
    stddevs[range(9,12)] = ( counter * stddevs[range(9,12)] + elem_product(ang_vel_err, ang_vel_err) ) / (counter + 1);
    
    stddevs[12] = ( counter * stddevs[12] + pos_err * pos_err ) / (counter + 1);
    stddevs[13] = ( counter * stddevs[13] + aa_err * aa_err ) / (counter + 1);
    stddevs[14] = ( counter * stddevs[14] + vel_err * vel_err ) / (counter + 1);
    stddevs[15] = ( counter * stddevs[15] + ang_vel_err * ang_vel_err ) / (counter + 1);
    
    const cov_matrix_type& P_xx = b.get_covariance().get_matrix();
    for(std::size_t l = 0; l < 12; ++l)
      stddevs[l + 16] = ( counter * stddevs[l + 16] + P_xx(l,l) ) / (counter + 1);
    
    ++counter;
  };
};




template <typename MeasureProvider, typename ResultLogger, typename Sat3DSystemType>
void batch_KF_on_timeseries(
    MeasureProvider meas_provider, 
    ResultLogger result_logger,
    const ReaK::ctrl::satellite_model_options& sat_options,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    typename Sat3DSystemType::state_belief_type b,
    typename Sat3DSystemType::input_belief_type b_u,
    typename Sat3DSystemType::output_belief_type b_z) {
  using namespace ReaK;
  
  result_logger.initialize();
  
  do {
    const sat3D_measurement_point& cur_meas = meas_provider.get_current_measurement();
    vect_n<double> z_vect(cur_meas.pose.size() + cur_meas.gyro.size() + cur_meas.IMU_a_m.size(), 0.0);
    z_vect[range(0,7)] = cur_meas.pose;
    if( cur_meas.gyro.size() ) {
      z_vect[range(7,10)] = cur_meas.gyro;
      if( cur_meas.IMU_a_m.size() )
        z_vect[range(10,16)] = cur_meas.IMU_a_m;
    };
    b_z.set_mean_state(z_vect);
    b_u.set_mean_state(cur_meas.u);
    
    if(sat_options.system_kind & ctrl::satellite_model_options::TSOSAKF) {
      ctrl::tsos_aug_inv_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    } else {
      ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    };
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() );
  
  result_logger.finalize();
};




template <typename MeasureProvider, typename ResultLogger, typename Sat3DSystemType>
void batch_KF_no_meas_predict(
    MeasureProvider meas_provider, 
    ResultLogger result_logger,
    const ReaK::ctrl::satellite_predictor_options& sat_options,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    double start_time,
    typename Sat3DSystemType::state_belief_type b,
    typename Sat3DSystemType::input_belief_type b_u,
    typename Sat3DSystemType::output_belief_type b_z) {
  using namespace ReaK;
  
  result_logger.initialize();
  
  // filtering phase:
  do {
    const sat3D_measurement_point& cur_meas = meas_provider.get_current_measurement();
    vect_n<double> z_vect(cur_meas.pose.size() + cur_meas.gyro.size() + cur_meas.IMU_a_m.size(), 0.0);
    z_vect[range(0,7)] = cur_meas.pose;
    if( cur_meas.gyro.size() ) {
      z_vect[range(7,10)] = cur_meas.gyro;
      if( cur_meas.IMU_a_m.size() )
        z_vect[range(10,16)] = cur_meas.IMU_a_m;
    };
    b_z.set_mean_state(z_vect);
    b_u.set_mean_state(cur_meas.u);
    
    if(sat_options.system_kind & ctrl::satellite_model_options::TSOSAKF) {
      ctrl::tsos_aug_inv_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    } else {
      ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    };
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() && ( meas_provider.get_current_time() < start_time ) );
  
  // prediction phase:
  do {
    b_u.set_mean_state(meas_provider.get_current_measurement().u);
    
    if(sat_options.system_kind & ctrl::satellite_model_options::TSOSAKF) {
      ctrl::tsos_aug_inv_kalman_predict(sat_sys, state_space, b, b_u, meas_provider.get_current_time());
    } else {
      ctrl::invariant_kalman_predict(sat_sys, state_space, b, b_u, meas_provider.get_current_time());
    };
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() );
  
  result_logger.finalize();
};




template <typename MeasureProvider, typename ResultLogger, typename Sat3DSystemType>
void batch_KF_ML_meas_predict(
    MeasureProvider meas_provider, 
    ResultLogger result_logger,
    const ReaK::ctrl::satellite_predictor_options& sat_options,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    double start_time,
    typename Sat3DSystemType::state_belief_type b,
    typename Sat3DSystemType::input_belief_type b_u,
    typename Sat3DSystemType::output_belief_type b_z) {
  using namespace ReaK;
  
  result_logger.initialize();
  
  // filtering phase:
  do {
    const sat3D_measurement_point& cur_meas = meas_provider.get_current_measurement();
    vect_n<double> z_vect(cur_meas.pose.size() + cur_meas.gyro.size() + cur_meas.IMU_a_m.size(), 0.0);
    z_vect[range(0,7)] = cur_meas.pose;
    if( cur_meas.gyro.size() ) {
      z_vect[range(7,10)] = cur_meas.gyro;
      if( cur_meas.IMU_a_m.size() )
        z_vect[range(10,16)] = cur_meas.IMU_a_m;
    };
    b_z.set_mean_state(z_vect);
    b_u.set_mean_state(cur_meas.u);
    
    ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() && ( meas_provider.get_current_time() < start_time ) );
  
  // prediction phase:
  do {
    b_u.set_mean_state(meas_provider.get_current_measurement().u);
    
    if(sat_options.system_kind & ctrl::satellite_model_options::TSOSAKF) {
      ctrl::tsos_aug_inv_kalman_predict(sat_sys, state_space, b, b_u, meas_provider.get_current_time());
    } else {
      ctrl::invariant_kalman_predict(sat_sys, state_space, b, b_u, meas_provider.get_current_time());
    };
    
    // apply ML assumption:
    b_z.set_mean_state(sat_sys.get_output(state_space, b.get_mean_state(), b_u.get_mean_state(), meas_provider.get_current_time()));
    
    if(sat_options.system_kind & ctrl::satellite_model_options::TSOSAKF) {
      ctrl::tsos_aug_inv_kalman_update(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    } else {
      ctrl::invariant_kalman_update(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    };
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() );
  
  result_logger.finalize();
};



template <typename Sat3DSystemType, typename MLTrajType>
typename boost::enable_if<
  ReaK::ctrl::is_augmented_ss_system<Sat3DSystemType>,
ReaK::shared_ptr< ReaK::pp::trajectory_base< ReaK::pp::temporal_space<ReaK::pp::se3_1st_order_topology<double>::type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> > > >::type
  construct_wrapped_trajectory(const ReaK::shared_ptr< MLTrajType >& ML_traj, double dt, double mtime) {
  using namespace ReaK;
  using namespace ctrl;
  using namespace pp;
  
  typedef se3_1st_order_topology<double>::type BaseSpaceType;
  typedef temporal_space<BaseSpaceType, time_poisson_topology, time_distance_only> TemporalBaseSpaceType;
  
#define RK_D_INF std::numeric_limits<double>::infinity()
  shared_ptr< TemporalBaseSpaceType > sat_base_temp_space(new TemporalBaseSpaceType(
    "satellite3D_temporal_space", 
    make_se3_space(
      "satellite3D_state_space",
      vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
      vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
      RK_D_INF, RK_D_INF),
    time_poisson_topology("satellite3D_time_space", dt, mtime)));
#undef RK_D_INF
  
  typedef transformed_trajectory<TemporalBaseSpaceType, MLTrajType, augmented_to_state_map> BaseTrajType;
  typedef trajectory_base<TemporalBaseSpaceType> StateTrajType;
  typedef trajectory_wrapper<BaseTrajType> WrappedStateTrajType;
  
  return shared_ptr< StateTrajType >(new WrappedStateTrajType("sat3D_predicted_traj", BaseTrajType(sat_base_temp_space, ML_traj)));
};

template <typename Sat3DSystemType, typename MLTrajType>
typename boost::enable_if<
  boost::mpl::not_< ReaK::ctrl::is_augmented_ss_system<Sat3DSystemType> >,
ReaK::shared_ptr< ReaK::pp::trajectory_base< ReaK::pp::temporal_space<ReaK::pp::se3_1st_order_topology<double>::type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> > > >::type
  construct_wrapped_trajectory(const ReaK::shared_ptr< MLTrajType >& ML_traj, double, double) {
  using namespace ReaK;
  using namespace pp;
  
  typedef ReaK::pp::trajectory_base< ReaK::pp::temporal_space<ReaK::pp::se3_1st_order_topology<double>::type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> > StateTrajType;
  typedef trajectory_wrapper<MLTrajType> WrappedStateTrajType;
  
  return shared_ptr< StateTrajType >(new WrappedStateTrajType("sat3D_predicted_traj", *ML_traj));
};


namespace {

template <typename Sat3DSystemType, typename MeasureProvider, typename ResultLogger>
struct prediction_updater {
  
  
  typedef typename Sat3DSystemType::state_space_type StateSpaceType;
  typedef typename Sat3DSystemType::temporal_state_space_type TempSpaceType;
  typedef typename Sat3DSystemType::belief_space_type BeliefSpaceType;
  typedef typename Sat3DSystemType::temporal_belief_space_type TempBeliefSpaceType;
  typedef typename Sat3DSystemType::covar_type CovarType;
  typedef typename CovarType::matrix_type CovarMatType;
  
  typedef typename Sat3DSystemType::point_type StateType;
  typedef typename Sat3DSystemType::state_belief_type StateBeliefType;
  typedef typename Sat3DSystemType::input_belief_type InputBeliefType;
  typedef typename Sat3DSystemType::output_belief_type OutputBeliefType;
  
  typedef typename ReaK::pp::topology_traits< TempBeliefSpaceType >::point_type TempBeliefPointType;
  
  typedef ReaK::pp::constant_trajectory< ReaK::pp::vector_topology< ReaK::vect_n<double> > > InputTrajType;
  
  typedef typename ReaK::ctrl::try_TSOSAIKF_belief_transfer_factory<Sat3DSystemType>::type PredFactoryType;
  typedef ReaK::ctrl::belief_predicted_trajectory<BeliefSpaceType, PredFactoryType, InputTrajType> BeliefPredTrajType;
  
  
  ReaK::shared_ptr< BeliefPredTrajType > predictor;
  
  ReaK::shared_ptr< Sat3DSystemType > satellite3D_system;
  ReaK::shared_ptr< TempSpaceType > sat_temp_space;
  
  MeasureProvider& meas_provider;
  ResultLogger& result_logger;
  
  StateBeliefType b;
  InputBeliefType b_u;
  OutputBeliefType b_z;
  
  double last_time;
  double diff_tolerance;
  
  ReaKaux::atomic<double>* current_target_anim_time;
  
  static ReaKaux::atomic<bool> should_stop;
  
  prediction_updater(
    ReaK::shared_ptr< BeliefPredTrajType > aPredictor,
    ReaK::shared_ptr< Sat3DSystemType > aSatSys,
    ReaK::shared_ptr< TempSpaceType > aSatTempSpace,
    MeasureProvider& aMeasProvider,
    ResultLogger& aResultLogger,
    StateBeliefType aB,
    InputBeliefType aBU,
    OutputBeliefType aBZ,
    double aLastTime, 
    double aDiffTolerance,
    ReaKaux::atomic<double>* aCurrentTargetAnimTime
  ) : 
    predictor(aPredictor),
    satellite3D_system(aSatSys),
    sat_temp_space(aSatTempSpace), 
    meas_provider(aMeasProvider),
    result_logger(aResultLogger),
    b(aB), b_u(aBU), b_z(aBZ),
    last_time(aLastTime), 
    diff_tolerance(aDiffTolerance),
    current_target_anim_time(aCurrentTargetAnimTime) { };
  
  int operator()() {
    
    typedef ReaK::pp::transformed_trajectory<TempSpaceType, BeliefPredTrajType, ReaK::ctrl::maximum_likelihood_map> MLTrajType;
    
    ReaK::shared_ptr< MLTrajType > ML_traj(new MLTrajType(sat_temp_space, predictor));
    
    typedef ReaK::pp::se3_1st_order_topology<double>::type BaseSpaceType;
    typedef ReaK::pp::temporal_space<BaseSpaceType, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> TemporalBaseSpaceType;
    typedef ReaK::pp::topology_traits<TemporalBaseSpaceType>::point_type BaseTempStateType;
    typedef ReaK::pp::topology_traits<BaseSpaceType>::point_type BaseStateType;
    typedef ReaK::ctrl::covariance_matrix< ReaK::vect_n<double> > CovarType;
    typedef ReaK::ctrl::gaussian_belief_state< BaseStateType, CovarType > BaseStateBeliefType;
    
    ReaK::shared_ptr< ReaK::pp::trajectory_base< TemporalBaseSpaceType > > pred_traj =
      construct_wrapped_trajectory<Sat3DSystemType>(ML_traj, 
                                                    sat_temp_space->get_time_topology().time_step, 
                                                    sat_temp_space->get_time_topology().mean_discrete_time);
    
    try {
      do {
        
        const sat3D_measurement_point& cur_meas = meas_provider.get_current_measurement();
        ReaK::vect_n<double> z_vect(cur_meas.pose.size() + cur_meas.gyro.size() + cur_meas.IMU_a_m.size(), 0.0);
        z_vect[ReaK::range(0,7)] = cur_meas.pose;
        if( cur_meas.gyro.size() ) {
          z_vect[ReaK::range(7,10)] = cur_meas.gyro;
          if( cur_meas.IMU_a_m.size() )
            z_vect[ReaK::range(10,16)] = cur_meas.IMU_a_m;
        };
        b_z.set_mean_state(z_vect);
        b_u.set_mean_state(cur_meas.u);
        
        ReaK::ctrl::tsos_aug_inv_kalman_filter_step(*satellite3D_system, 
                                                    sat_temp_space->get_space_topology(), 
                                                    b, b_u, b_z, last_time);
        
        last_time = meas_provider.get_current_time();
        (*current_target_anim_time) = last_time;
        
#if 1
        BaseTempStateType x_pred = pred_traj->get_point_at_time(last_time);
        result_logger.add_record(BaseStateBeliefType(x_pred.pt, CovarType(12)), b_u, b_z, meas_provider.get_current_time(), 
                                 meas_provider.get_current_gnd_truth_ptr());
#else
        TempBeliefPointType b_pred = predictor->get_point_at_time(last_time);
        result_logger.add_record(b_pred.pt, b_u, b_z, meas_provider.get_current_time(), 
                                 meas_provider.get_current_gnd_truth_ptr());
#endif
        
#if 0
        if( predictor->get_temporal_space().get_space_topology().distance(b, b_pred.pt) > diff_tolerance )
          predictor->set_initial_point(TempBeliefPointType(last_time, b));
#endif
      } while( !should_stop && meas_provider.step_once() );
    } catch (std::exception& e) {
      return 0;
    };
    return 0;
  };
  
  static ReaKaux::thread executer;
  
  static void stop_function() {
    should_stop = true;
    if(executer.joinable())
      executer.join();
  };
  
};


template <typename Sat3DSystemType, typename MeasureProvider, typename ResultLogger>
ReaKaux::atomic<bool> prediction_updater<Sat3DSystemType, MeasureProvider, ResultLogger>::should_stop = false;

template <typename Sat3DSystemType, typename MeasureProvider, typename ResultLogger>
ReaKaux::thread prediction_updater<Sat3DSystemType, MeasureProvider, ResultLogger>::executer = ReaKaux::thread();

boost::function<void()> pred_stop_function;

};





template <typename MeasureProvider, typename ResultLogger, typename Sat3DSystemType>
void batch_KF_meas_predict_with_predictor(
    MeasureProvider meas_provider, 
    ResultLogger result_logger,
    const ReaK::ctrl::satellite_predictor_options& sat_options,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    double start_time,
    typename Sat3DSystemType::state_belief_type b,
    typename Sat3DSystemType::input_belief_type b_u,
    typename Sat3DSystemType::output_belief_type b_z) {
  using namespace ReaK;
  using namespace ctrl;
  using namespace pp;
  
  result_logger.initialize();
  
  // filtering phase:
  do {
    const sat3D_measurement_point& cur_meas = meas_provider.get_current_measurement();
    vect_n<double> z_vect(cur_meas.pose.size() + cur_meas.gyro.size() + cur_meas.IMU_a_m.size(), 0.0);
    z_vect[range(0,7)] = cur_meas.pose;
    if( cur_meas.gyro.size() ) {
      z_vect[range(7,10)] = cur_meas.gyro;
      if( cur_meas.IMU_a_m.size() )
        z_vect[range(10,16)] = cur_meas.IMU_a_m;
    };
    b_z.set_mean_state(z_vect);
    b_u.set_mean_state(cur_meas.u);
    
    ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, meas_provider.get_current_time());
    
    result_logger.add_record(b, b_u, b_z, meas_provider.get_current_time(), 
                             meas_provider.get_current_gnd_truth_ptr());
    
  } while( meas_provider.step_once() && ( meas_provider.get_current_time() < start_time ) );
  

// Start a thread that will update the predictor as new data comes in.
  
  if( pred_stop_function ) {
    pred_stop_function();
    pred_stop_function.clear();
  };
  
  typedef constant_trajectory< vector_topology< vect_n<double> > > InputTrajType;
  typedef typename Sat3DSystemType::temporal_state_space_type TempSpaceType;
  typedef typename Sat3DSystemType::belief_space_type BeliefSpaceType;
  typedef typename Sat3DSystemType::temporal_belief_space_type TempBeliefSpaceType;
  typedef typename topology_traits< TempBeliefSpaceType >::point_type TempBeliefPointType;
  typedef typename Sat3DSystemType::covar_type CovarType;
  typedef typename CovarType::matrix_type CovarMatType;
  
  typedef typename try_TSOSAIKF_belief_transfer_factory<Sat3DSystemType>::type PredFactoryType;
  typedef belief_predicted_trajectory<BeliefSpaceType, PredFactoryType, InputTrajType> BeliefPredTrajType;
  
  typename BeliefPredTrajType::assumption pred_assumpt = BeliefPredTrajType::no_measurements;
  switch(sat_options.predict_assumption) {
    case satellite_predictor_options::no_measurements:  // No future measurements
      pred_assumpt = BeliefPredTrajType::no_measurements;
      break;
    case satellite_predictor_options::most_likely_measurements:  // Maximum-likelihood Measurements
      pred_assumpt = BeliefPredTrajType::most_likely_measurements;
      break;
    case satellite_predictor_options::full_certainty:  // Full certainty
      // FIXME: make this into what it really should be (what is that? ... no sure)
      pred_assumpt = BeliefPredTrajType::most_likely_measurements;
      break;
  };
  
  shared_ptr< Sat3DSystemType > sat_sys_ptr(&sat_sys, null_deleter());
  shared_ptr< TempSpaceType > sat_temp_space = sat_sys.get_temporal_state_space(0.0, sat_options.predict_time_horizon);
  shared_ptr< TempBeliefSpaceType > sat_temp_belief_space = sat_sys.get_temporal_belief_space(0.0, sat_options.predict_time_horizon);
  
  shared_ptr< BeliefPredTrajType > predictor(new BeliefPredTrajType(
    sat_temp_belief_space, 
    TempBeliefPointType(meas_provider.get_current_time(), b), 
    InputTrajType( vect_n<double>(0.0,0.0,0.0,0.0,0.0,0.0) ),
    PredFactoryType(sat_sys_ptr, 
                    CovarMatType(sat_options.input_disturbance), 
                    CovarMatType(sat_options.measurement_noise)),
    pred_assumpt
  ));
  
  
  ReaKaux::atomic<double> current_target_anim_time = meas_provider.get_current_time();
  
  prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::should_stop = false;
  prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::executer = ReaKaux::thread(
    prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>(
      predictor, sat_sys_ptr, sat_temp_space, meas_provider, result_logger,
      b, b_u, b_z, meas_provider.get_current_time(), 0.5, &current_target_anim_time));
  
  pred_stop_function = prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::stop_function;
  
  if(prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::executer->joinable())
    prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::executer->join();
  prediction_updater<Sat3DSystemType,MeasureProvider,ResultLogger>::executer.reset();
  pred_stop_function.clear();
  
  result_logger.finalize();
  
};








template <typename Sat3DSystemType>
void generate_timeseries(
    std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    typename Sat3DSystemType::point_type x,
    double start_time, double end_time,
    const cov_matrix_type& Qu,
    const cov_matrix_type& R,
    ReaK::shared_ptr< ReaK::recorder::data_recorder > stat_results = ReaK::shared_ptr< ReaK::recorder::data_recorder >() ) {
  using namespace ReaK;
  
  measurements.clear();
  ground_truth.clear();
  
  double time_step = sat_sys.get_time_step();
  vect_n< double > std_devs(R.get_row_count() + R.get_row_count() / 3, 0.0);
  for(double t = start_time; t < end_time; t += time_step) {
    vect_n<double> u = ctrl::sample_gaussian_point(vect_n<double>(6, 0.0), Qu);
    
    x = sat_sys.get_next_state(state_space, x, u, t);
    ground_truth.push_back(std::make_pair(t, get_sat3D_state(x)));
    
    vect_n<double> y       = sat_sys.get_output(state_space, x, u, t);
    vect_n<double> y_noise = ctrl::sample_gaussian_point(sat_sys.get_invariant_error(state_space, x, u, y, t), R);
    
    sat3D_measurement_point meas;
    meas.u = vect_n<double>(6, 0.0);
    meas.pose = vect_n<double>(7, 0.0);
    meas.pose[range(0,3)] = y[range(0,3)] + y_noise[range(0,3)];
    
    vect<double,3> aa_noise(y_noise[3],y_noise[4],y_noise[5]);
    quaternion<double> y_quat(y[range(3,7)]);
    y_quat *= axis_angle<double>(norm_2(aa_noise), aa_noise).getQuaternion();
    meas.pose[range(3,7)] = vect<double,4>(y_quat[0], y_quat[1], y_quat[2], y_quat[3]);
    
    std::size_t k = ground_truth.size();
    if( stat_results ) {
      std_devs[range(0,6)] = ( ( k - 1 ) * std_devs[range(0,6)] + elem_product(y_noise[range(0,6)], y_noise[range(0,6)]) ) / k;
      std_devs[6] = ( ( k - 1 ) * std_devs[6] + norm_2_sqr(y_noise[range(0,3)]) ) / k;
      std_devs[7] = ( ( k - 1 ) * std_devs[7] + norm_2_sqr(aa_noise) ) / k;
    };
    
    if( y.size() >= 10 ) {
      meas.gyro = y[range(7,10)] + y_noise[range(6,9)];
      if( stat_results ) {
        std_devs[range(8,11)]  = ( ( k - 1 ) * std_devs[range(8,11)] + elem_product(y_noise[range(6,9)], y_noise[range(6,9)]) ) / k;
        std_devs[11] = ( ( k - 1 ) * std_devs[11] + norm_2_sqr(y_noise[range(6,9)]) ) / k;
      };
      if( y.size() >= 16 ) {
        meas.IMU_a_m = y[range(10,16)] + y_noise[range(9,15)];
        if( stat_results ) {
          std_devs[range(12,15)]  = ( ( k - 1 ) * std_devs[range(12,15)] + elem_product(y_noise[range(9,12)], y_noise[range(9,12)]) ) / k;
          std_devs[15] = ( ( k - 1 ) * std_devs[15] + norm_2_sqr(y_noise[range(9,12)]) ) / k;
          std_devs[range(16,19)]  = ( ( k - 1 ) * std_devs[range(16,19)] + elem_product(y_noise[range(12,15)], y_noise[range(12,15)]) ) / k;
          std_devs[19] = ( ( k - 1 ) * std_devs[19] + norm_2_sqr(y_noise[range(12,15)]) ) / k;
        };
      };
    };
    measurements.push_back(std::make_pair(t, meas));
  };
  
  if( stat_results ) {
    for(std::size_t i = 0; i < std_devs.size(); ++i)
      (*stat_results) << std::sqrt(std_devs[i]);
    (*stat_results) << recorder::data_recorder::end_value_row << recorder::data_recorder::flush;
  };
};






template <typename Sat3DSystemType>
void do_online_run(
    ReaK::recorder::data_stream_options output_opt,
    const ReaK::ctrl::satellite_model_options& sat_options,
    const ReaK::shared_ptr< ReaK::recorder::data_extractor >& data_in,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    const typename Sat3DSystemType::state_belief_type& b,
    typename Sat3DSystemType::input_belief_type b_u,
    const typename Sat3DSystemType::output_belief_type& b_z) {
  
  std::stringstream ss;
  ss << "_" << std::setfill('0') << std::setw(4) << int(1000 * sat_options.time_step) << "_" << sat_options.get_kf_accronym() << ".";
  ReaK::recorder::data_stream_options cur_out_opt = output_opt;
  cur_out_opt.file_name += ss.str() + cur_out_opt.get_extension();
  sat_options.imbue_names_for_state_estimates(cur_out_opt);
  
  batch_KF_on_timeseries(
    sat3D_meas_true_from_extractor(data_in, sat_options),
    sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
    sat_options, sat_sys, state_space, b, b_u, b_z);
  
};

template <typename Sat3DSystemType>
void do_all_single_runs(
    ReaK::recorder::data_stream_options output_opt,
    const ReaK::ctrl::satellite_model_options& sat_options,
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    const typename Sat3DSystemType::state_belief_type& b,
    typename Sat3DSystemType::input_belief_type b_u,
    const typename Sat3DSystemType::output_belief_type& b_z,
    unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;
  
  cov_matrix_type Qu = b_u.get_covariance().get_matrix();
  
  for(unsigned int skips = min_skips; skips <= max_skips; ++skips) {
    
    sat_sys.set_time_step(skips * sat_options.time_step);
    
    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));
    
    std::stringstream ss;
    ss << "_" << std::setfill('0') << std::setw(4) << int(1000 * skips * sat_options.time_step) << "_" << sat_options.get_kf_accronym() << ".";
    recorder::data_stream_options cur_out_opt = output_opt;
    cur_out_opt.file_name += ss.str() + cur_out_opt.get_extension();
    sat_options.imbue_names_for_state_estimates(cur_out_opt);
    
    batch_KF_on_timeseries(
      sat3D_meas_true_from_vectors(&measurements, &ground_truth, skips),
      sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
      sat_options, sat_sys, state_space, b, b_u, b_z);
    
  };
  
  sat_sys.set_time_step(sat_options.time_step);
  
};



template <typename Sat3DSystemType>
void do_online_prediction(
    ReaK::recorder::data_stream_options output_opt,
    const ReaK::ctrl::satellite_predictor_options& sat_options,
    const ReaK::shared_ptr< ReaK::recorder::data_extractor >& data_in,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    const typename Sat3DSystemType::state_belief_type& b,
    typename Sat3DSystemType::input_belief_type b_u,
    const typename Sat3DSystemType::output_belief_type& b_z,
    double start_time) {
  using namespace ReaK;
  
  std::stringstream ss;
  ss << "_pred_" << std::setfill('0') << std::setw(5) << int(100 * start_time) << "_" << sat_options.get_kf_accronym() << ".";
  ReaK::recorder::data_stream_options cur_out_opt = output_opt;
  cur_out_opt.file_name += ss.str() + cur_out_opt.get_extension();
  sat_options.imbue_names_for_state_estimates(cur_out_opt);
  
#if 0
  batch_KF_meas_predict_with_predictor(
    sat3D_meas_true_from_extractor(data_in, sat_options),
    sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
    sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
#else
  if(sat_options.predict_assumption == ReaK::ctrl::satellite_predictor_options::no_measurements) {
    batch_KF_no_meas_predict(
      sat3D_meas_true_from_extractor(data_in, sat_options),
      sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
      sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
  } else {
    batch_KF_ML_meas_predict(
      sat3D_meas_true_from_extractor(data_in, sat_options),
      sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
      sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
  };
#endif
};


template <typename Sat3DSystemType>
void do_all_prediction_runs(
    ReaK::recorder::data_stream_options output_opt,
    const ReaK::ctrl::satellite_predictor_options& sat_options,
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    const typename Sat3DSystemType::state_belief_type& b,
    typename Sat3DSystemType::input_belief_type b_u,
    const typename Sat3DSystemType::output_belief_type& b_z,
    double start_intervals) {
  using namespace ReaK;
  
  if(measurements.size() == 0)
    return;
  
  double end_time = (measurements.end()-1)->first;
  
  for(double start_time = measurements.begin()->first + start_intervals; start_time < end_time; start_time += start_intervals) {
    
    std::stringstream ss;
    ss << "_pred_" << std::setfill('0') << std::setw(5) << int(100 * (start_time - measurements.begin()->first)) << "_" << sat_options.get_kf_accronym() << ".";
    recorder::data_stream_options cur_out_opt = output_opt;
    cur_out_opt.file_name += ss.str() + cur_out_opt.get_extension();
    sat_options.imbue_names_for_state_estimates(cur_out_opt);
    
#if 0
    batch_KF_meas_predict_with_predictor(
      sat3D_meas_true_from_vectors(&measurements, &ground_truth),
      sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
      sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
#else
    if(sat_options.predict_assumption == ctrl::satellite_predictor_options::no_measurements) {
      batch_KF_no_meas_predict(
        sat3D_meas_true_from_vectors(&measurements, &ground_truth),
        sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
        sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
    } else {
      batch_KF_ML_meas_predict(
        sat3D_meas_true_from_vectors(&measurements, &ground_truth),
        sat3D_estimate_result_to_recorder(cur_out_opt.create_recorder()), 
        sat_options, sat_sys, state_space, start_time, b, b_u, b_z);
    };
#endif
  };
  
};




template <typename Sat3DSystemType>
void do_single_monte_carlo_run(
    std::map< std::string, ReaK::shared_ptr< ReaK::recorder::data_recorder > >& results_map,
    ReaK::recorder::data_stream_options output_opt,
    const ReaK::ctrl::satellite_model_options& sat_options,
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const typename Sat3DSystemType::state_space_type& state_space,
    const typename Sat3DSystemType::state_belief_type& b,
    typename Sat3DSystemType::input_belief_type b_u,
    const typename Sat3DSystemType::output_belief_type& b_z,
    unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;
  
  cov_matrix_type Qu = b_u.get_covariance().get_matrix();
  
  for(unsigned int skips = min_skips; skips <= max_skips; ++skips) {
    
    sat_sys.set_time_step(skips * sat_options.time_step);
    
    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));
    
    std::stringstream ss;
    ss << "_" << std::setfill('0') << std::setw(4) << int(1000 * skips * sat_options.time_step) << "_" << sat_options.get_kf_accronym();
    std::string file_middle = ss.str();
    shared_ptr< recorder::data_recorder >& results = results_map[file_middle];
    if( !results ) {
      recorder::data_stream_options cur_out_opt = output_opt;
      cur_out_opt.file_name += file_middle + "_stddevs." + cur_out_opt.get_extension();
      sat_options.imbue_names_for_state_estimates_stddevs(cur_out_opt);
      results = cur_out_opt.create_recorder();
    };
    
    batch_KF_on_timeseries(
      sat3D_meas_true_from_vectors(&measurements, &ground_truth, skips),
      sat3D_collect_stddevs(results), 
      sat_options, sat_sys, state_space, b, b_u, b_z);
    
  };
  
  sat_sys.set_time_step(sat_options.time_step);
  
};


static
void get_timeseries_from_rec(
    const ReaK::shared_ptr< ReaK::recorder::data_extractor >& data_in,
    const std::vector<std::string>& names_in,
    const ReaK::ctrl::satellite_model_options& sat_options,
    std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    std::vector< std::pair< double, sat3D_state_type > >& ground_truth) {
  using namespace ReaK;
  
  measurements.clear();
  ground_truth.clear();
  
  std::cout << "Reading data file..." << std::endl;
  
  try {
    recorder::named_value_row nvr_in = data_in->getFreshNamedValueRow();
    while(true) {
      (*data_in) >> nvr_in;
      
      
      double t = nvr_in["time"];
      
      sat3D_measurement_point meas_actual;
      sat3D_measurement_point meas_noisy;
      
      meas_actual.pose.resize(7);
      meas_actual.pose[0] = nvr_in["p_x"];
      meas_actual.pose[1] = nvr_in["p_y"];
      meas_actual.pose[2] = nvr_in["p_z"];
      meas_actual.pose[3] = nvr_in["q_0"];
      meas_actual.pose[4] = nvr_in["q_1"];
      meas_actual.pose[5] = nvr_in["q_2"];
      meas_actual.pose[6] = nvr_in["q_3"];
      
      std::size_t merr_count = sat_options.get_meas_error_count();
      vect_n<double> added_noise = ctrl::sample_gaussian_point(vect_n<double>(sat_options.artificial_noise.get_row_count(), 0.0), sat_options.artificial_noise);
      if(sat_options.artificial_noise.get_row_count() < merr_count) 
        added_noise.resize(merr_count, 0.0);
      
      meas_noisy.pose.resize(7);
      meas_noisy.pose[range(0,3)] = meas_actual.pose[range(0,3)] + added_noise[range(0,3)];
      
      vect<double,3> aa_noise(added_noise[3], added_noise[4], added_noise[5]);
      quaternion<double> y_quat(meas_actual.pose[range(3,7)]);
      y_quat *= axis_angle<double>(norm_2(aa_noise), aa_noise).getQuaternion();
      meas_noisy.pose[range(3,7)] = vect<double,4>(y_quat[0], y_quat[1], y_quat[2], y_quat[3]);
      
      if( merr_count >= 9 ) {
        /* read off the IMU/gyro angular velocity measurements. */
        meas_actual.gyro.resize(3);
        meas_actual.gyro[0] = nvr_in["w_x"];
        meas_actual.gyro[1] = nvr_in["w_y"];
        meas_actual.gyro[2] = nvr_in["w_z"];
        meas_noisy.gyro = meas_actual.gyro + added_noise[range(6,9)];
        if( merr_count >= 15 ) {
          /* read off the IMU accel-mag measurements. */
          meas_actual.IMU_a_m.resize(6);
          meas_actual.IMU_a_m[0] = nvr_in["acc_x"];
          meas_actual.IMU_a_m[1] = nvr_in["acc_y"];
          meas_actual.IMU_a_m[2] = nvr_in["acc_z"];
          meas_actual.IMU_a_m[3] = nvr_in["mag_x"];
          meas_actual.IMU_a_m[4] = nvr_in["mag_y"];
          meas_actual.IMU_a_m[5] = nvr_in["mag_z"];
          meas_noisy.IMU_a_m = meas_actual.IMU_a_m + added_noise[range(9,15)];
        };
      };
      
      /* read off the input vector. */
      meas_actual.u.resize(6);
      meas_actual.u[0] = nvr_in["f_x"];
      meas_actual.u[1] = nvr_in["f_y"];
      meas_actual.u[2] = nvr_in["f_z"];
      meas_actual.u[3] = nvr_in["t_x"];
      meas_actual.u[4] = nvr_in["t_y"];
      meas_actual.u[5] = nvr_in["t_z"];
      meas_noisy.u  = meas_actual.u;
      
      /* now, the meas_actual and meas_noisy are fully formed. */
      measurements.push_back( std::make_pair(t, meas_noisy) );
      
      std::cout << "\r" << std::setw(10) << measurements.size() << std::flush;
      
      /* check if the file contains a ground-truth: */
      try {
        sat3D_state_type x;
        set_position(x, vect<double,3>(nvr_in["p_x_true"], nvr_in["p_y_true"], nvr_in["p_z_true"]));
        set_quaternion(x, unit_quat<double>(nvr_in["q_0_true"], nvr_in["q_1_true"], nvr_in["q_2_true"], nvr_in["q_3_true"]));
        set_velocity(x, vect<double,3>(nvr_in["v_x_true"], nvr_in["v_y_true"], nvr_in["v_z_true"]));
        set_ang_velocity(x, vect<double,3>(nvr_in["w_x_true"], nvr_in["w_y_true"], nvr_in["w_z_true"]));
        ground_truth.push_back( std::make_pair(t, x) );
      } catch(ReaK::recorder::out_of_bounds& e) { RK_UNUSED(e); 
        if( sat_options.artificial_noise.get_row_count() >= 6 ) {
          sat3D_state_type x;
          set_position(x, vect<double,3>(meas_actual.pose[0],meas_actual.pose[1],meas_actual.pose[2]));
          set_quaternion(x, unit_quat<double>(meas_actual.pose[3],meas_actual.pose[4],meas_actual.pose[5],meas_actual.pose[6]));
          set_velocity(x, vect<double,3>(0.0, 0.0, 0.0));
          set_ang_velocity(x, vect<double,3>(0.0, 0.0, 0.0));
          ground_truth.push_back( std::make_pair(t, x) );
        };
      };
      
    };
  } catch(recorder::end_of_record&) { };
  
  std::cout << "\nDone!" << std::endl;
  
};



template <typename Sat3DSystemType>
int do_required_tasks(ReaK::shared_ptr< Sat3DSystemType > satellite3D_system, 
                      const ReaK::ctrl::satellite_predictor_options& sat_options, 
                      boost::program_options::variables_map& vm,
                      ReaK::shared_ptr< ReaK::recorder::data_extractor > data_in,
                      const std::vector<std::string>& names_in,
                      const std::string& sys_output_stem_name,
                      const ReaK::recorder::data_stream_options& data_out_stem_opt) {
  using namespace ReaK;
  
  
  double start_time = vm["start-time"].as<double>();
  double end_time   = vm["end-time"].as<double>();
  
  unsigned int mc_runs    = vm["mc-runs"].as<unsigned int>();
  unsigned int min_skips  = vm["min-skips"].as<unsigned int>();
  unsigned int max_skips  = vm["max-skips"].as<unsigned int>();
  
  typedef typename Sat3DSystemType::temporal_state_space_type TempSpaceType;
  typedef typename Sat3DSystemType::covar_type CovarType;
  typedef typename CovarType::matrix_type CovarMatType;
  typedef typename Sat3DSystemType::point_type StateType;
  typedef typename Sat3DSystemType::state_belief_type StateBeliefType;
  typedef typename Sat3DSystemType::input_belief_type InputBeliefType;
  typedef typename Sat3DSystemType::output_belief_type OutputBeliefType;
  
  shared_ptr< TempSpaceType > sat_space = satellite3D_system->get_temporal_state_space(start_time, end_time);
  
  StateBeliefType b_init = satellite3D_system->get_zero_state_belief(10.0);
  if( ( b_init.get_covariance().get_matrix().get_row_count() > 12 ) &&
      ( sat_options.system_kind & ctrl::satellite_predictor_options::TSOSAKF ) &&
      ( sat_options.steady_param_covariance.get_row_count() + 12 == b_init.get_covariance().get_matrix().get_row_count() ) ) {
    mat<double,mat_structure::square> P(b_init.get_covariance().get_matrix());
    set_block(P, sat_options.steady_param_covariance, 12, 12);
    b_init.set_covariance(CovarType(CovarMatType(P)));
  };
  
  InputBeliefType b_u = satellite3D_system->get_zero_input_belief();
  b_u.set_covariance(CovarType(CovarMatType(sat_options.input_disturbance)));
  
  OutputBeliefType b_z = satellite3D_system->get_zero_output_belief();
  b_z.set_covariance(CovarType(CovarMatType(sat_options.measurement_noise)));
  
  std::vector< std::pair< double, sat3D_measurement_point > > measurements;
  std::vector< std::pair< double, sat3D_state_type > > ground_truth;
  
  typedef pp::discrete_point_trajectory< sat3D_temp_space_type > RecTrajType;
  shared_ptr< RecTrajType > traj_ptr;
  if( vm.count("generate-meas-file") && ( vm.count("xml") + vm.count("protobuf") + vm.count("binary") > 0 ) )
    traj_ptr = shared_ptr< RecTrajType >(new RecTrajType( 
      shared_ptr<sat3D_temp_space_type>(new sat3D_temp_space_type(
        "satellite3D_temporal_space",
         sat3D_state_space_type(),
         pp::time_poisson_topology("satellite3D_time_space", satellite3D_system->get_time_step(), (end_time - start_time) * 0.5)))));
  
  
  if( vm.count("generate-mdl-files") ) {
    try {
      *(serialization::open_oarchive(sys_output_stem_name + sat_options.get_sys_abbreviation() + "_mdl.rkx"))
        & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
    } catch(...) {
      RK_ERROR("An exception occurred during the saving the satellite system file!");
      return 14;
    };
  } else if( vm.count("online-run") ) {
    
    if( !data_in ) {
      RK_ERROR("Must have a defined input data-stream in order to run the estimator online!");
      return 15;
    };
    
    if(!vm.count("prediction-runs")) {
      
      do_online_run(data_out_stem_opt, sat_options, data_in, *satellite3D_system, 
                    sat_space->get_space_topology(), b_init, b_u, b_z);
      
    } else if( !vm.count("monte-carlo") ) {
      
      do_online_prediction(data_out_stem_opt, sat_options, data_in, *satellite3D_system, 
                            sat_space->get_space_topology(), b_init, b_u, b_z, vm["prediction-interval"].as<double>());
      
    };
    
  } else if( ! vm.count("monte-carlo") ) {
    
    if( data_in ) {
      get_timeseries_from_rec(data_in, names_in, sat_options, measurements, ground_truth);
    } else {
      // must generate the measurements and ground_truth vectors:
      StateType x_init;
      sat3D_state_type x_st;
      set_frame_3D(x_st, sat_options.initial_motion);
      set_sat3D_state(x_init, x_st);
      mat<double,mat_structure::diagonal> R_tot = sat_options.measurement_noise;
      if(sat_options.artificial_noise.get_row_count() == R_tot.get_row_count())
        R_tot += sat_options.artificial_noise;
      generate_timeseries(measurements, ground_truth, *satellite3D_system, sat_space->get_space_topology(),
                          x_init, start_time, end_time, CovarMatType(sat_options.input_disturbance), 
                          CovarMatType(R_tot));
    };
    
    // do a single run for each skips:
    
    if(!vm.count("generate-meas-file")) {
      
      if(!vm.count("prediction-runs")) {
        std::cout << "Running estimator on data series.." << std::flush;
        
        do_all_single_runs(data_out_stem_opt, sat_options, measurements, ground_truth, *satellite3D_system, 
                          sat_space->get_space_topology(), b_init, b_u, b_z, min_skips, max_skips);
        
        std::cout << "." << std::flush;
      } else {
        std::cout << "Running predictor on data series.." << std::flush;
        
        do_all_prediction_runs(data_out_stem_opt, sat_options, measurements, ground_truth, *satellite3D_system, 
                                sat_space->get_space_topology(), b_init, b_u, b_z, vm["prediction-interval"].as<double>());
        
        std::cout << "." << std::flush;
      };
      
    } else {
      
      std::cout << "Recording generated data series.." << std::flush;
      
      recorder::data_stream_options data_meas_opt = data_out_stem_opt;
      data_meas_opt.file_name += "_meas." + data_meas_opt.get_extension();
      sat_options.imbue_names_for_generated_meas(data_meas_opt);
      shared_ptr< recorder::data_recorder > data_meas = data_meas_opt.create_recorder();
      for(std::size_t i = 0; i < measurements.size(); ++i) {
        (*data_meas) << measurements[i].first;
        const sat3D_measurement_point& m = measurements[i].second;
        (*data_meas) << m.pose << m.gyro << m.IMU_a_m << m.u;  // if gyro-IMU not present, vectors will be zero-sized, not written to stream.
        if(ground_truth.size() == measurements.size()) {
          const sat3D_state_type& g = ground_truth[i].second;
          (*data_meas) << get_position(g) << get_quaternion(g) << get_velocity(g) << get_ang_velocity(g);
          (*data_meas) << recorder::data_recorder::end_value_row;
          if( traj_ptr )
            traj_ptr->push_back(sat3D_temp_point_type(ground_truth[i].first, g));
        };
      };
      (*data_meas) << recorder::data_recorder::flush;
      
      if( traj_ptr ) {
        std::cout << "\nSaving the generated trajectory.." << std::flush;
        std::cout << "." << std::flush;
        
        if( vm.count("xml") )
          *(serialization::open_oarchive(data_out_stem_opt.file_name + "_traj.rkx"))
            & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", *traj_ptr);
        
        std::cout << "." << std::flush;
        
        if( vm.count("protobuf") )
          *(serialization::open_oarchive(data_out_stem_opt.file_name + "_traj.pbuf"))
            & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", *traj_ptr);
        
        std::cout << "." << std::flush;
        
        if( vm.count("binary") )
          *(serialization::open_oarchive(data_out_stem_opt.file_name + "_traj.rkb"))
            & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", *traj_ptr);
        
        std::cout << "Finished!" << std::endl;
      };
      
    };
    
    std::cout << "Finished!" << std::endl;
    
  } else {
    // do monte-carlo runs:
    StateType x_init;
    sat3D_state_type x_st;
    set_frame_3D(x_st, sat_options.initial_motion);
    set_sat3D_state(x_init, x_st);
    mat<double,mat_structure::diagonal> R_tot = sat_options.measurement_noise;
    if(sat_options.artificial_noise.get_row_count() == R_tot.get_row_count())
      R_tot += sat_options.artificial_noise;
    
    recorder::data_stream_options data_stddev_opt = data_out_stem_opt;
    data_stddev_opt.file_name += "_meas_stddevs." + data_stddev_opt.get_extension();
    sat_options.imbue_names_for_meas_stddevs(data_stddev_opt);
    shared_ptr< recorder::data_recorder > data_stddev = data_stddev_opt.create_recorder();
    
    std::map< std::string, shared_ptr< recorder::data_recorder > > results_map;
    
    std::cout << "Running Monte-Carlo Simulations..." << std::endl;
    
    for(unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {
      
      std::cout << "\r" << std::setw(10) << mc_i << std::flush;
      
      generate_timeseries(measurements, ground_truth, *satellite3D_system, sat_space->get_space_topology(),
                          x_init, start_time, end_time, CovarMatType(sat_options.input_disturbance), 
                          CovarMatType(R_tot), data_stddev);
      
      std::cout << "." << std::flush;
      
      do_single_monte_carlo_run(
        results_map, data_out_stem_opt, sat_options, measurements, ground_truth, *satellite3D_system, 
        sat_space->get_space_topology(), b_init, b_u, b_z, min_skips, max_skips);
      
      std::cout << "." << std::flush;
      
    };
    
    std::cout << "Finished!" << std::endl;
    
  };
  
  
  return 0;
};




int main(int argc, char** argv) {
  using namespace ReaK;
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("generate-meas", "if set, the measurements used for the estimation will be generated from a simulation with the given initial conditions (default is not)")
    ("generate-meas-file,g", "if set, the measurement file will be generated from the output of a simulation with the given initial conditions (default is not)")
  ;
  
  po::options_description sim_options("Simulation options");
  sim_options.add_options()
    ("start-time,s", po::value< double >()->default_value(0.0), "start time of the estimation (default is 0.0)")
    ("end-time,e",   po::value< double >()->default_value(1.0), "end time of the estimation (default is 1.0)")
    ("monte-carlo",  "if set, will perform a Monte-Carlo set of randomized runs to gather estimation performance statistics")
    ("mc-runs",      po::value< unsigned int >()->default_value(1000), "number of Monte-Carlo runs to perform (default is 1000)")
    ("min-skips",    po::value< unsigned int >()->default_value(1), "minimum number of time-step skips between estimations when generating a series of Monte-Carlo statistics (default is 1, i.e., one estimation point per measurement point)")
    ("max-skips",    po::value< unsigned int >()->default_value(1), "maximum number of time-step skips between estimations when generating a series of Monte-Carlo statistics (default is 1, i.e., one estimation point per measurement point)")
    ("prediction-runs", "if set, will perform prediction runs instead of estimation runs")
    ("prediction-interval", po::value< double >()->default_value(1.0), "time interval between new trials of the predictor on the data")
    ("online-run", "if set, will perform estimation or prediction online (no pre-buffering of data stream)")
  ;
  
  po::options_description output_options("Output options (at least one must be set)");
  output_options.add_options()
    ("output-traj-file", "if set, output results in a trajectory file (not data-stream)")
    ("xml,x",      "if set, output results in XML format (rkx)")
    ("protobuf,p", "if set, output results in protobuf format (pbuf)")
    ("binary,b",   "if set, output results in binary format (rkb)")
  ;
  
  po::options_description model_options = ctrl::get_satellite_model_options_po_desc(true);
  
  po::options_description data_io_options = recorder::get_data_stream_options_po_desc(true, true);
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(sim_options).add(model_options).add(output_options).add(data_io_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  recorder::data_stream_options data_in_opt;
  shared_ptr< recorder::data_extractor > data_in;
  std::vector<std::string> names_in;
  if(!vm.count("generate-meas") && !vm.count("generate-mdl-files")) {
    try {
      data_in_opt  = recorder::get_data_stream_options_from_po(vm, false);
      boost::tie(data_in, names_in) = data_in_opt.create_extractor();
    } catch(std::invalid_argument& e) {
      std::cerr << "Error! Creation of input data-stream failed! Invalid argument: " << e.what() << std::endl;
      return 2;
    };
  };
  
  
  recorder::data_stream_options data_out_opt;
  std::string output_stem_name;
  try {
    data_out_opt = recorder::get_data_stream_options_from_po(vm, true);
    
    output_stem_name = data_out_opt.file_name;
    if(output_stem_name[output_stem_name.size()-1] == '/')
      output_stem_name += "output_record";
    else {
      std::size_t last_dot   = output_stem_name.find_last_of('.');
      last_dot = ( last_dot == std::string::npos ? 0 : last_dot );
      std::size_t last_slash = output_stem_name.find_last_of('/');
      last_slash = ( last_slash == std::string::npos ? 0 : last_slash );
      if( last_dot > last_slash ) 
        output_stem_name.erase(output_stem_name.begin() + last_dot, output_stem_name.end());
    };
    
  } catch(std::invalid_argument& e) {
    std::cerr << "Error! Creation of output data-stream failed! Invalid argument: " << e.what() << std::endl;
    return 1;
  };
  recorder::data_stream_options data_out_stem_opt = data_out_opt;
  data_out_stem_opt.file_name = output_stem_name;
  
  ctrl::satellite_predictor_options sat_options;
  try {
    sat_options = ctrl::get_satellite_predictor_options_from_po(vm);
  } catch(std::exception& e) {
    std::cerr << "Error! Creation of satellite modeling options failed! With exception: " << e.what() << std::endl;
    return 2;
  };
  
  std::string sys_output_stem_name;
  if( vm.count("system-output") )
    sys_output_stem_name = vm["system-output"].as<std::string>();
  if( !sys_output_stem_name.empty() ) {
    std::string sys_output_path_name = sys_output_stem_name;
    if( vm.count("generate-mdl-files") ) {
      if(sys_output_stem_name[sys_output_stem_name.size()-1] == '/')
        sys_output_stem_name += "satellite3D";
      else {
        std::size_t p = sys_output_path_name.find_last_of('/');
        if(p == std::string::npos)
          sys_output_path_name = "";
        else
          sys_output_path_name.erase(p);
      };
      while(sys_output_path_name[sys_output_path_name.length()-1] == '/') 
        sys_output_path_name.erase(sys_output_path_name.length()-1, 1);
      
      if(!sys_output_path_name.empty())
        fs::create_directory(sys_output_path_name.c_str());
    };
  };
  
  if( !vm.count("gyro") && !vm.count("IMU") ) {
    
    if( vm.count("imkf-em") ) {
      int errcode = do_required_tasks(sat_options.get_em_airship_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    } else if( vm.count("imkf-emd") ) {
      int errcode = do_required_tasks(sat_options.get_emd_airship_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    } else if( vm.count("imkf-emdJ") ) {
      int errcode = do_required_tasks(sat_options.get_emdJ_airship_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    } else {
      int errcode = do_required_tasks(sat_options.get_base_sat_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    };
    
  } else if( vm.count("gyro") && !vm.count("IMU") ) {
    
    if( vm.count("imkf-emd") ) {
      int errcode = do_required_tasks(sat_options.get_gyro_emd_airship_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    } else if( vm.count("imkf-emdJ") ) {
      int errcode = do_required_tasks(sat_options.get_gyro_emdJ_airship_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    } else {
      int errcode = do_required_tasks(sat_options.get_gyro_sat_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
      if(errcode)
        return errcode;
    };
    
  } else {
    
//     int errcode = do_required_tasks(sat_options.get_IMU_sat_system(), sat_options, vm, data_in, names_in, sys_output_stem_name, data_out_stem_opt);
//     if(errcode)
//       return errcode;
    
  };
  
  
  return 0;
};



