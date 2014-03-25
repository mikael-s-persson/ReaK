
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


#include "ss_systems/satellite_invar_models.hpp"
#include "ss_systems/satellite_modeling_po.hpp"

#include "ctrl_sys/kalman_filter.hpp"
#include "ctrl_sys/invariant_kalman_filter.hpp"

#include "ctrl_sys/gaussian_belief_state.hpp"
#include "ctrl_sys/covariance_matrix.hpp"

#include "serialization/archiver_factory.hpp"
#include "recorders/data_record_po.hpp"

#include "topologies/temporal_space.hpp"
#include "interpolation/discrete_point_trajectory.hpp"

#include <boost/random/normal_distribution.hpp>
#include "base/global_rng.hpp"

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



  
typedef ReaK::ctrl::satellite3D_lin_dt_system::state_space_type      sat3D_state_space_type;
typedef ReaK::ctrl::satellite3D_lin_dt_system::point_type            sat3D_state_type;
typedef ReaK::ctrl::satellite3D_lin_dt_system::point_difference_type sat3D_state_diff_type;
typedef ReaK::ctrl::satellite3D_lin_dt_system::input_type            sat3D_input_type;
typedef ReaK::ctrl::satellite3D_lin_dt_system::output_type           sat3D_output_type;

typedef ReaK::pp::temporal_space<sat3D_state_space_type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> sat3D_temp_space_type;
typedef ReaK::pp::topology_traits< sat3D_temp_space_type >::point_type sat3D_temp_point_type;

typedef ReaK::ctrl::covariance_matrix< ReaK::vect_n<double> > cov_type;
typedef cov_type::matrix_type cov_matrix_type;
typedef ReaK::ctrl::gaussian_belief_state< sat3D_state_type,  cov_type > sat3D_state_belief_type;
typedef ReaK::ctrl::gaussian_belief_state< sat3D_input_type,  cov_type > sat3D_input_belief_type;
typedef ReaK::ctrl::gaussian_belief_state< sat3D_output_type, cov_type > sat3D_output_belief_type;



template <typename Sat3DSystemType>
std::vector< std::pair< double, ReaK::vect_n<double> > > batch_KF_on_timeseries(
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    const Sat3DSystemType& sat_sys,
    const sat3D_state_space_type& state_space,
    sat3D_state_belief_type b,
    sat3D_input_belief_type b_u,
    sat3D_output_belief_type b_z,
    std::size_t skips = 1) {
  using namespace ReaK;
  
  std::vector< std::pair< double, vect_n<double> > > result_pts;
  
  std::vector< std::pair< double, sat3D_state_type > >::const_iterator it_orig = ground_truth.begin();
  std::vector< std::pair< double, sat3D_measurement_point > >::const_iterator it = measurements.begin();
  while( it != measurements.end() ) {
    vect_n<double> z_vect(it->second.pose.size() + it->second.gyro.size() + it->second.IMU_a_m.size(), 0.0);
    z_vect[range(0,6)] = it->second.pose;
    if( it->second.gyro.size() ) {
      z_vect[range(7,9)] = it->second.gyro;
      if( it->second.IMU_a_m.size() )
        z_vect[range(10,15)] = it->second.IMU_a_m;
    };
    b_z.set_mean_state(z_vect);
    
    b_u.set_mean_state(it->second.u);
    
    ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z, it->first);
    
    const sat3D_state_type& x_mean = b.get_mean_state();
    vect_n<double> result_vect(13 + 12 + 12, 0.0);
    result_vect[range(0,2)]   = get_position(x_mean);
    result_vect[range(3,6)]   = get_quaternion(x_mean);
    result_vect[range(7,9)]   = get_velocity(x_mean);
    result_vect[range(10,12)] = get_ang_velocity(x_mean);
    
    if( it_orig != ground_truth.end() ) {
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) * get_quaternion(it_orig->second).as_rotation());
      result_vect[range(13,15)] = get_position(x_mean) - get_position(it_orig->second);
      result_vect[range(16,18)] = aa_diff.angle() * aa_diff.axis();
      result_vect[range(19,21)] = get_velocity(x_mean) - get_velocity(it_orig->second);
      result_vect[range(22,24)] = get_ang_velocity(x_mean) - get_ang_velocity(it_orig->second);
    };
    
    const cov_matrix_type& P_xx = b.get_covariance().get_matrix();
    for(std::size_t l = 0; l < 12; ++l)
      result_vect[l + 25] = P_xx(l,l);
    
    result_pts.push_back(std::make_pair(it->first, result_vect));
    
    for(std::size_t i = 0; i < skips; ++i) {
      if( it != measurements.end() )
        ++it;
      if( it_orig != ground_truth.end() )
        ++it_orig;
    };
  };
  
  return result_pts;
};



template <typename Sat3DSystemType>
void generate_timeseries(
    std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    const Sat3DSystemType& sat_sys,
    const sat3D_state_space_type& state_space,
    sat3D_state_type x,
    double start_time, double end_time,
    const cov_matrix_type& Qu,
    const cov_matrix_type& R,
    ReaK::shared_ptr< ReaK::recorder::data_recorder > stat_results = ReaK::shared_ptr< ReaK::recorder::data_recorder >() ) {
  using namespace ReaK;
  
  measurements.clear();
  ground_truth.clear();
  
  double time_step = sat_sys.get_time_step();
  std::vector< double > std_devs(R.get_row_count() + R.get_row_count() / 3, 0.0);
  for(double t = start_time; t < end_time; t += time_step) {
    vect_n<double> u(6, 0.0);
    u[0] = var_rnd() * sqrt(Qu(0,0));
    u[1] = var_rnd() * sqrt(Qu(1,1));
    u[2] = var_rnd() * sqrt(Qu(2,2));
    u[3] = var_rnd() * sqrt(Qu(3,3));
    u[4] = var_rnd() * sqrt(Qu(4,4));
    u[5] = var_rnd() * sqrt(Qu(5,5));
    
    x = sat_sys.get_next_state(state_space, x, u, t);
    vect_n<double> y = sat_sys.get_output(state_space, x, u, t);
    
    ground_truth.push_back(std::make_pair(t, x));
    
    sat3D_measurement_point meas;
    meas.u = vect_n<double>(6, 0.0);
    meas.pose = vect_n<double>(7, 0.0);
    meas.pose[0] = y[0] + var_rnd() * sqrt(R(0,0));
    meas.pose[1] = y[1] + var_rnd() * sqrt(R(1,1));
    meas.pose[2] = y[2] + var_rnd() * sqrt(R(2,2));
    
    vect<double,3> aa_vect_noise(var_rnd() * sqrt(R(3,3)), var_rnd() * sqrt(R(4,4)), var_rnd() * sqrt(R(5,5)));
    axis_angle<double> aa_noise(norm_2(aa_vect_noise), aa_vect_noise);
    quaternion<double> y_quat(vect<double,4>(y[3],y[4],y[5],y[6]));
    y_quat *= aa_noise.getQuaternion();
    meas.pose[range(3,6)] = vect<double,4>(y_quat[0], y_quat[1], y_quat[2], y_quat[3]);
    
    std::size_t k = ground_truth.size();
    if( stat_results ) {
      std_devs[0] = ( ( k - 1 ) * std_devs[0] + (meas.pose[0] - y[0]) * (meas.pose[0] - y[0]) ) / k;
      std_devs[1] = ( ( k - 1 ) * std_devs[1] + (meas.pose[1] - y[1]) * (meas.pose[1] - y[1]) ) / k;
      std_devs[2] = ( ( k - 1 ) * std_devs[2] + (meas.pose[2] - y[2]) * (meas.pose[2] - y[2]) ) / k;
      std_devs[3] = ( ( k - 1 ) * std_devs[3] + aa_vect_noise[0] * aa_vect_noise[0] ) / k;
      std_devs[4] = ( ( k - 1 ) * std_devs[4] + aa_vect_noise[1] * aa_vect_noise[1] ) / k;
      std_devs[5] = ( ( k - 1 ) * std_devs[5] + aa_vect_noise[2] * aa_vect_noise[2] ) / k;
      
      std_devs[6] = ( ( k - 1 ) * std_devs[6] + ( (meas.pose[0] - y[0]) * (meas.pose[0] - y[0])
                                                + (meas.pose[1] - y[1]) * (meas.pose[1] - y[1])
                                                + (meas.pose[2] - y[2]) * (meas.pose[2] - y[2]) ) ) / k;
      std_devs[7] = ( ( k - 1 ) * std_devs[7] + norm_2_sqr(aa_vect_noise) ) / k;
    };
    
    if( y.size() >= 10 ) {
      meas.gyro = y[range(7,9)];
      meas.gyro[0] += var_rnd() * sqrt(R(6,6));
      meas.gyro[1] += var_rnd() * sqrt(R(7,7));
      meas.gyro[2] += var_rnd() * sqrt(R(8,8));
      if( stat_results ) {
        std_devs[8]  = ( ( k - 1 ) * std_devs[8] + (meas.gyro[0] - y[7]) * (meas.gyro[0] - y[7]) ) / k;
        std_devs[9]  = ( ( k - 1 ) * std_devs[9] + (meas.gyro[1] - y[8]) * (meas.gyro[1] - y[8]) ) / k;
        std_devs[10] = ( ( k - 1 ) * std_devs[10] + (meas.gyro[2] - y[9]) * (meas.gyro[2] - y[9]) ) / k;
        std_devs[11] = ( ( k - 1 ) * std_devs[11] + ( (meas.gyro[0] - y[7]) * (meas.gyro[0] - y[7])
                                                    + (meas.gyro[1] - y[8]) * (meas.gyro[1] - y[8])
                                                    + (meas.gyro[2] - y[9]) * (meas.gyro[2] - y[9]) ) ) / k;
      };
      if( y.size() >= 16 ) {
        meas.IMU_a_m = y[range(10,15)];
        meas.IMU_a_m[0] += var_rnd() * sqrt(R(9,9));
        meas.IMU_a_m[1] += var_rnd() * sqrt(R(10,10));
        meas.IMU_a_m[2] += var_rnd() * sqrt(R(11,11));
        meas.IMU_a_m[3] += var_rnd() * sqrt(R(12,12));
        meas.IMU_a_m[4] += var_rnd() * sqrt(R(13,13));
        meas.IMU_a_m[5] += var_rnd() * sqrt(R(14,14));
        if( stat_results ) {
          std_devs[12] = ( ( k - 1 ) * std_devs[12] + (meas.IMU_a_m[0] - y[10]) * (meas.IMU_a_m[0] - y[10]) ) / k;
          std_devs[13] = ( ( k - 1 ) * std_devs[13] + (meas.IMU_a_m[1] - y[11]) * (meas.IMU_a_m[1] - y[11]) ) / k;
          std_devs[14] = ( ( k - 1 ) * std_devs[14] + (meas.IMU_a_m[2] - y[12]) * (meas.IMU_a_m[2] - y[12]) ) / k;
          std_devs[16] = ( ( k - 1 ) * std_devs[16] + (meas.IMU_a_m[3] - y[13]) * (meas.IMU_a_m[3] - y[13]) ) / k;
          std_devs[17] = ( ( k - 1 ) * std_devs[17] + (meas.IMU_a_m[4] - y[14]) * (meas.IMU_a_m[4] - y[14]) ) / k;
          std_devs[18] = ( ( k - 1 ) * std_devs[18] + (meas.IMU_a_m[5] - y[15]) * (meas.IMU_a_m[5] - y[15]) ) / k;
          
          std_devs[15] = ( ( k - 1 ) * std_devs[15] + ( (meas.IMU_a_m[0] - y[10]) * (meas.IMU_a_m[0] - y[10])
                                                      + (meas.IMU_a_m[1] - y[11]) * (meas.IMU_a_m[1] - y[11])
                                                      + (meas.IMU_a_m[2] - y[12]) * (meas.IMU_a_m[2] - y[12]) ) ) / k;
          std_devs[19] = ( ( k - 1 ) * std_devs[19] + ( (meas.IMU_a_m[3] - y[13]) * (meas.IMU_a_m[3] - y[13])
                                                      + (meas.IMU_a_m[4] - y[14]) * (meas.IMU_a_m[4] - y[14])
                                                      + (meas.IMU_a_m[5] - y[15]) * (meas.IMU_a_m[5] - y[15]) ) ) / k;
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
void do_all_single_runs(
    ReaK::recorder::data_stream_options output_opt,
    const std::string& filter_name,
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const sat3D_state_space_type& state_space,
    const sat3D_state_belief_type& b,
    sat3D_input_belief_type b_u,
    const sat3D_output_belief_type& b_z,
    double time_step, unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;
  
  cov_matrix_type Qu = b_u.get_covariance().get_matrix();
  
  for(unsigned int skips = min_skips; skips <= max_skips; ++skips) {
    
    sat_sys.set_time_step(skips * time_step);
    
    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));
    
    std::vector< std::pair< double, vect_n<double> > > result_pts = 
      batch_KF_on_timeseries(measurements, ground_truth, sat_sys, state_space, b, b_u, b_z, skips);
    
    std::stringstream ss;
    ss << "_" << std::setfill('0') << std::setw(4) << int(1000 * skips * time_step) << "_" << filter_name << ".";
    
    recorder::data_stream_options cur_out_opt = output_opt;
    cur_out_opt.file_name += ss.str() + cur_out_opt.get_extension();
    cur_out_opt.names.clear();
    cur_out_opt
      .add_name("time").add_name("pos_x").add_name("pos_y").add_name("pos_z")
      .add_name("q0").add_name("q1").add_name("q2").add_name("q3")
      .add_name("vel_x").add_name("vel_y").add_name("vel_z")
      .add_name("avel_x").add_name("avel_y").add_name("avel_z")
      .add_name("ep_x").add_name("ep_y").add_name("ep_z")
      .add_name("ea_x").add_name("ea_y").add_name("ea_z")
      .add_name("ev_x").add_name("ev_y").add_name("ev_z")
      .add_name("ew_x").add_name("ew_y").add_name("ew_z")
      .add_name("P_xx").add_name("P_yy").add_name("P_zz")
      .add_name("P_aax").add_name("P_aay").add_name("P_aaz")
      .add_name("P_vvx").add_name("P_vvy").add_name("P_vvz")
      .add_name("P_wwx").add_name("P_wwy").add_name("P_wwz");
    shared_ptr< recorder::data_recorder > results = cur_out_opt.create_recorder();
    for(std::size_t i = 0; i < result_pts.size(); ++i) {
      (*results) << result_pts[i].first;
      for(std::size_t j = 0; j < result_pts[i].second.size(); ++j)
        (*results) << result_pts[i].second[j];
      (*results) << recorder::data_recorder::end_value_row;
    };
    (*results) << recorder::data_recorder::flush;
  };
  
  sat_sys.set_time_step(time_step);
  
};





template <typename Sat3DSystemType>
void do_single_monte_carlo_run(
    std::map< std::string, ReaK::shared_ptr< ReaK::recorder::data_recorder > >& results_map,
    ReaK::recorder::data_stream_options output_opt,
    const std::string& filter_name,
    const std::vector< std::pair< double, sat3D_measurement_point > >& measurements,
    const std::vector< std::pair< double, sat3D_state_type > >& ground_truth,
    Sat3DSystemType& sat_sys,
    const sat3D_state_space_type& state_space,
    const sat3D_state_belief_type& b,
    sat3D_input_belief_type b_u,
    const sat3D_output_belief_type& b_z,
    double time_step, unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;
  
  cov_matrix_type Qu = b_u.get_covariance().get_matrix();
  
  for(unsigned int skips = min_skips; skips <= max_skips; ++skips) {
    
    sat_sys.set_time_step(skips * time_step);
    
    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));
    
    std::vector< std::pair< double, vect_n<double> > > result_pts = 
      batch_KF_on_timeseries(measurements, ground_truth, sat_sys, state_space, b, b_u, b_z, skips);
    
    std::vector< double > std_devs(28, 0.0);
    for(std::size_t i = 0; i < result_pts.size(); ++i) {
      
      // state vector component errors:
      for(std::size_t j = 13; j < 25; ++j)
        std_devs[j - 13] = ( i * std_devs[j - 13] + result_pts[i].second[j] * result_pts[i].second[j] ) / (i + 1);
      
      // position distance error:
      std_devs[12] = ( i * std_devs[12] + ( result_pts[i].second[13] * result_pts[i].second[13]
                                          + result_pts[i].second[14] * result_pts[i].second[14]
                                          + result_pts[i].second[15] * result_pts[i].second[15] ) ) / (i + 1);
      
      // rotation angle error:
      std_devs[13] = ( i * std_devs[13] + ( result_pts[i].second[16] * result_pts[i].second[16]
                                          + result_pts[i].second[17] * result_pts[i].second[17]
                                          + result_pts[i].second[18] * result_pts[i].second[18] ) ) / (i + 1);
      
      // speed error:
      std_devs[14] = ( i * std_devs[14] + ( result_pts[i].second[19] * result_pts[i].second[19]
                                          + result_pts[i].second[20] * result_pts[i].second[20]
                                          + result_pts[i].second[21] * result_pts[i].second[21] ) ) / (i + 1);
      
      // angular speed error:
      std_devs[15] = ( i * std_devs[15] + ( result_pts[i].second[22] * result_pts[i].second[22]
                                          + result_pts[i].second[23] * result_pts[i].second[23]
                                          + result_pts[i].second[24] * result_pts[i].second[24] ) ) / (i + 1);
      
      // average estimated covariances:
      for(std::size_t j = 25; j < 37; ++j)
        std_devs[j - 25 + 16] = ( i * std_devs[j - 25 + 16] + result_pts[i].second[j] ) / (i + 1);
      
    };
    
    std::stringstream ss;
    ss << "_" << std::setfill('0') << std::setw(4) << int(1000 * skips * time_step) << "_" << filter_name;
    std::string file_middle = ss.str();
    shared_ptr< recorder::data_recorder >& results = results_map[file_middle];
    if( !results ) {
      recorder::data_stream_options cur_out_opt = output_opt;
      cur_out_opt.file_name += file_middle + "_stddevs." + cur_out_opt.get_extension();
      cur_out_opt.names.clear();
      cur_out_opt
        .add_name("ep_x").add_name("ep_y").add_name("ep_z")
        .add_name("ea_x").add_name("ea_y").add_name("ea_z")
        .add_name("ev_x").add_name("ev_y").add_name("ev_z")
        .add_name("ew_x").add_name("ew_y").add_name("ew_z")
        .add_name("ep_m").add_name("ea_m").add_name("ev_m").add_name("ew_m")
        .add_name("P_xx").add_name("P_yy").add_name("P_zz")
        .add_name("P_aax").add_name("P_aay").add_name("P_aaz")
        .add_name("P_vvx").add_name("P_vvy").add_name("P_vvz")
        .add_name("P_wwx").add_name("P_wwy").add_name("P_wwz");
      results = cur_out_opt.create_recorder();
    };
    
    for(std::size_t j = 0; j < 18; ++j)
      (*results) << std::sqrt(std_devs[j]); // turn variances into std-devs.
    (*results) << recorder::data_recorder::end_value_row << recorder::data_recorder::flush;
    
  };
  
  sat_sys.set_time_step(time_step);
  
  
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
  
  
  recorder::data_stream_options data_in_opt;
  shared_ptr< recorder::data_extractor > data_in;
  std::vector<std::string> names_in;
  if(!vm.count("generate-meas")) {
    try {
      data_in_opt  = recorder::get_data_stream_options_from_po(vm, false);
      boost::tie(data_in, names_in) = data_in_opt.create_extractor();
    } catch(std::invalid_argument& e) {
      std::cerr << "Error! Creation of input data-stream failed! Invalid argument: " << e.what() << std::endl;
      return 2;
    };
  };
  
  recorder::data_stream_options data_out_opt;
  shared_ptr< recorder::data_recorder > data_out;
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
      last_dot = ( last_slash == std::string::npos ? 0 : last_slash );
      if( last_dot > last_slash ) 
        output_stem_name.erase(output_stem_name.begin() + last_dot, output_stem_name.end());
    };
    
    if(!vm.count("output-traj-file")) {
      data_out = data_out_opt.create_recorder();
    };
  } catch(std::invalid_argument& e) {
    std::cerr << "Error! Creation of output data-stream failed! Invalid argument: " << e.what() << std::endl;
    return 1;
  };
  recorder::data_stream_options data_out_stem_opt = data_out_opt;
  data_out_stem_opt.file_name = output_stem_name;
  
  ctrl::satellite_model_options sat_options;
  try {
    sat_options = ctrl::get_satellite_model_options_from_po(vm);
  } catch(std::exception& e) {
    std::cerr << "Error! Creation of satellite modeling options failed! With exception: " << e.what() << std::endl;
    return 2;
  };
  
  std::string sys_output_stem_name = vm["system-output"].as<std::string>();
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
  
  double start_time = vm["start-time"].as<double>();
  double end_time   = vm["end-time"].as<double>();
  
  unsigned int mc_runs    = vm["mc-runs"].as<unsigned int>();
  unsigned int min_skips  = vm["min-skips"].as<unsigned int>();
  unsigned int max_skips  = vm["max-skips"].as<unsigned int>();
  
  double RAq0 = (sat_options.artificial_noise(3,3) + sat_options.artificial_noise(4,4) + sat_options.artificial_noise(5,5)) / 12.0;
  
  
  std::vector< std::pair< double, sat3D_measurement_point > > measurements;
  std::vector< std::pair< double, sat3D_state_type > > ground_truth;
  if( (!vm.count("monte-carlo")) && data_in ) {
    try {
      recorder::named_value_row nvr_in = data_in->getFreshNamedValueRow();
      while(true) {
        (*data_in) >> nvr_in;
        
        double t = nvr_in[ names_in[0] ];
        std::vector<double> meas;
        for(std::size_t i = 1; i < names_in.size(); ++i)
          meas.push_back( nvr_in[ names_in[i] ] );
        
        sat3D_measurement_point meas_actual, meas_noisy;
        
        /* read off the position-orientation measurements. */
        if(meas.size() < 7) {
          RK_ERROR("The measurement file does not even contain the position and quaternion measurements!");
          return 4;
        };
        meas_actual.pose = vect_n<double>(meas.begin(), meas.begin() + 7);
        meas_noisy.pose = meas_actual.pose;
        if( sat_options.artificial_noise.get_row_count() >= 6 ) {
          meas_noisy.pose += vect_n<double>(
            var_rnd() * sqrt(sat_options.artificial_noise(0,0)),
            var_rnd() * sqrt(sat_options.artificial_noise(1,1)),
            var_rnd() * sqrt(sat_options.artificial_noise(2,2)),
            var_rnd() * sqrt(RAq0),
            var_rnd() * sqrt(0.25 * sat_options.artificial_noise(3,3)),
            var_rnd() * sqrt(0.25 * sat_options.artificial_noise(4,4)),
            var_rnd() * sqrt(0.25 * sat_options.artificial_noise(5,5))
          );
        };
        meas.erase(meas.begin(), meas.begin() + 7);
        
        /* read off the IMU/gyro angular velocity measurements. */
        if( sat_options.get_meas_error_count() >= 9 ) {
          if(meas.size() < 3) {
            RK_ERROR("The measurement file does not contain the angular velocity measurements!");
            return 4;
          };
          meas_actual.gyro = vect_n<double>(meas.begin(), meas.begin() + 3);
          meas_noisy.gyro  = meas_actual.gyro;
          if( sat_options.artificial_noise.get_row_count() >= 9 ) {
            meas_noisy.gyro += vect_n<double>(
              var_rnd() * sqrt(sat_options.artificial_noise(6,6)),
              var_rnd() * sqrt(sat_options.artificial_noise(7,7)),
              var_rnd() * sqrt(sat_options.artificial_noise(8,8))
            );
          };
          meas.erase(meas.begin(), meas.begin() + 3);
        };
        
        /* read off the IMU accel-mag measurements. */
        if( sat_options.get_meas_error_count() >= 15 ) {
          if(meas.size() < 6) {
            RK_ERROR("The measurement file does not contain the accelerometer and magnetometer measurements!");
            return 4;
          };
          meas_actual.IMU_a_m = vect_n<double>(meas.begin(), meas.begin() + 6);
          meas_noisy.IMU_a_m  = meas_actual.IMU_a_m;
          if( sat_options.artificial_noise.get_row_count() >= 15 ) {
            meas_noisy.IMU_a_m += vect_n<double>(
              var_rnd() * sqrt(sat_options.artificial_noise( 9, 9)),
              var_rnd() * sqrt(sat_options.artificial_noise(10,10)),
              var_rnd() * sqrt(sat_options.artificial_noise(11,11)),
              var_rnd() * sqrt(sat_options.artificial_noise(12,12)),
              var_rnd() * sqrt(sat_options.artificial_noise(13,13)),
              var_rnd() * sqrt(sat_options.artificial_noise(14,14))
            );
          };
          meas.erase(meas.begin(), meas.begin() + 6);
        };
        
        /* read off the input vector. */
        if(meas.size() < 6) {
          RK_ERROR("The measurement file does not contain the input force-torque vector measurements!");
          return 4;
        };
        meas_actual.u = vect_n<double>(meas.begin(), meas.begin() + 6);
        meas_noisy.u  = meas_actual.u;
        meas.erase(meas.begin(), meas.begin() + 6);
        
        /* now, the meas_actual and meas_noisy are fully formed. */
        measurements.push_back( std::make_pair(t, meas_noisy) );
        
        /* check if the file contains a ground-truth: */
        if(meas.size() >= 7) {
          sat3D_state_type x;
          set_position(x, vect<double,3>(meas[0],meas[1],meas[2]));
          set_quaternion(x, unit_quat<double>(meas[3],meas[4],meas[5],meas[6]));
          if(meas.size() >= 13) {
            set_velocity(x, vect<double,3>(meas[7],meas[8],meas[9]));
            set_ang_velocity(x, vect<double,3>(meas[10],meas[11],meas[12]));
          };
          ground_truth.push_back( std::make_pair(t, x) );
          meas.erase(meas.begin(), meas.end());
        } else if( sat_options.artificial_noise.get_row_count() >= 6 ) {
          sat3D_state_type x;
          set_position(x, vect<double,3>(meas_actual.pose[0],meas_actual.pose[1],meas_actual.pose[2]));
          set_quaternion(x, unit_quat<double>(meas_actual.pose[3],meas_actual.pose[4],meas_actual.pose[5],meas_actual.pose[6]));
          ground_truth.push_back( std::make_pair(t, x) );
        };
        
      };
    } catch(recorder::end_of_record&) { };
  };
  
  
  
  
  sat3D_temp_space_type sat_space(
    "satellite3D_temporal_space",
    sat3D_state_space_type(),
    pp::time_poisson_topology("satellite3D_time_space", sat_options.time_step, (end_time - start_time) * 0.5));
  
  typedef pp::discrete_point_trajectory< sat3D_temp_space_type > sat3D_traj_type;
  
  shared_ptr< sat3D_traj_type > traj_ptr;
  if( vm.count("generate-meas-file") && vm.count("output-traj-file") )
    traj_ptr = shared_ptr< sat3D_traj_type >(new sat3D_traj_type( shared_ptr< sat3D_temp_space_type >(&sat_space, null_deleter()) ));
  
  
  sat3D_state_type x_init;
  set_position(x_init, vect<double,3>(0.0, 0.0, 0.0));
  set_velocity(x_init, vect<double,3>(0.0, 0.0, 0.0));
  set_quaternion(x_init, unit_quat<double>());
  set_ang_velocity(x_init, vect<double,3>(0.0, 0.0, 0.0));
  
  sat3D_state_belief_type b_init(x_init, 
                                 cov_type(cov_matrix_type(mat<double,mat_structure::diagonal>(12,10.0))));
  
  sat3D_input_belief_type b_u(sat3D_input_type(vect_n<double>(6, 0.0)),  
                              cov_type(cov_matrix_type(sat_options.input_disturbance)));
  
  if( !vm.count("gyro") && !vm.count("IMU") ) {
    
    // Create the set of satellite3D systems for when there is only pose measurements:
    
    ctrl::satellite3D_inv_dt_system sat3D_inv(
      "satellite3D_inv", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step);
    
    ctrl::satellite3D_imdt_sys sat3D_invmom(
      "satellite3D_invmom", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step);
    
    ctrl::satellite3D_imdt_sys sat3D_invmid(
      "satellite3D_invmid", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step, 2);
    
    sat3D_output_belief_type b_z(sat3D_output_type(vect_n<double>(7, 0.0)), 
                                 cov_type(cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise)));
    
    
    if( vm.count("generate-mdl-files") ) {
      try {
        shared_ptr< ctrl::satellite3D_inv_dt_system > satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_inv, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_inv_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
        
        satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmom, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmom_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
        
        satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmid, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmid_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch(...) {
        RK_ERROR("An exception occurred during the saving the satellite system file!");
        return 14;
      };
    } else if( ! vm.count("monte-carlo") ) {
      
      if( measurements.size() == 0 ) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, sat_options.initial_motion);
        generate_timeseries(measurements, ground_truth, sat3D_inv, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise));
        
        // and output those if asked for it:
        if( vm.count("generate-meas-file") ) {
          recorder::data_stream_options data_meas_opt = data_out_opt;
          data_meas_opt.file_name = output_stem_name + "_meas." + data_meas_opt.get_extension();
          data_meas_opt.names.clear();
          data_meas_opt
            .add_name("time").add_name("p_x").add_name("p_y").add_name("p_z")
            .add_name("q_0").add_name("q_1").add_name("q_2").add_name("q_3")
            .add_name("f_x").add_name("f_y").add_name("f_z")
            .add_name("t_x").add_name("t_y").add_name("t_z")
            .add_name("p_x_true").add_name("p_y_true").add_name("p_z_true")
            .add_name("q_0_true").add_name("q_1_true").add_name("q_2_true").add_name("q_3_true")
            .add_name("v_x_true").add_name("v_y_true").add_name("v_z_true")
            .add_name("w_x_true").add_name("w_y_true").add_name("w_z_true"); 
          shared_ptr< recorder::data_recorder > data_meas = data_meas_opt.create_recorder();
          for(std::size_t i = 0; i < measurements.size(); ++i) {
            (*data_meas) << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            (*data_meas) << m.pose << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            (*data_meas) << get_position(g) << get_quaternion(g) << get_velocity(g) << get_ang_velocity(g);
            (*data_meas) << recorder::data_recorder::end_value_row;
            if( vm.count("output-traj-file") ) {
              traj_ptr->push_back(sat3D_temp_point_type(ground_truth[i].first, g));
            };
          };
          (*data_meas) << recorder::data_recorder::flush;
        };
      };
      
      // do a single run for each skips:
      
      std::cout << "Running estimators on data series.." << std::flush;
      
      if( vm.count("iekf") ) {
        do_all_single_runs(data_out_stem_opt, "iekf", measurements, ground_truth, sat3D_inv, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkf") ) {
        do_all_single_runs(data_out_stem_opt, "imkf", measurements, ground_truth, sat3D_invmom, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkfv2") ) {
        do_all_single_runs(data_out_stem_opt, "imkfv2", measurements, ground_truth, sat3D_invmid, 
                           sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      std::cout << "Finished!" << std::endl;
      
    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, sat_options.initial_motion);
      
      recorder::data_stream_options data_stddev_opt = data_out_opt;
      data_stddev_opt.file_name = output_stem_name + "_meas_stddevs." + data_stddev_opt.get_extension();
      data_stddev_opt.names.clear();
      data_stddev_opt
        .add_name("ep_x").add_name("ep_y").add_name("ep_z")
        .add_name("ea_x").add_name("ea_y").add_name("ea_z").add_name("ep_m").add_name("ea_m"); 
      shared_ptr< recorder::data_recorder > data_stddev = data_stddev_opt.create_recorder();
      
      std::map< std::string, shared_ptr< recorder::data_recorder > > results_map;
      
      std::cout << "Running Monte-Carlo Simulations..." << std::endl;
      
      for(unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {
        
        std::cout << "\r" << std::setw(10) << mc_i << std::flush;
        
        generate_timeseries(measurements, ground_truth, sat3D_inv, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise), data_stddev);
        
        std::cout << "." << std::flush;
        
        if( vm.count("iekf") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "iekf", measurements, ground_truth, sat3D_inv, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkf") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkf", measurements, ground_truth, sat3D_invmom, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkfv2") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkfv2", measurements, ground_truth, sat3D_invmid, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
      };
      
      std::cout << "Finished!" << std::endl;
      
    };
    
    
  } else if( vm.count("gyro") && !vm.count("IMU") ) {
    
    // Create the set of satellite3D systems for when there is gyro measurements:
    
    ctrl::satellite3D_gyro_inv_dt_system sat3D_inv_gyro(
      "satellite3D_inv_with_gyros", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step);
    
    ctrl::satellite3D_gyro_imdt_sys sat3D_invmom_gyro(
      "satellite3D_invmom_with_gyros", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step);
    
    ctrl::satellite3D_gyro_imdt_sys sat3D_invmid_gyro(
      "satellite3D_invmid_with_gyros", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step, 2);
    
    sat3D_output_belief_type b_z(sat3D_output_type(vect_n<double>(10, 0.0)), 
                                 cov_type(cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise)));
    
    
    if( vm.count("generate-mdl-files") ) {
      try {
        shared_ptr< ctrl::satellite3D_inv_dt_system > satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_inv_gyro, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_inv_gyro_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
        
        satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmom_gyro, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmom_gyro_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
        
        satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmid_gyro, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmid_gyro_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch(...) {
        RK_ERROR("An exception occurred during the saving the satellite system file!");
        return 14;
      };
    } else if( ! vm.count("monte-carlo") ) {
      
      if( measurements.size() == 0 ) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, sat_options.initial_motion);
        generate_timeseries(measurements, ground_truth, sat3D_inv_gyro, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise));
        
        // and output those if asked for it:
        if( vm.count("generate-meas-file") ) {
          recorder::data_stream_options data_meas_opt = data_out_opt;
          data_meas_opt.file_name = output_stem_name + "_meas." + data_meas_opt.get_extension();
          data_meas_opt.names.clear();
          data_meas_opt
            .add_name("time").add_name("p_x").add_name("p_y").add_name("p_z")
            .add_name("q_0").add_name("q_1").add_name("q_2").add_name("q_3")
            .add_name("w_x").add_name("w_y").add_name("w_z")
            .add_name("f_x").add_name("f_y").add_name("f_z")
            .add_name("t_x").add_name("t_y").add_name("t_z")
            .add_name("p_x_true").add_name("p_y_true").add_name("p_z_true")
            .add_name("q_0_true").add_name("q_1_true").add_name("q_2_true").add_name("q_3_true")
            .add_name("v_x_true").add_name("v_y_true").add_name("v_z_true")
            .add_name("w_x_true").add_name("w_y_true").add_name("w_z_true"); 
          shared_ptr< recorder::data_recorder > data_meas = data_meas_opt.create_recorder();
          for(std::size_t i = 0; i < measurements.size(); ++i) {
            (*data_meas) << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            (*data_meas) << m.pose << m.gyro << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            (*data_meas) << get_position(g) << get_quaternion(g) << get_velocity(g) << get_ang_velocity(g);
            (*data_meas) << recorder::data_recorder::end_value_row;
            if( vm.count("output-traj-file") ) {
              traj_ptr->push_back(sat3D_temp_point_type(ground_truth[i].first, g));
            };
          };
          (*data_meas) << recorder::data_recorder::flush;
        };
      };
      
      // do a single run for each skips:
      
      std::cout << "Running estimators on data series.." << std::flush;
      
      if( vm.count("iekf") ) {
        do_all_single_runs(data_out_stem_opt, "iekf_gyro", measurements, ground_truth, sat3D_inv_gyro, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkf") ) {
        do_all_single_runs(data_out_stem_opt, "imkf_gyro", measurements, ground_truth, sat3D_invmom_gyro, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkfv2") ) {
        do_all_single_runs(data_out_stem_opt, "imkfv2_gyro", measurements, ground_truth, sat3D_invmid_gyro, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      std::cout << "Finished!" << std::endl;
      
    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, sat_options.initial_motion);
      
      recorder::data_stream_options data_stddev_opt = data_out_opt;
      data_stddev_opt.file_name = output_stem_name + "_meas_gyro_stddevs." + data_stddev_opt.get_extension();
      data_stddev_opt.names.clear();
      data_stddev_opt
        .add_name("ep_x").add_name("ep_y").add_name("ep_z")
        .add_name("ea_x").add_name("ea_y").add_name("ea_z")
        .add_name("ep_m").add_name("ea_m")
        .add_name("ew_x").add_name("ew_y").add_name("ew_z").add_name("ew_m");
      shared_ptr< recorder::data_recorder > data_stddev = data_stddev_opt.create_recorder();
      
      std::map< std::string, shared_ptr< recorder::data_recorder > > results_map;
      
      std::cout << "Running Monte-Carlo Simulations..." << std::endl;
      
      for(unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {
        
        std::cout << "\r" << std::setw(10) << mc_i << std::flush;
        
        generate_timeseries(measurements, ground_truth, sat3D_inv_gyro, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise), data_stddev);
        
        std::cout << "." << std::flush;
        
        if( vm.count("iekf") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "iekf_gyro", measurements, ground_truth, sat3D_inv_gyro, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkf") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkf_gyro", measurements, ground_truth, sat3D_invmom_gyro, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkfv2") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkfv2_gyro", measurements, ground_truth, sat3D_invmid_gyro, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
      };
      
      std::cout << "Finished!" << std::endl;
      
    };
    
  } else {
    
    // Create the set of satellite3D systems for when there is IMU measurements:
    
    ctrl::satellite3D_IMU_imdt_sys sat3D_invmom_IMU(
      "satellite3D_invmom_with_IMU", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step,
      sat_options.IMU_orientation, sat_options.IMU_location, sat_options.earth_orientation, sat_options.mag_field_direction);
    
    ctrl::satellite3D_IMU_imdt_sys sat3D_invmid_IMU(
      "satellite3D_invmid_with_IMU", sat_options.mass, sat_options.inertia_tensor, sat_options.time_step,
      sat_options.IMU_orientation, sat_options.IMU_location, sat_options.earth_orientation, sat_options.mag_field_direction, 2);
    
    sat3D_output_belief_type b_z(sat3D_output_type(vect_n<double>(16, 0.0)), 
                                 cov_type(cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise)));
    
    
    if( vm.count("generate-mdl-files") ) {
      try {
        shared_ptr< ctrl::satellite3D_inv_dt_system > satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmom_IMU, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmom_IMU_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
        
        satellite3D_system = 
          shared_ptr< ctrl::satellite3D_inv_dt_system >(&sat3D_invmid_IMU, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_invmid_IMU_mdl.rkx"))
          & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch(...) {
        RK_ERROR("An exception occurred during the saving the satellite system file!");
        return 14;
      };
    } else if( ! vm.count("monte-carlo") ) {
      
      if( measurements.size() == 0 ) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, sat_options.initial_motion);
        generate_timeseries(measurements, ground_truth, sat3D_invmom_IMU, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise));
        
        // and output those if asked for it:
        if( vm.count("generate-meas-file") ) {
          recorder::data_stream_options data_meas_opt = data_out_opt;
          data_meas_opt.file_name = output_stem_name + "_meas." + data_meas_opt.get_extension();
          data_meas_opt.names.clear();
          data_meas_opt
            .add_name("time").add_name("p_x").add_name("p_y").add_name("p_z")
            .add_name("q_0").add_name("q_1").add_name("q_2").add_name("q_3")
            .add_name("w_x").add_name("w_y").add_name("w_z")
            .add_name("acc_x").add_name("acc_y").add_name("acc_z")
            .add_name("mag_x").add_name("mag_y").add_name("mag_z")
            .add_name("f_x").add_name("f_y").add_name("f_z")
            .add_name("t_x").add_name("t_y").add_name("t_z")
            .add_name("p_x_true").add_name("p_y_true").add_name("p_z_true")
            .add_name("q_0_true").add_name("q_1_true").add_name("q_2_true").add_name("q_3_true")
            .add_name("v_x_true").add_name("v_y_true").add_name("v_z_true")
            .add_name("w_x_true").add_name("w_y_true").add_name("w_z_true"); 
          shared_ptr< recorder::data_recorder > data_meas = data_meas_opt.create_recorder();
          for(std::size_t i = 0; i < measurements.size(); ++i) {
            (*data_meas) << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            (*data_meas) << m.pose << m.gyro << m.IMU_a_m << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            (*data_meas) << get_position(g) << get_quaternion(g) << get_velocity(g) << get_ang_velocity(g);
            (*data_meas) << recorder::data_recorder::end_value_row;
            if( vm.count("output-traj-file") ) {
              traj_ptr->push_back(sat3D_temp_point_type(ground_truth[i].first, g));
            };
          };
          (*data_meas) << recorder::data_recorder::flush;
        };
      };
      
      // do a single run for each skips:
      
      std::cout << "Running estimators on data series.." << std::flush;
      
      if( vm.count("iekf") ) {
        std::cerr << "Warning: The invariant extended Kalman filter (IEKF) is not available for full IMU measurements!" << std::endl;
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkf") ) {
        do_all_single_runs(data_out_stem_opt, "imkf_IMU", measurements, ground_truth, sat3D_invmom_IMU, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      if( vm.count("imkfv2") ) {
        do_all_single_runs(data_out_stem_opt, "imkfv2_IMU", measurements, ground_truth, sat3D_invmid_IMU, 
                          sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
      };
      
      std::cout << "." << std::flush;
      
      std::cout << "Finished!" << std::endl;
      
    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, sat_options.initial_motion);
      
      recorder::data_stream_options data_stddev_opt = data_out_opt;
      data_stddev_opt.file_name = output_stem_name + "_meas_gyro_stddevs." + data_stddev_opt.get_extension();
      data_stddev_opt.names.clear();
      data_stddev_opt
        .add_name("ep_x").add_name("ep_y").add_name("ep_z")
        .add_name("ea_x").add_name("ea_y").add_name("ea_z")
        .add_name("ep_m").add_name("ea_m")
        .add_name("ew_x").add_name("ew_y").add_name("ew_z").add_name("ew_m")
        .add_name("eacc_x").add_name("eacc_y").add_name("eacc_z").add_name("eacc_m")
        .add_name("emag_x").add_name("emag_y").add_name("emag_z").add_name("emag_m");
      shared_ptr< recorder::data_recorder > data_stddev = data_stddev_opt.create_recorder();
      
      std::map< std::string, shared_ptr< recorder::data_recorder > > results_map;
      
      if( vm.count("iekf") ) {
        std::cerr << "Warning: The invariant extended Kalman filter (IEKF) is not available for full IMU measurements!" << std::endl;
      };
      
      std::cout << "Running Monte-Carlo Simulations..." << std::endl;
      
      for(unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {
        
        std::cout << "\r" << std::setw(10) << mc_i << std::flush;
        
        generate_timeseries(measurements, ground_truth, sat3D_invmom_IMU, sat_space.get_space_topology(),
                            x_init, start_time, end_time, cov_matrix_type(sat_options.input_disturbance), 
                            cov_matrix_type(sat_options.measurement_noise + sat_options.artificial_noise), data_stddev);
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkf") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkf_IMU", measurements, ground_truth, sat3D_invmom_IMU, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
        if( vm.count("imkfv2") ) {
          do_single_monte_carlo_run(
            results_map, data_out_stem_opt, "imkfv2_IMU", measurements, ground_truth, sat3D_invmid_IMU, 
            sat_space.get_space_topology(), b_init, b_u, b_z, sat_options.time_step, min_skips, max_skips);
        };
        
        std::cout << "." << std::flush;
        
      };
      
      std::cout << "Finished!" << std::endl;
      
    };
    
  };
  
  if( vm.count("generate-meas-file") && ( vm.count("xml") + vm.count("protobuf") + vm.count("binary") > 0 ) ) {
    std::cout << "Saving the generated trajectory.." << std::flush;
    
    std::cout << "." << std::flush;
    
    if( vm.count("xml") ) {
      
      *(serialization::open_oarchive(output_stem_name + "_traj.rkx"))
        & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
      
    };
    
    std::cout << "." << std::flush;
    
    if( vm.count("protobuf") ) {
      
      *(serialization::open_oarchive(output_stem_name + "_traj.pbuf"))
        & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
      
    };
    
    std::cout << "." << std::flush;
    
    if( vm.count("binary") ) {
      
      *(serialization::open_oarchive(output_stem_name + "_traj.rkb"))
        & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
      
    };
    
    std::cout << "Finished!" << std::endl;
  };
  
  
};



