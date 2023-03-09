
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include "ReaK/core/base/global_rng.h"

#include "ReaK/control/systems/satellite_invar_models.h"

#include "ReaK/control/estimators/invariant_kalman_filter.h"
#include "ReaK/control/estimators/kalman_filter.h"

#include "ReaK/control/estimators/covariance_matrix.h"
#include "ReaK/control/estimators/gaussian_belief_state.h"

#include "ReaK/core/recorders/ascii_recorder.h"
#include "ReaK/core/serialization/archiver_factory.h"

#include "ReaK/topologies/interpolation/discrete_point_trajectory.h"
#include "ReaK/topologies/spaces/temporal_space.h"

#include <memory>
#include <random>

#include <cctype>
#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/log.h"

// I/O options
ABSL_FLAG(std::string, measurements, "",
          "Specify the filename for the satellite's measurements.");
ABSL_FLAG(std::string, init, "models/satellite3D_init.rkx",
          "Specify the filename for the satellite's initial conditions, only "
          "used when Monte-Carlo simulations are done (default is "
          "'models/satellite3D_init.rkx').");
ABSL_FLAG(std::string, inertia, "models/satellite3D_inertia.rkx",
          "Specify the filename for the satellite's inertial data (default is "
          "'models/satellite3D_inertia.rkx').");
ABSL_FLAG(std::string, Q_matrix, "models/satellite3D_Q.rkx",
          "Specify the filename for the satellite's input disturbance "
          "covariance matrix (default is 'models/satellite3D_Q.rkx').");
ABSL_FLAG(std::string, R_matrix, "models/satellite3D_R.rkx",
          "Specify the filename for the satellite's measurement noise "
          "covariance matrix (default is 'models/satellite3D_R.rkx').");
ABSL_FLAG(std::string, R_added, "",
          "Specify the filename for the satellite's artificial measurement "
          "noise covariance matrix.");
ABSL_FLAG(std::string, IMU_config, "models/satellite3D_IMU_config.rkx",
          "Specify the filename for the satellite's IMU configuration data, "
          "specifying its placement on the satellite and the inertial / "
          "magnetic-field frame it is relative to (default is "
          "'models/satellite3D_IMU_config.rkx').");
ABSL_FLAG(std::string, output, "est_results/satellite3D/output_record",
          "Specify the filename stem (without extension) for the output of the "
          "results (default is 'sim_results/satellite3D/output_record').");
ABSL_FLAG(bool, generate_meas_file, false,
          "If set, the measurement file will be generated from the output of a "
          "simulation with the given initial conditions (default is not).");
ABSL_FLAG(bool, generate_mdl_files, false,
          "If set, the output will be the generation of all the modeling files "
          "(with default values).");
ABSL_FLAG(std::string, system_output, "models/satellite3D",
          "Specify the filename-stem for the output of the satellite system, "
          "when 'generate-files' is set (default is 'models/satellite3D').");

// Simulation options
ABSL_FLAG(double, start_time, 0.0,
          "Start time of the estimation (default is 0.0).");
ABSL_FLAG(double, end_time, 1.0,
          "End time of the estimation (default is 1.0).");
ABSL_FLAG(double, time_step, 0.01,
          "Time-step in the measurement files (default is 0.01).");
ABSL_FLAG(bool, monte_carlo, false,
          "If set, will perform a Monte-Carlo set of randomized runs to gather "
          "estimation performance statistics.");
ABSL_FLAG(int, mc_runs, 1000,
          "Number of Monte-Carlo runs to perform (default is 1000).");
ABSL_FLAG(int, min_skips, 1,
          "Minimum number of time-step skips between estimations when "
          "generating a series of Monte-Carlo statistics (default is 1, i.e., "
          "one estimation point per measurement point).");
ABSL_FLAG(int, max_skips, 1,
          "maximum number of time-step skips between estimations when "
          "generating a series of Monte-Carlo statistics (default is 1, i.e., "
          "one estimation point per measurement point).");

// Modeling options
ABSL_FLAG(bool, gyro, false,
          "If set, a set of gyros is added to the model (angular velocity "
          "measurements). This requires the 'R-matrix' file to contain a 9x9 "
          "matrix.");
ABSL_FLAG(bool, IMU, false,
          "If set, a set of gyros is added to the model (angular velocity, "
          "magnetic field, and accelerometer measurements). This requires the "
          "'R-matrix' file to contain a 15x15 matrix. This option also "
          "automatically implies the 'midpoint' option. This option will "
          "trigger the use of the 'IMU-config' file to obtain the information "
          "necessary about the IMU and the Earth's inertial frame.");
ABSL_FLAG(bool, iekf, false,
          "If set, results for the invariant extended Kalman filter (IEKF) "
          "will be generated.");
ABSL_FLAG(bool, imkf, false,
          "If set, results for the invariant momentum-tracking Kalman filter "
          "(IMKF) will be generated.");
ABSL_FLAG(bool, imkfv2, false,
          "If set, results for the invariant midpoint Kalman filter (IMKFv2) "
          "will be generated.");

// Output options (at least one must be set)
ABSL_FLAG(bool, xml, false, "If set, output results in XML format (rkx).");
ABSL_FLAG(bool, protobuf, false,
          "If set, output results in protobuf format (pb).");
ABSL_FLAG(bool, binary, false, "If set, output results in binary format (rkb)");
ABSL_FLAG(
    bool, ssv, false,
    "If set, output resulting trajectories as time-series in "
    "space-separated-values files (ssv) (easily loadable in matlab / octave).");

namespace fs = std::filesystem;

const auto var_rnd = []() {
  std::normal_distribution<double> dist;
  return dist(ReaK::get_global_rng());
};

struct sat3D_measurement_point {
  ReaK::vect_n<double> pose;
  ReaK::vect_n<double> gyro;
  ReaK::vect_n<double> IMU_a_m;
  ReaK::vect_n<double> u;
};

using sat3D_state_space_type =
    ReaK::ctrl::satellite3D_lin_dt_system::state_space_type;
using sat3D_state_type = ReaK::ctrl::satellite3D_lin_dt_system::point_type;
using sat3D_state_diff_type =
    ReaK::ctrl::satellite3D_lin_dt_system::point_difference_type;
using sat3D_input_type = ReaK::ctrl::satellite3D_lin_dt_system::input_type;
using sat3D_output_type = ReaK::ctrl::satellite3D_lin_dt_system::output_type;

using sat3D_temp_space_type =
    ReaK::pp::temporal_space<sat3D_state_space_type,
                             ReaK::pp::time_poisson_topology,
                             ReaK::pp::time_distance_only>;
using sat3D_temp_point_type =
    ReaK::pp::topology_traits<sat3D_temp_space_type>::point_type;

using cov_type = ReaK::ctrl::covariance_matrix<ReaK::vect_n<double>>;
using cov_matrix_type = cov_type::matrix_type;
using sat3D_state_belief_type =
    ReaK::ctrl::gaussian_belief_state<sat3D_state_type, cov_type>;
using sat3D_input_belief_type =
    ReaK::ctrl::gaussian_belief_state<sat3D_input_type, cov_type>;
using sat3D_output_belief_type =
    ReaK::ctrl::gaussian_belief_state<sat3D_output_type, cov_type>;

template <typename Sat3DSystemType>
std::vector<std::pair<double, ReaK::vect_n<double>>> batch_KF_on_timeseries(
    const std::vector<std::pair<double, sat3D_measurement_point>>& measurements,
    const std::vector<std::pair<double, sat3D_state_type>>& ground_truth,
    const Sat3DSystemType& sat_sys, const sat3D_state_space_type& state_space,
    sat3D_state_belief_type b, sat3D_input_belief_type b_u,
    sat3D_output_belief_type b_z, std::size_t skips = 1) {
  using namespace ReaK;

  std::vector<std::pair<double, vect_n<double>>> result_pts;

  auto it_orig = ground_truth.begin();
  auto it = measurements.begin();
  while (it != measurements.end()) {
    vect_n<double> z_vect(it->second.pose.size() + it->second.gyro.size() +
                              it->second.IMU_a_m.size(),
                          0.0);
    z_vect[range(0, 7)] = it->second.pose;
    if (!it->second.gyro.empty()) {
      z_vect[range(7, 10)] = it->second.gyro;
      if (!it->second.IMU_a_m.empty()) {
        z_vect[range(10, 16)] = it->second.IMU_a_m;
      }
    }
    b_z.set_mean_state(z_vect);

    b_u.set_mean_state(it->second.u);

    ctrl::invariant_kalman_filter_step(sat_sys, state_space, b, b_u, b_z,
                                       it->first);

    const sat3D_state_type& x_mean = b.get_mean_state();
    vect_n<double> result_vect(13 + 12 + 12, 0.0);
    result_vect[range(0, 3)] = get_position(x_mean);
    result_vect[range(3, 7)] = get_quaternion(x_mean);
    result_vect[range(7, 10)] = get_velocity(x_mean);
    result_vect[range(10, 13)] = get_ang_velocity(x_mean);

    if (it_orig != ground_truth.end()) {
      axis_angle<double> aa_diff(invert(get_quaternion(x_mean).as_rotation()) *
                                 get_quaternion(it_orig->second).as_rotation());
      result_vect[range(13, 16)] =
          get_position(x_mean) - get_position(it_orig->second);
      result_vect[range(16, 19)] = aa_diff.angle() * aa_diff.axis();
      result_vect[range(19, 22)] =
          get_velocity(x_mean) - get_velocity(it_orig->second);
      result_vect[range(22, 25)] =
          get_ang_velocity(x_mean) - get_ang_velocity(it_orig->second);
    }

    const cov_matrix_type& P_xx = b.get_covariance().get_matrix();
    for (std::size_t l = 0; l < 12; ++l) {
      result_vect[l + 25] = P_xx(l, l);
    }

    result_pts.emplace_back(it->first, result_vect);

    for (std::size_t i = 0; i < skips; ++i) {
      if (it != measurements.end()) {
        ++it;
      }
      if (it_orig != ground_truth.end()) {
        ++it_orig;
      }
    }
  }

  return result_pts;
}

template <typename Sat3DSystemType>
void generate_timeseries(
    std::vector<std::pair<double, sat3D_measurement_point>>& measurements,
    std::vector<std::pair<double, sat3D_state_type>>& ground_truth,
    const Sat3DSystemType& sat_sys, const sat3D_state_space_type& state_space,
    sat3D_state_type x, double start_time, double end_time,
    const cov_matrix_type& Qu, const cov_matrix_type& R,
    std::shared_ptr<ReaK::recorder::ascii_recorder> stat_results =
        std::shared_ptr<ReaK::recorder::ascii_recorder>()) {
  using namespace ReaK;

  measurements.clear();
  ground_truth.clear();

  double time_step = sat_sys.get_time_step();
  std::vector<double> std_devs(R.get_row_count() + R.get_row_count() / 3, 0.0);
  for (double t = start_time; t < end_time; t += time_step) {
    vect_n<double> u(6, 0.0);
    u[0] = var_rnd() * sqrt(Qu(0, 0));
    u[1] = var_rnd() * sqrt(Qu(1, 1));
    u[2] = var_rnd() * sqrt(Qu(2, 2));
    u[3] = var_rnd() * sqrt(Qu(3, 3));
    u[4] = var_rnd() * sqrt(Qu(4, 4));
    u[5] = var_rnd() * sqrt(Qu(5, 5));

    x = sat_sys.get_next_state(state_space, x, u, t);
    vect_n<double> y = sat_sys.get_output(state_space, x, u, t);

    ground_truth.emplace_back(t, x);

    sat3D_measurement_point meas;
    meas.u = vect_n<double>(6, 0.0);
    meas.pose = vect_n<double>(7, 0.0);
    meas.pose[0] = y[0] + var_rnd() * sqrt(R(0, 0));
    meas.pose[1] = y[1] + var_rnd() * sqrt(R(1, 1));
    meas.pose[2] = y[2] + var_rnd() * sqrt(R(2, 2));

    vect<double, 3> aa_vect_noise(var_rnd() * sqrt(R(3, 3)),
                                  var_rnd() * sqrt(R(4, 4)),
                                  var_rnd() * sqrt(R(5, 5)));
    axis_angle<double> aa_noise(norm_2(aa_vect_noise), aa_vect_noise);
    quaternion<double> y_quat(vect<double, 4>(y[3], y[4], y[5], y[6]));
    y_quat *= aa_noise.getQuaternion();
    meas.pose[range(3, 7)] =
        vect<double, 4>(y_quat[0], y_quat[1], y_quat[2], y_quat[3]);

    std::size_t k = ground_truth.size();
    if (stat_results) {
      std_devs[0] = ((k - 1) * std_devs[0] +
                     (meas.pose[0] - y[0]) * (meas.pose[0] - y[0])) /
                    k;
      std_devs[1] = ((k - 1) * std_devs[1] +
                     (meas.pose[1] - y[1]) * (meas.pose[1] - y[1])) /
                    k;
      std_devs[2] = ((k - 1) * std_devs[2] +
                     (meas.pose[2] - y[2]) * (meas.pose[2] - y[2])) /
                    k;
      std_devs[3] =
          ((k - 1) * std_devs[3] + aa_vect_noise[0] * aa_vect_noise[0]) / k;
      std_devs[4] =
          ((k - 1) * std_devs[4] + aa_vect_noise[1] * aa_vect_noise[1]) / k;
      std_devs[5] =
          ((k - 1) * std_devs[5] + aa_vect_noise[2] * aa_vect_noise[2]) / k;

      std_devs[6] = ((k - 1) * std_devs[6] +
                     ((meas.pose[0] - y[0]) * (meas.pose[0] - y[0]) +
                      (meas.pose[1] - y[1]) * (meas.pose[1] - y[1]) +
                      (meas.pose[2] - y[2]) * (meas.pose[2] - y[2]))) /
                    k;
      std_devs[7] = ((k - 1) * std_devs[7] + norm_2_sqr(aa_vect_noise)) / k;
    }

    if (y.size() >= 10) {
      meas.gyro = y[range(7, 10)];
      meas.gyro[0] += var_rnd() * sqrt(R(6, 6));
      meas.gyro[1] += var_rnd() * sqrt(R(7, 7));
      meas.gyro[2] += var_rnd() * sqrt(R(8, 8));
      if (stat_results) {
        std_devs[8] = ((k - 1) * std_devs[8] +
                       (meas.gyro[0] - y[7]) * (meas.gyro[0] - y[7])) /
                      k;
        std_devs[9] = ((k - 1) * std_devs[9] +
                       (meas.gyro[1] - y[8]) * (meas.gyro[1] - y[8])) /
                      k;
        std_devs[10] = ((k - 1) * std_devs[10] +
                        (meas.gyro[2] - y[9]) * (meas.gyro[2] - y[9])) /
                       k;
        std_devs[11] = ((k - 1) * std_devs[11] +
                        ((meas.gyro[0] - y[7]) * (meas.gyro[0] - y[7]) +
                         (meas.gyro[1] - y[8]) * (meas.gyro[1] - y[8]) +
                         (meas.gyro[2] - y[9]) * (meas.gyro[2] - y[9]))) /
                       k;
      }
      if (y.size() >= 16) {
        meas.IMU_a_m = y[range(10, 16)];
        meas.IMU_a_m[0] += var_rnd() * sqrt(R(9, 9));
        meas.IMU_a_m[1] += var_rnd() * sqrt(R(10, 10));
        meas.IMU_a_m[2] += var_rnd() * sqrt(R(11, 11));
        meas.IMU_a_m[3] += var_rnd() * sqrt(R(12, 12));
        meas.IMU_a_m[4] += var_rnd() * sqrt(R(13, 13));
        meas.IMU_a_m[5] += var_rnd() * sqrt(R(14, 14));
        if (stat_results) {
          std_devs[12] =
              ((k - 1) * std_devs[12] +
               (meas.IMU_a_m[0] - y[10]) * (meas.IMU_a_m[0] - y[10])) /
              k;
          std_devs[13] =
              ((k - 1) * std_devs[13] +
               (meas.IMU_a_m[1] - y[11]) * (meas.IMU_a_m[1] - y[11])) /
              k;
          std_devs[14] =
              ((k - 1) * std_devs[14] +
               (meas.IMU_a_m[2] - y[12]) * (meas.IMU_a_m[2] - y[12])) /
              k;
          std_devs[16] =
              ((k - 1) * std_devs[16] +
               (meas.IMU_a_m[3] - y[13]) * (meas.IMU_a_m[3] - y[13])) /
              k;
          std_devs[17] =
              ((k - 1) * std_devs[17] +
               (meas.IMU_a_m[4] - y[14]) * (meas.IMU_a_m[4] - y[14])) /
              k;
          std_devs[18] =
              ((k - 1) * std_devs[18] +
               (meas.IMU_a_m[5] - y[15]) * (meas.IMU_a_m[5] - y[15])) /
              k;

          std_devs[15] =
              ((k - 1) * std_devs[15] +
               ((meas.IMU_a_m[0] - y[10]) * (meas.IMU_a_m[0] - y[10]) +
                (meas.IMU_a_m[1] - y[11]) * (meas.IMU_a_m[1] - y[11]) +
                (meas.IMU_a_m[2] - y[12]) * (meas.IMU_a_m[2] - y[12]))) /
              k;
          std_devs[19] =
              ((k - 1) * std_devs[19] +
               ((meas.IMU_a_m[3] - y[13]) * (meas.IMU_a_m[3] - y[13]) +
                (meas.IMU_a_m[4] - y[14]) * (meas.IMU_a_m[4] - y[14]) +
                (meas.IMU_a_m[5] - y[15]) * (meas.IMU_a_m[5] - y[15]))) /
              k;
        }
      }
    }
    measurements.emplace_back(t, meas);
  }

  if (stat_results) {
    for (double std_dev : std_devs) {
      (*stat_results) << std::sqrt(std_dev);
    }
    (*stat_results) << recorder::data_recorder::end_value_row
                    << recorder::data_recorder::flush;
  }
}

template <typename Sat3DSystemType>
void do_all_single_runs(
    const std::string& output_stem_name, const std::string& filter_name,
    const std::vector<std::pair<double, sat3D_measurement_point>>& measurements,
    const std::vector<std::pair<double, sat3D_state_type>>& ground_truth,
    Sat3DSystemType& sat_sys, const sat3D_state_space_type& state_space,
    const sat3D_state_belief_type& b, sat3D_input_belief_type b_u,
    const sat3D_output_belief_type& b_z, double time_step,
    unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;

  cov_matrix_type Qu = b_u.get_covariance().get_matrix();

  for (unsigned int skips = min_skips; skips <= max_skips; ++skips) {

    sat_sys.set_time_step(skips * time_step);

    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));

    std::vector<std::pair<double, vect_n<double>>> result_pts =
        batch_KF_on_timeseries(measurements, ground_truth, sat_sys, state_space,
                               b, b_u, b_z, skips);

    std::stringstream ss;
    ss << output_stem_name << "_" << std::setfill('0') << std::setw(4)
       << int(1000 * skips * time_step) << "_" << filter_name << ".ssv";
    recorder::ascii_recorder results(ss.str());
    results << "time"
            << "pos_x"
            << "pos_y"
            << "pos_z"
            << "q0"
            << "q1"
            << "q2"
            << "q3"
            << "vel_x"
            << "vel_y"
            << "vel_z"
            << "avel_x"
            << "avel_y"
            << "avel_z"
            << "ep_x"
            << "ep_y"
            << "ep_z"
            << "ea_x"
            << "ea_y"
            << "ea_z"
            << "ev_x"
            << "ev_y"
            << "ev_z"
            << "ew_x"
            << "ew_y"
            << "ew_z"
            << "P_xx"
            << "P_yy"
            << "P_zz"
            << "P_aax"
            << "P_aay"
            << "P_aaz"
            << "P_vvx"
            << "P_vvy"
            << "P_vvz"
            << "P_wwx"
            << "P_wwy"
            << "P_wwz" << recorder::data_recorder::end_name_row;
    for (auto& result_pt : result_pts) {
      results << result_pt.first;
      for (double j : result_pt.second) {
        results << j;
      }
      results << recorder::data_recorder::end_value_row;
    }
    results << recorder::data_recorder::flush;
  }

  sat_sys.set_time_step(time_step);
}

template <typename Sat3DSystemType>
void do_single_monte_carlo_run(
    std::map<std::string, std::shared_ptr<ReaK::recorder::ascii_recorder>>&
        results_map,
    const std::string& output_stem_name, const std::string& filter_name,
    const std::vector<std::pair<double, sat3D_measurement_point>>& measurements,
    const std::vector<std::pair<double, sat3D_state_type>>& ground_truth,
    Sat3DSystemType& sat_sys, const sat3D_state_space_type& state_space,
    const sat3D_state_belief_type& b, sat3D_input_belief_type b_u,
    const sat3D_output_belief_type& b_z, double time_step,
    unsigned int min_skips, unsigned int max_skips) {
  using namespace ReaK;

  cov_matrix_type Qu = b_u.get_covariance().get_matrix();

  for (unsigned int skips = min_skips; skips <= max_skips; ++skips) {

    sat_sys.set_time_step(skips * time_step);

    b_u.set_covariance(cov_type(cov_matrix_type((1.0 / double(skips)) * Qu)));

    std::vector<std::pair<double, vect_n<double>>> result_pts =
        batch_KF_on_timeseries(measurements, ground_truth, sat_sys, state_space,
                               b, b_u, b_z, skips);

    std::vector<double> std_devs(28, 0.0);
    for (std::size_t i = 0; i < result_pts.size(); ++i) {

      // state vector component errors:
      for (std::size_t j = 13; j < 25; ++j) {
        std_devs[j - 13] = (i * std_devs[j - 13] +
                            result_pts[i].second[j] * result_pts[i].second[j]) /
                           (i + 1);
      }

      // position distance error:
      std_devs[12] = (i * std_devs[12] +
                      (result_pts[i].second[13] * result_pts[i].second[13] +
                       result_pts[i].second[14] * result_pts[i].second[14] +
                       result_pts[i].second[15] * result_pts[i].second[15])) /
                     (i + 1);

      // rotation angle error:
      std_devs[13] = (i * std_devs[13] +
                      (result_pts[i].second[16] * result_pts[i].second[16] +
                       result_pts[i].second[17] * result_pts[i].second[17] +
                       result_pts[i].second[18] * result_pts[i].second[18])) /
                     (i + 1);

      // speed error:
      std_devs[14] = (i * std_devs[14] +
                      (result_pts[i].second[19] * result_pts[i].second[19] +
                       result_pts[i].second[20] * result_pts[i].second[20] +
                       result_pts[i].second[21] * result_pts[i].second[21])) /
                     (i + 1);

      // angular speed error:
      std_devs[15] = (i * std_devs[15] +
                      (result_pts[i].second[22] * result_pts[i].second[22] +
                       result_pts[i].second[23] * result_pts[i].second[23] +
                       result_pts[i].second[24] * result_pts[i].second[24])) /
                     (i + 1);

      // average estimated covariances:
      for (std::size_t j = 25; j < 37; ++j) {
        std_devs[j - 25 + 16] =
            (i * std_devs[j - 25 + 16] + result_pts[i].second[j]) / (i + 1);
      }
    }

    std::stringstream ss;
    ss << "_" << std::setfill('0') << std::setw(4)
       << int(1000 * skips * time_step) << "_" << filter_name;
    std::string file_middle = ss.str();
    std::shared_ptr<recorder::ascii_recorder>& results =
        results_map[file_middle];
    if (!results) {
      results = std::make_shared<recorder::ascii_recorder>(
          output_stem_name + file_middle + "_stddevs.ssv");
      (*results) << "ep_x"
                 << "ep_y"
                 << "ep_z"
                 << "ea_x"
                 << "ea_y"
                 << "ea_z"
                 << "ev_x"
                 << "ev_y"
                 << "ev_z"
                 << "ew_x"
                 << "ew_y"
                 << "ew_z"
                 << "ep_m"
                 << "ea_m"
                 << "ev_m"
                 << "ew_m"
                 << "P_xx"
                 << "P_yy"
                 << "P_zz"
                 << "P_aax"
                 << "P_aay"
                 << "P_aaz"
                 << "P_vvx"
                 << "P_vvy"
                 << "P_vvz"
                 << "P_wwx"
                 << "P_wwy"
                 << "P_wwz" << recorder::data_recorder::end_name_row;
    }

    for (std::size_t j = 0; j < 18; ++j) {
      (*results) << std::sqrt(std_devs[j]);  // turn variances into std-devs.
    }
    (*results) << recorder::data_recorder::end_value_row
               << recorder::data_recorder::flush;
  }

  sat_sys.set_time_step(time_step);
}

static std::string strip_quotes(const std::string& s) {
  std::size_t first = 0;
  while ((first < s.length()) &&
         ((std::isspace(s[first]) != 0) || (s[first] == '"'))) {
    ++first;
  }
  std::size_t last = s.length();
  while ((last > 0) &&
         ((std::isspace(s[last - 1]) != 0) || (s[last - 1] == '"'))) {
    --last;
  }
  return s.substr(first, last - first);
}

int main(int argc, char** argv) {
  using namespace ReaK;

  absl::ParseCommandLine(argc, argv);

  std::string output_path_name = strip_quotes(absl::GetFlag(FLAGS_output));
  std::string output_stem_name = output_path_name;
  if (output_stem_name.back() == '/') {
    output_stem_name += "output_record";
  } else {
    std::size_t p = output_path_name.find_last_of('/');
    if (p == std::string::npos) {
      output_path_name = "";
    } else {
      output_path_name.erase(p);
    }
  }
  while (output_path_name.back() == '/') {
    output_path_name.erase(output_path_name.length() - 1, 1);
  }

  if (!output_path_name.empty()) {
    fs::create_directory(output_path_name.c_str());
  }

  std::string sys_output_path_name =
      strip_quotes(absl::GetFlag(FLAGS_system_output));
  std::string sys_output_stem_name = sys_output_path_name;
  if (absl::GetFlag(FLAGS_generate_mdl_files)) {
    if (sys_output_stem_name.back() == '/') {
      sys_output_stem_name += "satellite3D";
    } else {
      std::size_t p = sys_output_path_name.find_last_of('/');
      if (p == std::string::npos) {
        sys_output_path_name = "";
      } else {
        sys_output_path_name.erase(p);
      }
    }
    while (sys_output_path_name.back() == '/') {
      sys_output_path_name.erase(sys_output_path_name.length() - 1, 1);
    }

    if (!sys_output_path_name.empty()) {
      fs::create_directory(sys_output_path_name.c_str());
    }
  }

  double start_time = absl::GetFlag(FLAGS_start_time);
  double end_time = absl::GetFlag(FLAGS_end_time);
  double time_step = absl::GetFlag(FLAGS_time_step);

  unsigned int mc_runs = absl::GetFlag(FLAGS_mc_runs);
  unsigned int min_skips = absl::GetFlag(FLAGS_min_skips);
  unsigned int max_skips = absl::GetFlag(FLAGS_max_skips);

  /* initial states */
  frame_3D<double> initial_motion;
  if (!absl::GetFlag(FLAGS_init).empty()) {
    try {
      std::string init_filename = strip_quotes(absl::GetFlag(FLAGS_init));

      if (absl::GetFlag(FLAGS_generate_mdl_files)) {
        *(serialization::open_oarchive(init_filename)) &
            RK_SERIAL_SAVE_WITH_NAME(initial_motion);
      } else {
        if (!fs::exists(fs::path(init_filename))) {
          std::cout << "Initial-conditions file does not exist!" << std::endl;
          return 3;
        }

        *(serialization::open_iarchive(init_filename)) &
            RK_SERIAL_LOAD_WITH_NAME(initial_motion);
      }
    } catch (std::exception& e) {
      LOG(ERROR) << "An exception occurred during the loading of the initial "
                    "conditions! "
                    "what(): '"
                 << e.what() << "'.";
      return 11;
    }
  }

  /* inertial data */
  double mass = 1.0;
  ReaK::mat<double, ReaK::mat_structure::symmetric> inertia_tensor(
      1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
  try {

    std::string inertia_filename = strip_quotes(absl::GetFlag(FLAGS_inertia));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      *(serialization::open_oarchive(inertia_filename)) &
          RK_SERIAL_SAVE_WITH_NAME(mass) &
          RK_SERIAL_SAVE_WITH_NAME(inertia_tensor);
    } else {
      if (!fs::exists(fs::path(inertia_filename))) {
        std::cout << "Inertial-information file does not exist!" << std::endl;
        return 4;
      }

      *(serialization::open_iarchive(inertia_filename)) &
          RK_SERIAL_LOAD_WITH_NAME(mass) &
          RK_SERIAL_LOAD_WITH_NAME(inertia_tensor);
    }
  } catch (std::exception& e) {
    LOG(ERROR)
        << "An exception occurred during the loading of the inertia matrix! "
           "what(): '"
        << e.what() << "'.";
    return 12;
  }

  /* input disturbance */
  mat<double, mat_structure::diagonal> input_disturbance(6, true);
  try {

    std::string Qu_filename = strip_quotes(absl::GetFlag(FLAGS_Q_matrix));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      *(serialization::open_oarchive(Qu_filename)) &
          RK_SERIAL_SAVE_WITH_NAME(input_disturbance);
    } else {
      if (!fs::exists(fs::path(Qu_filename))) {
        std::cout << "Input disturbance covariance matrix file does not exist!"
                  << std::endl;
        return 5;
      }

      *(serialization::open_iarchive(Qu_filename)) &
          RK_SERIAL_LOAD_WITH_NAME(input_disturbance);
    }
  } catch (std::exception& e) {
    LOG(ERROR)
        << "An exception occurred during the loading of the input disturbance "
           "covariance matrix! what(): '"
        << e.what() << "'.";
    return 13;
  }

  /* measurement noise */
  std::size_t m_noise_size = 6;
  if (absl::GetFlag(FLAGS_gyro)) {
    m_noise_size = 9;
  }
  if (absl::GetFlag(FLAGS_IMU)) {
    m_noise_size = 15;
  }
  mat<double, mat_structure::diagonal> measurement_noise(m_noise_size, true);
  try {
    std::string R_filename = strip_quotes(absl::GetFlag(FLAGS_R_matrix));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      *(serialization::open_oarchive(R_filename)) &
          RK_SERIAL_SAVE_WITH_NAME(measurement_noise);
    } else {
      if (!fs::exists(fs::path(R_filename))) {
        std::cout << "Measurement noise covariance matrix file does not exist!"
                  << std::endl;
        return 6;
      }

      *(serialization::open_iarchive(R_filename)) &
          RK_SERIAL_LOAD_WITH_NAME(measurement_noise);
    }
  } catch (std::exception& e) {
    LOG(ERROR)
        << "An exception occurred during the loading of the measurement noise "
           "covariance matrix! what(): '"
        << e.what() << "'.";
    return 14;
  }

  /* artificial measurement noise */
  mat<double, mat_structure::diagonal> artificial_noise(m_noise_size, 0.0);
  if (!absl::GetFlag(FLAGS_R_added).empty()) {
    try {
      std::string R_added_filename = strip_quotes(absl::GetFlag(FLAGS_R_added));

      if (absl::GetFlag(FLAGS_generate_mdl_files)) {
        *(serialization::open_oarchive(R_added_filename)) &
            RK_SERIAL_SAVE_WITH_NAME(artificial_noise);
      } else {
        if (!fs::exists(fs::path(R_added_filename))) {
          std::cout << "Artificial noise covariance matrix file does not exist!"
                    << std::endl;
          return 6;
        }

        *(serialization::open_iarchive(R_added_filename)) &
            RK_SERIAL_LOAD_WITH_NAME(artificial_noise);
      }
    } catch (std::exception& e) {
      LOG(ERROR)
          << "An exception occurred during the loading of the artificial "
             "measurement noise covariance matrix! what(): '"
          << e.what() << "'.";
      return 3;
    }
  }
  double RAq0 = (artificial_noise(3, 3) + artificial_noise(4, 4) +
                 artificial_noise(5, 5)) /
                12.0;

  /* IMU configuration data */
  unit_quat<double> IMU_orientation;
  vect<double, 3> IMU_location;
  unit_quat<double> earth_orientation;
  vect<double, 3> mag_field_direction(1.0, 0.0, 0.0);
  if (absl::GetFlag(FLAGS_IMU)) {
    try {

      std::string IMUconf_filename =
          strip_quotes(absl::GetFlag(FLAGS_IMU_config));

      if (absl::GetFlag(FLAGS_generate_mdl_files)) {
        *(serialization::open_oarchive(IMUconf_filename)) &
            RK_SERIAL_SAVE_WITH_NAME(IMU_orientation) &
            RK_SERIAL_SAVE_WITH_NAME(IMU_location) &
            RK_SERIAL_SAVE_WITH_NAME(earth_orientation) &
            RK_SERIAL_SAVE_WITH_NAME(mag_field_direction);
      } else {
        if (!fs::exists(fs::path(IMUconf_filename))) {
          std::cout << "IMU configuration data file does not exist!"
                    << std::endl;
          return 6;
        }

        *(serialization::open_iarchive(IMUconf_filename)) &
            RK_SERIAL_LOAD_WITH_NAME(IMU_orientation) &
            RK_SERIAL_LOAD_WITH_NAME(IMU_location) &
            RK_SERIAL_LOAD_WITH_NAME(earth_orientation) &
            RK_SERIAL_LOAD_WITH_NAME(mag_field_direction);
      }
    } catch (std::exception& e) {
      LOG(ERROR) << "An exception occurred during the loading of the IMU "
                    "configuration "
                    "file! what(): '"
                 << e.what() << "'.";
      return 14;
    }
  }

  std::vector<std::pair<double, sat3D_measurement_point>> measurements;
  std::vector<std::pair<double, sat3D_state_type>> ground_truth;
  if ((!absl::GetFlag(FLAGS_monte_carlo)) &&
      !absl::GetFlag(FLAGS_measurements).empty() &&
      fs::exists(fs::path(strip_quotes(absl::GetFlag(FLAGS_measurements))))) {
    try {
      recorder::ascii_extractor measurements_file(
          strip_quotes(absl::GetFlag(FLAGS_measurements)));
      while (true) {
        double t = 0.0;
        measurements_file >> t;
        std::vector<double> meas;
        try {
          while (true) {
            double dummy = 0.0;
            measurements_file >> dummy;
            meas.push_back(dummy);
          }
        } catch (recorder::out_of_bounds&) {}
        measurements_file >> recorder::data_extractor::end_value_row;

        sat3D_measurement_point meas_actual;
        sat3D_measurement_point meas_noisy;

        /* read off the position-orientation measurements. */
        if (meas.size() < 7) {
          LOG(ERROR)
              << "The measurement file does not even contain the position and "
                 "quaternion measurements!";
          return 4;
        }
        meas_actual.pose = vect_n<double>(meas.begin(), meas.begin() + 7);
        meas_noisy.pose = meas_actual.pose;
        if (!absl::GetFlag(FLAGS_R_added).empty()) {
          meas_noisy.pose += vect_n<double>(
              var_rnd() * sqrt(artificial_noise(0, 0)),
              var_rnd() * sqrt(artificial_noise(1, 1)),
              var_rnd() * sqrt(artificial_noise(2, 2)), var_rnd() * sqrt(RAq0),
              var_rnd() * sqrt(0.25 * artificial_noise(3, 3)),
              var_rnd() * sqrt(0.25 * artificial_noise(4, 4)),
              var_rnd() * sqrt(0.25 * artificial_noise(5, 5)));
        }
        meas.erase(meas.begin(), meas.begin() + 7);

        /* read off the IMU/gyro angular velocity measurements. */
        if (absl::GetFlag(FLAGS_gyro) || absl::GetFlag(FLAGS_IMU)) {
          if (meas.size() < 3) {
            LOG(ERROR)
                << "The measurement file does not contain the angular velocity "
                   "measurements!";
            return 4;
          }
          meas_actual.gyro = vect_n<double>(meas.begin(), meas.begin() + 3);
          meas_noisy.gyro = meas_actual.gyro;
          if (!absl::GetFlag(FLAGS_R_added).empty() &&
              (artificial_noise.get_row_count() >= 9)) {
            meas_noisy.gyro +=
                vect_n<double>(var_rnd() * sqrt(artificial_noise(6, 6)),
                               var_rnd() * sqrt(artificial_noise(7, 7)),
                               var_rnd() * sqrt(artificial_noise(8, 8)));
          }
          meas.erase(meas.begin(), meas.begin() + 3);
        }

        /* read off the IMU accel-mag measurements. */
        if (absl::GetFlag(FLAGS_IMU)) {
          if (meas.size() < 6) {
            LOG(ERROR) << "The measurement file does not contain the "
                          "accelerometer and "
                          "magnetometer measurements!";
            return 4;
          }
          meas_actual.IMU_a_m = vect_n<double>(meas.begin(), meas.begin() + 6);
          meas_noisy.IMU_a_m = meas_actual.IMU_a_m;
          if (!absl::GetFlag(FLAGS_R_added).empty() &&
              (artificial_noise.get_row_count() >= 15)) {
            meas_noisy.IMU_a_m +=
                vect_n<double>(var_rnd() * sqrt(artificial_noise(9, 9)),
                               var_rnd() * sqrt(artificial_noise(10, 10)),
                               var_rnd() * sqrt(artificial_noise(11, 11)),
                               var_rnd() * sqrt(artificial_noise(12, 12)),
                               var_rnd() * sqrt(artificial_noise(13, 13)),
                               var_rnd() * sqrt(artificial_noise(14, 14)));
          }
          meas.erase(meas.begin(), meas.begin() + 6);
        }

        /* read off the input vector. */
        if (meas.size() < 6) {
          LOG(ERROR)
              << "The measurement file does not contain the input force-torque "
                 "vector measurements!";
          return 4;
        }
        meas_actual.u = vect_n<double>(meas.begin(), meas.begin() + 6);
        meas_noisy.u = meas_actual.u;
        meas.erase(meas.begin(), meas.begin() + 6);

        /* now, the meas_actual and meas_noisy are fully formed. */
        measurements.emplace_back(t, meas_noisy);

        /* check if the file contains a ground-truth: */
        if (meas.size() >= 7) {
          sat3D_state_type x;
          set_position(x, vect<double, 3>(meas[0], meas[1], meas[2]));
          set_quaternion(x,
                         unit_quat<double>(meas[3], meas[4], meas[5], meas[6]));
          if (meas.size() >= 13) {
            set_velocity(x, vect<double, 3>(meas[7], meas[8], meas[9]));
            set_ang_velocity(x, vect<double, 3>(meas[10], meas[11], meas[12]));
          }
          ground_truth.emplace_back(t, x);
          meas.erase(meas.begin(), meas.end());
        } else if (!absl::GetFlag(FLAGS_R_added).empty()) {
          sat3D_state_type x;
          set_position(x,
                       vect<double, 3>(meas_actual.pose[0], meas_actual.pose[1],
                                       meas_actual.pose[2]));
          set_quaternion(
              x, unit_quat<double>(meas_actual.pose[3], meas_actual.pose[4],
                                   meas_actual.pose[5], meas_actual.pose[6]));
          ground_truth.emplace_back(t, x);
        }
      }
    } catch (recorder::out_of_bounds&) {
      LOG(ERROR) << "The measurement file does not appear to have the required "
                    "number of "
                    "columns!";
      return 4;
    } catch (recorder::end_of_record&) {}
  }

  sat3D_temp_space_type sat_space(
      "satellite3D_temporal_space", sat3D_state_space_type(),
      pp::time_poisson_topology("satellite3D_time_space", time_step,
                                (end_time - start_time) * 0.5));

  using sat3D_traj_type = pp::discrete_point_trajectory<sat3D_temp_space_type>;

  std::shared_ptr<sat3D_traj_type> traj_ptr;
  if (absl::GetFlag(FLAGS_generate_meas_file) &&
      (absl::GetFlag(FLAGS_xml) || absl::GetFlag(FLAGS_protobuf) ||
       absl::GetFlag(FLAGS_binary))) {
    traj_ptr = std::make_shared<sat3D_traj_type>(
        std::shared_ptr<sat3D_temp_space_type>(&sat_space, null_deleter()));
  }

  sat3D_state_type x_init;
  set_position(x_init, vect<double, 3>(0.0, 0.0, 0.0));
  set_velocity(x_init, vect<double, 3>(0.0, 0.0, 0.0));
  set_quaternion(x_init, unit_quat<double>());
  set_ang_velocity(x_init, vect<double, 3>(0.0, 0.0, 0.0));

  sat3D_state_belief_type b_init(
      x_init, cov_type(cov_matrix_type(
                  mat<double, mat_structure::diagonal>(12, 10.0))));

  sat3D_input_belief_type b_u(sat3D_input_type(vect_n<double>(6, 0.0)),
                              cov_type(cov_matrix_type(input_disturbance)));

  if (!absl::GetFlag(FLAGS_gyro) && !absl::GetFlag(FLAGS_IMU)) {

    // Create the set of satellite3D systems for when there is only pose measurements:

    ctrl::satellite3D_inv_dt_system sat3D_inv("satellite3D_inv", mass,
                                              inertia_tensor, time_step);

    ctrl::satellite3D_imdt_sys sat3D_invmom("satellite3D_invmom", mass,
                                            inertia_tensor, time_step);

    ctrl::satellite3D_imdt_sys sat3D_invmid("satellite3D_invmid", mass,
                                            inertia_tensor, time_step, 2);

    sat3D_output_belief_type b_z(
        sat3D_output_type(vect_n<double>(7, 0.0)),
        cov_type(cov_matrix_type(measurement_noise + artificial_noise)));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      try {
        std::shared_ptr<ctrl::satellite3D_inv_dt_system> satellite3D_system =
            std::shared_ptr<ctrl::satellite3D_inv_dt_system>(&sat3D_inv,
                                                             null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name + "_inv_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);

        satellite3D_system = std::shared_ptr<ctrl::satellite3D_inv_dt_system>(
            &sat3D_invmom, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmom_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);

        satellite3D_system = std::shared_ptr<ctrl::satellite3D_inv_dt_system>(
            &sat3D_invmid, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmid_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch (...) {
        LOG(ERROR)
            << "An exception occurred during the saving the satellite system "
               "file!";
        return 14;
      }
    } else if (!absl::GetFlag(FLAGS_monte_carlo)) {

      if (measurements.empty()) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, initial_motion);
        generate_timeseries(
            measurements, ground_truth, sat3D_inv,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise));

        // and output those if asked for it:
        if (absl::GetFlag(FLAGS_generate_meas_file)) {
          recorder::ascii_recorder measurements_gen(output_stem_name +
                                                    "_meas.ssv");
          measurements_gen << "time"
                           << "p_x"
                           << "p_y"
                           << "p_z"
                           << "q_0"
                           << "q_1"
                           << "q_2"
                           << "q_2"
                           << "f_x"
                           << "f_y"
                           << "f_z"
                           << "t_x"
                           << "t_y"
                           << "t_z"
                           << "p_x_true"
                           << "p_y_true"
                           << "p_z_true"
                           << "q_0_true"
                           << "q_1_true"
                           << "q_2_true"
                           << "q_2_true"
                           << "v_x_true"
                           << "v_y_true"
                           << "v_z_true"
                           << "w_x_true"
                           << "w_y_true"
                           << "w_z_true"
                           << recorder::data_recorder::end_name_row;

          for (std::size_t i = 0; i < measurements.size(); ++i) {
            measurements_gen << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            measurements_gen << m.pose << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            measurements_gen << get_position(g) << get_quaternion(g)
                             << get_velocity(g) << get_ang_velocity(g);
            measurements_gen << recorder::data_recorder::end_value_row;

            if (absl::GetFlag(FLAGS_xml) || absl::GetFlag(FLAGS_protobuf) ||
                absl::GetFlag(FLAGS_binary)) {
              traj_ptr->push_back(
                  sat3D_temp_point_type(ground_truth[i].first, g));
            }
          }
          measurements_gen << recorder::data_recorder::flush;
        }
      }

      // do a single run for each skips:

      std::cout << "Running estimators on data series.." << std::flush;

      if (absl::GetFlag(FLAGS_iekf)) {
        do_all_single_runs(output_stem_name, "iekf", measurements, ground_truth,
                           sat3D_inv, sat_space.get_space_topology(), b_init,
                           b_u, b_z, time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkf)) {
        do_all_single_runs(output_stem_name, "imkf", measurements, ground_truth,
                           sat3D_invmom, sat_space.get_space_topology(), b_init,
                           b_u, b_z, time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkfv2)) {
        do_all_single_runs(output_stem_name, "imkfv2", measurements,
                           ground_truth, sat3D_invmid,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      std::cout << "Finished!" << std::endl;

    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, initial_motion);
      std::map<std::string, std::shared_ptr<recorder::ascii_recorder>>
          results_map;
      std::shared_ptr<recorder::ascii_recorder> ground_truth_stats(
          new recorder::ascii_recorder(output_stem_name + "_meas_stddevs.ssv"));
      (*ground_truth_stats) << "ep_x"
                            << "ep_y"
                            << "ep_z"
                            << "ea_x"
                            << "ea_y"
                            << "ea_z"
                            << "ep_m"
                            << "ea_m" << recorder::data_recorder::end_name_row;

      std::cout << "Running Monte-Carlo Simulations..." << std::endl;

      for (unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {

        std::cout << "\r" << std::setw(10) << mc_i << std::flush;

        generate_timeseries(
            measurements, ground_truth, sat3D_inv,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise),
            ground_truth_stats);

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_iekf)) {
          do_single_monte_carlo_run(results_map, output_stem_name, "iekf",
                                    measurements, ground_truth, sat3D_inv,
                                    sat_space.get_space_topology(), b_init, b_u,
                                    b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkf)) {
          do_single_monte_carlo_run(results_map, output_stem_name, "imkf",
                                    measurements, ground_truth, sat3D_invmom,
                                    sat_space.get_space_topology(), b_init, b_u,
                                    b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkfv2)) {
          do_single_monte_carlo_run(results_map, output_stem_name, "imkfv2",
                                    measurements, ground_truth, sat3D_invmid,
                                    sat_space.get_space_topology(), b_init, b_u,
                                    b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;
      }

      std::cout << "Finished!" << std::endl;
    }

  } else if (absl::GetFlag(FLAGS_gyro) && !absl::GetFlag(FLAGS_IMU)) {

    // Create the set of satellite3D systems for when there is gyro measurements:

    ctrl::satellite3D_gyro_inv_dt_system sat3D_inv_gyro(
        "satellite3D_inv_with_gyros", mass, inertia_tensor, time_step);

    ctrl::satellite3D_gyro_imdt_sys sat3D_invmom_gyro(
        "satellite3D_invmom_with_gyros", mass, inertia_tensor, time_step);

    ctrl::satellite3D_gyro_imdt_sys sat3D_invmid_gyro(
        "satellite3D_invmid_with_gyros", mass, inertia_tensor, time_step, 2);

    sat3D_output_belief_type b_z(
        sat3D_output_type(vect_n<double>(10, 0.0)),
        cov_type(cov_matrix_type(measurement_noise + artificial_noise)));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      try {
        std::shared_ptr<ctrl::satellite3D_inv_dt_system> satellite3D_system =
            std::shared_ptr<ctrl::satellite3D_inv_dt_system>(&sat3D_inv_gyro,
                                                             null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_inv_gyro_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);

        satellite3D_system = std::shared_ptr<ctrl::satellite3D_inv_dt_system>(
            &sat3D_invmom_gyro, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmom_gyro_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);

        satellite3D_system = std::shared_ptr<ctrl::satellite3D_inv_dt_system>(
            &sat3D_invmid_gyro, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmid_gyro_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch (...) {
        LOG(ERROR)
            << "An exception occurred during the saving the satellite system "
               "file!";
        return 14;
      }
    } else if (!absl::GetFlag(FLAGS_monte_carlo)) {

      if (measurements.empty()) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, initial_motion);
        generate_timeseries(
            measurements, ground_truth, sat3D_inv_gyro,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise));

        // and output those if asked for it:
        if (absl::GetFlag(FLAGS_generate_meas_file)) {
          recorder::ascii_recorder measurements_gen(output_stem_name +
                                                    "_meas.ssv");
          measurements_gen << "time"
                           << "p_x"
                           << "p_y"
                           << "p_z"
                           << "q_0"
                           << "q_1"
                           << "q_2"
                           << "q_2"
                           << "w_x"
                           << "w_y"
                           << "w_z"
                           << "f_x"
                           << "f_y"
                           << "f_z"
                           << "t_x"
                           << "t_y"
                           << "t_z"
                           << "p_x_true"
                           << "p_y_true"
                           << "p_z_true"
                           << "q_0_true"
                           << "q_1_true"
                           << "q_2_true"
                           << "q_2_true"
                           << "v_x_true"
                           << "v_y_true"
                           << "v_z_true"
                           << "w_x_true"
                           << "w_y_true"
                           << "w_z_true"
                           << recorder::data_recorder::end_name_row;

          for (std::size_t i = 0; i < measurements.size(); ++i) {
            measurements_gen << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            measurements_gen << m.pose << m.gyro << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            measurements_gen << get_position(g) << get_quaternion(g)
                             << get_velocity(g) << get_ang_velocity(g);
            measurements_gen << recorder::data_recorder::end_value_row;

            if (absl::GetFlag(FLAGS_xml) || absl::GetFlag(FLAGS_protobuf) ||
                absl::GetFlag(FLAGS_binary)) {
              traj_ptr->push_back(
                  sat3D_temp_point_type(ground_truth[i].first, g));
            }
          }
          measurements_gen << recorder::data_recorder::flush;
        }
      }

      // do a single run for each skips:

      std::cout << "Running estimators on data series.." << std::flush;

      if (absl::GetFlag(FLAGS_iekf)) {
        do_all_single_runs(output_stem_name, "iekf_gyro", measurements,
                           ground_truth, sat3D_inv_gyro,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkf)) {
        do_all_single_runs(output_stem_name, "imkf_gyro", measurements,
                           ground_truth, sat3D_invmom_gyro,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkfv2)) {
        do_all_single_runs(output_stem_name, "imkfv2_gyro", measurements,
                           ground_truth, sat3D_invmid_gyro,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      std::cout << "Finished!" << std::endl;

    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, initial_motion);
      std::map<std::string, std::shared_ptr<recorder::ascii_recorder>>
          results_map;
      std::shared_ptr<recorder::ascii_recorder> ground_truth_stats(
          new recorder::ascii_recorder(output_stem_name +
                                       "_meas_gyro_stddevs.ssv"));
      (*ground_truth_stats) << "ep_x"
                            << "ep_y"
                            << "ep_z"
                            << "ea_x"
                            << "ea_y"
                            << "ea_z"
                            << "ep_m"
                            << "ea_m"
                            << "ew_x"
                            << "ew_y"
                            << "ew_z"
                            << "ew_m" << recorder::data_recorder::end_name_row;

      std::cout << "Running Monte-Carlo Simulations..." << std::endl;

      for (unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {

        std::cout << "\r" << std::setw(10) << mc_i << std::flush;

        generate_timeseries(
            measurements, ground_truth, sat3D_inv_gyro,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise),
            ground_truth_stats);

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_iekf)) {
          do_single_monte_carlo_run(results_map, output_stem_name, "iekf_gyro",
                                    measurements, ground_truth, sat3D_inv_gyro,
                                    sat_space.get_space_topology(), b_init, b_u,
                                    b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkf)) {
          do_single_monte_carlo_run(
              results_map, output_stem_name, "imkf_gyro", measurements,
              ground_truth, sat3D_invmom_gyro, sat_space.get_space_topology(),
              b_init, b_u, b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkfv2)) {
          do_single_monte_carlo_run(
              results_map, output_stem_name, "imkfv2_gyro", measurements,
              ground_truth, sat3D_invmid_gyro, sat_space.get_space_topology(),
              b_init, b_u, b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;
      }

      std::cout << "Finished!" << std::endl;
    }

  } else {

    // Create the set of satellite3D systems for when there is IMU measurements:

    ctrl::satellite3D_IMU_imdt_sys sat3D_invmom_IMU(
        "satellite3D_invmom_with_IMU", mass, inertia_tensor, time_step,
        IMU_orientation, IMU_location, earth_orientation, mag_field_direction);

    ctrl::satellite3D_IMU_imdt_sys sat3D_invmid_IMU(
        "satellite3D_invmid_with_IMU", mass, inertia_tensor, time_step,
        IMU_orientation, IMU_location, earth_orientation, mag_field_direction,
        2);

    sat3D_output_belief_type b_z(
        sat3D_output_type(vect_n<double>(16, 0.0)),
        cov_type(cov_matrix_type(measurement_noise + artificial_noise)));

    if (absl::GetFlag(FLAGS_generate_mdl_files)) {
      try {
        std::shared_ptr<ctrl::satellite3D_inv_dt_system> satellite3D_system =
            std::shared_ptr<ctrl::satellite3D_inv_dt_system>(&sat3D_invmom_IMU,
                                                             null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmom_IMU_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);

        satellite3D_system = std::shared_ptr<ctrl::satellite3D_inv_dt_system>(
            &sat3D_invmid_IMU, null_deleter());
        *(serialization::open_oarchive(sys_output_stem_name +
                                       "_invmid_IMU_mdl.rkx")) &
            RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      } catch (...) {
        LOG(ERROR)
            << "An exception occurred during the saving the satellite system "
               "file!";
        return 14;
      }
    } else if (!absl::GetFlag(FLAGS_monte_carlo)) {

      if (measurements.empty()) {
        // must generate the measurements and ground_truth vectors:
        set_frame_3D(x_init, initial_motion);
        generate_timeseries(
            measurements, ground_truth, sat3D_invmom_IMU,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise));

        // and output those if asked for it:
        if (absl::GetFlag(FLAGS_generate_meas_file)) {
          recorder::ascii_recorder measurements_gen(output_stem_name +
                                                    "_meas.ssv");
          measurements_gen << "time"
                           << "p_x"
                           << "p_y"
                           << "p_z"
                           << "q_0"
                           << "q_1"
                           << "q_2"
                           << "q_2"
                           << "w_x"
                           << "w_y"
                           << "w_z"
                           << "acc_x"
                           << "acc_y"
                           << "acc_z"
                           << "mag_x"
                           << "mag_y"
                           << "mag_z"
                           << "f_x"
                           << "f_y"
                           << "f_z"
                           << "t_x"
                           << "t_y"
                           << "t_z"
                           << "p_x_true"
                           << "p_y_true"
                           << "p_z_true"
                           << "q_0_true"
                           << "q_1_true"
                           << "q_2_true"
                           << "q_2_true"
                           << "v_x_true"
                           << "v_y_true"
                           << "v_z_true"
                           << "w_x_true"
                           << "w_y_true"
                           << "w_z_true"
                           << recorder::data_recorder::end_name_row;

          for (std::size_t i = 0; i < measurements.size(); ++i) {
            measurements_gen << measurements[i].first;
            const sat3D_measurement_point& m = measurements[i].second;
            measurements_gen << m.pose << m.gyro << m.IMU_a_m << m.u;
            const sat3D_state_type& g = ground_truth[i].second;
            measurements_gen << get_position(g) << get_quaternion(g)
                             << get_velocity(g) << get_ang_velocity(g);
            measurements_gen << recorder::data_recorder::end_value_row;

            if (absl::GetFlag(FLAGS_xml) || absl::GetFlag(FLAGS_protobuf) ||
                absl::GetFlag(FLAGS_binary)) {
              traj_ptr->push_back(
                  sat3D_temp_point_type(ground_truth[i].first, g));
            }
          }
          measurements_gen << recorder::data_recorder::flush;
        }
      }

      // do a single run for each skips:

      std::cout << "Running estimators on data series.." << std::flush;

      if (absl::GetFlag(FLAGS_iekf)) {
        LOG(ERROR) << "Warning: The invariant extended Kalman filter (IEKF) is "
                      "not available for full IMU measurements!"
                   << std::endl;
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkf)) {
        do_all_single_runs(output_stem_name, "imkf_IMU", measurements,
                           ground_truth, sat3D_invmom_IMU,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      if (absl::GetFlag(FLAGS_imkfv2)) {
        do_all_single_runs(output_stem_name, "imkfv2_IMU", measurements,
                           ground_truth, sat3D_invmid_IMU,
                           sat_space.get_space_topology(), b_init, b_u, b_z,
                           time_step, min_skips, max_skips);
      }

      std::cout << "." << std::flush;

      std::cout << "Finished!" << std::endl;

    } else {
      // do monte-carlo runs:
      set_frame_3D(x_init, initial_motion);
      std::map<std::string, std::shared_ptr<recorder::ascii_recorder>>
          results_map;
      std::shared_ptr<recorder::ascii_recorder> ground_truth_stats(
          new recorder::ascii_recorder(output_stem_name +
                                       "_meas_IMU_stddevs.ssv"));
      (*ground_truth_stats)
          << "ep_x"
          << "ep_y"
          << "ep_z"
          << "ea_x"
          << "ea_y"
          << "ea_z"
          << "ep_m"
          << "ea_m"
          << "ew_x"
          << "ew_y"
          << "ew_z"
          << "ew_m"
          << "eacc_x"
          << "eacc_y"
          << "eacc_z"
          << "eacc_m"
          << "emag_y"
          << "emag_y"
          << "emag_z"
          << "emag_m" << recorder::data_recorder::end_name_row;

      if (absl::GetFlag(FLAGS_iekf)) {
        LOG(ERROR) << "Warning: The invariant extended Kalman filter (IEKF) is "
                      "not available for full IMU measurements!"
                   << std::endl;
      }

      std::cout << "Running Monte-Carlo Simulations..." << std::endl;

      for (unsigned int mc_i = 0; mc_i < mc_runs; ++mc_i) {

        std::cout << "\r" << std::setw(10) << mc_i << std::flush;

        generate_timeseries(
            measurements, ground_truth, sat3D_invmom_IMU,
            sat_space.get_space_topology(), x_init, start_time, end_time,
            cov_matrix_type(input_disturbance),
            cov_matrix_type(measurement_noise + artificial_noise),
            ground_truth_stats);

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkf)) {
          do_single_monte_carlo_run(
              results_map, output_stem_name, "imkf_IMU", measurements,
              ground_truth, sat3D_invmom_IMU, sat_space.get_space_topology(),
              b_init, b_u, b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;

        if (absl::GetFlag(FLAGS_imkfv2)) {
          do_single_monte_carlo_run(
              results_map, output_stem_name, "imkfv2_IMU", measurements,
              ground_truth, sat3D_invmid_IMU, sat_space.get_space_topology(),
              b_init, b_u, b_z, time_step, min_skips, max_skips);
        }

        std::cout << "." << std::flush;
      }

      std::cout << "Finished!" << std::endl;
    }
  }

  if (absl::GetFlag(FLAGS_generate_meas_file) &&
      (absl::GetFlag(FLAGS_xml) || absl::GetFlag(FLAGS_protobuf) ||
       absl::GetFlag(FLAGS_binary))) {
    std::cout << "Saving the generated trajectory.." << std::flush;

    std::cout << "." << std::flush;

    if (absl::GetFlag(FLAGS_xml)) {

      *(serialization::open_oarchive(output_stem_name + "_traj.rkx")) &
          RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    }

    std::cout << "." << std::flush;

    if (absl::GetFlag(FLAGS_protobuf)) {

      *(serialization::open_oarchive(output_stem_name + "_traj.pb")) &
          RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    }

    std::cout << "." << std::flush;

    if (absl::GetFlag(FLAGS_binary)) {

      *(serialization::open_oarchive(output_stem_name + "_traj.rkb")) &
          RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    }

    std::cout << "Finished!" << std::endl;
  }
}
