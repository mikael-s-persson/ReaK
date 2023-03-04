
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

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include <ReaK/topologies/spaces/differentiable_space.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>

// #define RK_ENABLE_TEST_LINEAR_INTERPOLATOR
// #define RK_ENABLE_TEST_CUBIC_INTERPOLATOR
// #define RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#define RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#define RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR

#include <ReaK/topologies/spaces/Ndof_spaces.hpp>
#include <ReaK/topologies/spaces/hyperbox_topology.hpp>

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#include <ReaK/topologies/interpolation/linear_interp.hpp>
#endif

#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#include <ReaK/topologies/interpolation/cubic_hermite_interp.hpp>
#endif

#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#include <ReaK/topologies/interpolation/quintic_hermite_interp.hpp>
#endif

#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#include <ReaK/topologies/interpolation/sustained_velocity_pulse_Ndof.hpp>
#include <ReaK/topologies/interpolation/svp_Ndof_reach_topologies.hpp>
#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#include <ReaK/topologies/interpolation/sap_Ndof_reach_topologies.hpp>
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof.hpp>
#endif

#include <ReaK/core/base/scope_guard.hpp>

#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

// I/O options
ABSL_FLAG(std::string, output_path, "test_interp_results",
          "Specify the output path (default is test_interp_results).");

// Monte-Carlo options
ABSL_FLAG(int, mc_runs, 100,
          "Number of monte-carlo runs to perform (default is 100).");

// Monte-Carlo options
ABSL_FLAG(int, space_dimensionality, 1,
          "Number of dimensions for the underlying space (default is 1).");
ABSL_FLAG(
    double, space_max_frequency, 10.0,
    "The maximum frequency of the sinusoidal curves (default is 10.0 Hz).");
ABSL_FLAG(double, interp_steps, 0.05,
          "The time-step between the interpolator's control-points, over a "
          "total curve-time of 1.0 second (default is 0.05 seconds)");

// Interpolator selection options
ABSL_FLAG(bool, all_interpolators, false,
          "Specify that all supported interpolators should be run (default if "
          "no particular interpolator is specified).");
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
ABSL_FLAG(bool, linear, false,
          "Specify that the linear interpolation should be run.");
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
ABSL_FLAG(bool, cubic, false,
          "Specify that the cubic interpolation should be run.");
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
ABSL_FLAG(bool, quintic, false,
          "Specify that the quintic interpolation should be run.");
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
ABSL_FLAG(
    bool, svp_Ndof, false,
    "Specify that the sustained-velocity-pulse interpolation should be run.");
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
ABSL_FLAG(bool, sap_Ndof, false,
          "Specify that the sustained-acceleration-pulse interpolation should "
          "be run.");
#endif

namespace fs = std::filesystem;

#define RK_TEST_IS_NAN(X) std::isnan(X[0])
#define RK_TEST_IS_INF(X) std::isinf(X[0])
#define RK_TEST_GET_VALUE(X) (X[0])

template <typename Vector>
bool vect_is_nan(const Vector& v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (std::isnan(v[i])) {
      return true;
    }
  }
  return false;
};

template <typename Vector>
bool vect_is_inf(const Vector& v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (std::isinf(v[i])) {
      return true;
    }
  }
  return false;
};

template <typename InterpTopoType, typename PtContainer>
void try_interpolation(const std::string& aMethodName,
                       const std::string& aTestData, std::size_t& succ_count,
                       std::size_t& graceful_failures,
                       const InterpTopoType& topo, const PtContainer& pts,
                       std::ostream& fail_reports) {

  using namespace ReaK;

  using PointType = typename pp::topology_traits<InterpTopoType>::point_type;
  using Iter = typename PtContainer::const_iterator;

  auto it = pts.begin();
  for (auto prev_it = it++; it != pts.end(); ++prev_it, ++it) {

    try {

      RK_SCOPE_EXIT_ROUTINE(report_dist_except) {
        fail_reports << aMethodName << " exception dist " << aTestData
                     << std::endl;
      };
      double d = get(pp::distance_metric, topo)(*prev_it, *it, topo);
      RK_SCOPE_EXIT_DISMISS(report_dist_except);

      if (std::isnan(d)) {
        fail_reports << aMethodName << " NaN dist " << aTestData << std::endl;
        throw std::domain_error("NaN condition encountered!");
      };

      if (std::isinf(d)) {
        ++graceful_failures;
        continue;
      };

      RK_SCOPE_EXIT_ROUTINE(report_interp_except) {
        fail_reports << aMethodName << " exception interp " << aTestData
                     << std::endl;
      };
      PointType p = topo.move_position_toward(*prev_it, 1.0, *it);
      RK_SCOPE_EXIT_DISMISS(report_interp_except);

      if (vect_is_nan((get<0>(p))) || vect_is_nan((get<1>(p))) ||
          vect_is_nan((get<2>(p)))) {
        fail_reports << aMethodName << " NaN interp " << aTestData << std::endl;
        throw std::domain_error("NaN condition encountered!");
      };

      if (vect_is_inf((get<0>(p))) || vect_is_inf((get<1>(p))) ||
          vect_is_inf((get<2>(p)))) {
        fail_reports << aMethodName << " INF interp " << aTestData << std::endl;
        throw std::domain_error("INF condition encountered!");
      };

      if (norm_2(get<0>(p) - get<0>(*it)) > 1e-3) {
        fail_reports << aMethodName << " pos_tol interp " << aTestData
                     << std::endl;
        throw std::domain_error("Position-tolerance exceeded!");
      };

      if (norm_2(get<1>(p) - get<1>(*it)) > 1e-3) {
        fail_reports << aMethodName << " vel_tol interp " << aTestData
                     << std::endl;
        throw std::domain_error("Velocity-tolerance exceeded!");
      };

      ++succ_count;

    } catch (std::exception& e) {
      RK_UNUSED(e);
    };
  };
};

template <std::size_t StaticSpDim>
struct interp_mc_test_space {
  using topo_type =
      typename ReaK::pp::Ndof_rl_space<double, StaticSpDim, 2>::type;
  using vector_type = ReaK::vect<double, StaticSpDim>;

  static vector_type default_vect(std::size_t /*unused*/) {
    return vector_type();
  };

  static topo_type create(const vector_type& lb, const vector_type& ub,
                          const vector_type& sb, const vector_type& ab,
                          const vector_type& jb) {
    return ReaK::pp::make_Ndof_rl_space<StaticSpDim>(lb, ub, sb, ab, jb);
  };
};

template <>
struct interp_mc_test_space<0> {
  using topo_type = ReaK::pp::Ndof_rl_space<double, 0, 2>::type;
  using vector_type = ReaK::vect_n<double>;

  static vector_type default_vect(std::size_t dyn_sp_size) {
    return vector_type(dyn_sp_size);
  };

  static topo_type create(const vector_type& lb, const vector_type& ub,
                          const vector_type& sb, const vector_type& ab,
                          const vector_type& jb) {
    return ReaK::pp::make_Ndof_rl_space(lb, ub, sb, ab, jb);
  };
};

template <std::size_t StaticSpDim>
void perform_mc_tests(std::size_t dyn_sp_dim) {

  using namespace ReaK;

  using Config = interp_mc_test_space<StaticSpDim>;
  using TopoType = typename Config::topo_type;
  using Vector = typename Config::vector_type;

  using PointType = typename pp::topology_traits<TopoType>::point_type;

  double max_freq = absl::GetFlag(FLAGS_space_max_frequency);
  double max_rad_freq = max_freq * 2.0 * M_PI;  // rad/s

  Vector lb = Config::default_vect(dyn_sp_dim);
  Vector ub = Config::default_vect(dyn_sp_dim);
  Vector sb = Config::default_vect(dyn_sp_dim);
  Vector ab = Config::default_vect(dyn_sp_dim);
  Vector jb = Config::default_vect(dyn_sp_dim);
  for (std::size_t i = 0; i < lb.size(); ++i) {
    lb[i] = -2.0;
    ub[i] = 2.0;
    sb[i] = 2.0 * max_rad_freq;
    ab[i] = 2.0 * max_rad_freq * max_rad_freq;
    jb[i] = 2.0 * max_rad_freq * max_rad_freq * max_rad_freq;
  };

  std::string output_path = absl::GetFlag(FLAGS_output_path);
  while (output_path[output_path.length() - 1] == '/') {
    output_path.erase(output_path.length() - 1, 1);
  }

  fs::create_directory(output_path.c_str());

  std::ofstream fail_reports((output_path + "/topo_fail_reports.txt").c_str());

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  std::size_t linear_succ_count = 0;
  std::size_t linear_graceful_fails = 0;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  std::size_t cubic_succ_count = 0;
  std::size_t cubic_graceful_fails = 0;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  std::size_t quintic_succ_count = 0;
  std::size_t quintic_graceful_fails = 0;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  std::size_t svp_Ndof_succ_count = 0;
  std::size_t svp_Ndof_graceful_fails = 0;
  pp::interpolated_topology<TopoType, pp::svp_Ndof_interpolation_tag>
      svp_Ndof_topo(Config::create(lb, ub, sb, ab, jb));
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  std::size_t sap_Ndof_succ_count = 0;
  std::size_t sap_Ndof_graceful_fails = 0;
  pp::interpolated_topology<TopoType, pp::sap_Ndof_interpolation_tag>
      sap_Ndof_topo(Config::create(lb, ub, sb, ab, jb));
#endif

  std::size_t mc_runs = absl::GetFlag(FLAGS_mc_runs);
  double interp_steps = absl::GetFlag(FLAGS_interp_steps);

  std::size_t segment_count = 0;
  for (double t = 0.0; t < 1.0 + 0.5 * interp_steps; t += interp_steps) {
    ++segment_count;
  }
  --segment_count;

  global_rng_type& gbl_rng = get_global_rng();

  for (std::size_t i = 0; i < mc_runs; ++i) {

    std::vector<PointType> pts;
    double curve_freq = double(gbl_rng() % 1000) * (max_rad_freq / 1000.0);
    Vector curve_ampl = Config::default_vect(dyn_sp_dim);
    Vector curve_phase = Config::default_vect(dyn_sp_dim);
    for (std::size_t j = 0; j < curve_ampl.size(); ++j) {
      curve_ampl[j] = double(gbl_rng() % 1000) * 0.001;
      curve_phase[j] = double(gbl_rng() % 1000) * (M_PI / 500.0);
    };

    for (double t = 0.0; t < 1.0 + 0.5 * interp_steps; t += interp_steps) {
      Vector pos = curve_ampl;
      Vector vel = curve_ampl;
      Vector acc = curve_ampl;

      for (std::size_t j = 0; j < pos.size(); ++j) {
        pos[j] *=
            std::sin(curve_freq * t + curve_phase[j]) / (2.0 * max_rad_freq);
        vel[j] *= curve_freq * std::cos(curve_freq * t + curve_phase[j]) /
                  (2.0 * max_rad_freq * max_rad_freq);
        acc[j] *= -curve_freq * curve_freq *
                  std::sin(curve_freq * t + curve_phase[j]) /
                  (2.0 * max_rad_freq * max_rad_freq * max_rad_freq);
      };

      pts.push_back(PointType(pos, vel, acc));
    };

    curve_ampl *= 0.5 / max_rad_freq;

    std::string test_data_str;
    {
      std::stringstream ss;
      ss << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps;
      test_data_str = ss.str();
    };

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_linear)) {
      try_interpolation("linear", test_data_str, linear_succ_count,
                        linear_graceful_fails, linear_topo, pts, fail_reports);
    };

#endif

#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_cubic)) {
      try_interpolation("cubic", test_data_str, cubic_succ_count,
                        cubic_graceful_fails, cubic_topo, pts, fail_reports);
    };

#endif

#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_quintic)) {
      try_interpolation("quintic", test_data_str, quintic_succ_count,
                        quintic_graceful_fails, quintic_topo, pts,
                        fail_reports);
    };

#endif

#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_svp_Ndof)) {
      try_interpolation("svp_Ndof", test_data_str, svp_Ndof_succ_count,
                        svp_Ndof_graceful_fails, svp_Ndof_topo, pts,
                        fail_reports);
    };

#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_sap_Ndof)) {
      try_interpolation("sap_Ndof", test_data_str, sap_Ndof_succ_count,
                        sap_Ndof_graceful_fails, sap_Ndof_topo, pts,
                        fail_reports);
    };

#endif
  };

  fail_reports.close();

  std::size_t total_segments = segment_count * mc_runs;
  std::ofstream succ_reports((output_path + "/topo_success_rates.txt").c_str());
  succ_reports << "Monte-Carlo runs of interpolation methods using random "
                  "sinusoidal curves.\n"
               << "  Number of runs: " << mc_runs << "\n"
               << "  Number of segments per run: " << segment_count << "\n"
               << "  Spatial dimensions: " << dyn_sp_dim << "\n"
               << "  Maximum frequency: " << max_freq << " Hz ( "
               << max_rad_freq << " rad/s )\n"
               << "  Interpolation steps of: " << interp_steps << std::endl;

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_linear))
    succ_reports << "Linear interp"
                 << "\n"
                 << "\t Successes: " << linear_succ_count << "\n"
                 << "\t Graceful Failures: " << linear_graceful_fails << "\n"
                 << "\t Total num. of trials: " << total_segments << "\n"
                 << "\t Success-rate: "
                 << (100.0 * double(linear_succ_count) / double(total_segments))
                 << "\%\n"
                 << "\t Graceful-fail-rate: "
                 << (100.0 * double(linear_graceful_fails) /
                     double(total_segments))
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << (100.0 *
                     (1.0 - double(linear_succ_count + linear_graceful_fails) /
                                double(total_segments)))
                 << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_cubic))
    succ_reports
        << "Cubic interp"
        << "\n"
        << "\t Successes: " << cubic_succ_count << "\n"
        << "\t Graceful Failures: " << cubic_graceful_fails << "\n"
        << "\t Total num. of trials: " << total_segments << "\n"
        << "\t Success-rate: "
        << (100.0 * double(cubic_succ_count) / double(total_segments)) << "\%\n"
        << "\t Graceful-fail-rate: "
        << (100.0 * double(cubic_graceful_fails) / double(total_segments))
        << "\%\n"
        << "\t Hard-fail-rate: "
        << (100.0 * (1.0 - double(cubic_succ_count + cubic_graceful_fails) /
                               double(total_segments)))
        << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_quintic))
    succ_reports
        << "Quintic interp"
        << "\n"
        << "\t Successes: " << quintic_succ_count << "\n"
        << "\t Graceful Failures: " << quintic_graceful_fails << "\n"
        << "\t Total num. of trials: " << total_segments << "\n"
        << "\t Success-rate: "
        << (100.0 * double(quintic_succ_count) / double(total_segments))
        << "\%\n"
        << "\t Graceful-fail-rate: "
        << (100.0 * double(quintic_graceful_fails) / double(total_segments))
        << "\%\n"
        << "\t Hard-fail-rate: "
        << (100.0 * (1.0 - double(quintic_succ_count + quintic_graceful_fails) /
                               double(total_segments)))
        << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_svp_Ndof)) {
    succ_reports
        << "SVP_Ndof interp"
        << "\n"
        << "\t Successes: " << svp_Ndof_succ_count << "\n"
        << "\t Graceful Failures: " << svp_Ndof_graceful_fails << "\n"
        << "\t Total num. of trials: " << total_segments << "\n"
        << "\t Success-rate: "
        << (100.0 * double(svp_Ndof_succ_count) / double(total_segments))
        << "%\n"
        << "\t Graceful-fail-rate: "
        << (100.0 * double(svp_Ndof_graceful_fails) / double(total_segments))
        << "%\n"
        << "\t Hard-fail-rate: "
        << (100.0 *
            (1.0 - double(svp_Ndof_succ_count + svp_Ndof_graceful_fails) /
                       double(total_segments)))
        << "%." << std::endl;
  }
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_sap_Ndof)) {
    succ_reports
        << "SAP_Ndof interp"
        << "\n"
        << "\t Successes: " << sap_Ndof_succ_count << "\n"
        << "\t Graceful Failures: " << sap_Ndof_graceful_fails << "\n"
        << "\t Total num. of trials: " << total_segments << "\n"
        << "\t Success-rate: "
        << (100.0 * double(sap_Ndof_succ_count) / double(total_segments))
        << "%\n"
        << "\t Graceful-fail-rate: "
        << (100.0 * double(sap_Ndof_graceful_fails) / double(total_segments))
        << "%\n"
        << "\t Hard-fail-rate: "
        << (100.0 *
            (1.0 - double(sap_Ndof_succ_count + sap_Ndof_graceful_fails) /
                       double(total_segments)))
        << "%." << std::endl;
  }
#endif

  succ_reports.close();
};

int main(int argc, char** argv) {

  using namespace ReaK;

  absl::ParseCommandLine(argc, argv);

  std::size_t sp_dim = absl::GetFlag(FLAGS_space_dimensionality);

  switch (sp_dim) {
    case 1:
      perform_mc_tests<1>(sp_dim);
      break;
    case 2:
      perform_mc_tests<2>(sp_dim);
      break;
    case 3:
      perform_mc_tests<3>(sp_dim);
      break;
    case 4:
      perform_mc_tests<4>(sp_dim);
      break;
    case 5:
      perform_mc_tests<5>(sp_dim);
      break;
    case 6:
      perform_mc_tests<6>(sp_dim);
      break;
    case 7:
      perform_mc_tests<7>(sp_dim);
      break;
    case 8:
      perform_mc_tests<8>(sp_dim);
      break;
    case 9:
      perform_mc_tests<9>(sp_dim);
      break;
    default:
      perform_mc_tests<0>(sp_dim);
      break;
  };

  return 0;
};
