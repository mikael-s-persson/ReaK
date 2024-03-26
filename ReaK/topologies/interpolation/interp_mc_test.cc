
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

#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/temporal_space.h"
#include "ReaK/topologies/spaces/time_poisson_topology.h"
#include "ReaK/topologies/spaces/time_topology.h"

#define RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#define RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#define RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#define RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#define RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR

#include "ReaK/topologies/spaces/hyperbox_topology.h"
#include "ReaK/topologies/spaces/ndof_spaces.h"

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#include "ReaK/topologies/interpolation/linear_interp.h"
#endif

#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#include "ReaK/topologies/interpolation/cubic_hermite_interp.h"
#endif

#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#include "ReaK/topologies/interpolation/quintic_hermite_interp.h"
#endif

#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#include "ReaK/topologies/interpolation/sustained_velocity_pulse_ndof.h"
#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#include "ReaK/topologies/interpolation/sustained_acceleration_pulse_ndof.h"
#endif

#include "ReaK/core/base/scope_guard.h"

#include <filesystem>
#include <memory>

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
}

template <typename Vector>
bool vect_is_inf(const Vector& v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (std::isinf(v[i])) {
      return true;
    }
  }
  return false;
}

template <typename InterpTrajType, typename Vector, typename TempTopoType,
          typename PtContainer>
void try_interpolation(const std::string& aMethodName, std::size_t& succ_count,
                       const Vector& curve_ampl, const Vector& curve_phase,
                       double curve_freq, double interp_steps,
                       const std::shared_ptr<TempTopoType>& topo,
                       const PtContainer& pts, std::ostream& fail_reports) {

  using namespace ReaK;

  using TempPointType = typename pp::topology_traits<TempTopoType>::point_type;

  try {
    RK_SCOPE_EXIT_ROUTINE(report_construct_except) {
      fail_reports << aMethodName << " exception construct 0 "
                   << norm_2(curve_ampl) << " " << curve_freq << " "
                   << interp_steps << std::endl;
    };
    InterpTrajType interp(pts.begin(), pts.end(), topo);
    RK_SCOPE_EXIT_DISMISS(report_construct_except);

    for (double t = 0.0; t < 1.0 - 0.5 * interp_steps;) {
      TempPointType p = pts[0];
      for (std::size_t j = 0; j < 100; ++j) {
        t += 0.01 * interp_steps;
        RK_SCOPE_EXIT_ROUTINE(report_interp_except) {
          fail_reports << aMethodName << " exception interp " << t << " "
                       << norm_2(curve_ampl) << " " << curve_freq << " "
                       << interp_steps << std::endl;
        };
        p = interp.get_point_at_time(t);
        RK_SCOPE_EXIT_DISMISS(report_interp_except);

        if (vect_is_nan((get<0>(p.pt))) || vect_is_nan((get<1>(p.pt))) ||
            vect_is_nan((get<2>(p.pt)))) {
          fail_reports << aMethodName << " NaN interp " << t << " "
                       << norm_2(curve_ampl) << " " << curve_freq << " "
                       << interp_steps << std::endl;
          throw std::domain_error("NaN condition encountered!");
        }

        if (vect_is_inf((get<0>(p.pt))) || vect_is_inf((get<1>(p.pt))) ||
            vect_is_inf((get<2>(p.pt)))) {
          fail_reports << aMethodName << " INF interp " << t << " "
                       << norm_2(curve_ampl) << " " << curve_freq << " "
                       << interp_steps << std::endl;
          throw std::domain_error("INF condition encountered!");
        }
      }

      Vector ref_pos = curve_ampl;
      for (std::size_t j = 0; j < ref_pos.size(); ++j) {
        ref_pos[j] *= std::sin(curve_freq * t + curve_phase[j]);
      }
      if (norm_2(get<0>(p.pt) - ref_pos) > 1e-3) {
        fail_reports << aMethodName << " pos_tol interp " << t << " "
                     << norm_2(curve_ampl) << " " << curve_freq << " "
                     << interp_steps << std::endl;
        throw std::domain_error("Position-tolerance exceeded!");
      }
    }

    ++succ_count;

  } catch (std::exception& e) {
    RK_UNUSED(e);
  }
}

template <std::size_t StaticSpDim>
struct interp_mc_test_space {
  using topo_type =
      typename ReaK::pp::Ndof_rl_space<double, StaticSpDim, 2>::type;
  using temp_topo_type =
      ReaK::pp::temporal_space<topo_type, ReaK::pp::time_poisson_topology>;
  using vector_type = ReaK::vect<double, StaticSpDim>;

  static vector_type default_vect(std::size_t /*unused*/) {
    return vector_type();
  };

  static std::shared_ptr<temp_topo_type> create(const vector_type& lb,
                                                const vector_type& ub,
                                                const vector_type& sb,
                                                const vector_type& ab,
                                                const vector_type& jb) {
    return std::shared_ptr<temp_topo_type>(new temp_topo_type(
        "temporal_space",
        ReaK::pp::make_Ndof_rl_space<StaticSpDim>(lb, ub, sb, ab, jb)));
  }
};

template <>
struct interp_mc_test_space<0> {
  using topo_type = ReaK::pp::Ndof_rl_space<double, 0, 2>::type;
  using temp_topo_type =
      ReaK::pp::temporal_space<topo_type, ReaK::pp::time_poisson_topology>;
  using vector_type = ReaK::vect_n<double>;

  static vector_type default_vect(std::size_t dyn_sp_size) {
    return vector_type(dyn_sp_size);
  }

  static std::shared_ptr<temp_topo_type> create(const vector_type& lb,
                                                const vector_type& ub,
                                                const vector_type& sb,
                                                const vector_type& ab,
                                                const vector_type& jb) {
    return std::make_shared<temp_topo_type>(
        "temporal_space", ReaK::pp::make_Ndof_rl_space(lb, ub, sb, ab, jb));
  }
};

template <std::size_t StaticSpDim>
void perform_mc_tests(std::size_t dyn_sp_dim) {

  using namespace ReaK;

  using Config = interp_mc_test_space<StaticSpDim>;
  using TopoType = typename Config::topo_type;
  using TempTopoType = typename Config::temp_topo_type;
  using Vector = typename Config::vector_type;

  using PointType = typename pp::topology_traits<TopoType>::point_type;
  using TempPointType = typename pp::topology_traits<TempTopoType>::point_type;

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
  }

  std::shared_ptr<TempTopoType> topo = Config::create(lb, ub, sb, ab, jb);

  std::string output_path = absl::GetFlag(FLAGS_output_path);
  while (output_path.back() == '/') {
    output_path.erase(output_path.length() - 1, 1);
  }

  fs::create_directory(output_path.c_str());

  std::ofstream fail_reports((output_path + "/mc_fail_reports.txt").c_str());

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  std::size_t linear_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  std::size_t cubic_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  std::size_t quintic_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  std::size_t svp_Ndof_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  std::size_t sap_Ndof_succ_count = 0;
#endif

  std::size_t mc_runs = absl::GetFlag(FLAGS_mc_runs);
  double interp_steps = absl::GetFlag(FLAGS_interp_steps);

  global_rng_type& gbl_rng = get_global_rng();

  for (std::size_t i = 0; i < mc_runs; ++i) {

    std::vector<TempPointType> pts;
    double curve_freq = double(gbl_rng() % 1000) * (max_rad_freq / 1000.0);

    Vector curve_ampl = Config::default_vect(dyn_sp_dim);
    Vector curve_phase = Config::default_vect(dyn_sp_dim);
    for (std::size_t j = 0; j < curve_ampl.size(); ++j) {
      curve_ampl[j] = double(gbl_rng() % 1000) * 0.001;
      curve_phase[j] = double(gbl_rng() % 1000) * (M_PI / 500.0);
    }

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
      }

      pts.push_back(TempPointType(t, PointType(pos, vel, acc)));
    }

    curve_ampl *= 0.5 / max_rad_freq;

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_linear)) {
      try_interpolation<pp::linear_interp_traj<TempTopoType>>(
          "linear", linear_succ_count, curve_ampl, curve_phase, curve_freq,
          interp_steps, topo, pts, fail_reports);
    }

#endif

#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_cubic)) {
      try_interpolation<pp::cubic_hermite_interp_traj<TempTopoType>>(
          "cubic", cubic_succ_count, curve_ampl, curve_phase, curve_freq,
          interp_steps, topo, pts, fail_reports);
    }

#endif

#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_quintic)) {
      try_interpolation<pp::quintic_hermite_interp_traj<TempTopoType>>(
          "quintic", quintic_succ_count, curve_ampl, curve_phase, curve_freq,
          interp_steps, topo, pts, fail_reports);
    }

#endif

#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_svp_Ndof)) {
      try_interpolation<pp::svp_Ndof_interp_traj<TempTopoType>>(
          "svp_ndof", svp_Ndof_succ_count, curve_ampl, curve_phase, curve_freq,
          interp_steps, topo, pts, fail_reports);
    }

#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR

    if (absl::GetFlag(FLAGS_all_interpolators) ||
        absl::GetFlag(FLAGS_sap_Ndof)) {
      try_interpolation<pp::sap_Ndof_interp_traj<TempTopoType>>(
          "sap_ndof", sap_Ndof_succ_count, curve_ampl, curve_phase, curve_freq,
          interp_steps, topo, pts, fail_reports);
    }

#endif
  }

  fail_reports.close();

  std::ofstream succ_reports((output_path + "/mc_success_rates.txt").c_str());
  succ_reports << "Monte-Carlo runs of interpolation methods using random "
                  "sinusoidal curves.\n"
               << "  Number of runs: " << mc_runs << "\n"
               << "  Spatial dimensions: " << dyn_sp_dim << "\n"
               << "  Maximum frequency: " << max_freq << " Hz ( "
               << max_rad_freq << " rad/s )\n"
               << "  Interpolation steps of: " << interp_steps << std::endl;

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_linear)) {
    succ_reports << "Linear interp succeeded " << linear_succ_count
                 << " out of " << mc_runs << " which is "
                 << (100.0 * double(linear_succ_count) / double(mc_runs))
                 << "% success-rate." << std::endl;
  }
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_cubic)) {
    succ_reports << "Cubic interp succeeded " << cubic_succ_count << " out of "
                 << mc_runs << " which is "
                 << (100.0 * double(cubic_succ_count) / double(mc_runs))
                 << "% success-rate." << std::endl;
  }
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_quintic)) {
    succ_reports << "Quintic interp succeeded " << quintic_succ_count
                 << " out of " << mc_runs << " which is "
                 << (100.0 * double(quintic_succ_count) / double(mc_runs))
                 << "% success-rate." << std::endl;
  }
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_svp_Ndof)) {
    succ_reports << "SVP_Ndof interp succeeded " << svp_Ndof_succ_count
                 << " out of " << mc_runs << " which is "
                 << (100.0 * double(svp_Ndof_succ_count) / double(mc_runs))
                 << "% success-rate." << std::endl;
  }
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  if (absl::GetFlag(FLAGS_all_interpolators) || absl::GetFlag(FLAGS_sap_Ndof)) {
    succ_reports << "SAP_Ndof interp succeeded " << sap_Ndof_succ_count
                 << " out of " << mc_runs << " which is "
                 << (100.0 * double(sap_Ndof_succ_count) / double(mc_runs))
                 << "% success-rate." << std::endl;
  }
#endif

  succ_reports.close();
}

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
  }

  return 0;
}
