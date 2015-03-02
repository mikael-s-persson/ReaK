
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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <ReaK/topologies/spaces/time_topology.hpp>
#include <ReaK/topologies/spaces/time_poisson_topology.hpp>
#include <ReaK/topologies/spaces/differentiable_space.hpp>
#include <ReaK/topologies/spaces/temporal_space.hpp>



#define RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#define RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#define RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#define RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#define RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR



#include <ReaK/topologies/spaces/hyperbox_topology.hpp>
#include <ReaK/topologies/spaces/Ndof_spaces.hpp>


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
#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof.hpp>
#endif


#include <ReaK/core/base/scope_guard.hpp>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;



#define RK_TEST_IS_NAN(X) std::isnan(X[0])
#define RK_TEST_IS_INF(X) std::isinf(X[0])
#define RK_TEST_GET_VALUE(X) (X[0])


template <typename Vector>
bool vect_is_nan(const Vector& v) {
  for(std::size_t i = 0; i < v.size(); ++i) 
    if( std::isnan(v[i]) )
      return true;
  return false;
};

template <typename Vector>
bool vect_is_inf(const Vector& v) {
  for(std::size_t i = 0; i < v.size(); ++i) 
    if( std::isinf(v[i]) )
      return true;
  return false;
};


template <typename InterpTrajType, typename Vector, typename TempTopoType, typename PtContainer >
void try_interpolation(const std::string& aMethodName, std::size_t& succ_count,
                       const Vector& curve_ampl, const Vector& curve_phase, 
                       double curve_freq, double interp_steps,
                       const ReaK::shared_ptr< TempTopoType >& topo, const PtContainer& pts, std::ostream& fail_reports) {
  
  using namespace ReaK;
  
  typedef typename pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  try {
    RK_SCOPE_EXIT_ROUTINE(report_construct_except) {
      fail_reports << aMethodName << " exception construct 0 " << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps << std::endl;
    };
    InterpTrajType interp(pts.begin(), pts.end(), topo);
    RK_SCOPE_EXIT_DISMISS(report_construct_except);
    
    for(double t = 0.0; t < 1.0 - 0.5 * interp_steps; ) {
      TempPointType p = pts[0];
      for(std::size_t j = 0; j < 100; ++j) {
        t += 0.01 * interp_steps;
        RK_SCOPE_EXIT_ROUTINE(report_interp_except) {
          fail_reports << aMethodName << " exception interp " << t << " " << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps << std::endl;
        };
        p = interp.get_point_at_time(t);
        RK_SCOPE_EXIT_DISMISS(report_interp_except);
        
        if( vect_is_nan((get<0>(p.pt))) || 
            vect_is_nan((get<1>(p.pt))) || 
            vect_is_nan((get<2>(p.pt))) ) {
          fail_reports << aMethodName << " NaN interp " << t << " " << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps << std::endl;
          throw std::domain_error("NaN condition encountered!");
        };
        
        if( vect_is_inf((get<0>(p.pt))) || 
            vect_is_inf((get<1>(p.pt))) || 
            vect_is_inf((get<2>(p.pt))) ) {
          fail_reports << aMethodName << " INF interp " << t << " " << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps << std::endl;
          throw std::domain_error("INF condition encountered!");
        };
        
      };
      
      Vector ref_pos = curve_ampl;
      for(std::size_t j = 0; j < ref_pos.size(); ++j) {
        ref_pos[j] *= std::sin(curve_freq * t + curve_phase[j]);
      };
      if( norm_2(get<0>(p.pt) - ref_pos) > 1e-3) {
        fail_reports << aMethodName << " pos_tol interp " << t << " " << norm_2(curve_ampl) << " " << curve_freq << " " << interp_steps << std::endl;
        throw std::domain_error("Position-tolerance exceeded!");
      };
      
    };
    
    ++succ_count;
    
  } catch (std::exception& e) { };
};




template <std::size_t StaticSpDim>
struct interp_mc_test_space {
  typedef typename ReaK::pp::Ndof_rl_space< double, StaticSpDim, 2>::type topo_type;
  typedef ReaK::pp::temporal_space< topo_type, ReaK::pp::time_poisson_topology> temp_topo_type;
  typedef ReaK::vect<double, StaticSpDim> vector_type;
  
  static vector_type default_vect(std::size_t ) {
    return vector_type();
  };
  
  static ReaK::shared_ptr< temp_topo_type > create(const vector_type& lb, const vector_type& ub, const vector_type& sb, 
                                                   const vector_type& ab, const vector_type& jb) {
    return ReaK::shared_ptr< temp_topo_type >( new temp_topo_type( "temporal_space",
      ReaK::pp::make_Ndof_rl_space<StaticSpDim>( lb, ub, sb, ab, jb)));
  };
};

template <>
struct interp_mc_test_space<0> {
  typedef ReaK::pp::Ndof_rl_space< double, 0, 2>::type topo_type;
  typedef ReaK::pp::temporal_space< topo_type, ReaK::pp::time_poisson_topology> temp_topo_type;
  typedef ReaK::vect_n<double> vector_type;
  
  static vector_type default_vect(std::size_t dyn_sp_size) {
    return vector_type(dyn_sp_size);
  };
  
  static ReaK::shared_ptr< temp_topo_type > create(const vector_type& lb, const vector_type& ub, const vector_type& sb, 
                                                   const vector_type& ab, const vector_type& jb) {
    return ReaK::shared_ptr< temp_topo_type >( new temp_topo_type( "temporal_space",
      ReaK::pp::make_Ndof_rl_space( lb, ub, sb, ab, jb)));
  };
};


template <std::size_t StaticSpDim>
void perform_mc_tests(const po::variables_map& vm, std::size_t dyn_sp_dim) {
  
  using namespace ReaK;
  
  typedef interp_mc_test_space< StaticSpDim > Config;
  typedef typename Config::topo_type      TopoType;
  typedef typename Config::temp_topo_type TempTopoType;
  typedef typename Config::vector_type    Vector;
  
  typedef typename pp::topology_traits<TopoType>::point_type     PointType;
  typedef typename pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  
  double max_freq     = vm["space-max-frequency"].as<double>();
  double max_rad_freq = max_freq * 2.0 * M_PI;  // rad/s
  
  Vector lb = Config::default_vect(dyn_sp_dim);
  Vector ub = Config::default_vect(dyn_sp_dim);
  Vector sb = Config::default_vect(dyn_sp_dim);
  Vector ab = Config::default_vect(dyn_sp_dim);
  Vector jb = Config::default_vect(dyn_sp_dim);
  for(std::size_t i = 0; i < lb.size(); ++i) {
    lb[i] = -2.0;
    ub[i] =  2.0;
    sb[i] =  2.0 * max_rad_freq;
    ab[i] =  2.0 * max_rad_freq * max_rad_freq;
    jb[i] =  2.0 * max_rad_freq * max_rad_freq * max_rad_freq;
  };
  
  shared_ptr< TempTopoType > topo = Config::create(lb, ub, sb, ab, jb);
  
  
  std::string output_path = vm["output-path"].as<std::string>();
  while(output_path[output_path.length()-1] == '/') 
    output_path.erase(output_path.length()-1, 1);
  
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
  
  std::size_t mc_runs = vm["mc-runs"].as<std::size_t>();
  double interp_steps = vm["interp-steps"].as<double>();
  
  for(std::size_t i = 0; i < mc_runs; ++i) {
    
    std::vector< TempPointType > pts;
    double curve_freq  = double(std::rand() % 1000) * (max_rad_freq / 1000.0);
    
    Vector curve_ampl = Config::default_vect(dyn_sp_dim);
    Vector curve_phase = Config::default_vect(dyn_sp_dim);
    for(std::size_t j = 0; j < curve_ampl.size(); ++j) {
      curve_ampl[j]  = double(std::rand() % 1000) * 0.001;
      curve_phase[j] = double(std::rand() % 1000) * (M_PI / 500.0);
    };
    
    for(double t = 0.0; t < 1.0 + 0.5 * interp_steps; t += interp_steps) {
      Vector pos = curve_ampl;
      Vector vel = curve_ampl;
      Vector acc = curve_ampl;
      
      for(std::size_t j = 0; j < pos.size(); ++j) {
        pos[j] *=  std::sin(curve_freq * t + curve_phase[j]) / (2.0 * max_rad_freq);
        vel[j] *=  curve_freq * std::cos(curve_freq * t + curve_phase[j]) / (2.0 * max_rad_freq * max_rad_freq);
        acc[j] *= -curve_freq * curve_freq * std::sin(curve_freq * t + curve_phase[j]) / (2.0 * max_rad_freq * max_rad_freq * max_rad_freq);
      };
      
      pts.push_back( TempPointType( t, PointType(pos, vel, acc) ) );
    };
    
    
    curve_ampl *= 0.5 / max_rad_freq;
    
    
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
    
    if(vm.count("all-interpolators") || vm.count("linear")) {
      try_interpolation< pp::linear_interp_traj<TempTopoType> >(
        "linear", linear_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    };
    
#endif
    
    
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
    
    if(vm.count("all-interpolators") || vm.count("cubic")) {
      try_interpolation< pp::cubic_hermite_interp_traj<TempTopoType> >(
        "cubic", cubic_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    };
    
#endif
    
    
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
    
    if(vm.count("all-interpolators") || vm.count("quintic")) {
      try_interpolation< pp::quintic_hermite_interp_traj<TempTopoType> >(
        "quintic", quintic_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    };
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
    
    if(vm.count("all-interpolators") || vm.count("svp-Ndof")) {
      try_interpolation< pp::svp_Ndof_interp_traj<TempTopoType> >(
        "svp_Ndof", svp_Ndof_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    };
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
    
    if(vm.count("all-interpolators") || vm.count("sap-Ndof")) {
      try_interpolation< pp::sap_Ndof_interp_traj<TempTopoType> >(
        "sap_Ndof", sap_Ndof_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    };
    
#endif
    
  };
  
  fail_reports.close();
  
  std::ofstream succ_reports((output_path + "/mc_success_rates.txt").c_str());
  succ_reports << "Monte-Carlo runs of interpolation methods using random sinusoidal curves.\n"
               << "  Number of runs: " << mc_runs << "\n"
               << "  Spatial dimensions: " << dyn_sp_dim << "\n"
               << "  Maximum frequency: " << max_freq << " Hz ( " << max_rad_freq << " rad/s )\n"
               << "  Interpolation steps of: " << interp_steps << std::endl;
  
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  if(vm.count("all-interpolators") || vm.count("linear"))
    succ_reports << "Linear interp succeeded " << linear_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(linear_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  if(vm.count("all-interpolators") || vm.count("cubic"))
    succ_reports << "Cubic interp succeeded " << cubic_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(cubic_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  if(vm.count("all-interpolators") || vm.count("quintic"))
    succ_reports << "Quintic interp succeeded " << quintic_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(quintic_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  if(vm.count("all-interpolators") || vm.count("svp-Ndof"))
    succ_reports << "SVP_Ndof interp succeeded " << svp_Ndof_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(svp_Ndof_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  if(vm.count("all-interpolators") || vm.count("sap-Ndof"))
    succ_reports << "SAP_Ndof interp succeeded " << sap_Ndof_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(sap_Ndof_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
  
  succ_reports.close();
  
  
};





int main(int argc, char** argv) {
  
  std::srand( static_cast<unsigned int>( std::time(NULL) ) );
  
  using namespace ReaK;
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("output-path,o", po::value< std::string >()->default_value("test_interp_results"), "specify the output path (default is test_interp_results)")
  ;
  
  po::options_description mc_options("Monte-Carlo options");
  mc_options.add_options()
    ("mc-runs", po::value< std::size_t >()->default_value(100), "number of monte-carlo runs to perform (default is 100)")
  ;
  
  po::options_description space_options("Monte-Carlo options");
  mc_options.add_options()
    ("space-dimensionality", po::value< std::size_t >()->default_value(1), "number of dimensions for the underlying space (default is 1)")
    ("space-max-frequency", po::value< double >()->default_value(10.0), "the maximum frequency of the sinusoidal curves (default is 10.0 Hz)")
    ("interp-steps", po::value< double >()->default_value(0.05), "the time-step between the interpolator's control-points, over a total curve-time of 1.0 second (default is 0.05 seconds)")
  ;
  
  po::options_description interp_select_options("Interpolator selection options");
  interp_select_options.add_options()
    ("all-interpolators,a", "specify that all supported interpolators should be run (default if no particular interpolator is specified)")
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
    ("linear", "specify that the uni-directional RRT algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
    ("cubic", "specify that the bi-directional RRT algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
    ("quintic", "specify that the RRT* algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
    ("svp-Ndof", "specify that the PRM algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
    ("sap-Ndof", "specify that the FADPRM algorithm should be run")
#endif
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(mc_options).add(space_options).add(interp_select_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::size_t sp_dim  = vm["space-dimensionality"].as<std::size_t>();
  
  
  switch(sp_dim) {
    case 1:
      perform_mc_tests<1>(vm, sp_dim);
      break;
    case 2:
      perform_mc_tests<2>(vm, sp_dim);
      break;
    case 3:
      perform_mc_tests<3>(vm, sp_dim);
      break;
    case 4:
      perform_mc_tests<4>(vm, sp_dim);
      break;
    case 5:
      perform_mc_tests<5>(vm, sp_dim);
      break;
    case 6:
      perform_mc_tests<6>(vm, sp_dim);
      break;
    case 7:
      perform_mc_tests<7>(vm, sp_dim);
      break;
    case 8:
      perform_mc_tests<8>(vm, sp_dim);
      break;
    case 9:
      perform_mc_tests<9>(vm, sp_dim);
      break;
    default:
      perform_mc_tests<0>(vm, sp_dim);
      break;
  };
  
  
  return 0;
};















