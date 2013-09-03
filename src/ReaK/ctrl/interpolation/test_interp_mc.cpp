
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

#include "topologies/time_topology.hpp"
#include "topologies/time_poisson_topology.hpp"
#include "topologies/differentiable_space.hpp"
#include "topologies/temporal_space.hpp"



#define RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#define RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#define RK_ENABLE_TEST_QUINTIC_INTERPOLATOR


#define RK_INTERP_MC_TEST_USE_HYPERBOX
#define RK_INTERP_MC_TEST_USE_RL_NDOF_SPACE


#ifndef RK_INTERP_MC_TEST_USE_HYPERBOX
#define RK_ENABLE_TEST_SVP_INTERPOLATOR
#define RK_ENABLE_TEST_SAP_INTERPOLATOR
#else
#define RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#define RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#endif



#ifndef RK_INTERP_MC_TEST_USE_HYPERBOX
#include "topologies/line_topology.hpp"
#else
#include "topologies/hyperbox_topology.hpp"
#include "topologies/Ndof_spaces.hpp"
#endif


#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
#include "linear_interp.hpp"
#endif

#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
#include "cubic_hermite_interp.hpp"
#endif

#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
#include "quintic_hermite_interp.hpp"
#endif

#ifdef RK_ENABLE_TEST_SVP_INTERPOLATOR
#include "sustained_velocity_pulse.hpp"
#endif

#ifdef RK_ENABLE_TEST_SAP_INTERPOLATOR
#include "sustained_acceleration_pulse.hpp"
#endif

#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
#include "sustained_velocity_pulse_Ndof.hpp"
#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#include "sustained_acceleration_pulse_Ndof.hpp"
#endif


#include "base/scope_guard.hpp"


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;



#ifndef RK_INTERP_MC_TEST_USE_HYPERBOX
#define RK_TEST_IS_NAN(X) std::isnan(X)
#define RK_TEST_IS_INF(X) std::isinf(X)
#define RK_TEST_GET_VALUE(X) (X)
#else
#define RK_TEST_IS_NAN(X) std::isnan(X[0])
#define RK_TEST_IS_INF(X) std::isinf(X[0])
#define RK_TEST_GET_VALUE(X) (X[0])
#endif



template <typename InterpTrajType, typename TempTopoType, typename PtContainer >
void try_interpolation(const std::string& aMethodName, std::size_t& succ_count,
                       double curve_ampl, double curve_phase, double curve_freq, double interp_steps,
                       const ReaK::shared_ptr< TempTopoType >& topo, const PtContainer& pts, std::ostream& fail_reports) {
  
  using namespace ReaK;
  
  typedef typename pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  try {
    RK_SCOPE_EXIT_ROUTINE(report_construct_except) {
      fail_reports << aMethodName << " exception construct 0 " << curve_ampl << " " << curve_phase << " " << curve_freq << " " << interp_steps << std::endl;
    };
    InterpTrajType interp(pts.begin(), pts.end(), topo);
    RK_SCOPE_EXIT_DISMISS(report_construct_except);
    
    for(double t = 0.0; t < 1.0 - 0.5 * interp_steps; ) {
      TempPointType p = pts[0];
      for(std::size_t j = 0; j < 100; ++j) {
        t += 0.01 * interp_steps;
        RK_SCOPE_EXIT_ROUTINE(report_interp_except) {
          fail_reports << aMethodName << " exception interp " << t << " " << curve_ampl << " " << curve_phase << " " << curve_freq << " " << interp_steps << std::endl;
        };
        p = interp.get_point_at_time(t);
        RK_SCOPE_EXIT_DISMISS(report_interp_except);
        
        if( RK_TEST_IS_NAN((get<0>(p.pt))) || 
            RK_TEST_IS_NAN((get<1>(p.pt))) || 
            RK_TEST_IS_NAN((get<2>(p.pt))) ) {
          fail_reports << aMethodName << " NaN interp " << t << " " << curve_ampl << " " << curve_phase << " " << curve_freq << " " << interp_steps << std::endl;
          throw std::domain_error("NaN condition encountered!");
        };
        
        if( RK_TEST_IS_INF((get<0>(p.pt))) || 
            RK_TEST_IS_INF((get<1>(p.pt))) || 
            RK_TEST_IS_INF((get<2>(p.pt))) ) {
          fail_reports << aMethodName << " INF interp " << t << " " << curve_ampl << " " << curve_phase << " " << curve_freq << " " << interp_steps << std::endl;
          throw std::domain_error("INF condition encountered!");
        };
        
      };
      
      if(std::fabs( RK_TEST_GET_VALUE(get<0>(p.pt)) - curve_ampl * std::sin(curve_freq * t + curve_phase)) > 1e-3) {
//         std::cout << aMethodName << " failed tolerance at time " << t 
//                   << " with value " << get<0>(p.pt) << " instead of " << (curve_ampl * std::sin(curve_freq * t + curve_phase)) << std::endl;
        fail_reports << aMethodName << " pos_tol interp " << t << " " << curve_ampl << " " << curve_phase << " " << curve_freq << " " << interp_steps << std::endl;
        throw std::domain_error("Position-tolerance exceeded!");
      };
      
    };
    
    ++succ_count;
    
  } catch (std::exception& e) { };
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
#ifdef RK_ENABLE_TEST_SVP_INTERPOLATOR
    ("svp", "specify that the PRM algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_SAP_INTERPOLATOR
    ("sap", "specify that the FADPRM algorithm should be run")
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
  
  std::string output_path = vm["output-path"].as<std::string>();
  while(output_path[output_path.length()-1] == '/') 
    output_path.erase(output_path.length()-1, 1);
  
  fs::create_directory(output_path.c_str());
  
  
  std::size_t mc_runs = vm["mc-runs"].as<std::size_t>();
  
  std::size_t sp_dim  = vm["space-dimensionality"].as<std::size_t>();
  // NOTE: temporary..
  if(sp_dim > 1) {
    std::cout << "Sorry, only dimensionality of 1 is supported for the moment.." << std::endl;
    return 1;
  };
  
  double max_freq     = vm["space-max-frequency"].as<double>();
  double max_rad_freq = max_freq * 2.0 * M_PI;  // rad/s
  double interp_steps = vm["interp-steps"].as<double>();
  
  
  std::ofstream fail_reports((output_path + "/mc_fail_reports.txt").c_str());
  
#ifndef RK_INTERP_MC_TEST_USE_HYPERBOX
  typedef arithmetic_tuple< 
            pp::line_segment_topology<double>, 
            pp::line_segment_topology<double>, 
            pp::line_segment_topology<double>, 
            pp::line_segment_topology<double> 
          > SpaceTupleType;

  typedef pp::differentiable_space< pp::time_poisson_topology, SpaceTupleType > TopoType;
  typedef pp::topology_traits<TopoType>::point_type PointType;
  typedef pp::temporal_space< TopoType, pp::time_poisson_topology> TempTopoType;
  typedef pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  shared_ptr< TempTopoType > topo = 
    shared_ptr< TempTopoType >( new TempTopoType( "temporal_space",
      SpaceTupleType(pp::line_segment_topology<double>("pos_topo", -2.0, 2.0),
                     pp::line_segment_topology<double>("vel_topo", -2.0 * max_rad_freq, 2.0 * max_rad_freq),
                     pp::line_segment_topology<double>("acc_topo", -2.0 * max_rad_freq * max_rad_freq, 2.0 * max_rad_freq * max_rad_freq),
                     pp::line_segment_topology<double>("jerk_topo",-2.0 * max_rad_freq * max_rad_freq * max_rad_freq, 2.0 * max_rad_freq * max_rad_freq * max_rad_freq))));
  
#else
  
#ifndef RK_INTERP_MC_TEST_USE_RL_NDOF_SPACE
  
  typedef arithmetic_tuple< 
            pp::hyperbox_topology< vect<double, 1> >, 
            pp::hyperbox_topology< vect<double, 1> >, 
            pp::hyperbox_topology< vect<double, 1> >, 
            pp::hyperbox_topology< vect<double, 1> > 
          > SpaceTupleType;

  typedef pp::differentiable_space< pp::time_poisson_topology, SpaceTupleType > TopoType;
  typedef pp::topology_traits<TopoType>::point_type PointType;
  typedef pp::temporal_space< TopoType, pp::time_poisson_topology> TempTopoType;
  typedef pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  shared_ptr< TempTopoType > topo = 
    shared_ptr< TempTopoType >( new TempTopoType( "temporal_space",
      SpaceTupleType(pp::hyperbox_topology< vect<double, 1> >("pos_topo", vect<double, 1>(-2.0), vect<double, 1>(2.0)),
                     pp::hyperbox_topology< vect<double, 1> >("vel_topo", vect<double, 1>(-2.0 * max_rad_freq), vect<double, 1>(2.0 * max_rad_freq)),
                     pp::hyperbox_topology< vect<double, 1> >("acc_topo", vect<double, 1>(-2.0 * max_rad_freq * max_rad_freq), vect<double, 1>(2.0 * max_rad_freq * max_rad_freq)),
                     pp::hyperbox_topology< vect<double, 1> >("jerk_topo",vect<double, 1>(-2.0 * max_rad_freq * max_rad_freq * max_rad_freq), vect<double, 1>(2.0 * max_rad_freq * max_rad_freq * max_rad_freq)))));
  
#else
  
  typedef pp::Ndof_rl_space< double, 1, 2>::type TopoType;
  typedef pp::topology_traits<TopoType>::point_type PointType;
  typedef pp::temporal_space< TopoType, pp::time_poisson_topology> TempTopoType;
  typedef pp::topology_traits<TempTopoType>::point_type TempPointType;
  
  shared_ptr< TempTopoType > topo = 
    shared_ptr< TempTopoType >( new TempTopoType( "temporal_space",
      pp::make_Ndof_rl_space<1, vect<double, 1> >(
        vect<double, 1>(-2.0), vect<double, 1>(2.0), 
        vect<double, 1>(2.0 * max_rad_freq),
        vect<double, 1>(2.0 * max_rad_freq * max_rad_freq),
        vect<double, 1>(2.0 * max_rad_freq * max_rad_freq * max_rad_freq))));
  
#endif
  
#endif
  
  
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  std::size_t linear_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  std::size_t cubic_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  std::size_t quintic_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SVP_INTERPOLATOR
  std::size_t svp_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SAP_INTERPOLATOR
  std::size_t sap_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  std::size_t svp_Ndof_succ_count = 0;
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  std::size_t sap_Ndof_succ_count = 0;
#endif
  
  for(std::size_t i = 0; i < mc_runs; ++i) {
    
    std::vector< TempPointType > pts;
    double curve_ampl  = double(std::rand() % 1000) * 0.001;
    double curve_phase = double(std::rand() % 1000) * (M_PI / 500.0);
    double curve_freq  = double(std::rand() % 1000) * (max_rad_freq / 1000.0);
    
    for(double t = 0.0; t < 1.0 + 0.5 * interp_steps; t += interp_steps) {
#ifndef RK_INTERP_MC_TEST_USE_HYPERBOX
      pts.push_back( TempPointType( t, PointType(
         curve_ampl * std::sin(curve_freq * t + curve_phase), 
         curve_ampl * curve_freq * std::cos(curve_freq * t + curve_phase), 
        -curve_ampl * curve_freq * curve_freq * std::sin(curve_freq * t + curve_phase), 
        -curve_ampl * curve_freq * curve_freq * curve_freq * std::cos(curve_freq * t + curve_phase))) );
#else
#ifndef RK_INTERP_MC_TEST_USE_RL_NDOF_SPACE
      pts.push_back( TempPointType( t, PointType(
        vect<double, 1>( curve_ampl * std::sin(curve_freq * t + curve_phase)), 
        vect<double, 1>( curve_ampl * curve_freq * std::cos(curve_freq * t + curve_phase)), 
        vect<double, 1>(-curve_ampl * curve_freq * curve_freq * std::sin(curve_freq * t + curve_phase)), 
        vect<double, 1>(-curve_ampl * curve_freq * curve_freq * curve_freq * std::cos(curve_freq * t + curve_phase)))) );
#else
      pts.push_back( TempPointType( t, PointType(
        vect<double, 1>( curve_ampl / (2.0 * max_rad_freq) * std::sin(curve_freq * t + curve_phase)), 
        vect<double, 1>( curve_ampl * curve_freq / (2.0 * max_rad_freq * max_rad_freq) * std::cos(curve_freq * t + curve_phase)), 
        vect<double, 1>(-curve_ampl * curve_freq * curve_freq / (2.0 * max_rad_freq * max_rad_freq * max_rad_freq) * std::sin(curve_freq * t + curve_phase)))) );
#endif
#endif
    };
    
    
#ifdef RK_INTERP_MC_TEST_USE_RL_NDOF_SPACE
    curve_ampl /= 2.0 * max_rad_freq;
#endif
    
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
    
    try_interpolation< pp::linear_interp_traj<TempTopoType> >(
      "linear", linear_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
    
    try_interpolation< pp::cubic_hermite_interp_traj<TempTopoType> >(
      "cubic", cubic_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
    
    try_interpolation< pp::quintic_hermite_interp_traj<TempTopoType> >(
      "quintic", quintic_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SVP_INTERPOLATOR
    
    try_interpolation< pp::svp_interp_traj<TempTopoType> >(
      "svp", svp_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SAP_INTERPOLATOR
    
    try_interpolation< pp::sap_interp_traj<TempTopoType> >(
      "sap", sap_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
    
    try_interpolation< pp::svp_Ndof_interp_traj<TempTopoType> >(
      "svp_Ndof", svp_Ndof_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
    
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
    
    try_interpolation< pp::sap_Ndof_interp_traj<TempTopoType> >(
      "sap_Ndof", sap_Ndof_succ_count, curve_ampl, curve_phase, curve_freq, interp_steps, topo, pts, fail_reports);
    
#endif
    
  };
  
  fail_reports.close();
  
  std::ofstream succ_reports((output_path + "/mc_success_rates.txt").c_str());
  succ_reports << "Monte-Carlo runs of interpolation methods using random sinusoidal curves.\n"
               << "  Number of runs: " << mc_runs << "\n"
               << "  Spatial dimensions: " << sp_dim << "\n"
               << "  Maximum frequency: " << max_freq << " Hz ( " << max_rad_freq << " rad/s )\n"
               << "  Interpolation steps of: " << interp_steps << std::endl;
  
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  succ_reports << "Linear interp succeeded " << linear_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(linear_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  succ_reports << "Cubic interp succeeded " << cubic_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(cubic_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  succ_reports << "Quintic interp succeeded " << quintic_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(quintic_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SVP_INTERPOLATOR
  succ_reports << "SVP interp succeeded " << svp_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(svp_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SAP_INTERPOLATOR
  succ_reports << "SAP interp succeeded " << sap_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(sap_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  succ_reports << "SVP_Ndof interp succeeded " << svp_Ndof_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(svp_Ndof_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  succ_reports << "SAP_Ndof interp succeeded " << sap_Ndof_succ_count << " out of " << mc_runs << " which is " << (100.0 * double(sap_Ndof_succ_count) / double(mc_runs)) << "\% success-rate." << std::endl;
#endif
  
  succ_reports.close();

  return 0;
};















