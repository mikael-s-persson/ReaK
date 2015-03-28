
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
#include <ReaK/topologies/spaces/differentiable_space.hpp>


// #define RK_ENABLE_TEST_LINEAR_INTERPOLATOR
// #define RK_ENABLE_TEST_CUBIC_INTERPOLATOR
// #define RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
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
#include <ReaK/topologies/interpolation/svp_Ndof_reach_topologies.hpp>
#endif

#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof.hpp>
#include <ReaK/topologies/interpolation/sap_Ndof_reach_topologies.hpp>
#endif


#include <ReaK/core/base/scope_guard.hpp>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


#define RK_TEST_IS_NAN( X ) std::isnan( X[0] )
#define RK_TEST_IS_INF( X ) std::isinf( X[0] )
#define RK_TEST_GET_VALUE( X ) ( X[0] )


template < typename Vector >
bool vect_is_nan( const Vector& v ) {
  for( std::size_t i = 0; i < v.size(); ++i )
    if( std::isnan( v[i] ) )
      return true;
  return false;
};

template < typename Vector >
bool vect_is_inf( const Vector& v ) {
  for( std::size_t i = 0; i < v.size(); ++i )
    if( std::isinf( v[i] ) )
      return true;
  return false;
};


template < typename InterpTopoType >
void try_interpolation( const std::string& aMethodName, std::size_t mc_runs, std::size_t& succ_count,
                        std::size_t& graceful_failures, const InterpTopoType& topo, std::ostream& fail_reports ) {

  using namespace ReaK;

  typedef typename pp::topology_traits< InterpTopoType >::point_type PointType;

  for( std::size_t i = 0; i < mc_runs; ++i ) {

    PointType p1 = get( pp::random_sampler, topo )( topo );
    PointType p2 = get( pp::random_sampler, topo )( topo );

    try {

      RK_SCOPE_EXIT_ROUTINE( report_dist_except ) { fail_reports << aMethodName << " exception dist" << std::endl; };
      double d = get( pp::distance_metric, topo )( p1, p2, topo );
      RK_SCOPE_EXIT_DISMISS( report_dist_except );

      if( std::isnan( d ) ) {
        fail_reports << aMethodName << " NaN dist" << std::endl;
        throw std::domain_error( "NaN condition encountered!" );
      };

      if( std::isinf( d ) ) {
        ++graceful_failures;
        continue;
      };

      RK_SCOPE_EXIT_ROUTINE( report_interp_except ) {
        fail_reports << aMethodName << " exception interp" << std::endl;
      };
      PointType p = topo.move_position_toward( p1, 1.0, p2 );
      RK_SCOPE_EXIT_DISMISS( report_interp_except );

      if( vect_is_nan( ( get< 0 >( p ) ) ) || vect_is_nan( ( get< 1 >( p ) ) ) || vect_is_nan( ( get< 2 >( p ) ) ) ) {
        fail_reports << aMethodName << " NaN interp" << std::endl;
        throw std::domain_error( "NaN condition encountered!" );
      };

      if( vect_is_inf( ( get< 0 >( p ) ) ) || vect_is_inf( ( get< 1 >( p ) ) ) || vect_is_inf( ( get< 2 >( p ) ) ) ) {
        fail_reports << aMethodName << " INF interp" << std::endl;
        throw std::domain_error( "INF condition encountered!" );
      };

      if( norm_2( get< 0 >( p ) - get< 0 >( p2 ) ) > 1e-3 ) {
        fail_reports << aMethodName << " pos_tol interp " << norm_2( get< 0 >( p ) - get< 0 >( p2 ) ) << std::endl;
        throw std::domain_error( "Position-tolerance exceeded!" );
      };

      if( norm_2( get< 1 >( p ) - get< 1 >( p2 ) ) > 1e-3 ) {
        fail_reports << aMethodName << " vel_tol interp " << norm_2( get< 1 >( p ) - get< 1 >( p2 ) ) << std::endl;
        throw std::domain_error( "Velocity-tolerance exceeded!" );
      };

      ++succ_count;

    } catch( std::exception& e ) {
    };
  };
};


template < std::size_t StaticSpDim >
struct interp_mc_test_space {
  typedef typename ReaK::pp::Ndof_rl_space< double, StaticSpDim, 2 >::type topo_type;
  typedef ReaK::vect< double, StaticSpDim > vector_type;

  static vector_type default_vect( std::size_t ) { return vector_type(); };

  static topo_type create( const vector_type& lb, const vector_type& ub, const vector_type& sb, const vector_type& ab,
                           const vector_type& jb ) {
    return ReaK::pp::make_Ndof_rl_space< StaticSpDim >( lb, ub, sb, ab, jb );
  };
};

template <>
struct interp_mc_test_space< 0 > {
  typedef ReaK::pp::Ndof_rl_space< double, 0, 2 >::type topo_type;
  typedef ReaK::vect_n< double > vector_type;

  static vector_type default_vect( std::size_t dyn_sp_size ) { return vector_type( dyn_sp_size ); };

  static topo_type create( const vector_type& lb, const vector_type& ub, const vector_type& sb, const vector_type& ab,
                           const vector_type& jb ) {
    return ReaK::pp::make_Ndof_rl_space( lb, ub, sb, ab, jb );
  };
};


template < std::size_t StaticSpDim >
void perform_mc_tests( const po::variables_map& vm, std::size_t dyn_sp_dim ) {

  using namespace ReaK;

  typedef interp_mc_test_space< StaticSpDim > Config;
  typedef typename Config::topo_type TopoType;
  typedef typename Config::vector_type Vector;

  Vector lb = Config::default_vect( dyn_sp_dim );
  Vector ub = Config::default_vect( dyn_sp_dim );
  Vector sb = Config::default_vect( dyn_sp_dim );
  Vector ab = Config::default_vect( dyn_sp_dim );
  Vector jb = Config::default_vect( dyn_sp_dim );
  for( std::size_t i = 0; i < lb.size(); ++i ) {
    lb[i] = -1.0;
    ub[i] = 1.0;
    sb[i] = 1.0;
    ab[i] = 1.0;
    jb[i] = 1.0;
  };

  std::string output_path = vm["output-path"].as< std::string >();
  while( output_path[output_path.length() - 1] == '/' )
    output_path.erase( output_path.length() - 1, 1 );

  fs::create_directory( output_path.c_str() );

  std::ofstream fail_reports( ( output_path + "/topo_rand_fail_reports.txt" ).c_str() );


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
  pp::interpolated_topology< TopoType, pp::svp_Ndof_interpolation_tag > svp_Ndof_topo(
    Config::create( lb, ub, sb, ab, jb ) );
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  std::size_t sap_Ndof_succ_count = 0;
  std::size_t sap_Ndof_graceful_fails = 0;
  pp::interpolated_topology< TopoType, pp::sap_Ndof_interpolation_tag > sap_Ndof_topo(
    Config::create( lb, ub, sb, ab, jb ) );
#endif

  std::size_t mc_runs = vm["mc-runs"].as< std::size_t >();

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR

  if( vm.count( "all-interpolators" ) || vm.count( "linear" ) ) {
    try_interpolation( "linear", mc_runs, linear_succ_count, linear_graceful_fails, linear_topo, fail_reports );
  };

#endif


#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR

  if( vm.count( "all-interpolators" ) || vm.count( "cubic" ) ) {
    try_interpolation( "cubic", mc_runs, cubic_succ_count, cubic_graceful_fails, cubic_topo, fail_reports );
  };

#endif


#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR

  if( vm.count( "all-interpolators" ) || vm.count( "quintic" ) ) {
    try_interpolation( "quintic", mc_runs, quintic_succ_count, quintic_graceful_fails, quintic_topo, fail_reports );
  };

#endif


#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR

  if( vm.count( "all-interpolators" ) || vm.count( "svp-Ndof" ) ) {
    try_interpolation( "svp_Ndof", mc_runs, svp_Ndof_succ_count, svp_Ndof_graceful_fails, svp_Ndof_topo, fail_reports );
  };

#endif


#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR

  if( vm.count( "all-interpolators" ) || vm.count( "sap-Ndof" ) ) {
    try_interpolation( "sap_Ndof", mc_runs, sap_Ndof_succ_count, sap_Ndof_graceful_fails, sap_Ndof_topo, fail_reports );
  };

#endif

  fail_reports.close();

  std::ofstream succ_reports( ( output_path + "/topo_rand_success_rates.txt" ).c_str() );
  succ_reports << "Monte-Carlo runs of interpolation methods using random samples.\n"
               << "  Number of runs: " << mc_runs << "\n"
               << "  Spatial dimensions: " << dyn_sp_dim << std::endl;

#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
  if( vm.count( "all-interpolators" ) || vm.count( "linear" ) )
    succ_reports << "Linear interp"
                 << "\n"
                 << "\t Successes: " << linear_succ_count << "\n"
                 << "\t Graceful Failures: " << linear_graceful_fails << "\n"
                 << "\t Total num. of trials: " << mc_runs << "\n"
                 << "\t Success-rate: " << ( 100.0 * double( linear_succ_count ) / double( mc_runs ) ) << "\%\n"
                 << "\t Graceful-fail-rate: " << ( 100.0 * double( linear_graceful_fails ) / double( mc_runs ) )
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << ( 100.0 * ( 1.0 - double( linear_succ_count + linear_graceful_fails ) / double( mc_runs ) ) )
                 << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
  if( vm.count( "all-interpolators" ) || vm.count( "cubic" ) )
    succ_reports << "Cubic interp"
                 << "\n"
                 << "\t Successes: " << cubic_succ_count << "\n"
                 << "\t Graceful Failures: " << cubic_graceful_fails << "\n"
                 << "\t Total num. of trials: " << mc_runs << "\n"
                 << "\t Success-rate: " << ( 100.0 * double( cubic_succ_count ) / double( mc_runs ) ) << "\%\n"
                 << "\t Graceful-fail-rate: " << ( 100.0 * double( cubic_graceful_fails ) / double( mc_runs ) )
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << ( 100.0 * ( 1.0 - double( cubic_succ_count + cubic_graceful_fails ) / double( mc_runs ) ) ) << "\%."
                 << std::endl;
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
  if( vm.count( "all-interpolators" ) || vm.count( "quintic" ) )
    succ_reports << "Quintic interp"
                 << "\n"
                 << "\t Successes: " << quintic_succ_count << "\n"
                 << "\t Graceful Failures: " << quintic_graceful_fails << "\n"
                 << "\t Total num. of trials: " << mc_runs << "\n"
                 << "\t Success-rate: " << ( 100.0 * double( quintic_succ_count ) / double( mc_runs ) ) << "\%\n"
                 << "\t Graceful-fail-rate: " << ( 100.0 * double( quintic_graceful_fails ) / double( mc_runs ) )
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << ( 100.0 * ( 1.0 - double( quintic_succ_count + quintic_graceful_fails ) / double( mc_runs ) ) )
                 << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
  if( vm.count( "all-interpolators" ) || vm.count( "svp-Ndof" ) )
    succ_reports << "SVP_Ndof interp"
                 << "\n"
                 << "\t Successes: " << svp_Ndof_succ_count << "\n"
                 << "\t Graceful Failures: " << svp_Ndof_graceful_fails << "\n"
                 << "\t Total num. of trials: " << mc_runs << "\n"
                 << "\t Success-rate: " << ( 100.0 * double( svp_Ndof_succ_count ) / double( mc_runs ) ) << "\%\n"
                 << "\t Graceful-fail-rate: " << ( 100.0 * double( svp_Ndof_graceful_fails ) / double( mc_runs ) )
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << ( 100.0 * ( 1.0 - double( svp_Ndof_succ_count + svp_Ndof_graceful_fails ) / double( mc_runs ) ) )
                 << "\%." << std::endl;
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
  if( vm.count( "all-interpolators" ) || vm.count( "sap-Ndof" ) )
    succ_reports << "SAP_Ndof interp"
                 << "\n"
                 << "\t Successes: " << sap_Ndof_succ_count << "\n"
                 << "\t Graceful Failures: " << sap_Ndof_graceful_fails << "\n"
                 << "\t Total num. of trials: " << mc_runs << "\n"
                 << "\t Success-rate: " << ( 100.0 * double( sap_Ndof_succ_count ) / double( mc_runs ) ) << "\%\n"
                 << "\t Graceful-fail-rate: " << ( 100.0 * double( sap_Ndof_graceful_fails ) / double( mc_runs ) )
                 << "\%\n"
                 << "\t Hard-fail-rate: "
                 << ( 100.0 * ( 1.0 - double( sap_Ndof_succ_count + sap_Ndof_graceful_fails ) / double( mc_runs ) ) )
                 << "\%." << std::endl;
#endif

  succ_reports.close();
};


int main( int argc, char** argv ) {

  std::srand( static_cast< unsigned int >( std::time( NULL ) ) );

  using namespace ReaK;


  po::options_description generic_options( "Generic options" );
  generic_options.add_options()( "help,h", "produce this help message." );

  po::options_description io_options( "I/O options" );
  io_options.add_options()( "output-path,o", po::value< std::string >()->default_value( "test_interp_results" ),
                            "specify the output path (default is test_interp_results)" );

  po::options_description mc_options( "Monte-Carlo options" );
  mc_options.add_options()( "mc-runs", po::value< std::size_t >()->default_value( 100 ),
                            "number of monte-carlo runs to perform (default is 100)" );

  po::options_description space_options( "Monte-Carlo options" );
  mc_options.add_options()( "space-dimensionality", po::value< std::size_t >()->default_value( 1 ),
                            "number of dimensions for the underlying space (default is 1)" );

  po::options_description interp_select_options( "Interpolator selection options" );
  interp_select_options.add_options()(
    "all-interpolators,a",
    "specify that all supported interpolators should be run (default if no particular interpolator is specified)" )
#ifdef RK_ENABLE_TEST_LINEAR_INTERPOLATOR
    ( "linear", "specify that the uni-directional RRT algorithm should be run" )
#endif
#ifdef RK_ENABLE_TEST_CUBIC_INTERPOLATOR
    ( "cubic", "specify that the bi-directional RRT algorithm should be run" )
#endif
#ifdef RK_ENABLE_TEST_QUINTIC_INTERPOLATOR
    ( "quintic", "specify that the RRT* algorithm should be run" )
#endif
#ifdef RK_ENABLE_TEST_SVP_NDOF_INTERPOLATOR
    ( "svp-Ndof", "specify that the PRM algorithm should be run" )
#endif
#ifdef RK_ENABLE_TEST_SAP_NDOF_INTERPOLATOR
    ( "sap-Ndof", "specify that the FADPRM algorithm should be run" )
#endif
    ;

  po::options_description cmdline_options;
  cmdline_options.add( generic_options ).add( io_options ).add( mc_options ).add( space_options ).add(
    interp_select_options );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, cmdline_options ), vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };

  std::size_t sp_dim = vm["space-dimensionality"].as< std::size_t >();


  switch( sp_dim ) {
    case 1:
      perform_mc_tests< 1 >( vm, sp_dim );
      break;
    case 2:
      perform_mc_tests< 2 >( vm, sp_dim );
      break;
    case 3:
      perform_mc_tests< 3 >( vm, sp_dim );
      break;
    case 4:
      perform_mc_tests< 4 >( vm, sp_dim );
      break;
    case 5:
      perform_mc_tests< 5 >( vm, sp_dim );
      break;
    case 6:
      perform_mc_tests< 6 >( vm, sp_dim );
      break;
    case 7:
      perform_mc_tests< 7 >( vm, sp_dim );
      break;
    case 8:
      perform_mc_tests< 8 >( vm, sp_dim );
      break;
    case 9:
      perform_mc_tests< 9 >( vm, sp_dim );
      break;
    default:
      perform_mc_tests< 0 >( vm, sp_dim );
      break;
  };


  return 0;
};
