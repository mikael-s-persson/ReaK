
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


#include <ReaK/mbd/models/manip_P3R3R_arm.hpp>

#include <ReaK/mbd/models/joint_space_limits.hpp>

#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/shapes/sphere.hpp>
#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/core/recorders/data_record_po.hpp>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


int main( int argc, char** argv ) {

  using namespace ReaK;
  using namespace geom;
  using namespace kte;


  po::options_description generic_options( "Generic options" );
  generic_options.add_options()( "help,h", "produce this help message." );

  po::options_description mdl_options( "Model I/O options" );
  mdl_options.add_options()(
    "CRS-kin-model", po::value< std::string >()->default_value( "models/CRS_A465.model.rkx" ),
    "specify the input file for the kinematic model of the CRS robot (default is 'models/CRS_A465.model.rkx')" );

  po::options_description io_options = recorder::get_data_stream_options_po_desc( true, true );


  po::options_description cmdline_options;
  cmdline_options.add( generic_options ).add( mdl_options ).add( io_options );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, cmdline_options ), vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };

  std::string CRS_mdl_fname = vm["CRS-kin-model"].as< std::string >();

  shared_ptr< frame_3D< double > > base_frame;
  shared_ptr< manip_P3R3R_kinematics > CRS_kte_model;
  shared_ptr< joint_limits_collection< double > > joint_rate_limits;

  try {
    shared_ptr< colored_model_3D > mdl_geom;
    shared_ptr< proxy_query_model_3D > mdl_prox;
    ( *serialization::open_iarchive( CRS_mdl_fname ) ) >> base_frame >> CRS_kte_model >> joint_rate_limits >> mdl_geom
      >> mdl_prox;
  } catch( std::exception& e ) {
    std::cerr << "An error occurred while trying to load the CRS robot's kinematics model! Got exception: '" << e.what()
              << "'." << std::endl;
    return 2;
  };


  recorder::data_stream_options data_in_opt = recorder::get_data_stream_options_from_po( vm );

  shared_ptr< recorder::data_extractor > data_in;
  std::vector< std::string > data_in_names;
  try {
    std::tie( data_in, data_in_names ) = data_in_opt.create_extractor();
  } catch( std::exception& e ) {
    std::cerr << "An error occurred while trying to open the input data stream! Got exception: '" << e.what() << "'."
              << std::endl;
    return 3;
  };

  recorder::named_value_row nvr_in = data_in->getFreshNamedValueRow();

  try {
    nvr_in["time"];
    nvr_in["q_0"];
    nvr_in["q_1"];
    nvr_in["q_2"];
    nvr_in["q_3"];
    nvr_in["q_4"];
    nvr_in["q_5"];
    nvr_in["q_6"];
  } catch( recorder::out_of_bounds& e ) {
    RK_UNUSED( e );
    std::cerr << "Could not recognize the input data fields!" << std::endl;
    return 4;
  };

  recorder::data_stream_options data_out_opt = recorder::get_data_stream_options_from_po( vm, true );
  data_out_opt.names.clear();
  data_out_opt.add_name( "time" )
    .add_name( "pos_x" )
    .add_name( "pos_y" )
    .add_name( "pos_z" )
    .add_name( "q0" )
    .add_name( "q1" )
    .add_name( "q2" )
    .add_name( "q3" )
    .add_name( "vel_x" )
    .add_name( "vel_y" )
    .add_name( "vel_z" )
    .add_name( "avel_x" )
    .add_name( "avel_y" )
    .add_name( "avel_z" );
  shared_ptr< recorder::data_recorder > data_out = data_out_opt.create_recorder();

  try {
    double prev_time = 0.0;
    vect_n< double > prev_jtctrl( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    while( true ) {
      ( *data_in ) >> nvr_in;

      vect_n< double > cur_jtctrl( nvr_in["q_0"], nvr_in["q_1"], nvr_in["q_2"], nvr_in["q_3"], nvr_in["q_4"],
                                   nvr_in["q_5"], nvr_in["q_6"] );
      CRS_kte_model->setJointPositions( cur_jtctrl );
      if( prev_time > 1e-6 ) {
        double dt = nvr_in["time"] - prev_time;
        CRS_kte_model->setJointVelocities( ( cur_jtctrl - prev_jtctrl ) * ( 1.0 / dt ) );
      } else {
        CRS_kte_model->setJointVelocities( vect_n< double >( 7, 0.0 ) );
      };
      CRS_kte_model->doDirectMotion();
      frame_3D< double > cur_EE = CRS_kte_model->getDependentFrame3D( 0 )->mFrame->getGlobalFrame();

      ( *data_out ) << nvr_in["time"] << cur_EE.Position << cur_EE.Quat[0] << cur_EE.Quat[1] << cur_EE.Quat[2]
                    << cur_EE.Quat[3] << cur_EE.Velocity << cur_EE.AngVelocity;
      ( *data_out ) << recorder::data_recorder::end_value_row;
    };
  } catch( recorder::end_of_record& e ) {
    RK_UNUSED( e );
  };

  ( *data_out ) << recorder::data_recorder::flush;

  return 0;
};
