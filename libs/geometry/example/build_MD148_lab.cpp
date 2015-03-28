/**
 * \file build_MD148_lab.cpp
 *
 * This application constructs a geometric model of the relevant objects in the MD148 laboratory
 * environment for the experimental tests with the CRS robot and the airship.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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


#include <ReaK/geometry/shapes/plane.hpp>
#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>
#include <ReaK/geometry/shapes/capped_cylinder.hpp>
#include <ReaK/geometry/shapes/colored_model.hpp>

#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


int main( int argc, char** argv ) {

  po::options_description generic_options( "Generic options" );
  generic_options.add_options()( "help,h", "produce this help message." );

  po::options_description io_options( "I/O options" );
  io_options.add_options()( "output-path,p", po::value< std::string >()->default_value( "models" ),
                            "specify the output path (default is 'models')" )(
    "output-name,o", po::value< std::string >()->default_value( "MD148_lab" ),
    "specify the output base-name (default is 'MD148_lab')" )(
    "format", po::value< std::string >()->default_value( "xml" ),
    "specify the format that should be outputted (default is 'xml', but can also be 'bin' or 'protobuf')" );

  po::options_description cmdline_options;
  cmdline_options.add( generic_options ).add( io_options );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, cmdline_options ), vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };

  std::string output_base_name = vm["output-name"].as< std::string >();

  std::string output_path_name = vm["output-path"].as< std::string >();
  while( output_path_name[output_path_name.length() - 1] == '/' )
    output_path_name.erase( output_path_name.length() - 1, 1 );

  fs::create_directory( output_path_name.c_str() );

  std::string output_extension = ".rkx";
  if( vm["format"].as< std::string >() == "bin" ) {
    output_extension = ".rkb";
  } else if( vm["format"].as< std::string >() == "protobuf" ) {
    output_extension = ".pbuf";
  };


  using namespace ReaK;
  using namespace geom;

  shared_ptr< colored_model_3D > MD148_basic_lab( new colored_model_3D( "MD148_basic_lab_render" ) );
  shared_ptr< proxy_query_model_3D > MD148_lab_proxy( new proxy_query_model_3D( "MD148_basic_lab_proxy" ) );


  shared_ptr< plane > lab_floor( new plane(
    "MD148_floor", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( -0.8, -1.0, 0.0 ), quaternion< double >() ),
    vect< double, 2 >( 4.0, 6.0 ) ) );

  shared_ptr< plane > lab_n_wall( new plane(
    "MD148_north_wall", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( 1.2, -1.0, 1.5 ),
                       axis_angle< double >( M_PI * 0.5, vect< double, 3 >( 0.0, -1.0, 0.0 ) ).getQuaternion() ),
    vect< double, 2 >( 3.0, 6.0 ) ) );

  shared_ptr< plane > lab_w_wall( new plane(
    "MD148_west_wall", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( -0.8, 2.0, 1.5 ),
                       axis_angle< double >( M_PI * 0.5, vect< double, 3 >( 1.0, 0.0, 0.0 ) ).getQuaternion() ),
    vect< double, 2 >( 4.0, 3.0 ) ) );

  shared_ptr< box > lab_robot_track( new box(
    "MD148_robot_track", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( 0.0, -1.71, 0.15 ), quaternion< double >() ),
    vect< double, 3 >( 0.4, 3.42, 0.3 ) ) );

  shared_ptr< plane > lab_operator_wall( new plane(
    "MD148_operator_wall", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( -0.8, -3.5, 1.5 ),
                       axis_angle< double >( M_PI * 0.5, vect< double, 3 >( -1.0, 0.0, 0.0 ) ).getQuaternion() ),
    vect< double, 2 >( 4.0, 3.0 ) ) );

  shared_ptr< capped_cylinder > lab_robot_track_left( new capped_cylinder(
    "MD148_robot_track_left", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( 0.1, -1.71, 0.15 ),
                       axis_angle< double >( M_PI * 0.5, vect< double, 3 >( 1.0, 0.0, 0.0 ) ).getQuaternion() ),
    3.42, 0.18 ) );

  shared_ptr< capped_cylinder > lab_robot_track_right( new capped_cylinder(
    "MD148_robot_track_right", shared_ptr< pose_3D< double > >(),
    pose_3D< double >( weak_ptr< pose_3D< double > >(), vect< double, 3 >( -0.1, -1.71, 0.15 ),
                       axis_angle< double >( M_PI * 0.5, vect< double, 3 >( 1.0, 0.0, 0.0 ) ).getQuaternion() ),
    3.42, 0.18 ) );

  shared_ptr< coord_arrows_3D > lab_global_arrows(
    new coord_arrows_3D( "global_frame_arrows", shared_ptr< pose_3D< double > >(), pose_3D< double >(), 1.0 ) );

  ( *MD148_basic_lab )
    .addElement( color( 0, 0, 0 ), lab_global_arrows )
    .addElement( color( 0.3, 0.3, 0.3 ), lab_floor )
    //     .addElement(color(0.8,0.8,0.8),lab_n_wall)
    //     .addElement(color(0.8,0.8,0.8),lab_w_wall)
    //     .addElement(color(0.8,0.8,0.8),lab_operator_wall)
    .addElement( color( 0.35, 0.35, 0.35 ), lab_robot_track );

  ( *MD148_lab_proxy )
    .addShape( lab_floor )
    .addShape( lab_n_wall )
    .addShape( lab_w_wall )
    .addShape( lab_operator_wall )
    .addShape( lab_robot_track_left )
    .addShape( lab_robot_track_right );

  ( *serialization::open_oarchive( output_path_name + "/" + output_base_name + ".geom" + output_extension ) )
    << MD148_basic_lab << MD148_lab_proxy;


  return 0;
};
