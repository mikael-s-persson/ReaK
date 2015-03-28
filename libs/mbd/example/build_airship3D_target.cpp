
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


#include <ReaK/mbd/kte/inertia.hpp>
#include <ReaK/mbd/kte/driving_actuator.hpp>
#include <ReaK/mbd/kte/state_measures.hpp>
#include <ReaK/mbd/kte/free_joints.hpp>
#include <ReaK/mbd/kte/kte_map_chain.hpp>

#include <ReaK/mbd/models/manip_dynamics_model.hpp>

#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/shapes/sphere.hpp>
#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>
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
    "output-name,o", po::value< std::string >()->default_value( "airship3D" ),
    "specify the output base-name (default is 'airship3D')" )(
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
  using namespace kte;

  shared_ptr< frame_3D< double > > global_base( new frame_3D< double >() );

  //  Output from calibration program:
  //   (Position = (0.176438; -0.356764; -0.0521845); Quaternion = (0.999986; 0.00301881; -0.00257751; 0.00343278))
  global_base->Position = vect< double, 3 >( 0.176438, -0.356764, -0.0521845 );
  global_base->Quat = quaternion< double >( vect< double, 4 >( 0.999986, 0.00301881, -0.00257751, 0.00343278 ) );


  shared_ptr< frame_3D< double > > airship3D_frame( new frame_3D< double >() );

  shared_ptr< frame_3D< double > > airship3D_output_frame( new frame_3D< double >() );

  shared_ptr< jacobian_3D_3D< double > > airship3D_joint_jac( new jacobian_3D_3D< double >() );

  shared_ptr< free_joint_3D > airship3D_joint(
    new free_joint_3D( "airship3D_joint", airship3D_frame, global_base, airship3D_output_frame, airship3D_joint_jac ) );

  shared_ptr< joint_dependent_frame_3D > airship3D_dep_frame( new joint_dependent_frame_3D( airship3D_output_frame ) );
  airship3D_dep_frame->add_joint( airship3D_frame, airship3D_joint_jac );

  shared_ptr< inertia_3D > airship3D_inertia(
    new inertia_3D( "airship3D_inertia", airship3D_dep_frame, 1.0,
                    mat< double, mat_structure::symmetric >( mat< double, mat_structure::identity >( 3 ) ) ) );

  shared_ptr< driving_actuator_3D > airship3D_actuator(
    new driving_actuator_3D( "airship3D_actuator", airship3D_frame, airship3D_joint ) );

  shared_ptr< position_measure_3D > airship3D_position(
    new position_measure_3D( "airship3D_position", airship3D_frame ) );

  shared_ptr< rotation_measure_3D > airship3D_rotation(
    new rotation_measure_3D( "airship3D_rotation", airship3D_frame ) );

  shared_ptr< kte_map_chain > airship3D_model( new kte_map_chain( "airship3D_model" ) );

  ( *airship3D_model ) << airship3D_position << airship3D_rotation << airship3D_actuator << airship3D_joint
                       << airship3D_inertia;


  shared_ptr< manipulator_dynamics_model > airship3D_dyn_model(
    new manipulator_dynamics_model( "airship3D_dyn_model" ) );
  airship3D_dyn_model->setModel( airship3D_model );
  ( *airship3D_dyn_model ) << airship3D_frame;
  ( *airship3D_dyn_model ) << airship3D_inertia; // also adds the dependent frame 'airship3D_dep_frame'

  ( *airship3D_dyn_model ) << airship3D_actuator;
  ( *airship3D_dyn_model ) << airship3D_position;
  ( *airship3D_dyn_model ) << airship3D_rotation;

  /* Old values (none-sense, I think):
  shared_ptr< frame_3D<double> >
    airship3D_grasp_frame( new frame_3D<double>(airship3D_output_frame,
      vect<double,3>(0.97 * std::sin(0.2 / 0.93),0.0,0.97 * std::cos(0.2 / 0.93)),
      axis_angle<double>(0.2 / 0.93 / 2.0,vect<double,3>(0.0,1.0,0.0)).getQuaternion()
      * quaternion<double>::yrot(M_PI) * quaternion<double>::zrot(0.5 * M_PI),
      vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0),
      vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0),
      vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0)
    ));
  airship3D_grasp_frame->Position += airship3D_grasp_frame->Quat * (-0.3 * vect_k);
  */

  // Original position (near prop A)
  //   double gr_radius = 0.93;
  //   double gr_arc_length = 0.225;

  // New position (lower, during bottom-heavy tests)
  double gr_radius = 0.93;
  double gr_arc_length = 0.584;

  double gr_beta = gr_arc_length / gr_radius;

  shared_ptr< frame_3D< double > > airship3D_grasp_frame( new frame_3D< double >(
    airship3D_output_frame, vect< double, 3 >( gr_radius * std::cos( gr_beta ), 0.0, gr_radius * std::sin( gr_beta ) ),
    quaternion< double >::yrot( -0.5 * M_PI - gr_beta ).getQuaternion(), vect< double, 3 >( 0.0, 0.0, 0.0 ),
    vect< double, 3 >( 0.0, 0.0, 0.0 ), vect< double, 3 >( 0.0, 0.0, 0.0 ), vect< double, 3 >( 0.0, 0.0, 0.0 ),
    vect< double, 3 >( 0.0, 0.0, 0.0 ), vect< double, 3 >( 0.0, 0.0, 0.0 ) ) );
  airship3D_grasp_frame->Position += airship3D_grasp_frame->Quat * ( -0.3 * vect_k );

  shared_ptr< joint_dependent_frame_3D > airship3D_dep_grasp_frame(
    new joint_dependent_frame_3D( airship3D_grasp_frame ) );
  airship3D_dep_grasp_frame->add_joint( airship3D_frame, airship3D_joint_jac );


  shared_ptr< sphere > hull( new sphere( "airship3D_hull", airship3D_output_frame, pose_3D< double >(), 0.93 ) );

  shared_ptr< box > grapple(
    new box( "airship3D_grapple", airship3D_output_frame,
             pose_3D< double >( shared_ptr< pose_3D< double > >(),
                                vect< double, 3 >( ( gr_radius + 0.05 ) * std::cos( gr_beta ), 0.0,
                                                   ( gr_radius + 0.05 ) * std::sin( gr_beta ) ),
                                quaternion< double >::yrot( -0.5 * M_PI - gr_beta ).getQuaternion() ),
             vect< double, 3 >( 0.08, 0.003, 0.1 ) ) );

  shared_ptr< coord_arrows_3D > grapple_arrows(
    new coord_arrows_3D( "airship3D_grapple_arrows", airship3D_grasp_frame, pose_3D< double >(), 0.3 ) );


  shared_ptr< colored_model_3D > geom_mdl
    = shared_ptr< colored_model_3D >( new colored_model_3D( "airship3D_geom_render" ) );
  ( *geom_mdl )
    .addAnchor( airship3D_output_frame )
    .addElement( color( 0, 0, 0 ), grapple_arrows )
    .addElement( color( 1, 1, 1 ), hull )
    .addElement( color( 0.2, 0.2, 0.2 ), grapple );

  shared_ptr< proxy_query_model_3D > proxy_mdl
    = shared_ptr< proxy_query_model_3D >( new proxy_query_model_3D( "airship3D_geom_render" ) );
  ( *proxy_mdl ).addShape( hull );


  shared_ptr< manipulator_kinematics_model > airship3D_kin_model(
    new manipulator_kinematics_model( "airship3D_kin_model" ) );
  airship3D_kin_model->setModel( airship3D_model );
  ( *airship3D_kin_model ) << airship3D_frame;
  ( *airship3D_kin_model ) << airship3D_dep_grasp_frame;


  ( *serialization::open_oarchive( output_path_name + "/" + output_base_name + ".model" + output_extension ) )
    << global_base << airship3D_kin_model << airship3D_grasp_frame << geom_mdl << proxy_mdl;
};
