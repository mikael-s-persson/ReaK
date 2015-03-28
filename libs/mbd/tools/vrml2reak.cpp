/**
 * \file vrml2reak.cpp
 *
 * This application converts a VRML / X3D / OpenInventor 3D model into any serialized
 * format of ReaK geometries.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2013
 */


#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/mbd/coin3D/oi_reader.hpp>

#include <ReaK/core/serialization/xml_archiver.hpp>
#include <ReaK/core/serialization/bin_archiver.hpp>
#include <ReaK/core/serialization/protobuf_archiver.hpp>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace po = boost::program_options;

int main( int argc, char** argv ) {

  using namespace ReaK;
  using namespace geom;


  po::options_description generic_options( "Generic options" );
  generic_options.add_options()( "help,h", "produce this help message." );

  po::options_description io_options( "I/O options" );
  io_options.add_options()( "input-file,i", po::value< std::string >(), "specify the input file (default is stdin)" )(
    "output-file,o", po::value< std::string >(), "specify the output file (default is stdout)" );

  po::options_description format_options( "Format options" );
  format_options.add_options()(
    "xml,x", "output an XML serialization format for the ReaK geometry (this is the default format)" )(
    "bin,b", "output a binary serialization format for the ReaK geometry" )(
    "proto-buf,p", "output a proto-buf serialization format for the ReaK geometry" );

  po::options_description geom_options( "Geometry options" );
  geom_options.add_options()( "geometry,g",
                              "output the colored model of the geometry (for rendering) (this is the default output)" )(
    "proxy-query,q", "output the proximity model of the geometry (for proximity queries)" );

  po::options_description cmdline_options;
  cmdline_options.add( generic_options ).add( io_options ).add( format_options ).add( geom_options );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, cmdline_options ), vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };

  shared_ptr< std::istream > file_source( &std::cin, null_deleter() );
  if( vm.count( "input-file" ) ) {
    file_source = shared_ptr< std::istream >( new std::ifstream( vm["input-file"].as< std::string >().c_str() ) );
    if( file_source->fail() ) {
      std::cout << "Fatal Error: Input file couldn't not be opened!" << std::endl;
      return 2;
    };
  } else {
    shared_ptr< std::stringstream > tmp_ss( new std::stringstream() );
    ( *tmp_ss ) << std::cin.rdbuf();
    file_source = tmp_ss;
  };

  shared_ptr< std::ostream > file_dest( &std::cout, null_deleter() );
  if( vm.count( "output-file" ) ) {
    file_dest = shared_ptr< std::ostream >( new std::ofstream( vm["output-file"].as< std::string >().c_str() ) );
    if( file_dest->fail() ) {
      std::cout << "Fatal Error: Output file couldn't not be opened!" << std::endl;
      return 3;
    };
  };

  if( vm.count( "bin" ) + vm.count( "proto-buf" ) + vm.count( "xml" ) > 1 ) {
    std::cout << "Fatal Error: More than one output format was specified!" << std::endl;
    return 4;
  };

  shared_ptr< colored_model_3D > geom_model;
  shared_ptr< proxy_query_model_3D > proxy_model;
  if( vm.count( "geometry" ) ) {
    geom_model = shared_ptr< colored_model_3D >( new colored_model_3D( "geometric_model" ) );
    if( vm.count( "proxy-query" ) )
      proxy_model = shared_ptr< proxy_query_model_3D >( new proxy_query_model_3D( "proximity_query_model" ) );
  } else {
    if( vm.count( "proxy-query" ) )
      proxy_model = shared_ptr< proxy_query_model_3D >( new proxy_query_model_3D( "proximity_query_model" ) );
    else // must still create a geometric-model (default):
      geom_model = shared_ptr< colored_model_3D >( new colored_model_3D( "geometric_model" ) );
  };

  oi_reader sg_in( *file_source );

  if( geom_model && proxy_model ) {
    sg_in.translate_into( *geom_model, *proxy_model );
  } else if( geom_model ) {
    sg_in >> *geom_model;
  } else if( proxy_model ) {
    sg_in >> *proxy_model;
  };

  shared_ptr< serialization::oarchive > serial_dest;
  if( vm.count( "bin" ) ) {
    serial_dest = shared_ptr< serialization::oarchive >( new serialization::bin_oarchive( *file_dest ) );
  } else if( vm.count( "proto-buf" ) ) {
    serial_dest = shared_ptr< serialization::oarchive >( new serialization::protobuf_oarchive( *file_dest ) );
  } else {
    serial_dest = shared_ptr< serialization::oarchive >( new serialization::xml_oarchive( *file_dest ) );
  };

  if( geom_model )
    ( *serial_dest ) << geom_model;
  if( proxy_model )
    ( *serial_dest ) << proxy_model;

  return 0;
};
