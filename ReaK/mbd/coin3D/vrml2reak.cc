/**
 * \file vrml2reak.cpp
 *
 * This application converts a VRML / X3D / OpenInventor 3D model into any serialized
 * format of ReaK geometries.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2013
 */

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/geometry/shapes/colored_model.h"

#include "ReaK/mbd/coin3D/oi_reader.h"

#include "ReaK/core/serialization/bin_archiver.h"
#include "ReaK/core/serialization/protobuf_archiver.h"
#include "ReaK/core/serialization/xml_archiver.h"

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <fstream>
#include <iostream>
#include <sstream>

// I/O options
ABSL_FLAG(std::string, input_file, "stdin",
          "Specify the input file (default is stdin).");
ABSL_FLAG(std::string, output_file, "stdout",
          "Specify the output file (default is stdout).");

// Format options
ABSL_FLAG(bool, xml, false,
          "Output an XML serialization format for the ReaK geometry (this is "
          "the default format).");
ABSL_FLAG(bool, bin, false,
          "Output a binary serialization format for the ReaK geometry.");
ABSL_FLAG(bool, protobuf, false,
          "Output a protobuf serialization format for the ReaK geometry.");

// Geometry options
ABSL_FLAG(bool, geometry, false,
          "Output the colored model of the geometry (for rendering) (this is "
          "the default output).");
ABSL_FLAG(
    bool, proxy_query, false,
    "Output the proximity model of the geometry (for proximity queries).");

int main(int argc, char** argv) {

  using namespace ReaK;
  using namespace geom;

  absl::ParseCommandLine(argc, argv);

  std::shared_ptr<std::istream> file_source(&std::cin, null_deleter());
  if (absl::GetFlag(FLAGS_input_file) != "stdin") {
    file_source = std::make_shared<std::ifstream>(
        absl::GetFlag(FLAGS_input_file).c_str());
    if (file_source->fail()) {
      std::cout << "Fatal Error: Input file couldn't not be opened!"
                << std::endl;
      return 2;
    }
  } else {
    std::shared_ptr<std::stringstream> tmp_ss(new std::stringstream());
    (*tmp_ss) << std::cin.rdbuf();
    file_source = tmp_ss;
  }

  std::shared_ptr<std::ostream> file_dest(&std::cout, null_deleter());
  if (absl::GetFlag(FLAGS_output_file) != "stdout") {
    file_dest = std::make_shared<std::ofstream>(
        absl::GetFlag(FLAGS_output_file).c_str());
    if (file_dest->fail()) {
      std::cout << "Fatal Error: Output file couldn't not be opened!"
                << std::endl;
      return 3;
    }
  }

  if (static_cast<int>(absl::GetFlag(FLAGS_bin)) +
          static_cast<int>(absl::GetFlag(FLAGS_protobuf)) +
          static_cast<int>(absl::GetFlag(FLAGS_xml)) >
      1) {
    std::cout << "Fatal Error: More than one output format was specified!"
              << std::endl;
    return 4;
  }

  std::shared_ptr<colored_model_3D> geom_model;
  std::shared_ptr<proxy_query_model_3D> proxy_model;
  if (absl::GetFlag(FLAGS_geometry)) {
    geom_model = std::make_shared<colored_model_3D>("geometric_model");
    if (absl::GetFlag(FLAGS_proxy_query)) {
      proxy_model =
          std::make_shared<proxy_query_model_3D>("proximity_query_model");
    }
  } else {
    if (absl::GetFlag(FLAGS_proxy_query)) {
      proxy_model =
          std::make_shared<proxy_query_model_3D>("proximity_query_model");
    } else {
      geom_model = std::make_shared<colored_model_3D>("geometric_model");
    }
  };

  oi_reader sg_in(*file_source);

  if (geom_model && proxy_model) {
    sg_in.translate_into(*geom_model, *proxy_model);
  } else if (geom_model) {
    sg_in >> *geom_model;
  } else if (proxy_model) {
    sg_in >> *proxy_model;
  }

  std::shared_ptr<serialization::oarchive> serial_dest;
  if (absl::GetFlag(FLAGS_bin)) {
    serial_dest = std::make_shared<serialization::bin_oarchive>(*file_dest);
  } else if (absl::GetFlag(FLAGS_protobuf)) {
    serial_dest =
        std::make_shared<serialization::protobuf_oarchive>(*file_dest);
  } else {
    serial_dest = std::make_shared<serialization::xml_oarchive>(*file_dest);
  }

  if (geom_model) {
    (*serial_dest) << geom_model;
  }
  if (proxy_model) {
    (*serial_dest) << proxy_model;
  }

  return 0;
};
