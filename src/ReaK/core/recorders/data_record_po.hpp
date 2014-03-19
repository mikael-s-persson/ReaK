/**
 * \file data_record_po.hpp
 *
 * This library declares utility functions for creating and dealing with program-options related 
 * to a data recording or extraction stream (see data_record.hpp). Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2014
 */

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

#ifndef REAK_DATA_RECORD_PO_HPP
#define REAK_DATA_RECORD_PO_HPP

#include "data_record.hpp"
#include "data_record_options.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Data Recorders and Extractors */
namespace recorder {

namespace detail {
  
  static std::string strip_quotes(const std::string& s) {
    std::size_t first = 0;
    while( ( first < s.length() ) && ( std::isspace(s[first]) || (s[first] == '"') ) )
      ++first;
    std::size_t last = s.length();
    while( ( last > 0 ) && ( std::isspace(s[last-1]) || (s[last-1] == '"') ) )
      --last;
    return s.substr(first, last - first);
  };
  
};

/**
 * This function constructs a Boost.Program-Options descriptor for data-streaming options.
 * This function can either construct options for input, output, or both.
 */
boost::program_options::options_description get_data_stream_options_po_desc(bool aInput, bool aOutput = false) {
  
  
  po::options_description io_options("I/O options");
  if( aInput ) {
    io_options.add_options()
      ("input",         po::value< std::string >(), "specify the filename for the input data-stream")
      ("input-ip",      po::value< std::string >(), "specify the IP-address of the data-server")
      ("input-port",    po::value< unsigned int >()->default_value(17017), "specify the IP-port to connect to the data-server (default is 17017)")
      ("input-tcp",     "if set, will try to listen to an input TCP data-stream")
      ("input-udp",     "if set, will try to listen to an input UDP data-stream")
      ("input-raw-udp", "if set, will try to listen to an input RAW UDP data-stream, for this to work, you must specify the list of columns via the 'keep-columns' option")
      ("input-format",  po::value< std::string >(), "specify the format for the input file (default is to use the file-extension of input-file (ssv, tsv, bin, etc.), or if piped, use 'ssv')")
    ;
  };
  if( aOutput ) {
    io_options.add_options()
      ("output",        po::value< std::string >(), "specify the filename for the output file for the datastream (default is to output to the console 'stdout', for piping the output)")
      ("output-ip",     po::value< std::string >(), "specify the IP-address of the RAW UDP data-stream, this is because a raw udp stream is connection-less and requires a pre-defined output IP-address")
      ("output-port",   po::value< unsigned int >()->default_value(17017), "specify the IP-port for the data-server that will be created (default is 17017)")
      ("output-tcp",    "if set, will output a TCP data-stream")
      ("output-udp",    "if set, will output a UDP data-stream")
      ("output-raw-udp","if set, will output a RAW UDP data-stream")
      ("output-format", po::value< std::string >(), "specify the format for the output file (default is to use the file-extension of input-file (ssv, tsv, bin, etc.), or if piped, use 'ssv')")
    ;
  };
  
  po::options_description filt_options("Filtering options");
  filt_options.add_options()
    ("time-column-sync", po::value< std::string >(), "name of the time-column from the datastream to use as a timed-output, i.e., for a network streaming, this will deliver the rows at each correct time")
    ("keep-columns",     po::value< std::string >(), "specify a semi-colon-separated list of the columns to output (keep) in the datastream")
  ;
  
  po::options_description result;
  result.add(io_options).add(filt_options);
  
  return result;
};


/**
 * This function constructs a data-streaming options object from a Boost.Program-Options variable-map.
 * This function can either construct options for input or output.
 */
data_stream_options get_data_stream_options_from_po(boost::program_options::variables_map& vm, bool aForOutput = false) {
  
  data_stream_options result;
  
  if( !aForOutput ) {
    // load options for input.
    
    std::string input_extension = "ssv";
    if(vm.count("input")) {
      result.file_name = detail::strip_quotes(vm["input"].as<std::string>());
      std::size_t p = result.file_name.find_last_of('.');
      input_extension = result.file_name.substr(p+1);
    } else {
      result.file_name = "stdin";
    };
    if(vm.count("input-format"))
      input_extension = detail::strip_quotes(vm["input-format"].as<std::string>());
    
    if( vm.count("input-tcp") + vm.count("input-udp") + vm.count("input-raw-udp") > 0 ) {
      std::stringstream ss;
      ss << vm["input-ip"].as<std::string>() << ":" << vm["input-port"].as<unsigned int>();
      result.file_name = ss.str();
      if( vm.count("input-udp") ) {
        result.kind = data_stream_options::udp_stream;
      } else if( vm.count("input-raw-udp") ) {
        result.kind = data_stream_options::raw_udp_stream;
      } else {
        result.kind = data_stream_options::tcp_stream;
      };
    } else {
      if(input_extension == "ssv") {
        result.kind = data_stream_options::space_separated;
      } else if(input_extension == "tsv") {
        result.kind = data_stream_options::tab_separated;
      } else if(input_extension == "bin") {
        result.kind = data_stream_options::binary;
      };
    };
    
  } else {
    // load options for output.
    
    std::string output_extension = "ssv";
    if(vm.count("output-format"))
      output_extension = detail::strip_quotes( vm["output-format"].as< std::string >() );
    
    // take care of 'file_name' and 'kind':
    if(vm.count("output-tcp") + vm.count("output-udp") + vm.count("output-raw-udp") == 0) {
      if(vm.count("output")) {
        result.file_name = detail::strip_quotes(vm["output"].as<std::string>());
        
        std::string output_path_name = result.file_name;
        std::size_t p = output_path_name.find_last_of('/');
        if(p == std::string::npos)
          output_path_name = "";
        else
          output_path_name.erase(p);
        while(output_path_name[output_path_name.length()-1] == '/') 
          output_path_name.erase(output_path_name.length()-1, 1);
        if(!output_path_name.empty())
          boost::filesystem::create_directory(output_path_name.c_str());
        
        std::size_t p_dot = output_file_name.find_last_of('.');
        output_extension = output_file_name.substr(p_dot + 1);
        
      } else {
        result.file_name = "stdout";
      };
      
      if(output_extension == "ssv") {
        result.kind = data_stream_options::space_separated;
      } else if(output_extension == "tsv") {
        result.kind = data_stream_options::tab_separated;
      } else if(output_extension == "bin") {
        result.kind = data_stream_options::binary;
      };
    } else {
      std::stringstream ss;
      if( vm.count("output-raw-udp") ) {
        ss << vm["output-ip"].as<std::string>() << ":" << vm["output-port"].as<unsigned int>();
        result.kind = data_stream_options::raw_udp_stream;
      } else {
        ss << vm["output-port"].as<unsigned int>();
        if( vm.count("output-udp") )
          result.kind = data_stream_options::udp_stream;
        else
          result.kind = data_stream_options::tcp_stream;
      };
      result.file_name = ss.str();
    };
    
  };
  
  if(vm.count("keep-columns")) {
    std::string tmp_keep_columns = detail::strip_quotes(vm["keep-columns"].as<std::string>());
    for(std::string::iterator it = tmp_keep_columns.begin(); it != tmp_keep_columns.end(); ++it) {
      std::string::iterator it_next = std::find(it, tmp_keep_columns.end(), ';');
      std::string new_name(it, it_next);
      result.names.push_back(new_name);
      if(it_next == tmp_keep_columns.end())
        break;
      it = it_next;
    };
  };
  
  if(vm.count("time-column-sync"))
    result.time_sync_name = detail::strip_quotes(vm["time-column-sync"].as<std::string>());
  
  return result;
};


};

};

#endif // RK_DATA_RECORD_PO_HPP












