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

/**
 * This function constructs a Boost.Program-Options descriptor for data-streaming options.
 * This function can either construct options for input, output, or both.
 */
boost::program_options::options_description get_data_stream_options_po_desc(bool aInput, bool aOutput = false) {
  using boost::program_options::options_description;
  using boost::program_options::value;
  
  options_description result;
  
  if( aInput ) {
    options_description input_options("Data-stream Input Options");
    input_options.add_options()
      ("input",         value< std::string >(), "specify the filename for the input data-stream")
      ("input-ip",      value< std::string >(), "specify the IP-address of the data-server")
      ("input-port",    value< unsigned int >()->default_value(17017), "specify the IP-port to connect to the data-server (default is 17017)")
      ("input-tcp",     "if set, will try to listen to an input TCP data-stream")
      ("input-udp",     "if set, will try to listen to an input UDP data-stream")
      ("input-raw-udp", "if set, will try to listen to an input RAW UDP data-stream, for this to work, you must specify the list of columns via the 'keep-columns' option")
      ("input-format",  value< std::string >(), "specify the format for the input file (default is to use the file-extension of input-file (ssv, tsv, csv, bin, etc.), or if piped, use 'ssv')")
      ("input-flush-freq",  value< unsigned int >()->default_value(50), "specify the flushing frequency of the input datastream (default is 50Hz, if set to 0, input will be immediate (as soon as available))")
      ("input-buffer-size", value< unsigned int >()->default_value(500), "specify the desired buffer size of the input datastream, without guarantee that the buffer will be limited to that amount (default is 500 bytes, if set to 0, input will not be buffered, i.e., read as you go)")
      ("input-immediate-mode", "if set, will make the input datastream without buffering or regular flushing, i.e., the operations are immediate (note: this option overrides the flush-freq and buffer-size options)")
    ;
    result.add(input_options);
  };
  if( aOutput ) {
    options_description output_options("Data-stream Output Options");
    output_options.add_options()
      ("output",        value< std::string >(), "specify the filename for the output file for the datastream (default is to output to the console 'stdout', for piping the output)")
      ("output-ip",     value< std::string >(), "specify the IP-address of the RAW UDP data-stream, this is because a raw udp stream is connection-less and requires a pre-defined output IP-address")
      ("output-port",   value< unsigned int >()->default_value(17017), "specify the IP-port for the data-server that will be created (default is 17017)")
      ("output-tcp",    "if set, will output a TCP data-stream")
      ("output-udp",    "if set, will output a UDP data-stream")
      ("output-raw-udp","if set, will output a RAW UDP data-stream")
      ("output-format", value< std::string >(), "specify the format for the output file (default is to use the file-extension of input-file (ssv, tsv, csv, bin, etc.), or if piped, use 'ssv')")
      ("output-flush-freq",  value< unsigned int >()->default_value(50), "specify the flushing frequency of the output datastream (default is 50Hz, if set to 0, output will be immediate)")
      ("output-buffer-size", value< unsigned int >()->default_value(500), "specify the desired buffer size of the output datastream, without guarantee that the buffer will be limited to that amount (default is 500 bytes, if set to 0, output will not be buffered, i.e., sent as you go)")
      ("output-immediate-mode", "if set, will make the output datastream without buffering or regular flushing, i.e., the operations are immediate (note: this option overrides the flush-freq and buffer-size options)")
    ;
    result.add(output_options);
  };
  
  options_description filt_options("Data-stream Filtering Options");
  filt_options.add_options()
    ("time-column-sync", value< std::string >(), "name of the time-column from the datastream to use as a timed-output, i.e., for a network streaming, this will deliver the rows at each correct time")
    ("keep-columns",     value< std::string >(), "specify a semi-colon-separated list of the columns to output (keep) in the datastream")
  ;
  
  result.add(filt_options);
  
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
      result.file_name = vm["input"].as<std::string>();
      std::size_t p = result.file_name.find_last_of('.');
      input_extension = result.file_name.substr(p+1);
    } else {
      result.file_name = "stdin";
    };
    if(vm.count("input-format"))
      input_extension = vm["input-format"].as<std::string>();
    
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
      } else if(input_extension == "csv") {
        result.kind = data_stream_options::comma_separated;
      } else if(input_extension == "bin") {
        result.kind = data_stream_options::binary;
      };
    };
    
    if(vm.count("input-immediate-mode")) {
      result.set_unbuffered();
    } else {
      result.flush_rate = vm["input-flush-freq"].as<unsigned int>();
      result.buffer_size = vm["input-buffer-size"].as<unsigned int>();      
    };
    
  } else {
    // load options for output.
    
    // take care of 'file_name' and 'kind':
    if(vm.count("output-tcp") + vm.count("output-udp") + vm.count("output-raw-udp") == 0) {
      std::string output_extension = "ssv";
      if(vm.count("output")) {
        result.file_name = vm["output"].as<std::string>();
        
        boost::filesystem::create_directories(boost::filesystem::path(result.file_name).parent_path());
        
        std::size_t p_dot = result.file_name.find_last_of('.');
        output_extension = result.file_name.substr(p_dot + 1);
        
      } else {
        result.file_name = "stdout";
      };
      
      if(vm.count("output-format"))
        output_extension = vm["output-format"].as< std::string >();
      
      if(output_extension == "ssv") {
        result.kind = data_stream_options::space_separated;
      } else if(output_extension == "tsv") {
        result.kind = data_stream_options::tab_separated;
      } else if(output_extension == "csv") {
        result.kind = data_stream_options::comma_separated;
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
    
    if(vm.count("output-immediate-mode")) {
      result.set_unbuffered();
    } else {
      result.flush_rate = vm["output-flush-freq"].as<unsigned int>();
      result.buffer_size = vm["output-buffer-size"].as<unsigned int>();      
    };
    
  };
  
  if(vm.count("keep-columns")) {
    std::string tmp_keep_columns = vm["keep-columns"].as<std::string>();
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
    result.time_sync_name = vm["time-column-sync"].as<std::string>();
  
  return result;
};


};

};

#endif // RK_DATA_RECORD_PO_HPP












