
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


#include "ssv_recorder.hpp"
#include "tsv_recorder.hpp"
#include "bin_recorder.hpp"
#include "tcp_recorder.hpp"
#include "udp_recorder.hpp"
#include "raw_udp_recorder.hpp"

#include <sstream>

#include "base/chrono_incl.hpp"
#include "base/thread_incl.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;





static std::string strip_quotes(const std::string& s) {
  std::size_t first = 0;
  while( ( first < s.length() ) && ( std::isspace(s[first]) || (s[first] == '"') ) )
    ++first;
  std::size_t last = s.length();
  while( ( last > 0 ) && ( std::isspace(s[last-1]) || (s[last-1] == '"') ) )
    --last;
  return s.substr(first, last - first);
};


namespace ch = ReaKaux::chrono;
typedef ch::steady_clock stc;

int main(int argc, char** argv) {
  using namespace ReaK;
  using namespace recorder;
  
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("input,i",         po::value< std::string >(), "specify the filename for the input data-stream")
    ("input-ip,I",      po::value< std::string >(), "specify the IP-address of the data-server")
    ("input-port,P",    po::value< unsigned int >()->default_value(17017), "specify the IP-port to connect to the data-server (default is 17017)")
    ("input-udp,U",     "if set, will try to listen to an input UDP data-stream")
    ("input-raw-udp,R", "if set, will try to listen to an input RAW UDP data-stream, for this to work, you must specify the list of columns via the 'keep-columns' option")
    ("input-format,F",  po::value< std::string >(), "specify the format for the input file (default is to use the file-extension of input-file (ssv, tsv, bin, etc.), or if piped, use 'ssv')")
    ("output,o",        po::value< std::string >(), "specify the filename for the output file for the datastream (default is to output to the console 'stdout', for piping the output)")
    ("output-port,p",   po::value< unsigned int >()->default_value(17017), "specify the IP-port for the data-server that will be created (default is 17017)")
    ("output-tcp,t",    "if set, will output a TCP data-stream")
    ("output-udp,u",    "if set, will output a UDP data-stream")
    ("output-raw-udp,r","if set, will output a RAW UDP data-stream")
    ("output-ip",       po::value< std::string >(), "specify the IP-address of the RAW UDP data-stream, this is because a raw udp stream is connection-less and requires a pre-defined output IP-address")
    ("output-format,f", po::value< std::string >(), "specify the format for the output file (default is to use the file-extension of input-file (ssv, tsv, bin, etc.), or if piped, use 'ssv')")
  ;
  
  po::options_description filt_options("Filtering options");
  filt_options.add_options()
    ("time-column-sync", po::value< std::string >()->default_value("t"), "name of the time-column from the datastream to use as a timed-output, i.e., for a network streaming, this will deliver the rows at each correct time (default is 't')")
    ("keep-columns",     po::value< std::string >(), "specify a semi-colon-separated list of the columns to output (keep) in the datastream")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(filt_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  
  shared_ptr< data_recorder > data_out;
  std::string output_extension = "ssv";
  std::string output_file_name;
  if(vm.count("output")) {
    // make sure that the output directory exists:
    output_file_name = strip_quotes(vm["output"].as<std::string>());
    std::string output_path_name = output_file_name;
    std::size_t p = output_path_name.find_last_of('/');
    if(p == std::string::npos)
      output_path_name = "";
    else
      output_path_name.erase(p);
    while(output_path_name[output_path_name.length()-1] == '/') 
      output_path_name.erase(output_path_name.length()-1, 1);
    if(!output_path_name.empty())
      fs::create_directory(output_path_name.c_str());
    
    std::size_t p_dot = output_file_name.find_last_of('.');
    output_extension = output_file_name.substr(p_dot + 1);
  };
  
  if(vm.count("output-tcp")) {
    std::stringstream ss;
    ss << vm["output-port"].as<unsigned int>();
    data_out = shared_ptr< data_recorder >(new tcp_recorder(ss.str()));
  } else if(vm.count("output-udp")) {
    std::stringstream ss;
    ss << vm["output-port"].as<unsigned int>();
    data_out = shared_ptr< data_recorder >(new udp_recorder(ss.str()));
  } else if(vm.count("output-raw-udp")) {
    std::stringstream ss;
    ss << vm["output-ip"].as<std::string>() << ":" << vm["output-port"].as<unsigned int>();
    data_out = shared_ptr< data_recorder >(new raw_udp_recorder(ss.str()));
  } else {
    if(vm.count("output-format"))
      output_extension = strip_quotes(vm["output-format"].as<std::string>());
    if(output_extension == "ssv") {
      data_out = shared_ptr< data_recorder >(new ssv_recorder());
    } else if(output_extension == "tsv") {
      data_out = shared_ptr< data_recorder >(new tsv_recorder());
    } else if(output_extension == "bin") {
      data_out = shared_ptr< data_recorder >(new bin_recorder());
    } else {
      std::cerr << "Error! Unsupported output format!" << std::endl;
      return 1;
    };
    
    if(vm.count("output")) {
      data_out->setFileName(output_file_name);
    } else {
      data_out->setStream(std::cout);
    };
  };
  
  shared_ptr< data_extractor > data_in;
  std::string input_extension = "ssv";
  std::string input_file_name;
  if(vm.count("input")) {
    input_file_name = strip_quotes(vm["input"].as<std::string>());
    std::size_t p = input_file_name.find_last_of('.');
    input_extension = input_file_name.substr(p+1);
  };
  
  if(vm.count("input-ip") && !vm.count("input-udp") && !vm.count("input-raw-udp")) {
    std::stringstream ss;
    ss << vm["input-ip"].as<std::string>() << ":" << vm["input-port"].as<unsigned int>();
    data_in = shared_ptr< data_extractor >(new tcp_extractor(ss.str()));
  } else if(vm.count("input-ip") && !vm.count("input-raw-udp")) {
    std::stringstream ss;
    ss << vm["input-ip"].as<std::string>() << ":" << vm["input-port"].as<unsigned int>();
    data_in = shared_ptr< data_extractor >(new udp_extractor(ss.str()));
  } else if(!vm.count("input-raw-udp")) {
    if(vm.count("input-format"))
      input_extension = strip_quotes(vm["input-format"].as<std::string>());
    if(input_extension == "ssv") {
      data_in = shared_ptr< data_extractor >(new ssv_extractor());
    } else if(input_extension == "tsv") {
      data_in = shared_ptr< data_extractor >(new tsv_extractor());
    } else if(input_extension == "bin") {
      data_in = shared_ptr< data_extractor >(new bin_extractor());
    } else {
      std::cerr << "Error! Unsupported input format!" << std::endl;
      return 1;
    };
    
    if(vm.count("input")) {
      data_in->setFileName(input_file_name);
    } else {
      data_in->setStream(std::cin);
    };
  };
  
  
  std::vector<std::string> names_in;
  std::vector<std::string> names_out;
  std::string time_name;
  if(!vm.count("input-raw-udp")) {
    names_in.resize(data_in->getColCount(), "");
    for(std::size_t i = 0; i < names_in.size(); ++i)
      (*data_in) >> names_in[i];
    
    if(vm.count("keep-columns")) {
      std::string tmp_keep_columns = strip_quotes(vm["keep-columns"].as<std::string>());
      for(std::string::iterator it = tmp_keep_columns.begin(); it != tmp_keep_columns.end(); ++it) {
        std::string::iterator it_next = std::find(it, tmp_keep_columns.end(), ';');
        std::string new_name(it, it_next);
        if( std::find(names_in.begin(), names_in.end(), new_name) != names_in.end() )
          names_out.push_back(new_name);
        if(it_next == tmp_keep_columns.end())
          break;
        it = it_next;
      };
    } else {
      names_out = names_in;
    };
    
  } else if(vm.count("input-raw-udp") && vm.count("keep-columns")) {
    std::stringstream ss;
    ss << vm["input-ip"].as<std::string>() << ":" << vm["input-port"].as<unsigned int>();
    shared_ptr< raw_udp_extractor > data_in_tmp(new raw_udp_extractor());
    data_in = data_in_tmp;
    
    std::string tmp_keep_columns = strip_quotes(vm["keep-columns"].as<std::string>());
    for(std::string::iterator it = tmp_keep_columns.begin(); it != tmp_keep_columns.end(); ++it) {
      std::string::iterator it_next = std::find(it, tmp_keep_columns.end(), ';');
      std::string new_name(it, it_next);
      data_in_tmp->addName(new_name);
      names_in.push_back(new_name);
      names_out.push_back(new_name);
      if(it_next == tmp_keep_columns.end())
        break;
      it = it_next;
    };
    
  } else {
    std::cerr << "Error! Infeasible input specifications! Either you requested raw-udp stream by mistake, or you forgot to specify the 'keep-columns' names for the raw-udp input stream." << std::endl;
    return 2;
  };
  
  if(vm.count("time-column-sync")) {
    if( std::find(names_in.begin(), names_in.end(), strip_quotes(vm["time-column-sync"].as<std::string>())) != names_in.end() ) {
      time_name = strip_quotes(vm["time-column-sync"].as<std::string>());
    };
  };
  
  for(std::size_t i = 0; i < names_out.size(); ++i)
    (*data_out) << names_out[i];
  (*data_out) << data_recorder::end_name_row;
  
  named_value_row nvr_in  = data_in->getFreshNamedValueRow();
  named_value_row nvr_out = data_out->getFreshNamedValueRow();
  
  try {
    stc::time_point t_0 = stc::now();
    while(true) {
      (*data_in) >> nvr_in;
      for(std::size_t i = 0; i < names_out.size(); ++i)
        nvr_out[ names_out[i] ] = nvr_in[ names_out[i] ];
      if( time_name != "" ) {
        // wait until the proper time to output the value.
        stc::time_point t_to_reach = t_0 + ch::duration_cast<stc::duration>(ch::duration<double, ReaKaux::ratio<1,1000> >(nvr_in[time_name]));
        ReaKaux::this_thread::sleep_until( t_to_reach );
      };
      (*data_out) << nvr_out;
    };
  } catch(end_of_record& e) { RK_UNUSED(e); };
  
  
  return 0;
};





















