
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

#include <ReaK/core/recorders/data_record_po.hpp>

#include <sstream>
#include <iomanip>

#include <ReaK/core/base/chrono_incl.hpp>
#include <ReaK/core/base/thread_incl.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace ch = ReaKaux::chrono;
typedef ch::steady_clock stc;

int main(int argc, char** argv) {
  using namespace ReaK;
  using namespace recorder;
  
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
    ("echo", "echo all the output to the terminal (do not use with std-out streaming).")
  ;
  
  po::options_description io_options = get_data_stream_options_po_desc(true, true);
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  try {
    data_stream_options data_in_opt  = get_data_stream_options_from_po(vm, false);
    data_stream_options data_out_opt = get_data_stream_options_from_po(vm, true);
    
    shared_ptr< data_extractor > data_in;
    std::vector<std::string> names_in;
    boost::tie(data_in, names_in) = data_in_opt.create_extractor();
    
    if(data_out_opt.names.size() == 0)
      data_out_opt.names = names_in;
    else
      names_in = data_out_opt.names;
    shared_ptr< data_recorder > data_out = data_out_opt.create_recorder();
    
    named_value_row nvr_in  = data_in->getFreshNamedValueRow();
    named_value_row nvr_out = data_out->getFreshNamedValueRow();
    
    stc::time_point t_0 = stc::now();
    while(true) {
      
      (*data_in) >> nvr_in;
      
      if(vm.count("echo"))
        std::cout << "\n\n\n\n";
      
      for(std::size_t i = 0; i < names_in.size(); ++i) {
        try {
          nvr_out[ names_in[i] ] = nvr_in[ names_in[i] ];
        } catch(out_of_bounds& e) { RK_UNUSED(e);
          nvr_out[ names_in[i] ] = 0.0;
        };
        if(vm.count("echo"))
          std::cout << names_in[i] << '\t' << std::setw(16) << nvr_out[ names_in[i] ] << '\n';
      };
      if(vm.count("echo"))
        std::cout << std::flush;
      if( !data_out_opt.time_sync_name.empty() ) {
        // wait until the proper time to output the value.
        stc::time_point t_to_reach = t_0 + ch::duration_cast<stc::duration>(ch::duration<double, ReaKaux::ratio<1,1> >(nvr_in[data_out_opt.time_sync_name]));
        ReaKaux::this_thread::sleep_until( t_to_reach );
      };
      (*data_out) << nvr_out;
    };
    
  } catch(std::invalid_argument& e) {
    std::cerr << "Error! Creation of data-streams failed! Invalid argument: " << e.what() << std::endl;
    return 1;
  } catch(end_of_record& e) { 
    RK_UNUSED(e); 
  };
  
  return 0;
};





















