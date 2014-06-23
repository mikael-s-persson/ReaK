
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

#include <ReaK/core/recorders/data_record_options.hpp>

#include <ReaK/core/recorders/bin_recorder.hpp>
#include <ReaK/core/recorders/ssv_recorder.hpp>
#include <ReaK/core/recorders/tsv_recorder.hpp>
#include <ReaK/core/recorders/network_recorder.hpp>
#include <ReaK/core/recorders/vector_recorder.hpp>


namespace ReaK {

namespace recorder {


shared_ptr< data_recorder > data_stream_options::create_recorder() const {
  
  shared_ptr< data_recorder > result;
  switch(kind) {
    case binary:
      result = shared_ptr< data_recorder >(new bin_recorder());
      break;
    case space_separated:
      result = shared_ptr< data_recorder >(new ssv_recorder());
      break;
    case tab_separated:
      result = shared_ptr< data_recorder >(new tsv_recorder());
      break;
    case tcp_stream:
    case udp_stream: 
    case raw_udp_stream: 
      result = shared_ptr< data_recorder >(new network_recorder());
      break;
    case vector_stream:
      result = shared_ptr< data_recorder >(new vector_recorder());
      break;
  };
  
  if(file_name != "stdout") {
    switch(kind) {
      case tcp_stream:
        result->setFileName("tcp:" + file_name);
        break;
      case udp_stream: 
        result->setFileName("udp:" + file_name);
        break;
      case raw_udp_stream: 
        result->setFileName("raw_udp:" + file_name);
        break;
      default:
        result->setFileName(file_name);
        break;
    };
  } else {
    result->setStream(std::cout);
  };
  
  if( names.size() ) {
    for(std::vector< std::string >::const_iterator it = names.begin(), it_end = names.end(); it != it_end; ++it) {
      (*result) << *it;
    };
    (*result) << data_recorder::end_name_row;
  };
  
  return result;
};


std::pair< shared_ptr< data_extractor >, std::vector< std::string > > data_stream_options::create_extractor() const {
  
  std::pair< shared_ptr< data_extractor >, std::vector< std::string > > result;
  switch(kind) {
    case binary:
      result.first = shared_ptr< data_extractor >(new bin_extractor());
      break;
    case space_separated:
      result.first = shared_ptr< data_extractor >(new ssv_extractor());
      break;
    case tab_separated:
      result.first = shared_ptr< data_extractor >(new tsv_extractor());
      break;
    case tcp_stream:
    case udp_stream:
    case raw_udp_stream:
      result.first = shared_ptr< data_extractor >(new network_extractor());
      break;
    case vector_stream:
      result.first = shared_ptr< data_extractor >(new vector_extractor());
      break;
  };
  
  if(file_name != "stdin") {
    switch(kind) {
      case tcp_stream:
        result.first->setFileName("tcp:" + file_name);
        break;
      case udp_stream: 
        result.first->setFileName("udp:" + file_name);
        break;
      case raw_udp_stream: 
        result.first->setFileName("raw_udp:" + file_name);
        break;
      default:
        result.first->setFileName(file_name);
        break;
    };
  } else {
    result.first->setStream(std::cin);
  };
  
  if( ( kind != raw_udp_stream ) && ( kind != vector_stream ) ) {
    result.second.resize(result.first->getColCount(), "");
    for(std::size_t i = 0; i < result.second.size(); ++i)
      (*result.first) >> result.second[i];
    
    if( names.empty() )
      names = result.second;
    
    if( !time_sync_name.empty() &&
        ( std::find(result.second.begin(), result.second.end(), time_sync_name) == result.second.end() ) )
      throw std::invalid_argument(time_sync_name + " as time-sync column-name");
  } else if( kind == raw_udp_stream ) {
    if( names.empty() )
      throw std::invalid_argument("empty names for a raw-udp-extractor");
    result.second = names;
    
    network_extractor* data_in_tmp = static_cast< network_extractor* >(result.first.get());
    for(std::vector< std::string >::const_iterator it = names.begin(), it_end = names.end(); it != it_end; ++it)
      data_in_tmp->addName(*it);
  } else {
    if( names.empty() )
      throw std::invalid_argument("empty names for a vector-extractor");
    result.second = names;
    
    vector_extractor* data_in_tmp = static_cast< vector_extractor* >(result.first.get());
    for(std::vector< std::string >::const_iterator it = names.begin(), it_end = names.end(); it != it_end; ++it)
      data_in_tmp->addName(*it);
  };
  
  return result;
};


};

};








