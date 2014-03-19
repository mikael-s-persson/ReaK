/**
 * \file data_record_options.hpp
 *
 * This library declares utility functions for building option-sets related 
 * to a data recording or extraction stream (see data_record.hpp), and using those options 
 * in a factory-function to create an abstract stream. Here, "data" is meant as
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

#ifndef REAK_DATA_RECORD_OPTIONS_HPP
#define REAK_DATA_RECORD_OPTIONS_HPP

#include "data_record.hpp"

#include <string>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Data Recorders and Extractors */
namespace recorder {


/**
 * This class stores a number of options related to data-streaming with either data-recorders 
 * or with data-extractors. Most of the options are to be set before calling any "create" 
 * functions.
 * \note This class is mainly intended to be used with data_record_po (program-options).
 */
struct data_stream_options {
  
  enum stream_type {
    binary = 0,
    space_separated,
    tab_separated,
    tcp_stream,
    udp_stream,
    raw_udp_stream
  } kind; ///< Stores the kind of stream (format) to use.
  
  /// Stores the file-name (or the ip:port name) for the stream.
  std::string file_name;
  
  /// Stores the names to put on a recorder or to keep from an extractor, will be ignored if empty.
  std::vector< std::string > names;
  
  /// Stores the name of time-sync column for when streaming is paced with time.
  std::string time_sync_name;
  
  /**
   * Default constructor.
   */
  data_stream_options() : kind(space_separated) { };
  
  /**
   * This function creates a data-recorder from the options that were set in
   * this options object.
   * \return A data-recorder corresponding to the options that were set in this object.
   */
  shared_ptr< data_recorder > create_recorder() const;
  
  /**
   * This function creates a data-extractor (with a list of names of columns from it) 
   * from the options that were set in this options object.
   * \return A data-extractor (with a list of names of columns from it) corresponding to 
   *         the options that were set in this object.
   */
  std::pair< shared_ptr< data_extractor >, std::vector< std::string > > create_extractor() const;
  
};



};

};

#endif // RK_DATA_RECORD_PO_HPP












