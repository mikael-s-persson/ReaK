/**
 * \file data_record_options.h
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

#ifndef REAK_CORE_RECORDERS_DATA_RECORD_OPTIONS_H_
#define REAK_CORE_RECORDERS_DATA_RECORD_OPTIONS_H_

#include "ReaK/core/recorders/data_record.h"

#include <string>

/** Main namespace for ReaK */
namespace ReaK::recorder {

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
    comma_separated,
    tcp_stream,
    udp_stream,
    raw_udp_stream,
    vector_stream
  } kind{space_separated};  ///< Stores the kind of stream (format) to use.

  /**
   * This function returns the file-extension for the data-stream file (if it's a file-based stream).
   * \return The file-extension for the data-stream file (if it's a file-based stream).
   */
  std::string get_extension() const {
    switch (kind) {
      case binary:
        return "bin";
      case space_separated:
        return "ssv";
      case tab_separated:
        return "tsv";
      case comma_separated:
        return "csv";
      default:
        return "";
    }
  }

  /// Stores the file-name (or the ip:port name) for the stream.
  std::string file_name;

  /// Stores the names to put on a recorder or to keep from an extractor, will be ignored if empty.
  mutable std::vector<std::string> names;

  /**
   * This function adds a name to the 'names' vector.
   * \param aName The name to be added to the 'names' vector.
   * \return A reference to this object.
   */
  data_stream_options& add_name(const std::string& aName) {
    names.push_back(aName);
    return *this;
  }

  /// Stores the name of time-sync column for when streaming is paced with time.
  std::string time_sync_name;

  /// Stores the frequency (Hz) of data flushes, if 0, then data is sent / received immediately.
  unsigned int flush_rate{50};

  /// Stores the desired average size of the data buffer, if 0, then data is sent / received immediately without
  /// buffering.
  unsigned int buffer_size{500};

  void set_unbuffered() {
    flush_rate = 0;
    buffer_size = 0;
  }
  void set_buffered(unsigned int aFlushRate = 50,
                    unsigned int aBufferSize = 500) {
    flush_rate = aFlushRate;
    buffer_size = aBufferSize;
  }

  std::string get_uri() const;
  void set_from_uri(const std::string& aURI);

  /**
   * Default constructor.
   */
  data_stream_options() = default;

  /**
   * This function creates a data-recorder from the options that were set in
   * this options object.
   * \return A data-recorder corresponding to the options that were set in this object.
   */
  std::shared_ptr<data_recorder> create_recorder() const;

  /**
   * This function creates a data-extractor (with a list of names of columns from it)
   * from the options that were set in this options object.
   * \return A data-extractor (with a list of names of columns from it) corresponding to
   *         the options that were set in this object.
   */
  std::pair<std::shared_ptr<data_extractor>, std::vector<std::string>>
  create_extractor() const;

  /**
   * Loads all the configurations from the given file-name (ReaK archive).
   * \param a_file_name The file-name of the ReaK archive from which to load all the configurations.
   */
  void load_all_configs(const std::string& a_file_name);

  /**
   * Saves all the configurations to the given file-name (ReaK archive).
   * \param a_file_name The file-name of the ReaK archive to which to save all the configurations.
   */
  void save_all_configs(const std::string& a_file_name) const;
};

}  // namespace ReaK::recorder

#endif  // REAK_CORE_RECORDERS_DATA_RECORD_OPTIONS_H_
