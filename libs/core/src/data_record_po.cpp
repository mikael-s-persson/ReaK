/**
 * \file data_record_po.cpp
 *
 * This library declares utility functions for creating and dealing with program-options related
 * to a data recording or extraction stream (see data_record.hpp). Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2023
 */

/*
 *    Copyright 2023 Sven Mikael Persson
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

#include "ReaK/core/recorders/data_record_po.hpp"
#include "ReaK/core/recorders/data_record_options.hpp"

#include "absl/flags/flag.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"

// Input options.
ABSL_FLAG(std::string, recorder_input, "stdin",
          "Specify the filename for the input data-stream.");
ABSL_FLAG(std::string, recorder_input_ip, "",
          "Specify the IP-address of the data-server.");
ABSL_FLAG(
    int, recorder_input_port, 17017,
    "Specify the IP-port to connect to the data-server (default is 17017).");
ABSL_FLAG(bool, recorder_input_tcp, false,
          "If set, will try to listen to an input TCP data-stream.");
ABSL_FLAG(bool, recorder_input_udp, false,
          "If set, will try to listen to an input UDP data-stream.");
ABSL_FLAG(bool, recorder_input_raw_udp, false,
          "If set, will try to listen to an input RAW UDP data-stream, for "
          "this to work, you must specify the list of columns via the "
          "'recorder_keep_columns' option.");
ABSL_FLAG(std::string, recorder_input_format, "",
          "Specify the format for the input file (default is to use the "
          "file-extension of input-file (ssv, tsv, csv, bin, etc.), or if "
          "piped, use 'ssv').");
ABSL_FLAG(
    int, recorder_input_flush_freq, 50,
    "Specify the flushing frequency of the input data-stream (default is 50Hz, "
    "if set to 0, input will be immediate (as soon as available)).");
ABSL_FLAG(
    int, recorder_input_buffer_size, 512,
    "Specify the desired buffer size of the input data-stream, without "
    "guarantee that the buffer will be limited to that amount (default is 512 "
    "bytes, if set to 0, input will not be buffered, i.e., read as you go).");
ABSL_FLAG(bool, recorder_input_immediate_mode, false,
          "If set, will make the input data-stream without buffering or "
          "regular flushing, i.e., the operations are immediate (note: this "
          "option overrides the flush-freq and buffer-size options).");

// Output options.
ABSL_FLAG(
    std::string, recorder_output, "stdout",
    "Specify the filename for the output file for the output data-stream "
    "(default is to output to the console 'stdout', for piping the output).");
ABSL_FLAG(std::string, recorder_output_ip, "",
          "Specify the IP-address of the RAW UDP output data-stream, this is "
          "because a raw udp stream is connection-less and requires a "
          "pre-defined output IP-address.");
ABSL_FLAG(int, recorder_output_port, 17017,
          "Specify the IP-port for the data-server that will be created "
          "(default is 17017).");
ABSL_FLAG(bool, recorder_output_tcp, false,
          "If set, will output a TCP data-stream.");
ABSL_FLAG(bool, recorder_output_udp, false,
          "If set, will output a UDP data-stream.");
ABSL_FLAG(bool, recorder_output_raw_udp, false,
          "If set, will output a RAW UDP data-stream.");
ABSL_FLAG(std::string, recorder_output_format, "",
          "Specify the format for the output file (default is to use the "
          "file-extension of input-file (ssv, tsv, csv, bin, etc.), or if "
          "piped, use 'ssv').");
ABSL_FLAG(int, recorder_output_flush_freq, 50,
          "Specify the flushing frequency of the output data-stream (default "
          "is 50Hz, if set to 0, output will be immediate (blocking)).");
ABSL_FLAG(
    int, recorder_output_buffer_size, 512,
    "Specify the desired buffer size of the output data-stream, without "
    "guarantee that the buffer will be limited to that amount (default is 512 "
    "bytes, if set to 0, output will not be buffered, i.e., sent as you go).");
ABSL_FLAG(bool, recorder_output_immediate_mode, false,
          "If set, will make the output data-stream without buffering or "
          "regular flushing, i.e., the operations are immediate (note: this "
          "option overrides the flush-freq and buffer-size options).");

// Filtering options.
ABSL_FLAG(std::string, recorder_time_column_sync, "",
          "Name of the time-column from the data-stream to use as a "
          "timed-output, i.e., for a network streaming, this will deliver the "
          "rows at each correct time.");
ABSL_FLAG(std::string, recorder_keep_columns, "",
          "Specify a semi-colon-separated list of the columns to output (keep) "
          "in the data-stream.");

namespace ReaK::recorder {

data_stream_options get_data_stream_options_from_flags(bool aForOutput) {

  data_stream_options result;

  if (!aForOutput) {
    // load options for input.

    std::string input_extension = "ssv";
    if (absl::GetFlag(FLAGS_recorder_input) != "stdin") {
      result.file_name = absl::GetFlag(FLAGS_recorder_input);
      std::vector<std::string_view> substrs =
          absl::StrSplit(result.file_name, '.');
      if (!substrs.empty()) {
        input_extension = std::string(substrs.back());
      }
    } else {
      result.file_name = "stdin";
    }
    if (!absl::GetFlag(FLAGS_recorder_input_format).empty()) {
      input_extension = absl::GetFlag(FLAGS_recorder_input_format);
    }

    if (absl::GetFlag(FLAGS_recorder_input_tcp) ||
        absl::GetFlag(FLAGS_recorder_input_udp) ||
        absl::GetFlag(FLAGS_recorder_input_raw_udp)) {
      result.file_name =
          absl::StrCat(absl::GetFlag(FLAGS_recorder_input_ip), ":",
                       absl::GetFlag(FLAGS_recorder_input_port));
      if (absl::GetFlag(FLAGS_recorder_input_udp)) {
        result.kind = data_stream_options::udp_stream;
      } else if (absl::GetFlag(FLAGS_recorder_input_raw_udp)) {
        result.kind = data_stream_options::raw_udp_stream;
      } else {
        result.kind = data_stream_options::tcp_stream;
      }
    } else {
      if (input_extension == "ssv") {
        result.kind = data_stream_options::space_separated;
      } else if (input_extension == "tsv") {
        result.kind = data_stream_options::tab_separated;
      } else if (input_extension == "csv") {
        result.kind = data_stream_options::comma_separated;
      } else if (input_extension == "bin") {
        result.kind = data_stream_options::binary;
      }
    }

    if (absl::GetFlag(FLAGS_recorder_input_immediate_mode)) {
      result.set_unbuffered();
    } else {
      result.flush_rate = absl::GetFlag(FLAGS_recorder_input_flush_freq);
      result.buffer_size = absl::GetFlag(FLAGS_recorder_input_buffer_size);
    }

  } else {
    // load options for output.

    // take care of 'file_name' and 'kind':
    if (absl::GetFlag(FLAGS_recorder_output_tcp) ||
        absl::GetFlag(FLAGS_recorder_output_udp) ||
        absl::GetFlag(FLAGS_recorder_output_raw_udp)) {
      if (absl::GetFlag(FLAGS_recorder_output_raw_udp)) {
        result.file_name =
            absl::StrCat(absl::GetFlag(FLAGS_recorder_output_ip), ":",
                         absl::GetFlag(FLAGS_recorder_output_port));
        result.kind = data_stream_options::raw_udp_stream;
      } else {
        result.file_name =
            absl::StrCat(absl::GetFlag(FLAGS_recorder_output_port));
        if (absl::GetFlag(FLAGS_recorder_output_udp)) {
          result.kind = data_stream_options::udp_stream;
        } else {
          result.kind = data_stream_options::tcp_stream;
        }
      }
    } else {
      std::string output_extension = "ssv";
      if (absl::GetFlag(FLAGS_recorder_output) != "stdout") {
        result.file_name = absl::GetFlag(FLAGS_recorder_output);

        std::filesystem::create_directories(
            std::filesystem::path(result.file_name).parent_path());

        std::vector<std::string_view> substrs =
            absl::StrSplit(result.file_name, '.');
        if (!substrs.empty()) {
          output_extension = std::string(substrs.back());
        }
      } else {
        result.file_name = "stdout";
      }

      if (!absl::GetFlag(FLAGS_recorder_output_format).empty()) {
        output_extension = absl::GetFlag(FLAGS_recorder_output_format);
      }

      if (output_extension == "ssv") {
        result.kind = data_stream_options::space_separated;
      } else if (output_extension == "tsv") {
        result.kind = data_stream_options::tab_separated;
      } else if (output_extension == "csv") {
        result.kind = data_stream_options::comma_separated;
      } else if (output_extension == "bin") {
        result.kind = data_stream_options::binary;
      }
    }

    if (absl::GetFlag(FLAGS_recorder_output_immediate_mode)) {
      result.set_unbuffered();
    } else {
      result.flush_rate = absl::GetFlag(FLAGS_recorder_output_flush_freq);
      result.buffer_size = absl::GetFlag(FLAGS_recorder_output_buffer_size);
    }
  }

  // Load filtering options.
  result.names = absl::StrSplit(absl::GetFlag(FLAGS_recorder_keep_columns), ';',
                                absl::SkipWhitespace());
  result.time_sync_name = absl::GetFlag(FLAGS_recorder_time_column_sync);

  return result;
}

}  // namespace ReaK::recorder
