
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

#include <iomanip>
#include <sstream>

#include <chrono>
#include <thread>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

namespace ch = std::chrono;
typedef ch::steady_clock stc;
typedef ch::high_resolution_clock hrc;

ABSL_FLAG(
    bool, echo, false,
    "Echo all the output to the terminal (do not use with stdout sreaming).");
ABSL_FLAG(std::string, add_relative_time, "",
          "Add a relative timing of the received data to the output stream, as "
          "first column, with the name given.");

int main(int argc, char** argv) {
  using namespace ReaK;
  using namespace recorder;

  absl::ParseCommandLine(argc, argv);

  try {
    data_stream_options data_in_opt = get_data_stream_options_from_flags(false);
    data_stream_options data_out_opt = get_data_stream_options_from_flags(true);

    std::shared_ptr<data_extractor> data_in;
    std::vector<std::string> names_in;
    std::tie(data_in, names_in) = data_in_opt.create_extractor();

    if (data_out_opt.names.size() == 0) {
      data_out_opt.names = names_in;
    } else {
      names_in = data_out_opt.names;
    }
    if (!absl::GetFlag(FLAGS_add_relative_time).empty()) {
      data_out_opt.names.insert(data_out_opt.names.begin(),
                                absl::GetFlag(FLAGS_add_relative_time));
    }
    std::shared_ptr<data_recorder> data_out = data_out_opt.create_recorder();

    named_value_row nvr_in = data_in->getFreshNamedValueRow();
    named_value_row nvr_out = data_out->getFreshNamedValueRow();

    stc::time_point t_0 = stc::now();
    hrc::time_point hrt_0 = hrc::now();
    double in_time_0 = std::numeric_limits<double>::infinity();
    while (true) {

      (*data_in) >> nvr_in;
      hrc::time_point hrt_1 = hrc::now();

      for (std::size_t i = 0; i < names_in.size(); ++i) {
        try {
          nvr_out[names_in[i]] = nvr_in[names_in[i]];
        } catch (out_of_bounds& e) {
          RK_UNUSED(e);
          nvr_out[names_in[i]] = 0.0;
        }
      }
      if (!data_out_opt.time_sync_name.empty()) {
        if (in_time_0 == std::numeric_limits<double>::infinity()) {
          in_time_0 = nvr_in[data_out_opt.time_sync_name];
        }
        // wait until the proper time to output the value.
        stc::time_point t_to_reach =
            t_0 + ch::duration_cast<stc::duration>(
                      ch::duration<double, std::ratio<1, 1>>(
                          nvr_in[data_out_opt.time_sync_name] - in_time_0));
        std::this_thread::sleep_until(t_to_reach);
      } else if (!absl::GetFlag(FLAGS_add_relative_time).empty()) {
        nvr_out[absl::GetFlag(FLAGS_add_relative_time)] =
            ch::duration_cast<ch::duration<double>>(hrt_1 - hrt_0).count();
      }
      if (absl::GetFlag(FLAGS_echo)) {
        std::cout << "\n\n\n\n";
        if (!absl::GetFlag(FLAGS_add_relative_time).empty()) {
          std::cout << absl::GetFlag(FLAGS_add_relative_time) << '\t'
                    << std::setw(16)
                    << nvr_out[absl::GetFlag(FLAGS_add_relative_time)] << '\n';
        }
        for (std::size_t i = 0; i < names_in.size(); ++i) {
          std::cout << names_in[i] << '\t' << std::setw(16)
                    << nvr_out[names_in[i]] << '\n';
        }
        std::cout << std::flush;
      }
      (*data_out) << nvr_out;
    }

  } catch (std::invalid_argument& e) {
    std::cerr << "Error! Creation of data-streams failed! Invalid argument: "
              << e.what() << std::endl;
    return 1;
  } catch (end_of_record& e) {
    RK_UNUSED(e);
  }

  return 0;
}
