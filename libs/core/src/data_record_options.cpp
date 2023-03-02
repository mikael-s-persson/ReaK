
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

#include <ReaK/core/recorders/ascii_recorder.hpp>
#include <ReaK/core/recorders/bin_recorder.hpp>
#include <ReaK/core/recorders/network_recorder.hpp>
#include <ReaK/core/recorders/vector_recorder.hpp>

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <algorithm>
#include <boost/graph/graph_concepts.hpp>

#include <ctime>

namespace ReaK::recorder {

std::string data_stream_options::get_URI() const {
  std::string result;

  switch (kind) {
    case binary:
    case space_separated:
    case tab_separated:
    case comma_separated:
      result += "file:";
      break;
    case tcp_stream:
      result += "tcp:";
      break;
    case udp_stream:
      result += "udp:";
      break;
    case raw_udp_stream:
      result += "raw_udp:";
      break;
    case vector_stream:
      result = "memory:";
      break;
  }

  result += file_name;
  std::size_t dot = result.find_last_of('.');
  if ((dot != std::string::npos) && (result.size() - dot == 4) &&
      ((result.substr(dot, 4) == ".bin") || (result.substr(dot, 4) == ".ssv") ||
       (result.substr(dot, 4) == ".tsv") ||
       (result.substr(dot, 4) == ".csv"))) {
    result.erase(dot, 4);
  }
  switch (kind) {
    case binary:
      result += ".bin";
      break;
    case space_separated:
      result += ".ssv";
      break;
    case tab_separated:
      result += ".tsv";
      break;
    case comma_separated:
      result += ".csv";
      break;
    default:
      break;
  }

  if (!time_sync_name.empty()) {
    result += "|time=" + time_sync_name;
  }

  if (flush_rate != 50) {
    std::stringstream ss;
    ss << "|freq=" << flush_rate;
    result += ss.str();
  }

  if (buffer_size != 500) {
    std::stringstream ss;
    ss << "|buf=" << buffer_size;
    result += ss.str();
  }

  return result;
}

void data_stream_options::set_from_URI(const std::string& aURI) {
  std::size_t cur_i = 0;
  std::size_t cur_sep = aURI.find(':');
  std::string hdr_name = aURI.substr(0, cur_sep);
  if (cur_sep != std::string::npos) {
    cur_i = cur_sep + 1;
  } else {
    cur_i = std::string::npos;
  }

  if (hdr_name == "memory") {
    kind = vector_stream;
    cur_sep = aURI.find('|', cur_i);
    cur_i = cur_sep;
  } else if (hdr_name == "file") {
    cur_sep = aURI.find('|', cur_i);
    file_name = aURI.substr(cur_i, cur_sep - cur_i);
    cur_i = cur_sep;
    std::size_t dot = file_name.find_last_of('.');
    if ((dot == std::string::npos) || (file_name.size() - dot != 4)) {
      kind = space_separated;
    } else if (file_name.substr(dot, 4) == ".bin") {
      kind = binary;
    } else if (file_name.substr(dot, 4) == ".tsv") {
      kind = tab_separated;
    } else if (file_name.substr(dot, 4) == ".csv") {
      kind = comma_separated;
    } else {
      kind = space_separated;
    }
    if (file_name.substr(0, dot) == "stdout") {
      file_name = "stdout";
    }
    if (file_name.substr(0, dot) == "stdin") {
      file_name = "stdin";
    }
  } else {
    if (hdr_name == "tcp") {
      kind = tcp_stream;
    } else if (hdr_name == "udp") {
      kind = udp_stream;
    } else if (hdr_name == "raw_udp") {
      kind = raw_udp_stream;
    }
    cur_sep = aURI.find('|', cur_i);
    file_name = aURI.substr(cur_i, cur_sep - cur_i);
    cur_i = cur_sep;
  }

  // Default values:
  time_sync_name = "";
  flush_rate = 50;
  buffer_size = 500;
  while (cur_i != std::string::npos) {
    cur_sep = aURI.find('|', cur_i);
    std::string cur_sec = aURI.substr(cur_i, cur_sep - cur_i);
    cur_i = cur_sep;
    std::size_t eq_sign = cur_sec.find('=');
    std::string cur_sec_hdr = cur_sec.substr(0, eq_sign);
    if (cur_sec_hdr == "|time") {
      time_sync_name = cur_sec.substr(eq_sign + 1);
    } else if (cur_sec_hdr == "|freq") {
      std::string cur_sec_val = cur_sec.substr(eq_sign + 1);
      std::stringstream ss(cur_sec_val);
      ss >> flush_rate;
    } else if (cur_sec_hdr == "|buf") {
      std::string cur_sec_val = cur_sec.substr(eq_sign + 1);
      std::stringstream ss(cur_sec_val);
      ss >> buffer_size;
    }
  }
}

std::shared_ptr<data_recorder> data_stream_options::create_recorder() const {

  std::shared_ptr<data_recorder> result;
  switch (kind) {
    case binary:
      result = std::make_shared<bin_recorder>();
      break;
    case space_separated:
    case tab_separated:
    case comma_separated:
      result = std::make_shared<ascii_recorder>();
      break;
    case tcp_stream:
    case udp_stream:
    case raw_udp_stream:
      result = std::make_shared<network_recorder>();
      break;
    case vector_stream:
      result = std::make_shared<vector_recorder>();
      break;
  }

  result->setFlushSampleRate(flush_rate);
  result->setMaxBufferSize(buffer_size);

  switch (kind) {
    case space_separated:
      static_cast<ascii_recorder*>(result.get())->delimiter = " ";
      break;
    case tab_separated:
      static_cast<ascii_recorder*>(result.get())->delimiter = "\t";
      break;
    case comma_separated:
      static_cast<ascii_recorder*>(result.get())->delimiter = ",";
      break;
    default:
      break;
  }

  if (file_name != "stdout") {
    std::string new_fname;

    // Check if the file name contains a wildcard for a date or time stamp, or both:
    std::size_t d = file_name.find("$d");
    std::size_t t = file_name.find("$t");
    if ((d != std::string::npos) || (t != std::string::npos)) {
      // Record the current date and time:
      std::time_t t_ctime = std::time(nullptr);
      std::array<char, 16> cdate_as_str;
      if (std::strftime(cdate_as_str.data(), cdate_as_str.size(), "%Y%m%d",
                        std::localtime(&t_ctime)) == 0) {
        cdate_as_str[0] = '\0';
      }
      std::array<char, 16> ctime_as_str;
      if (std::strftime(ctime_as_str.data(), ctime_as_str.size(), "%H%M%S",
                        std::localtime(&t_ctime)) == 0) {
        ctime_as_str[0] = '\0';
      }

      if (d < t) {
        new_fname.append(file_name.begin(), file_name.begin() + d);
        new_fname.append(cdate_as_str.data());
        if (t != std::string::npos) {
          new_fname.append(file_name.begin() + d + 2, file_name.begin() + t);
          new_fname.append(ctime_as_str.data());
          new_fname.append(file_name.begin() + t + 2, file_name.end());
        } else {
          new_fname.append(file_name.begin() + d + 2, file_name.end());
        }
      } else if (t < d) {
        new_fname.append(file_name.begin(), file_name.begin() + t);
        new_fname.append(ctime_as_str.data());
        if (d != std::string::npos) {
          new_fname.append(file_name.begin() + t + 2, file_name.begin() + d);
          new_fname.append(cdate_as_str.data());
          new_fname.append(file_name.begin() + d + 2, file_name.end());
        } else {
          new_fname.append(file_name.begin() + t + 2, file_name.end());
        }
      }
    } else {
      new_fname = file_name;
    }

    switch (kind) {
      case tcp_stream:
        result->setFileName("tcp:" + new_fname);
        break;
      case udp_stream:
        result->setFileName("udp:" + new_fname);
        break;
      case raw_udp_stream:
        result->setFileName("raw_udp:" + new_fname);
        break;
      default:
        result->setFileName(new_fname);
        break;
    }
  } else {
    result->setStream(std::cout);
  }

  if (static_cast<unsigned int>(!names.empty()) != 0U) {
    for (auto& name : names) {
      (*result) << name;
    }
    (*result) << data_recorder::end_name_row;
  }

  return result;
}

std::pair<std::shared_ptr<data_extractor>, std::vector<std::string>>
data_stream_options::create_extractor() const {

  std::pair<std::shared_ptr<data_extractor>, std::vector<std::string>> result;
  switch (kind) {
    case binary:
      result.first = std::make_shared<bin_extractor>();
      break;
    case space_separated:
    case tab_separated:
    case comma_separated:
      result.first = std::make_shared<ascii_extractor>();
      break;
    case tcp_stream:
    case udp_stream:
    case raw_udp_stream:
      result.first = std::make_shared<network_extractor>();
      break;
    case vector_stream:
      result.first = std::make_shared<vector_extractor>();
      break;
  }

  result.first->setFlushSampleRate(flush_rate);
  result.first->setMinBufferSize(buffer_size);

  switch (kind) {
    case space_separated:
      static_cast<ascii_extractor*>(result.first.get())->delimiter = " ";
      break;
    case tab_separated:
      static_cast<ascii_extractor*>(result.first.get())->delimiter = "\t";
      break;
    case comma_separated:
      static_cast<ascii_extractor*>(result.first.get())->delimiter = ",";
      break;
    default:
      break;
  }

  if (kind == raw_udp_stream) {
    if (names.empty()) {
      throw std::invalid_argument("empty names for a raw-udp-extractor");
    }
    result.second = names;

    auto* data_in_tmp = static_cast<network_extractor*>(result.first.get());
    for (auto& name : names) {
      data_in_tmp->addName(name);
    }
  } else if (kind == vector_stream) {
    if (names.empty()) {
      throw std::invalid_argument("empty names for a vector-extractor");
    }
    result.second = names;

    auto* data_in_tmp = static_cast<vector_extractor*>(result.first.get());
    for (auto& name : names) {
      data_in_tmp->addName(name);
    }
  }

  if (file_name != "stdin") {
    switch (kind) {
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
    }
  } else {
    result.first->setStream(std::cin);
  }

  if ((kind != raw_udp_stream) && (kind != vector_stream)) {
    result.second.resize(result.first->getColCount(), "");
    for (std::size_t i = 0; i < result.second.size(); ++i) {
      (*result.first) >> result.second[i];
    }

    if (names.empty()) {
      names = result.second;
    }

    if (!time_sync_name.empty() &&
        (std::find(result.second.begin(), result.second.end(),
                   time_sync_name) == result.second.end())) {
      throw std::invalid_argument(time_sync_name + " as time-sync column-name");
    }
  }

  return result;
}

void data_stream_options::load_all_configs(const std::string& aFileName) {

  std::shared_ptr<serialization::iarchive> p_ia =
      serialization::open_iarchive(aFileName);

  std::string stream_kind;

  (*p_ia) & RK_SERIAL_LOAD_WITH_NAME(stream_kind) &
      RK_SERIAL_LOAD_WITH_NAME(file_name) & RK_SERIAL_LOAD_WITH_NAME(names) &
      RK_SERIAL_LOAD_WITH_NAME(time_sync_name) &
      RK_SERIAL_LOAD_WITH_NAME(flush_rate) &
      RK_SERIAL_LOAD_WITH_NAME(buffer_size);

  if (stream_kind == "binary") {
    kind = binary;
  } else if (stream_kind == "space_separated") {
    kind = space_separated;
  } else if (stream_kind == "tab_separated") {
    kind = tab_separated;
  } else if (stream_kind == "comma_separated") {
    kind = comma_separated;
  } else if (stream_kind == "tcp_stream") {
    kind = tcp_stream;
  } else if (stream_kind == "udp_stream") {
    kind = udp_stream;
  } else if (stream_kind == "raw_udp_stream") {
    kind = raw_udp_stream;
  } else {  // if(stream_kind == "vector_stream")
    kind = vector_stream;
  }
}

void data_stream_options::save_all_configs(const std::string& aFileName) const {

  std::shared_ptr<serialization::oarchive> p_oa =
      serialization::open_oarchive(aFileName);

  std::string stream_kind;

  switch (kind) {
    case binary:
      stream_kind = "binary";
      break;
    case space_separated:
      stream_kind = "space_separated";
      break;
    case tab_separated:
      stream_kind = "tab_separated";
      break;
    case comma_separated:
      stream_kind = "comma_separated";
      break;
    case tcp_stream:
      stream_kind = "tcp_stream";
      break;
    case udp_stream:
      stream_kind = "udp_stream";
      break;
    case raw_udp_stream:
      stream_kind = "raw_udp_stream";
      break;
    default:
      stream_kind = "vector_stream";
      break;
  }

  (*p_oa) & RK_SERIAL_SAVE_WITH_NAME(stream_kind) &
      RK_SERIAL_SAVE_WITH_NAME(file_name) & RK_SERIAL_SAVE_WITH_NAME(names) &
      RK_SERIAL_SAVE_WITH_NAME(time_sync_name) &
      RK_SERIAL_SAVE_WITH_NAME(flush_rate) &
      RK_SERIAL_SAVE_WITH_NAME(buffer_size);
}

}  // namespace ReaK::recorder
