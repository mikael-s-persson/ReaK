
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

#include <cmath>

#include <ReaK/core/rpc/detail/rpc_headers.hpp>
#include <ReaK/core/rpc/rpc_exceptions.hpp>

#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ReaK::rpc::detail {

static const std::array<const char*, rpc_header_type_count>
    rpc_header_type_to_str = {"CALL", "RETURN", "EXCEPTION", "UNRECOGNIZED",
                              "STARTPORTSERVER"};

std::string generate_rpc_header(rpc_header_type h, msg_format fmt,
                                std::size_t msg_size) {
  std::stringstream ss;
  ss << "REAK.RPC " << std::setprecision(2) << rpc_protocol_version << " "
     << rpc_header_type_to_str[h] << " " << msg_format_to_str[fmt] << " "
     << msg_size << "\r\n"
     << std::flush;
  return ss.str();
}

std::tuple<rpc_header_type, msg_format, std::size_t> parse_rpc_header(
    std::istream& hdr_in) {

  std::string hdr_lead;
  if (!(hdr_in >> hdr_lead) || (hdr_lead != "REAK.RPC")) {
    throw communication_error(
        "The header received does not start with 'REAK.RPC!'");
  }

  double hdr_ver = NAN;
  if (!(hdr_in >> hdr_ver) || (hdr_ver > rpc_protocol_version + 0.001)) {
    throw communication_error(
        "The header received demands a version of ReaK.RPC which is not "
        "supported by this implementation!");
  }

  std::tuple<rpc_header_type, msg_format, std::size_t> result;

  std::string hdr_type;
  if (!(hdr_in >> hdr_type)) {
    throw communication_error(
        "The header received does not have a valid type identifier!");
  }
  if (hdr_type == rpc_header_type_to_str[rpc_call]) {
    std::get<0>(result) = rpc_call;
  } else if (hdr_type == rpc_header_type_to_str[rpc_return]) {
    std::get<0>(result) = rpc_return;
  } else if (hdr_type == rpc_header_type_to_str[rpc_exception]) {
    std::get<0>(result) = rpc_exception;
  } else if (hdr_type == rpc_header_type_to_str[rpc_unrecognized]) {
    std::get<0>(result) = rpc_unrecognized;
  } else if (hdr_type == rpc_header_type_to_str[rpc_startportserver]) {
    std::get<0>(result) = rpc_startportserver;
  }

  std::string hdr_format;
  if (!(hdr_in >> hdr_format)) {
    throw communication_error(
        "The header received does not have a valid format identifier!");
  }
  if (hdr_format == msg_format_to_str[binary_format]) {
    std::get<1>(result) = binary_format;
  } else if (hdr_format == msg_format_to_str[xml_format]) {
    std::get<1>(result) = xml_format;
  } else if (hdr_format == msg_format_to_str[protobuf_format]) {
    std::get<1>(result) = protobuf_format;
  }

  if (!(hdr_in >> std::get<2>(result))) {
    throw communication_error(
        "The header received does not have a valid format identifier!");
  }

  while (true) {
    int p = hdr_in.peek();
    if ((p == '\r') || (p == '\n')) {
      hdr_in.get();
      continue;
    }
    break;
  }

  return result;
}

}  // namespace ReaK::rpc::detail
