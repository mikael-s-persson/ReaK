/**
 * \file version.hpp
 *
 * This library declares helpers to the ReaK.RPC protocol.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2014
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

#ifndef REAK_RPC_VERSION_HPP
#define REAK_RPC_VERSION_HPP

#include <array>
#include <string>

namespace ReaK::rpc {

extern const double rpc_protocol_version;

enum msg_format {
  binary_format = 0,
  xml_format,
  protobuf_format,
  msg_format_count
};

extern const std::array<const char*, msg_format_count> msg_format_to_str;

}  // namespace ReaK::rpc

#endif  // REAK_RPC_VERSION_HPP
