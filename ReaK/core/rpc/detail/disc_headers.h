/**
 * \file disc_headers.h
 *
 * This library declares helpers to generate and parse the headers of the ReaK.DISC protocol.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2015
 */

/*
 *    Copyright 2015 Sven Mikael Persson
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

#ifndef REAK_CORE_RPC_DETAIL_DISC_HEADERS_H_
#define REAK_CORE_RPC_DETAIL_DISC_HEADERS_H_

#include "ReaK/core/serialization/archiver.h"
#include "ReaK/core/serialization/serializable.h"

#include "ReaK/core/rpc/version.h"

#include <iostream>
#include <string>
#include <tuple>
#include <utility>

namespace ReaK {

namespace rpc::detail {

enum disc_header_type {
  disc_newport = 0,
  disc_giveport,
  disc_advertise,
  disc_unadvertise,
  disc_echo,
  disc_header_count
};

std::string generate_disc_header(disc_header_type h, msg_format fmt,
                                 std::size_t msg_size);

std::tuple<disc_header_type, msg_format, std::size_t> parse_disc_header(
    std::istream& hdr_in);

struct disc_basic_header : serializable {
  std::string name;
  std::string addr;

  explicit disc_basic_header(std::string aName = "",
                             std::string aAddr = "127.0.0.0")
      : name(std::move(aName)), addr(std::move(aAddr)) {}

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(name) & RK_SERIAL_SAVE_WITH_NAME(addr);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(name) & RK_SERIAL_LOAD_WITH_NAME(addr);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(disc_basic_header, 1, serializable)
};

inline bool operator<(const disc_basic_header& lhs,
                      const disc_basic_header& rhs) {
  if (lhs.name < rhs.name) {
    return true;
  }
  if (lhs.name == rhs.name) {
    if (lhs.addr < rhs.addr) {
      return true;
    }
  }
  return false;
}

struct disc_with_port_header : disc_basic_header {
  unsigned int port;

  explicit disc_with_port_header(const std::string& aName = "",
                                 const std::string& aAddr = "127.0.0.0",
                                 unsigned int aPort = 0)
      : disc_basic_header(aName, aAddr), port(aPort) {}

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    disc_basic_header::save(A, 1);
    A& RK_SERIAL_SAVE_WITH_NAME(port);
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    disc_basic_header::load(A, 1);
    A& RK_SERIAL_LOAD_WITH_NAME(port);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(disc_with_port_header, 1, disc_basic_header)
};

}  // namespace rpc::detail

namespace rtti {

#define RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS(VARIABLE, CLASSID)        \
  template <>                                                                 \
  struct get_type_id<::ReaK::rpc::detail::VARIABLE> {                         \
    static constexpr unsigned int ID = CLASSID;                               \
    static constexpr auto type_name = std::string_view{#VARIABLE};            \
    static construct_ptr CreatePtr() noexcept { return nullptr; }             \
                                                                              \
    typedef const serializable& save_type;                                    \
    typedef serializable& load_type;                                          \
  };                                                                          \
                                                                              \
  template <typename Tail>                                                    \
  struct get_type_info<::ReaK::rpc::detail::VARIABLE, Tail> {                 \
    typedef type_id<::ReaK::rpc::detail::VARIABLE, typename Tail::type> type; \
    static constexpr auto type_name =                                         \
        ct_concat_v<get_type_id<std::string>::type_name,                      \
                    get_type_name_tail<Tail>::value>;                         \
  };

RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS(disc_basic_header, 0xC1600001)
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS(disc_with_port_header, 0xC1600002)

#undef RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_CORE_RPC_DETAIL_DISC_HEADERS_H_
