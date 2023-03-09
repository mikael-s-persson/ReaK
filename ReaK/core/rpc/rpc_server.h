/**
 * \file rpc_server.h
 *
 * This library declares a singleton class to act as the server to publish and subscribe to
 * remote procedure call (rpc) functions.
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

#ifndef REAK_CORE_RPC_RPC_SERVER_H_
#define REAK_CORE_RPC_RPC_SERVER_H_

#include "detail/remote_function.h"

#include <future>
#include <thread>

namespace ReaK::rpc {

// Singleton
class server {
 private:
  server();
  ~server();

  struct impl;
  impl* pimpl;

 public:
  static server& instance();

  static void set_name(const std::string& aName);

  friend class detail::remote_function;

 private:
  void publish_function(detail::remote_function* aFunc);
  void unpublish_function(detail::remote_function* aFunc);

  detail::call_preparations prepare_call(detail::remote_function* aPFunc) const;
  std::future<detail::call_results> make_remote_call(
      detail::remote_function* aPFunc,
      detail::call_preparations&& pre_data) const;
};

}  // namespace ReaK::rpc

#endif  // REAK_CORE_RPC_RPC_SERVER_H_
