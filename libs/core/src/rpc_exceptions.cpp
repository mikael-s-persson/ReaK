
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

#include "ReaK/core/rpc/rpc_exceptions.hpp"

#include "ReaK/core/rpc/detail/remote_function.hpp"

#include <iomanip>
#include <sstream>

namespace ReaK::rpc {

publishing_mismatch::publishing_mismatch(const std::string& aAction,
                                         detail::remote_function* aGiven,
                                         detail::remote_function* aExisting) {
  std::stringstream ss;
  ss << "Could not " << aAction
     << " function because of a mismatch with existing RPC function! Tried to "
     << aAction << " function '" << aGiven->name << ":"
     << aGiven->get_params_hash() << "' at address 0x" << std::hex << aGiven
     << ", but already had function '" << aExisting->name << ":"
     << aExisting->get_params_hash() << "' at address 0x" << std::hex
     << aExisting << "!" << std::flush;
  msg = ss.str();
};

unrecognized_function::unrecognized_function(detail::remote_function* aPFunc) {
  std::stringstream ss;
  ss << "Unrecognized function '" << aPFunc->name << ":"
     << aPFunc->get_params_hash()
     << "' from host '" + aPFunc->get_host() + "'.";
  set_message(ss.str());
};

marshalling_error::marshalling_error(detail::remote_function* aPFunc,
                                     const std::string& aSerialExceptMsg) {
  std::stringstream ss;
  ss << "Mashalling error with function '" << aPFunc->name << ":"
     << aPFunc->get_params_hash() << "' with marshalling error message: '"
     << aSerialExceptMsg << "'.";
  set_message(ss.str());
};

}  // namespace ReaK::rpc
