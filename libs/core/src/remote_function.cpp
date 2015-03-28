
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


#include <ReaK/core/rpc/detail/remote_function.hpp>

#include <ReaK/core/rpc/rpc_server.hpp>
#include <ReaK/core/base/thread_incl.hpp>

namespace ReaK {

namespace rpc {

namespace detail {


remote_function::remote_function( const std::string& aName ) : name( aName ){};

remote_function::~remote_function(){};

void remote_function::publish() { server::instance().publish_function( this ); };

void remote_function::unpublish() { server::instance().unpublish_function( this ); };


call_preparations remote_function::prepare_for_call() { return server::instance().prepare_call( this ); };

call_results remote_function::do_remote_call( call_preparations&& pre_data ) {
  ReaKaux::future< call_results > fp_ari = server::instance().make_remote_call( this, std::move( pre_data ) );
  return fp_ari.get();
};
};
};
};
