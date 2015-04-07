/*
 * \file remote_function.hpp
 *
 * This library declares and defines some underlying details for the implementation of the
 * remote procedure call (rpc) library from ReaK. None of the elements within this header
 * are meant to be used in user-side code, only library code.
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


#ifndef REAK_DETAIL_REMOTE_FUNCTION_HPP
#define REAK_DETAIL_REMOTE_FUNCTION_HPP

#include <ReaK/core/rpc/version.hpp>

#include <ReaK/core/serialization/archiver.hpp>
#include <sstream>
#include <string>

namespace ReaK {

namespace rpc {

namespace detail {


struct call_preparations {
  msg_format fmt;
  std::unique_ptr< std::stringstream > ss;
  std::unique_ptr< serialization::oarchive > p_aro;
  std::size_t call_seq;

  call_preparations() BOOST_NOEXCEPT : fmt(), ss(), p_aro(), call_seq(){};

  call_preparations( call_preparations&& rhs ) BOOST_NOEXCEPT : fmt( rhs.fmt ),
                                                                ss( std::move( rhs.ss ) ),
                                                                p_aro( std::move( rhs.p_aro ) ),
                                                                call_seq( rhs.call_seq ) {
    rhs.call_seq = 0;
  };
  call_preparations& operator=( call_preparations&& rhs ) BOOST_NOEXCEPT {
    fmt = rhs.fmt;
    ss = std::move( rhs.ss );
    p_aro = std::move( rhs.p_aro );
    call_seq = rhs.call_seq;
    rhs.call_seq = 0;
    return *this;
  };
};

struct call_results {
  std::unique_ptr< std::stringstream > ss;
  std::unique_ptr< serialization::iarchive > p_ari;

  call_results() BOOST_NOEXCEPT : ss(), p_ari(){};

  call_results( call_results&& rhs ) BOOST_NOEXCEPT : ss( std::move( rhs.ss ) ), p_ari( std::move( rhs.p_ari ) ){};
  call_results& operator=( call_results&& rhs ) BOOST_NOEXCEPT {
    ss = std::move( rhs.ss );
    p_ari = std::move( rhs.p_ari );
    return *this;
  };
};

class remote_function {
public:
  std::string name;

  virtual void execute( serialization::iarchive& inputs, serialization::oarchive& outputs ) = 0;

  remote_function( const std::string& aName = "" );
  virtual ~remote_function();

  virtual std::size_t get_params_hash() const { return 0; };
  virtual std::string get_host() const { return "localhost"; };

protected:
  void publish();
  void unpublish();

  call_preparations prepare_for_call();
  call_results do_remote_call( call_preparations&& pre_data );
};

class dummy_remote_function : public remote_function {
public:
  std::size_t params_hash;

  void execute( serialization::iarchive& inputs, serialization::oarchive& outputs ){};

  dummy_remote_function( const std::string& aName = "", std::size_t aParamsHash = 0 )
      : remote_function( aName ), params_hash( aParamsHash ){};

  std::size_t get_params_hash() const { return params_hash; };
};
};
};
};


#endif
