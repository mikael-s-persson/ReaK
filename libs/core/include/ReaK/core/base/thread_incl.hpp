/**
 * \file thread_incl.hpp
 *
 * This library contains a few useful macros and inclusions to handle the use of threads.
 * This library should be included instead of either Boost.Thread libraries or C++11 standard
 * thread libraries. This header takes care of figuring out which thread library is appropriate
 * and imports the relevant objects into the ReaKaux namespace.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date March 2012
 */

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

#ifndef REAK_THREAD_INCL_HPP
#define REAK_THREAD_INCL_HPP

#include "defs.hpp"

#include "function_incl.hpp"

//  C++11 thread features in GCC 4.7.0 and later
//  NOTE: This does not include futures and the notify-all-at-exit, because they are not fully supported.
// #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || !defined(__GXX_EXPERIMENTAL_CXX0X__)

#if !defined( BOOST_NO_CXX11_HDR_THREAD )

// This is a bit nasty...
#if __GNUC__
#ifndef _GLIBCXX_USE_SCHED_YIELD
#define _GLIBCXX_USE_SCHED_YIELD
#endif
#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif
#endif

// must use standard thread library because most supported versions of boost, up to 1.48 are broken for C++11 under GCC
// 4.7 or higher.
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <exception>


namespace ReaKaux {

// thread:

using ::std::thread;
// using ::std::swap;

namespace this_thread {

using ::std::this_thread::get_id;
using ::std::this_thread::yield;
using ::std::this_thread::sleep_until;
using ::std::this_thread::sleep_for;
};

// mutex:

using ::std::mutex;
using ::std::recursive_mutex;
using ::std::timed_mutex;
using ::std::recursive_timed_mutex;

using ::std::defer_lock_t;
using ::std::try_to_lock_t;
using ::std::adopt_lock_t;

using ::std::defer_lock;
using ::std::try_to_lock;
using ::std::adopt_lock;

using ::std::lock_guard;
using ::std::unique_lock;

// using ::std::swap;

//   using ::std::try_lock;
//   using ::std::lock;

using ::std::once_flag;

using ::std::call_once;

// condition_variable:

using ::std::condition_variable;
using ::std::condition_variable_any;

// using ::std::notify_all_at_thread_exit;

using ::std::exception_ptr;
using ::std::make_exception_ptr;
};

#else


// must use the Boost.Thread library, because there was some indication that the gnu implementation is broken.
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/once.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/exception_ptr.hpp>
#include <boost/version.hpp>


namespace ReaKaux {

// thread:

using ::boost::thread;
// using std::swap;

namespace this_thread {

using ::boost::this_thread::get_id;
using ::boost::this_thread::yield;
#if( BOOST_VERSION >= 105000 )
using ::boost::this_thread::sleep_until;
using ::boost::this_thread::sleep_for;
#endif
};

// mutex:

using ::boost::mutex;
using ::boost::recursive_mutex;
using ::boost::timed_mutex;
using ::boost::recursive_timed_mutex;

using ::boost::defer_lock_t;
using ::boost::try_to_lock_t;
using ::boost::adopt_lock_t;

using ::boost::defer_lock;
using ::boost::try_to_lock;
using ::boost::adopt_lock;

using ::boost::lock_guard;
using ::boost::unique_lock;

// using ::std::swap;

//   using ::boost::try_lock;
//   using ::boost::lock;

using ::boost::once_flag;

using ::boost::call_once;

// condition_variable:

using ::boost::condition_variable;
using ::boost::condition_variable_any;

// using ::boost::notify_all_at_thread_exit;

using ::boost::exception_ptr;

//   using ::boost::make_exception_ptr;
template < class T >
exception_ptr make_exception_ptr( T const& e ) {
  return ::boost::copy_exception( e );
};
};

#endif


#ifndef BOOST_NO_CXX11_HDR_FUTURE

// must use standard thread library because most supported versions of boost, up to 1.48 are broken for C++11 under GCC
// 4.7 or higher.
#include <future>

namespace ReaKaux {

using ::std::future;
using ::std::future_error;
using ::std::promise;
};

#else


// must use the Boost library.
#define BOOST_THREAD_PROVIDES_FUTURE
#include <boost/thread/future.hpp>

namespace ReaKaux {

using ::boost::future;
using ::boost::future_error;
using ::boost::promise;
};

#endif


#endif
