/**
 * \file atomic_incl.hpp
 *
 * This library contains imports of inclusions of atomics classes and operations depending on 
 * whether the standard library supports it or not (fallback is Boost).
 * This library should be included instead of either Boost.Atomics libraries or C++11 standard
 * atomic libraries. This header takes care of figuring out which atomic library is appropriate
 * and imports the relevant objects into the ReaKaux namespace.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date August 2014
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

#ifndef REAK_ATOMIC_INCL_HPP
#define REAK_ATOMIC_INCL_HPP

#include "defs.hpp"

//  NOTE: Currently, there is no definite and simple test for whether the <atomic> header exists (what's up with that Boost devs!!!)
#if 0

#include <atomic>

namespace ReaKaux {
  
  // atomics:
  
  using ::std::atomic;
  
//   using ::std::atomic_is_lock_free;
//   using ::std::atomic_store;
//   using ::std::atomic_store_explicit;
//   using ::std::atomic_load;
//   using ::std::atomic_load_explicit;
//   using ::std::atomic_exchange;
//   using ::std::atomic_exchange_explicit;
//   using ::std::atomic_compare_exchange_weak;
//   using ::std::atomic_compare_exchange_weak_explicit;
//   using ::std::atomic_compare_exchange_strong;
//   using ::std::atomic_compare_exchange_strong_explicit;
//   using ::std::atomic_fetch_add;
//   using ::std::atomic_fetch_add_explicit;
//   using ::std::atomic_fetch_sub;
//   using ::std::atomic_fetch_sub_explicit;
//   using ::std::atomic_fetch_and;
//   using ::std::atomic_fetch_and_explicit;
//   using ::std::atomic_fetch_or;
//   using ::std::atomic_fetch_or_explicit;
//   using ::std::atomic_fetch_xor;
//   using ::std::atomic_fetch_xor_explicit;
  
  using ::std::atomic_flag;
//   using ::std::atomic_flag_test_and_set;
//   using ::std::atomic_flag_test_and_set_explicit;
//   using ::std::atomic_flag_clear;
//   using ::std::atomic_flag_clear_explicit;
  
//   using ::std::atomic_init;
  
  using ::std::memory_order;
//   using ::std::kill_dependency;
  using ::std::atomic_thread_fence;
  using ::std::atomic_signal_fence;
  
};

#else

#include <boost/atomic.hpp>


namespace ReaKaux {
  
  // atomics:
  
  using ::boost::atomic;
  
//   using ::boost::atomic_is_lock_free;
//   using ::boost::atomic_store;
//   using ::boost::atomic_store_explicit;
//   using ::boost::atomic_load;
//   using ::boost::atomic_load_explicit;
//   using ::boost::atomic_exchange;
//   using ::boost::atomic_exchange_explicit;
//   using ::boost::atomic_compare_exchange_weak;
//   using ::boost::atomic_compare_exchange_weak_explicit;
//   using ::boost::atomic_compare_exchange_strong;
//   using ::boost::atomic_compare_exchange_strong_explicit;
//   using ::boost::atomic_fetch_add;
//   using ::boost::atomic_fetch_add_explicit;
//   using ::boost::atomic_fetch_sub;
//   using ::boost::atomic_fetch_sub_explicit;
//   using ::boost::atomic_fetch_and;
//   using ::boost::atomic_fetch_and_explicit;
//   using ::boost::atomic_fetch_or;
//   using ::boost::atomic_fetch_or_explicit;
//   using ::boost::atomic_fetch_xor;
//   using ::boost::atomic_fetch_xor_explicit;
  
  using ::boost::atomic_flag;
//   using ::boost::atomic_flag_test_and_set;
//   using ::boost::atomic_flag_test_and_set_explicit;
//   using ::boost::atomic_flag_clear;
//   using ::boost::atomic_flag_clear_explicit;
  
//   using ::boost::atomic_init;
  
  using ::boost::memory_order;
//   using ::boost::kill_dependency;
  using ::boost::atomic_thread_fence;
  using ::boost::atomic_signal_fence;
  
};

#endif




#endif





