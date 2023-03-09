/**
 * \file global_rng.h
 *
 * This library provides defines the global random-number generator object.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_CORE_BASE_GLOBAL_RNG_H_
#define REAK_CORE_BASE_GLOBAL_RNG_H_

#include <random>

namespace ReaK {

/// This is the type of the global random-number generator used in ReaK's algorithms.
using global_rng_type = std::mt19937;

/**
 * This function returns the global (static) instance of the random-number generator (seeded at first use).
 * \return A reference to the global (static) instance of the random-number generator (seeded at first use).
 */
inline global_rng_type& get_global_rng() {
  static std::random_device rd;
  static auto instance = global_rng_type(rd());
  return instance;
}

}  // namespace ReaK

#endif  // REAK_CORE_BASE_GLOBAL_RNG_H_
