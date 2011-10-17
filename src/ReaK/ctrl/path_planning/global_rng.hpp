/**
 * \file global_rng.hpp
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

#ifndef REAK_GLOBAL_RNG_HPP
#define REAK_GLOBAL_RNG_HPP


#include "base/defs.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>

#include <ctime>

namespace ReaK {

namespace pp {


typedef boost::minstd_rand global_rng_type;

inline global_rng_type& get_global_rng() {
  static global_rng_type instance = global_rng_type(static_cast<unsigned int>(std::time(NULL)));
  return instance;
};


};

};

#endif








