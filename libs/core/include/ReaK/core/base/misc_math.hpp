/**
 * \file misc_math.hpp
 *
 * This library contains a few useful basic math functions.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date May 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_MISC_MATH_HPP
#define REAK_MISC_MATH_HPP

#include <cmath>


namespace ReaK {

namespace math {


inline std::size_t highest_pow2(std::size_t N) {
  std::size_t temp;
  for(std::size_t shift = 1; ((shift < 8 * sizeof(std::size_t)) && (temp = N >> shift)); shift <<= 1)
    N |= temp;
  return N & ~(N >> 1);
};

inline std::size_t highest_set_bit(std::size_t N) {
  std::size_t temp = 0;
  for(std::size_t shift = sizeof(std::size_t) * 4; (shift && (N != 1)); shift >>= 1) {
    if(N >> shift) {
      temp |= shift;
      N >>= shift;
    };
  };
  return temp;
};




};

};


#endif





