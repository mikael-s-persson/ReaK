/**
 * \file defs.h
 *
 * This library contains a few useful macros for the ReaK platform. Mainly including debugging outputs.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date february 2010
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

#ifndef REAK_CORE_BASE_DEFS_H_
#define REAK_CORE_BASE_DEFS_H_

#include <cmath>
#include <iostream>
#include <string>

#ifndef M_PI
#define M_PI 3.1415926535898
#endif

#define RK_ORDER_LITTLE_ENDIAN 1
#define RK_ORDER_BIG_ENDIAN 2
#define RK_ORDER_PDP_ENDIAN 3

#ifdef __GNUC__

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define RK_BYTE_ORDER RK_ORDER_LITTLE_ENDIAN
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define RK_BYTE_ORDER RK_ORDER_BIG_ENDIAN
#else
#define RK_BYTE_ORDER RK_ORDER_PDP_ENDIAN
#endif

#endif  // __GNUC__

/**
 * This MACRO is used to signify that a declared variable is not used, intentionally.
 * This creates a no-op that uses the variable X, and thus, avoid annoying compiler warnings
 * such as "parameter X is never used" or "variable X is set but never used".
 */
#define RK_UNUSED(X) \
  { (void)X; }

#endif  // REAK_CORE_BASE_DEFS_H_
