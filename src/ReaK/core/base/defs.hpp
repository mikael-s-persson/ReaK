/**
 * \file defs.hpp
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

#ifndef rk_defs_H
#define rk_defs_H

#ifdef WIN32
#if defined(_M_X64) || defined(__x86_64__)
#define RK_CALL
#else
#define RK_CALL __stdcall
#endif
#else
#ifdef __x86_64__
#define RK_CALL
#else
#define RK_CALL __attribute__((__stdcall__))
#endif
#endif
#define RK_EXTERN extern "C"

#ifndef RK_VERBOSITY
#define RK_VERBOSITY 5
#endif

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#define RK_ENABLE_CXX0X_FEATURES
#endif

#include <string>
#include <iostream>

/**
 * This function will convert an absolute path that goes through a folder named "ReaK" and
 * return a the path relative to that ReaK trunk folder.
 *
 * \param S Some string containing a full file/folder path.
 * \return string containing the file/folder path relative to the ReaK trunk folder.
 */
inline std::string RK_RELATIVE_PATH(const std::string& S) {
  for(unsigned int i=0;i + 5 < S.size();++i)
    if(S.substr(i,4) == "ReaK")
      return S.substr(i+4);
  return S;
};

/**
 * This MACRO expands into an output of at string as "ReaK/.../current_file.hpp:42" if put in 
 * file "current_file.hpp" at line 42.
 */
#define RK_HERE __FILE__ << ":" << __LINE__

/**
 * This MACRO outputs to std::cout an error message containing the filename, line number, and message X.
 */
#define RK_ERROR(X) std::cout << RK_RELATIVE_PATH(__FILE__) << ":" << __LINE__ << " Error: " << X << std::endl

/**
 * This MACRO outputs to std::cout a warning message containing the filename, line number, and message X.
 */
#define RK_WARNING(X) std::cout << RK_RELATIVE_PATH(__FILE__) << ":" << __LINE__ << " Warning: " << X << std::endl

/**
 * This MACRO outputs to std::cout a notice message containing the filename, line number, and message Y, 
 * only if the RK_VERBOSITY is set to lower higher than X.
 */
#define RK_NOTICE(X,Y) if(X <= RK_VERBOSITY)  std::cout << RK_RELATIVE_PATH(__FILE__) << ":" << __LINE__ << " " << Y << std::endl;

/**
 * This MACRO is used to signify that a declared variable is not used, intentionally. 
 * This creates a no-op that uses the variable X, and thus, avoid annoying compiler warnings 
 * such as "parameter X is never used" or "variable X is set but never used".
 */
#define RK_UNUSED(X) { (void)X; }


namespace ReaK {

template<bool> struct CompileTimeChecker
{
   CompileTimeChecker(...);
};
template<> struct CompileTimeChecker<false> { };
};
#define RK_CT_ASSERT(expr, msg) \
   {\
       class ERROR_##msg {}; \
       (void)sizeof(CompileTimeChecker<(expr) != 0>((ERROR_##msg())));\
   }


#endif





