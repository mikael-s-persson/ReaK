/**
 *\file py_fixes.hpp
 *
 * This header file declares a few fixes required to make Boost.Python work in ReaK.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date June 2012
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

#ifndef REAK_PY_FIXES_HPP
#define REAK_PY_FIXES_HPP


#include "defs.hpp"

#include <boost/version.hpp>
#include <boost/python.hpp>


#ifndef BOOST_NO_CXX11_SMART_PTR

#if (BOOST_VERSION < 105300)

namespace std {
  
  template <typename T>
  T* get_pointer(const std::shared_ptr<T>& s) {
    return s.get();
  };
  
};

#endif

#endif


#endif


