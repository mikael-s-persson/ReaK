/**
 * \file chrono_incl.hpp
 *
 * This library contains a few useful macros and inclusions to handle the use of chrono library.
 * This library should be included instead of either Boost.Chrono libraries or C++11 standard
 * chrono libraries. This header takes care of figuring out which chrono library is appropriate
 * and imports the relevant objects into the ReaKaux namespace.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date April 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_CHRONO_INCL_HPP
#define REAK_CHRONO_INCL_HPP

#include "defs.hpp"



#ifndef BOOST_NO_CXX11_HDR_CHRONO

// must use standard chrono library because most supported versions of boost, up to 1.48 are broken for C++11 under GCC 4.7 or higher.
#include <chrono>
#include <ratio>

namespace ReaKaux {
  
  using std::ratio;
  
  namespace chrono {
    
    using std::chrono::duration;
    using std::chrono::time_point;
    using std::chrono::duration_cast;
    using std::chrono::time_point_cast;
    
    using std::chrono::nanoseconds;
    using std::chrono::microseconds;
    using std::chrono::milliseconds;
    using std::chrono::seconds;
    using std::chrono::minutes;
    using std::chrono::hours;
    
    using std::chrono::system_clock;
    using std::chrono::steady_clock;
    using std::chrono::high_resolution_clock;
    
  };

};

#else

#include <boost/version.hpp>

#if (BOOST_VERSION >= 104700)

// must use the Boost.Chrono library, because there was some indication that the gnu implementation is broken.
#include <boost/chrono/chrono.hpp>
#include <boost/ratio/ratio.hpp>

namespace ReaKaux {
  
  using boost::ratio;
  
  namespace chrono {
    
    using boost::chrono::duration;
    using boost::chrono::time_point;
    using boost::chrono::duration_cast;
    using boost::chrono::time_point_cast;
    
    using boost::chrono::nanoseconds;
    using boost::chrono::microseconds;
    using boost::chrono::milliseconds;
    using boost::chrono::seconds;
    using boost::chrono::minutes;
    using boost::chrono::hours;
    
    using boost::chrono::system_clock;
    using boost::chrono::steady_clock;
    using boost::chrono::high_resolution_clock;
    
  };

};

#endif

#endif




#endif





