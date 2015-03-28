/**
 * \file function_incl.hpp
 *
 * This library contains a inclusions to handle the use of std::function / boost::function.
 * This library should be included instead of either Boost.Function libraries or C++11 standard
 * functional library. This header takes care of figuring out which function library is appropriate
 * and imports the relevant objects into the ReaKaux namespace.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date July 2014
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

#ifndef REAK_FUNCTION_INCL_HPP
#define REAK_FUNCTION_INCL_HPP

#include "defs.hpp"


#ifndef BOOST_NO_CXX11_HDR_FUNCTIONAL

// Use standard functional library things:

#include <functional>

namespace ReaKaux {

using ::std::bind;
using ::std::mem_fn;
using ::std::is_bind_expression;
using ::std::is_placeholder;

namespace placeholders = ::std::placeholders;

using ::std::function;
using ::std::bad_function_call;

using ::std::reference_wrapper;
using ::std::cref;
using ::std::ref;

using ::std::hash;
};

#else

// Use Boost functional library things:

#include <boost/ref.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/functional/hash.hpp>

namespace ReaKaux {

using ::boost::bind;
using ::boost::mem_fn;
using ::boost::is_bind_expression;
using ::boost::is_placeholder;

namespace placeholders {
using ::_1;
using ::_2;
using ::_3;
using ::_4;
using ::_5;
using ::_6;
using ::_7;
using ::_8;
using ::_9;
};

using ::boost::function;
using ::boost::bad_function_call;

using ::boost::reference_wrapper;
using ::boost::cref;
using ::boost::ref;

using ::boost::hash;
};

#endif


#endif
