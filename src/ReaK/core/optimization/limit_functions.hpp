/**
 * \file limit_functions.hpp
 *
 * The following library provides a number of limit functions for vector-spaces.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_LIMIT_FUNCTIONS_HPP
#define REAK_LIMIT_FUNCTIONS_HPP

#include "base/defs.hpp"

#include "lin_alg/vect_alg.hpp"


namespace ReaK {
  
  
namespace optim {


/**
 * This limit functor is a no-op limiter which imposes no limits on the candidate step.
 * TEST PASSED
 */
struct no_limit_functor {
  /**
   * This function imposes no limits on the step proposed over the given vector.
   */
  template <typename Vector>
  void operator()(const Vector&, Vector&) const { };
};


/**
 * This function limits a proposed step size such that the resulting vector is within 
 * the bounds of a hyperbox defined by lower and upper bound vectors.
 * TEST PASSED
 * \tparam Vector A vector type.
 * \param x The current vector.
 * \param dx The proposed step.
 * \param l The lower bound vector.
 * \param u The upper bound vector.
 */
template <typename Vector>
typename boost::enable_if<
  is_writable_vector<Vector>,
void >::type box_limit_function(const Vector& x, Vector& dx, const Vector& l, const Vector& u) {
  typedef typename vect_traits<Vector>::size_type SizeType;
  for(SizeType i = 0; i < x.size(); ++i) {
    if( x[i] + dx[i] < l[i] )
      dx[i] = l[i] - x[i];
    if( x[i] + dx[i] > u[i] )
      dx[i] = u[i] - x[i];
  };
};


};

};


#endif

