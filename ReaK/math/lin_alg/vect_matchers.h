/// \file vect_matchers.h
///
/// This library provides several googletest matchers for the various vector comparisons.
///
/// \author Sven Mikael Persson <mikael.s.persson@gmail.com>
/// \date March 2023

/*
 *    Copyright 2023 Sven Mikael Persson
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

#ifndef REAK_MATH_LIN_ALG_VECT_MATCHERS_H_
#define REAK_MATH_LIN_ALG_VECT_MATCHERS_H_

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_traits.h"
#include "gmock/gmock.h"

namespace ReaK {
namespace testing {

/// This function template computes the element-wise comparison of two vectors.
MATCHER_P2(VectorIsNear, expected, tolerance, "") {
  if (arg.size() != expected.size()) {
    *result_listener << "\nExpected vector =\n  " << expected
                     << "\nActual vector =\n  " << arg
                     << "\nSizes do not match, expected " << expected.size()
                     << " but got " << arg.size();
    return false;
  }
  for (int i = 0; i < expected.size(); ++i) {
    if (!::testing::ExplainMatchResult(
            ::testing::DoubleNear(expected[i], tolerance), arg[i],
            result_listener)) {
      *result_listener << "\nExpected vector =\n  " << expected
                       << "\nActual vector =\n  " << arg << "\nDifference =\n"
                       << (arg - expected) << "\nTolerance = " << tolerance
                       << "\n";
      return false;
    }
  }
  return true;
}

/// This function template computes the element-wise comparison of two vectors.
MATCHER_P(VectorIsZero, tolerance, "") {
  for (int i = 0; i < arg.size(); ++i) {
    if (!::testing::ExplainMatchResult(::testing::DoubleNear(0.0, tolerance),
                                       arg[i], result_listener)) {
      *result_listener << "\nExpected zero vector"
                       << "\nActual vector =\n  " << arg
                       << "\nTolerance = " << tolerance << "\n";
      return false;
    }
  }
  return true;
}

/// This function template computes the element-wise comparison of two vectors.
MATCHER_P2(VectorIsGreater, expected, tolerance, "") {
  if (arg.size() != expected.size()) {
    *result_listener << "\nComparing to vector =\n  " << expected
                     << "\nActual vector =\n  " << arg
                     << "\nSizes do not match, expected " << expected.size()
                     << " but got " << arg.size();
    return false;
  }
  for (int i = 0; i < expected.size(); ++i) {
    if (!::testing::ExplainMatchResult(::testing::Gt(expected[i] - tolerance),
                                       arg[i], result_listener)) {
      *result_listener << "\nComparing to vector =\n  " << expected
                       << "\nActual vector =\n  " << arg << "\nDifference =\n"
                       << (arg - expected) << "\nTolerance = " << tolerance
                       << "\n";
      return false;
    }
  }
  return true;
}

/// This function template computes the element-wise comparison of two vectors.
MATCHER_P2(VectorIsLess, expected, tolerance, "") {
  if (arg.size() != expected.size()) {
    *result_listener << "\nComparing to vector =\n  " << expected
                     << "\nActual vector =\n  " << arg
                     << "\nSizes do not match, expected " << expected.size()
                     << " but got " << arg.size();
    return false;
  }
  for (int i = 0; i < expected.size(); ++i) {
    if (!::testing::ExplainMatchResult(::testing::Lt(expected[i] + tolerance),
                                       arg[i], result_listener)) {
      *result_listener << "\nComparing to vector =\n  " << expected
                       << "\nActual vector =\n  " << arg << "\nDifference =\n"
                       << (arg - expected) << "\nTolerance = " << tolerance
                       << "\n";
      return false;
    }
  }
  return true;
}

}  // namespace testing
}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_VECT_MATCHERS_H_
