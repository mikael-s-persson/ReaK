
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

#include "ReaK/math/optimization/mehrotra_method.h"
#include "ReaK/math/optimization/simplex_method.h"

#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

namespace ReaK::optim {
namespace {

// Switch from standard form of:
//  max c*x
//  A*x <= b
//  l <= x <= u
// to augmented (general / slack) form of:
//  max c*x
//  [A I] * [x s] = b
//  l <= x <= u
//  0 <= s
void simplex_method_on_augmented_form(
    const mat<double, mat_structure::rectangular>& A, const vect_n<double>& b,
    const vect_n<double>& c, vect_n<double>& x, const vect_n<double>& l,
    const vect_n<double>& u, double tolerance) {
  mat<double, mat_structure::rectangular> A_aug(
      A.get_row_count(), A.get_col_count() + A.get_row_count(), 0.0);
  sub(A_aug)(range(0, A.get_row_count()), range(0, A.get_col_count())) = A;
  sub(A_aug)(range(0, A.get_row_count()),
             range(A.get_col_count(), A.get_col_count() + A.get_row_count())) =
      mat_ident<double>(A.get_row_count());
  vect_n<double> x_aug(x.size() + b.size(), 0.0);
  sub(x_aug)[range(0, x.size())] = x;
  vect_n<double> l_aug(b.size() + x.size(), 0.0);
  sub(l_aug)[range(0, x.size())] = l;
  vect_n<double> u_aug(b.size() + x.size(),
                       std::numeric_limits<double>::infinity());
  sub(u_aug)[range(0, x.size())] = u;
  vect_n<double> c_aug(b.size() + x.size(), 0.0);
  sub(c_aug)[range(0, x.size())] = c;

  simplex_method(A_aug, b, c_aug, x_aug, l_aug, u_aug, tolerance);

  x = sub(x_aug)[range(0, x.size())];
}

// Switch from standard form of:
//  max c*x
//  A*x <= b
//  0 <= x
// to augmented (general / slack) form of:
//  max c*x
//  [A I] * [x s] = b
//  0 <= x
//  0 <= s
void mehrotra_method_on_augmented_form(
    const mat<double, mat_structure::rectangular>& A, const vect_n<double>& b,
    const vect_n<double>& c, vect_n<double>& x, int max_iter,
    double tolerance) {
  mat<double, mat_structure::rectangular> A_aug(
      A.get_row_count(), A.get_col_count() + A.get_row_count(), 0.0);
  sub(A_aug)(range(0, A.get_row_count()), range(0, A.get_col_count())) = A;
  sub(A_aug)(range(0, A.get_row_count()),
             range(A.get_col_count(), A.get_col_count() + A.get_row_count())) =
      mat_ident<double>(A.get_row_count());
  vect_n<double> x_aug(x.size() + b.size(), 0.0);
  sub(x_aug)[range(0, x.size())] = x;
  vect_n<double> c_aug(b.size() + x.size(), 0.0);
  sub(c_aug)[range(0, x.size())] = c;

  mehrotra_method(A_aug, b, c_aug, x_aug, max_iter, tolerance);

  x = sub(x_aug)[range(0, x.size())];
}

TEST(SimplexMethod, BasicTests) {
  const double tolerance = 1e-6;

  vect_n<double> c(2);
  mat<double, mat_structure::rectangular> A(1, 2);
  vect_n<double> b(1);
  vect_n<double> x(c.size(), 0.2);
  vect_n<double> l(c.size(), 0.0);
  vect_n<double> u(c.size(), 0.75);
  // Favor maximum x coord.
  c[0] = 1.0;
  c[1] = 0.001;  // tie-break
  // Within triangle near origin.
  A(0, 0) = 1.0;
  A(0, 1) = 1.0;
  b[0] = 1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(
      simplex_method_on_augmented_form(A, b, c, x, l, u, tolerance));
  EXPECT_NEAR(x[0], 0.75, tolerance);
  EXPECT_NEAR(x[1], 0.25, tolerance);

  // Favor maximum y coord.
  c[0] = 0.001;  // tie-break
  c[1] = 1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(
      simplex_method_on_augmented_form(A, b, c, x, l, u, tolerance));
  EXPECT_NEAR(x[0], 0.25, tolerance);
  EXPECT_NEAR(x[1], 0.75, tolerance);

  // Oppose constraint.
  c[0] = -1.0;
  c[1] = -1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(
      simplex_method_on_augmented_form(A, b, c, x, l, u, tolerance));
  EXPECT_NEAR(x[0], 0.0, tolerance);
  EXPECT_NEAR(x[1], 0.0, tolerance);

  // Aligned to constraint.
  c[0] = 1.0;
  c[1] = 1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(
      simplex_method_on_augmented_form(A, b, c, x, l, u, tolerance));
  EXPECT_NEAR(x * c, 1.0, tolerance);
}

TEST(MehrotraMethod, BasicTests) {
  const double tolerance = 1e-6;

  vect_n<double> c(2);
  mat<double, mat_structure::rectangular> A(1, 2);
  vect_n<double> b(1);
  vect_n<double> x(c.size(), 0.2);
  // Favor maximum x coord.
  c[0] = 1.0;
  c[1] = 0.0;
  // Within triangle near origin.
  A(0, 0) = 1.0;
  A(0, 1) = 1.0;
  b[0] = 1.0;

  EXPECT_NO_THROW(mehrotra_method_on_augmented_form(A, b, -c, x, 100,
                                                    tolerance * tolerance));
  EXPECT_NEAR(x[0], 1.0, std::sqrt(tolerance));
  EXPECT_NEAR(x[1], 0.0, std::sqrt(tolerance));

  // Favor maximum y coord.
  c[0] = 0.0;
  c[1] = 1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(mehrotra_method_on_augmented_form(A, b, -c, x, 100,
                                                    tolerance * tolerance));
  EXPECT_NEAR(x[0], 0.0, std::sqrt(tolerance));
  EXPECT_NEAR(x[1], 1.0, std::sqrt(tolerance));

#if 0
  // TODO: MEHROTRA IS FAILING CASES BELOW 
  // Oppose constraint.
  c[0] = -1.0;
  c[1] = -1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(mehrotra_method_on_augmented_form(A, b, -c, x, 100, tolerance*tolerance));
  EXPECT_NEAR(x[0], 0.0, std::sqrt(tolerance));
  EXPECT_NEAR(x[1], 0.0, std::sqrt(tolerance));

  // Aligned to constraint.
  c[0] = 1.0;
  c[1] = 1.0;

  x[0] = 0.2;
  x[1] = 0.2;
  EXPECT_NO_THROW(mehrotra_method_on_augmented_form(A, b, -c, x, 100, tolerance*tolerance));
  EXPECT_NEAR(x * c, 1.0, std::sqrt(tolerance));
#endif
}

}  // namespace
}  // namespace ReaK::optim
