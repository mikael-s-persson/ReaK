
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

#include "ReaK/math/root_finders/bisection_method.h"
#include "ReaK/math/root_finders/secant_method.h"

#include <cmath>
#include <vector>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

double func1(double x) {
  return 4.0 * std::cos(x) - std::exp(x);
}

double func2(double x) {
  double result = 0.0;
  for (int i = 1; i < 11; ++i) {
    double t = 0.1 * i;
    result += std::exp(x * t) - std::exp(5.0 * t);
  }
  return result;
}

double func3(double x) {
  return 2.0 * x * std::exp(-20.0) + 1.0 - 2.0 * std::exp(-20.0 * x);
}

double func4(double x) {
  return std::exp(1.0 / x - 25.0) - 1.0;
}

TEST(RootFinding, Tests) {
  using func_ptr = double (*)(double);

  std::vector<func_ptr> funcs;
  std::vector<double> lows;
  std::vector<double> his;
  std::vector<double> sols;

  funcs.push_back(func1);
  lows.push_back(0.0);
  his.push_back(1.5);
  sols.push_back(0.904788);
  funcs.push_back(func2);
  lows.push_back(4.0);
  his.push_back(6.5);
  sols.push_back(5.0);
  funcs.push_back(func3);
  lows.push_back(0.0);
  his.push_back(1.0);
  sols.push_back(0.03465736);
  funcs.push_back(func4);
  lows.push_back(0.03);
  his.push_back(0.09);
  sols.push_back(0.04);

  for (unsigned int i = 0; i < funcs.size(); ++i) {
    double rel_tol = 1e-10 / (his[i] - lows[i]);

    double l = lows[i];
    double h = his[i];
    EXPECT_NO_THROW(bisection_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);

    l = lows[i];
    h = his[i];
    EXPECT_NO_THROW(l = secant_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);

    l = lows[i];
    h = his[i];
    EXPECT_NO_THROW(illinois_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);

    l = lows[i];
    h = his[i];
    EXPECT_NO_THROW(ford3_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);

    l = lows[i];
    h = his[i];
    EXPECT_NO_THROW(brent_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);

    l = lows[i];
    h = his[i];
    EXPECT_NO_THROW(ridders_method(l, h, funcs[i], rel_tol));
    EXPECT_NEAR(l, sols[i], 5e-5);
  }
}

}  // namespace
}  // namespace ReaK
