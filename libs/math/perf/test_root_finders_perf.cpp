
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

#include <cmath>
#include "ReaK/math/root_finders/bisection_method.hpp"
#include "ReaK/math/root_finders/secant_method.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

using std::cos;
using std::exp;
using std::log;
using std::pow;
using std::sin;

int iteration_count;

double func1(double x) {
  ++iteration_count;
  return 4.0 * cos(x) - exp(x);
};

double func2(double x) {
  ++iteration_count;
  double result = 0.0;
  for (int i = 1; i < 11; ++i) {
    double t = 0.1 * i;
    result += exp(x * t) - exp(5.0 * t);
  };
  return result;
};

double func3(double x) {
  ++iteration_count;
  return 2.0 * x * exp(-20.0) + 1.0 - 2.0 * exp(-20.0 * x);
};

double func4(double x) {
  ++iteration_count;
  return exp(1.0 / x - 25.0) - 1.0;
};

double func5(double x) {
  ++iteration_count;
  double result = 1e-8 * (x - 1.0);
  for (int i = 1; i < 11; ++i) {
    result *= x * x + x + i;
  };
  return result;
};

double func6(double x) {
  ++iteration_count;
  return 1e10 * pow(x, 1.0 / x) - 1.0;
};

double func7(double x) {
  ++iteration_count;
  return pow(x, 20.0) - 1.0;
};

double func8(double x) {
  ++iteration_count;
  return exp(21000.0 / x) / (1.11e11 * x * x) - 1.0;
};

double func9(double x) {
  ++iteration_count;
  return 1.0 / x + log(x) - 100.0;
};

double func10(double x) {
  ++iteration_count;
  return exp(exp(x)) - exp(exp(1.0));
};

double func11(double x) {
  ++iteration_count;
  return sin(0.01 / x) - 0.01;
};

int main() {

  using func_ptr = double (*)(double);

  std::vector<func_ptr> funcs;
  std::vector<double> lows;
  std::vector<double> his;

  funcs.push_back(func1);
  lows.push_back(0.0);
  his.push_back(1.5);
  funcs.push_back(func1);
  lows.push_back(-1.0);
  his.push_back(3.0);
  funcs.push_back(func1);
  lows.push_back(-1.5);
  his.push_back(6.0);

  funcs.push_back(func2);
  lows.push_back(4.0);
  his.push_back(6.5);
  funcs.push_back(func2);
  lows.push_back(2.0);
  his.push_back(8.0);
  funcs.push_back(func2);
  lows.push_back(0.0);
  his.push_back(15.0);
  funcs.push_back(func2);
  lows.push_back(-5.0);
  his.push_back(25.0);

  funcs.push_back(func3);
  lows.push_back(0.0);
  his.push_back(1.0);
  funcs.push_back(func3);
  lows.push_back(-0.1);
  his.push_back(1.5);
  funcs.push_back(func3);
  lows.push_back(-0.5);
  his.push_back(2.0);
  funcs.push_back(func3);
  lows.push_back(-1.0);
  his.push_back(4.0);

  funcs.push_back(func4);
  lows.push_back(0.035);
  his.push_back(0.05);
  funcs.push_back(func4);
  lows.push_back(0.03);
  his.push_back(0.09);
  funcs.push_back(func4);
  lows.push_back(0.025);
  his.push_back(0.5);
  funcs.push_back(func4);
  lows.push_back(0.02);
  his.push_back(1.0);

  funcs.push_back(func5);
  lows.push_back(0.9);
  his.push_back(1.1);
  funcs.push_back(func5);
  lows.push_back(0.5);
  his.push_back(1.5);
  funcs.push_back(func5);
  lows.push_back(-5.0);
  his.push_back(5.0);
  funcs.push_back(func5);
  lows.push_back(5.0);
  his.push_back(10.0);

  funcs.push_back(func6);
  lows.push_back(0.095);
  his.push_back(1.0);
  funcs.push_back(func6);
  lows.push_back(0.075);
  his.push_back(0.15);
  funcs.push_back(func6);
  lows.push_back(0.08);
  his.push_back(0.5);
  funcs.push_back(func6);
  lows.push_back(0.05);
  his.push_back(0.2);

  funcs.push_back(func7);
  lows.push_back(0.9);
  his.push_back(1.05);
  funcs.push_back(func7);
  lows.push_back(0.7);
  his.push_back(1.2);
  funcs.push_back(func7);
  lows.push_back(0.0);
  his.push_back(2.5);
  funcs.push_back(func7);
  lows.push_back(-0.5);
  his.push_back(5.0);

  funcs.push_back(func8);
  lows.push_back(550.0);
  his.push_back(560.0);
  funcs.push_back(func8);
  lows.push_back(525.0);
  his.push_back(590.0);
  funcs.push_back(func8);
  lows.push_back(400.0);
  his.push_back(600.0);
  funcs.push_back(func8);
  lows.push_back(350.0);
  his.push_back(850.0);

  funcs.push_back(func9);
  lows.push_back(0.005);
  his.push_back(0.02);
  funcs.push_back(func9);
  lows.push_back(0.001);
  his.push_back(0.05);
  funcs.push_back(func9);
  lows.push_back(0.0001);
  his.push_back(0.1);
  funcs.push_back(func9);
  lows.push_back(0.001);
  his.push_back(100.0);

  funcs.push_back(func10);
  lows.push_back(0.0);
  his.push_back(2.0);
  funcs.push_back(func10);
  lows.push_back(-4.0);
  his.push_back(2.0);
  funcs.push_back(func10);
  lows.push_back(-10.0);
  his.push_back(3.0);
  funcs.push_back(func10);
  lows.push_back(0.5);
  his.push_back(3.5);

  funcs.push_back(func11);
  lows.push_back(0.5);
  his.push_back(2.0);
  funcs.push_back(func11);
  lows.push_back(0.2);
  his.push_back(6.0);
  funcs.push_back(func11);
  lows.push_back(0.05);
  his.push_back(20.0);
  funcs.push_back(func11);
  lows.push_back(0.004);
  his.push_back(200.0);

  for (unsigned int i = 0; i < funcs.size(); ++i) {
    std::cout << " " << std::setw(5) << i << " | ";

    double rel_tol = 1e-10 / (his[i] - lows[i]);

    double l = lows[i];
    double h = his[i];
    iteration_count = 0;
    ReaK::bisection_method(l, h, funcs[i], rel_tol);
    std::cout << "Bisection found solution = " << l << " to " << h << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    l = lows[i];
    h = his[i];
    iteration_count = 0;
    double sec_sol = ReaK::secant_method(l, h, funcs[i], rel_tol);
    std::cout << "Secant found solution = " << sec_sol << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    l = lows[i];
    h = his[i];
    iteration_count = 0;
    ReaK::illinois_method(l, h, funcs[i], rel_tol);
    std::cout << "Illinois found solution = " << l << " to " << h << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    l = lows[i];
    h = his[i];
    iteration_count = 0;
    ReaK::ford3_method(l, h, funcs[i], rel_tol);
    std::cout << "Ford3 found solution = " << l << " to " << h << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    l = lows[i];
    h = his[i];
    iteration_count = 0;
    ReaK::brent_method(l, h, funcs[i], rel_tol);
    std::cout << "Brent found solution = " << l << " to " << h << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    l = lows[i];
    h = his[i];
    iteration_count = 0;
    ReaK::ridders_method(l, h, funcs[i], rel_tol);
    std::cout << "Ridders found solution = " << l << " to " << h << std::endl;
    std::cout << std::setw(6) << iteration_count << std::endl;

    std::cout << std::endl;
  };
};
