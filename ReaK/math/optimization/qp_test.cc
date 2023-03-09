
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

#include "ReaK/math/optimization/line_search.h"

#include "ReaK/math/optimization/finite_diff_jacobians.h"

#include "ReaK/math/optimization/mehrotra_method.h"
#include "ReaK/math/optimization/quadratic_programs.h"

#include "ReaK/math/lin_alg/mat_svd_method.h"

#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

using namespace ReaK;

int main() {

  std::vector<mat<double, mat_structure::symmetric>> Gs;
  std::vector<mat<double, mat_structure::rectangular>> As;
  std::vector<vect_n<double>> bs;
  std::vector<vect_n<double>> cs;
  std::vector<vect_n<double>> xs;

  Gs.emplace_back(4.0, 1.5, 5.0, 8.0, 2.0, 12.0);
  cs.emplace_back(1.0, 1.0, 0.5);
  bs.emplace_back(1, 0.0);
  mat<double, mat_structure::rectangular> A(1, 3);
  A(0, 0) = -2.0;
  A(0, 1) = 1.0;
  A(0, 2) = 4.0;
  As.push_back(A);
  xs.emplace_back(1.0, 1.0, 1.0);

  Gs.emplace_back(6.0, 2.0, 1.0, 5.0, 2.0, 4.0);
  cs.emplace_back(-8.0, -3.0, -3.0);
  bs.emplace_back(vect<double, 2>(3.0, 0.0));
  A = mat<double, mat_structure::rectangular>(2, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 0.0;
  A(0, 2) = 1.0;
  A(1, 0) = 0.0;
  A(1, 1) = 1.0;
  A(1, 2) = 1.0;
  As.push_back(A);
  xs.emplace_back(1.0, 1.0, 1.0);

  vect_n<double> x;
  for (std::size_t i = 0; i < Gs.size(); ++i) {
    std::cout << "*************************************************************"
                 "***************************"
              << std::endl;
    std::cout << "Solving the quadratic problem #" << i
              << " with:" << std::endl;
    std::cout << "  G = " << Gs[i] << std::endl;
    std::cout << "  c = " << cs[i] << std::endl;
    std::cout << "  A = " << As[i] << std::endl;
    std::cout << "  b = " << bs[i] << std::endl << std::endl;

    x = xs[i];
    try {
      optim::null_space_QP_method(As[i], bs[i], Gs[i], cs[i], x, 1e-8);
      std::cout << "  Null-Space QP method gives:\n"
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    } catch (std::exception& e) {
      std::cout << "  Null-Space QP method failed with error: " << e.what()
                << std::endl
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    };

    x = xs[i];
    try {
      optim::projected_CG_method(As[i], bs[i], Gs[i], cs[i], x, 100, 1e-8);
      std::cout << "  Projected CG method gives:\n"
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    } catch (std::exception& e) {
      std::cout << "  Projected CG method failed with error: " << e.what()
                << std::endl
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    };

    x = xs[i];
    try {
      optim::mehrotra_QP_method(mat<double, mat_structure::rectangular>(0, 3),
                                vect_n<double>(0), Gs[i], cs[i], As[i], bs[i],
                                x, 100, 1e-8);
      std::cout << "  Mehrotra QP method gives:\n"
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    } catch (std::exception& e) {
      std::cout << "  Mehrotra QP method failed with error: " << e.what()
                << std::endl
                << "    x = " << x
                << " with |Ax - b| = " << norm_2(As[i] * x - bs[i]) << "\n"
                << "    and with xGx + cx = "
                << (0.5 * (x * Gs[i]) * x + cs[i] * x)
                << " so A(Gx + c) = " << (As[i] * (Gs[i] * x + cs[i]))
                << std::endl;
    };

    std::cout << "*************************************************************"
                 "***************************"
              << std::endl;
  };

  return 0;
};
