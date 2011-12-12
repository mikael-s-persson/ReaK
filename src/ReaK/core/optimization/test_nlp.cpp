
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
#include "line_search.hpp"

#include "finite_diff_jacobians.hpp"

#include "augmented_lagrangian_methods.hpp"
#include "sequential_qp_methods.hpp"
#include "nl_interior_points_methods.hpp"

#include "lin_alg/mat_svd_method.hpp"

#include <iostream>
#include <cmath>

#ifndef M_PI 
#define M_PI 3.14159265359
#endif

static int evalCount;
static int gradCount;

using namespace ReaK;

typedef double (*FunctionPtr)(const vect_n<double>&);
typedef vect_n<double> (*GradFunctionPtr)(const vect_n<double>&);
typedef void (*HessFunctionPtr)(mat<double,mat_structure::symmetric>&,const vect_n<double>&,double,const vect_n<double>&);
typedef vect_n<double> (*GFunctionPtr)(const vect_n<double>&);
typedef void (*GJacFunctionPtr)(mat<double,mat_structure::rectangular>&,const vect_n<double>&,const vect_n<double>&);
typedef vect_n<double> (*HFunctionPtr)(const vect_n<double>&);
typedef void (*HJacFunctionPtr)(mat<double,mat_structure::rectangular>&,const vect_n<double>&,const vect_n<double>&);



double p01_f(const vect_n<double>& x) {
  ++evalCount;
  return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0]);
};

vect_n<double> p01_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
  result[1] = 200.0 * (x[1] - x[0] * x[0]);
  return result;
};

void p01_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0;
  H(0,1) = -400.0 * x[0];
  H(1,1) = 200.0;
};

vect_n<double> p01_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p01_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p01_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = x[1] + 1.5;
  return result;
};

void p01_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = 0.0; J(0,1) = 1.0;
};

vect_n<double> p01_sol = vect_n<double>(vect<double,2>(1.0,1.0));
vect_n<double> p01_start = vect_n<double>(vect<double,2>(-2.0,1.0));
vect_n<double> p01_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p01_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));




FunctionPtr p02_f = p01_f;
GradFunctionPtr p02_grad = p01_grad;
HessFunctionPtr p02_H = p01_H;
GFunctionPtr p02_g = p01_g;
GJacFunctionPtr p02_g_jac = p01_g_jac;

vect_n<double> p02_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = x[1] - 1.5;
  return result;
};

void p02_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = 0.0; J(0,1) = 1.0;
};

vect_n<double> p02_sol = vect_n<double>(vect<double,2>(-1.2210274,1.5));
vect_n<double> p02_start = vect_n<double>(vect<double,2>(-2.0,1.0));
vect_n<double> p02_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p02_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));



double p03_f(const vect_n<double>& x) {
  ++evalCount;
  return x[1] + 1.0e-5 * (x[1] - x[0]) * (x[1] - x[0]);
};

vect_n<double> p03_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = -2.0e-5 * (x[1] - x[0]);
  result[1] = 1.0 - result[0];
  return result;
};

void p03_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 2.0e-5;
  H(0,1) = -2.0e-5;
  H(1,1) = 2.0;
};

vect_n<double> p03_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p03_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p03_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = x[1];
  return result;
};

void p03_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = 0.0; J(0,1) = 1.0;
};

vect_n<double> p03_sol = vect_n<double>(vect<double,2>(0.0,0.0));
vect_n<double> p03_start = vect_n<double>(vect<double,2>(10.0,1.0));
vect_n<double> p03_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p03_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));




double p04_f(const vect_n<double>& x) {
  ++evalCount;
  return (x[0] + 1.0) * (x[0] + 1.0) * (x[0] + 1.0) / 3.0 + x[1]; 
};

vect_n<double> p04_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = (x[0] + 1.0) * (x[0] + 1.0);
  result[1] = 1.0;
  return result;
};

void p04_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 2.0 * (x[0] + 1.0);
  H(0,1) = 0.0;
  H(1,1) = 0.0;
};

vect_n<double> p04_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p04_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p04_h(const vect_n<double>& x) {
  vect_n<double> result(2);
  result[0] = x[0] - 1.0;
  result[1] = x[1];
  return result;
};

void p04_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(2);
  J(0,0) = 1.0; J(0,1) = 0.0;
  J(1,0) = 0.0; J(1,1) = 1.0;
};

vect_n<double> p04_sol = vect_n<double>(vect<double,2>(1.0,0.0));
vect_n<double> p04_start = vect_n<double>(vect<double,2>(1.125,0.125));
vect_n<double> p04_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p04_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));




double p05_f(const vect_n<double>& x) {
  ++evalCount;
  return sin(x[0] + x[1]) + (x[0] - x[1]) * (x[0] - x[1]) - 1.5 * x[0] + 2.5 * x[1] + 1.0;      
};

vect_n<double> p05_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  double v1 = cos(x[0] + x[1]);
  double v2 = 2.0 * (x[0] - x[1]);
  result[0] = v1 + v2 - 1.5; 
  result[1] = v1 - v2 + 2.5; 
  return result;
};

void p05_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  double v1 = -sin(x[0] + x[1]);
  H(0,0) = v1 + 2.0;
  H(0,1) = v1 - 2.0;
  H(1,1) = v1 + 2.0;
};

vect_n<double> p05_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p05_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p05_h(const vect_n<double>& x) {
  vect_n<double> result(4);
  result[0] = x[0] + 1.5;
  result[1] = x[1] + 3.0;
  result[2] = 4.0 - x[0];
  result[3] = 3.0 - x[1];
  return result;
};

void p05_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(4);
  J(0,0) = 1.0; J(0,1) = 0.0;
  J(1,0) = 0.0; J(1,1) = 1.0;
  J(2,0) = -1.0; J(2,1) = 0.0;
  J(3,0) = 0.0;  J(3,1) = -1.0;
};

vect_n<double> p05_sol = vect_n<double>(vect<double,2>(0.5 - 4.0 * atan(1.0) / 3.0,
						       0.5 - 4.0 * atan(1.0) / 3.0 - 1.0));
vect_n<double> p05_start = vect_n<double>(vect<double,2>(0.0,0.0));
vect_n<double> p05_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p05_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));




double p06_f(const vect_n<double>& x) {
  ++evalCount;
  return (1.0 - x[0]) * (1.0 - x[0]);
};

vect_n<double> p06_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = -2.0 * (1.0 - x[0]); 
  result[1] = 0.0; 
  return result;
};

void p06_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 2.0;
  H(0,1) = 0.0;
  H(1,1) = 0.0;
};

vect_n<double> p06_g(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = 10.0 * (x[1] - x[0] * x[0]);
  return result;
};

void p06_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -20.0 * x[0];
  J(0,1) = 10.0;
};

vect_n<double> p06_h(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p06_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p06_sol = vect_n<double>(vect<double,2>(1.0,
						       1.0));
vect_n<double> p06_start = vect_n<double>(vect<double,2>(-1.2,1.0));
vect_n<double> p06_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p06_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));






double p07_f(const vect_n<double>& x) {
  ++evalCount;
  return log(1.0 + x[0] * x[0]) - x[1];
};

vect_n<double> p07_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = 2.0 * x[0] / (1.0 + x[0] * x[0]);
  result[1] = -1.0; 
  return result;
};

void p07_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  double x2 = x[0] * x[0];
  H(0,0) = 2.0 / (1.0 + x2) - 4.0 * x2 / ((1.0 + x2) * (1.0 + x2));
  H(0,1) = 0.0;
  H(1,1) = 0.0;
};

vect_n<double> p07_g(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = (1.0 + x[0] * x[0]) * (1.0 + x[0] * x[0]) + x[1] * x[1] - 4.0;
  return result;
};

void p07_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = 4.0 * x[0] * (1.0 + x[0] * x[0]);
  J(0,1) = 2.0 * x[1];
};

vect_n<double> p07_h(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p07_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p07_sol = vect_n<double>(vect<double,2>(0.0,
						       sqrt(3.0)));
vect_n<double> p07_start = vect_n<double>(vect<double,2>(2.0,2.0));
vect_n<double> p07_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p07_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));







double p10_f(const vect_n<double>& x) {
  ++evalCount;
  return x[0] - x[1];
};

vect_n<double> p10_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = 1.0;
  result[1] = -1.0; 
  return result;
};

void p10_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 0.0;
  H(0,1) = 0.0;
  H(1,1) = 0.0;
};

vect_n<double> p10_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p10_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p10_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = -3.0 * x[0] * x[0] + 2.0 * x[0] * x[1] - x[1] * x[1] + 1.0;
  return result;
};

void p10_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -6.0 * x[0] + 2.0 * x[1];
  J(0,1) = 2.0 * (x[0] - x[1]);
};

vect_n<double> p10_sol = vect_n<double>(vect<double,2>(0.0,
						       1.0));
vect_n<double> p10_start = vect_n<double>(vect<double,2>(-10.0,10.0));
vect_n<double> p10_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p10_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));





double p11_f(const vect_n<double>& x) {
  ++evalCount;
  return (x[0] - 5.0) * (x[0] - 5.0) + x[1] * x[1] - 25.0;
};

vect_n<double> p11_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = 2.0 * (x[0] - 5.0);
  result[1] = 2.0 * x[1];
  return result;
};

void p11_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>&, double, const vect_n<double>&) {
  H.set_col_count(2);
  H(0,0) = 2.0;
  H(0,1) = 0.0;
  H(1,1) = 2.0;
};

vect_n<double> p11_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p11_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p11_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = x[1] - x[0] * x[0];
  return result;
};

void p11_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -2.0 * x[0];
  J(0,1) = 1.0;
};

vect_n<double> p11_get_sol() {
  vect_n<double> result(2);
  double AEX = 7.5 * sqrt(6.0);
  double AW = pow(sqrt(AEX * AEX + 1.0) + AEX,1.0 / 3.0);  
  double QAW = AW * AW;  
  result[0] = (AW - 1.0 / AW) / sqrt(6.0);   
  result[1] = (QAW - 2.0 + 1.0 / QAW) / 6.0;
  return result;
};


vect_n<double> p11_sol = p11_get_sol();
vect_n<double> p11_start = vect_n<double>(vect<double,2>(4.9,0.1));
vect_n<double> p11_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p11_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));






double p12_f(const vect_n<double>& x) {
  ++evalCount;
  return 0.5 * x[0] * x[0] + x[1] * x[1] - x[0] * x[1] - 7.0 * x[0] - 7.0 * x[1];
};

vect_n<double> p12_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = x[0] - x[1] - 7.0;
  result[1] = 2.0 * x[1] - x[0] - 7.0; 
  return result;
};

void p12_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 1.0;
  H(0,1) = -1.0;
  H(1,1) = 2.0;
};

vect_n<double> p12_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p12_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p12_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = 25.0 - 4.0 * x[0] * x[0] - x[1] * x[1];
  return result;
};

void p12_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -8.0 * x[0];
  J(0,1) = -2.0 * x[1];
};

vect_n<double> p12_sol = vect_n<double>(vect<double,2>(2.0,
						       3.0));
vect_n<double> p12_start = vect_n<double>(vect<double,2>(0.0,0.0));
vect_n<double> p12_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p12_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));




double p14_f(const vect_n<double>& x) {
  ++evalCount;
  return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 1.0) * (x[1] - 1.0);
};

vect_n<double> p14_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = 2.0 * (x[0] - 2.0);
  result[1] = 2.0 * (x[1] - 1.0); 
  return result;
};

void p14_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 2.0;
  H(0,1) = 0.0;
  H(1,1) = 2.0;
};

vect_n<double> p14_g(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = x[0] - 2.0 * x[1] + 1.0;
  return result;
};

void p14_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = 1.0;
  J(0,1) = -2.0;
};

vect_n<double> p14_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = 1.0 - x[0] * x[0] * 0.25 - x[1] * x[1];
  return result;
};

void p14_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -0.5 * x[0];
  J(0,1) = -2.0 * x[1];
};

vect_n<double> p14_sol = vect_n<double>(vect<double,2>((sqrt(7.0) - 1.0) * 0.5,
						       (sqrt(7.0) + 1.0) * 0.25));
vect_n<double> p14_start = vect_n<double>(vect<double,2>(2.0,2.0));
vect_n<double> p14_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p14_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));





double p15_f(const vect_n<double>& x) {
  ++evalCount;
  return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 1.0) * (x[1] - 1.0);
};

vect_n<double> p15_grad(const vect_n<double>& x) {
  ++gradCount;
  vect_n<double> result(2);
  result[0] = 2.0 * (x[0] - 2.0);
  result[1] = 2.0 * (x[1] - 1.0); 
  return result;
};

void p15_H(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) {
  H.set_col_count(2);
  H(0,0) = 2.0;
  H(0,1) = 0.0;
  H(1,1) = 2.0;
};

vect_n<double> p15_g(const vect_n<double>&) {
  return vect_n<double>(0);
};

void p15_g_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>&, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(0);
};

vect_n<double> p15_h(const vect_n<double>& x) {
  vect_n<double> result(1);
  result[0] = 1.0 - x[0] * x[0] * 0.25 - x[1] * x[1];
  return result;
};

void p15_h_jac(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  J.set_col_count(2);
  J.set_row_count(1);
  J(0,0) = -0.5 * x[0];
  J(0,1) = -2.0 * x[1];
};

vect_n<double> p15_sol = vect_n<double>(vect<double,2>(1.66497,
						       0.55405));
vect_n<double> p15_start = vect_n<double>(vect<double,2>(2.0,2.0));
vect_n<double> p15_lower = vect_n<double>(vect<double,2>(-std::numeric_limits<double>::infinity(),
							 -std::numeric_limits<double>::infinity()));
vect_n<double> p15_upper = vect_n<double>(vect<double,2>(std::numeric_limits<double>::infinity(),
							 std::numeric_limits<double>::infinity()));





int main() {
  
  std::vector< FunctionPtr > funcs_f;
  std::vector< GradFunctionPtr > funcs_grad;
  std::vector< HessFunctionPtr > funcs_hess;
  std::vector< GFunctionPtr > funcs_g;
  std::vector< GJacFunctionPtr > funcs_g_jac;
  std::vector< HFunctionPtr > funcs_h;
  std::vector< HJacFunctionPtr > funcs_h_jac;
  std::vector< vect_n<double> > funcs_sol;
  std::vector< vect_n<double> > funcs_start;
  std::vector< vect_n<double> > funcs_lower;
  std::vector< vect_n<double> > funcs_upper;
  
  funcs_f.push_back(p01_f);
  funcs_grad.push_back(p01_grad);
  funcs_hess.push_back(p01_H);
  funcs_g.push_back(p01_g);
  funcs_g_jac.push_back(p01_g_jac);
  funcs_h.push_back(p01_h);
  funcs_h_jac.push_back(p01_h_jac);
  funcs_sol.push_back(p01_sol);
  funcs_start.push_back(p01_start);
  funcs_lower.push_back(p01_lower);
  funcs_upper.push_back(p01_upper);
  
  funcs_f.push_back(p02_f);
  funcs_grad.push_back(p02_grad);
  funcs_hess.push_back(p02_H);
  funcs_g.push_back(p02_g);
  funcs_g_jac.push_back(p02_g_jac);
  funcs_h.push_back(p02_h);
  funcs_h_jac.push_back(p02_h_jac);
  funcs_sol.push_back(p02_sol);
  funcs_start.push_back(p02_start);
  funcs_lower.push_back(p02_lower);
  funcs_upper.push_back(p02_upper);
  
  funcs_f.push_back(p03_f);
  funcs_grad.push_back(p03_grad);
  funcs_hess.push_back(p03_H);
  funcs_g.push_back(p03_g);
  funcs_g_jac.push_back(p03_g_jac);
  funcs_h.push_back(p03_h);
  funcs_h_jac.push_back(p03_h_jac);
  funcs_sol.push_back(p03_sol);
  funcs_start.push_back(p03_start);
  funcs_lower.push_back(p03_lower);
  funcs_upper.push_back(p03_upper);
  
  funcs_f.push_back(p04_f);
  funcs_grad.push_back(p04_grad);
  funcs_hess.push_back(p04_H);
  funcs_g.push_back(p04_g);
  funcs_g_jac.push_back(p04_g_jac);
  funcs_h.push_back(p04_h);
  funcs_h_jac.push_back(p04_h_jac);
  funcs_sol.push_back(p04_sol);
  funcs_start.push_back(p04_start);
  funcs_lower.push_back(p04_lower);
  funcs_upper.push_back(p04_upper);
  
  funcs_f.push_back(p05_f);
  funcs_grad.push_back(p05_grad);
  funcs_hess.push_back(p05_H);
  funcs_g.push_back(p05_g);
  funcs_g_jac.push_back(p05_g_jac);
  funcs_h.push_back(p05_h);
  funcs_h_jac.push_back(p05_h_jac);
  funcs_sol.push_back(p05_sol);
  funcs_start.push_back(p05_start);
  funcs_lower.push_back(p05_lower);
  funcs_upper.push_back(p05_upper);
  
  funcs_f.push_back(p06_f);
  funcs_grad.push_back(p06_grad);
  funcs_hess.push_back(p06_H);
  funcs_g.push_back(p06_g);
  funcs_g_jac.push_back(p06_g_jac);
  funcs_h.push_back(p06_h);
  funcs_h_jac.push_back(p06_h_jac);
  funcs_sol.push_back(p06_sol);
  funcs_start.push_back(p06_start);
  funcs_lower.push_back(p06_lower);
  funcs_upper.push_back(p06_upper);
  
  funcs_f.push_back(p07_f);
  funcs_grad.push_back(p07_grad);
  funcs_hess.push_back(p07_H);
  funcs_g.push_back(p07_g);
  funcs_g_jac.push_back(p07_g_jac);
  funcs_h.push_back(p07_h);
  funcs_h_jac.push_back(p07_h_jac);
  funcs_sol.push_back(p07_sol);
  funcs_start.push_back(p07_start);
  funcs_lower.push_back(p07_lower);
  funcs_upper.push_back(p07_upper);
  
  funcs_f.push_back(p10_f);
  funcs_grad.push_back(p10_grad);
  funcs_hess.push_back(p10_H);
  funcs_g.push_back(p10_g);
  funcs_g_jac.push_back(p10_g_jac);
  funcs_h.push_back(p10_h);
  funcs_h_jac.push_back(p10_h_jac);
  funcs_sol.push_back(p10_sol);
  funcs_start.push_back(p10_start);
  funcs_lower.push_back(p10_lower);
  funcs_upper.push_back(p10_upper);
  
  funcs_f.push_back(p11_f);
  funcs_grad.push_back(p11_grad);
  funcs_hess.push_back(p11_H);
  funcs_g.push_back(p11_g);
  funcs_g_jac.push_back(p11_g_jac);
  funcs_h.push_back(p11_h);
  funcs_h_jac.push_back(p11_h_jac);
  funcs_sol.push_back(p11_sol);
  funcs_start.push_back(p11_start);
  funcs_lower.push_back(p11_lower);
  funcs_upper.push_back(p11_upper);
  
  funcs_f.push_back(p12_f);
  funcs_grad.push_back(p12_grad);
  funcs_hess.push_back(p12_H);
  funcs_g.push_back(p12_g);
  funcs_g_jac.push_back(p12_g_jac);
  funcs_h.push_back(p12_h);
  funcs_h_jac.push_back(p12_h_jac);
  funcs_sol.push_back(p12_sol);
  funcs_start.push_back(p12_start);
  funcs_lower.push_back(p12_lower);
  funcs_upper.push_back(p12_upper);
  
  funcs_f.push_back(p14_f);
  funcs_grad.push_back(p14_grad);
  funcs_hess.push_back(p14_H);
  funcs_g.push_back(p14_g);
  funcs_g_jac.push_back(p14_g_jac);
  funcs_h.push_back(p14_h);
  funcs_h_jac.push_back(p14_h_jac);
  funcs_sol.push_back(p14_sol);
  funcs_start.push_back(p14_start);
  funcs_lower.push_back(p14_lower);
  funcs_upper.push_back(p14_upper);
  
  funcs_f.push_back(p15_f);
  funcs_grad.push_back(p15_grad);
  funcs_hess.push_back(p15_H);
  funcs_g.push_back(p15_g);
  funcs_g_jac.push_back(p15_g_jac);
  funcs_h.push_back(p15_h);
  funcs_h_jac.push_back(p15_h_jac);
  funcs_sol.push_back(p15_sol);
  funcs_start.push_back(p15_start);
  funcs_lower.push_back(p15_lower);
  funcs_upper.push_back(p15_upper);
  
  
  vect_n<double> x;
  for(std::size_t i = 0; i < funcs_f.size(); ++i) {
    std::cout << "****************************************************************************************" << std::endl;
    x = funcs_start[i];
    std::size_t M = funcs_g[i](x).size();
    std::size_t K = funcs_h[i](x).size();
    std::cout << "Searching for minimum of function #" << i << " which has N = " << x.size() << " M = " << M << " K = " << K << std::endl;
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_newton_method_tr(funcs_f[i],funcs_grad[i],funcs_hess[i],2.0,300,1e-6,1e-6,1e-3)
        (x);
      std::cout << "  Newton method (unconstrained) gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Newton method (unconstrained) failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_constraint_newton_method_tr(funcs_f[i],funcs_grad[i],funcs_hess[i],2.0,300,1e-6,1e-6,1e-3)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	.set_ineq_constraints(funcs_h[i],funcs_h_jac[i])
	(x);
      std::cout << "  Augmented Lagrangian method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Augmented Lagrangian method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_constraint_newton_method_tr(funcs_f[i],funcs_grad[i],funcs_hess[i],2.0,300,1e-6,1e-6,1e-3)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	.set_ineq_constraints(funcs_h[i],funcs_h_jac[i])
	.regularize(1e-8)
	(x);
      std::cout << "  Regularized Augmented Lagrangian method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Regularized Augmented Lagrangian method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_nlip_newton_tr(funcs_f[i],funcs_grad[i],funcs_hess[i],2.0,0.1,300,1e-6,1e-3,0.99,0.8)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	.set_ineq_constraints(funcs_h[i],funcs_h_jac[i])
	(x);
      std::cout << "  NL Interior-point method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  NL Interior-point method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_nlip_quasi_newton_tr(funcs_f[i],funcs_grad[i],2.0,0.1,300,1e-6,1e-3,0.99,0.8)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	.set_ineq_constraints(funcs_h[i],funcs_h_jac[i])
	(x);
      std::cout << "  NL Interior-point Quasi-Newton method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  NL Interior-point Quasi-Newton method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    
    if(K == 0) {
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_bosqp_newton_tr(funcs_f[i],funcs_grad[i],funcs_hess[i],2.0,300,1e-6,1e-3,0.8)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	(x);
      std::cout << "  Byrd-Omojokun SQP method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Byrd-Omojokun SQP method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    x = funcs_start[i];
    evalCount = 0; gradCount = 0;
    try {
      optim::make_bosqp_quasi_newton_tr(funcs_f[i],funcs_grad[i],2.0,300,1e-6,1e-3,0.8)
        .set_limiter(boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,funcs_lower[i],funcs_upper[i]))
	.set_eq_constraints(funcs_g[i],funcs_g_jac[i])
	(x);
      std::cout << "  Byrd-Omojokun SQP Quasi-Newton method gives:\n"
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Byrd-Omojokun SQP Quasi-Newton method failed with error: " << e.what() << std::endl
                << "    x = " << x << " with error = " << norm_2(x - funcs_sol[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << "\n"
                << "    f(x) = " << funcs_f[i](x) << " with f(x_opt) = " << funcs_f[i](funcs_sol[i]) << " |f(x) - f(x_opt)| = " << fabs(funcs_f[i](x) - funcs_f[i](funcs_sol[i])) << std::endl;
    };
    
    };
    
    std::cout << "****************************************************************************************" << std::endl;
  };
  
  
  return 0;
};








