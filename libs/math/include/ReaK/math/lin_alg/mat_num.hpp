/**
 * \file mat_num.hpp
 *
 * ReaK Matrix Numerical Methods Library
 *
 * This library declares all functions for matrix numerical methods. Mainly includes inversion, determinants,
 * singular value decomposition, eigen values / vectors, pseudo-inverse, linear system solution, etc.
 *
 * Note: All matrix memory is organized, by default, such that columns are concatenated. This
 *       was found to be a more efficient representation since columns often have
 *       more significances than rows (representing basis vectors for example).
 *
 * \author Mikael Persson, B.Eng.
 * \date march 2010
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

#ifndef REAK_MAT_NUM_HPP
#define REAK_MAT_NUM_HPP

#include "mat_alg.hpp"

#include "mat_are_solver.hpp"
#include "mat_balance.hpp"
#include "mat_cholesky.hpp"
#include "mat_damped_matrix.hpp"
#include "mat_exp_methods.hpp"
#include "mat_gaussian_elim.hpp"
#include "mat_givens_rot.hpp"
#include "mat_hess_decomp.hpp"
#include "mat_householder.hpp"
#include "mat_jacobi_method.hpp"
#include "mat_norms.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_schur_decomp.hpp"
#include "mat_svd_method.hpp"

#endif
