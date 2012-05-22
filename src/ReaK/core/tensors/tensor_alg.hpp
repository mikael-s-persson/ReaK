/**
 * \file tensor_alg.hpp
 * 
 * This library declares (mostly by inclusions) all types of tensors currently available 
 * in the ReaK platform. These are all STL-vector based storage tensors or special tensors
 * such as nil or identity, or tensor adaptors, views, compositions and slices.
 * 
 * \todo Implement a suite of statically-sized matrix classes.
 * \todo Port the code related to the upper-triangular and lower-triangular matrices.
 * \todo Implement expression templates to optimize compound matrix expressions.
 * \todo Implement additional matrix views, for example, transposed view.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_TENSOR_ALG_HPP
#define REAK_TENSOR_ALG_HPP

#include "mat_alg_general.hpp"
#include "mat_operators.hpp"
#include "mat_comparisons.hpp"
#include "mat_alg_rectangular.hpp"
#include "mat_alg_square.hpp"
#include "mat_alg_nil.hpp"
#include "mat_alg_identity.hpp"
#include "mat_alg_scalar.hpp"
#include "mat_alg_symmetric.hpp"
#include "mat_alg_skew_symmetric.hpp"
#include "mat_alg_diagonal.hpp"
#include "mat_alg_orthogonal.hpp"
#include "mat_alg_lower_triangular.hpp"
#include "mat_alg_upper_triangular.hpp"
#include "mat_alg_permutation.hpp"

#include "mat_vector_adaptor.hpp"
#include "mat_views.hpp"
#include "mat_transpose_view.hpp"
#include "mat_slices.hpp"
#include "mat_composite_adaptor.hpp"

/** Main namespace for ReaK */
namespace ReaK {
  




};



#endif









