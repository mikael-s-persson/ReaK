/**
 * \file mat_alg_orthogonal.hpp
 *
 * This library is not done yet. Not sure if there is any point in making it.
 *
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

#ifndef REAK_MAT_ALG_ORTHOGONAL_HPP
#define REAK_MAT_ALG_ORTHOGONAL_HPP

#include "ReaK/math/lin_alg/mat_alg_general.hpp"

namespace ReaK {

template <>
struct mat_indexer<mat_structure::orthogonal, mat_alignment::column_major> {
  int rowCount;
  explicit mat_indexer<mat_structure::orthogonal, mat_alignment::column_major>(
      int aRowCount)
      : rowCount(aRowCount) {}
  int operator()(int i, int j) const { return j * rowCount + i; }
};

template <>
struct mat_indexer<mat_structure::orthogonal, mat_alignment::row_major> {
  int rowCount;
  explicit mat_indexer<mat_structure::orthogonal, mat_alignment::row_major>(
      int aRowCount)
      : rowCount(aRowCount) {}
  int operator()(int i, int j) const { return i * rowCount + j; }
};

}  // namespace ReaK

#endif
