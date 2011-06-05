
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

#ifndef MAT_NORMS_HPP
#define MAT_NORMS_HPP

#include "mat_traits.hpp"
#include "mat_concepts.hpp"

namespace ReaK {



template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename Matrix::value_type >::type norm_1(const Matrix& M) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  ValueType max = ValueType();
  for(SizeType j = 0; j < M.get_col_count(); ++j) {
    ValueType sum = ValueType();
    for(SizeType i = 0; i < M.get_row_count(); ++i)
      sum += fabs(M(i,j));
    if(sum > max)
      max = sum;
  };
  return max;
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename Matrix::value_type >::type norm_inf(const Matrix& M) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  ValueType max = ValueType();
  for(SizeType i = 0; i < M.get_row_count(); ++i) {
    ValueType sum = ValueType();
    for(SizeType j = 0; j < M.get_col_count(); ++j)
      sum += fabs(M(i,j));
    if(sum > max)
      max = sum;
  };
  return max;
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename Matrix::value_type >::type elem_norm_2(const Matrix& M) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::sqrt;
  
  ValueType sum = ValueType();
  for(SizeType i = 0; i < M.get_row_count(); ++i)
    for(SizeType j = 0; j < M.get_col_count(); ++j)
      sum += M(i,j) * M(i,j); 
  return sqrt(sum);
};


template <typename Matrix>
typename Matrix::value_type frobenius_norm(const Matrix& M) {
  return elem_norm_2(M);
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename Matrix::value_type >::type elem_norm_max(const Matrix& M) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  ValueType max = ValueType();
  for(SizeType i = 0; i < M.get_row_count(); ++i)
    for(SizeType j = 0; j < M.get_col_count(); ++j)
      if(max < fabs(M(i,j)))
        max = fabs(M(i,j));
  return max;
};




};

#endif











