
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

#ifndef MAT_BALANCE_HPP
#define MAT_BALANCE_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

namespace ReaK {
  
/**
 * Performs matrix balancing of a square matrix such that the sum of absolute values of row elements and 
 * column elements match in order of magnitude. This algorithm produces a diagonal matrix that can
 * scale the matrix, as A_balanced = D^-1 A D
 *
 * \param A square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param D holds as output, the diagonal matrix D which balances A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 * 
 * Taken from Golub & Van Loan, "Matrix Computations" (3rd ed).
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix1>::value && 
                             is_writable_matrix<Matrix2>::value, 
void >::type balance(Matrix1& A, Matrix2& D)
{  
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("Upper-Hessenberg decomposition is only possible on a square matrix!");

  using std::fabs;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::value_type SizeType;
  
  SizeType N = A.get_row_count();
  D = mat<typename mat_traits<Matrix2>::value_type, mat_structure::identity>(N);

  bool keep_going = true;

  while (keep_going) {
    size_t i, j;
    double g, f, s;

    keep_going = false;

    for (SizeType i=0; i<N; ++i) {
      ValueType row_mag = ValueType();
      ValueType col_mag = ValueType();
      
      for (SizeType j=0; j<N; ++j) {
	if (j != i) {
	  col_mag += fabs(A(j,i));
          row_mag += fabs(A(i,j));
	};
      };

      if ((col_mag < std::numeric_limits<ValueType>::epsilon()) 
	|| (row_mag < std::numeric_limits<ValueType>::epsilon())) 
	continue;

      ValueType g = row_mag / ValueType(2.0);
      ValueType f = ValueType(1.0);
      ValueType s = col_mag + row_mag;

      while (col_mag < g) {
	f *= ValueType(2.0);
        col_mag *= ValueType(4.0);
      };

      g = row_mag * ValueType(2.0);
      while (col_mag > g) {
	f /= ValueType(2.0);
        col_mag /= ValueType(4.0);
      };

      if ((row_mag + col_mag) < ValueType(0.95) * s * f) {
        keep_going = true;
        g = ValueType(1.0) / f;
        mat_row_slice<Matrix1> row_i(A,i,0,N); row_i *= g;
	mat_col_slice<Matrix1> col_i(A,i,0,N); col_i *= f;
	D(i,i) *= f;
      };
    };
  };
};


};

#endif















