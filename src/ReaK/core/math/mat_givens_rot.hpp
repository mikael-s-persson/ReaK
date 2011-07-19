
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

#ifndef MAT_GIVENS_ROT_HPP
#define MAT_GIVENS_ROT_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

namespace ReaK {


template <typename T>
class givens_rot_matrix {
  public:
    typedef givens_rot_matrix<T> self;
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef void allocator_type;
    
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* iterator;
    typedef const T* const_iterator;
    
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::orthogonal);
    
    value_type c;
    value_type s;
    
  private:
    void calculate_givens(const value_type& NumTol) {      
      using std::fabs;
      using std::sqrt;
      if(fabs(s) < NumTol * fabs(c)) {
	c = value_type(1.0);
	s = value_type(0.0);
	return;
      };
      
      if(fabs(c) < fabs(s)) {
	value_type tau = - c / s;
	s = value_type(1.0) / sqrt(value_type(1.0) + tau * tau);
	c = s * tau;
      } else {
	value_type tau = - s / c;
	c = value_type(1.0) / sqrt(value_type(1.0) + tau * tau);
	s = c * tau;
      };
    };
  public:
    void set(const value_type& aA, const value_type& aB, const value_type& NumTol = value_type(1E-8)) {
      c = aA;
      s = aB;
      calculate_givens(NumTol);
    };
    
    explicit givens_rot_matrix(const value_type& aA, const value_type& aB, const value_type& NumTol = value_type(1E-8)) :
                               c(aA), s(aB) {
      calculate_givens(NumTol);
    };
    
    givens_rot_matrix(const self& rhs) : c(rhs.c), s(rhs.s) { };
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.c,rhs.c);
      swap(lhs.s,rhs.s);
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    size_type get_row_count() const { return 2; };
    size_type get_col_count() const { return 2; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(2,2); };
    
    allocator_type get_allocator() const { return; };
    
    value_type operator()(size_type i,size_type j) { return (i == j ? c : (i - j) * s); };
    
    friend self transpose(self M) {
      M.s = -M.s;
      return M;
    };
    friend self transpose_move(self& M) {
      self result;
      swap(M,result);
      result.s = - result.s;
      return result;
    };
    
    friend value_type trace(const self& M) {
      return M.c + M.c;
    };
    
};



template <typename T>
struct is_readable_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct is_writable_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct is_resizable_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct has_allocator_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct mat_product_priority< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<mat_structure::diagonal>::value+1);
};

template <typename T>
struct mat_addition_priority< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<mat_structure::square>::value);
};

template <typename T>
struct is_square_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true);
  typedef is_square_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct is_symmetric_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_symmetric_matrix< givens_rot_matrix<T> > type;
};

template <typename T>
struct is_diagonal_matrix< givens_rot_matrix<T> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_diagonal_matrix< givens_rot_matrix<T> > type;
};





template <typename Matrix, typename T>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
void >::type givens_rot_prod(Matrix& A, const givens_rot_matrix<T>& G,
                             typename mat_traits<Matrix>::size_type j = 0,
                             typename mat_traits<Matrix>::size_type k = 1) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  if(fabs(G.s) < std::numeric_limits<typename givens_rot_matrix<T>::value_type>::epsilon())
    return;
  for(SizeType i=0;i<A.get_row_count();++i) {
    ValueType t0 = A(i,j);
    ValueType t1 = A(i,k);
    A(i,j) = G.c * t0 + G.s * t1;
    A(i,k) = G.c * t1 - G.s * t0;
  };
  return;
};

template <typename Matrix, typename T>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
void >::type givens_rot_prod(const givens_rot_matrix<T>& G, Matrix& A,
                             typename mat_traits<Matrix>::size_type j = 0,
                             typename mat_traits<Matrix>::size_type k = 1) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  if(fabs(G.s) < std::numeric_limits<typename givens_rot_matrix<T>::value_type>::epsilon())
    return;
  for(SizeType i=0;i<A.get_col_count();++i) {
    ValueType t0 = A(j,i);
    ValueType t1 = A(k,i);
    A(j,i) = G.c * t0 - G.s * t1;
    A(k,i) = G.s * t0 + G.c * t1;
  };
  
  return;
};




template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
Matrix >::type operator *(const Matrix& A, const givens_rot_matrix<Vector>& G) {
  Matrix result(A);
  givens_rot_matrix(result,G);
  return result;
};

template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
Matrix >::type operator *(const givens_rot_matrix<Vector>& G, const Matrix& A) {
  Matrix result(A);
  givens_rot_matrix(G,result);
  return result;
};




};

#endif









