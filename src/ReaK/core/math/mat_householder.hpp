/**
 * \file mat_householder.hpp
 * 
 * This library provides the means to perform Householder reflections on matrices, efficiently. 
 * This library provides a class template that can be constructed from the relevant elements 
 * to be reflected by the Householder reflection, and it creates a representation of a Householder reflection 
 * that can then be used to multiply a matrix (or sub-matrix) efficiently (with minimal calculations).
 * 
 * Householder reflection are useful in a number of matrix numerical methods. It is a fundamental construct.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef MAT_HOUSEHOLDER_HPP
#define MAT_HOUSEHOLDER_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

namespace ReaK {


template <typename Vector>
class householder_matrix {
  public:
    typedef householder_matrix<Vector> self;
    typedef typename vect_traits<Vector>::value_type value_type;
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::difference_type difference_type;
    typedef typename vect_traits<Vector>::allocator_type allocator_type;
    
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef typename vect_traits<Vector>::const_pointer const_pointer;
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::const_reference const_reference;
    typedef typename vect_traits<Vector>::iterator iterator;
    typedef typename vect_traits<Vector>::const_iterator const_iterator;
    
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::symmetric);
    
    value_type beta;
    Vector v;
    
  private:
    void calculate_hhv(const value_type& NumTol) {      
      if(v.size() == 0)
	return;
      
      using std::sqrt;
      /*
      value_type alpha = value_type(0.0);
      for(size_type i=0;i<v.size();++i)
	alpha += v[i] * v[i];
      alpha = sqrt(alpha);
      if(v[0] > 0)
	alpha = -alpha;
      
      value_type rsq = alpha * alpha - alpha * v[0];
      beta = 2.0 / rsq;
      v[0] -= alpha;
      */
      
      value_type sigma = value_type(0.0);
      for(size_type i=1;i<v.size();++i)
	sigma += v[i] * v[i];
      if(sigma < NumTol) {
	v[0] = value_type(1.0);
	beta = value_type(0.0);
	return;
      };
      value_type mu = sqrt(sigma + v[0]*v[0]);
      if(v[0] < NumTol) {
	v[0] -= mu;
      } else {
	v[0] = -sigma / (v[0] + mu);
      };
      beta = 2.0 / (sigma + v[0]*v[0]);
    };
  public:
    template <typename Vector2>
    typename boost::enable_if_c< is_readable_vector<Vector2>::value, 
    void >::type set(const Vector2& aE, const value_type& NumTol = value_type(1E-8)) {
      beta = value_type(0.0);
      v = aE;
      calculate_hhv(NumTol);
    };
    
    explicit householder_matrix(const value_type& aBeta = 0.0, 
				const Vector& aV = Vector()) : 
				beta(aBeta),
				v(aV) { };

    template <typename Vector2>
    explicit householder_matrix(const Vector2& aE, const value_type& NumTol = value_type(1E-8), typename boost::enable_if_c< is_readable_vector<Vector2>::value, void* >::type dummy = NULL) :
                                beta(0.0), v(aE) {
      calculate_hhv(NumTol);
    };
    
    householder_matrix(const self& rhs) : beta(rhs.beta), v(rhs.v) { };
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.beta,rhs.beta);
      swap(lhs.v,rhs.v);
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    size_type get_row_count() const { return v.size(); };
    size_type get_col_count() const { return v.size(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(v.size(),v.size()); };
    
    allocator_type get_allocator() const { return v.get_allocator(); };
    
    value_type operator()(size_type i,size_type j) { return (i == j ? value_type(1.0) : value_type(0.0)) - beta * v[i] * v[j]; };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    friend self transpose_move(self& M) {
      self result;
      swap(M,result);
      return result;
    };
    
    friend value_type trace(const self& M) {
      return value_type(M.v.size()) - M.beta * (M.v * M.v);
    };
    
};



template <typename Vector>
struct is_readable_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct is_writable_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct is_resizable_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct has_allocator_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector<Vector>::value );
  typedef has_allocator_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct mat_product_priority< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<mat_structure::diagonal>::value+1);
};

template <typename Vector>
struct mat_addition_priority< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<mat_structure::square>::value);
};

template <typename Vector>
struct is_square_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = true);
  typedef is_square_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct is_symmetric_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = true);
  typedef is_symmetric_matrix< householder_matrix<Vector> > type;
};

template <typename Vector>
struct is_diagonal_matrix< householder_matrix<Vector> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_diagonal_matrix< householder_matrix<Vector> > type;
};






template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
void >::type householder_prod(Matrix& A, const householder_matrix<Vector>& P) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  if(fabs(P.beta) < std::numeric_limits<typename householder_matrix<Vector>::value_type>::epsilon())
    return;
  for(SizeType i=0;i<A.get_row_count();++i) {
    ValueType temp = ValueType(0);
    for(SizeType j=0;j<P.v.size();++j)
      temp += P.beta * A(i,j) * P.v[j];
    for(SizeType j=0;j<P.v.size();++j)
      A(i,j) -= temp * P.v[j];
  };
  
  return;
};

template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
void >::type householder_prod(const householder_matrix<Vector>& P, Matrix& A) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  if(fabs(P.beta) < std::numeric_limits<typename householder_matrix<Vector>::value_type>::epsilon())
    return;
  for(SizeType i=0;i<A.get_col_count();++i) {
    ValueType temp = ValueType(0);
    for(SizeType j=0;j<P.v.size();++j)
      temp += P.beta * A(j,i) * P.v[j];
    for(SizeType j=0;j<P.v.size();++j)
      A(j,i) -= temp * P.v[j];
  };
  
  return;
};




template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
Matrix >::type operator *(const Matrix& A, const householder_matrix<Vector>& P) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  Matrix result(A);
  householder_prod(result,P);
  return result;
};

template <typename Matrix, typename Vector>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
Matrix >::type operator *(const householder_matrix<Vector>& P, const Matrix& A) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  
  Matrix result(A);
  householder_prod(P,result);
  return result;
};






};

#endif









