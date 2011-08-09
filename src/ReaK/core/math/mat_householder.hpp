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

/**
 * This class represents a Householder reflection as a NxN matrix. It can be constructed 
 * by providing the column-vector u such that the resulting Householder reflection H has
 * the effect of zero-ing all but the first element: H * u = (a, 0) (where u is a column-vector 
 * and a is some scalar value). To apply a Householder reflection efficiently, it is preferred 
 * to use the householder_prod function templates, although this Householder reflection matrix models 
 * a readable matrix concept and can thus be involved in other matrix arithmetic (although rarely 
 * useful). In order to apply the Householder reflection to a matrix with larger dimensions, which 
 * is most often the case, one can use the classes provided in mat_views.hpp and 
 * mat_composite_adaptor.hpp to create such sub-matrices and composites of them to extract only 
 * the rows or columns affected by the reflection (the resulting sub-matrix (or matrix-view) 
 * should have the appropriate dimensions). This class also provides friend functions to 
 * perform the transposition of the Householder reflection (and effectively, its inverse).
 * 
 * Models: ReadableMatrixConcept.
 * 
 * \tparam Vector The vector type to store the householder vector.
 */
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
    
    value_type beta; ///< Holds the householder coefficient.
    Vector v; ///< Holds the householder vector.
    
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
    /**
     * Set the Householder reflection by providing a column-vector u such that the 
     * resulting Householder reflection H has the effect of zero-ing the second element: H * u = (a, 0).
     * \tparam Vector2 A readable vector type.
     * \param aE The vector from which to calculate the Householder reflection (vector u).
     * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
     */
    template <typename Vector2>
    typename boost::enable_if_c< is_readable_vector<Vector2>::value, 
    void >::type set(const Vector2& aE, const value_type& NumTol = value_type(1E-8)) {
      beta = value_type(0.0);
      v = aE;
      calculate_hhv(NumTol);
    };
    /**
     * Default constructor.
     */
    explicit householder_matrix(const value_type& aBeta = 0.0, 
				const Vector& aV = Vector()) : 
				beta(aBeta),
				v(aV) { };

    /**
     * Constructs the Householder reflection by providing a column-vector u such that the 
     * resulting Householder reflection H has the effect of zero-ing the second element: H * u = (a, 0).
     * \tparam Vector2 A readable vector type.
     * \param aE The vector from which to calculate the Householder reflection (vector u).
     * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
     */
    template <typename Vector2>
    explicit householder_matrix(const Vector2& aE, const value_type& NumTol = value_type(1E-8), typename boost::enable_if_c< is_readable_vector<Vector2>::value, void* >::type dummy = NULL) :
                                beta(0.0), v(aE) {
      calculate_hhv(NumTol);
    };
    
    /**
     * Standard copy-constructor.
     */
    householder_matrix(const self& rhs) : beta(rhs.beta), v(rhs.v) { };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.beta,rhs.beta);
      swap(lhs.v,rhs.v);
    };
    
    /**
     * Standard assignment operator (copy-and-swap and move-and-swap).
     */
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const { return v.size(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const { return v.size(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(v.size(),v.size()); };
    
    /**
     * Returns the allocator object of the underlying container, which is none at all in this case.
     * \return the allocator object of the underlying container, which is none at all in this case.
     */
    allocator_type get_allocator() const { return v.get_allocator(); };
    
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) { return (i == j ? value_type(1.0) : value_type(0.0)) - beta * v[i] * v[j]; };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend const self& transpose(const self& M) {
      return M;
    };
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend self transpose_move(self& M) {
      self result;
      swap(M,result);
      return result;
    };
    
    /**
     * Returns the trace of the matrix M.
     * \param M The matrix.
     * \return The trace of M.
     */
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





/**
 * This function template allows for efficient post-multiplication of a matrix with a 
 * Householder reflection matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam Vector A vector-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param A The matrix to be multiplied by the Householder reflection, stores, as output, the resulting matrix.
 * \param P The Householder reflection which will post-multiply A.
 */
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

/**
 * This function template allows for efficient pre-multiplication of a matrix with a 
 * Householder reflection matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam Vector A vector-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param A The matrix to be multiplied by the Householder reflection, stores, as output, the resulting matrix.
 * \param P The Householder reflection which will pre-multiply A.
 */
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









