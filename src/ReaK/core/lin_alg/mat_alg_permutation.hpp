/**
 * \file mat_alg_permutation.hpp
 * 
 * This library declares matrix specializations for representing and manipulating permutation matrices.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_MAT_ALG_PERMUTATION_HPP
#define REAK_MAT_ALG_PERMUTATION_HPP

#include "mat_alg_general.hpp"
#include "mat_alg_identity.hpp"

#include <boost/mpl/less.hpp>
#include <boost/mpl/and.hpp>


namespace ReaK {


  
  


  
/**
 * This class implements a place-holder or interface-implementation to represent
 * a identity matrix (all entries zero except diagonal is all unity). This is useful to build for example a
 * block-matrix with some identity-matrix blocks, and, of course, requires minimal storage space.
 * 
 * Models: ReadableMatrixConcept and ResizableMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T,mat_structure::permutation, Alignment, Allocator> : public serializable {
  public:    
    
    typedef mat<T,mat_structure::permutation, Alignment, Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef unsigned int size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::vector<size_type> container_type;
    
    typedef void reference;
    typedef T const_reference;
    typedef void pointer;
    typedef void const_pointer;
  
    typedef void row_iterator;
    typedef void const_row_iterator;
    typedef void col_iterator;
    typedef void const_col_iterator;
    
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::permutation);
    
  private:
    container_type idx;
    size_type rowCount; ///< Row Count.
  public:
    /**
     * Default constructor. Sets dimensions to zero.
     */
    mat() : idx(), rowCount(0) { };
    /**
     * Constructs an identity matrix to the given dimensions.
     */
    mat(size_type aRowCount) : idx(aRowCount), rowCount(aRowCount) { 
      for(size_type i = 0; i < rowCount; ++i)
        idx[i] = i;
    };
    
    /**
     * Standard swap function (works with ADL).
     */
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
      lhs.idx.swap(rhs.idx);
      swap(lhs.rowCount,rhs.rowCount);
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    const_reference operator()(size_type i,size_type j) const { 
      if(idx[i] == j) 
        return value_type(1); 
      else 
        return value_type(0); 
    };
    
    /**
     * Sub-matrix operator, accessor for read only.
     * \test PASSED
     */
    mat_const_sub_block<self> operator()(const std::pair<size_type,size_type>& r, const std::pair<size_type,size_type>& c) const {
      return sub(*this)(r,c);
    };
    
    /**
     * Sub-matrix operator, accessor for read only.
     * \test PASSED
     */
    mat_const_col_slice<self> operator()(size_type r, const std::pair<size_type,size_type>& c) const {
      return slice(*this)(r,c);
    };
    
    /**
     * Sub-matrix operator, accessor for read only.
     * \test PASSED
     */
    mat_const_row_slice<self> operator()(const std::pair<size_type,size_type>& r, size_type c) const {
      return slice(*this)(r,c);
    };
    
    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const { return rowCount; };
    
    /**
     * Sets the row-count (number of rows) of the matrix.
     * \param aRowCount new number of rows for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_row_count(size_type aRowCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      rowCount = aRowCount;
      idx.resize(rowCount);
    };
    
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const { return rowCount; };
    
    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_col_count(size_type aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      rowCount = aColCount;
      idx.resize(rowCount);
    };
    
    /**
     * Get the allocator object of the underlying container.
     * \return The allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return idx.get_allocator(); };
    
    /**
     * Get the index of the row that should be in-place of the given row index.
     * Note that for column swap this mapping can simply be applied in reverse, that is, 
     * the returned index is the index of the original column that should appear at the 
     * given destination column index.
     * \param i The index of the original row (or destination column).
     * \return the index of the destination row (or original column).
     */
    size_type operator[](size_type i) const { return idx[i]; };
    
    /**
     * Append a row swap to this permutation matrix.
     * Essentially equivalent to permute(i,j) * (this).
     * \param i The first row involved in the swap.
     * \param j The second row involved in the swap.
     */
    void add_row_swap(size_type i, size_type j) {
      if(i == j)
        return;
      using std::swap;
      swap(idx[i],idx[j]);
    };
    
    /**
     * Append a column swap to this permutation matrix.
     * Essentially equivalent to (this) * permute(i,j).
     * Note also that it is generally more efficient to add row swaps (for example, if 
     * you want to accumulated many column swaps, it is more efficient to accumulated them 
     * as row swaps and then invert (or transpose) the final permutation matrix).
     * \param i The first column involved in the swap.
     * \param j The second column involved in the swap.
     */
    void add_column_swap(size_type i, size_type j) {
      if(i == j)
        return;
      using std::swap;
      size_type p_i = 0;
      size_type p_j = 0;
      for(size_type k = 0; k < rowCount; ++k) {
        if(idx[k] == i)
          p_i = k;
        if(idx[k] == j)
          p_j = k;
      };
      swap(idx[p_i],idx[p_j]);
    };
    
    /**
     * Negate the matrix.
     * \return This matrix, by constant reference.
     * \test PASSED
     */
    mat<value_type,mat_structure::scalar> operator -() const {
      return mat<value_type,mat_structure::scalar>(rowCount,value_type(-1));
    };
    
    friend self operator*(const self& M1, const self& M2) {
      if( M1.get_row_count() != M2.get_row_count() )
        throw std::range_error("Matrix dimensions mismatch!");
      self result(M1.get_row_count());
      for(size_type i = 0; i < M1.get_row_count(); ++i)
        result.idx[i] = M2.idx[M1.idx[i]];
      return result;
    };
    
    /**
     * Transpose the matrix.
     * \param rhs the matrix to be transposed.
     * \return The transpose matrix, by value.
     * \test PASSED
     */
    friend self transpose(const self& rhs) {
      self result(rhs.rowCount);
      for(size_type i = 0; i < rhs.rowCount; ++i)
        result.idx[rhs.idx[i]] = i;
      return result;
    };
    
    /**
     * Transpose and move the matrix.
     * \param rhs the matrix to be transposed and moved (emptied).
     * \return The transpose matrix, by value.
     * \test PASSED
     */
    friend self transpose_move(const self& rhs) {
      return transpose(rhs);
    };
    
    /**
     * Invert the matrix.
     * \param rhs the matrix to be inverted.
     * \return The inverse matrix, by value.
     * \test PASSED
     */
    friend self invert(const self& rhs) {
      return transpose(rhs);
    };
    
    /**
     * Returns the trace of a matrix.
     * \param M A matrix.
     * \return the trace of matrix M.
     * \test PASSED
     */
    friend value_type trace(const self& M) {
      value_type result(0);
      for(size_type i = 0; i < M.rowCount; ++i)
        if(M.idx[i] == i)
          result += value_type(1);
      return result;
    };
        
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(idx)
        & std::pair<std::string, size_type>("rowCount",rowCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(idx)
        & std::pair<std::string, size_type&>("rowCount",rowCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serializable)

};


template <typename T>
struct mat_permutation {
  typedef mat<T,mat_structure::permutation> type;
};



template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_square_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_square_matrix< mat<T,mat_structure::permutation,Alignment,Allocator> > type;
};


/**
 * Column-vector multiplication, always results in a null vector.
 * \param M some matrix.
 * \param V some vector.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if< is_readable_vector<Vector>, 
vect_copy<Vector> >::type::type operator *(const mat<T,mat_structure::permutation,Alignment,Allocator>& M, 
                                           const Vector& V) {
  if(V.size() != M.get_col_count())
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename vect_copy<Vector>::type result_type;
  typedef typename vect_traits<Vector>::size_type size_type;
  
  result_type result(V);
  for(size_type i = 0; i < result.size(); ++i)
    result[i] = V[M[i]];
  return result;
};
    
/**
 * Row-vector multiplication with null matrix, always results in a null vector.
 * \param V some row-vector.
 * \param M a null-matrix.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if< is_readable_vector<Vector>, 
vect_copy<Vector> >::type::type operator *(const Vector& V,
                                           const mat<T,mat_structure::permutation,Alignment,Allocator>& M) {
  if(V.size() != M.get_row_count())
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename vect_copy<Vector>::type result_type;
  typedef typename vect_traits<Vector>::size_type size_type;
  
  result_type result(V);
  for(size_type i = 0; i < result.size(); ++i)
    result[M[i]] = V[i];
  return result;
};

/**
 * Column-vector multiplication, always results in a null vector.
 * \param M some matrix.
 * \param S some scalar.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if< 
  boost::mpl::and_<
    boost::mpl::not_< is_readable_vector<T> >,
    boost::mpl::not_< is_readable_matrix<T> >
  >, 
mat<T,mat_structure::square,Alignment,Allocator> >::type 
 operator *(const mat<T,mat_structure::permutation,Alignment,Allocator>& M, const T& S) {
  typedef typename mat_traits< mat<T,mat_structure::square,Alignment,Allocator> >::size_type size_type;
  mat<T,mat_structure::square,Alignment,Allocator> result(M.get_row_count());
  for(size_type i = 0; i < M.get_row_count(); ++i)
    result(i,M[i]) = S;
  return result;
};
    
/**
 * Row-vector multiplication with null matrix, always results in a null vector.
 * \param V some row-vector.
 * \param M a null-matrix.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if< 
  boost::mpl::and_<
    boost::mpl::not_< is_readable_vector<T> >,
    boost::mpl::not_< is_readable_matrix<T> >
  >, 
mat<T,mat_structure::square,Alignment,Allocator> >::type 
 operator *(const T& S,const mat<T,mat_structure::permutation,Alignment,Allocator>& M) {
  typedef typename mat_traits< mat<T,mat_structure::square,Alignment,Allocator> >::size_type size_type;
  mat<T,mat_structure::square,Alignment,Allocator> result(M.get_row_count());
  for(size_type i = 0; i < M.get_row_count(); ++i)
    result(i,M[i]) = S;
  return result;
};

/**
 * Matrix multiplication with permutation matrix (column permutation).
 * \param M1 some matrix.
 * \param M2 a permutation-matrix.
 * \return A permuted matrix.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_matrix<Matrix>,
    boost::mpl::less<
      mat_product_priority< Matrix >,
      mat_product_priority< mat<T,mat_structure::permutation,Alignment,Allocator> >
    >
  >,
mat<T,mat_structure::rectangular,Alignment,Allocator> >::type
 operator *(const Matrix& M1, const mat<T,mat_structure::permutation,Alignment,Allocator>& M2) {
  if(M1.get_col_count() != M2.get_row_count())
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename mat_traits< mat<T,mat_structure::rectangular,Alignment,Allocator> >::size_type SizeType;
  mat<T,mat_structure::rectangular,Alignment,Allocator> result(M1.get_row_count(),M1.get_col_count());
  for(SizeType i = 0; i < result.get_col_count(); ++i)
    for(SizeType j = 0; j < result.get_row_count(); ++j)
      result(j,M2[i]) = M1(j,i);
  return result;
};

/**
 * Matrix multiplication with permutation matrix (row permutation).
 * \param M1 a permutation-matrix.
 * \param M2 some matrix.
 * \return A permuted matrix.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_matrix<Matrix>,
    boost::mpl::less<
      mat_product_priority< Matrix >,
      mat_product_priority< mat<T,mat_structure::permutation,Alignment,Allocator> >
    >
  >,
mat<T,mat_structure::rectangular,Alignment,Allocator> >::type
 operator *(const mat<T,mat_structure::permutation,Alignment,Allocator>& M1, const Matrix& M2) {
  if(M1.get_col_count() != M2.get_row_count())
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename mat_traits< mat<T,mat_structure::rectangular,Alignment,Allocator> >::size_type SizeType;
  mat<T,mat_structure::rectangular,Alignment,Allocator> result(M2.get_row_count(),M2.get_col_count());
  for(SizeType j = 0; j < result.get_row_count(); ++j)
    for(SizeType i = 0; i < result.get_col_count(); ++i)
      result(j,i) = M2(M1[j],i);
  return result;
};




#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template class mat<double, mat_structure::permutation>;
extern template class mat<float, mat_structure::permutation>;


extern template vect<double,2> operator *(const mat<double,mat_structure::permutation>& M, const vect<double,2>& V);
extern template vect<double,3> operator *(const mat<double,mat_structure::permutation>& M, const vect<double,3>& V);
extern template vect<double,4> operator *(const mat<double,mat_structure::permutation>& M, const vect<double,4>& V);
extern template vect<double,6> operator *(const mat<double,mat_structure::permutation>& M, const vect<double,6>& V);
extern template vect_n<double> operator *(const mat<double,mat_structure::permutation>& M, const vect_n<double>& V);

extern template vect<double,2> operator *(const vect<double,2>& V, const mat<double,mat_structure::permutation>& M);
extern template vect<double,3> operator *(const vect<double,3>& V, const mat<double,mat_structure::permutation>& M);
extern template vect<double,4> operator *(const vect<double,4>& V, const mat<double,mat_structure::permutation>& M);
extern template vect<double,6> operator *(const vect<double,6>& V, const mat<double,mat_structure::permutation>& M);
extern template vect_n<double> operator *(const vect_n<double>& V, const mat<double,mat_structure::permutation>& M);

extern template mat<double,mat_structure::square> operator *(const mat<double,mat_structure::permutation>& M, const double& S);
extern template mat<double,mat_structure::square> operator *(const double& S,const mat<double,mat_structure::permutation>& M);

extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::rectangular>& M1, const mat<double,mat_structure::permutation>& M2);
extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::square>& M1, const mat<double,mat_structure::permutation>& M2);

extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::permutation>& M1, const mat<double,mat_structure::rectangular>& M2);
extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::permutation>& M1, const mat<double,mat_structure::square>& M2);


extern template vect<float,2> operator *(const mat<float,mat_structure::permutation>& M, const vect<float,2>& V);
extern template vect<float,3> operator *(const mat<float,mat_structure::permutation>& M, const vect<float,3>& V);
extern template vect<float,4> operator *(const mat<float,mat_structure::permutation>& M, const vect<float,4>& V);
extern template vect<float,6> operator *(const mat<float,mat_structure::permutation>& M, const vect<float,6>& V);
extern template vect_n<float> operator *(const mat<float,mat_structure::permutation>& M, const vect_n<float>& V);

extern template vect<float,2> operator *(const vect<float,2>& V, const mat<float,mat_structure::permutation>& M);
extern template vect<float,3> operator *(const vect<float,3>& V, const mat<float,mat_structure::permutation>& M);
extern template vect<float,4> operator *(const vect<float,4>& V, const mat<float,mat_structure::permutation>& M);
extern template vect<float,6> operator *(const vect<float,6>& V, const mat<float,mat_structure::permutation>& M);
extern template vect_n<float> operator *(const vect_n<float>& V, const mat<float,mat_structure::permutation>& M);

extern template mat<float,mat_structure::square> operator *(const mat<float,mat_structure::permutation>& M, const float& S);
extern template mat<float,mat_structure::square> operator *(const float& S,const mat<float,mat_structure::permutation>& M);

extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::rectangular>& M1, const mat<float,mat_structure::permutation>& M2);
extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::square>& M1, const mat<float,mat_structure::permutation>& M2);

extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::permutation>& M1, const mat<float,mat_structure::rectangular>& M2);
extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::permutation>& M1, const mat<float,mat_structure::square>& M2);


#endif





};


#endif










