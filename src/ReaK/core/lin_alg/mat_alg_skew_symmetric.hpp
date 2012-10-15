/**
 * \file mat_alg_skew_symmetric.hpp
 * 
 * This library implements the specialization of the mat<> template for a 
 * skew-symmetric matrix (dynamic dimension). This matrix type fulfills the matrix 
 * concepts of Readable, Writable, Resizable and DynAlloc, but not FullyWritable.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011 (originally february 2010)
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

#ifndef REAK_MAT_ALG_SKEW_SYMMETRIC_HPP
#define REAK_MAT_ALG_SKEW_SYMMETRIC_HPP

#include "mat_alg_general.hpp"


namespace ReaK {




//....not good, the following: (ill-defined)
//NOT USED. DEPRECATED
template <mat_alignment::tag Alignment>
struct mat_indexer<mat_structure::skew_symmetric,Alignment> {
  std::size_t rowCount;
  mat_indexer<mat_structure::skew_symmetric, Alignment>(std::size_t aRowCount) : rowCount(aRowCount) { };
  std::size_t mat_triangular_size(std::size_t Size) {
    return (Size * (Size - 1)) / 2 + Size;
  };
  std::size_t operator()(std::size_t i,std::size_t j) const {
    if(i > j)
      return mat_triangular_size(i) + j;
    else
      return mat_triangular_size(j) + i;      
 /*   if(i > j)
      return -q[mat_triangular_size(i-1) + j];
    else if(i < j)
      return q[mat_triangular_size(j-1) + i];
    else
      return T(0.0);*/
  };
};




/**
 * This class holds a skew-symmetric matrix. This class will hold only the strict upper-triangular part
 * since the lower part is assumed to be equal to the negative of the upper one.
 * 
 * Models: ReadableMatrixConcept, WritableMatrixConcept, ResizableMatrixConcept, and DynAllocMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T,mat_structure::skew_symmetric,Alignment,Allocator> : public serialization::serializable {
  public:    
    
    typedef mat<T,mat_structure::skew_symmetric,Alignment,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
  
    typedef void col_iterator;
    typedef void const_col_iterator;
    typedef void row_iterator;
    typedef void const_row_iterator;
  
    typedef unsigned int size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::skew_symmetric);
    
  private:
    std::vector<T,Allocator> q; ///< Holds the array of scalar entries.
    size_type rowCount; ///< Holds the dimension, both row and column count are equal to size.
    
    static size_type mat_triangular_size(size_type Size) {
      return (Size * (Size - 1)) / 2 + Size;
    };
    
  public:

   
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
    /**
     * Default constructor: sets all to zero.
     * \test PASSED
     */
    mat(const allocator_type& aAlloc = allocator_type()) :
        q(0,value_type(0),aAlloc),
	rowCount(0) { };

    /**
     * Constructor for a sized matrix.
     * \test PASSED
     */
    mat(size_type aRowCount, value_type aFill = 0, const allocator_type& aAlloc = allocator_type()) :
        q(mat_triangular_size(aRowCount-1),aFill,aAlloc),
	rowCount(aRowCount) { };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(const self& M) : q(M.q), rowCount(M.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(self&& M) : q(std::move(M.q)), rowCount(std::move(M.rowCount)) { };
#endif

    /**
     * Explicit constructor from any type of matrix. The "(M - M.transpose) / 2" is applied to guarantee skew-symmetry.
     * \test PASSED
     */
    template <typename Matrix>
    explicit mat(const Matrix&  M, const allocator_type& aAlloc = allocator_type(), 
		                   typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                               !(boost::is_same<Matrix,self>::value) , void* >::type dummy = NULL) :
             q(mat_triangular_size((M.get_row_count() > M.get_col_count() ? M.get_row_count() : M.get_col_count()) - 1),value_type(0),aAlloc),
	     rowCount((M.get_row_count() > M.get_col_count() ? M.get_row_count() : M.get_col_count())) {
      size_type k=0;
      size_type i=1;
      size_type min_size = (M.get_row_count() > M.get_col_count() ? M.get_col_count() : M.get_row_count());
      for(;i<min_size;k += i++) {
	for(size_type j=0;j<i;++j) {
	  q[k+j] = value_type(0.5) * (M(j,i) - M(i,j));
	};
      };
      if(M.get_row_count() > M.get_col_count()) {
        for(;i<rowCount;k += i++) {
   	  for(size_type j=0;j<min_size;++j)
	    q[k+j] = value_type(-0.5) * M(i,j);
        };
      } else {
	for(;i<rowCount;k += i++) {
   	  for(size_type j=0;j<min_size;++j)
	    q[k+j] = value_type(0.5) * M(j,i);
        };
      };
    };
    
    /**
     * Constructor from a symmetric matrix, i.e., takes the skew-symmetric part of a symmetric matrix, which is null.
     * \test PASSED
     */
    template <mat_alignment::tag Align2, typename Alloc2>
    mat(const mat<value_type,mat_structure::symmetric,Align2,Alloc2>& M, 
	const allocator_type& aAlloc = allocator_type()) : 
	q(mat_triangular_size(M.get_row_count()-1),value_type(0),aAlloc), 
	rowCount(M.get_row_count()) { };

    /**
     * Destructor.
     * \test PASSED
     */
    ~mat() { };

    /**
     * Constructs a 2x2 skew-symmetric matrix from one element.
     * \test PASSED
     */
    explicit mat(const_reference a12, 
		 const allocator_type& aAlloc = allocator_type()) : q(1,value_type(0),aAlloc), rowCount(2) {
      q[0] = a12;
    };

    /**
     * Constructs a 3x3 skew-symmetric matrix from 3 elements.
     * \test PASSED
     */
    mat(const_reference a12,const_reference a13,
	                    const_reference a23, 
	const allocator_type& aAlloc = allocator_type()) : q(3,value_type(0),aAlloc), rowCount(3) {
      q[0] = a12;
      q[1] = a13;
      q[2] = a23;
    };

    /**
     * Constructs a 4x4 skew-symmetric matrix from six elements.
     * \test PASSED
     */
    mat(const_reference a12,const_reference a13,const_reference a14,
	                    const_reference a23,const_reference a24,
			                        const_reference a34, 
	const allocator_type& aAlloc = allocator_type()) : q(6,value_type(0),aAlloc), rowCount(4) {
      q[0] = a12;
      q[1] = a13;
      q[2] = a23;
      q[3] = a14;
      q[4] = a24;
      q[5] = a34;
    };

    /**
     * Explicit constructor of a skew-symmetric matrix from a 3D vector (cross-product matrix).
     * \throw std::range_error if the size of the vector is not 3.
     * \test PASSED
     */
    explicit mat(const vect_n<T>& V, const allocator_type& aAlloc = allocator_type()) : q(3,value_type(0),aAlloc), rowCount(3) {
      if(V.size() != 3)
	throw std::range_error("To construct a skew-matrix from a vector, that vector must have dimension 3");
      q[0] = -V[2];
      q[1] = V[1];
      q[2] = -V[0];
    };

    /**
     * Explicit constructor of a skew-symmetric matrix from a 3D vector (cross-product matrix).
     * \test PASSED
     */
    explicit mat(const vect<T,3>& V, const allocator_type& aAlloc = allocator_type()) : q(3,value_type(0),aAlloc), rowCount(3) {
      q[0] = -V[2];
      q[1] = V[1];
      q[2] = -V[0];
    };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.q,rhs.q);
      swap(lhs.rowCount,rhs.rowCount);
    };


/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    
    /**
     * Matrix indexing accessor for read-write access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \throw std::range_error if the element accessed cannot be written to (diagonal elements).
     * \test PASSED
     */
    reference operator()(size_type i,size_type j) {
      if(i > j)
	return q[mat_triangular_size(i-1) + j];
      else if(i < j)
	return q[mat_triangular_size(j-1) + i];
      else
        throw std::range_error("Cannot set the elements of the diagonal of a skew-symmetric matrix!");
    };

    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) const {
      if(i > j)
	return -q[mat_triangular_size(i-1) + j];
      else if(i < j)
        return q[mat_triangular_size(j-1) + i];
      else
	return value_type(0.0);
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const {
      return rowCount;
    };

    /**
     * Sets the row-count (number of rows) of the matrix.
     * \param aRowCount new number of rows for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_row_count(size_type aRowCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(mat_triangular_size(aRowCount-1),value_type(0));
      rowCount = aRowCount;
    };

    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const {
      return rowCount;
    };

    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_col_count(unsigned int aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(mat_triangular_size(aColCount-1),value_type(0));
      rowCount = aColCount;
    };

    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };
    
    /**
     * Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
     * \param sz new dimensions for the matrix.
     * \test PASSED
     */
    void resize(const std::pair<size_type,size_type>& sz) {
      set_row_count(sz.first,true);
    };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard Assignment operator with a symmetric matrix.
     * \test PASSED
     */
    self& operator =(self M) {
      swap(*this,M);
      return *this;
    };
    
    /**
     * Standard Assignment operator with a symmetric matrix.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator =(const Matrix& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /**
     * Add-and-store operator with standard semantics.
     * \param M the other matrix to be added to this.
     * \return this matrix by reference.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    self& operator +=(const self& M) {
      if(M.rowCount != rowCount)
	throw std::range_error("Matrix dimension mismatch.");
      typename container_type::iterator it = q.begin();
      for(typename container_type::const_iterator cit = M.q.begin(); it != q.end(); ++it, ++cit)
	*it += *cit;
      return *this;
    };

    /**
     * Sub-and-store operator with standard semantics.
     * \param M the other matrix to be substracted from this.
     * \return this matrix by reference.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    self& operator -=(const self& M) {
      if(M.rowCount != rowCount)
	throw std::range_error("Matrix dimension mismatch.");
      typename container_type::iterator it = q.begin();
      for(typename container_type::const_iterator cit = M.q.begin(); it != q.end(); ++it, ++cit)
	*it -= *cit;
      return *this;
    };

    /**
     * Scalar-multiply-and-store operator with standard semantics.
     * \param S the scalar to be multiplied to this.
     * \return this matrix by reference.
     * \test PASSED
     */
    self& operator *=(const T& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };


/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Addition operator with standard semantics.
     * \param M the other matrix to be added to this.
     * \return the matrix sum.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    friend self operator +(self M1, const self& M2) {
      M1 += M2;
      return M1;
    };

    /**
     * Negation operator with standard semantics.
     * \return the negative of this matrix.
     * \test PASSED
     */
    self operator -() const {
      self result(rowCount,value_type(0),q.get_allocator());
      typename container_type::iterator it = result.q.begin();
      for(typename container_type::const_iterator cit = q.begin(); cit != q.end(); ++it, ++cit)
	*it -= *cit;
      return result;
    };

    /**
     * Substraction operator with standard semantics.
     * \param M the other matrix to be substracted from this.
     * \return the matrix difference.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    friend self operator -(self M1, const self& M2) {
      M1 -= M2;
      return M1;
    };

    /**
     * Multiplication operator with standard semantics.
     * \param M1 the first matrix (the skew-symmetric one).
     * \param M2 the other matrix.
     * \return the matrix multiplication result, this * M.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    template <typename Matrix>
    friend 
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                (mat_product_priority<Matrix>::value < mat_product_priority<self>::value), 
    mat<value_type,mat_structure::rectangular,Alignment,Allocator> >::type 
     operator *(const self& M1,const Matrix& M2) {
      if(M1.rowCount != M2.get_row_count())
        throw std::range_error("Matrix dimension mismatch.");
      mat<value_type,mat_structure::rectangular,Alignment,Allocator> result(M1.rowCount,M2.get_col_count(),value_type(0),M1.get_allocator());
      size_type k=0; size_type i=1;
      for(;i<M1.rowCount;k += i++) {
	for(size_type l=0;l<M2.get_col_count();++l) {
	  for(size_type j=0;j<i;++j) {
	    result(j,l) += M1.q[k+j] * M2(i,l);
	    result(i,l) -= M1.q[k+j] * M2(j,l);
	  };
	};
      };
      return result;
    };
    
    /**
     * Multiplication operator with standard semantics.
     * \param M1 the first matrix.
     * \param M2 the other matrix (the skew-symmetric one).
     * \return the matrix multiplication result.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    template <typename Matrix>
    friend 
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                (mat_product_priority<Matrix>::value < mat_product_priority<self>::value), 
    mat<value_type,mat_structure::rectangular,Alignment,Allocator> >::type 
     operator *(const Matrix& M1,const self& M2) {
      if(M2.rowCount != M1.get_col_count())
        throw std::range_error("Matrix dimension mismatch.");
      mat<value_type,mat_structure::rectangular,Alignment,Allocator> result(M1.get_row_count(),M2.rowCount,value_type(0),M2.get_allocator());
      size_type k=0; size_type i=1;
      for(;i<M2.rowCount;k += i++) {
	for(size_type l=0;l<M1.get_row_count();++l) {
	  for(size_type j=0;j<i;++j) {
	    result(l,j) -= M2.q[k+j] * M1(l,i);
	    result(l,i) += M2.q[k+j] * M1(l,j);
	  };
	};
      };
      return result;
    };

    /**
     * Multiplication operator with a skew-symmetric matrix.
     * \param M the other skew-symmetric matrix to be multiplied by this.
     * \return the matrix multiplication result, this * M.
     * \throw std::range_error if the matrix dimensions don't match.
     * \test PASSED
     */
    friend 
    mat<value_type,mat_structure::square,alignment,allocator_type> operator *(const self& M1, const self& M2) {
      if(M1.rowCount != M2.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      mat<value_type,mat_structure::square,alignment,allocator_type> result(M1.rowCount,value_type(0),M1.get_allocator());
      size_type k=0; size_type i=1;
      for(;i<M1.rowCount;k += i++) {
	size_type h=0; size_type l=0;
	for(;l<=i;h += l++) {
	  size_type m=0; size_type j=0;
	  for(;j<l;m += j++) {
	    result(j,l) -= M1.q[k+j] * M2.q[k+l];
	    result(i,l) -= M1.q[k+j] * M2.q[h+j];
	  };
	  result(j,l) -= M1.q[k+j] * M2.q[k+l];
	  for(m += j++;j<i;m += j++) {
	    result(j,l) -= M1.q[k+j] * M2.q[k+l];
	    result(i,l) += M1.q[k+j] * M2.q[m+l];
	  };
	};
	for(;l<M2.rowCount;h += l++) {
	  size_type m=0; size_type j=0;
	  for(;j<i;m += j++) {
	    result(j,l) += M1.q[k+j] * M2.q[h+i];
	    result(i,l) -= M1.q[k+j] * M2.q[h+j];
	  };
	};
      };
      return result;
    };   


    /**
     * Multiplication by a column-vector.
     * \param M the matrix.
     * \param V the column vector.
     * \return column-vector equal to this * V.
     * \throw std::range_error if this matrix and the vector dimensions don't match.
     */
    template <unsigned int Size>
    friend vect<T,Size> operator *(const self& M,const vect<T,Size>& V) {
      if(M.rowCount != Size)
        throw std::range_error("Matrix dimension mismatch.");
      vect<T,Size> result;
      size_type k=0; size_type i=1;
      for(;i<Size;k += i++) {
	for(size_type j=0;j<i;++j) {
	  result[j] += M.q[k+j] * V[i];
	  result[i] -= M.q[k+j] * V[j];
	};
      };
      return result;
    };
    
    /**
     * Multiplication by a row-vector.
     * \param M the matrix.
     * \param V the row vector.
     * \return column-vector equal to this * V.
     * \throw std::range_error if this matrix and the vector dimensions don't match.
     */
    template <unsigned int Size>
    friend vect<T,Size> operator *(const vect<T,Size>& V,const self& M) {
      if(M.rowCount != Size)
        throw std::range_error("Matrix dimension mismatch.");
      vect<T,Size> result;
      size_type k=0; size_type i=1;
      for(;i<Size;k += i++) {
	for(size_type j=0;j<i;++j) {
	  result[j] -= M.q[k+j] * V[i];
	  result[i] += M.q[k+j] * V[j];
	};
      };
      return result;
    };
    
    /**
     * Multiplication by a column-vector.
     * \param M the matrix.
     * \param V the column vector.
     * \return column-vector equal to this * V.
     * \throw std::range_error if this matrix and the vector dimensions don't match.
     */
    template <typename Vector>
    friend 
    typename boost::enable_if_c< is_writable_vector<Vector>::value,
    Vector>::type operator *(const self& M,const Vector& V) {
      if(M.rowCount != V.size())
        throw std::range_error("Matrix dimension mismatch.");
      Vector result(V.size(),value_type(0),V.get_allocator());
      size_type k=0; size_type i=1;
      for(;i<V.size();k += i++) {
	for(size_type j=0;j<i;++j) {
	  result[j] += M.q[k+j] * V[i];
	  result[i] -= M.q[k+j] * V[j];
	};
      };
      return result;
    };
    
    /**
     * Multiplication by a row-vector.
     * \param M the matrix.
     * \param V the column vector.
     * \return column-vector equal to this * V.
     * \throw std::range_error if this matrix and the vector dimensions don't match.
     */
    template <typename Vector>
    friend 
    typename boost::enable_if_c< is_writable_vector<Vector>::value,
    Vector>::type operator *(const Vector& V,const self& M) {
      if(M.rowCount != V.size())
        throw std::range_error("Matrix dimension mismatch.");
      Vector result(V.size(),value_type(0),V.get_allocator());
      size_type k=0; size_type i=1;
      for(;i<V.size();k += i++) {
	for(size_type j=0;j<i;++j) {
	  result[j] -= M.q[k+j] * V[i];
	  result[i] += M.q[k+j] * V[j];
	};
      };
      return result;
    };
    
    
    
    
    
/*******************************************************************************
                         Special Methods
*******************************************************************************/



    /**
     * Extracts a sub-matrix from this matrix.
     * \param aRowOffset Number of rows before the start of the sub-matrix rows.
     * \param aColOffset Number of columns before the start of the sub-matrix columns.
     * \param aRowCountOut Number of rows of the sub-matrix.
     * \param aColCountOut Number of columns of the sub-matrix.
     * \return The sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     *//* This can just use the generic overload instead.
    mat_cm<T> getSubMat(unsigned int aRowOffset,unsigned int aColOffset,unsigned int aRowCountOut,unsigned int aColCountOut) const throw(std::range_error) {
      if((aRowOffset + aRowCountOut > size) || (aColOffset + aColCountOut > size))
        throw std::range_error("Matrix dimension mismatch.");
      mat_cm<T> result(aRowCountOut,aColCountOut);
      for(unsigned int j=0;j<aColCountOut;++j)
        for(unsigned int i=0;i < aRowCountOut;++i)
	  result.q[j*aRowCountOut+i] = this->operator()(i+aRowOffset,j+aColOffset);
      return result;
    };*/

    /**
     * Extracts a skew-symmetric sub-matrix from this matrix.
     * \param M The skew-symmetric matrix from which the sub-block is obtained.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \param aSizeOut Number of rows/columns of the sub-matrix.
     * \return The skew-symmetric sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend self get_block(const self& M,size_type aDiagOffset,size_type aSizeOut) {
      if(aDiagOffset + aSizeOut > M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      self result(aSizeOut,value_type(0),M.get_allocator());
      size_type k=mat_triangular_size(aDiagOffset);
      size_type k_out=0;
      for(size_type i=1;i<aSizeOut; k += (i + aDiagOffset), k_out += i++) 
        for(size_type j=0;j<i;++j)
          result.q[k_out+j] = M.q[k+j+aDiagOffset];
      return result;
    };


    /** Sets the sub-part of this matrix to a skew-symmetric sub-matrix M.
     * \param M A skew-symmetric sub-matrix that will be written in the sub-part of this matrix.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \return This matrix, by reference.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend self& set_block(self& M, const self& subM,size_type aDiagOffset) {
      if(aDiagOffset + subM.rowCount > M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      size_type k=mat_triangular_size(aDiagOffset);
      size_type k_in=0;
      for(size_type i=1;i<subM.rowCount;k += (i + aDiagOffset), k_in += i++) 
        for(size_type j=0;j<=i;++j) 
           M.q[k+j+aDiagOffset] = subM.q[k_in+j];
      return M;
    };
    
    
    /**
     * Appends the matrix 'rhs' to the end of the matrix 'lhs'.
     * \param lhs The matrix to which to append the other.
     * \param rhs The matrix to be appended to 'lhs'.
     */
    friend void append_block_diag(self& lhs, const self& rhs) {
      size_type oldCount = lhs.get_col_count();
      lhs.set_col_count(oldCount + rhs.get_col_count(),true);
      set_block(lhs,rhs,oldCount);
    };
    
        
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(const self& M) {
      return -M;
    };
    
    /**
     * Transposes and moves the matrix M.
     * \param M The matrix to be transposed and moved.
     * \return The transpose of M.
     */
    friend self transpose_move(self& M) {
      self result;
      swap(result,M);
      for(typename container_type::iterator it = result.q.begin(); it != result.q.end(); ++it)
	*it = -(*it);
      return result;
    };
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Transposes and moves the matrix M.
     * \param M The matrix to be transposed and moved.
     * \return The transpose of M.
     */
    friend self transpose(self&& M) {
      self result(std::move(M));
      for(typename container_type::iterator it = result.q.begin(); it != result.q.end(); ++it)
	*it = -(*it);
      return result;
    };
#endif
    
    /**
     * Returns the trace of matrix M.
     * \param M A diagonal matrix.
     * \return the trace of matrix M.
     */
    friend value_type trace(const self& M) {
      return value_type(0);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("rowCount",rowCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("rowCount",rowCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
};


#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

extern template class mat<double, mat_structure::skew_symmetric>;
extern template class mat<float, mat_structure::skew_symmetric>;


extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);

extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::skew_symmetric>& M2);
extern template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::skew_symmetric>& M2);

extern template vect<double,2> operator *(const mat<double,mat_structure::skew_symmetric>& M,const vect<double,2>& V);
extern template vect<double,3> operator *(const mat<double,mat_structure::skew_symmetric>& M,const vect<double,3>& V);
extern template vect<double,4> operator *(const mat<double,mat_structure::skew_symmetric>& M,const vect<double,4>& V);
extern template vect<double,6> operator *(const mat<double,mat_structure::skew_symmetric>& M,const vect<double,6>& V);

extern template vect<double,2> operator *(const vect<double,2>& V,const mat<double,mat_structure::skew_symmetric>& M);
extern template vect<double,3> operator *(const vect<double,3>& V,const mat<double,mat_structure::skew_symmetric>& M);
extern template vect<double,4> operator *(const vect<double,4>& V,const mat<double,mat_structure::skew_symmetric>& M);
extern template vect<double,6> operator *(const vect<double,6>& V,const mat<double,mat_structure::skew_symmetric>& M);

extern template vect_n<double> operator *(const mat<double,mat_structure::skew_symmetric>& M,const vect_n<double>& V);
extern template vect_n<double> operator *(const vect_n<double>& V,const mat<double,mat_structure::skew_symmetric>& M);


extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);

extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::skew_symmetric>& M2);
extern template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::skew_symmetric>& M2);

extern template vect<float,2> operator *(const mat<float,mat_structure::skew_symmetric>& M,const vect<float,2>& V);
extern template vect<float,3> operator *(const mat<float,mat_structure::skew_symmetric>& M,const vect<float,3>& V);
extern template vect<float,4> operator *(const mat<float,mat_structure::skew_symmetric>& M,const vect<float,4>& V);
extern template vect<float,6> operator *(const mat<float,mat_structure::skew_symmetric>& M,const vect<float,6>& V);

extern template vect<float,2> operator *(const vect<float,2>& V,const mat<float,mat_structure::skew_symmetric>& M);
extern template vect<float,3> operator *(const vect<float,3>& V,const mat<float,mat_structure::skew_symmetric>& M);
extern template vect<float,4> operator *(const vect<float,4>& V,const mat<float,mat_structure::skew_symmetric>& M);
extern template vect<float,6> operator *(const vect<float,6>& V,const mat<float,mat_structure::skew_symmetric>& M);

extern template vect_n<float> operator *(const mat<float,mat_structure::skew_symmetric>& M,const vect_n<float>& V);
extern template vect_n<float> operator *(const vect_n<float>& V,const mat<float,mat_structure::skew_symmetric>& M);


#endif




};


#endif










