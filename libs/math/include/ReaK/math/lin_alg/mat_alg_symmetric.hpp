/**
 * \file mat_alg_symmetric.hpp
 * 
 * This library implements the specialization of the mat<> template for a 
 * general symmetric matrix (dynamic dimension). This matrix type fulfills the matrix 
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

#ifndef REAK_MAT_ALG_SYMMETRIC_HPP
#define REAK_MAT_ALG_SYMMETRIC_HPP

#include "mat_alg_general.hpp"

namespace ReaK {


//....not good, the following: (ill-defined)

template <mat_alignment::tag Alignment>
struct mat_indexer<mat_structure::symmetric,Alignment> {
  std::size_t rowCount;
  mat_indexer<mat_structure::symmetric, Alignment>(std::size_t aRowCount) : rowCount(aRowCount) { };
  std::size_t mat_triangular_size(std::size_t Size) {
    return (Size * (Size - 1)) / 2 + Size;
  };
  std::size_t operator()(std::size_t i,std::size_t j) const {
    if(i > j)
      return mat_triangular_size(i) + j;
    else
      return mat_triangular_size(j) + i;
  };
};




/**
 * This class holds a symmetric matrix. This class will hold only the upper-triangular part
 * since the lower part is assumed to be equal to the upper one.
 * 
 * Models: ReadableMatrixConcept, WritableMatrixConcept, ResizableMatrixConcept, and DynAllocMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T,mat_structure::symmetric,Alignment,Allocator> : public serializable {
  public:    
    
    typedef mat<T,mat_structure::symmetric,Alignment,Allocator> self;
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
  
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::symmetric);
    
  private:
    container_type q; ///< Holds the array of scalar entries.
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
        q(0,value_type(),aAlloc),
        rowCount(0) { };

    /**
     * Constructor for a sized matrix.
     * \test PASSED
     */
    mat(size_type aRowCount, T aFill = 0,const allocator_type& aAlloc = allocator_type()) :
        q(mat_triangular_size(aRowCount),aFill,aAlloc),
        rowCount(aRowCount) { };

    /**
     * Constructor for an identity matrix.
     * \test PASSED
     */
    mat(size_type aRowCount, bool aIdentity,const allocator_type& aAlloc = allocator_type()) :
        q(mat_triangular_size(aRowCount),T(0.0),aAlloc),
        rowCount(aRowCount) {
      if(aIdentity) {
        size_type k=0;
        for(size_type i=0; i < rowCount; k += ++i)
          q[k+i] = 1.0;
      };
    };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(const self& M) : q(M.q), rowCount(M.rowCount) { };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES     
    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(self&& M) : q(std::move(M.q)), rowCount(std::move(M.rowCount)) { };
#endif

    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& M1, self& M2) throw() {
      using std::swap;
      swap(M1.q,M2.q);
      swap(M1.rowCount,M2.rowCount);
    };

    /**
     * Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
     * \test PASSED
     */    
    template <typename Matrix>
    explicit mat(const Matrix&  M, typename boost::enable_if< 
                                              boost::mpl::and_<
                                                is_readable_matrix<Matrix>,
                                                boost::mpl::not_< is_symmetric_matrix<Matrix> >,
                                                boost::mpl::not_< boost::is_same<Matrix,self> >
                                              >, void* >::type dummy = NULL) :
                 q(mat_triangular_size((M.get_row_count() > M.get_col_count() ? M.get_row_count() : M.get_col_count())),T(0.0)),
                 rowCount((M.get_row_count() > M.get_col_count() ? M.get_row_count() : M.get_col_count())) {
      size_type k=0;
      size_type i=0;
      size_type min_size = (M.get_row_count() > M.get_col_count() ? M.get_col_count() : M.get_row_count());
      for(;i<min_size;k += ++i) {
        for(size_type j=0;j<i;++j) {
          q[k+j] = value_type(0.5) * (M(j,i) + M(i,j));
        };
        q[k+i] = M(i,i);
      };
      if(M.get_row_count() > M.get_col_count()) {
        for(;i<rowCount;k += ++i) {
          for(size_type j=0;j<min_size;++j)
            q[k+j] = value_type(0.5) * M(i,j);
        };
      } else {
        for(;i<rowCount;k += ++i) {
          for(size_type j=0;j<min_size;++j)
            q[k+j] = value_type(0.5) * M(j,i);
        };
      };
    };
    
    /**
     * Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
     * \test PASSED
     */    
    template <typename Matrix>
    explicit mat(const Matrix&  M, typename boost::enable_if< 
                                              boost::mpl::and_<
                                                is_readable_matrix<Matrix>,
                                                is_symmetric_matrix<Matrix>,
                                                boost::mpl::not_< boost::is_same<Matrix,self> >
                                              >, void* >::type dummy = NULL) :
                 q(mat_triangular_size(M.get_row_count()),T(0.0)),
                 rowCount(M.get_row_count()) {
      size_type k=0;
      size_type i=0;
      for(;i<rowCount;k += ++i) {
        for(size_type j=0;j<i;++j)
          q[k+j] = M(i,j);
        q[k+i] = M(i,i);
      };
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~mat() { };

    /**
     * Constructs a 2x2 symmetric matrix from three elements.
     * \test PASSED
     */
    mat(const_reference a11,const_reference a12,
                            const_reference a22) : q(3), rowCount(2) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a22;
    };

    /**
     * Constructs a 3x3 symmetric matrix from six elements.
     * \test PASSED
     */
    mat(const_reference a11,const_reference a12,const_reference a13,
                            const_reference a22,const_reference a23,
                                                const_reference a33) : q(6), rowCount(3) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a22;
      q[3] = a13;
      q[4] = a23;
      q[5] = a33;
    };

    /**
     * Constructs a 4x4 symmetric matrix from ten elements.
     * \test PASSED
     */
    mat(const_reference a11,const_reference a12,const_reference a13,const_reference a14,
                            const_reference a22,const_reference a23,const_reference a24,
                                                const_reference a33,const_reference a34,
                                                                    const_reference a44) : q(10), rowCount(4) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a22;
      q[3] = a13;
      q[4] = a23;
      q[5] = a33;
      q[6] = a14;
      q[7] = a24;
      q[8] = a34;
      q[9] = a44;
    };


/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/


    /**
     * Matrix indexing accessor for read-write access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    reference operator()(size_type i,size_type j) {
      if(i > j)
        return q[mat_triangular_size(i) + j];
      else
        return q[mat_triangular_size(j) + i];
    };

    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    const_reference operator()(size_type i,size_type j) const {
      if(i > j)
        return q[mat_triangular_size(i) + j];
      else
        return q[mat_triangular_size(j) + i];
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
      q.resize(mat_triangular_size(aRowCount),T(0.0));
      rowCount = aRowCount;
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
      q.resize(mat_triangular_size(aColCount),T(0.0));
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
      swap(*this, M);
      return *this;
    };

    /**
     * Standard Assignment operator with a matrix of any type. The "(M + M.transpose) / 2" formula is applied to guarantee symmetry.
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
     * \test PASSED
     */
    template <mat_alignment::tag Align2, typename Alloc2>
    self& operator +=(const mat<T,mat_structure::symmetric,Align2,Alloc2>& M) {
      if(M.get_row_count() != rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      size_type k=0;
      for(size_type i=0; i < rowCount; k += ++i)
        for(size_type j=0;j <= i;++j)
          q[k+j] += M(i,j);
      return *this;
    };
    
    /**
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <mat_alignment::tag Align2, typename Alloc2>
    self& operator +=(const mat<T,mat_structure::diagonal,Align2,Alloc2>& M) {
      if(M.get_row_count() != rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      size_type k=0;
      for(size_type i=0; i < rowCount; k += ++i)
        q[k+i] += M(i,i);
      return *this;
    };

    /**
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <mat_alignment::tag Align2, typename Alloc2>
    self& operator -=(const mat<T,mat_structure::symmetric,Align2,Alloc2>& M) {
      if(M.get_row_count() != rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      unsigned int k=0;
      for(size_type i=0; i < rowCount; k += ++i)
        for(size_type j=0;j <= i;++j)
          q[k+j] -= M(i,j);
      return *this;
    };
    
    /**
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <mat_alignment::tag Align2, typename Alloc2>
    self& operator -=(const mat<T,mat_structure::diagonal,Align2,Alloc2>& M) {
      if(M.get_row_count() != rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      size_type k=0;
      for(size_type i=0; i < rowCount; k += ++i)
        q[k+i] -= M(i,i);
      return *this;
    };

    /**
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const T& S) {
      for(typename container_type::iterator it = q.begin(); it != q.end(); ++it) 
        *it *= S;
      return *this;
    };
    
    /** 
     * Negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return Symmetric matrix.
     * \test PASSED
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin();it != q.end();++it,++itr)
        *itr = -(*it);
      return result;
    };


/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Add two matrices.
     * \param M1 the first matrix (first operand).
     * \param M2 the other matrix (second operand).
     * \return symmetric matrix sum of M1 and M2.
     * \throw std::range_error if the two matrix dimensions don't match.
     * \test PASSED
     */
    friend self operator +(const self& M1,const self& M2) {
      if(M1.rowCount != M2.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      self result(M1.rowCount);
      size_type k=0;
      for(size_type i=0; i < M1.rowCount; k += ++i)
        for(size_type j=0;j <= i;++j)
          result.q[k+j] = M1.q[k+j] + M2.q[k+j];
      return result;
    };

    /**
     * Sub two matrices.
     * \param M1 the first matrix (first operand).
     * \param M2 the other matrix (second operand).
     * \return symmetric matrix difference of M1 and M2.
     * \throw std::range_error if the two matrix dimensions don't match.
     * \test PASSED
     */
    friend self operator -(const self& M1,const self& M2) {
      if(M1.rowCount != M2.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      self result(M1.rowCount);
      size_type k=0;
      for(size_type i=0; i < M1.rowCount; k += ++i)
        for(size_type j=0;j <= i;++j)
          result.q[k+j] = M1.q[k+j] - M2.q[k+j];
      return result;
    };

    /**
     * General Matrix multiplication.
     * \param M the other matrix (second operand).
     * \return general matrix multiplication result of this and M.
     * \throw std::range_error if the two matrix dimensions don't match.
     */
    template <typename Matrix>
    friend 
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                (mat_product_priority<Matrix>::value < mat_product_priority<self>::value), 
    mat<value_type,mat_structure::rectangular,Alignment,Allocator> >::type 
     operator *(const self& M1, const Matrix& M2) {
      typedef mat<value_type,mat_structure::rectangular,Alignment,Allocator> result_type;
      if(M1.rowCount != M2.get_row_count())
        throw std::range_error("Matrix dimension mismatch.");
      result_type result(M1.rowCount,M2.get_col_count(),value_type(0),M1.get_allocator());
      size_type k=0; size_type i=0;
      for(;i<M1.rowCount;k += ++i) {
        for(size_type l=0;l<M2.get_col_count();++l) {
          for(size_type j=0;j<i;++j) {
            result(j,l) += M1.q[k+j] * M2(i,l);
            result(i,l) += M1.q[k+j] * M2(j,l);
          };
          result(i,l) += M1.q[k+i] * M2(i,l);
        };
      };
      return result;
    };
    
    /**
     * General Matrix multiplication.
     * \param M the other matrix (second operand).
     * \return general matrix multiplication result of this and M.
     * \throw std::range_error if the two matrix dimensions don't match.
     */
    template <typename Matrix>
    friend 
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                (mat_product_priority<Matrix>::value < mat_product_priority<self>::value), 
    mat<value_type,mat_structure::rectangular,alignment,allocator_type> >::type 
     operator *(const Matrix& M1, const self& M2) {
      typedef mat<typename mat_traits<Matrix>::value_type,mat_structure::rectangular,alignment,allocator_type> result_type;
      if(M2.rowCount != M1.get_col_count())
        throw std::range_error("Matrix dimension mismatch.");
      result_type result(M1.get_row_count(),M2.rowCount,value_type(0),M2.get_allocator());
      size_type k=0; size_type i=0;
      for(;i<M2.rowCount;k += ++i) {
        for(size_type l=0;l<M1.get_row_count();++l) {
          for(size_type j=0;j<i;++j) {
            result(l,i) += M2.q[k+j] * M1(l,j);
            result(l,j) += M2.q[k+j] * M1(l,i);
          };
          result(l,i) += M2.q[k+i] * M1(l,i);
        };
      };
      return result;
    };

    /**
     * Symmetric Matrix multiplication.
     * \param M1 the first symmetric matrix (first operand).
     * \param M2 the other symmetric matrix (second operand).
     * \return square matrix, result of M1 times M2.
     * \throw std::range_error if the two matrix dimensions don't match.
     */
    friend mat<T,mat_structure::square,alignment,allocator_type> 
      operator *(const self& M1,const self& M2) {
      typedef mat<T,mat_structure::square,alignment,allocator_type> result_type;
      if(M1.rowCount != M2.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      result_type result(M1.rowCount,T(0),M1.get_allocator());
      size_type k=0; size_type i=0;
      for(;i<M1.rowCount;k += ++i) {
        size_type h=0; size_type l=0;
        for(;l<=i;h += ++l) {
          size_type m=0; size_type j=0;
          for(;j<l;m += ++j) {
            result(j,l) += M1.q[k+j] * M2.q[k+l];
            result(i,l) += M1.q[k+j] * M2.q[h+j];
          };
          for(;j<i;m += ++j) {
            result(j,l) += M1.q[k+j] * M2.q[k+l];
            result(i,l) += M1.q[k+j] * M2.q[m+l];
          };
          result(i,l) += M1.q[k+i] * M2.q[k+l];
        };
        for(;l<M2.rowCount;h += ++l) {
          size_type m=0; size_type j=0;
          for(;j<i;m += ++j) {
            result(j,l) += M1.q[k+j] * M2.q[h+i];
            result(i,l) += M1.q[k+j] * M2.q[h+j];
          };
          result(i,l) += M1.q[k+i] * M2.q[h+i];
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
      if(Size != M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      vect<T,Size> result;
      size_type k=0; size_type i=0;
      for(;i<Size;k += ++i) {
        for(size_type j=0;j<i;++j) {
          result[i] += M.q[k+j] * V[j];
          result[j] += M.q[k+j] * V[i];
        };
        result[i] += M.q[k+i] * V[i];
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
      return M * V;
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
      if(V.size() != M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      Vector result(V.size(),value_type(0),V.get_allocator());
      size_type k=0; size_type i=0;
      for(;i<V.size();k += ++i) {
        for(size_type j=0;j<i;++j) {
          result[i] += M.q[k+j] * V[j];
          result[j] += M.q[k+j] * V[i];
        };
        result[i] += M.q[k+i] * V[i];
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
      return M * V;
    };
    
/*******************************************************************************
                         Special Methods
*******************************************************************************/
    /**
     * Extracts a sub-matrix from this matrix.
     * \param M A symmetric matrix.
     * \param aRowOffset Number of rows before the start of the sub-matrix rows.
     * \param aColOffset Number of columns before the start of the sub-matrix columns.
     * \param aRowCountOut Number of rows of the sub-matrix.
     * \param aColCountOut Number of columns of the sub-matrix.
     * \return The sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend 
    mat<T,mat_structure::rectangular,Alignment,Allocator> get_block(const self& M,
                                                                    size_type aRowOffset,size_type aColOffset,
                                                                    size_type aRowCountOut,size_type aColCountOut) {
      if((aRowOffset + aRowCountOut > M.rowCount) || (aColOffset + aColCountOut > M.rowCount))
        throw std::range_error("Matrix dimension mismatch.");
      mat<T,mat_structure::rectangular,Alignment,Allocator> result(aRowCountOut,aColCountOut,T(0),M.get_allocator());
      size_type k=mat_triangular_size(aColOffset);
      for(size_type j=0;j<aColCountOut;k += (++j + aColOffset)) {
        size_type h=mat_triangular_size(aRowOffset); size_type i=0;
        for(;((i<aRowCountOut) && (i+aRowOffset <= j+aColOffset));h += (++i + aRowOffset))
          result(i,j) = M.q[k+i+aRowOffset];
        for(;i<aRowCountOut;h += (++i + aRowOffset)) {
          result(i,j) = M.q[h+j+aColOffset];
        };
      };
      return result;
    };
    
    /**
     * Extracts a sub-matrix from this matrix.
     * \param M A symmetric matrix.
     * \param aRowOffset Number of rows before the start of the sub-matrix rows.
     * \param aColOffset Number of columns before the start of the sub-matrix columns.
     * \param aSizeOut Number of rows and columns of the sub-matrix.
     * \return The sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend 
    mat<T,mat_structure::square,Alignment,Allocator> get_block(const self& M,
                                                               size_type aRowOffset,size_type aColOffset,
                                                               size_type aSizeOut) {
      if((aRowOffset + aSizeOut > M.rowCount) || (aColOffset + aSizeOut > M.rowCount))
        throw std::range_error("Matrix dimension mismatch.");
      mat<T,mat_structure::square,Alignment,Allocator> result(aSizeOut,T(0),M.get_allocator());
      size_type k=mat_triangular_size(aColOffset);
      for(size_type j=0;j<aSizeOut;k += (++j + aColOffset)) {
        size_type h=mat_triangular_size(aRowOffset); size_type i=0;
        for(;((i<aSizeOut) && (i+aRowOffset <= j+aColOffset));h += (++i + aRowOffset))
          result(i,j) = M.q[k+i+aRowOffset];
        for(;i<aSizeOut;h += (++i + aRowOffset)) {
          result(i,j) = M.q[h+j+aColOffset];
        };
      };
      return result;
    };

    /**
     * Extracts a symmetric sub-matrix from this matrix.
     * \param M A symmetric matrix.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \param aSizeOut Number of rows/columns of the sub-matrix.
     * \return The symmetric sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend 
    self get_block(const self& M, size_type aDiagOffset,size_type aSizeOut) {
      if(aDiagOffset + aSizeOut > M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      self result(aSizeOut);
      size_type k=mat_triangular_size(aDiagOffset);
      size_type k_out=0;
      for(size_type i=0;i<aSizeOut;k_out += ++i, k += (i + aDiagOffset))
        for(size_type j=0;j<=i;++j)
          result.q[k_out+j] = M.q[k+j+aDiagOffset];
      return result;
    };


    /** Sets the sub-part of this matrix to a symmetric sub-matrix M.
     * \param M A symmetric sub-matrix that will be written in the sub-part of this matrix.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \return This matrix, by reference.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    friend
    self& set_block(self& M, const self& subM,size_type aDiagOffset) {
      if(aDiagOffset + subM.rowCount > M.rowCount)
        throw std::range_error("Matrix dimension mismatch.");
      size_type k=mat_triangular_size(aDiagOffset);
      size_type k_in=0;
      for(size_type i=0;i<subM.rowCount;k_in += ++i, k += (i + aDiagOffset))
        for(size_type j=0;j<=i;++j)
          M.q[k+j+aDiagOffset] = subM.q[k_in+j];
      return M;
    };
    
    /**
     * Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both symmetric matrices.
     * \param lhs The symmetric matrix to which to append the other.
     * \param rhs The symmetric matrix to be appended to 'lhs'.
     */
    friend void append_block_diag(self& lhs, const self& rhs) {
      size_type oldCount = lhs.get_col_count();
      lhs.set_col_count(oldCount + rhs.get_col_count(),true);
      set_block(lhs,rhs,oldCount);
    };
    
    /**
     * Transposes the matrix M (which has no effect since M is symmetric, simply copies it).
     * \param M The symmetric matrix to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(const self& M) {
      return M;
    };
    
    /**
     * Transposes the matrix M in a potentially destructive way (move-semantics, pre-C++0x).
     * \param M The symmetric matrix to be transposed and moved.
     * \return The transpose of M.
     */
    friend self transpose_move(self& M) {
      self result;
      swap(result,M);
      return result;
    };
    
    /**
     * Returns the trace of matrix M.
     * \param M A symmetric matrix.
     * \return the trace of matrix M.
     */
    friend value_type trace(const self& M) {
      value_type sum = value_type(0);
      size_type k = 0;
      for(size_type i = 0; i < M.rowCount; k += ++i)
        sum += M.q[k+i];
      return sum;
    };

  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("rowCount",rowCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int temp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("rowCount",temp);
      rowCount = temp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serializable)
  
};



#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template class mat<double, mat_structure::symmetric>;
extern template class mat<float, mat_structure::symmetric>;


extern template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator +=(const mat<double,mat_structure::symmetric>& M);
extern template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator +=(const mat<double,mat_structure::diagonal>& M);
extern template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator -=(const mat<double,mat_structure::symmetric>& M);
extern template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator -=(const mat<double,mat_structure::diagonal>& M);

extern template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator +=(const mat<float,mat_structure::symmetric>& M);
extern template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator +=(const mat<float,mat_structure::diagonal>& M);
extern template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator -=(const mat<float,mat_structure::symmetric>& M);
extern template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator -=(const mat<float,mat_structure::diagonal>& M);


#endif



};


#endif











