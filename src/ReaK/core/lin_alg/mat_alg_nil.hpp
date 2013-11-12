/**
 * \file mat_alg_nil.hpp
 * 
 * This library declares matrix specializations for representing and manipulating nil matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented 
 * if specialized for the nil matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the nil matrix class is simply a partial specialization of the "ReaK::mat" 
 * class template, so, the burden on the user is minimal.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
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

#ifndef REAK_MAT_ALG_NIL_HPP
#define REAK_MAT_ALG_NIL_HPP

#include "mat_alg_general.hpp"


namespace ReaK {




  

/**
 * This class implements a place-holder or interface-implementation to represent
 * a nil matrix (all entries zero). This is useful to build for example a
 * block-matrix with some zero-matrix blocks, and, of course, the storage is minimal.
 * 
 * Models: ReadableMatrixConcept and ResizableMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T,mat_structure::nil,Alignment,Allocator> : public serialization::serializable {
  public:    
    
    typedef mat<T,mat_structure::nil,Alignment,Allocator> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef void reference;
    typedef T const_reference;
    typedef void pointer;
    typedef void const_pointer;
  
    typedef void row_iterator;
    typedef void const_row_iterator;
    typedef void col_iterator;
    typedef void const_col_iterator;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::nil);
    
  private:
    size_type rowCount; ///< Row Count.
    size_type colCount; ///< Column Count.
  public:
    /**
     * Default constructor. Sets dimensions to zero.
     */
    mat() : rowCount(0), colCount(0) { };
    /**
     * Constructs a null matrix to the given dimensions.
     */
    mat(size_type aRowCount,size_type aColCount) : rowCount(aRowCount), colCount(aColCount) { };
    
    mat(const self& rhs) : rowCount(rhs.rowCount), colCount(rhs.colCount) { };
    /**
     * Default destructor.
     */
    ~mat() { };
    
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
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
    const_reference operator()(size_type i,size_type j) const { return T(0.0); };

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
    };

    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const { return colCount; };

    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_col_count(size_type aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      colCount = aColCount;
    };
    
    
    /**
     * Negate the matrix, has no effect of course.
     * \return This matrix, by constant reference.
     * \test PASSED
     */
    const self& operator -() const {
      return *this;
    };
    
    /**
     * Transposes the matrix M.
     * \param rhs The nil matrix to be transposed.
     * \return The transpose of rhs.
     */
    friend self transpose(self rhs) {
      using std::swap;
      swap(rhs.colCount,rhs.rowCount);
      return rhs;
    };
    
    /**
     * Transposes the matrix M.
     * \param rhs The nil matrix to be transposed.
     * \return The transpose of rhs.
     */
    friend self transpose_move(self& rhs) {
      self result(rhs.colCount,rhs.rowCount);
      return result;
    };
    
    /**
     * Returns the trace of the matrix.
     * \return the trace of the matrix.
     */
    friend value_type trace(const self&) {
      return value_type(0);
    };
    
    /**
     * Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both nil matrices.
     * \param lhs The nil matrix to which to append the other.
     * \param rhs The nil matrix to be appended to 'lhs'.
     */
    friend void append_block_diag(self& lhs, const self& rhs) {
      lhs.colCount += rhs.colCount;
      lhs.rowCount += rhs.rowCount;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(rowCount)
	& RK_SERIAL_SAVE_WITH_NAME(colCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(rowCount)
	& RK_SERIAL_LOAD_WITH_NAME(colCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
    
};

template <typename T>
struct mat_null {
  typedef mat<T,mat_structure::nil> type;
};


template <typename T>
mat<T,mat_structure::nil> mat_nil(typename mat<T,mat_structure::nil>::size_type aRowCount,
				  typename mat<T,mat_structure::nil>::size_type aColCount) {
  return mat<T,mat_structure::nil>(aRowCount,aColCount);
};


template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_fully_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};


/**
 * Column-vector multiplication, always results in a null vector.
 * \param M some matrix.
 * \param V some vector.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_vector<Vector>::value, 
 vect_n< T, Allocator > >::type
  operator *(const mat<T,mat_structure::nil,Alignment,Allocator>& M, const Vector& V) {
    if(V.size() != M.get_col_count())
      throw std::range_error("Matrix dimension mismatch.");
    return vect_n< T, Allocator >(M.get_row_count());
  };
    
/**
 * Row-vector multiplication with null matrix, always results in a null vector.
 * \param V some row-vector.
 * \param M a null-matrix.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_vector<Vector>::value, 
 vect_n< T, Allocator > >::type 
  operator *(const Vector& V,const mat<T,mat_structure::nil,Alignment,Allocator>& M) {
    if(V.size() != M.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    return vect_n< T, Allocator >(M.get_col_count());
  };
  
/**
 * Multiplication by a column-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return column-vector equal to M * V.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, unsigned int Size, mat_alignment::tag Alignment, typename Allocator>
vect<T,Size> operator *(const mat<T,mat_structure::nil,Alignment,Allocator>& M,const vect<T,Size>& V) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  return vect<T,Size>();
};


/**
 * Multiplication by a row-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return row-vector equal to V * M.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, unsigned int Size, mat_alignment::tag Alignment, typename Allocator>
vect<T,Size> operator *(const vect<T,Size>& V,const mat<T,mat_structure::nil,Alignment,Allocator>& M) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  return vect<T,Size>();
};

#if 0
/**
 * Matrix multiplication with null matrix, always results in a null matrix.
 * \param M1 some matrix.
 * \param M2 a null-matrix.
 * \return A null matrix, by value.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, typename Matrix, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
 mat<T,mat_structure::nil,Alignment,Allocator> >::type
  operator *(const Matrix& M1, const mat<T,mat_structure::nil,Alignment,Allocator>& M2) {
    if(M1.get_col_count() != M2.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    return mat<T,mat_structure::nil,Alignment,Allocator>(M1.get_row_count(),M2.get_col_count());
  };

/**
 * Matrix multiplication, always results in a null matrix.
 * \param M some matrix.
 * \return A null matrix, by value.
 * \throw std::range_error if matrix dimensions are not proper for multiplication.
 */
template <typename T, typename Matrix, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
 mat<T,mat_structure::nil,Alignment,Allocator> >::type 
  operator *(const mat<T,mat_structure::nil,Alignment,Allocator>& M1, const Matrix& M2) {
    if(M2.get_row_count() != M1.get_col_count())
      throw std::range_error("Matrix dimension mismatch.");
    return mat<T,mat_structure::nil,Alignment,Allocator>(M1.get_row_count(),M2.get_col_count());
  };
#endif

#if 0
  
/**
 * Add two matrices, just returns the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The sum of the matrices (actually just returns M), by constant reference.
 * \throw std::range_error if matrix dimensions are not proper for addition.
 * \test PASSED
 */
template <typename T, typename Matrix, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, const Matrix& >::type
  operator +(const mat<T,mat_structure::nil,Alignment,Allocator>& M1, const Matrix& M2) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    return M2;
  };
    
/**
 * Add two matrices, just returns the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The sum of the matrices (actually just returns M), by constant reference.
 * \throw std::range_error if matrix dimensions are not proper for addition.
 * \test PASSED
 */
template <typename Matrix, typename T, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, const Matrix& >::type
  operator +(const Matrix& M2, const mat<T,mat_structure::nil,Alignment,Allocator>& M1) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    return M2;
  };
  
  
/**
 * Sub two matrices, just returns the negative of the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The substraction of the matrices (actually just returns -M), by value.
 * \throw std::range_error if matrix dimensions are not proper for substraction.
 * \test PASSED
 */
template <typename T, typename Matrix, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, Matrix >::type
  operator -(const mat<T,mat_structure::nil,Alignment,Allocator>& M1, const Matrix& M2) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    return -M2;
  };
  
/**
 * Sub two matrices, just returns the negative of the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The substraction of the matrices (actually just returns -M), by value.
 * \throw std::range_error if matrix dimensions are not proper for substraction.
 * \test PASSED
 */
template <typename T, typename Matrix, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, const Matrix& >::type
  operator -(const Matrix& M2, const mat<T,mat_structure::nil,Alignment,Allocator>& M1) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    return M2;
  };
  

#endif



#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template class mat<double, mat_structure::nil>;
extern template class mat<float, mat_structure::nil>;

extern template mat<double,mat_structure::nil> mat_nil<double>(mat<double,mat_structure::nil>::size_type aRowCount, mat<double,mat_structure::nil>::size_type aColCount);
extern template mat<float,mat_structure::nil> mat_nil<float>(mat<float,mat_structure::nil>::size_type aRowCount, mat<float,mat_structure::nil>::size_type aColCount);

extern template vect<double,2> operator *(const mat<double,mat_structure::nil>& M,const vect<double,2>& V);
extern template vect<double,3> operator *(const mat<double,mat_structure::nil>& M,const vect<double,3>& V);
extern template vect<double,4> operator *(const mat<double,mat_structure::nil>& M,const vect<double,4>& V);
extern template vect<double,6> operator *(const mat<double,mat_structure::nil>& M,const vect<double,6>& V);
extern template vect_n<double> operator *(const mat<double,mat_structure::nil>& M, const vect_n<double>& V);

extern template vect<double,2> operator *(const vect<double,2>& V,const mat<double,mat_structure::nil>& M);
extern template vect<double,3> operator *(const vect<double,3>& V,const mat<double,mat_structure::nil>& M);
extern template vect<double,4> operator *(const vect<double,4>& V,const mat<double,mat_structure::nil>& M);
extern template vect<double,6> operator *(const vect<double,6>& V,const mat<double,mat_structure::nil>& M);
extern template vect_n<double> operator *(const vect_n<double>& V,const mat<double,mat_structure::nil>& M);

extern template vect<float,2> operator *(const mat<float,mat_structure::nil>& M,const vect<float,2>& V);
extern template vect<float,3> operator *(const mat<float,mat_structure::nil>& M,const vect<float,3>& V);
extern template vect<float,4> operator *(const mat<float,mat_structure::nil>& M,const vect<float,4>& V);
extern template vect<float,6> operator *(const mat<float,mat_structure::nil>& M,const vect<float,6>& V);
extern template vect_n<float> operator *(const mat<float,mat_structure::nil>& M, const vect_n<float>& V);

extern template vect<float,2> operator *(const vect<float,2>& V,const mat<float,mat_structure::nil>& M);
extern template vect<float,3> operator *(const vect<float,3>& V,const mat<float,mat_structure::nil>& M);
extern template vect<float,4> operator *(const vect<float,4>& V,const mat<float,mat_structure::nil>& M);
extern template vect<float,6> operator *(const vect<float,6>& V,const mat<float,mat_structure::nil>& M);
extern template vect_n<float> operator *(const vect_n<float>& V,const mat<float,mat_structure::nil>& M);

#endif



};

#endif











