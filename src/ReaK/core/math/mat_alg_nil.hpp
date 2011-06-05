
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

#ifndef MAT_ALG_NIL_HPP
#define MAT_ALG_NIL_HPP

#include "mat_alg_general.hpp"


namespace ReaK {




  

/**
 * This class implements a place-holder or interface-implementation to represent
 * a null matrix (all entries zero). This is useful to build for example a
 * block-matrix with some zero-matrix blocks.
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
  
    typedef unsigned int size_type;
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

    const_reference operator()(size_type i,size_type j) const { return T(0.0); };

    size_type get_row_count() const { return rowCount; };
    size_type get_col_count() const { return colCount; };
    
    /**
     * Negate the matrix, has no effect of course.
     * \return This matrix, by constant reference.
     * \test PASSED
     */
    const self& operator -() const {
      return *this;
    };
    
    
    friend self transpose(self rhs) {
      using std::swap;
      swap(rhs.colCount,rhs.rowCount);
      return rhs;
    };
    friend self transpose_move(self& rhs) {
      self result(rhs.colCount,rhs.rowCount);
      return result;
    };
    
    friend value_type trace(const self&) {
      return value_type(0);
    };
    
    friend void append_block_diag(self& lhs, const self& rhs) {
      lhs.colCount += rhs.colCount;
      lhs.rowCount += rhs.rowCount;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, unsigned int>("rowCount",rowCount)
	& std::pair<std::string, unsigned int>("colCount",colCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & std::pair<std::string, unsigned int&>("rowCount",rowCount)
	& std::pair<std::string, unsigned int&>("colCount",colCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
    
};

template <typename T>
struct mat_null {
  typedef mat<T,mat_structure::nil> type;
};


template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix< mat<T,mat_structure::nil, Alignment, Allocator> > {
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
  






};

#endif











