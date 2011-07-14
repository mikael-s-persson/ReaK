
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

#ifndef MAT_ALG_IDENTITY_HPP
#define MAT_ALG_IDENTITY_HPP

#include "mat_alg_general.hpp"
#include "mat_alg_nil.hpp"


namespace ReaK {


  
  


  
/**
 * This class implements a place-holder or interface-implementation to represent
 * a identity matrix (all entries zero except diagonal). This is useful to build for example a
 * block-matrix with some identity-matrix blocks.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T,mat_structure::identity, Alignment, Allocator> : public serialization::serializable {
  public:    
    
    typedef mat<T,mat_structure::identity, Alignment, Allocator> self;
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
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::identity);
    
  private:
    size_type rowCount; ///< Row Count.
  public:
    /**
     * Default constructor. Sets dimensions to zero.
     */
    mat() : rowCount(0) { };
    /**
     * Constructs a null matrix to the given dimensions.
     */
    mat(size_type aRowCount) : rowCount(aRowCount) { };
    
    mat(const self& rhs) : rowCount(rhs.rowCount) { };
    /**
     * Default destructor.
     */
    ~mat() { };
    
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
      swap(lhs.rowCount,rhs.rowCount);
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    const_reference operator()(size_type i,size_type j) const { if(i == j) return value_type(1); else return value_type(0); };

    size_type get_row_count() const { return rowCount; };
    size_type get_col_count() const { return rowCount; };
    
    /**
     * Negate the matrix.
     * \return This matrix, by constant reference.
     * \test PASSED
     */
    mat<value_type,mat_structure::diagonal> operator -() const {
      return mat<value_type,mat_structure::diagonal>(rowCount,value_type(-1));
    };
    
    
    friend self transpose(self rhs) {
      return rhs;
    };
    
    friend self transpose_move(self& rhs) {
      self result; 
      swap(result, rhs);
      return result;
    };
    
    friend value_type trace(const self& M) {
      return M.rowCount;
    };
    
    friend void append_block_diag(self& lhs, const self& rhs) {
      lhs.rowCount += rhs.rowCount;
    };
        
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, unsigned int>("rowCount",rowCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & std::pair<std::string, unsigned int&>("rowCount",rowCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)

};


template <typename T>
struct mat_identity {
  typedef mat<T,mat_structure::identity> type;
};

template <typename T>
mat<T,mat_structure::identity> mat_ident(typename mat<T,mat_structure::identity>::size_type aRowCount) {
  return mat<T,mat_structure::identity>(aRowCount);
};



template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
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
 Vector >::type
  operator *(const mat<T,mat_structure::identity,Alignment,Allocator>& M, const Vector& V) {
    if(V.size() != M.get_col_count())
      throw std::range_error("Matrix dimension mismatch.");
    return V;
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
 Vector >::type 
  operator *(const Vector& V,const mat<T,mat_structure::identity,Alignment,Allocator>& M) {
    if(V.size() != M.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    return V;
  };

/**
 * Column-vector multiplication, always results in a null vector.
 * \param M some matrix.
 * \param S some scalar.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< !is_readable_vector<T>::value &&
                             !is_readable_matrix<T>::value, 
 mat<T,mat_structure::diagonal,Alignment,Allocator> >::type
  operator *(const mat<T,mat_structure::identity,Alignment,Allocator>& M, const T& S) {
    return mat<T,mat_structure::diagonal,Alignment,Allocator>(M.get_row_count(),S);
  };
    
/**
 * Row-vector multiplication with null matrix, always results in a null vector.
 * \param V some row-vector.
 * \param M a null-matrix.
 * \return A null vector, by value.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
typename boost::enable_if_c< !is_readable_vector<T>::value &&
                             !is_readable_matrix<T>::value, 
 mat<T,mat_structure::diagonal,Alignment,Allocator> >::type 
  operator *(const T& S,const mat<T,mat_structure::identity,Alignment,Allocator>& M) {
    return mat<T,mat_structure::diagonal,Alignment,Allocator>(M.get_row_count(),S);
  };

/**
 * Matrix multiplication with identity matrix.
 * \param M1 some matrix.
 * \param M2 a identity-matrix.
 * \return A matrix.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< mat_product_priority< mat<T,Structure,Alignment,Allocator> >::value <= mat_product_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value,
 const mat<T,Structure,Alignment,Allocator>& >::type
  operator *(const mat<T,Structure,Alignment,Allocator>& M1, const mat<T,mat_structure::identity,Alignment2,Allocator2>& M2) {
    if(M1.get_col_count() != M2.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    return M1;
  };

/**
 * Matrix multiplication, always results in a identity matrix.
 * \param M1 a identity-matrix.
 * \param M2 some matrix.
 * \return A matrix.
 * \throw std::range_error if matrix dimensions are not proper for multiplication.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< mat_product_priority< mat<T,Structure,Alignment,Allocator> >::value < mat_product_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value,
 const mat<T,Structure,Alignment,Allocator>& >::type 
  operator *(const mat<T,mat_structure::identity,Alignment2,Allocator2>& M1, const mat<T,Structure,Alignment,Allocator>& M2) {
    if(M2.get_row_count() != M1.get_col_count())
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
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< is_writable_matrix< mat<T,Structure,Alignment,Allocator> >::value &&
                             (mat_addition_priority< mat<T,Structure,Alignment,Allocator> >::value < mat_addition_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value), 
 mat<T,Structure,Alignment,Allocator> >::type
  operator +(const mat<T,mat_structure::identity,Alignment2,Allocator2>& M1, mat<T,Structure,Alignment,Allocator> M2) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    for(std::size_t i=0;i<M1.get_row_count();++i)
      M2(i,i) += T(1);
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
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< is_writable_matrix< mat<T,Structure,Alignment,Allocator> >::value && 
                             (mat_addition_priority< mat<T,Structure,Alignment,Allocator> >::value < mat_addition_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value), 
 mat<T,Structure,Alignment,Allocator> >::type
  operator +(mat<T,Structure,Alignment,Allocator> M2, const mat<T,mat_structure::identity,Alignment2,Allocator2>& M1) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    for(std::size_t i=0;i<M1.get_row_count();++i)
      M2(i,i) += T(1);
    return M2;
  };
  
/**
 * Add two identity matrices together.
 * \param M1 some identity matrix.
 * \param M2 some identity matrix.
 * \return diagonal matrix.
 * \throw std::range_error if the matrix dimensions are not proper for addition.
 * 
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
mat<T,mat_structure::diagonal,Alignment,Allocator> operator+(const mat<T,mat_structure::identity,Alignment,Allocator>& M1,const mat<T,mat_structure::identity,Alignment2,Allocator2>& M2) {
  if(M1.get_row_count() != M2.get_row_count())
    throw std::range_error("Matrix dimension mismatch.");
  return mat<T,mat_structure::diagonal,Alignment,Allocator>(M1.get_row_count(),T(2));
};
  
  
/**
 * Sub two matrices, just returns the negative of the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The substraction of the matrices (actually just returns -M), by value.
 * \throw std::range_error if matrix dimensions are not proper for substraction.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< is_writable_matrix< mat<T,Structure,Alignment,Allocator> >::value && 
                             (mat_addition_priority< mat<T,Structure,Alignment,Allocator> >::value < mat_addition_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value), 
 mat<T,Structure,Alignment,Allocator> >::type
  operator -(const mat<T,mat_structure::identity,Alignment2,Allocator2>& M1, const mat<T,Structure,Alignment,Allocator>& M2) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    mat<T,Structure,Alignment,Allocator> result = -M2;
    for(std::size_t i=0;i<M1.get_row_count();++i)
      result(i,i) += T(1);
    return result;
  };
  
/**
 * Sub two matrices, just returns the negative of the parameter.
 * \param M1 some matrix.
 * \param M2 some matrix.
 * \return The substraction of the matrices (actually just returns -M), by value.
 * \throw std::range_error if matrix dimensions are not proper for substraction.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
typename boost::enable_if_c< is_writable_matrix< mat<T,Structure,Alignment,Allocator> >::value && 
                             (mat_addition_priority< mat<T,Structure,Alignment,Allocator> >::value < mat_addition_priority< mat<T,mat_structure::identity,Alignment2,Allocator2> >::value), 
 mat<T,Structure,Alignment,Allocator> >::type
  operator -(mat<T,Structure,Alignment,Allocator> M2, const mat<T,mat_structure::identity,Alignment2,Allocator2>& M1) {
    if((M2.get_row_count() != M1.get_row_count()) || (M2.get_col_count() != M1.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    for(std::size_t i=0;i<M1.get_row_count();++i)
      M2(i,i) -= T(1);
    return M2;
  };
  
/**
 * Sub two identity matrices together.
 * \param M1 some identity matrix.
 * \param M2 some identity matrix.
 * \return null matrix.
 * \throw std::range_error if the matrix dimensions are not proper for addition.
 * 
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator, mat_alignment::tag Alignment2, typename Allocator2>
 mat<T,mat_structure::nil,Alignment,Allocator> operator-(const mat<T,mat_structure::identity,Alignment,Allocator>& M1,const mat<T,mat_structure::identity,Alignment2,Allocator2>& M2) {
  if(M1.get_row_count() != M2.get_row_count())
    throw std::range_error("Matrix dimension mismatch.");
  return mat<T,mat_structure::nil,Alignment,Allocator>(M1.get_row_count(),M1.get_row_count());
};
  
  
  







};


#endif










