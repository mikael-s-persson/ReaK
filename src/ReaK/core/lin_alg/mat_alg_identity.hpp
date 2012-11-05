/**
 * \file mat_alg_identity.hpp
 * 
 * This library declares matrix specializations for representing and manipulating identity matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented 
 * if specialized for the identity matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the identity matrix class is simply a partial specialization of the "ReaK::mat" 
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

#ifndef REAK_MAT_ALG_IDENTITY_HPP
#define REAK_MAT_ALG_IDENTITY_HPP

#include "mat_alg_general.hpp"
#include "mat_alg_nil.hpp"


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
  
    typedef std::size_t size_type;
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
     * Constructs an identity matrix to the given dimensions.
     */
    mat(size_type aRowCount) : rowCount(aRowCount) { };
    
    mat(const self& rhs) : rowCount(rhs.rowCount) { };
    /**
     * Default destructor.
     */
    ~mat() { };
    
    /**
     * Standard swap function (works with ADL).
     */
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
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
    const_reference operator()(size_type i,size_type j) const { if(i == j) return value_type(1); else return value_type(0); };
    
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
    size_type get_col_count() const { return rowCount; };
    
    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_col_count(size_type aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      rowCount = aColCount;
    };
    
    
    /**
     * Negate the matrix.
     * \return This matrix, by constant reference.
     * \test PASSED
     */
    mat<value_type,mat_structure::scalar> operator -() const {
      return mat<value_type,mat_structure::scalar>(rowCount,value_type(-1));
    };
    
    /**
     * Transpose the matrix.
     * \param rhs the matrix to be transposed.
     * \return The rhs matrix, by value.
     * \test PASSED
     */
    friend self transpose(const self& rhs) {
      return rhs;
    };
    
    /**
     * Transpose and move the matrix.
     * \param rhs the matrix to be transposed and moved (emptied).
     * \return The rhs matrix, by value.
     * \test PASSED
     */
    friend self transpose_move(self& rhs) {
      self result; 
      swap(result, rhs);
      return result;
    };
    
    /**
     * Returns the trace of a matrix.
     * \param M A matrix.
     * \return the trace of matrix M.
     * \test PASSED
     */
    friend value_type trace(const self& M) {
      return M.rowCount;
    };
    
    /**
     * Appends a matrix to another.
     * \param lhs the matrix to which 'rhs' will be appended to.
     * \param rhs the matrix to append to 'lhs'.
     * \test PASSED
     */
    friend void append_block_diag(self& lhs, const self& rhs) {
      lhs.rowCount += rhs.rowCount;
    };
        
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(rowCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(rowCount);
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
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > type;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix< mat<T,mat_structure::identity,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
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
 mat<T,mat_structure::scalar,Alignment,Allocator> >::type
  operator *(const mat<T,mat_structure::identity,Alignment,Allocator>& M, const T& S) {
    return mat<T,mat_structure::scalar,Alignment,Allocator>(M.get_row_count(),S);
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
 mat<T,mat_structure::scalar,Alignment,Allocator> >::type 
  operator *(const T& S,const mat<T,mat_structure::identity,Alignment,Allocator>& M) {
    return mat<T,mat_structure::scalar,Alignment,Allocator>(M.get_row_count(),S);
  };

#if 0

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
  
#endif
  



#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

extern template class mat<double, mat_structure::identity>;
extern template class mat<float, mat_structure::identity>;

extern template mat<double,mat_structure::identity> mat_ident<double>(mat<double,mat_structure::identity>::size_type aRowCount);
extern template mat<float,mat_structure::identity> mat_ident<float>(mat<float,mat_structure::identity>::size_type aRowCount);

extern template vect<double,2> operator *<double, vect<double,2>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,2>& V);
extern template vect<double,3> operator *<double, vect<double,3>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,3>& V);
extern template vect<double,4> operator *<double, vect<double,4>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,4>& V);
extern template vect<double,6> operator *<double, vect<double,6>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,6>& V);
extern template vect_n<double> operator *<double, vect_n<double>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M, const vect_n<double>& V);

extern template vect<double,2> operator *<double, vect<double,2>, mat_alignment::column_major, std::allocator<double> >(const vect<double,2>& V,const mat<double,mat_structure::identity>& M);
extern template vect<double,3> operator *<double, vect<double,3>, mat_alignment::column_major, std::allocator<double> >(const vect<double,3>& V,const mat<double,mat_structure::identity>& M);
extern template vect<double,4> operator *<double, vect<double,4>, mat_alignment::column_major, std::allocator<double> >(const vect<double,4>& V,const mat<double,mat_structure::identity>& M);
extern template vect<double,6> operator *<double, vect<double,6>, mat_alignment::column_major, std::allocator<double> >(const vect<double,6>& V,const mat<double,mat_structure::identity>& M);
extern template vect_n<double> operator *<double, vect_n<double>, mat_alignment::column_major, std::allocator<double> >(const vect_n<double>& V,const mat<double,mat_structure::identity>& M);

extern template vect<float,2> operator *<float, vect<float,2>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,2>& V);
extern template vect<float,3> operator *<float, vect<float,3>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,3>& V);
extern template vect<float,4> operator *<float, vect<float,4>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,4>& V);
extern template vect<float,6> operator *<float, vect<float,6>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,6>& V);
extern template vect_n<float> operator *<float, vect_n<float>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M, const vect_n<float>& V);

extern template vect<float,2> operator *<float, vect<float,2>, mat_alignment::column_major, std::allocator<float> >(const vect<float,2>& V,const mat<float,mat_structure::identity>& M);
extern template vect<float,3> operator *<float, vect<float,3>, mat_alignment::column_major, std::allocator<float> >(const vect<float,3>& V,const mat<float,mat_structure::identity>& M);
extern template vect<float,4> operator *<float, vect<float,4>, mat_alignment::column_major, std::allocator<float> >(const vect<float,4>& V,const mat<float,mat_structure::identity>& M);
extern template vect<float,6> operator *<float, vect<float,6>, mat_alignment::column_major, std::allocator<float> >(const vect<float,6>& V,const mat<float,mat_structure::identity>& M);
extern template vect_n<float> operator *<float, vect_n<float>, mat_alignment::column_major, std::allocator<float> >(const vect_n<float>& V,const mat<float,mat_structure::identity>& M);


extern template mat<double,mat_structure::scalar> operator *(const mat<double,mat_structure::identity>& M, const double& S);
extern template mat<double,mat_structure::scalar> operator *(const double& S, const mat<double,mat_structure::identity>& M);

extern template mat<float,mat_structure::scalar> operator *(const mat<float,mat_structure::identity>& M, const float& S);
extern template mat<float,mat_structure::scalar> operator *(const float& S, const mat<float,mat_structure::identity>& M);


#endif






};


#endif










