/**
 * \file mat_alg_general_hpp
 * 
 * This library implements the general versions of many meta-functions (templates), 
 * functions, and operators. These are meant to be used when no more-specialized 
 * implementations exist for the matrix types involved.
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

#ifndef MAT_ALG_GENERAL_HPP
#define MAT_ALG_GENERAL_HPP

#include "vect_alg.hpp"
#include "vect_concepts.hpp"
#include "mat_concepts.hpp"
#include "mat_traits.hpp"
#include "stride_iterator.hpp"

#include <boost/concept_check.hpp>

#include "base/serializable.hpp"
#include "rtti/so_register_type.hpp"

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>


namespace ReaK {


  
template <typename T, 
          mat_structure::tag Structure = mat_structure::rectangular,
	  mat_alignment::tag Alignment = mat_alignment::column_major,
	  typename Allocator = std::allocator<T> >
class mat { char this_specialization_is_not_available_or_possible[0]; };


template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_readable_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_writable_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_resizable_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct has_allocator_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct mat_product_priority< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<Structure>::value);
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct mat_addition_priority< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<Structure>::value);
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_square_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = ((Structure != mat_structure::rectangular) && (Structure != mat_structure::nil)));
  typedef is_square_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_symmetric_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::symmetric) || (Structure == mat_structure::diagonal) || (Structure == mat_structure::tridiagonal) || (Structure == mat_structure::identity)));
  typedef is_symmetric_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_diagonal_matrix< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::diagonal) || (Structure == mat_structure::identity)));
  typedef is_diagonal_matrix< mat<T,Structure,Alignment,Allocator> > type;
};





template <typename T, 
          mat_structure::tag Structure = mat_structure::rectangular,
	  unsigned int RowCount = 1,
	  unsigned int ColCount = RowCount,
	  mat_alignment::tag Alignment = mat_alignment::column_major >
class mat_fix { char this_specialization_is_not_available_or_possible[0]; };


template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_readable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_writable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_resizable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct has_allocator_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct mat_product_priority< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<Structure>::value);
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct mat_addition_priority< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<Structure>::value);
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_square_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = (RowCount == ColCount));
  typedef is_square_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_symmetric_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::symmetric) || (Structure == mat_structure::diagonal) || (Structure == mat_structure::tridiagonal) || (Structure == mat_structure::identity)));
  typedef is_symmetric_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_diagonal_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::diagonal) || (Structure == mat_structure::identity)));
  typedef is_diagonal_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};







template <typename Matrix1, typename Matrix2>
struct mat_product_result {
  typedef Matrix1 type;
};

template <typename Matrix1, typename Matrix2>
struct mat_addition_result {
  typedef Matrix1 type;
};



template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix2>
struct mat_product_result< mat<T,Structure,Alignment,Allocator>, Matrix2> {
  typedef mat<T,mat_structure::rectangular,Alignment,Allocator> type;
};




/**
 * Prints a matrix to a standard output stream (<<) as "((a11; a12); (a21; a22))". \test PASSED
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, std::ostream& >::type
  operator <<(std::ostream& out_stream,const Matrix& M) {
  out_stream << "(";
  if((M.get_row_count() != 0) && (M.get_col_count() != 0)) {
    for(unsigned int i=0;i<M.get_row_count();++i) {
      out_stream << "(" << M(i,0);
      for(unsigned int j=1;j<M.get_col_count();++j) {
        out_stream << "; " << M(i,j);
      };
      out_stream << ")";
      if(i != M.get_row_count()-1)
        out_stream << "; ";
    };
  };
  return (out_stream << ")");
};



template <mat_structure::tag Structure, mat_alignment::tag Alignment>
struct mat_indexer { };



namespace rtti {

template <typename T,mat_structure::tag Structure, mat_alignment::tag Alignment,typename Allocator>
struct get_type_id< mat<T,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000012);
  static std::string type_name() { return "mat"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const serialization::serializable& save_type;
  typedef serialization::serializable& load_type;
};

template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Tail>
struct get_type_info< mat<T,Structure,Alignment,Allocator>, Tail > {
  typedef detail::type_id< mat<T,Structure,Alignment,Allocator> , typename get_type_info<T,
                                                                           get_type_info< boost::mpl::integral_c<mat_structure::tag,Structure>,
									   get_type_info< boost::mpl::integral_c<mat_alignment::tag,Alignment>, Tail> > >::type > type;
  static std::string type_name() { return get_type_id< mat<T,Structure,Alignment,Allocator> >::type_name() + "<" + get_type_id<T>::type_name() + "," 
                                                                                                                 + get_type_id< boost::mpl::integral_c<mat_structure::tag,Structure> >::type_name() + "," 
														 + get_type_id< boost::mpl::integral_c<mat_alignment::tag,Alignment> >::type_name() + ">" + "," + Tail::type_name(); };
};

template <typename T,
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct get_type_id< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000013);
  static std::string type_name() { return "mat_fix"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const serialization::serializable& save_type;
  typedef serialization::serializable& load_type;
};

template <typename T, mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount, mat_alignment::tag Alignment, typename Tail>
struct get_type_info< mat_fix<T,Structure,RowCount,ColCount,Alignment>, Tail > {
  typedef detail::type_id< mat_fix<T,Structure,RowCount,ColCount,Alignment> , typename get_type_info<T,
                                                                                       get_type_info< boost::mpl::integral_c<mat_structure::tag,Structure>,
									               get_type_info< boost::mpl::integral_c<unsigned int,RowCount>,
									               get_type_info< boost::mpl::integral_c<unsigned int,ColCount>,
									               get_type_info< boost::mpl::integral_c<mat_alignment::tag,Alignment>, Tail> > > > >::type > type;
  static std::string type_name() { return get_type_id< mat_fix<T,Structure,RowCount,ColCount,Alignment> >::type_name() + "<" + get_type_id<T>::type_name() + "," 
                                                                                                                             + get_type_id< boost::mpl::integral_c<mat_structure::tag,Structure> >::type_name() + "," 
                                                                                                                             + get_type_id< boost::mpl::integral_c<unsigned int,RowCount> >::type_name() + "," 
                                                                                                                             + get_type_id< boost::mpl::integral_c<unsigned int,ColCount> >::type_name() + "," 
														             + get_type_id< boost::mpl::integral_c<mat_alignment::tag,Alignment> >::type_name() + ">" + "," + Tail::type_name(); };
};

};






/*******************************************************************************
                         Basic Operators
*******************************************************************************/

/**
 * General (least-specialized) addition operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
 mat< T, mat_structure::rectangular, Alignment, Allocator > >::type
  operator +(const mat<T,Structure,Alignment,Allocator>& M1,const Matrix2& M2) {
    typedef mat< T, mat_structure::rectangular, Alignment, Allocator > result_type; 
    if((M1.get_row_count() != M2.get_row_count()) || (M1.get_col_count() != M2.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    result_type result(M1.get_row_count(),M1.get_col_count(),T(0),M1.get_allocator());
    for(unsigned int j=0;j<M1.get_col_count();++j)
      for(unsigned int i=0;i<M1.get_row_count();++i)
        result(i,j) = M1(i,j) + M2(i,j);
    return result;
  };



/**
 * General substraction operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
 mat< T, mat_structure::rectangular, Alignment, Allocator > >::type
  operator -(const mat<T,Structure,Alignment,Allocator>& M1,const Matrix2& M2) {
    typedef mat< T, mat_structure::rectangular, Alignment, Allocator > result_type;
    if((M1.get_row_count() != M2.get_row_count()) || (M1.get_col_count() != M2.get_col_count()))
      throw std::range_error("Matrix dimension mismatch.");
    result_type result(M1.get_row_count(),M1.get_col_count(),T(0),M1.get_allocator());
    for(unsigned int j=0;j<M1.get_col_count();++j)
      for(unsigned int i=0;i<M1.get_row_count();++i)
        result(i,j) = M1(i,j) - M2(i,j);
    return result;
  };
  


/**
 * General multiplication operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix2>::value &&
                             ((mat_product_priority< Matrix1 >::value < detail::product_priority<mat_structure::diagonal>::value) &&
                              (mat_product_priority<Matrix2>::value < detail::product_priority<mat_structure::diagonal>::value)), 
 typename mat_product_result<Matrix1,Matrix2>::type >::type
  operator *(const Matrix1& M1,const Matrix2& M2) {
    typedef typename mat_product_result<Matrix1,Matrix2>::type result_type;
    if(M1.get_col_count() != M2.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    typedef typename mat_traits<result_type>::value_type ValueType;
    typedef typename mat_traits<result_type>::value_type SizeType;
    result_type result(M1.get_row_count(),M2.get_col_count(),ValueType(0),M1.get_allocator());
    for(SizeType i=0;i<M1.get_row_count();++i) {
      for(SizeType jj=0;jj<M2.get_col_count();++jj) {
        for(SizeType j=0;j<M1.get_col_count();++j) {
          result(i,jj) += M1(i,j) * M2(j,jj);
        };
      };
    };
    return result;
  };

#if 0
/**
 * General multiplication operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix2>::value &&
                             ((mat_product_priority< mat<T,Structure,Alignment,Allocator> >::value < detail::product_priority<mat_structure::diagonal>::value) &&
                              (mat_product_priority<Matrix2>::value < detail::product_priority<mat_structure::diagonal>::value)), 
 mat< T, mat_structure::rectangular, Alignment, Allocator > >::type
  operator *(const Matrix2& M1, const mat<T,Structure,Alignment,Allocator>& M2) {
    typedef mat< T, mat_structure::rectangular, Alignment, Allocator > result_type;
    if(M1.get_col_count() != M2.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    result_type result(M1.get_row_count(),M2.get_col_count(),T(0),M2.get_allocator());
    for(unsigned int i=0;i<M1.get_row_count();++i) {
      for(unsigned int jj=0;jj<M2.get_col_count();++jj) {
        for(unsigned int j=0;j<M1.get_col_count();++j) {
          result(i,jj) += M1(i,j) * M2(j,jj);
        };
      };
    };
    return result;
  };
#endif
/**
 * Matrix multiplication operator with a column vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param V The column vector (second operand).
 * \return Column vector equal to M * V.
 * \throw std::range_error if the matrix column count does not correspond to the vector dimension.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Vector>
typename boost::enable_if_c< is_readable_vector<Vector>::value, 
 vect_n< T, Allocator > >::type
  operator *(const mat<T,Structure,Alignment,Allocator>& M, const Vector& V) {
    typedef vect_n< T, Allocator > result_type;
    if(V.size() != M.get_col_count())
      throw std::range_error("Matrix dimension mismatch.");
    result_type result(M.get_row_count(),T(0),M.get_allocator());
    for(unsigned int i=0;i<M.get_row_count();++i)
      for(unsigned int j=0;j<M.get_col_count();++j)
        result[i] += M(i,j) * V[j];
    return result;
  };
  
/**
 * Matrix multiplication operator with a row vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param V The row vector (first operand).
 * \param M The matrix (second operand).
 * \return Row vector equal to V * M.
 * \throw std::range_error if the matrix row count does not correspond to the vector dimension.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Vector>
typename boost::enable_if_c< is_readable_vector<Vector>::value, 
 vect_n< T, Allocator > >::type
  operator *(const Vector& V,const mat<T,Structure,Alignment,Allocator>& M) {
    typedef vect_n< T, Allocator > result_type;
    if(V.size() != M.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    result_type result(M.get_col_count(),T(0),V.get_allocator());
    for(unsigned int j=0;j<M.get_col_count();++j) {
      result[j] = 0;
      for(unsigned int i=0;i<M.get_row_count();++i) {
        result[j] += M(i,j) * V[i];
      };
    };
    return result;
  };

  
/**
 * Multiplication by a column-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return column-vector equal to M * V.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, unsigned int Size>
vect<T,Size> operator *(const mat<T,Structure,Alignment,Allocator>& M,const vect<T,Size>& V) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  vect<T,Size> result;
  for(std::size_t i=0;i<Size;++i) {
    result[i] = 0;
    for(std::size_t j=0;j<Size;++j)
      result[i] += M(i,j) * V[j];
  };
  return result;
};


/**
 * Multiplication by a row-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return row-vector equal to V * M.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, unsigned int Size>
vect<T,Size> operator *(const vect<T,Size>& V,const mat<T,Structure,Alignment,Allocator>& M) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  vect<T,Size> result;
  for(std::size_t j=0;j<Size;++j) {
    result[j] = 0;
    for(std::size_t i=0;i<Size;++i)
      result[j] += M(i,j) * V[i];
  };
  return result;
};
  

/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param S The scalar (second operand).
 * \return Column-major matrix equal to M * S.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
mat<T,Structure,Alignment,Allocator> operator *(mat<T,Structure,Alignment,Allocator> M, const T& S) {
  M *= S;
  return M;
};

/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param S The scalar (first operand).
 * \param M The matrix (second operand).
 * \return Column-major matrix equal to S * M.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
mat<T,Structure,Alignment,Allocator> operator *(const T& S, mat<T,Structure,Alignment,Allocator> M) {
  M *= S;
  return M;
};





/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && is_readable_matrix<Matrix2>::value, bool >::type
  operator ==(const Matrix1& M1,const Matrix2& M2) {
    if((M1.get_row_count() != M2.get_row_count()) || (M1.get_col_count() != M2.get_col_count()))
      return false;
    for(unsigned int j=0;j<M1.get_col_count();++j)
      for(unsigned int i=0;i<M1.get_row_count();++i)
        if(M1(i,j) != M2(i,j))
          return false;
    return true;
  };

/**
 * Inequality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && is_readable_matrix<Matrix2>::value, bool >::type
  operator !=(const Matrix1& M1,const Matrix2& M2) {
    if((M1.get_row_count() != M2.get_row_count()) || (M1.get_col_count() != M2.get_col_count()))
      return true;
    for(unsigned int j=0;j<M1.get_col_count();++j)
      for(unsigned int i=0;i<M1.get_row_count();++i)
        if(M1(i,j) != M2(i,j))
          return true;
    return false;
  };  





/**
 * Verifies that the matrix A is diagonal up to a tolerance.
 */
template <class Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
 bool >::type is_diagonal(const Matrix& A, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  using std::fabs;
  if(A.get_row_count() != A.get_col_count())
    return false;
  for(typename mat_traits<Matrix>::size_type i=0;i<A.get_row_count();i++)
    for(typename mat_traits<Matrix>::size_type j=0;j<i;j++)
      if((fabs(A(j,i)) > NumTol) || (fabs(A(i,j)) > NumTol))
	return false;
  return true;
};

/**
 * Verifies that the matrix A is symmetric up to a tolerance.
 */
template <class Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
bool >::type is_symmetric(const Matrix& A, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  using std::fabs;
  if(A.get_row_count() != A.get_col_count())
    return false;
  for(typename mat_traits<Matrix>::size_type i=0;i<A.get_row_count();i++)
    for(typename mat_traits<Matrix>::size_type j=0;j<i;j++)
      if(fabs(A(j,i) - A(i,j)) > NumTol)
	return false;
  return true;
};






template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix1>::value &&
                             is_resizable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             (((mat_traits<Matrix1>::structure == mat_structure::square) && (is_square_matrix<Matrix2>::value)) ||
                             (mat_traits<Matrix1>::structure == mat_structure::rectangular)),
void >::type append_block_diag(Matrix1& A, const Matrix2& B) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType oldRowCount = A.get_row_count();
  SizeType oldColCount = A.get_col_count();
  A.set_col_count(oldColCount + B.get_col_count(),true);
  A.set_row_count(oldRowCount + B.get_row_count(),true);
  for(SizeType i = 0; i < B.get_row_count(); ++i)
    for(SizeType j = 0; j < B.get_col_count(); ++j)
      A(i + oldRowCount,j + oldColCount) = B(i,j);
};












};


#endif










