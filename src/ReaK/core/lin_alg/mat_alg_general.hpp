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

#ifndef REAK_MAT_ALG_GENERAL_HPP
#define REAK_MAT_ALG_GENERAL_HPP

#include "base/defs.hpp"
#include "vect_alg.hpp"
#include "vect_concepts.hpp"
#include "mat_concepts.hpp"
#include "mat_traits.hpp"
#include "stride_iterator.hpp"

#include "base/serializable.hpp"
#include "rtti/so_register_type.hpp"

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>


namespace ReaK {


/**
 * This class is the general class template for all matrix classes in the ReaK linear algebra
 * libraries. The general template itself should never be used and will cause a compilation 
 * error if it is, all useful matrix class templates are, in fact, partial specializations of
 * this general class template.
 * 
 * Models: ReadableMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Structure Enum which defines the structure of the matrix, see mat_structure::tag.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */  
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
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_writable_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_resizable_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct has_allocator_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct mat_product_priority< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<Structure>::value);
  typedef detail::product_priority<Structure> type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct mat_addition_priority< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<Structure>::value);
  typedef detail::addition_priority<Structure> type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_square_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure != mat_structure::rectangular) && (Structure != mat_structure::nil)));
  typedef is_square_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_symmetric_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::symmetric) || (Structure == mat_structure::diagonal) || (Structure == mat_structure::tridiagonal) || (Structure == mat_structure::identity)));
  typedef is_symmetric_matrix< mat<T,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  mat_alignment::tag Alignment,
	  typename Allocator>
struct is_diagonal_matrix< mat<T,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
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
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_writable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_resizable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct has_allocator_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct mat_product_priority< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::product_priority<Structure>::value);
  typedef detail::product_priority<Structure> type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct mat_addition_priority< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = detail::addition_priority<Structure>::value);
  typedef detail::addition_priority<Structure> type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_square_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = (RowCount == ColCount));
  typedef is_square_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_symmetric_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::symmetric) || (Structure == mat_structure::diagonal) || (Structure == mat_structure::tridiagonal) || (Structure == mat_structure::identity)));
  typedef is_symmetric_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
};

template <typename T, 
          mat_structure::tag Structure,
	  unsigned int RowCount, unsigned int ColCount,
	  mat_alignment::tag Alignment>
struct is_diagonal_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == mat_structure::diagonal) || (Structure == mat_structure::identity)));
  typedef is_diagonal_matrix< mat_fix<T,Structure,RowCount,ColCount,Alignment> > type;
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



/*
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< 
  is_fully_writable_matrix<Matrix1>::value &&
  is_resizable_matrix<Matrix1>::value &&
  is_readable_matrix<Matrix2>::value &&
  boost::mpl::or_< boost::mpl::equal_to< boost::mpl::integral_c< mat_structure::tag, mat_traits<Matrix1>::structure>,
                                         boost::mpl::integral_c< mat_structure::tag, mat_structure::square> >,
                   boost::mpl::equal_to< boost::mpl::integral_c< mat_structure::tag, mat_traits<Matrix2>::structure>,
                                         boost::mpl::integral_c< mat_structure::tag, mat_structure::rectangular> > >::value,
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
*/











};


#endif










