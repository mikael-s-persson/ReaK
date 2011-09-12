/**
 * \file mat_traits.hpp
 * 
 * This header declares a number of standard matrix type traits. These include 
 * matrix alignment tags (column-major and row-major), matrix structure tags 
 * (rectangular, square, symmetric, diagonal, skew-symmetric, etc.), some ReaK::rtti
 * template specializations for these tags, the main mat_traits<> template which
 * defines the nested typedefs related to a matrix class (and required by an implementation
 * of a matrix class), and, finally, a series of meta-functions to compute the product-priority
 * additive-priority of matrix classes based on their structural tags. 
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

#ifndef REAK_MAT_TRAITS_HPP
#define REAK_MAT_TRAITS_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <boost/type_traits.hpp>

#include "rtti/so_type.hpp"
#include "vect_alg.hpp"

namespace ReaK {
  
  
namespace mat_alignment {
  
  /**
   * These tags are used to mark a matrix as having a column-major or row-major alignment
   * in memory (note that array-of-array is not an option, due to the obvious performance problems 
   * with such a memory model for matrix representations).
   */
  enum tag {
    column_major = 1,
    row_major = 2
  };
};

namespace mat_structure {
  
  /**
   * These tags are used to identify a matrix class by its underlying structure. The most 
   * general structure is, of course, the rectangular matrix (any size, any shape, and densely 
   * populated). Then, the frequently used tags are square, symmetric, skew_symmetric, and 
   * diagonal. Additionally, the nil and identity tags signify that a matrix type is constraint
   * to never be anything other than nil or identity (obviously, these types will be read-only
   * and do not require storage other than the size information).
   */
  enum tag {
    rectangular = 1,
    square = 2,
    symmetric = 3,
    skew_symmetric = 4,
    diagonal = 5,
    upper_triangular = 6,
    lower_triangular = 7,
    orthogonal = 8,
    tridiagonal = 9,
    nil = 10,
    identity = 11,
    scalar = 12
  };
};
  
/*
 * The following are ReaK::rtti class templates which are used to associate 
 * type information to the matrix tags (mat_alignment::tag and mat_structure::tag).
 * 
 * Please ignore, unless interested in the inner-workings of the ReaK::rtti system.
 */
namespace rtti {
  
template <>
struct get_type_id< boost::mpl::integral_c<mat_alignment::tag, mat_alignment::column_major > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
  static std::string type_name() { return "column_major"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_alignment::tag, mat_alignment::row_major > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
  static std::string type_name() { return "row_major"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <mat_alignment::tag U, typename Tail>
struct get_type_info< boost::mpl::integral_c<mat_alignment::tag, U >, Tail > {
  typedef detail::type_id<boost::mpl::integral_c<mat_alignment::tag, U >, typename Tail::type> type;
  static std::string type_name() { return get_type_id< boost::mpl::integral_c<mat_alignment::tag, U > >::type_name() + "," + Tail::type_name(); };
};

  
template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::rectangular > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
  static std::string type_name() { return "rectangular"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::square > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
  static std::string type_name() { return "square"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::symmetric > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 3);
  static std::string type_name() { return "symmetric"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::skew_symmetric > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 4);
  static std::string type_name() { return "skew_symmetric"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::diagonal > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 5);
  static std::string type_name() { return "diagonal"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::upper_triangular > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 6);
  static std::string type_name() { return "upper_triangular"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::lower_triangular > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 7);
  static std::string type_name() { return "lower_triangular"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::orthogonal > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 8);
  static std::string type_name() { return "orthogonal"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::tridiagonal > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 9);
  static std::string type_name() { return "tridiagonal"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::nil > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 10);
  static std::string type_name() { return "nil"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::identity > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 11);
  static std::string type_name() { return "identity"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<mat_structure::tag, mat_structure::scalar > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 12);
  static std::string type_name() { return "scalar"; };
  static construct_ptr CreatePtr() { return NULL; };
};


template <mat_structure::tag U, typename Tail>
struct get_type_info< boost::mpl::integral_c<mat_structure::tag, U >, Tail > {
  typedef detail::type_id<boost::mpl::integral_c<mat_structure::tag, U >, typename Tail::type> type;
  static std::string type_name() { return get_type_id< boost::mpl::integral_c<mat_structure::tag, U > >::type_name() + "," + Tail::type_name(); };
};

  
};


/**
 * This type-traits definition describes the nested typedefs that are expected 
 * from an implementation of a matrix class. They are mostly inspired from 
 * STL-containers' traits. The simplest why to implement a new matrix class 
 * that conforms with these mat_traits requirements is to provide all those 
 * nested typedefs and static constants. The other alternative, which is non-intrusive,
 * is to define a specialization for mat_traits for the new matrix class, providing
 * all the public members as seen in this general trait template.
 */
template <typename Matrix>
struct mat_traits {
  /// The type of the elements of the matrix.
  typedef typename Matrix::value_type value_type;
  /// The type for a reference to an element of the matrix
  typedef typename Matrix::reference reference;
  /// The type for a const reference to an element of the matrix.
  typedef typename Matrix::const_reference const_reference;
  /// The type for a pointer to an element of the matrix.
  typedef typename Matrix::pointer pointer;
  /// The type for a const pointer to an element of the matrix.
  typedef typename Matrix::const_pointer const_pointer;
  /// The type of the allocator for the matrix (can be void if the matrix does not have an allocator).
  typedef typename Matrix::allocator_type allocator_type;
  
  /// The type of the row-iterator for the matrix (a row-iterator goes from one row to another (on the same column)).
  typedef typename Matrix::row_iterator row_iterator;
  /// The type of the const row-iterator for the matrix (a row-iterator goes from one row to another (on the same column)).
  typedef typename Matrix::const_row_iterator const_row_iterator;
  /// The type of the column-iterator for the matrix (a column-iterator goes from one column to another (on the same row)).
  typedef typename Matrix::col_iterator col_iterator;
  /// The type of the const column-iterator for the matrix (a column-iterator goes from one column to another (on the same row)).
  typedef typename Matrix::const_col_iterator const_col_iterator;
  
  /// The type of the size descriptors (or index descriptors) of the matrix class.
  typedef typename Matrix::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the matrix class.
  typedef typename Matrix::difference_type difference_type;
  
  /// The static row count. Should be 0 if the row count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_row_count = Matrix::static_row_count);
  /// The static column count. Should be 0 if the column count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_col_count = Matrix::static_col_count);
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = Matrix::alignment);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(mat_structure::tag, structure = Matrix::structure);
  
};



/*
 * This detail section includes several template specializations which are used 
 * to define the priority of certain matrix structures when it comes to doing 
 * additions (or subtractions) and multiplications. These meta-functions are 
 * used to perform Sfinae-based switching between overloaded versions of operators 
 * and functions based on the types. For example, when multiplying a symmetric 
 * matrix with a diagonal matrix, surely, the general matrix multiplication operator
 * should be turned off by sfinae, additionally, the diagonal matrix multiplication 
 * operator should be selected over the symmetric matrix multiplication operator. These
 * meta-functions allow this kind of situations to be detected and the correct overload 
 * to be activated (by turning off the others (via Sfinae) which have lower priority).
 *
 * Please ignore this section unless interested in the inner-workings of Sfinae-based 
 * switching of operator overload selection for the matrix classes. 
 */
namespace detail {
  
  template <mat_structure::tag Structure>
  struct product_priority {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 0);
    typedef product_priority<Structure> type;
  };
  
  template <>
  struct product_priority<mat_structure::rectangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 1);
    typedef product_priority<mat_structure::rectangular> type;
  };
  
  template <>
  struct product_priority<mat_structure::square> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 2);
    typedef product_priority<mat_structure::square> type;
  };
  
  template <>
  struct product_priority<mat_structure::symmetric> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 10);
    typedef product_priority<mat_structure::symmetric> type;
  };
  
  template <>
  struct product_priority<mat_structure::skew_symmetric> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 11);
    typedef product_priority<mat_structure::skew_symmetric> type;
  };
  
  template <>
  struct product_priority<mat_structure::diagonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 40);
    typedef product_priority<mat_structure::diagonal> type;
  };
  
  template <>
  struct product_priority<mat_structure::scalar> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 41);
    typedef product_priority<mat_structure::scalar> type;
  };
  
  template <>
  struct product_priority<mat_structure::upper_triangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 20);
    typedef product_priority<mat_structure::upper_triangular> type;
  };
  
  template <>
  struct product_priority<mat_structure::lower_triangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 21);
    typedef product_priority<mat_structure::lower_triangular> type;
  };
  
  template <>
  struct product_priority<mat_structure::orthogonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 3);
    typedef product_priority<mat_structure::orthogonal> type;
  };
  
  template <>
  struct product_priority<mat_structure::tridiagonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 30);
    typedef product_priority<mat_structure::tridiagonal> type;
  };
  
  template <>
  struct product_priority<mat_structure::nil> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 50);
    typedef product_priority<mat_structure::nil> type;
  };
  
  template <>
  struct product_priority<mat_structure::identity> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 49);
    typedef product_priority<mat_structure::identity> type;
  };
  
  
  
  
  template <mat_structure::tag Structure>
  struct addition_priority {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 0);
    typedef addition_priority<Structure> type;
  };
  
  template <>
  struct addition_priority<mat_structure::rectangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 1);
    typedef addition_priority<mat_structure::rectangular> type;
  };
  
  template <>
  struct addition_priority<mat_structure::square> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 2);
    typedef addition_priority<mat_structure::square> type;
  };
  
  template <>
  struct addition_priority<mat_structure::symmetric> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 10);
    typedef addition_priority<mat_structure::symmetric> type;
  };
  
  template <>
  struct addition_priority<mat_structure::skew_symmetric> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 11);
    typedef addition_priority<mat_structure::skew_symmetric> type;
  };
  
  template <>
  struct addition_priority<mat_structure::diagonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 40);
    typedef addition_priority<mat_structure::diagonal> type;
  };
  
  template <>
  struct addition_priority<mat_structure::scalar> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 40);
    typedef addition_priority<mat_structure::scalar> type;
  };
  
  template <>
  struct addition_priority<mat_structure::identity> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 40);
    typedef addition_priority<mat_structure::identity> type;
  };
  
  template <>
  struct addition_priority<mat_structure::upper_triangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 20);
    typedef addition_priority<mat_structure::upper_triangular> type;
  };
  
  template <>
  struct addition_priority<mat_structure::lower_triangular> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 21);
    typedef addition_priority<mat_structure::lower_triangular> type;
  };
  
  template <>
  struct addition_priority<mat_structure::orthogonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 3);
    typedef addition_priority<mat_structure::orthogonal> type;
  };
  
  template <>
  struct addition_priority<mat_structure::tridiagonal> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 30);
    typedef addition_priority<mat_structure::tridiagonal> type;
  };
  
  template <>
  struct addition_priority<mat_structure::nil> {
    typedef boost::mpl::integral_c_tag tag;
    typedef std::size_t value_type;
    BOOST_STATIC_CONSTANT(std::size_t, value = 50);
    typedef addition_priority<mat_structure::nil> type;
  };
  
};


template <typename Matrix>
struct mat_product_priority {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = 0);
  typedef mat_product_priority<Matrix> type;
};

template <typename Matrix>
struct mat_addition_priority {
  typedef boost::mpl::integral_c_tag tag;
  typedef std::size_t value_type;
  BOOST_STATIC_CONSTANT(std::size_t, value = 0);
  typedef mat_addition_priority<Matrix> type;
};






};




#endif









