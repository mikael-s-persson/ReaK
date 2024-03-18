/**
 * \file mat_composite_adaptor.h
 *
 * This library provides a number of adaptor classes that can be used to concatenate matrices together.
 * It is sometimes desirable to concatenate matrices, for various reasons. You might have a few matrices
 * stored separately but that actually are parts of an imaginary block-structured super-matrix, if most
 * of the time these matrices are manipulated separately but once in a while an operation on the super-matrix
 * is required, than the classes included in this file can come in handy. Another use is to be able to make
 * a block-structured matrix where different blocks are better represented (or stored) by different types
 * of matrices (e.g. some symmetric blocks, some nil-blocks, or some identity blocks, etc.). The matrix
 * composition classes provided by this library allow you to have heterogeneous block-structured matrices in
 * addition to allowing you to obtain const and non-const views over a collection of matrices as a super-matrix.
 *
 * Additionally, this library is best used with C++0x compatibility because it can use overloading based
 * on rvalue-references which will increase the ease of use. Also, this library provides operator overloads
 * to allow for the concatenation of the matrices to have a very neat syntax (and these operators operate
 * better under C++0x support for rvalue-references), as so:
 *  A = ( A_11 & A_12 |
 *        A_21 & A_22 );
 *
 * Also note that all concatenations are based on concatenating two matrices horizontally or vertically. So,
 * for the concatenation of several matrices, you have to use several nested concatenations (or recursive).
 * Fortunately, factory function templates and operator overloads hide all that nasty syntax away. Moreover,
 * C++0x's feature of type inference (like 'auto' keyword and decltype()) also greatly simplify the syntax
 * on the caller's side.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_COMPOSITE_ADAPTOR_H_
#define REAK_MATH_LIN_ALG_MAT_COMPOSITE_ADAPTOR_H_

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/math/lin_alg/mat_views.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_views.h"

#include <type_traits>

namespace ReaK {

/// This class template forms the horizontal concatenation of two matrices, which it stores by value.
/// This class makes the concatenation of the two matrices look as if it was just one matrix (and so,
/// this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both LeftMatrix and RightMatrix,
/// except for ResizableMatrix.
///
/// \tparam LeftMatrix Matrix type for the left matrix.
/// \tparam RightMatrix Matrix type for the right matrix.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
class mat_horiz_cat {
 public:
  using self = mat_horiz_cat<LeftMatrix, RightMatrix>;

  using value_type = mat_value_type_t<LeftMatrix>;

  using reference = typename mat_traits<LeftMatrix>::reference;
  using const_reference = typename mat_traits<LeftMatrix>::const_reference;
  using pointer = typename mat_traits<LeftMatrix>::pointer;
  using const_pointer = typename mat_traits<LeftMatrix>::const_pointer;

  using col_iterator = typename mat_traits<LeftMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<LeftMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<LeftMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<LeftMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<LeftMatrix>;
  using difference_type = typename mat_traits<LeftMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfExpectedEqual(mat_traits<LeftMatrix>::static_row_count,
                                   mat_traits<RightMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfConcat(mat_traits<LeftMatrix>::static_col_count,
                            mat_traits<RightMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<LeftMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  LeftMatrix ml;   ///< Holds the left part of the matrix.
  RightMatrix mr;  ///< Holds the right part of the matrix.
 public:
  /// Default constructor.
  mat_horiz_cat() : ml(), mr() {}

  /// Parametrized constructor.
  /// \param aML Matrix to fill the left part of the matrix.
  /// \param aMR Matrix to fill the right part of the matrix.
  mat_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR)
      : ml(aML), mr(aMR) {
    if (ml.get_row_count() != mr.get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Standard copy-constructor.
  /// \param aObj Right-hand-side of the copy.
  mat_horiz_cat(const self& aObj) = default;

  /// Parametrized move-constructor.
  /// \param aML Matrix to fill and be moved into the left part of the matrix.
  /// \param aMR Matrix to fill and be moved into the right part of the matrix.
  mat_horiz_cat(LeftMatrix&& aML, RightMatrix&& aMR)
      : ml(std::move(aML)), mr(std::move(aMR)) {
    if (ml.get_row_count() != mr.get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Standard move-constructor.
  /// \param aObj Right-hand-side of the move.
  mat_horiz_cat(self&& aObj) noexcept = default;

  /// Standard swap function.
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.ml, rhs.ml);
    swap(lhs.mr, rhs.mr);
  }

  /// Standard assignment operator.
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /// Templated assignment operator to assign the content of the matrix with the content
  /// of a matrix of another type (Matrix)
  /// \param rhs Right-hand-side of the assignment.
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != ml.get_row_count()) ||
        (rhs.get_col_count() != ml.get_col_count() + mr.get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    ml = sub(rhs)(range(0, ml.get_row_count()), range(0, ml.get_col_count()));
    mr = sub(rhs)(
        range(0, mr.get_row_count()),
        range(ml.get_col_count(), ml.get_col_count() + mr.get_col_count()));
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) {
    if (j < ml.get_col_count()) {
      return ml(i, j);
    }
    return mr(i, j - ml.get_col_count());
  }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (j < ml.get_col_count()) {
      return ml(i, j);
    }
    return mr(i, j - ml.get_col_count());
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept { return ml.get_row_count(); }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept {
    return ml.get_col_count() + mr.get_col_count();
  }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(ml.get_row_count(),
                          ml.get_col_count() + mr.get_col_count());
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_row_count() != ml.get_row_count()) ||
        (M.get_col_count() != ml.get_col_count() + mr.get_col_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    ml += sub(M)(range(0, ml.get_row_count()), range(0, ml.get_col_count()));
    mr += sub(M)(
        range(0, mr.get_row_count()),
        range(ml.get_col_count(), ml.get_col_count() + mr.get_col_count()));
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_row_count() != ml.get_row_count()) ||
        (M.get_col_count() != ml.get_col_count() + mr.get_col_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    ml -= sub(M)(range(0, ml.get_row_count()), range(0, ml.get_col_count()));
    mr -= sub(M)(
        range(0, mr.get_row_count()),
        range(ml.get_col_count(), ml.get_col_count() + mr.get_col_count()));
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  self& operator*=(const value_type& S) {
    ml *= S;
    mr *= S;
    return *this;
  }

  /// General Matrix multiplication.
  template <ReadableMatrix Matrix>
  self& operator*=(const Matrix& M) {
    if ((M.get_col_count() != get_col_count()) ||
        (M.get_row_count() != get_col_count())) {
      throw std::range_error("Matrix Dimension Mismatch.");
    }
    *this = *this * M;
    return *this;
  }
};

template <typename LeftMatrix, typename RightMatrix>
static constexpr bool is_fully_writable_matrix_v<mat_horiz_cat<LeftMatrix, RightMatrix>> = is_fully_writable_matrix_v<LeftMatrix> &&
                                is_fully_writable_matrix_v<RightMatrix>;

/// This class template forms the horizontal concatenation of two matrices, which it takes by reference
/// (and stores by pointer, to be copyable). This class makes the concatenation of the two matrices
/// look as if it was just one matrix (and so, this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both LeftMatrix and RightMatrix,
/// except for ResizableMatrix.
///
/// \tparam LeftMatrix Matrix type for the left matrix.
/// \tparam RightMatrix Matrix type for the right matrix.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
class mat_ref_horiz_cat {
 public:
  using self = mat_ref_horiz_cat<LeftMatrix, RightMatrix>;

  using value_type = mat_value_type_t<LeftMatrix>;

  using reference = typename mat_traits<LeftMatrix>::reference;
  using const_reference = typename mat_traits<LeftMatrix>::const_reference;
  using pointer = typename mat_traits<LeftMatrix>::pointer;
  using const_pointer = typename mat_traits<LeftMatrix>::const_pointer;

  using col_iterator = typename mat_traits<LeftMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<LeftMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<LeftMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<LeftMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<LeftMatrix>;
  using difference_type = typename mat_traits<LeftMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfExpectedEqual(mat_traits<LeftMatrix>::static_row_count,
                                   mat_traits<RightMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfConcat(mat_traits<LeftMatrix>::static_col_count,
                            mat_traits<RightMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<LeftMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  LeftMatrix* ml;   ///< Holds the reference to the left part of the matrix.
  RightMatrix* mr;  ///< Holds the reference to the right part of the matrix.
 public:
  /// Default constructor.
  mat_ref_horiz_cat() : ml(), mr() {}

  /// Parametrized constructor.
  /// \param aML Matrix to become the left part of the matrix.
  /// \param aMR Matrix to become the right part of the matrix.
  mat_ref_horiz_cat(LeftMatrix& aML, RightMatrix& aMR) : ml(&aML), mr(&aMR) {
    if (ml->get_row_count() != mr->get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Copy-constructor (shallow-copy).
  mat_ref_horiz_cat(const self& aObj) = default;

  /// Move-constructor (shallow-move).
  mat_ref_horiz_cat(self&& aObj) noexcept = default;

  /// Standard swap function (shallow).
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.ml, rhs.ml);
    swap(lhs.mr, rhs.mr);
  }

  /// Templated assignment operator to assign the content of the matrix with the content
  /// of a matrix of another type (Matrix)
  /// \param rhs Right-hand-side of the assignment.
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != ml->get_row_count()) ||
        (rhs.get_col_count() != ml->get_col_count() + mr->get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    (*ml) =
        sub(rhs)(range(0, ml->get_row_count()), range(0, ml->get_col_count()));
    (*mr) = sub(rhs)(
        range(0, mr->get_row_count()),
        range(ml->get_col_count(), ml->get_col_count() + mr->get_col_count()));
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) {
    if (j < ml->get_col_count()) {
      return (*ml)(i, j);
    }
    return (*mr)(i, j - ml->get_col_count());
  }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (j < ml->get_col_count()) {
      return (*ml)(i, j);
    }
    return (*mr)(i, j - ml->get_col_count());
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept { return ml->get_row_count(); }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept {
    return ml->get_col_count() + mr->get_col_count();
  }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(ml->get_row_count(),
                          ml->get_col_count() + mr->get_col_count());
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_row_count() != ml->get_row_count()) ||
        (M.get_col_count() != ml->get_col_count() + mr->get_col_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    (*ml) +=
        sub(M)(range(0, ml->get_row_count()), range(0, ml->get_col_count()));
    (*mr) += sub(M)(
        range(0, mr->get_row_count()),
        range(ml->get_col_count(), ml->get_col_count() + mr->get_col_count()));
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_row_count() != ml->get_row_count()) ||
        (M.get_col_count() != ml->get_col_count() + mr->get_col_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    (*ml) -=
        sub(M)(range(0, ml->get_row_count()), range(0, ml->get_col_count()));
    (*mr) -= sub(M)(
        range(0, mr->get_row_count()),
        range(ml->get_col_count(), ml->get_col_count() + mr->get_col_count()));
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  self& operator*=(const value_type& S) {
    (*ml) *= S;
    (*mr) *= S;
    return *this;
  }

  /// General Matrix multiplication.
  template <ReadableMatrix Matrix>
  self& operator*=(const Matrix& M) {
    if ((M.get_col_count() != get_col_count()) ||
        (M.get_row_count() != get_col_count())) {
      throw std::range_error("Matrix Dimension Mismatch.");
    }
    *this = *this * M;
    return *this;
  }
};

template <typename LeftMatrix, typename RightMatrix>
static constexpr bool is_fully_writable_matrix_v<mat_ref_horiz_cat<LeftMatrix, RightMatrix>> = is_fully_writable_matrix_v<LeftMatrix> &&
                                is_fully_writable_matrix_v<RightMatrix>;

/// This class template forms the horizontal concatenation of two matrices, which it takes by const-reference
/// (and stores by const-pointer, to be copyable). This class makes the concatenation of the two matrices
/// look as if it was just one matrix (and so, this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both LeftMatrix and RightMatrix,
/// except for ResizableMatrix.
///
/// \tparam LeftMatrix Matrix type for the left matrix.
/// \tparam RightMatrix Matrix type for the right matrix.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
class mat_const_ref_horiz_cat {
 public:
  using self = mat_const_ref_horiz_cat<LeftMatrix, RightMatrix>;

  using value_type = mat_value_type_t<LeftMatrix>;

  using reference = typename mat_traits<LeftMatrix>::reference;
  using const_reference = typename mat_traits<LeftMatrix>::const_reference;
  using pointer = typename mat_traits<LeftMatrix>::pointer;
  using const_pointer = typename mat_traits<LeftMatrix>::const_pointer;

  using col_iterator = typename mat_traits<LeftMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<LeftMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<LeftMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<LeftMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<LeftMatrix>;
  using difference_type = typename mat_traits<LeftMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfExpectedEqual(mat_traits<LeftMatrix>::static_row_count,
                                   mat_traits<RightMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfConcat(mat_traits<LeftMatrix>::static_col_count,
                            mat_traits<RightMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<LeftMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  const LeftMatrix* ml;
  const RightMatrix* mr;

  self& operator=(const self&);
  mat_const_ref_horiz_cat(LeftMatrix&&, RightMatrix&&);

 public:
  /// Parametrized constructor.
  /// \param aML Matrix to fill the left part of the matrix.
  /// \param aMR Matrix to fill the right part of the matrix.
  mat_const_ref_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR)
      : ml(&aML), mr(&aMR) {
    if (ml->get_row_count() != mr->get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Standard copy-constructor (shallow).
  /// \param aObj Right-hand-side of the copy.
  mat_const_ref_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) {}

  /// Standard move-constructor.
  /// \param aObj Right-hand-side of the move.
  mat_const_ref_horiz_cat(self&& aObj) noexcept : ml(aObj.ml), mr(aObj.mr) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (j < ml->get_col_count()) {
      return (*ml)(i, j);
    }
    return (*mr)(i, j - ml->get_col_count());
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept { return ml->get_row_count(); }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept {
    return ml->get_col_count() + mr->get_col_count();
  }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(ml->get_row_count(),
                          ml->get_col_count() + mr->get_col_count());
  }
};

/// This function template will horizontally concatenate two matrices, by copying them into a
/// composite matrix.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices, by copy.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat_copy(const LeftMatrix& ML, const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, RightMatrix>(ML, MR);
}

/// This function template will horizontally concatenate two non-const matrices, by reference to them in a
/// composite matrix.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices, by reference.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(LeftMatrix& ML, RightMatrix& MR) {
  return mat_ref_horiz_cat<LeftMatrix, RightMatrix>(ML, MR);
}

/// This function template will horizontally concatenate two const matrices, by reference to them in a
/// composite matrix.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices, by const-reference.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(const LeftMatrix& ML, const RightMatrix& MR) {
  return mat_const_ref_horiz_cat<LeftMatrix, RightMatrix>(ML, MR);
}

/// This function template will horizontally concatenate two rvalue matrices, by moving them into a
/// composite matrix. This is an overload that will be selected when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices, by copy.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(LeftMatrix&& ML, RightMatrix&& MR) {
  return mat_horiz_cat<LeftMatrix, RightMatrix>(std::move(ML), std::move(MR));
}

/// This function template will horizontally concatenate one lvalue matrix and one rvalue matrix, by referring
/// to the former and moving the latter into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(LeftMatrix& ML, RightMatrix&& MR) {
  return mat_horiz_cat<mat_sub_block<LeftMatrix>, RightMatrix>(
      mat_sub_block<LeftMatrix>(ML), std::move(MR));
}

/// This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(LeftMatrix&& ML, RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix>>(
      std::move(ML), mat_sub_block<RightMatrix>(MR));
}

/// This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(const LeftMatrix& ML, RightMatrix&& MR) {
  return mat_horiz_cat<mat_const_sub_block<LeftMatrix>, RightMatrix>(
      mat_const_sub_block<LeftMatrix>(ML), std::move(MR));
}

/// This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto hcat(LeftMatrix&& ML, const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix>>(
      std::move(ML), mat_const_sub_block<RightMatrix>(MR));
}

/// This operator overload template will horizontally concatenate two rvalue matrices, by moving them into a
/// composite matrix. This is an overload that will be selected when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices, by copy.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(LeftMatrix&& ML, RightMatrix&& MR) {
  return mat_horiz_cat<LeftMatrix, RightMatrix>(std::move(ML), std::move(MR));
}

/// This operator overload template will horizontally concatenate one lvalue matrix and one rvalue matrix, by referring
/// to the former and moving the latter into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(LeftMatrix& ML, RightMatrix&& MR) {
  return mat_horiz_cat<mat_sub_block<LeftMatrix>, RightMatrix>(
      mat_sub_block<LeftMatrix>(ML), std::move(MR));
}

/// This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(LeftMatrix&& ML, RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix>>(
      std::move(ML), mat_sub_block<RightMatrix>(MR));
}

/// This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The value of the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(const LeftMatrix& ML, RightMatrix&& MR) {
  return mat_horiz_cat<mat_const_sub_block<LeftMatrix>, RightMatrix>(
      mat_const_sub_block<LeftMatrix>(ML), std::move(MR));
}

/// This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The value of the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(LeftMatrix&& ML, const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix>>(
      std::move(ML), mat_const_sub_block<RightMatrix>(MR));
}

/// This operator overload template will horizontally concatenate two lvalue matrices, by referring
/// to them in a composite matrix. This is an overload that will be selected
/// when given non-const lvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(LeftMatrix& ML, RightMatrix& MR) {
  return mat_ref_horiz_cat<LeftMatrix, RightMatrix>(ML, MR);
}

/// This operator overload template will horizontally concatenate two lvalue const matrices, by referring
/// to them in a composite matrix. This is an overload that will be selected
/// when given const lvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam LeftMatrix Matrix type for the left part of the composite matrix.
/// \tparam RightMatrix Matrix type for the right part of the composite matrix.
/// \param MU The matrix storing the left part of the composite matrix.
/// \param ML The matrix storing the right part of the composite matrix.
/// \return The composite matrix that horizontally concatenates the two given matrices.
template <ReadableMatrix LeftMatrix, ReadableMatrix RightMatrix>
auto operator&(const LeftMatrix& ML, const RightMatrix& MR) {
  return mat_const_ref_horiz_cat<LeftMatrix, RightMatrix>(ML, MR);
}

/// This class template forms the vertical concatenation of two matrices, which it stores by value.
/// This class makes the concatenation of the two matrices look as if it was just one matrix (and so,
/// this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both LeftMatrix and RightMatrix,
/// except for ResizableMatrix.
///
/// \tparam UpperMatrix Matrix type for the upper matrix.
/// \tparam LowerMatrix Matrix type for the lower matrix.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
class mat_vert_cat {
 public:
  using self = mat_vert_cat<UpperMatrix, LowerMatrix>;

  using value_type = mat_value_type_t<UpperMatrix>;

  using reference = typename mat_traits<UpperMatrix>::reference;
  using const_reference = typename mat_traits<UpperMatrix>::const_reference;
  using pointer = typename mat_traits<UpperMatrix>::pointer;
  using const_pointer = typename mat_traits<UpperMatrix>::const_pointer;

  using col_iterator = typename mat_traits<UpperMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<UpperMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<UpperMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<UpperMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<UpperMatrix>;
  using difference_type = typename mat_traits<UpperMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfConcat(mat_traits<UpperMatrix>::static_row_count,
                            mat_traits<LowerMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfExpectedEqual(mat_traits<UpperMatrix>::static_col_count,
                                   mat_traits<LowerMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<UpperMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  UpperMatrix mu;  ///< Holds the upper part of the matrix.
  LowerMatrix ml;  ///< Holds the lower part of the matrix.
 public:
  /// Default constructor.
  mat_vert_cat() : mu(), ml() {}

  /// Parametrized constructor.
  /// \param aMU Matrix to fill the upper part of the matrix.
  /// \param aML Matrix to fill the lower part of the matrix.
  mat_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML)
      : mu(aMU), ml(aML) {
    if (ml.get_col_count() != mu.get_col_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Copy-constructor.
  mat_vert_cat(const self& aObj) = default;

  /// Parametrized Move-constructor.
  /// \param aMU Matrix to fill the upper part of the matrix.
  /// \param aML Matrix to fill the lower part of the matrix.
  mat_vert_cat(UpperMatrix&& aMU, LowerMatrix&& aML)
      : mu(std::move(aMU)), ml(std::move(aML)) {
    if (ml.get_col_count() != mu.get_col_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Move-constructor.
  mat_vert_cat(self&& aObj) noexcept = default;

  /// Standard swap function.
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.mu, rhs.mu);
    swap(lhs.ml, rhs.ml);
  }

  /// Standard assignment operator.
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /// Templated assignment operator to assign the content of the matrix with the content
  /// of a matrix of another type (Matrix)
  /// \param rhs Right-hand-side of the assignment.
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != mu.get_row_count() + ml.get_row_count()) ||
        (rhs.get_col_count() != mu.get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    mu = sub(rhs)(range(0, mu.get_row_count()), range(0, mu.get_col_count()));
    ml = sub(rhs)(
        range(mu.get_row_count(), mu.get_row_count() + ml.get_row_count()),
        range(0, ml.get_col_count()));
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) {
    if (i < mu.get_row_count()) {
      return mu(i, j);
    }
    return ml(i - mu.get_row_count(), j);
  }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (i < mu.get_row_count()) {
      return mu(i, j);
    }
    return ml(i - mu.get_row_count(), j);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept {
    return mu.get_row_count() + ml.get_row_count();
  }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept { return mu.get_col_count(); }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(mu.get_row_count() + ml.get_row_count(),
                          mu.get_col_count());
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_row_count() != mu.get_row_count() + ml.get_row_count()) ||
        (M.get_col_count() != mu.get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    mu += sub(M)(range(0, mu.get_row_count()), range(0, mu.get_col_count()));
    ml += sub(M)(
        range(mu.get_row_count(), mu.get_row_count() + ml.get_row_count()),
        range(0, ml.get_col_count()));
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_row_count() != mu.get_row_count() + ml.get_row_count()) ||
        (M.get_col_count() != mu.get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    mu -= sub(M)(range(0, mu.get_row_count()), range(0, mu.get_col_count()));
    ml -= sub(M)(
        range(mu.get_row_count(), mu.get_row_count() + ml.get_row_count()),
        range(0, ml.get_col_count()));
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  self& operator*=(const value_type& S) {
    mu *= S;
    ml *= S;
    return *this;
  }

  /// General Matrix multiplication.
  template <ReadableMatrix Matrix>
  self& operator*=(const Matrix& M) {
    if constexpr (!ReadableMatrix<Matrix>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != get_col_count()) ||
          (M.get_row_count() != get_col_count())) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }
};

template <typename UpperMatrix, typename LowerMatrix>
static constexpr bool is_fully_writable_matrix_v<mat_vert_cat<UpperMatrix, LowerMatrix>> = is_fully_writable_matrix_v<UpperMatrix> &&
                                is_fully_writable_matrix_v<LowerMatrix>;

/// This class template forms the vertical concatenation of two matrices, which it takes by reference
/// (and stores by pointer, to be copyable). This class makes the concatenation of the two matrices
/// look as if it was just one matrix (and so, this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both UpperMatrix and LowerMatrix,
/// except for ResizableMatrix.
///
/// \tparam UpperMatrix Matrix type for the left matrix.
/// \tparam LowerMatrix Matrix type for the right matrix.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
class mat_ref_vert_cat {
 public:
  using self = mat_ref_vert_cat<UpperMatrix, LowerMatrix>;

  using value_type = mat_value_type_t<UpperMatrix>;

  using reference = typename mat_traits<UpperMatrix>::reference;
  using const_reference = typename mat_traits<UpperMatrix>::const_reference;
  using pointer = typename mat_traits<UpperMatrix>::pointer;
  using const_pointer = typename mat_traits<UpperMatrix>::const_pointer;

  using col_iterator = typename mat_traits<UpperMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<UpperMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<UpperMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<UpperMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<UpperMatrix>;
  using difference_type = typename mat_traits<UpperMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfConcat(mat_traits<UpperMatrix>::static_row_count,
                            mat_traits<LowerMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfExpectedEqual(mat_traits<UpperMatrix>::static_col_count,
                                   mat_traits<LowerMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<UpperMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  UpperMatrix* mu;  ///< Refers to the upper part of the matrix.
  LowerMatrix* ml;  ///< Refers to the lower part of the matrix.
 public:
  /// Parametrized constructor.
  /// \param aMU Matrix to become the upper part of the matrix.
  /// \param aML Matrix to become the lower part of the matrix.
  mat_ref_vert_cat(UpperMatrix& aMU, LowerMatrix& aML) : mu(&aMU), ml(&aML) {
    if (ml->get_col_count() != mu->get_col_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Standard copy-constructor (shallow-copy).
  mat_ref_vert_cat(const self& aObj) = default;

  /// Standard move-constructor (shallow-move).
  mat_ref_vert_cat(self&& aObj) noexcept = default;

  /// Templated assignment operator to assign the content of the matrix with the content
  /// of a matrix of another type (Matrix)
  /// \param rhs Right-hand-side of the assignment.
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != mu->get_row_count() + ml->get_row_count()) ||
        (rhs.get_col_count() != mu->get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    *mu =
        sub(rhs)(range(0, mu->get_row_count()), range(0, mu->get_col_count()));
    *ml = sub(rhs)(
        range(mu->get_row_count(), mu->get_row_count() + ml->get_row_count()),
        range(0, ml->get_col_count()));
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) {
    if (i < mu->get_row_count()) {
      return (*mu)(i, j);
    }
    return (*ml)(i - mu->get_row_count(), j);
  }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (i < mu->get_row_count()) {
      return (*mu)(i, j);
    }
    return (*ml)(i - mu->get_row_count(), j);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept {
    return mu->get_row_count() + ml->get_row_count();
  }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept { return mu->get_col_count(); }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(mu->get_row_count() + ml->get_row_count(),
                          mu->get_col_count());
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_row_count() != mu->get_row_count() + ml->get_row_count()) ||
        (M.get_col_count() != mu->get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    *mu += sub(M)(range(0, mu->get_row_count()), range(0, mu->get_col_count()));
    *ml += sub(M)(
        range(mu->get_row_count(), mu->get_row_count() + ml->get_row_count()),
        range(0, ml->get_col_count()));
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_row_count() != mu->get_row_count() + ml->get_row_count()) ||
        (M.get_col_count() != mu->get_col_count())) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    *mu -= sub(M)(range(0, mu->get_row_count()), range(0, mu->get_col_count()));
    *ml -= sub(M)(
        range(mu->get_row_count(), mu->get_row_count() + ml->get_row_count()),
        range(0, ml->get_col_count()));
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  self& operator*=(const value_type& S) {
    mu *= S;
    ml *= S;
    return *this;
  }

  /// General Matrix multiplication.
  template <typename Matrix>
  self& operator*=(const Matrix& M) {
    if constexpr (!ReadableMatrix<Matrix>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != get_col_count()) ||
          (M.get_row_count() != get_col_count())) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }
};

template <typename UpperMatrix, typename LowerMatrix>
static constexpr bool is_fully_writable_matrix_v<mat_ref_vert_cat<UpperMatrix, LowerMatrix>> = is_fully_writable_matrix_v<UpperMatrix> &&
                                is_fully_writable_matrix_v<LowerMatrix>;

/// This class template forms the vertical concatenation of two matrices, which it takes by const-reference
/// (and stores by const-pointer, to be copyable). This class makes the concatenation of the two matrices
/// look as if it was just one matrix (and so, this class is an adaptor).
///
/// Models: ReadableMatrix and all matrix concepts modeled by both UpperMatrix and LowerMatrix,
/// except for ResizableMatrix, WritableMatrix, and FullyWritableMatrix.
///
/// \tparam UpperMatrix Matrix type for the left matrix.
/// \tparam LowerMatrix Matrix type for the right matrix.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
class mat_const_ref_vert_cat {
 public:
  using self = mat_const_ref_vert_cat<UpperMatrix, LowerMatrix>;

  using value_type = mat_value_type_t<UpperMatrix>;

  using reference = typename mat_traits<UpperMatrix>::reference;
  using const_reference = typename mat_traits<UpperMatrix>::const_reference;
  using pointer = typename mat_traits<UpperMatrix>::pointer;
  using const_pointer = typename mat_traits<UpperMatrix>::const_pointer;

  using col_iterator = typename mat_traits<UpperMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<UpperMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<UpperMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<UpperMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<UpperMatrix>;
  using difference_type = typename mat_traits<UpperMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfConcat(mat_traits<UpperMatrix>::static_row_count,
                            mat_traits<LowerMatrix>::static_row_count);
  static constexpr unsigned int static_col_count =
      MatStaticSizeIfExpectedEqual(mat_traits<UpperMatrix>::static_col_count,
                                   mat_traits<LowerMatrix>::static_col_count);
  static constexpr mat_alignment::tag alignment =
      mat_traits<UpperMatrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  const UpperMatrix* mu;
  const LowerMatrix* ml;

 public:
  /// Parametrized constructor.
  /// \param aMU The matrix which will become the upper part of this matrix.
  /// \param aML The matrix which will become the lower part of this matrix.
  mat_const_ref_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML)
      : mu(&aMU), ml(&aML) {
    if (ml->get_col_count() != mu->get_col_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
  }

  /// Standard copy-constructor (shallow-copy).
  mat_const_ref_vert_cat(const self& aObj) = default;

  /// Standard move-constructor (shallow-move).
  mat_const_ref_vert_cat(self&& aObj) noexcept = default;

  self& operator=(const self&) = delete;

  mat_const_ref_vert_cat(UpperMatrix&&, LowerMatrix&&) = delete;

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (i < mu->get_row_count()) {
      return (*mu)(i, j);
    }
    return (*ml)(i - mu->get_row_count(), j);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept {
    return mu->get_row_count() + ml->get_row_count();
  }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept { return mu->get_col_count(); }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(mu->get_row_count() + ml->get_row_count(),
                          mu->get_col_count());
  }
};

/// This function template will vertically concatenate two matrices, by copying them into a
/// composite matrix.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices, by copy.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat_copy(const UpperMatrix& MU, const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, LowerMatrix>(MU, ML);
}

/// This function template will vertically concatenate two non-const matrices, by reference to them in a
/// composite matrix.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices, by reference.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(UpperMatrix& MU, LowerMatrix& ML) {
  return mat_ref_vert_cat<UpperMatrix, LowerMatrix>(MU, ML);
}

/// This function template will vertically concatenate two const matrices, by reference to them in a
/// composite matrix.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices, by const-reference.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(const UpperMatrix& MU, const LowerMatrix& ML) {
  return mat_const_ref_vert_cat<UpperMatrix, LowerMatrix>(MU, ML);
}

/// This function template will vertically concatenate two rvalue matrices, by moving them into a
/// composite matrix. This is an overload that will be selected when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices, by copy.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(UpperMatrix&& MU, LowerMatrix&& ML) {
  return mat_vert_cat<UpperMatrix, LowerMatrix>(std::move(MU), std::move(ML));
}

/// This function template will vertically concatenate one lvalue matrix and one rvalue matrix, by referring
/// to the former and moving the latter into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(UpperMatrix& MU, LowerMatrix&& ML) {
  return mat_vert_cat<mat_sub_block<UpperMatrix>, LowerMatrix>(
      mat_sub_block<UpperMatrix>(MU), std::move(ML));
}

/// This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(UpperMatrix&& MU, LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix>>(
      std::move(MU), mat_sub_block<LowerMatrix>(ML));
}

/// This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(const UpperMatrix& MU, LowerMatrix&& ML) {
  return mat_vert_cat<mat_const_sub_block<UpperMatrix>, LowerMatrix>(
      mat_const_sub_block<UpperMatrix>(MU), std::move(ML));
}

/// This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto vcat(UpperMatrix&& MU, const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix>>(
      std::move(MU), mat_const_sub_block<LowerMatrix>(ML));
}

/// This operator overload template will vertically concatenate two rvalue matrices, by moving them into a
/// composite matrix. This is an overload that will be selected when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices, by copy.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(UpperMatrix&& MU, LowerMatrix&& ML) {
  return mat_vert_cat<UpperMatrix, LowerMatrix>(std::move(MU), std::move(ML));
}

/// This operator overload template will vertically concatenate one lvalue matrix and one rvalue matrix, by referring
/// to the former and moving the latter into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(UpperMatrix& MU, LowerMatrix&& ML) {
  return mat_vert_cat<mat_sub_block<UpperMatrix>, LowerMatrix>(
      mat_sub_block<UpperMatrix>(MU), std::move(ML));
}

/// This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(UpperMatrix&& MU, LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix>>(
      std::move(MU), mat_sub_block<LowerMatrix>(ML));
}

/// This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The value of the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(const UpperMatrix& MU, LowerMatrix&& ML) {
  return mat_vert_cat<mat_const_sub_block<UpperMatrix>, LowerMatrix>(
      mat_const_sub_block<UpperMatrix>(MU), std::move(ML));
}

/// This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring
/// to the latter and moving the former into a composite matrix. This is an overload that will be selected
/// when given rvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The value of the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(UpperMatrix&& MU, const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix>>(
      std::move(MU), mat_const_sub_block<LowerMatrix>(ML));
}

/// This operator overload template will vertically concatenate two lvalue matrices, by referring
/// to them in a composite matrix. This is an overload that will be selected
/// when given non-const lvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(UpperMatrix& MU, LowerMatrix& ML) {
  return mat_ref_vert_cat<UpperMatrix, LowerMatrix>(MU, ML);
}

/// This operator overload template will vertically concatenate two lvalue const matrices, by referring
/// to them in a composite matrix. This is an overload that will be selected
/// when given const lvalues.
/// \note Requires C++0x support for rvalue-references and move-semantics.
/// \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
/// \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
/// \param MU The matrix storing the upper part of the composite matrix.
/// \param ML The matrix storing the lower part of the composite matrix.
/// \return The composite matrix that vertically concatenates the two given matrices.
template <ReadableMatrix UpperMatrix, ReadableMatrix LowerMatrix>
auto operator|(const UpperMatrix& MU, const LowerMatrix& ML) {
  return mat_const_ref_vert_cat<UpperMatrix, LowerMatrix>(MU, ML);
}
}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_COMPOSITE_ADAPTOR_H_
