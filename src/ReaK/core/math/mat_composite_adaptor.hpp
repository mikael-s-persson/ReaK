/**
 * \file mat_composite_adaptor.hpp
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

#ifndef MAT_COMPOSITE_ADAPTOR_HPP
#define MAT_COMPOSITE_ADAPTOR_HPP


#include "mat_alg_general.hpp"
#include "mat_views.hpp"

#include <boost/static_assert.hpp>

namespace ReaK {

/**
 * This class template forms the horizontal concatenation of two matrices, which it stores by value.
 * This class makes the concatenation of the two matrices look as if it was just one matrix (and so,
 * this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both LeftMatrix and RightMatrix, 
 * except for ResizableMatrixConcept.
 * 
 * \tparam LeftMatrix Matrix type for the left matrix.
 * \tparam RightMatrix Matrix type for the right matrix.
 */
template <typename LeftMatrix, typename RightMatrix>
class mat_horiz_cat {
  public:
    typedef mat_horiz_cat<LeftMatrix,RightMatrix> self;
    typedef typename mat_traits<LeftMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<LeftMatrix>::value_type value_type;
    
    typedef typename mat_traits<LeftMatrix>::reference reference;
    typedef typename mat_traits<LeftMatrix>::const_reference const_reference;
    typedef typename mat_traits<LeftMatrix>::pointer pointer;
    typedef typename mat_traits<LeftMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<LeftMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<LeftMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<LeftMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<LeftMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<LeftMatrix>::size_type size_type;
    typedef typename mat_traits<LeftMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<LeftMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    LeftMatrix ml; ///< Holds the left part of the matrix.
    RightMatrix mr; ///< Holds the right part of the matrix.
  public:
    /**
     * Default constructor.
     */
    mat_horiz_cat() : ml(), mr() { };
    
    /**
     * Parametrized constructor.
     * \param aML Matrix to fill the left part of the matrix.
     * \param aMR Matrix to fill the right part of the matrix.
     */
    mat_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR) : ml(aML), mr(aMR) { 
      if(ml.get_row_count() != mr.get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Standard copy-constructor.
     * \param aObj Right-hand-side of the copy.
     */
    mat_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Parametrized move-constructor.
     * \param aML Matrix to fill and be moved into the left part of the matrix.
     * \param aMR Matrix to fill and be moved into the right part of the matrix.
     */
    mat_horiz_cat(LeftMatrix&& aML, RightMatrix&& aMR) : ml(std::move(aML)), mr(std::move(aMR)) { 
      if(ml.get_row_count() != mr.get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Standard move-constructor.
     * \param aObj Right-hand-side of the move.
     */
    mat_horiz_cat(self&& aObj) : ml(std::move(aObj.ml)), mr(std::move(aObj.mr)) { };
#endif
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.ml,rhs.ml);
      swap(lhs.mr,rhs.mr);
    };
    
    /**
     * Standard assignment operator, implemented by copy-and-swap (and move-and-swap).
     */
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    /**
     * Templated assignment operator to assign the content of the matrix with the content 
     * of a matrix of another type (Matrix)
     * \tparam Matrix A readable matrix (models ReadableMatrixConcept)
     * \param rhs Right-hand-side of the assignment.
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != ml.get_row_count()) || (rhs.get_col_count() != ml.get_col_count() + mr.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      ml = sub(rhs)(range(0,ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      mr = sub(rhs)(range(0,mr.get_row_count()-1),range(ml.get_col_count(),ml.get_col_count() + mr.get_col_count()-1));
      return *this;
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
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) const { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return ml.get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return ml.get_col_count() + mr.get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml.get_row_count(),ml.get_col_count() + mr.get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return ml.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != ml.get_row_count()) || (M.get_col_count() != ml.get_col_count() + mr.get_col_count()))
	throw std::range_error("Matrix dimension mismatch.");
      ml += sub(M)(range(0,ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      mr += sub(M)(range(0,mr.get_row_count()-1),range(ml.get_col_count(),ml.get_col_count() + mr.get_col_count()-1));
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != ml.get_row_count()) || (M.get_col_count() != ml.get_col_count() + mr.get_col_count()))
	throw std::range_error("Matrix dimension mismatch.");
      ml -= sub(M)(range(0,ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      mr -= sub(M)(range(0,mr.get_row_count()-1),range(ml.get_col_count(),ml.get_col_count() + mr.get_col_count()-1));
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      ml *= S;
      mr *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
      if((M.get_col_count() != get_col_count()) || (M.get_row_count() != get_col_count()))
	throw std::range_error("Matrix Dimension Mismatch.");
      *this = *this * M;
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename LeftMatrix, typename RightMatrix>
struct is_readable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<LeftMatrix>::value && is_readable_matrix<RightMatrix>::value );
  typedef is_readable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_writable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<LeftMatrix>::value && is_writable_matrix<RightMatrix>::value);
  typedef is_writable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_fully_writable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<LeftMatrix>::value && is_fully_writable_matrix<RightMatrix>::value );
  typedef is_fully_writable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_resizable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct has_allocator_matrix< mat_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<LeftMatrix>::value );
  typedef has_allocator_matrix<LeftMatrix> type;
};




/**
 * This class template forms the horizontal concatenation of two matrices, which it takes by reference 
 * (and stores by pointer, to be copyable). This class makes the concatenation of the two matrices 
 * look as if it was just one matrix (and so, this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both LeftMatrix and RightMatrix, 
 * except for ResizableMatrixConcept.
 * 
 * \tparam LeftMatrix Matrix type for the left matrix.
 * \tparam RightMatrix Matrix type for the right matrix.
 */
template <typename LeftMatrix, typename RightMatrix>
class mat_ref_horiz_cat {
  public:
    typedef mat_ref_horiz_cat<LeftMatrix,RightMatrix> self;
    typedef typename mat_traits<LeftMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<LeftMatrix>::value_type value_type;
    
    typedef typename mat_traits<LeftMatrix>::reference reference;
    typedef typename mat_traits<LeftMatrix>::const_reference const_reference;
    typedef typename mat_traits<LeftMatrix>::pointer pointer;
    typedef typename mat_traits<LeftMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<LeftMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<LeftMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<LeftMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<LeftMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<LeftMatrix>::size_type size_type;
    typedef typename mat_traits<LeftMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<LeftMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    LeftMatrix* ml; ///< Holds the reference to the left part of the matrix.
    RightMatrix* mr; ///< Holds the reference to the right part of the matrix.
  public:
    /**
     * Default constructor.
     */
    mat_ref_horiz_cat() : ml(NULL), mr(NULL) { };
    
    /**
     * Parametrized constructor.
     * \param aML Matrix to become the left part of the matrix.
     * \param aMR Matrix to become the right part of the matrix.
     */
    mat_ref_horiz_cat(LeftMatrix& aML, RightMatrix& aMR) : ml(&aML), mr(&aMR) { 
      if(ml->get_row_count() != mr->get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Copy-constructor (shallow-copy).
     */
    mat_ref_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Move-constructor (shallow-move).
     */
    mat_ref_horiz_cat(self&& aObj) : ml(aObj.ml), mr(aObj.mr) { };
#endif
    
    /**
     * Standard swap function (shallow).
     */
    friend void swap(self& lhs,self& rhs) {
      using std::swap;
      swap(lhs.ml,rhs.ml);
      swap(lhs.mr,rhs.mr);
      return;
    };
    
    /**
     * Templated assignment operator to assign the content of the matrix with the content 
     * of a matrix of another type (Matrix)
     * \tparam Matrix A readable matrix (models ReadableMatrixConcept)
     * \param rhs Right-hand-side of the assignment.
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != ml->get_row_count()) || (rhs.get_col_count() != ml->get_col_count() + mr->get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      (*ml) = sub(rhs)(range(0,ml->get_row_count()-1),range(0,ml->get_col_count()-1));
      (*mr) = sub(rhs)(range(0,mr->get_row_count()-1),range(ml->get_col_count(),ml->get_col_count() + mr->get_col_count()-1));
      return *this;
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
      if(j < ml->get_col_count())
	return (*ml)(i,j);
      else
	return (*mr)(i,j - ml->get_col_count());
    };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) const { 
      if(j < ml->get_col_count())
	return (*ml)(i,j);
      else
	return (*mr)(i,j - ml->get_col_count());
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return ml->get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return ml->get_col_count() + mr->get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml->get_row_count(),ml->get_col_count() + mr->get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return ml->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != ml->get_row_count()) || (M.get_col_count() != ml->get_col_count() + mr->get_col_count()))
	throw std::range_error("Matrix dimension mismatch.");
      (*ml) += sub(M)(range(0,ml->get_row_count()-1),range(0,ml->get_col_count()-1));
      (*mr) += sub(M)(range(0,mr->get_row_count()-1),range(ml->get_col_count(),ml->get_col_count() + mr->get_col_count()-1));
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != ml->get_row_count()) || (M.get_col_count() != ml->get_col_count() + mr->get_col_count()))
	throw std::range_error("Matrix dimension mismatch.");
      (*ml) -= sub(M)(range(0,ml->get_row_count()-1),range(0,ml->get_col_count()-1));
      (*mr) -= sub(M)(range(0,mr->get_row_count()-1),range(ml->get_col_count(),ml->get_col_count() + mr->get_col_count()-1));
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      (*ml) *= S;
      (*mr) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
      if((M.get_col_count() != get_col_count()) || (M.get_row_count() != get_col_count()))
	throw std::range_error("Matrix Dimension Mismatch.");
      *this = *this * M;
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename LeftMatrix, typename RightMatrix>
struct is_readable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<LeftMatrix>::value && is_readable_matrix<RightMatrix>::value );
  typedef is_readable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_writable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<LeftMatrix>::value && is_writable_matrix<RightMatrix>::value);
  typedef is_writable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_fully_writable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<LeftMatrix>::value && is_fully_writable_matrix<RightMatrix>::value );
  typedef is_fully_writable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_resizable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct has_allocator_matrix< mat_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<LeftMatrix>::value );
  typedef has_allocator_matrix<LeftMatrix> type;
};







/**
 * This class template forms the horizontal concatenation of two matrices, which it takes by const-reference 
 * (and stores by const-pointer, to be copyable). This class makes the concatenation of the two matrices 
 * look as if it was just one matrix (and so, this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both LeftMatrix and RightMatrix, 
 * except for ResizableMatrixConcept.
 * 
 * \tparam LeftMatrix Matrix type for the left matrix.
 * \tparam RightMatrix Matrix type for the right matrix.
 */
template <typename LeftMatrix, typename RightMatrix>
class mat_const_ref_horiz_cat {
  public:
    typedef mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> self;
    typedef typename mat_traits<LeftMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<LeftMatrix>::value_type value_type;
    
    typedef typename mat_traits<LeftMatrix>::reference reference;
    typedef typename mat_traits<LeftMatrix>::const_reference const_reference;
    typedef typename mat_traits<LeftMatrix>::pointer pointer;
    typedef typename mat_traits<LeftMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<LeftMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<LeftMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<LeftMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<LeftMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<LeftMatrix>::size_type size_type;
    typedef typename mat_traits<LeftMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<LeftMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    const LeftMatrix* ml;
    const RightMatrix* mr;
    
    self& operator=(const self&);
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_ref_horiz_cat(LeftMatrix&&, RightMatrix&&);
#endif
  public:
    /**
     * Parametrized constructor.
     * \param aML Matrix to fill the left part of the matrix.
     * \param aMR Matrix to fill the right part of the matrix.
     */
    mat_const_ref_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR) : ml(&aML), mr(&aMR) { 
      if(ml->get_row_count() != mr->get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Standard copy-constructor (shallow).
     * \param aObj Right-hand-side of the copy.
     */
    mat_const_ref_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES    
    /**
     * Standard move-constructor.
     * \param aObj Right-hand-side of the move.
     */
    mat_const_ref_horiz_cat(self&& aObj) : ml(aObj.ml), mr(aObj.mr) { };
#endif
    
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
    value_type operator()(size_type i,size_type j) const { 
      if(j < ml->get_col_count())
	return (*ml)(i,j);
      else
	return (*mr)(i,j - ml->get_col_count());
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return ml->get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return ml->get_col_count() + mr->get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml->get_row_count(),ml->get_col_count() + mr->get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return ml->get_allocator(); };
    
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename LeftMatrix, typename RightMatrix>
struct is_readable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<LeftMatrix>::value && is_readable_matrix<RightMatrix>::value );
  typedef is_readable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_writable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_writable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_fully_writable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_fully_writable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct is_resizable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false);
  typedef is_resizable_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > type;
};

template <typename LeftMatrix, typename RightMatrix>
struct has_allocator_matrix< mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<LeftMatrix>::value );
  typedef has_allocator_matrix<LeftMatrix> type;
};




/**
 * This function template will horizontally concatenate two matrices, by copying them into a 
 * composite matrix.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by copy.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix,RightMatrix> >::type hcat_copy(const LeftMatrix& ML,const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

/**
 * This function template will horizontally concatenate two non-const matrices, by reference to them in a 
 * composite matrix.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by reference.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_ref_horiz_cat<LeftMatrix,RightMatrix> >::type hcat(LeftMatrix& ML,RightMatrix& MR) {
  return mat_ref_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

/**
 * This function template will horizontally concatenate two const matrices, by reference to them in a 
 * composite matrix.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by const-reference.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> >::type hcat(const LeftMatrix& ML,const RightMatrix& MR) {
  return mat_const_ref_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

#ifndef RK_ENABLE_CXX0X_FEATURES

/**
 * This operator template will horizontally concatenate two non-const matrices, by copying them into a 
 * composite matrix. Making a copy is the only safe option in C++03 because it is not safe to assume
 * that a const-reference is anything but an rvalue (unless explicitly implied by the use of the named
 * function templates (vcat)).
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by copy.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix,RightMatrix> >::type operator&(const LeftMatrix& ML,const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

#else

/**
 * This function template will horizontally concatenate two rvalue matrices, by moving them into a 
 * composite matrix. This is an overload that will be selected when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by copy.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix,RightMatrix> >::type hcat(LeftMatrix&& ML,RightMatrix&& MR) {
  return mat_horiz_cat<LeftMatrix,RightMatrix>(std::move(ML),std::move(MR));
};

/**
 * This function template will horizontally concatenate one lvalue matrix and one rvalue matrix, by referring 
 * to the former and moving the latter into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<mat_sub_block<LeftMatrix>,RightMatrix> >::type hcat(LeftMatrix& ML,RightMatrix&& MR) {
  return mat_horiz_cat<mat_sub_block<LeftMatrix>,RightMatrix>(mat_sub_block<LeftMatrix>(ML),std::move(MR));
};

/**
 * This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix> > >::type hcat(LeftMatrix&& ML,RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix> >(std::move(ML),mat_sub_block<RightMatrix>(MR));
};

/**
 * This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<mat_const_sub_block<LeftMatrix>,RightMatrix> >::type hcat(const LeftMatrix& ML,RightMatrix&& MR) {
  return mat_horiz_cat<mat_const_sub_block<LeftMatrix>,RightMatrix>(mat_const_sub_block<LeftMatrix>(ML),std::move(MR));
};

/**
 * This function template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix> > >::type hcat(LeftMatrix&& ML,const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix> >(std::move(ML),mat_const_sub_block<RightMatrix>(MR));
};



/**
 * This operator overload template will horizontally concatenate two rvalue matrices, by moving them into a 
 * composite matrix. This is an overload that will be selected when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices, by copy.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix,RightMatrix> >::type operator&(LeftMatrix&& ML,RightMatrix&& MR) {
  return mat_horiz_cat<LeftMatrix,RightMatrix>(std::move(ML),std::move(MR));
};

/**
 * This operator overload template will horizontally concatenate one lvalue matrix and one rvalue matrix, by referring 
 * to the former and moving the latter into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<mat_sub_block<LeftMatrix>,RightMatrix> >::type operator&(LeftMatrix& ML,RightMatrix&& MR) {
  return mat_horiz_cat<mat_sub_block<LeftMatrix>,RightMatrix>(mat_sub_block<LeftMatrix>(ML),std::move(MR));
};

/**
 * This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix> > >::type operator&(LeftMatrix&& ML,RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_sub_block<RightMatrix> >(std::move(ML),mat_sub_block<RightMatrix>(MR));
};

/**
 * This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The value of the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<mat_const_sub_block<LeftMatrix>,RightMatrix> >::type operator&(const LeftMatrix& ML,RightMatrix&& MR) {
  return mat_horiz_cat<mat_const_sub_block<LeftMatrix>,RightMatrix>(mat_const_sub_block<LeftMatrix>(ML),std::move(MR));
};

/**
 * This operator overload template will horizontally concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The value of the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix> > >::type operator&(LeftMatrix&& ML,const RightMatrix& MR) {
  return mat_horiz_cat<LeftMatrix, mat_const_sub_block<RightMatrix> >(std::move(ML),mat_const_sub_block<RightMatrix>(MR));
};


/**
 * This operator overload template will horizontally concatenate two lvalue matrices, by referring 
 * to them in a composite matrix. This is an overload that will be selected 
 * when given non-const lvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_ref_horiz_cat<LeftMatrix,RightMatrix> >::type operator&(LeftMatrix& ML,RightMatrix& MR) {
  return mat_ref_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

/**
 * This operator overload template will horizontally concatenate two lvalue const matrices, by referring 
 * to them in a composite matrix. This is an overload that will be selected 
 * when given const lvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam LeftMatrix Matrix type for the left part of the composite matrix.
 * \tparam RightMatrix Matrix type for the right part of the composite matrix.
 * \param MU The matrix storing the left part of the composite matrix.
 * \param ML The matrix storing the right part of the composite matrix.
 * \return The composite matrix that horizontally concatenates the two given matrices.
 */
template <typename LeftMatrix, typename RightMatrix>
typename boost::enable_if_c< is_readable_matrix<LeftMatrix>::value &&
                             is_readable_matrix<RightMatrix>::value,
mat_const_ref_horiz_cat<LeftMatrix,RightMatrix> >::type operator&(const LeftMatrix& ML,const RightMatrix& MR) {
  return mat_const_ref_horiz_cat<LeftMatrix,RightMatrix>(ML,MR);
};

#endif











/**
 * This class template forms the vertical concatenation of two matrices, which it stores by value.
 * This class makes the concatenation of the two matrices look as if it was just one matrix (and so,
 * this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both LeftMatrix and RightMatrix, 
 * except for ResizableMatrixConcept.
 * 
 * \tparam UpperMatrix Matrix type for the upper matrix.
 * \tparam LowerMatrix Matrix type for the lower matrix.
 */
template <typename UpperMatrix, typename LowerMatrix>
class mat_vert_cat {
  public:
    typedef mat_vert_cat<UpperMatrix,LowerMatrix> self;
    typedef typename mat_traits<UpperMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<UpperMatrix>::value_type value_type;
    
    typedef typename mat_traits<UpperMatrix>::reference reference;
    typedef typename mat_traits<UpperMatrix>::const_reference const_reference;
    typedef typename mat_traits<UpperMatrix>::pointer pointer;
    typedef typename mat_traits<UpperMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<UpperMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<UpperMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<UpperMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<UpperMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<UpperMatrix>::size_type size_type;
    typedef typename mat_traits<UpperMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<UpperMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    UpperMatrix mu; ///< Holds the upper part of the matrix.
    LowerMatrix ml; ///< Holds the lower part of the matrix.
  public:
    /**
     * Default constructor.
     */
    mat_vert_cat() : mu(), ml() { };
    
    /**
     * Parametrized constructor.
     * \param aMU Matrix to fill the upper part of the matrix.
     * \param aML Matrix to fill the lower part of the matrix.
     */
    mat_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML) : mu(aMU), ml(aML) { 
      if(ml.get_col_count() != mu.get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Copy-constructor.
     */
    mat_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Parametrized Move-constructor.
     * \param aMU Matrix to fill the upper part of the matrix.
     * \param aML Matrix to fill the lower part of the matrix.
     */
    mat_vert_cat(UpperMatrix&& aMU, LowerMatrix&& aML) : mu(std::move(aMU)), ml(std::move(aML)) { 
      if(ml.get_col_count() != mu.get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Move-constructor.
     */
    mat_vert_cat(self&& aObj) : mu(std::move(aObj.mu)), ml(std::move(aObj.ml)) { };
#endif
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mu,rhs.mu);
      swap(lhs.ml,rhs.ml);
    };
    
    /**
     * Standard assignment operator (copy-and-swap and move-and-swap).
     */
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    /**
     * Templated assignment operator to assign the content of the matrix with the content 
     * of a matrix of another type (Matrix)
     * \tparam Matrix A readable matrix (models ReadableMatrixConcept)
     * \param rhs Right-hand-side of the assignment.
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != mu.get_row_count() + ml.get_row_count()) || (rhs.get_col_count() != mu.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      mu = sub(rhs)(range(0,mu.get_row_count()-1),range(0,mu.get_col_count()-1));
      ml = sub(rhs)(range(mu.get_row_count(),mu.get_row_count() + ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      return *this;
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
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) const {  
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return mu.get_row_count() + ml.get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return mu.get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu.get_row_count() + ml.get_row_count(),mu.get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return mu.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != mu.get_row_count() + ml.get_row_count()) || (M.get_col_count() != mu.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      mu += sub(M)(range(0,mu.get_row_count()-1),range(0,mu.get_col_count()-1));
      ml += sub(M)(range(mu.get_row_count(),mu.get_row_count() + ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != mu.get_row_count() + ml.get_row_count()) || (M.get_col_count() != mu.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      mu -= sub(M)(range(0,mu.get_row_count()-1),range(0,mu.get_col_count()-1));
      ml -= sub(M)(range(mu.get_row_count(),mu.get_row_count() + ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      mu *= S;
      ml *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
      if((M.get_col_count() != get_col_count()) || (M.get_row_count() != get_col_count()))
	throw std::range_error("Matrix Dimension Mismatch.");
      *this = *this * M;
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename UpperMatrix, typename LowerMatrix>
struct is_readable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<UpperMatrix>::value && is_readable_matrix<LowerMatrix>::value );
  typedef is_readable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_writable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<UpperMatrix>::value && is_writable_matrix<LowerMatrix>::value);
  typedef is_writable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_fully_writable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<UpperMatrix>::value && is_fully_writable_matrix<LowerMatrix>::value );
  typedef is_fully_writable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_resizable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct has_allocator_matrix< mat_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<UpperMatrix>::value );
  typedef has_allocator_matrix<UpperMatrix> type;
};




/**
 * This class template forms the vertical concatenation of two matrices, which it takes by reference 
 * (and stores by pointer, to be copyable). This class makes the concatenation of the two matrices 
 * look as if it was just one matrix (and so, this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both UpperMatrix and LowerMatrix, 
 * except for ResizableMatrixConcept.
 * 
 * \tparam UpperMatrix Matrix type for the left matrix.
 * \tparam LowerMatrix Matrix type for the right matrix.
 */
template <typename UpperMatrix, typename LowerMatrix>
class mat_ref_vert_cat {
  public:
    typedef mat_ref_vert_cat<UpperMatrix,LowerMatrix> self;
    typedef typename mat_traits<UpperMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<UpperMatrix>::value_type value_type;
    
    typedef typename mat_traits<UpperMatrix>::reference reference;
    typedef typename mat_traits<UpperMatrix>::const_reference const_reference;
    typedef typename mat_traits<UpperMatrix>::pointer pointer;
    typedef typename mat_traits<UpperMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<UpperMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<UpperMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<UpperMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<UpperMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<UpperMatrix>::size_type size_type;
    typedef typename mat_traits<UpperMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<UpperMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    UpperMatrix* mu; ///< Refers to the upper part of the matrix.
    LowerMatrix* ml; ///< Refers to the lower part of the matrix.
  public:
    /**
     * Parametrized constructor.
     * \param aMU Matrix to become the upper part of the matrix.
     * \param aML Matrix to become the lower part of the matrix.
     */
    mat_ref_vert_cat(UpperMatrix& aMU, LowerMatrix& aML) : mu(&aMU), ml(&aML) { 
      if(ml->get_col_count() != mu->get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Standard copy-constructor (shallow-copy).
     */
    mat_ref_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Standard move-constructor (shallow-move).
     */
    mat_ref_vert_cat(self&& aObj) : mu(aObj.mu), ml(aObj.ml) { };
#endif
    
    /**
     * Templated assignment operator to assign the content of the matrix with the content 
     * of a matrix of another type (Matrix)
     * \tparam Matrix A readable matrix (models ReadableMatrixConcept)
     * \param rhs Right-hand-side of the assignment.
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != mu->get_row_count() + ml->get_row_count()) || (rhs.get_col_count() != mu->get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      *mu = sub(rhs)(range(0,mu->get_row_count()-1),range(0,mu->get_col_count()-1));
      *ml = sub(rhs)(range(mu->get_row_count(),mu->get_row_count() + ml->get_row_count()-1),range(0,ml->get_col_count()-1));
      return *this;
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
      if(i < mu->get_row_count())
	return (*mu)(i,j);
      else
	return (*ml)(i - mu->get_row_count(),j);
    };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    value_type operator()(size_type i,size_type j) const {  
      if(i < mu->get_row_count())
	return (*mu)(i,j);
      else
	return (*ml)(i - mu->get_row_count(),j);
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return mu->get_row_count() + ml->get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return mu->get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu->get_row_count() + ml->get_row_count(),mu->get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return mu->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != mu->get_row_count() + ml->get_row_count()) || (M.get_col_count() != mu->get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      *mu += sub(M)(range(0,mu.get_row_count()-1),range(0,mu.get_col_count()-1));
      *ml += sub(M)(range(mu.get_row_count(),mu.get_row_count() + ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_row_count() != mu->get_row_count() + ml->get_row_count()) || (M.get_col_count() != mu->get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      *mu -= sub(M)(range(0,mu.get_row_count()-1),range(0,mu.get_col_count()-1));
      *ml -= sub(M)(range(mu.get_row_count(),mu.get_row_count() + ml.get_row_count()-1),range(0,ml.get_col_count()-1));
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      mu *= S;
      ml *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
      if((M.get_col_count() != get_col_count()) || (M.get_row_count() != get_col_count()))
	throw std::range_error("Matrix Dimension Mismatch.");
      *this = *this * M;
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename UpperMatrix, typename LowerMatrix>
struct is_readable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<UpperMatrix>::value && is_readable_matrix<LowerMatrix>::value );
  typedef is_readable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_writable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<UpperMatrix>::value && is_writable_matrix<LowerMatrix>::value);
  typedef is_writable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_fully_writable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<UpperMatrix>::value && is_fully_writable_matrix<LowerMatrix>::value );
  typedef is_fully_writable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_resizable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct has_allocator_matrix< mat_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<UpperMatrix>::value );
  typedef has_allocator_matrix<UpperMatrix> type;
};




/**
 * This class template forms the vertical concatenation of two matrices, which it takes by const-reference 
 * (and stores by const-pointer, to be copyable). This class makes the concatenation of the two matrices 
 * look as if it was just one matrix (and so, this class is an adaptor).
 * 
 * Models: ReadableMatrixConcept and all matrix concepts modeled by both UpperMatrix and LowerMatrix, 
 * except for ResizableMatrixConcept, WritableMatrixConcept, and FullyWritableMatrixConcept.
 * 
 * \tparam UpperMatrix Matrix type for the left matrix.
 * \tparam LowerMatrix Matrix type for the right matrix.
 */
template <typename UpperMatrix, typename LowerMatrix>
class mat_const_ref_vert_cat {
  public:
    typedef mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> self;
    typedef typename mat_traits<UpperMatrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<UpperMatrix>::value_type value_type;
    
    typedef typename mat_traits<UpperMatrix>::reference reference;
    typedef typename mat_traits<UpperMatrix>::const_reference const_reference;
    typedef typename mat_traits<UpperMatrix>::pointer pointer;
    typedef typename mat_traits<UpperMatrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<UpperMatrix>::col_iterator col_iterator;
    typedef typename mat_traits<UpperMatrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<UpperMatrix>::row_iterator row_iterator;
    typedef typename mat_traits<UpperMatrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<UpperMatrix>::size_type size_type;
    typedef typename mat_traits<UpperMatrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<UpperMatrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    const UpperMatrix* mu;
    const LowerMatrix* ml;
    
    self& operator=(const self&);
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_ref_vert_cat(UpperMatrix&&, LowerMatrix&&);
#endif
  public:
    /**
     * Parametrized constructor.
     * \param aMU The matrix which will become the upper part of this matrix.
     * \param aML The matrix which will become the lower part of this matrix.
     */
    mat_const_ref_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML) : mu(&aMU), ml(&aML) { 
      if(ml->get_col_count() != mu->get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    /**
     * Standard copy-constructor (shallow-copy).
     */
    mat_const_ref_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Standard move-constructor (shallow-move).
     */
    mat_const_ref_vert_cat(self&& aObj) : mu(aObj.mu), ml(aObj.ml) { };
#endif
    
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
    value_type operator()(size_type i,size_type j) const {  
      if(i < mu->get_row_count())
	return (*mu)(i,j);
      else
	return (*ml)(i - mu->get_row_count(),j);
    };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return mu->get_row_count() + ml->get_row_count(); };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return mu->get_col_count(); };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu->get_row_count() + ml->get_row_count(),mu->get_col_count()); };

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return mu->get_allocator(); };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::rectangular> operator -() const {
      mat<value_type,mat_structure::rectangular> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    /**
     * Transposes the matrix M.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename UpperMatrix, typename LowerMatrix>
struct is_readable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<UpperMatrix>::value && is_readable_matrix<LowerMatrix>::value );
  typedef is_readable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_writable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_fully_writable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct is_resizable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > type;
};

template <typename UpperMatrix, typename LowerMatrix>
struct has_allocator_matrix< mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<UpperMatrix>::value );
  typedef has_allocator_matrix<UpperMatrix> type;
};





/**
 * This function template will vertically concatenate two matrices, by copying them into a 
 * composite matrix.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by copy.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix,LowerMatrix> >::type vcat_copy(const UpperMatrix& MU,const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

/**
 * This function template will vertically concatenate two non-const matrices, by reference to them in a 
 * composite matrix.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by reference.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_ref_vert_cat<UpperMatrix,LowerMatrix> >::type vcat(UpperMatrix& MU,LowerMatrix& ML) {
  return mat_ref_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

/**
 * This function template will vertically concatenate two const matrices, by reference to them in a 
 * composite matrix.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by const-reference.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> >::type vcat(const UpperMatrix& MU,const LowerMatrix& ML) {
  return mat_const_ref_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

#ifndef RK_ENABLE_CXX0X_FEATURES

/**
 * This operator template will vertically concatenate two non-const matrices, by copying them into a 
 * composite matrix. Making a copy is the only safe option in C++03 because it is not safe to assume
 * that a const-reference is anything but an rvalue (unless explicitly implied by the use of the named
 * function templates (vcat)).
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by copy.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix,LowerMatrix> >::type operator|(const UpperMatrix& MU,const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

#else

/**
 * This function template will vertically concatenate two rvalue matrices, by moving them into a 
 * composite matrix. This is an overload that will be selected when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by copy.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix,LowerMatrix> >::type vcat(UpperMatrix&& MU,LowerMatrix&& ML) {
  return mat_vert_cat<UpperMatrix,LowerMatrix>(std::move(MU),std::move(ML));
};

/**
 * This function template will vertically concatenate one lvalue matrix and one rvalue matrix, by referring 
 * to the former and moving the latter into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat< mat_sub_block<UpperMatrix> ,LowerMatrix> >::type vcat(UpperMatrix& MU,LowerMatrix&& ML) {
  return mat_vert_cat<mat_sub_block<UpperMatrix>,LowerMatrix>(mat_sub_block<UpperMatrix>(MU),std::move(ML));
};

/**
 * This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix> > >::type vcat(UpperMatrix&& MU,LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix> >(std::move(MU),mat_sub_block<LowerMatrix>(ML));
};

/**
 * This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat< mat_const_sub_block<UpperMatrix> ,LowerMatrix> >::type vcat(const UpperMatrix& MU,LowerMatrix&& ML) {
  return mat_vert_cat<mat_const_sub_block<UpperMatrix>,LowerMatrix>(mat_const_sub_block<UpperMatrix>(MU),std::move(ML));
};

/**
 * This function template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix> > >::type vcat(UpperMatrix&& MU,const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix> >(std::move(MU),mat_const_sub_block<LowerMatrix>(ML));
};




/**
 * This operator overload template will vertically concatenate two rvalue matrices, by moving them into a 
 * composite matrix. This is an overload that will be selected when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices, by copy.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix,LowerMatrix> >::type operator|(UpperMatrix&& MU,LowerMatrix&& ML) {
  return mat_vert_cat<UpperMatrix,LowerMatrix>(std::move(MU),std::move(ML));
};

/**
 * This operator overload template will vertically concatenate one lvalue matrix and one rvalue matrix, by referring 
 * to the former and moving the latter into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<mat_sub_block<UpperMatrix>,LowerMatrix> >::type operator|(UpperMatrix& MU,LowerMatrix&& ML) {
  return mat_vert_cat<mat_sub_block<UpperMatrix>,LowerMatrix>(mat_sub_block<UpperMatrix>(MU),std::move(ML));
};

/**
 * This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix> > >::type operator|(UpperMatrix&& MU,LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_sub_block<LowerMatrix> >(std::move(MU),mat_sub_block<LowerMatrix>(ML));
};

/**
 * This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The value of the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<mat_const_sub_block<UpperMatrix>,LowerMatrix> >::type operator|(const UpperMatrix& MU,LowerMatrix&& ML) {
  return mat_vert_cat<mat_const_sub_block<UpperMatrix>,LowerMatrix>(mat_const_sub_block<UpperMatrix>(MU),std::move(ML));
};

/**
 * This operator overload template will vertically concatenate one rvalue matrix and one lvalue matrix, by referring 
 * to the latter and moving the former into a composite matrix. This is an overload that will be selected 
 * when given rvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The value of the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix> > >::type operator|(UpperMatrix&& MU,const LowerMatrix& ML) {
  return mat_vert_cat<UpperMatrix, mat_const_sub_block<LowerMatrix> >(std::move(MU),mat_const_sub_block<LowerMatrix>(ML));
};





/**
 * This operator overload template will vertically concatenate two lvalue matrices, by referring 
 * to them in a composite matrix. This is an overload that will be selected 
 * when given non-const lvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_ref_vert_cat<UpperMatrix,LowerMatrix> >::type operator|(UpperMatrix& MU,LowerMatrix& ML) {
  return mat_ref_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

/**
 * This operator overload template will vertically concatenate two lvalue const matrices, by referring 
 * to them in a composite matrix. This is an overload that will be selected 
 * when given const lvalues.
 * \note Requires C++0x support for rvalue-references and move-semantics.
 * \tparam UpperMatrix Matrix type for the upper part of the composite matrix.
 * \tparam LowerMatrix Matrix type for the lower part of the composite matrix.
 * \param MU The matrix storing the upper part of the composite matrix.
 * \param ML The matrix storing the lower part of the composite matrix.
 * \return The composite matrix that vertically concatenates the two given matrices.
 */
template <typename UpperMatrix, typename LowerMatrix>
typename boost::enable_if_c< is_readable_matrix<UpperMatrix>::value &&
                             is_readable_matrix<LowerMatrix>::value,
mat_const_ref_vert_cat<UpperMatrix,LowerMatrix> >::type operator|(const UpperMatrix& MU,const LowerMatrix& ML) {
  return mat_const_ref_vert_cat<UpperMatrix,LowerMatrix>(MU,ML);
};

#endif










};

#endif






















