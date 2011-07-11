
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

#ifndef MAT_VIEWS_HPP
#define MAT_VIEWS_HPP

#include "mat_alg_general.hpp"
#include <boost/static_assert.hpp>

namespace ReaK {


template <typename T>
std::pair<T,T> range(T aFirst,T aLast) {
  return std::pair<T,T>(aFirst,aLast);
};

  
  
  

template <typename Matrix>
class mat_copy_sub_block {
  public:    
    
    typedef mat_copy_sub_block<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    Matrix m;
    size_type rowOffset;
    size_type colOffset;
    size_type rowCount;
    size_type colCount;
  public:
    mat_copy_sub_block() : m(), rowOffset(0), colOffset(0), rowCount(0), colCount(0) { };
    
    explicit mat_copy_sub_block(const Matrix& aM) : m(aM), rowOffset(0), colOffset(0), rowCount(aM.get_row_count()), colCount(aM.get_col_count()) { };
    
    mat_copy_sub_block(const Matrix& aM, 
		       size_type aRowCount, 
		       size_type aColCount,
		       size_type aRowOffset = 0,
                       size_type aColOffset = 0) : m(aM), rowOffset(aRowOffset), colOffset(aColOffset), rowCount(aRowCount), colCount(aColCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    explicit mat_copy_sub_block(Matrix&& aM) : m(std::move(aM)), rowOffset(0), colOffset(0), rowCount(0), colCount(0) { rowCount = m.get_row_count(); colCount = m.get_col_count(); };
    
    mat_copy_sub_block(Matrix&& aM, 
		       size_type aRowCount, 
		       size_type aColCount,
		       size_type aRowOffset = 0,
                       size_type aColOffset = 0) : m(std::move(aM)), rowOffset(aRowOffset), colOffset(aColOffset), rowCount(aRowCount), colCount(aColCount) { };
#endif

    mat_copy_sub_block(const self& aObj) : m(aObj.m), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_copy_sub_block(self&& aObj) : m(std::move(aObj.m)), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.rowOffset,rhs.rowOffset);
      swap(lhs.colOffset,rhs.colOffset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value &&
                                 !boost::is_same<Matrix2,self>::value ,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  m(rowOffset + i,colOffset + j) = rhs(i,j);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return m(rowOffset + i, colOffset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return m(rowOffset + i, colOffset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return m.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  m(rowOffset + i, colOffset + j) += M(i,j);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  m(rowOffset + i, colOffset + j) -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  m(rowOffset + i, colOffset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != colCount) || (M.get_row_count() != colCount))
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename Matrix>
struct is_readable_matrix< mat_copy_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_writable_matrix< mat_copy_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_fully_writable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_fully_writable_matrix< mat_copy_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_fully_writable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_resizable_matrix< mat_copy_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_copy_sub_block<Matrix> > type;
};

template <typename Matrix>
struct has_allocator_matrix< mat_copy_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};

  
  
  


template <typename Matrix>
class mat_sub_block {
  public:    
    
    typedef mat_sub_block<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    Matrix* m;
    size_type rowOffset;
    size_type colOffset;
    size_type rowCount;
    size_type colCount;
  public:
    explicit mat_sub_block(Matrix& aM) : m(&aM), rowOffset(0), colOffset(0), rowCount(aM.get_row_count()), colCount(aM.get_col_count()) { };
    
    mat_sub_block(Matrix& aM, 
		  size_type aRowCount, 
		  size_type aColCount,
		  size_type aRowOffset = 0,
                  size_type aColOffset = 0) : m(&aM), rowOffset(aRowOffset), colOffset(aColOffset), rowCount(aRowCount), colCount(aColCount) { };
   
    mat_sub_block(const self& aObj) : m(aObj.m), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
   
    mat_sub_block(self&& aObj) : m(aObj.m), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.rowOffset,rhs.rowOffset);
      swap(lhs.colOffset,rhs.colOffset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
        
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value ,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  (*m)(rowOffset + i,colOffset + j) = rhs(i,j);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*m)(rowOffset + i, colOffset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return (*m)(rowOffset + i, colOffset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  (*m)(rowOffset + i, colOffset + j) += M(i,j);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  (*m)(rowOffset + i, colOffset + j) -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i)
	  (*m)(rowOffset + i, colOffset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != colCount) || (M.get_row_count() != colCount))
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename Matrix>
struct is_readable_matrix< mat_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_writable_matrix< mat_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_fully_writable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_fully_writable_matrix< mat_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_fully_writable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_resizable_matrix< mat_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_sub_block<Matrix> > type;
};

template <typename Matrix>
struct has_allocator_matrix< mat_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};





template <typename Matrix>
class mat_const_sub_block {
  public:    
    
    typedef mat_const_sub_block<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    const Matrix* m;
    size_type rowOffset;
    size_type colOffset;
    size_type rowCount;
    size_type colCount;
    
    self& operator=(const self&);  
#ifdef RK_ENABLE_CXX0X_FEATURES
    explicit mat_const_sub_block(Matrix&&);
 
    mat_const_sub_block(Matrix&&, size_type, size_type, size_type aRowOffset = 0, size_type aColOffset = 0);
#endif
    
  public:
    explicit mat_const_sub_block(const Matrix& aM) : m(&aM), rowOffset(0), colOffset(0), rowCount(aM.get_row_count()), colCount(aM.get_col_count()) { };
    
    mat_const_sub_block(const Matrix& aM, 
			size_type aRowCount, 
			size_type aColCount,
			size_type aRowOffset = 0,
			size_type aColOffset = 0) : m(&aM), rowOffset(aRowOffset), colOffset(aColOffset), rowCount(aRowCount), colCount(aColCount) { };
    
    mat_const_sub_block(const self& aObj) : m(aObj.m), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_block(self&& aObj) : m(aObj.m), rowOffset(aObj.rowOffset), colOffset(aObj.colOffset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.rowOffset,rhs.rowOffset);
      swap(lhs.colOffset,rhs.colOffset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const { 
      return (*m)(rowOffset + i, colOffset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return (*m).get_allocator(); };
    
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
    friend mat<value_type,mat_structure::rectangular> transpose_move(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.colCount, M.rowCount);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
};


template <typename Matrix>
struct is_readable_matrix< mat_const_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix>
struct is_writable_matrix< mat_const_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_sub_block<Matrix> > type;
};

template <typename Matrix>
struct is_fully_writable_matrix< mat_const_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_sub_block<Matrix> > type;
};

template <typename Matrix>
struct is_resizable_matrix< mat_const_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_sub_block<Matrix> > type;
};

template <typename Matrix>
struct has_allocator_matrix< mat_const_sub_block<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};






template <typename Matrix>
struct mat_copy_sub_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;

#ifndef RK_ENABLE_CXX0X_FEATURES
  Matrix m;
  mat_copy_sub_block_factory(const Matrix& aM) : m(aM) { };
  mat_copy_sub_block<Matrix> operator()(const std::pair<size_type,size_type>& rows,
					const std::pair<size_type,size_type>& cols) {
    return mat_copy_sub_block<Matrix>(m,rows.second - rows.first + 1,cols.second - cols.first + 1,rows.first,cols.first);
  };
#else
  Matrix m;
  mat_copy_sub_block_factory(Matrix&& aM) : m(std::move(aM)) { };
  mat_copy_sub_block<Matrix> operator()(const std::pair<size_type,size_type>& rows,
					const std::pair<size_type,size_type>& cols) {
    return mat_copy_sub_block<Matrix>(std::move(m),rows.second - rows.first + 1,cols.second - cols.first + 1,rows.first,cols.first);
  };
#endif
};

template <typename Matrix>
struct mat_sub_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  Matrix& m;
  mat_sub_block_factory(Matrix& aM) : m(aM) { };
  mat_sub_block<Matrix> operator()(const std::pair<size_type,size_type>& rows,
				   const std::pair<size_type,size_type>& cols) {
    return mat_sub_block<Matrix>(m,rows.second - rows.first + 1,cols.second - cols.first + 1,rows.first,cols.first);
  };
};

template <typename Matrix>
struct mat_const_sub_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  const Matrix& m;
  mat_const_sub_block_factory(const Matrix& aM) : m(aM) { };
  mat_const_sub_block<Matrix> operator()(const std::pair<size_type,size_type>& rows,
				         const std::pair<size_type,size_type>& cols) {
    return mat_const_sub_block<Matrix>(m,rows.second - rows.first + 1,cols.second - cols.first + 1,rows.first,cols.first);
  };
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_sub_block_factory<Matrix> >::type sub(Matrix& M) { return mat_sub_block_factory<Matrix>(M); };

template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_const_sub_block_factory<Matrix> >::type sub(const Matrix& M) { return mat_const_sub_block_factory<Matrix>(M); };

template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_copy_sub_block_factory<Matrix> >::type sub_copy(const Matrix& M) { return mat_copy_sub_block_factory<Matrix>(M); };

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_copy_sub_block_factory<Matrix> >::type sub(Matrix&& M) { return mat_copy_sub_block_factory<Matrix>(std::move(M)); };
#endif



















template <typename Matrix, mat_structure::tag Structure = mat_traits<Matrix>::structure >
class mat_copy_sub_sym_block { char invalid_matrix_type[0]; };

template <typename Matrix, mat_structure::tag Structure = mat_traits<Matrix>::structure >
class mat_sub_sym_block { char invalid_matrix_type[0]; };

template <typename Matrix, mat_structure::tag Structure = mat_traits<Matrix>::structure >
class mat_const_sub_sym_block { char invalid_matrix_type[0]; };


template <typename Matrix>
class mat_copy_sub_sym_block<Matrix,mat_structure::symmetric> {
  public:
    typedef mat_sub_sym_block<Matrix,mat_structure::symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::symmetric);
    
  private:
    Matrix m;
    size_type offset;
    size_type rowCount;
  public:
    mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) { };
    
    explicit mat_copy_sub_sym_block(const Matrix& aM) : m(aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_copy_sub_sym_block(const Matrix& aM, 
		      size_type aSize,
		      size_type aOffset = 0) : m(aM), offset(aOffset), rowCount(aSize) { };
		      
    mat_copy_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    explicit mat_copy_sub_sym_block(Matrix&& aM) : m(std::move(aM)), offset(0), rowCount(0) { rowCount = m.get_row_count(); };
    
    mat_copy_sub_sym_block(Matrix&& aM, 
		           size_type aSize,
		           size_type aOffset = 0) : m(std::move(aM)), offset(aOffset), rowCount(aSize) { };
		      
    mat_copy_sub_sym_block(self&& aObj) : m(std::move(aObj.m)), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
       
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i,offset + j) = (rhs(i,j) + rhs(j,i)) * value_type(0.5);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return m(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return m(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) += (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) -= (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::symmetric> operator -() const {
      mat<value_type,mat_structure::symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = j; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend self&& transpose(self&& M) {
      return std::move(M);
    };
#endif
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};


template <typename Matrix>
class mat_sub_sym_block<Matrix,mat_structure::symmetric> {
  public:
    typedef mat_sub_sym_block<Matrix,mat_structure::symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::symmetric);
    
  private:
    Matrix* m;
    size_type offset;
    size_type rowCount;
  public:
    explicit mat_sub_sym_block(Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_sub_sym_block(Matrix& aM, 
		      size_type aSize,
		      size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
    mat_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
   
    mat_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  (*m)(offset + i,offset + j) = (rhs(i,j) + rhs(j,i)) * value_type(0.5);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*m)(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  (*m)(offset + i, offset + j) += (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  (*m)(offset + i, offset + j) -= (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  (*m)(offset + i, offset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::symmetric> operator -() const {
      mat<value_type,mat_structure::symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = j; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};



template <typename Matrix>
class mat_const_sub_sym_block<Matrix,mat_structure::symmetric> {
  public:
    typedef mat_const_sub_sym_block<Matrix,mat_structure::symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::symmetric);
    
  private:
    const Matrix* m;
    size_type offset;
    size_type rowCount;
    
    self& operator=(const self&);
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(Matrix&&);
    
    mat_const_sub_sym_block(Matrix&&, size_type, size_type aOffset = 0);
#endif
  public:
    
    mat_const_sub_sym_block(const Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_const_sub_sym_block(const Matrix& aM, 
			    size_type aSize,
			    size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
			    
    mat_const_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    

#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> operator -() const {
      mat<value_type,mat_structure::symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = j; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend self&& transpose(self&& M) {
      return std::move(M);
    };
#endif
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};








template <typename Matrix>
class mat_copy_sub_sym_block<Matrix,mat_structure::skew_symmetric> {
  public:
    typedef mat_sub_sym_block<Matrix,mat_structure::skew_symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::skew_symmetric);
    
  private:
    Matrix m;
    size_type offset;
    size_type rowCount;
  public:
    mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) { };
    
    explicit mat_copy_sub_sym_block(const Matrix& aM) : m(aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_copy_sub_sym_block(const Matrix& aM, 
		      size_type aSize,
		      size_type aOffset = 0) : m(aM), offset(aOffset), rowCount(aSize) { };
		      
    mat_copy_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    explicit mat_copy_sub_sym_block(Matrix&& aM) : m(std::move(aM)), offset(0), rowCount(0) { rowCount = m.get_row_count(); };
    
    mat_copy_sub_sym_block(Matrix&& aM, 
		           size_type aSize,
		           size_type aOffset = 0) : m(std::move(aM)), offset(aOffset), rowCount(aSize) { };
		      
    mat_copy_sub_sym_block(self&& aObj) : m(std::move(aObj.m)), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
       
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i,offset + j) = (rhs(i,j) + rhs(j,i)) * value_type(0.5);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return m(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return m(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) += (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) -= (M(i,j) + M(j,i)) * value_type(0.5);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=j;i<rowCount;++i)
	  m(offset + i, offset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::skew_symmetric> operator -() const {
      mat<value_type,mat_structure::skew_symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = j; i < result.get_row_count(); ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose(const self& M) {
      return -M;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose_move(const self& M) {
      return -M;
    };
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};


template <typename Matrix>
class mat_sub_sym_block<Matrix,mat_structure::skew_symmetric> {
  public:
    typedef mat_sub_sym_block<Matrix,mat_structure::skew_symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::skew_symmetric);
    
  private:
    Matrix* m;
    size_type offset;
    size_type rowCount;
  public:
    
    mat_sub_sym_block(Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_sub_sym_block(Matrix& aM, 
		      size_type aSize,
		      size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
		      
    mat_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type j=1;j<rowCount;++j)
        for(size_type i=0;i<j;++i)
	  (*m)(offset + i,offset + j) = (rhs(i,j) - rhs(j,i)) * value_type(0.5);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*m)(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=0;i<j;++i)
	  (*m)(offset + i, offset + j) += (M(i,j) - M(j,i)) * value_type(0.5);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=0;i<j;++i)
	  (*m)(offset + i, offset + j) -= (M(i,j) - M(j,i)) * value_type(0.5);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type j=0;j<rowCount;++j)
        for(size_type i=0;i<j;++i)
	  (*m)(offset + i, offset + j) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::skew_symmetric> operator -() const {
      mat<value_type,mat_structure::skew_symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < j; ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose(const self& M) {
      return -M;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose_move(const self& M) {
      return -M;
    };
    
    friend value_type trace(const self& M) {
      return value_type(0.0);
    };
    
};




template <typename Matrix>
class mat_const_sub_sym_block<Matrix,mat_structure::skew_symmetric> {
  public:
    typedef mat_const_sub_sym_block<Matrix,mat_structure::skew_symmetric> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::skew_symmetric);
    
  private:
    const Matrix* m;
    size_type offset;
    size_type rowCount;
    
    self& operator=(const self&);
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(Matrix&& aM);
    
    mat_const_sub_sym_block(Matrix&& aM, size_type aSize, size_type aOffset = 0);
#endif
  public:
    
    mat_const_sub_sym_block(const Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_const_sub_sym_block(const Matrix& aM, 
			    size_type aSize,
			    size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
    
    mat_const_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> operator -() const {
      mat<value_type,mat_structure::skew_symmetric> result(*this);
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < j; ++i)
	  result(i,j) = -result(i,j);
      return result;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose(const self& M) {
      return -M;
    };
    
    friend mat<value_type,mat_structure::skew_symmetric> transpose_move(const self& M) {
      return -M;
    };
    
    friend value_type trace(const self& M) {
      return value_type(0.0);
    };
    
};








template <typename Matrix>
class mat_copy_sub_sym_block<Matrix,mat_structure::diagonal> {
  public:
    typedef mat_copy_sub_sym_block<Matrix,mat_structure::diagonal> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::diagonal);
    
  private:
    Matrix m;
    size_type offset;
    size_type rowCount;
  public:
    mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) { rowCount = m.get_row_count(); };
    
    mat_copy_sub_sym_block(const Matrix& aM) : m(aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_copy_sub_sym_block(const Matrix& aM, 
		           size_type aSize,
		           size_type aOffset = 0) : m(aM), offset(aOffset), rowCount(aSize) { };
			   
    mat_copy_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_copy_sub_sym_block(Matrix&& aM) : m(std::move(aM)), offset(0), rowCount(0) { rowCount = m.get_row_count(); };
    
    mat_copy_sub_sym_block(Matrix&& aM, 
		           size_type aSize,
		           size_type aOffset = 0) : m(std::move(aM)), offset(aOffset), rowCount(aSize) { };
			   
    mat_copy_sub_sym_block(self&& aObj) : m(std::move(aObj.m)), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
       
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type i=0;i<rowCount;++i)
        m(offset + i,offset + i) = rhs(i,i);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return m(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return m(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m.get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type i=0;i<rowCount;++i)
        m(offset + i, offset + i) += M(i,i);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type i=0;i<rowCount;++i)
	m(offset + i, offset + i) -= M(i,i);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type i=0;i<rowCount;++i)
	m(offset + i, offset + i) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::diagonal> operator -() const {
      mat<value_type,mat_structure::diagonal> result(*this);
      for(size_type i = 0; i < result.get_col_count(); ++i)
	result(i,i) = -result(i,i);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend self&& transpose(self&& M) {
      return std::move(M);
    };
#endif
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};




template <typename Matrix>
class mat_sub_sym_block<Matrix,mat_structure::diagonal> {
  public:
    typedef mat_sub_sym_block<Matrix,mat_structure::diagonal> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::diagonal);
    
  private:
    Matrix* m;
    size_type offset;
    size_type rowCount;
  public:
    mat_sub_sym_block(Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_sub_sym_block(Matrix& aM, 
		      size_type aSize,
		      size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
    mat_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
       
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value,
    self& >::type operator=(const Matrix2& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != rowCount))
	throw std::range_error("Matrix dimensions mismatch.");
      for(size_type i=0;i<rowCount;++i)
        (*m)(offset + i,offset + i) = rhs(i,i);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*m)(offset + i, offset + j);
    };
    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator +=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type i=0;i<rowCount;++i)
        (*m)(offset + i, offset + i) += M(i,i);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix2>
    self& operator -=(const Matrix2& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix2> >();
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      for(size_type i=0;i<rowCount;++i)
	(*m)(offset + i, offset + i) -= M(i,i);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type i=0;i<rowCount;++i)
	(*m)(offset + i, offset + i) *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix2>
    typename boost::enable_if_c< is_readable_matrix<Matrix2>::value, 
     self&>::type operator *=(const Matrix2& M) {
      if((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount))
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
    mat<value_type,mat_structure::diagonal> operator -() const {
      mat<value_type,mat_structure::diagonal> result(*this);
      for(size_type i = 0; i < result.get_col_count(); ++i)
	result(i,i) = -result(i,i);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend self&& transpose(self&& M) {
      return std::move(M);
    };
#endif
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};




template <typename Matrix>
class mat_const_sub_sym_block<Matrix,mat_structure::diagonal> {
  public:
    typedef mat_const_sub_sym_block<Matrix,mat_structure::diagonal> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator col_iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_col_iterator;
    typedef typename mat_traits<Matrix>::row_iterator row_iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_row_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_traits<Matrix>::alignment);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::diagonal);
    
  private:
    const Matrix* m;
    size_type offset;
    size_type rowCount;
    
    self& operator=(const self&);
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(Matrix&&);
    mat_const_sub_sym_block(Matrix&&, size_type, size_type aOffset = 0);
#endif
  public:
    mat_const_sub_sym_block(const Matrix& aM) : m(&aM), offset(0), rowCount(aM.get_row_count()) { };
    
    mat_const_sub_sym_block(const Matrix& aM, 
			    size_type aSize,
			    size_type aOffset = 0) : m(&aM), offset(aOffset), rowCount(aSize) { };
    
    mat_const_sub_sym_block(const self& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
    
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_sub_sym_block(self&& aObj) : m(aObj.m), offset(aObj.offset), rowCount(aObj.rowCount) { };
#endif
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const { 
      return (*m)(offset + i, offset + j); 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return rowCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,rowCount); };

    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::diagonal> operator -() const {
      mat<value_type,mat_structure::diagonal> result(*this);
      for(size_type i = 0; i < result.get_col_count(); ++i)
	result(i,i) = -result(i,i);
      return result;
    };
    
    friend const self& transpose(const self& M) {
      return M;
    };
    
    friend const self& transpose_move(const self& M) {
      return M;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend self&& transpose(self&& M) {
      return std::move(M);
    };
#endif
    
    friend value_type trace(const self& M) {
      value_type result(0.0);
      for(size_type i = 0; i < M.rowCount; ++i)
	result += M(i,i);
      return result;
    };
    
};







template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<Matrix>::value );
  typedef is_writable_matrix<Matrix> type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_resizable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_copy_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct has_allocator_matrix< mat_copy_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};



template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix< mat_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix< mat_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix<Matrix>::value );
  typedef is_writable_matrix<Matrix> type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix< mat_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_resizable_matrix< mat_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct has_allocator_matrix< mat_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};



template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix< mat_const_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix<Matrix>::value );
  typedef is_readable_matrix<Matrix> type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix< mat_const_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix< mat_const_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_resizable_matrix< mat_const_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_sub_sym_block<Matrix, Structure> > type;
};

template <typename Matrix, mat_structure::tag Structure>
struct has_allocator_matrix< mat_const_sub_sym_block<Matrix, Structure> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_matrix<Matrix> type;
};






template <typename Matrix>
struct mat_copy_sub_sym_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;

#ifndef RK_ENABLE_CXX0X_FEATURES
  Matrix m;
  mat_copy_sub_sym_block_factory(const Matrix& aM) : m(aM) { };
  mat_copy_sub_sym_block<Matrix> operator()(const std::pair<size_type,size_type>& rows) {
    return mat_copy_sub_sym_block<Matrix>(m,rows.second - rows.first + 1,rows.first);
  };
#else
  Matrix m;
  mat_copy_sub_sym_block_factory(Matrix&& aM) : m(std::move(aM)) { };
  mat_copy_sub_sym_block<Matrix> operator()(const std::pair<size_type,size_type>& rows) {
    return mat_copy_sub_sym_block<Matrix>(std::move(m),rows.second - rows.first + 1,rows.first);
  };
#endif
};

template <typename Matrix>
struct mat_sub_sym_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  Matrix& m;
  mat_sub_sym_block_factory(Matrix& aM) : m(aM) { };
  mat_sub_sym_block<Matrix> operator()(const std::pair<size_type,size_type>& rows) {
    return mat_sub_sym_block<Matrix>(m,rows.second - rows.first + 1,rows.first);
  };
};

template <typename Matrix>
struct mat_const_sub_sym_block_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  const Matrix& m;
  mat_const_sub_sym_block_factory(const Matrix& aM) : m(aM) { };
  mat_const_sub_sym_block<Matrix> operator()(const std::pair<size_type,size_type>& rows) {
    return mat_const_sub_sym_block<Matrix>(m,rows.second - rows.first + 1,rows.first);
  };
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_sub_sym_block_factory<Matrix> >::type sub_sym(Matrix& M) { return mat_sub_sym_block_factory<Matrix>(M); };

template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_const_sub_sym_block_factory<Matrix> >::type sub_sym(const Matrix& M) { return mat_const_sub_sym_block_factory<Matrix>(M); };

template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_copy_sub_sym_block_factory<Matrix> >::type sub_sym_copy(const Matrix& M) { return mat_copy_sub_sym_block_factory<Matrix>(M); };

#ifdef RK_ENABLE_CXX0X_FEATURES
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_copy_sub_sym_block_factory<Matrix> >::type sub_sym(Matrix&& M) { return mat_copy_sub_sym_block_factory<Matrix>(std::move(M)); };
#endif











};

#endif



















