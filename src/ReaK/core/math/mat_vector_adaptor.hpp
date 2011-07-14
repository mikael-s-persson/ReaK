
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

#ifndef MAT_VECTOR_ADAPTOR_HPP
#define MAT_VECTOR_ADAPTOR_HPP

#include "mat_alg_general.hpp"

namespace ReaK {


template <typename Vector, mat_alignment::tag Alignment = mat_alignment::column_major>
class mat_vect_adaptor { };

template <typename Vector>
class mat_vect_adaptor<Vector,mat_alignment::column_major> {
  public:    
    
    typedef mat_vect_adaptor<Vector,mat_alignment::column_major> self;
    typedef typename vect_traits<Vector>::allocator_type allocator_type;
    
    typedef typename vect_traits<Vector>::value_type value_type;
    
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::const_reference const_reference;
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef typename vect_traits<Vector>::const_pointer const_pointer;
  
    typedef typename vect_traits<Vector>::iterator col_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_col_iterator;
    typedef typename vect_traits<Vector>::iterator row_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_row_iterator;
  
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    Vector* v;
    size_type offset;
    size_type rowCount;
    size_type colCount;
  public:
    mat_vect_adaptor(Vector& aV) : v(&aV), offset(0), rowCount(aV.size()), colCount(1) { };
    
    mat_vect_adaptor(Vector& aV, 
		     size_type aRowCount, 
		     size_type aColCount,
		     size_type aOffset = 0) : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) { };
    mat_vect_adaptor(const self& aObj) : v(aObj.v), offset(aObj.offset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
    friend void swap(self& lhs, self& rhs) {
      using std::swap;
      swap(lhs.v,rhs.v);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      size_type it = offset;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  (*v)[it] = rhs(i,j);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*v)[offset + j*rowCount + i]; 
    };
    const_reference operator()(size_type i,size_type j) const { 
      return (*v)[offset + j*rowCount + i]; 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return v->get_allocator(); };
    
    
        /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      size_type it = offset;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  (*v)[it] += M(i,j);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      size_type it = offset;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  (*v)[it] -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type it = offset;it < offset + rowCount*colCount;++it)
        (*v)[it] *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
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
    
    friend mat_vect_adaptor<Vector,mat_alignment::row_major> transpose(const self& M) {
      return mat_vect_adaptor<Vector,mat_alignment::row_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
    friend mat_vect_adaptor<Vector,mat_alignment::row_major> transpose_move(self& M) {
      return mat_vect_adaptor<Vector,mat_alignment::row_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
};




template <typename Vector>
class mat_vect_adaptor<Vector,mat_alignment::row_major> {
  public:    
    
    typedef mat_vect_adaptor<Vector,mat_alignment::row_major> self;
    typedef typename vect_traits<Vector>::allocator_type allocator_type;
    
    typedef typename vect_traits<Vector>::value_type value_type;
    
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::const_reference const_reference;
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef typename vect_traits<Vector>::const_pointer const_pointer;
  
    typedef typename vect_traits<Vector>::iterator col_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_col_iterator;
    typedef typename vect_traits<Vector>::iterator row_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_row_iterator;
  
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::row_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    Vector* v;
    size_type offset;
    size_type rowCount;
    size_type colCount;
  public:
    mat_vect_adaptor(Vector& aV) : v(&aV), offset(0), rowCount(1), colCount(aV.size()) { };
    
    mat_vect_adaptor(Vector& aV, 
		     size_type aRowCount, 
		     size_type aColCount,
		     size_type aOffset = 0) : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) { };
    mat_vect_adaptor(const self& aObj) : v(aObj.v), offset(aObj.offset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
    friend void swap(self& lhs, self& rhs) {
      using std::swap;
      swap(lhs.v,rhs.v);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != rowCount) || (rhs.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      size_type it = offset;
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  (*v)[it] = rhs(i,j);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      return (*v)[offset + i*colCount + j]; 
    };
    const_reference operator()(size_type i,size_type j) const { 
      return (*v)[offset + i*colCount + j]; 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return v->get_allocator(); };
    
    
        /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      size_type it = offset;
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  (*v)[it] += M(i,j);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      size_type it = offset;
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  (*v)[it] -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(size_type it = offset;it < offset + rowCount*colCount;++it)
        (*v)[it] *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
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
    
    friend mat_vect_adaptor<Vector,mat_alignment::column_major> transpose(const self& M) {
      return mat_vect_adaptor<Vector,mat_alignment::column_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
    friend mat_vect_adaptor<Vector,mat_alignment::column_major> transpose_move(self& M) {
      return mat_vect_adaptor<Vector,mat_alignment::column_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
};





template <typename Vector, mat_alignment::tag Alignment>
struct is_readable_matrix< mat_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_vector< Vector >::value );
  typedef is_readable_vector< Vector > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_writable_matrix< mat_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_vector<Vector>::value );
  typedef is_writable_vector<Vector> type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_fully_writable_matrix< mat_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_vector<Vector>::value );
  typedef is_writable_vector<Vector> type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_resizable_matrix< mat_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_vect_adaptor<Vector,Alignment> > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct has_allocator_matrix< mat_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector<Vector>::value );
  typedef has_allocator_vector<Vector> type;
};








template <typename Vector, mat_alignment::tag Alignment = mat_alignment::column_major>
class mat_const_vect_adaptor { };

template <typename Vector>
class mat_const_vect_adaptor<Vector,mat_alignment::column_major> {
  public:    
    
    typedef mat_const_vect_adaptor<Vector,mat_alignment::column_major> self;
    typedef typename vect_traits<Vector>::allocator_type allocator_type;
    
    typedef typename vect_traits<Vector>::value_type value_type;
    
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::const_reference const_reference;
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef typename vect_traits<Vector>::const_pointer const_pointer;
  
    typedef typename vect_traits<Vector>::iterator col_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_col_iterator;
    typedef typename vect_traits<Vector>::iterator row_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_row_iterator;
  
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    const Vector* v;
    size_type offset;
    size_type rowCount;
    size_type colCount;
    
    self& operator=(const self&); //non-assignable.
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_vect_adaptor(Vector&&);
    mat_const_vect_adaptor(Vector&&, size_type, size_type, size_type aOffset = 0);
#endif
  public:
    mat_const_vect_adaptor(const Vector& aV) : v(&aV), offset(0), rowCount(aV.size()), colCount(1) { };
    
    mat_const_vect_adaptor(const Vector& aV, 
		     size_type aRowCount, 
		     size_type aColCount,
		     size_type aOffset = 0) : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) { };
    mat_const_vect_adaptor(const self& aObj) : v(aObj.v), offset(aObj.offset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
    friend void swap(self& lhs, self& rhs) {
      using std::swap;
      swap(lhs.v,rhs.v);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    const_reference operator()(size_type i,size_type j) const { 
      return (*v)[offset + j*rowCount + i]; 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return v->get_allocator(); };
    
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
    
    friend mat_const_vect_adaptor<Vector,mat_alignment::row_major> transpose(const self& M) {
      return mat_const_vect_adaptor<Vector,mat_alignment::row_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
    friend mat_const_vect_adaptor<Vector,mat_alignment::row_major> transpose_move(self& M) {
      return mat_const_vect_adaptor<Vector,mat_alignment::row_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
};




template <typename Vector>
class mat_const_vect_adaptor<Vector,mat_alignment::row_major> {
  public:    
    
    typedef mat_const_vect_adaptor<Vector,mat_alignment::row_major> self;
    typedef typename vect_traits<Vector>::allocator_type allocator_type;
    
    typedef typename vect_traits<Vector>::value_type value_type;
    
    typedef typename vect_traits<Vector>::reference reference;
    typedef typename vect_traits<Vector>::const_reference const_reference;
    typedef typename vect_traits<Vector>::pointer pointer;
    typedef typename vect_traits<Vector>::const_pointer const_pointer;
  
    typedef typename vect_traits<Vector>::iterator col_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_col_iterator;
    typedef typename vect_traits<Vector>::iterator row_iterator;
    typedef typename vect_traits<Vector>::const_iterator const_row_iterator;
  
    typedef typename vect_traits<Vector>::size_type size_type;
    typedef typename vect_traits<Vector>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::row_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  private:
    const Vector* v;
    size_type offset;
    size_type rowCount;
    size_type colCount;
    
    self& operator=(const self&); //non-assignable.
    
        
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_vect_adaptor(Vector&&);
    mat_const_vect_adaptor(Vector&&, size_type, size_type, size_type aOffset = 0);
#endif
  public:
    mat_const_vect_adaptor(const Vector& aV) : v(&aV), offset(0), rowCount(1), colCount(aV.size()) { };
    
    mat_const_vect_adaptor(const Vector& aV, 
			   size_type aRowCount, 
			   size_type aColCount,
			   size_type aOffset = 0) : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) { };
    mat_const_vect_adaptor(const self& aObj) : v(aObj.v), offset(aObj.offset), rowCount(aObj.rowCount), colCount(aObj.colCount) { };
    
    friend void swap(self& lhs, self& rhs) {
      using std::swap;
      swap(lhs.v,rhs.v);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowCount,rhs.rowCount);
      swap(lhs.colCount,rhs.colCount);
      return;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    const_reference operator()(size_type i,size_type j) const { 
      return (*v)[offset + i*colCount + j]; 
    };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };

    allocator_type get_allocator() const { return v->get_allocator(); };
    
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
    
    friend mat_const_vect_adaptor<Vector,mat_alignment::column_major> transpose(const self& M) {
      return mat_const_vect_adaptor<Vector,mat_alignment::column_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
    friend mat_const_vect_adaptor<Vector,mat_alignment::column_major> transpose_move(self& M) {
      return mat_const_vect_adaptor<Vector,mat_alignment::column_major>(*M.v,M.offset,M.colCount,M.rowCount);
    };
    
};





template <typename Vector, mat_alignment::tag Alignment>
struct is_readable_matrix< mat_const_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_vector< Vector >::value );
  typedef is_readable_vector< Vector > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_writable_matrix< mat_const_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_vect_adaptor<Vector,Alignment> > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_fully_writable_matrix< mat_const_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_vect_adaptor<Vector,Alignment> > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct is_resizable_matrix< mat_const_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_vect_adaptor<Vector,Alignment> > type;
};

template <typename Vector, mat_alignment::tag Alignment>
struct has_allocator_matrix< mat_const_vect_adaptor<Vector,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector<Vector>::value );
  typedef has_allocator_vector<Vector> type;
};


template <typename Vector, mat_alignment::tag Alignment, typename Matrix2>
struct mat_product_result<mat_const_vect_adaptor<Vector,Alignment>,Matrix2> {
  typedef typename vect_traits<Vector>::value_type value_type;
  typedef typename mat_product_result<mat<value_type,mat_structure::rectangular>,Matrix2>::type type;
};

template <typename Vector, mat_alignment::tag Alignment, typename Matrix2>
struct mat_addition_result<mat_const_vect_adaptor<Vector,Alignment>,Matrix2> {
  typedef typename vect_traits<Vector>::value_type value_type;
  typedef typename mat_addition_result<mat<value_type,mat_structure::rectangular>,Matrix2>::type type;
};

template <typename Vector, mat_alignment::tag Alignment, typename Matrix1>
struct mat_product_result< Matrix1,mat_const_vect_adaptor<Vector,Alignment> > {
  typedef typename vect_traits<Vector>::value_type value_type;
  typedef typename mat_product_result<Matrix1,mat<value_type,mat_structure::rectangular> >::type type;
};

template <typename Vector, mat_alignment::tag Alignment, typename Matrix1>
struct mat_addition_result< Matrix1,mat_const_vect_adaptor<Vector,Alignment> > {
  typedef typename vect_traits<Vector>::value_type value_type;
  typedef typename mat_addition_result<Matrix1,mat<value_type,mat_structure::rectangular> >::type type;
};




template <typename Vector>
struct mat_vect_adaptor_factory {
  typedef typename vect_traits<Vector>::size_type size_type;
  
  Vector& v;
  mat_vect_adaptor_factory(Vector& aV) : v(aV) { };
  mat_vect_adaptor<Vector,mat_alignment::row_major> operator()(size_type rowCount,
				      const std::pair<size_type,size_type>& cols) {
    return mat_vect_adaptor<Vector,mat_alignment::row_major>(v,rowCount,cols.second - cols.first + 1,cols.first);
  };
  mat_vect_adaptor<Vector,mat_alignment::column_major> operator()(const std::pair<size_type,size_type>& rows,
				   size_type colCount) {
    return mat_vect_adaptor<Vector,mat_alignment::column_major>(v,rows.second - rows.first + 1,colCount,rows.first);
  };
};

template <typename Vector>
struct mat_const_vect_adaptor_factory {
  typedef typename vect_traits<Vector>::size_type size_type;
  
  const Vector& v;
  mat_const_vect_adaptor_factory(const Vector& aV) : v(aV) { };
  mat_const_vect_adaptor<Vector,mat_alignment::row_major> operator()(size_type rowCount,
				   const std::pair<size_type,size_type>& cols) {
    return mat_const_vect_adaptor<Vector,mat_alignment::row_major>(v,rowCount,cols.second - cols.first + 1,cols.first);
  };
  mat_const_vect_adaptor<Vector,mat_alignment::column_major> operator()(const std::pair<size_type,size_type>& rows,
				   size_type colCount) {
    return mat_const_vect_adaptor<Vector,mat_alignment::column_major>(v,rows.second - rows.first + 1,colCount,rows.first);
  };
};


template <typename Vector>
typename boost::enable_if_c< is_readable_vector<Vector>::value,
mat_vect_adaptor_factory<Vector> >::type make_mat(Vector& V) { return mat_vect_adaptor_factory<Vector>(V); };

template <typename Vector>
typename boost::enable_if_c< is_readable_vector<Vector>::value,
mat_const_vect_adaptor_factory<Vector> >::type make_mat(const Vector& V) { return mat_const_vect_adaptor_factory<Vector>(V); };











};

#endif













