
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
    LeftMatrix ml;
    RightMatrix mr;
  public:
    mat_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR) : ml(aML), mr(aMR) { 
      if(ml.get_row_count() != mr.get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.ml,rhs.ml);
      swap(lhs.mr,rhs.mr);
    };
    
    self& operator=(const self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != ml.get_row_count()) || (rhs.get_col_count() != ml.get_col_count() + mr.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      ml = mat_const_sub_block< Matrix >(rhs,ml.get_row_count(),ml.get_col_count(),0,0);
      mr = mat_const_sub_block< Matrix >(rhs,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };
    value_type operator()(size_type i,size_type j) const { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };

    size_type get_row_count() const throw() { return ml.get_row_count(); };
    size_type get_col_count() const throw() { return ml.get_col_count() + mr.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml.get_row_count(),ml.get_col_count() + mr.get_col_count()); };

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
      ml += mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),0,0);
      mr += mat_const_sub_block< Matrix >(M,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
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
      ml -= mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),0,0);
      mr -= mat_const_sub_block< Matrix >(M,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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
    LeftMatrix& ml;
    RightMatrix& mr;
  public:
    mat_ref_horiz_cat(LeftMatrix& aML, RightMatrix& aMR) : ml(aML), mr(aMR) { 
      if(ml.get_row_count() != mr.get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_ref_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != ml.get_row_count()) || (rhs.get_col_count() != ml.get_col_count() + mr.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      ml = mat_const_sub_block< Matrix >(rhs,ml.get_row_count(),ml.get_col_count(),0,0);
      mr = mat_const_sub_block< Matrix >(rhs,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };
    value_type operator()(size_type i,size_type j) const { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };

    size_type get_row_count() const throw() { return ml.get_row_count(); };
    size_type get_col_count() const throw() { return ml.get_col_count() + mr.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml.get_row_count(),ml.get_col_count() + mr.get_col_count()); };

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
      ml += mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),0,0);
      mr += mat_const_sub_block< Matrix >(M,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
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
      ml -= mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),0,0);
      mr -= mat_const_sub_block< Matrix >(M,mr.get_row_count(),mr.get_col_count(),0,ml.get_col_count());
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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
    const LeftMatrix& ml;
    const RightMatrix& mr;
    
    self& operator=(const self&);
  public:
    mat_const_ref_horiz_cat(const LeftMatrix& aML, const RightMatrix& aMR) : ml(aML), mr(aMR) { 
      if(ml.get_row_count() != mr.get_row_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_const_ref_horiz_cat(const self& aObj) : ml(aObj.ml), mr(aObj.mr) { };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const { 
      if(j < ml.get_col_count())
	return ml(i,j);
      else
	return mr(i,j - ml.get_col_count());
    };

    size_type get_row_count() const throw() { return ml.get_row_count(); };
    size_type get_col_count() const throw() { return ml.get_col_count() + mr.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(ml.get_row_count(),ml.get_col_count() + mr.get_col_count()); };

    allocator_type get_allocator() const { return ml.get_allocator(); };
    
    
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
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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
    UpperMatrix mu;
    LowerMatrix ml;
  public:
    mat_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML) : mu(aMU), ml(aML) { 
      if(ml.get_col_count() != mu.get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mu,rhs.mu);
      swap(lhs.ml,rhs.ml);
    };
    
    self& operator=(const self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != mu.get_row_count() + ml.get_row_count()) || (rhs.get_col_count() != mu.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      mu = mat_const_sub_block< Matrix >(rhs,mu.get_row_count(),mu.get_col_count(),0,0);
      ml = mat_const_sub_block< Matrix >(rhs,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };
    value_type operator()(size_type i,size_type j) const {  
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };

    size_type get_row_count() const throw() { return mu.get_row_count() + ml.get_row_count(); };
    size_type get_col_count() const throw() { return mu.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu.get_row_count() + ml.get_row_count(),mu.get_col_count()); };

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
      mu += mat_const_sub_block< Matrix >(M,mu.get_row_count(),mu.get_col_count(),0,0);
      ml += mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
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
      mu -= mat_const_sub_block< Matrix >(M,mu.get_row_count(),mu.get_col_count(),0,0);
      ml -= mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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
    UpperMatrix& mu;
    LowerMatrix& ml;
  public:
    mat_ref_vert_cat(UpperMatrix& aMU, LowerMatrix& aML) : mu(aMU), ml(aML) { 
      if(ml.get_col_count() != mu.get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_ref_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    self& >::type operator=(const Matrix& rhs) {
      if((rhs.get_row_count() != mu.get_row_count() + ml.get_row_count()) || (rhs.get_col_count() != mu.get_col_count()))
	throw std::range_error("Matrix dimensions mismatch.");
      mu = mat_const_sub_block< Matrix >(rhs,mu.get_row_count(),mu.get_col_count(),0,0);
      ml = mat_const_sub_block< Matrix >(rhs,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
      return *this;
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { 
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };
    value_type operator()(size_type i,size_type j) const {  
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };

    size_type get_row_count() const throw() { return mu.get_row_count() + ml.get_row_count(); };
    size_type get_col_count() const throw() { return mu.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu.get_row_count() + ml.get_row_count(),mu.get_col_count()); };

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
      mu += mat_const_sub_block< Matrix >(M,mu.get_row_count(),mu.get_col_count(),0,0);
      ml += mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
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
      mu -= mat_const_sub_block< Matrix >(M,mu.get_row_count(),mu.get_col_count(),0,0);
      ml -= mat_const_sub_block< Matrix >(M,ml.get_row_count(),ml.get_col_count(),mu.get_row_count(),0);
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
    
    friend mat<value_type,mat_structure::rectangular> transpose(const self& M) {
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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
    const UpperMatrix& mu;
    const LowerMatrix& ml;
    
    self& operator=(const self&);
  public:
    mat_const_ref_vert_cat(const UpperMatrix& aMU, const LowerMatrix& aML) : mu(aMU), ml(aML) { 
      if(ml.get_col_count() != mu.get_col_count())
	throw std::range_error("Matrix dimensions mismatch.");
    };
    
    mat_const_ref_vert_cat(const self& aObj) : mu(aObj.mu), ml(aObj.ml) { };
    
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    value_type operator()(size_type i,size_type j) const {  
      if(i < mu.get_row_count())
	return mu(i,j);
      else
	return ml(i - mu.get_row_count(),j);
    };

    size_type get_row_count() const throw() { return mu.get_row_count() + ml.get_row_count(); };
    size_type get_col_count() const throw() { return mu.get_col_count(); };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(mu.get_row_count() + ml.get_row_count(),mu.get_col_count()); };

    allocator_type get_allocator() const { return mu.get_allocator(); };
    
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
      mat<value_type,mat_structure::rectangular> result(M.get_col_count(), M.get_row_count());
      for(size_type j = 0; j < result.get_col_count(); ++j)
	for(size_type i = 0; i < result.get_row_count(); ++i)
	  result(i,j) = M(j,i);
      return result;
    };
    
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








};

#endif






















