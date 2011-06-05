/**
 * \file mat_alg_rectangular_fixed.hpp
 * 
 * This library implements the specialization of the mat<> template for a 
 * general rectangular matrix (fixed dimensions) of both column-major and 
 * row-major alignment. This matrix type fulfills all the general matrix 
 * concepts (Readable, Writable, and Fully-Writable).
 * 
 * This library also implements transposition of matrices via alignment 
 * switching (switching from column-major to row-major, or vice versa). This 
 * is very efficient and can even avoid copies completely (via transpose_move()) on an
 * optimizing compiler.
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

#ifndef MAT_ALG_RECTANGULAR_FIXED_HPP
#define MAT_ALG_RECTANGULAR_FIXED_HPP

#include "mat_alg_general.hpp"

namespace ReaK {
  
  

template <typename T,
          unsigned int RowCount, unsigned int ColCount,
          mat_alignment::tag Alignment>
struct is_fully_writable_matrix< mat_fixed<T,mat_structure::rectangular,RowCount,ColCount,Alignment> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_fully_writable_matrix< mat_fixed<T,mat_structure::rectangular,RowCount,ColCount,Alignment> > type;
};

  
/**
 * This class template specialization implements a matrix with rectangular structure
 * and column-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This matrix type is dynamically resizable.
 */
template <typename T,
          unsigned int rowCount,
	  unsigned int colCount>
class mat_fixed<T,mat_structure::rectangular,rowCount,colCount,mat_alignment::column_major> : public serialization::serializable {
  public:    
    
    typedef mat_fixed<T,mat_structure::rectangular,rowCount,colCount,mat_alignment::column_major> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef stride_iterator< T* > col_iterator;
    typedef stride_iterator< const T* > const_col_iterator;
    typedef T* row_iterator;
    typedef const T* const_row_iterator;
  
    typedef unsigned int size_type;
    typedef int difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = rowCount);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = colCount);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  
  private:
    value_type q[rowCount * colCount]; ///< Array which holds all the values of the matrix (dimension: rowCount x colCount).
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor.
     * \test PASSED
     */
    explicit mat_fixed(const value_type& aFill = value_type()) {  
      for(pointer it = q; it != q + rowCount * colCount; ++it)
	*it = aFill;
    };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat_fixed(const self& M) { 
      for(size_type i = 0; i < rowCount * colCount; ++i)
	q[i] = M.q[i];
    };

    /**
     * Explicit constructor from a any type of matrix.
     * \test PASSED
     */
    template <typename Matrix>
    explicit mat_fixed(const Matrix&  M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                                      !(boost::is_same<Matrix,self>::value) , void* >::type dummy = NULL) {
      if((M.get_row_count() != rowCount) || (M.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      pointer it = q;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  *it = M(i,j);
    };

    /**
     * Constructor from a vector of column major values.
     */
    explicit mat_fixed(const_pointer Q) {
      for(size_type i = 0; i < rowCount * colCount; ++i)
	q[i] = Q[i];
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~mat_fixed() { };

    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      for(size_type i = 0; i < rowCount * colCount; ++i)
	swap(m1.q[i],m2.q[i]);
    };
    
    friend void swap(self& m1, pointer q2) throw() {
      using std::swap;
      for(size_type i = 0; i < rowCount * colCount; ++i)
	swap(m1.q[i],q2[i]);
    };
    
    
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { return q[j*rowCount + i]; };
    const_reference operator()(size_type i,size_type j) const { return q[j*rowCount + i]; };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };
    void resize(const std::pair<size_type,size_type>& sz) { };

    void set_row_count(size_type,bool aPreserveData = false) { RK_UNUSED(aPreserveData); };

    void set_col_count(size_type,bool aPreserveData = false) { RK_UNUSED(aPreserveData); };
    
    row_iterator first_row() { 
      return q; 
    };
    const_row_iterator first_row() const { 
      return q;
    };
    row_iterator last_row() { 
      return q + rowCount; 
    };
    const_row_iterator last_row() const {
      return q + rowCount; 
    };
    row_iterator first_row(col_iterator cit) {
      size_type diff = cit.base() - q;
      return q + ((diff / rowCount) * rowCount);
    };
    const_row_iterator first_row(const_col_iterator cit) const {
      size_type diff = cit.base() - q;
      return q + ((diff / rowCount) * rowCount);
    };
    row_iterator last_row(col_iterator cit) {
      size_type diff = cit.base() - q;
      return q + ((diff / rowCount + 1) * rowCount);
    };
    const_row_iterator last_row(const_col_iterator cit) const {
      size_type diff = cit.base() - q;
      return q + ((diff / rowCount + 1) * rowCount);
    };
    std::pair<row_iterator,row_iterator> rows() { 
      return std::make_pair(first_row(),last_row());
    };
    std::pair<const_row_iterator,const_row_iterator> rows() const { 
      return std::make_pair(first_row(),last_row());
    };
    std::pair<row_iterator,row_iterator> rows(col_iterator cit) { 
      return std::make_pair(first_row(cit),last_row(cit));
    };
    std::pair<const_row_iterator,const_row_iterator> rows(const_col_iterator cit) const { 
      return std::make_pair(first_row(cit),last_row(cit));
    };
    
    col_iterator first_col() {
      return col_iterator(q,rowCount);
    };
    const_col_iterator first_col() const {
      return const_col_iterator(q,rowCount);
    };
    col_iterator last_col() {
      return col_iterator(q + colCount * rowCount,rowCount);
    };
    const_col_iterator last_col() const {
      return const_col_iterator(q + colCount * rowCount,rowCount);
    };
    col_iterator first_col(row_iterator rit) {
      return col_iterator(q + ((rit - q) % rowCount),rowCount);
    };
    const_col_iterator first_col(const_row_iterator rit) const {
      return const_col_iterator(q + ((rit - q) % rowCount),rowCount);
    };
    col_iterator last_col(row_iterator rit) {
      return col_iterator(q + ((rit - q) % rowCount) + colCount * rowCount,rowCount);
    };
    const_col_iterator last_col(const_row_iterator rit) const {
      return const_col_iterator(q + ((rit - q) % rowCount) + colCount * rowCount,rowCount);
    };
    std::pair<col_iterator,col_iterator> cols() { 
      return std::make_pair(first_col(),last_col());
    };
    std::pair<const_col_iterator,const_col_iterator> cols() const { 
      return std::make_pair(first_col(),last_col());
    };
    std::pair<col_iterator,col_iterator> cols(row_iterator rit) { 
      return std::make_pair(first_col(rit),last_col(rit));
    };
    std::pair<const_col_iterator,const_col_iterator> cols(const_row_iterator rit) const { 
      return std::make_pair(first_col(rit),last_col(rit));
    };
    

    allocator_type get_allocator() const { };
    
    

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** COL-MAJOR ONLY
     * Standard Assignment operator with standard semantics.
     * Strong exception-safety.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator =(const Matrix& M) {
      if((M.get_row_count() != rowCount) || (M.get_col_count() != colCount))
	throw std::range_error("Matrix dimensions mismatch.");
      pointer it = q;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  *it = M(i,j);
      return *this;
    };

    /** COL-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      pointer it = q;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  *it += M(i,j);
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
      pointer it = q;
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
	  *it -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(pointer it = q;it != q + rowCount * colCount;++it)
        *it *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
     self&>::type operator *=(const Matrix& M) {
      *this = *this * M;
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin();it != q.end();++it,++itr)
        *itr = -(*it);
      return result;
    };
    
    
    
    friend mat_fixed<T,mat_structure::rectangular,colCount,rowCount,mat_alignment::row_major> transpose(const self& M) {
      return mat_fixed<T,mat_structure::rectangular,colCount,rowCount,mat_alignment::row_major>(M.q);
    };
    
    friend mat_fixed<T,mat_structure::rectangular,colCount,rowCount,mat_alignment::row_major> transpose_move(self& M) {
      mat_fixed<T,mat_structure::rectangular,colCount,rowCount,mat_alignment::row_major> result;
      using std::swap;
      swap(result,M.q);
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      for(size_type i = 0; i < rowCount * colCount; ++i)
        A & RK_SERIAL_SAVE_WITH_NAME(q[i]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      for(size_type i = 0; i < rowCount * colCount; ++i)
        A & RK_SERIAL_LOAD_WITH_NAME(q[i]);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
    
  
};


/**
 * This class template specialization implements a matrix with rectangular structure
 * and row-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This matrix type is dynamically resizable.
 */
template <typename T,
	  typename Allocator>
class mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
  
    typedef stride_iterator< typename container_type::iterator > row_iterator;
    typedef stride_iterator< typename container_type::const_iterator > const_row_iterator;
    typedef typename container_type::iterator col_iterator;
    typedef typename container_type::const_iterator const_col_iterator;
  
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::row_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::rectangular);
    
  
  private:
    std::vector<value_type,allocator_type> q; ///< Array which holds all the values of the matrix (dimension: rowCount x colCount).
    size_type rowCount; ///< Row Count.
    size_type colCount; ///< Column Count.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor: sets all to zero.
     * \test PASSED
     */
    mat(const allocator_type& aAlloc = allocator_type()) :
             q(0,value_type(0),aAlloc),
	     rowCount(0),
	     colCount(0) { };

    /**
     * Constructor for a sized matrix.
     * \test PASSED
     */
    mat(size_type aRowCount,
        size_type aColCount, const value_type& aFill = value_type(),const allocator_type& aAlloc = allocator_type()) :
        q(aRowCount * aColCount,aFill,aAlloc),
	rowCount(aRowCount),
	colCount(aColCount) { };

    /**
     * Constructor for an identity matrix.
     * \test PASSED
     */
    mat(size_type aRowCount,size_type aColCount, bool aIdentity,const allocator_type& aAlloc = allocator_type()) :
             q(aRowCount * aColCount,0,aAlloc),
	     rowCount(aRowCount),
	     colCount(aColCount) {
      if(aIdentity) {
	size_type minN = (colCount < rowCount ? colCount : rowCount) * (colCount + 1);
	for(size_type i=0;i < minN;i += colCount + 1)
	  q[i] = 1;
      };
    };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(const self& M) :
             q(M.q),
	     rowCount(M.rowCount),
	     colCount(M.colCount) { };

	     
    /**
     * Explicit constructor from a any type of matrix.
     * \test PASSED
     */
    template <typename Matrix>
    explicit mat(const Matrix&  M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                               !(boost::is_same<Matrix,self>::value) , void* >::type dummy = NULL) :
             q(M.get_row_count()*M.get_col_count(),T(0.0)),
	     rowCount(M.get_row_count()),
	     colCount(M.get_col_count()) {
      typename container_type::iterator it = q.begin();
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  *it = M(i,j);
    };

    /**
     * Constructor from a vector of column major values.
     */
    mat(const container_type& Q,size_type aRowCount,size_type aColCount) :
             q(Q), rowCount(aRowCount), colCount(aColCount) { };

    /**
     * Destructor.
     * \test PASSED
     */
    ~mat() { };
    
    /**
     * Constructs a 2x2 matrix from four elements.
     * \test PASSED
     */
    mat(const value_type& a11,const value_type& a12,
	const value_type& a21,const value_type& a22) : q(4), rowCount(2), colCount(2) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a21;
      q[3] = a22;
    };

    /**
     * Constructs a 3x3 matrix from nine elements.
     * \test PASSED
     */
    mat(const value_type& a11,const value_type& a12,const value_type& a13,
	const value_type& a21,const value_type& a22,const value_type& a23,
	const value_type& a31,const value_type& a32,const value_type& a33) : 
	q(9), rowCount(3), colCount(3) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a13;
      q[3] = a21;
      q[4] = a22;
      q[5] = a23;
      q[6] = a31;
      q[7] = a32;
      q[8] = a33;
    };

    /**
     * Constructs a 4x4 matrix from sixteen elements.
     * \test PASSED
     */
    mat(const value_type& a11,const value_type& a12,const value_type& a13,const value_type& a14,
	const value_type& a21,const value_type& a22,const value_type& a23,const value_type& a24,
	const value_type& a31,const value_type& a32,const value_type& a33,const value_type& a34,
	const value_type& a41,const value_type& a42,const value_type& a43,const value_type& a44) : 
	q(16), rowCount(4), colCount(4) {
      q[0] = a11;
      q[1] = a12;
      q[2] = a13;
      q[3] = a14;
      q[4] = a21;
      q[5] = a22;
      q[6] = a23;
      q[7] = a24;
      q[8] = a31;
      q[9] = a32;
      q[10] = a33;
      q[11] = a34;
      q[12] = a41;
      q[13] = a42;
      q[14] = a43;
      q[15] = a44;
    };

    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.rowCount,m2.rowCount);
      swap(m1.colCount,m2.colCount);
    };
    
    friend void swap(self& m1, container_type& q2, size_type& rowCount2, size_type& colCount2) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.rowCount,rowCount2);
      swap(m1.colCount,colCount2);
    };
    
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    reference operator()(size_type i,size_type j) { return q[i*colCount + j]; };
    const_reference operator()(size_type i,size_type j) const { return q[i*colCount + j]; };

    size_type get_row_count() const throw() { return rowCount; };
    size_type get_col_count() const throw() { return colCount; };
    
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };
    void resize(const std::pair<size_type,size_type>& sz) {
      set_col_count(sz.second,true);
      set_row_count(sz.first,true);
    };

    void set_col_count(size_type aColCount,bool aPreserveData = false) {
      if(aPreserveData) {
	if(aColCount > colCount)
	  for(size_type i=rowCount;i!=0;--i)
	    q.insert(q.begin() + i*colCount,aColCount - colCount,value_type(0.0));
	else if(aColCount < colCount)
	  for(size_type i=rowCount;i!=0;--i)
	    q.erase(q.begin() + (i-1)*colCount + aColCount, q.begin() + i*colCount);
      } else
	q.resize(rowCount * aColCount,value_type(0.0));
      colCount = aColCount;
    };

    void set_row_count(size_type aRowCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(aRowCount * colCount,value_type(0.0));
      rowCount = aRowCount;
    };
    
    row_iterator first_row() { 
      return row_iterator(q.begin(),colCount); 
    };
    const_row_iterator first_row() const { 
      return const_row_iterator(q.begin(),colCount);
    };
    row_iterator last_row() { 
      return row_iterator(q.begin() + colCount * rowCount,colCount);
    };
    const_row_iterator last_row() const {
      return const_row_iterator(q.begin() + colCount * rowCount,colCount);
    };
    row_iterator first_row(col_iterator cit) {
      return row_iterator(q.begin() + ((cit - q.begin()) % colCount),colCount);
    };
    const_row_iterator first_row(const_col_iterator cit) const {
      return const_row_iterator(q.begin() + ((cit - q.begin()) % colCount),colCount);
    };
    row_iterator last_row(col_iterator cit) {
      return row_iterator(q.begin() + ((cit - q.begin()) % colCount) + colCount * rowCount,colCount);
    };
    const_row_iterator last_row(const_col_iterator cit) const {
      return const_row_iterator(q.begin() + ((cit - q.begin()) % colCount) + colCount * rowCount,colCount);
    };
    std::pair<row_iterator,row_iterator> rows() { 
      return std::make_pair(first_row(),last_row());
    };
    std::pair<const_row_iterator,const_row_iterator> rows() const { 
      return std::make_pair(first_row(),last_row());
    };
    std::pair<row_iterator,row_iterator> rows(col_iterator cit) { 
      return std::make_pair(first_row(cit),last_row(cit));
    };
    std::pair<const_row_iterator,const_row_iterator> rows(const_col_iterator cit) const { 
      return std::make_pair(first_row(cit),last_row(cit));
    };
    
    col_iterator first_col() {
      return col_iterator(q.begin(),rowCount);
    };
    const_col_iterator first_col() const {
      return const_col_iterator(q.begin(),rowCount);
    };
    col_iterator last_col() {
      return col_iterator(q.begin() + colCount * rowCount,rowCount);
    };
    const_col_iterator last_col() const {
      return const_col_iterator(q.begin() + colCount * rowCount,rowCount);
    };
    col_iterator first_col(row_iterator rit) {
      size_type diff = rit.base() - q.begin();
      return q.begin() + ((diff / colCount) * colCount);
    };
    const_col_iterator first_col(const_row_iterator rit) const {
      size_type diff = rit.base() - q.begin();
      return q.begin() + ((diff / colCount) * colCount);
    };
    col_iterator last_col(row_iterator rit) {
      size_type diff = rit.base() - q.begin();
      return q.begin() + ((diff / colCount + 1) * colCount);
    };
    const_col_iterator last_col(const_row_iterator rit) const {
      size_type diff = rit.base() - q.begin();
      return q.begin() + ((diff / colCount + 1) * colCount);
    };
    std::pair<col_iterator,col_iterator> cols() { 
      return std::make_pair(first_col(),last_col());
    };
    std::pair<const_col_iterator,const_col_iterator> cols() const { 
      return std::make_pair(first_col(),last_col());
    };
    std::pair<col_iterator,col_iterator> cols(row_iterator rit) { 
      return std::make_pair(first_col(rit),last_col(rit));
    };
    std::pair<const_col_iterator,const_col_iterator> cols(const_row_iterator rit) const { 
      return std::make_pair(first_col(rit),last_col(rit));
    };
    

    allocator_type get_allocator() const { return q.get_allocator(); };

    
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** ROW-MAJOR ONLY
     * Standard Assignment operator with standard semantics.
     * Strong exception safety.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator =(const Matrix& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /** ROW-MAJOR ONLY
     * Add-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator +=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  *it += M(i,j);
      return *this;
    };

    /** ROW-MAJOR ONLY
     * Sub-and-store operator with standard semantics.
     * \test PASSED
     */
    template <typename Matrix>
    self& operator -=(const Matrix& M) {
      boost::function_requires< ReadableMatrixConcept<Matrix> >();
      if((M.get_col_count() != colCount) || (M.get_row_count() != rowCount))
	throw std::range_error("Matrix dimension mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i=0;i<rowCount;++i)
        for(size_type j=0;j<colCount;++j,++it)
	  *it -= M(i,j);
      return *this;
    };

    /** WORKS FOR ALL
     * Scalar-multiply-and-store operator with standard semantics.
     * \test PASSED
     */
    self& operator *=(const value_type& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };

    /** WORKS FOR ALL
     * General Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
     self&>::type operator *=(const Matrix& M) {
      self result = *this * M;
      swap(*this,result);
      return *this;
    };
    
    /** WORKS FOR ALL
     * General negation operator for any type of matrices. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General column-major matrix.
     * \test PASSED
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin();it != q.end();++it,++itr)
        *itr = -(*it);
      return result;
    };
    
    
    
    friend mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> transpose(const self& M) {
      return mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator>(M.q,M.colCount,M.rowCount);
    };
    friend mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> transpose_move(self& M) {
      mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.colCount,M.rowCount);
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("rowCount",rowCount)
	& std::pair<std::string, unsigned int>("colCount",colCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("rowCount",rowCount)
	& std::pair<std::string, unsigned int&>("colCount",colCount);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
};






/**
 * Extracts a sub-matrix from this matrix.
 * \param M Matrix from which the sub-matrix is obtained.
 * \param aRowOffset Number of rows before the start of the sub-matrix rows.
 * \param aColOffset Number of columns before the start of the sub-matrix columns.
 * \param aRowCountOut Number of rows of the sub-matrix.
 * \param aColCountOut Number of columns of the sub-matrix.
 * \return The sub-matrix contained in this matrix.
 * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
mat<T,mat_structure::rectangular,Alignment,Allocator> 
  get_block(const mat<T,Structure,Alignment,Allocator>& M, 
	    std::size_t aRowOffset,std::size_t aColOffset,
	    std::size_t aRowCountOut,std::size_t aColCountOut) {
  if((aRowOffset + aRowCountOut > M.get_row_count()) || (aColOffset + aColCountOut > M.get_col_count()))
    throw std::range_error("Matrix dimension mismatch.");
  mat<T,mat_structure::rectangular,Alignment,Allocator> result(aRowCountOut,aColCountOut);
  for(unsigned int i=0;i<aRowCountOut;++i)
    for(unsigned int j=0;j<aColCountOut;++j)
      result(i,j) = M(i+aRowOffset,j+aColOffset);
  return result;
};

/** Sets the sub-part of this matrix to a sub-matrix.
 * \param M The matrix into which the sub-matrix should be set.
 * \param subM A sub-matrix of any type that will be written in the sub-part of this matrix.
 * \param aRowOffset Number of rows before the start of the sub-matrix rows.
 * \param aColOffset Number of columns before the start of the sub-matrix columns.
 * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix>
typename boost::enable_if_c< (((Structure == mat_structure::rectangular) || 
                               (Structure == mat_structure::square)) &&
                               is_readable_matrix<Matrix>::value) >::type
  set_block(mat<T,Structure,Alignment,Allocator>& M,const Matrix& subM,unsigned int aRowOffset,unsigned int aColOffset) {
  if((aRowOffset + subM.get_row_count() > M.get_row_count()) || (aColOffset + subM.get_col_count() > M.get_col_count()))
    throw std::range_error("Matrix dimension mismatch.");
  for(unsigned int i=0;i<subM.get_row_count();++i)
    for(unsigned int j=0;j<subM.get_col_count();++j)
      M(i+aRowOffset,j+aColOffset) = subM(i,j);
};





  
};

#endif







