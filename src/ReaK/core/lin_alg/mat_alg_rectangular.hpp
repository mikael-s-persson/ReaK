/**
 * \file mat_alg_rectangular.hpp
 * 
 * This library implements the specialization of the mat<> template for a 
 * general rectangular matrix (dynamic dimensions) of both column-major and 
 * row-major alignment. This matrix type fulfills all the general matrix 
 * concepts (Readable, Writable, Fully-Writable, Resizable and DynAlloc).
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

#ifndef REAK_MAT_ALG_RECTANGULAR_HPP
#define REAK_MAT_ALG_RECTANGULAR_HPP

#include "mat_alg_general.hpp"

namespace ReaK {
  
  
  
  #if 0
// NOT USED, DEPRECATED
template <>
struct mat_indexer<mat_structure::rectangular, mat_alignment::column_major> {
  std::size_t rowCount, colCount;
  mat_indexer<mat_structure::rectangular, mat_alignment::column_major>(std::size_t aRowCount, std::size_t aColCount) : rowCount(aRowCount), colCount(aColCount) { };
  std::size_t operator()(std::size_t i, std::size_t j) const { return j * rowCount + i; };
};

// NOT USED, DEPRECATED
template <>
struct mat_indexer<mat_structure::rectangular, mat_alignment::row_major> {
  std::size_t rowCount, colCount;
  mat_indexer<mat_structure::rectangular, mat_alignment::row_major>(std::size_t aRowCount, std::size_t aColCount) : rowCount(aRowCount), colCount(aColCount) { };
  std::size_t operator()(std::size_t i, std::size_t j) const { return i * colCount + j; };
};
#endif


template <typename T,
          mat_alignment::tag Alignment,
	  typename Allocator>
struct is_fully_writable_matrix< mat<T,mat_structure::rectangular,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_fully_writable_matrix< mat<T,mat_structure::rectangular,Alignment,Allocator> > type;
};

  
/**
 * This class template specialization implements a matrix with rectangular structure
 * and column-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This matrix type is dynamically resizable.
 * 
 * Models: ReadableMatrixConcept, WritableMatrixConcept, FullyWritableMatrixConcept, 
 * ResizableMatrixConcept, and DynAllocMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
	  typename Allocator>
class mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
  
    typedef stride_iterator< typename container_type::iterator > col_iterator;
    typedef stride_iterator< typename container_type::const_iterator > const_col_iterator;
    typedef typename container_type::iterator row_iterator;
    typedef typename container_type::const_iterator const_row_iterator;
  
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 0);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 0);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
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
             q(0,value_type(),aAlloc),
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
	size_type minN = (colCount < rowCount ? colCount : rowCount) * (rowCount + 1);
	for(size_type i=0;i < minN;i += rowCount + 1)
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
	     
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(self&& M) :
             q(std::move(M.q)),
	     rowCount(std::move(M.rowCount)),
	     colCount(std::move(M.colCount)) { };
#endif

    /**
     * Explicit constructor from a any type of matrix.
     * \test PASSED
     */
    template <typename Matrix>
    explicit mat(const Matrix&  M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                                boost::mpl::not_< boost::is_same<Matrix,self> >::value , void* >::type dummy = NULL) :
             q(M.get_row_count()*M.get_col_count(),T(0.0)),
	     rowCount(M.get_row_count()),
	     colCount(M.get_col_count()) {
      typename container_type::iterator it = q.begin();
      for(size_type j=0;j<colCount;++j)
        for(size_type i=0;i<rowCount;++i,++it)
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
      q[1] = a21;
      q[2] = a12;
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
      q[1] = a21;
      q[2] = a31;
      q[3] = a12;
      q[4] = a22;
      q[5] = a32;
      q[6] = a13;
      q[7] = a23;
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
      q[1] = a21;
      q[2] = a31;
      q[3] = a41;
      q[4] = a12;
      q[5] = a22;
      q[6] = a32;
      q[7] = a42;
      q[8] = a13;
      q[9] = a23;
      q[10] = a33;
      q[11] = a43;
      q[12] = a14;
      q[13] = a24;
      q[14] = a34;
      q[15] = a44;
    };

    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.rowCount,m2.rowCount);
      swap(m1.colCount,m2.colCount);
    };
    
    /**
     * A swap function to swap the matrix with a container of values to fill the matrix.
     * \param m1 The matrix to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param rowCount2 The row-count corresponding to q2's data.
     * \param colCount2 The column-count corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& rowCount2, size_type& colCount2) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.rowCount,rowCount2);
      swap(m1.colCount,colCount2);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
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
    reference operator()(size_type i,size_type j) { return q[j*rowCount + i]; };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    const_reference operator()(size_type i,size_type j) const { return q[j*rowCount + i]; };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return rowCount; };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return colCount; };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };
    /**
     * Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
     * \param sz new dimensions for the matrix.
     * \test PASSED
     */
    void resize(const std::pair<size_type,size_type>& sz) {
      set_row_count(sz.first,true);
      set_col_count(sz.second,true);
    };

    /**
     * Sets the row-count (number of rows) of the matrix.
     * \param aRowCount new number of rows for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_row_count(size_type aRowCount,bool aPreserveData = false) {
      if(aPreserveData) {
	if(aRowCount > rowCount)
	  for(size_type i=0;i<colCount;++i)
	    q.insert(q.begin() + i*aRowCount,aRowCount - rowCount,value_type(0.0));
	else if(aRowCount < rowCount)
	  for(size_type i=1;i<colCount;++i)
	    for(size_type j=0;j<aRowCount;++j)
	      q[i*aRowCount+j] = q[i*rowCount+j];
      } else
	q.resize(colCount * aRowCount,value_type(0.0));
      rowCount = aRowCount;
    };

    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
    void set_col_count(size_type aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(aColCount * rowCount,value_type(0.0));
      colCount = aColCount;
    };
    
    row_iterator first_row() { 
      return q.begin(); 
    };
    const_row_iterator first_row() const { 
      return q.begin();
    };
    row_iterator last_row() { 
      return q.begin() + rowCount; 
    };
    const_row_iterator last_row() const {
      return q.begin() + rowCount; 
    };
    row_iterator first_row(col_iterator cit) {
      size_type diff = cit.base() - q.begin();
      return q.begin() + ((diff / rowCount) * rowCount);
    };
    const_row_iterator first_row(const_col_iterator cit) const {
      size_type diff = cit.base() - q.begin();
      return q.begin() + ((diff / rowCount) * rowCount);
    };
    row_iterator last_row(col_iterator cit) {
      size_type diff = cit.base() - q.begin();
      return q.begin() + ((diff / rowCount + 1) * rowCount);
    };
    const_row_iterator last_row(const_col_iterator cit) const {
      size_type diff = cit.base() - q.begin();
      return q.begin() + ((diff / rowCount + 1) * rowCount);
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
      return col_iterator(q.begin() + ((rit - q.begin()) % rowCount),rowCount);
    };
    const_col_iterator first_col(const_row_iterator rit) const {
      return const_col_iterator(q.begin() + ((rit - q.begin()) % rowCount),rowCount);
    };
    col_iterator last_col(row_iterator rit) {
      return col_iterator(q.begin() + ((rit - q.begin()) % rowCount) + colCount * rowCount,rowCount);
    };
    const_col_iterator last_col(const_row_iterator rit) const {
      return const_col_iterator(q.begin() + ((rit - q.begin()) % rowCount) + colCount * rowCount,rowCount);
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
    
    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };
    
    

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
      self tmp(M);
      swap(*this,tmp);
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
      typename container_type::iterator it = q.begin();
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
      typename container_type::iterator it = q.begin();
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
    
    
    
    /**
     * Transposes the matrix M by a simple copy with a change of alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> transpose(const self& M) {
      return mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator>(M.q,M.colCount,M.rowCount);
    };
    
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> transpose_move(self& M) {
      mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.colCount,M.rowCount);
      return result;
    };

#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> transpose(self&& M) {
      mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.colCount,M.rowCount);
      return result;
    };    
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("rowCount",rowCount)
	& std::pair<std::string, unsigned int>("colCount",colCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("rowCount",tmp);
      rowCount = tmp;
      A & std::pair<std::string, unsigned int&>("colCount",tmp);
      colCount = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
    
  
};


/**
 * This class template specialization implements a matrix with rectangular structure
 * and row-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This matrix type is dynamically resizable.
 * 
 * Models: ReadableMatrixConcept, WritableMatrixConcept, FullyWritableMatrixConcept, 
 * ResizableMatrixConcept, and DynAllocMatrixConcept.
 * 
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
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
	     
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    mat(self&& M) :
             q(std::move(M.q)),
	     rowCount(std::move(M.rowCount)),
	     colCount(std::move(M.colCount)) { };
#endif
	     
    /**
     * Explicit constructor from a any type of matrix.
     * \test PASSED
     */
    template <typename Matrix>
    explicit mat(const Matrix&  M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                                boost::mpl::not_< boost::is_same<Matrix,self> >::value , void* >::type dummy = NULL) :
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
    
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.rowCount,m2.rowCount);
      swap(m1.colCount,m2.colCount);
    };
    
    /**
     * A swap function to swap the matrix with a container of values to fill the matrix.
     * \param m1 The matrix to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param rowCount2 The row-count corresponding to q2's data.
     * \param colCount2 The column-count corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& rowCount2, size_type& colCount2) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.rowCount,rowCount2);
      swap(m1.colCount,colCount2);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
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
    reference operator()(size_type i,size_type j) { return q[i*colCount + j]; };
    /**
     * Matrix indexing accessor for read-only access.
     * \param i Row index.
     * \param j Column index.
     * \return the element at the given position.
     * \test PASSED
     */
    const_reference operator()(size_type i,size_type j) const { return q[i*colCount + j]; };

    /**
     * Gets the row-count (number of rows) of the matrix.
     * \return number of rows of the matrix.
     * \test PASSED
     */
    size_type get_row_count() const throw() { return rowCount; };
    /**
     * Gets the column-count (number of columns) of the matrix.
     * \return number of columns of the matrix.
     * \test PASSED
     */
    size_type get_col_count() const throw() { return colCount; };
    
    /**
     * Gets the row-count and column-count of the matrix, as a std::pair of values.
     * \return the row-count and column-count of the matrix, as a std::pair of values.
     * \test PASSED
     */
    std::pair<size_type,size_type> size() const throw() { return std::make_pair(rowCount,colCount); };
    /**
     * Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
     * \param sz new dimensions for the matrix.
     * \test PASSED
     */
    void resize(const std::pair<size_type,size_type>& sz) {
      set_col_count(sz.second,true);
      set_row_count(sz.first,true);
    };

    /**
     * Sets the column-count (number of columns) of the matrix.
     * \param aColCount new number of columns for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
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

    /**
     * Sets the row-count (number of rows) of the matrix.
     * \param aRowCount new number of rows for the matrix.
     * \param aPreserveData If true, the resizing will preserve all the data it can.
     * \test PASSED
     */
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
    

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
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
    
    
    
    /**
     * Transposes the matrix M by a simple copy with a change of alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> transpose(const self& M) {
      return mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator>(M.q,M.colCount,M.rowCount);
    };
    
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> transpose_move(self& M) {
      mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.colCount,M.rowCount);
      return result;
    };

#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> transpose(self&& M) {
      mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.colCount,M.rowCount);
      return result;
    };
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("rowCount",rowCount)
	& std::pair<std::string, unsigned int>("colCount",colCount);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("rowCount",tmp);
      rowCount = tmp;
      A & std::pair<std::string, unsigned int&>("colCount",tmp);
      colCount = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
};

/*
//NOT USED. DEPRECATED (see friend functions)
template <typename T, typename Allocator>
mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> 
  transpose(mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> M) {
  using std::swap;
  mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> result;
  swap(result.q,M.q);
  swap(result.rowCount, M.colCount);
  swap(result.colCount, M.rowCount);
  return result;
};

//NOT USED. DEPRECATED (see friend functions)
template <typename T, typename Allocator>
mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> 
  transpose(mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> M) {
  using std::swap;
  mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> result;
  swap(result.q,M.q);
  swap(result.rowCount, M.colCount);
  swap(result.colCount, M.rowCount);
  return result;
};


//NOT USED. DEPRECATED (see friend functions)
template <typename T, typename Allocator>
mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> 
  transpose_move(mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator>& M) {
  using std::swap;
  mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator> result;
  swap(result.q,M.q);
  swap(result.rowCount, M.colCount);
  swap(result.colCount, M.rowCount);
  return result;
};

//NOT USED. DEPRECATED (see friend functions)
template <typename T, typename Allocator>
mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> 
  transpose_move(mat<T,mat_structure::rectangular,mat_alignment::column_major,Allocator>& M) {
  using std::swap;
  mat<T,mat_structure::rectangular,mat_alignment::row_major,Allocator> result;
  swap(result.q,M.q);
  swap(result.rowCount, M.colCount);
  swap(result.colCount, M.rowCount);
  return result;
};
*/




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
typename boost::enable_if_c< 
  boost::mpl::and_< boost::mpl::or_< boost::mpl::equal_to< boost::mpl::integral_c< mat_structure::tag, Structure>,
                                                           boost::mpl::integral_c< mat_structure::tag, mat_structure::rectangular> >,
                                     boost::mpl::equal_to< boost::mpl::integral_c< mat_structure::tag, Structure>,
                                                           boost::mpl::integral_c< mat_structure::tag, mat_structure::square> > >,
                    is_readable_matrix<Matrix> >::value, 
void >::type
  set_block(mat<T,Structure,Alignment,Allocator>& M,const Matrix& subM,unsigned int aRowOffset,unsigned int aColOffset) {
  if((aRowOffset + subM.get_row_count() > M.get_row_count()) || (aColOffset + subM.get_col_count() > M.get_col_count()))
    throw std::range_error("Matrix dimension mismatch.");
  for(unsigned int i=0;i<subM.get_row_count();++i)
    for(unsigned int j=0;j<subM.get_col_count();++j)
      M(i+aRowOffset,j+aColOffset) = subM(i,j);
};





  
};

#endif







