
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

#ifndef MAT_SLICES_HPP
#define MAT_SLICES_HPP

#include "mat_alg_general.hpp"

namespace ReaK {


template <typename Matrix>
class mat_row_slice {
  public:
    typedef mat_row_slice<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::row_iterator iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  private:
    Matrix* m;
    size_type offset;
    size_type colIndex;
    size_type count;
  public:
    
    mat_row_slice(Matrix& aM) : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) { };
    
    mat_row_slice(Matrix& aM, 
		  size_type aColIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) { };
		  
    mat_row_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), colIndex(rhs.colIndex), count(rhs.count) { };
    
    
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.colIndex,rhs.colIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator=(const Vector& v) {
      if(v.size() != count) 
	throw std::range_error("Vector dimensions mismatch.");
      for(size_type i = 0; i < count; ++i)
	(*m)(offset + i,colIndex) = v[i];
      return *this;
    };
    
    
    size_type size() const { return count; };
    size_type max_size() const { return count; };
    size_type capacity() const { return count; };
    void resize(size_type sz, const_reference c = value_type()) const { };
    bool empty() const { return false; };
    void reserve(size_type sz) const { };
    
    iterator begin() { return m->first_row(m->first_col() + colIndex) + offset; };
    const_iterator begin() const { return m->first_row(m->first_col() + colIndex) + offset; };
    iterator end() { return m->first_row(m->first_col() + colIndex) + offset + count; };
    const_iterator end() const { return m->first_row(m->first_col() + colIndex) + offset + count; };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator [](size_type i) {
      return (*m)(offset + i,colIndex);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return (*m)(offset + i,colIndex);
    };

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator ()(size_type i) {
      return (*m)(offset + i,colIndex);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(offset + i,colIndex);
    };
    
    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/


    /**
     * Standard add-and-store operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator +=(const Vector& V) {
      if(count != V.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<count;++i)
	(*m)(offset + i,colIndex) += V[i];
      return *this;
    };

    /**
     * Standard sub-and-store operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator -=(const Vector& V) {
      if(count != V.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<count;++i)
	(*m)(offset + i,colIndex) -= V[i];
      return *this;
    };

    /**
     * Scalar multiply-and-store operator for gain.
     * \test PASSED
     */
    self& operator *=(const_reference S) {
      for(size_type i=0;i<count;++i)
	(*m)(offset + i,colIndex) *= S;
      return *this;
    };

    /**
     * Scalar divide-and-store operator for gain.
     * \test PASSED
     */
    self& operator /=(const_reference S) {
      for(size_type i=0;i<count;++i)
	(*m)(offset + i,colIndex) /= S;
      return *this;
    };
  
};




template <typename Matrix>
struct is_readable_vector< mat_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< mat_row_slice<Matrix> > type;
};

template <typename Matrix>
struct is_writable_vector< mat_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_writable_vector< mat_row_slice<Matrix> > type;
};

template <typename Matrix>
struct is_resizable_vector< mat_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< mat_row_slice<Matrix> > type;
};


template <typename Matrix>
struct has_allocator_vector< mat_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_vector< mat_row_slice<Matrix> > type;
};







template <typename Matrix>
class mat_const_row_slice {
  public:
    typedef mat_const_row_slice<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::row_iterator iterator;
    typedef typename mat_traits<Matrix>::const_row_iterator const_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  private:
    const Matrix* m;
    size_type offset;
    size_type colIndex;
    size_type count;
    
    self& operator=(const self&); //non-assignable.

#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_row_slice(Matrix&&);
    mat_const_row_slice(Matrix&&, size_type, size_type, size_type);
#endif
  public:
    
    mat_const_row_slice(const Matrix& aM) : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) { };
    
    mat_const_row_slice(const Matrix& aM, 
		  size_type aColIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) { };
		  
    mat_const_row_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), colIndex(rhs.colIndex), count(rhs.count) { };
    
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.colIndex,rhs.colIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    size_type size() const { return count; };
    size_type max_size() const { return count; };
    size_type capacity() const { return count; };
    void resize(size_type sz, const_reference c = value_type()) const { };
    bool empty() const { return false; };
    void reserve(size_type sz) const { };
    
    const_iterator begin() const { return m->first_row(m->first_col() + colIndex) + offset; };
    const_iterator end() const { return m->first_row(m->first_col() + colIndex) + offset + count; };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return (*m)(offset + i,colIndex);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(offset + i,colIndex);
    };
    
    allocator_type get_allocator() const { return m->get_allocator(); };

  
};




template <typename Matrix>
struct is_readable_vector< mat_const_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< mat_const_row_slice<Matrix> > type;
};

template <typename Matrix>
struct is_writable_vector< mat_const_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector< mat_const_row_slice<Matrix> > type;
};

template <typename Matrix>
struct is_resizable_vector< mat_const_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< mat_const_row_slice<Matrix> > type;
};


template <typename Matrix>
struct has_allocator_vector< mat_const_row_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_vector< mat_const_row_slice<Matrix> > type;
};






















template <typename Matrix>
class mat_col_slice {
  public:
    typedef mat_col_slice<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  private:
    Matrix* m;
    size_type offset;
    size_type rowIndex;
    size_type count;
  public:
    
    mat_col_slice(Matrix& aM) : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) { };
    
    mat_col_slice(Matrix& aM, 
		  size_type aRowIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) { };
		  
    mat_col_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), rowIndex(rhs.rowIndex), count(rhs.count) { };
    
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowIndex,rhs.rowIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator=(const Vector& v) {
      if(v.size() != count) 
	throw std::range_error("Vector dimensions mismatch.");
      for(size_type i = 0; i < count; ++i)
	(*m)(rowIndex,offset + i) = v[i];
      return *this;
    };
    
    
    size_type size() const { return count; };
    size_type max_size() const { return count; };
    size_type capacity() const { return count; };
    void resize(size_type sz, const_reference c = value_type()) const { };
    bool empty() const { return false; };
    void reserve(size_type sz) const { };
    
    iterator begin() { return m->first_col(m->first_row() + rowIndex) + offset; };
    const_iterator begin() const { return m->first_col(m->first_row() + rowIndex) + offset; };
    iterator end() { return m->first_col(m->first_row() + rowIndex) + offset + count; };
    const_iterator end() const { return m->first_col(m->first_row() + rowIndex) + offset + count; };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator [](size_type i) {
      return (*m)(rowIndex,offset + i);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return (*m)(rowIndex,offset + i);
    };

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator ()(size_type i) {
      return (*m)(rowIndex,offset + i);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(rowIndex,offset + i);
    };
    
    allocator_type get_allocator() const { return m->get_allocator(); };
    
    
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/


    /**
     * Standard add-and-store operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator +=(const Vector& V) {
      if(count != V.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<count;++i)
	(*m)(rowIndex,offset + i) += V[i];
      return *this;
    };

    /**
     * Standard sub-and-store operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator -=(const Vector& V) {
      if(count != V.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<count;++i)
	(*m)(rowIndex,offset + i) -= V[i];
      return *this;
    };

    /**
     * Scalar multiply-and-store operator for gain.
     * \test PASSED
     */
    self& operator *=(const_reference S) {
      for(size_type i=0;i<count;++i)
	(*m)(rowIndex,offset + i) *= S;
      return *this;
    };

    /**
     * Scalar divide-and-store operator for gain.
     * \test PASSED
     */
    self& operator /=(const_reference S) {
      for(size_type i=0;i<count;++i)
	(*m)(rowIndex,offset + i) /= S;
      return *this;
    };
  
};




template <typename Matrix>
struct is_readable_vector< mat_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< mat_col_slice<Matrix> > type;
};

template <typename Matrix>
struct is_writable_vector< mat_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix<Matrix>::value );
  typedef is_writable_vector< mat_col_slice<Matrix> > type;
};

template <typename Matrix>
struct is_resizable_vector< mat_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< mat_col_slice<Matrix> > type;
};


template <typename Matrix>
struct has_allocator_vector< mat_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_vector< mat_col_slice<Matrix> > type;
};











template <typename Matrix>
class mat_const_col_slice {
  public:
    typedef mat_const_col_slice<Matrix> self;
    typedef typename mat_traits<Matrix>::allocator_type allocator_type;
    
    typedef typename mat_traits<Matrix>::value_type value_type;
    
    typedef typename mat_traits<Matrix>::reference reference;
    typedef typename mat_traits<Matrix>::const_reference const_reference;
    typedef typename mat_traits<Matrix>::pointer pointer;
    typedef typename mat_traits<Matrix>::const_pointer const_pointer;
  
    typedef typename mat_traits<Matrix>::col_iterator iterator;
    typedef typename mat_traits<Matrix>::const_col_iterator const_iterator;
  
    typedef typename mat_traits<Matrix>::size_type size_type;
    typedef typename mat_traits<Matrix>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  private:
    const Matrix* m;
    size_type offset;
    size_type rowIndex;
    size_type count;
    
    self& operator=(const self&); //non-assignable.
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_col_slice(Matrix&&);
    mat_const_col_slice(Matrix&&, size_type, size_type, size_type);
#endif
  public:
    
    mat_const_col_slice(const Matrix& aM) : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) { };
    
    mat_const_col_slice(const Matrix& aM, 
		  size_type aRowIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) { };
		  
    mat_const_col_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), rowIndex(rhs.rowIndex), count(rhs.count) { };
    
    
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowIndex,rhs.rowIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    size_type size() const { return count; };
    size_type max_size() const { return count; };
    size_type capacity() const { return count; };
    void resize(size_type sz, const_reference c = value_type()) const { };
    bool empty() const { return false; };
    void reserve(size_type sz) const { };
    
    const_iterator begin() const { return m->first_col(m->first_row() + rowIndex) + offset; };
    const_iterator end() const { return m->first_col(m->first_row() + rowIndex) + offset + count; };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return (*m)(rowIndex,offset + i);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(rowIndex,offset + i);
    };
    
    allocator_type get_allocator() const { return m->get_allocator(); };
    
};




template <typename Matrix>
struct is_readable_vector< mat_const_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< mat_const_col_slice<Matrix> > type;
};

template <typename Matrix>
struct is_writable_vector< mat_const_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector< mat_const_col_slice<Matrix> > type;
};

template <typename Matrix>
struct is_resizable_vector< mat_const_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< mat_const_col_slice<Matrix> > type;
};


template <typename Matrix>
struct has_allocator_vector< mat_const_col_slice<Matrix> > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix<Matrix>::value );
  typedef has_allocator_vector< mat_const_col_slice<Matrix> > type;
};










template <typename Matrix>
struct mat_slice_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  Matrix& m;
  mat_slice_factory(Matrix& aM) : m(aM) { };
  mat_col_slice<Matrix> operator()(size_type row,
				   const std::pair<size_type,size_type>& cols) {
    return mat_col_slice<Matrix>(m,row,cols.first,cols.second - cols.first + 1);
  };
  mat_row_slice<Matrix> operator()(const std::pair<size_type,size_type>& rows,
				   size_type col) {
    return mat_row_slice<Matrix>(m,col,rows.first,rows.second - rows.first + 1);
  };
};

template <typename Matrix>
struct mat_const_slice_factory {
  typedef typename mat_traits<Matrix>::size_type size_type;
  
  const Matrix& m;
  mat_const_slice_factory(const Matrix& aM) : m(aM) { };
  mat_const_col_slice<Matrix> operator()(size_type row,
				   const std::pair<size_type,size_type>& cols) {
    return mat_const_col_slice<Matrix>(m,row,cols.first,cols.second - cols.first + 1);
  };
  mat_const_row_slice<Matrix> operator()(const std::pair<size_type,size_type>& rows,
				   size_type col) {
    return mat_const_row_slice<Matrix>(m,col,rows.first,rows.second - rows.first + 1);
  };
};


template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_slice_factory<Matrix> >::type slice(Matrix& M) { return mat_slice_factory<Matrix>(M); };

template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_const_slice_factory<Matrix> >::type slice(const Matrix& M) { return mat_const_slice_factory<Matrix>(M); };














};

#endif


















