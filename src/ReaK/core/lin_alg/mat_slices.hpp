/**
 * \file mat_slices.hpp
 * 
 * This library provides a number of classes and functions to take a slice of a matrix.
 * A slice is simply a reduction of the order (i.e. a matrix is a 2nd order tensor, 
 * and a vector is a 1st order tensor). So, this library allows one to extract a row 
 * or column from a matrix and have this row or column appear as a vector. In other words,
 * this library provides adaptors from matrices to vectors.
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

#ifndef REAK_MAT_SLICES_HPP
#define REAK_MAT_SLICES_HPP

#include "mat_alg_general.hpp"

#include "vect_alg.hpp"

namespace ReaK {

/**
 * This class template can be used to view a column of a matrix as a vector (i.e. a row-slice). It 
 * takes the matrix object by reference (internally by pointer). 
 * 
 * Models: ReadableVectorConcept and WritableVectorConcept (if the matrix is writable)
 * \tparam Matrix A matrix type.
 */
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
    Matrix* m; ///< Holds the reference to the matrix object.
    size_type offset; ///< Holds the offset from the start of the column.
    size_type colIndex; ///< Holds the index of column of the slice.
    size_type count; ///< Holds the number of elements of the column to take.
  public:
    
    /**
     * Constructs the row-slice from a matrix M, taking the entire first column.
     * \param aM The matrix from which to take the slice.
     */
    mat_row_slice(Matrix& aM) : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) { };
    
    /**
     * Constructs the row-slice from a matrix M, taking the aColIndex column 
     * from aOffset with aCount elements.
     * \param aM The matrix from which to take the slice.
     * \param aColIndex The column to use for the slice.
     * \param aOffset The offset into the column used for the slice.
     * \param aCount The number of elements to take in the slice.
     */
    mat_row_slice(Matrix& aM, 
		  size_type aColIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) { };
		  
    /**
     * Standard copy-constructor (shallow).
     */
    mat_row_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), colIndex(rhs.colIndex), count(rhs.count) { };
    
    
    /**
     * Standard swap function (shallow)
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.colIndex,rhs.colIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    /**
     * Standard assignment operator for any readable vector type.
     * \tparam Vector A readable vector type.
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator=(const Vector& v) {
      if(v.size() != count) 
	throw std::range_error("Vector dimensions mismatch.");
      for(size_type i = 0; i < count; ++i)
	(*m)(offset + i,colIndex) = v[i];
      return *this;
    };
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return count; };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return count; };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return count; };
    /**
     * Resizes the vector.
     * \param sz The new size for the vector.
     * \param c The value to fill any additional elements to the vector.
     */
    void resize(size_type sz, const_reference c = value_type()) const { RK_UNUSED(sz); RK_UNUSED(c); };
    /**
     * Checks whether the vector is empty.
     */
    bool empty() const { return false; };
    /**
     * Reserve a certain amount of capacity for future additions.
     * \param sz The new capacity for the vector.
     */
    void reserve(size_type sz) const { RK_UNUSED(sz); };
    
    /**
     * Returns an iterator to the first element of the vector.
     */
    iterator begin() { return m->first_row(m->first_col() + colIndex) + offset; };
    /**
     * Returns an const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return m->first_row(m->first_col() + colIndex) + offset; };
    /**
     * Returns an iterator to the one-passed-last element of the vector.
     */
    iterator end() { return m->first_row(m->first_col() + colIndex) + offset + count; };
    /**
     * Returns an const-iterator to the one-passed-last element of the vector.
     */
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
     * Sub-vector operator, accessor for read/write.
     * \test PASSED
     */
    vect_ref_view<self> operator[](const std::pair<size_type,size_type>& r) {
      return sub(*this)[r];
    };

    /**
     * Sub-vector operator, accessor for read only.
     * \test PASSED
     */
    vect_const_ref_view<self> operator[](const std::pair<size_type,size_type>& r) const {
      return sub(*this)[r];
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
    
    /**
     * Returns the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return m->get_allocator(); };
    
  
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
struct vect_copy< mat_row_slice<Matrix> > {
  typedef vect_n< typename vect_traits< mat_row_slice<Matrix> >::value_type > type;
};





/**
 * This class template can be used to view a column of a matrix as a const-vector (i.e. a row-slice). It 
 * takes the matrix object by const-reference (internally by pointer). 
 * 
 * Models: ReadableVectorConcept
 * \tparam Matrix A matrix type.
 */
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
    const Matrix* m; ///< Holds the reference to the matrix object.
    size_type offset; ///< Holds the offset from the start of the column.
    size_type colIndex; ///< Holds the index of column of the slice.
    size_type count; ///< Holds the number of elements of the column to take.
    
    self& operator=(const self&); //non-assignable.

#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_row_slice(Matrix&&);
    mat_const_row_slice(Matrix&&, size_type, size_type, size_type);
#endif
  public:
    
    /**
     * Constructs the row-slice from a matrix M, taking the entire first column.
     * \param aM The matrix from which to take the slice.
     */
    mat_const_row_slice(const Matrix& aM) : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) { };
    
    /**
     * Constructs the row-slice from a matrix M, taking the aColIndex column 
     * from aOffset with aCount elements.
     * \param aM The matrix from which to take the slice.
     * \param aColIndex The column to use for the slice.
     * \param aOffset The offset into the column used for the slice.
     * \param aCount The number of elements to take in the slice.
     */
    mat_const_row_slice(const Matrix& aM, 
		  size_type aColIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) { };
		  
    /**
     * Standard copy-constructor (shallow).
     */
    mat_const_row_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), colIndex(rhs.colIndex), count(rhs.count) { };
    
    
    /**
     * Standard swap function (shallow)
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.colIndex,rhs.colIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return count; };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return count; };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return count; };
    /**
     * Resizes the vector.
     * \param sz The new size for the vector.
     * \param c The value to fill any additional elements to the vector.
     */
    void resize(size_type sz, const_reference c = value_type()) const { };
    /**
     * Checks whether the vector is empty.
     */
    bool empty() const { return false; };
    /**
     * Reserve a certain amount of capacity for future additions.
     * \param sz The new capacity for the vector.
     */
    void reserve(size_type sz) const { };
    
    /**
     * Returns an const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return m->first_row(m->first_col() + colIndex) + offset; };
    /**
     * Returns an const-iterator to the one-passed-last element of the vector.
     */
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
     * Sub-vector operator, accessor for read only.
     * \test PASSED
     */
    vect_const_ref_view<self> operator[](const std::pair<size_type,size_type>& r) const {
      return sub(*this)[r];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(offset + i,colIndex);
    };
    
    /**
     * Returns the allocator object of the underlying container.
     */
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
struct vect_copy< mat_const_row_slice<Matrix> > {
  typedef vect_n< typename vect_traits< mat_const_row_slice<Matrix> >::value_type > type;
};




















/**
 * This class template can be used to view a row of a matrix as a vector (i.e. a column-slice). It 
 * takes the matrix object by reference (internally by pointer). 
 * 
 * Models: ReadableVectorConcept and WritableVectorConcept (if the matrix is writable)
 * \tparam Matrix A matrix type.
 */
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
    Matrix* m; ///< Holds the reference to the matrix object.
    size_type offset; ///< Holds the offset from the start of the row.
    size_type rowIndex; ///< Holds the index of row of the slice.
    size_type count; ///< Holds the number of elements of the row to take.
  public:
    
    /**
     * Constructs the column-slice from a matrix M, taking the entire first row.
     * \param aM The matrix from which to take the slice.
     */
    mat_col_slice(Matrix& aM) : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) { };
    
    /**
     * Constructs the column-slice from a matrix M, taking the aRowIndex row 
     * from aOffset with aCount elements.
     * \param aM The matrix from which to take the slice.
     * \param aRowIndex The row to use for the slice.
     * \param aOffset The offset into the row used for the slice.
     * \param aCount The number of elements to take in the slice.
     */
    mat_col_slice(Matrix& aM, 
		  size_type aRowIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) { };
		  
    /**
     * Standard copy-constructor (shallow).
     */
    mat_col_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), rowIndex(rhs.rowIndex), count(rhs.count) { };
    
    
    /**
     * Standard swap function (shallow)
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowIndex,rhs.rowIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    /**
     * Standard assignment operator for any readable vector type.
     * \tparam Vector A readable vector type.
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value,
    self& >::type operator=(const Vector& v) {
      if(v.size() != count) 
	throw std::range_error("Vector dimensions mismatch.");
      for(size_type i = 0; i < count; ++i)
	(*m)(rowIndex,offset + i) = v[i];
      return *this;
    };
    
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return count; };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return count; };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return count; };
    /**
     * Resizes the vector.
     * \param sz The new size for the vector.
     * \param c The value to fill any additional elements to the vector.
     */
    void resize(size_type sz, const_reference c = value_type()) const { };
    /**
     * Checks whether the vector is empty.
     */
    bool empty() const { return false; };
    /**
     * Reserve a certain amount of capacity for future additions.
     * \param sz The new capacity for the vector.
     */
    void reserve(size_type sz) const { };
    
    /**
     * Returns an iterator to the first element of the vector.
     */
    iterator begin() { return m->first_col(m->first_row() + rowIndex) + offset; };
    /**
     * Returns an const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return m->first_col(m->first_row() + rowIndex) + offset; };
    /**
     * Returns an iterator to the one-passed-last element of the vector.
     */
    iterator end() { return m->first_col(m->first_row() + rowIndex) + offset + count; };
    /**
     * Returns an const-iterator to the one-passed-last element of the vector.
     */
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
     * Sub-vector operator, accessor for read/write.
     * \test PASSED
     */
    vect_ref_view<self> operator[](const std::pair<size_type,size_type>& r) {
      return sub(*this)[r];
    };

    /**
     * Sub-vector operator, accessor for read only.
     * \test PASSED
     */
    vect_const_ref_view<self> operator[](const std::pair<size_type,size_type>& r) const {
      return sub(*this)[r];
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
    
    /**
     * Returns the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return m->get_allocator(); };
    

  
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
struct vect_copy< mat_col_slice<Matrix> > {
  typedef vect_n< typename vect_traits< mat_col_slice<Matrix> >::value_type > type;
};











/**
 * This class template can be used to view a row of a matrix as a vector (i.e. a column-slice). It 
 * takes the matrix object by const-reference (internally by const-pointer). 
 * 
 * Models: ReadableVectorConcept
 * \tparam Matrix A matrix type.
 */
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
    const Matrix* m; ///< Holds the reference to the matrix object.
    size_type offset; ///< Holds the offset from the start of the row.
    size_type rowIndex; ///< Holds the index of row of the slice.
    size_type count; ///< Holds the number of elements of the row to take.
    
    self& operator=(const self&); //non-assignable.
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    mat_const_col_slice(Matrix&&);
    mat_const_col_slice(Matrix&&, size_type, size_type, size_type);
#endif
  public:
    
    /**
     * Constructs the column-slice from a matrix M, taking the entire first row.
     * \param aM The matrix from which to take the slice.
     */
    mat_const_col_slice(const Matrix& aM) : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) { };
    
    /**
     * Constructs the column-slice from a matrix M, taking the aRowIndex row 
     * from aOffset with aCount elements.
     * \param aM The matrix from which to take the slice.
     * \param aRowIndex The row to use for the slice.
     * \param aOffset The offset into the row used for the slice.
     * \param aCount The number of elements to take in the slice.
     */
    mat_const_col_slice(const Matrix& aM, 
		  size_type aRowIndex, 
		  size_type aOffset, 
		  size_type aCount) : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) { };
		  
    /**
     * Standard copy-constructor (shallow).
     */
    mat_const_col_slice(const self& rhs) : m(rhs.m), offset(rhs.offset), rowIndex(rhs.rowIndex), count(rhs.count) { };
    
    
    /**
     * Standard swap function (shallow)
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.m,rhs.m);
      swap(lhs.offset,rhs.offset);
      swap(lhs.rowIndex,rhs.rowIndex);
      swap(lhs.count,rhs.count);
      return;
    };
    
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return count; };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return count; };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return count; };
    /**
     * Resizes the vector.
     * \param sz The new size for the vector.
     * \param c The value to fill any additional elements to the vector.
     */
    void resize(size_type sz, const_reference c = value_type()) const { };
    /**
     * Checks whether the vector is empty.
     */
    bool empty() const { return false; };
    /**
     * Reserve a certain amount of capacity for future additions.
     * \param sz The new capacity for the vector.
     */
    void reserve(size_type sz) const { };
    
    /**
     * Returns an const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return m->first_col(m->first_row() + rowIndex) + offset; };
    /**
     * Returns an const-iterator to the one-passed-last element of the vector.
     */
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
     * Sub-vector operator, accessor for read only.
     * \test PASSED
     */
    vect_const_ref_view<self> operator[](const std::pair<size_type,size_type>& r) const {
      return sub(*this)[r];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return (*m)(rowIndex,offset + i);
    };
    
    /**
     * Returns the allocator object of the underlying container.
     */
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
struct vect_copy< mat_const_col_slice<Matrix> > {
  typedef vect_n< typename vect_traits< mat_const_col_slice<Matrix> >::value_type > type;
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

/**
 * This function template prepares a matrix to be sliced via an intermediate factor class.
 * Once the factor class has been created with this function, one can use the range() function
 * to define to row or column range to use in the slice. Given matrix M, one can create 
 * a slice as follow: slice(M)(range(0,M.get_row_count()-1), 0) which creates a vector view 
 * on the entire first column.
 * \tparam Matrix A readable matrix type.
 * \param M the matrix to slice.
 * \return A factory class to create a slice from a row or column range and a row or column index.
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_slice_factory<Matrix> >::type slice(Matrix& M) { return mat_slice_factory<Matrix>(M); };

/**
 * This function template prepares a matrix to be sliced via an intermediate factor class.
 * Once the factor class has been created with this function, one can use the range() function
 * to define to row or column range to use in the slice. Given matrix M, one can create 
 * a slice as follow: slice(M)(range(0,M.get_row_count()-1), 0) which creates a vector view 
 * on the entire first column.
 * \tparam Matrix A readable matrix type.
 * \param M the matrix to slice.
 * \return A factory class to create a slice from a row or column range and a row or column index.
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
mat_const_slice_factory<Matrix> >::type slice(const Matrix& M) { return mat_const_slice_factory<Matrix>(M); };














};

#endif


















