/**
 * \file rotations_2D.hpp
 *
 * This library declares all geometric 2D rotation classes for fixed (2,3) and variable dimensions.
 *
 * Note: All matrix memory is organized, by default, such that columns are concatenated. This
 *       was found to be a more efficient representation since columns often have
 *       more significances than rows (representing basis vectors for example).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef ROTATIONS_2D_HPP
#define ROTATIONS_2D_HPP

#include "mat_alg.hpp"


namespace ReaK {


//forward declaration, for friend declaration.
template <class T> class trans_mat_2D;

/**
 * This class is a rotation matrix (proper orthogonal) of dimension 2 by 2.
 * \test All unit test for this class have been passed!
 */
template <typename T>
class rot_mat_2D : public serialization::serializable {
  public:
    typedef rot_mat_2D<T> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef void col_iterator;
    typedef void const_col_iterator;
    typedef void row_iterator;
    typedef void const_row_iterator;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 2);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 2);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::orthogonal);
    
    
  private:
    T q[2];

    explicit rot_mat_2D(const_reference cos_a,const_reference sin_a) { q[0] = cos_a; q[1] = sin_a; };

  public:
    friend class trans_mat_2D<T>;

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default Constructor with no rotation.
     * \test PASSED
     */
    rot_mat_2D() { q[0] = 1.0; q[1] = 0.0; };

    /**
     * Explicit constructor of a rotation matrix from a rotation angle.
     * \test PASSED
     */
    explicit rot_mat_2D(const_reference Angle) {
      q[0] = cos(Angle);
      q[1] = sin(Angle);
    };
    
    
    /**
     * Explicit constructor of a rotation matrix from a cosine and sine of an angle.
     * \test PASSED
     */
    explicit rot_mat_2D(vect<value_type,2> v) {
      v = unit(v);
      q[0] = v[0];
      q[1] = v[1];
    };
    
    template <typename Matrix>
    explicit rot_mat_2D(const Matrix& R, typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                                                      !boost::is_same<self, Matrix>::value, void* >::type dummy = NULL) {
      if((R.get_col_count() != 2) || (R.get_row_count() != 2))
	throw std::range_error("Right-hand-side of 2D rotation matrix assignment is not a 2x2 matrix!");
      vect<value_type,2> v = unit(vect<value_type,2>(R(0,0),R(1,0)));
      q[0] = v[0]; q[1] = v[1];
    };

    //Default copy constructor is fine for this class.

    /**
     * Destructor.
     * \test PASSED
     */
    ~rot_mat_2D() { };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides a copy of the rotation matrix as an ordinary 2x2 matrix.
     * \return this rotation matrix as a normal column-major matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      return mat<T,mat_structure::square>(q[0],-q[1],q[1],q[0]);
    };

    /**
     * Returns the angle (-pi .. pi) of the rotation matrix.
     * \test PASSED
     */
    value_type getAngle() const {
      using std::atan2;
      return atan2(q[1],q[0]);
    };

    /**
     * Sets the angle (in radians) of the rotation matrix.
     * \test PASSED
     */
    void setAngle(const_reference Angle) {
      using std::sin;
      using std::cos;
      q[0] = cos(Angle);
      q[1] = sin(Angle);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    value_type operator [](size_type i) const {
      if(i >= 4)
	throw std::range_error("Matrix index out of range.");
      if(i == 3)
	return q[0];
      if(i == 2)
	return -q[1];
      return q[i];
    };

    /**
     * Array double-indexing operator, ith row and jth column, accessor for read only.
     * \test PASSED
     */
    value_type operator ()(size_type i,size_type j) const {
      if((i >= 2) || (j >= 2))
	throw std::range_error("Matrix index out of range.");
      if((j == 1) && (i == 0))
	return -q[1];
      if((j == 1) && (i == 1))
	return q[0];
      return q[i];
    };
    
    size_type get_row_count() const { return 2; };
    size_type get_col_count() const { return 2; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Assignment operator.
     * \test PASSED
     */
    self& operator =(const self& R) {
      q[0] = R.q[0];
      q[1] = R.q[1];
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<self, Matrix>::value,
    self& >::type operator =(const Matrix& R) {
      if((R.get_col_count() != 2) || (R.get_row_count() != 2))
	throw std::range_error("Right-hand-side of 2D rotation matrix assignment is not a 2x2 matrix!");
      vect<value_type,2> v = unit(vect<value_type,2>(R(0,0),R(1,0)));
      q[0] = v[0]; q[1] = v[1];
      return *this;
    };

    /**
     * Multiply-and-store operator.
     * \test PASSED
     */
    self& operator *=(const self& R) {
      value_type tmp = q[0] * R.q[0] - q[1] * R.q[1];
      q[1] = q[1] * R.q[0] + q[0] * R.q[1];
      q[0] = tmp;
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<self, Matrix>::value,
    self& >::type operator *=(const Matrix& R) {
      if((R.get_col_count() != 2) || (R.get_row_count() != 2))
	throw std::range_error("Right-hand-side of 2D rotation matrix assignment is not a 2x2 matrix!");
      vect<value_type,2> v = unit(vect<value_type,2>(R(0,0),R(1,0)));
      value_type tmp = q[0] * v[0] - q[1] * v[1];
      q[1] = q[1] * v[0] + q[0] * v[1];
      q[0] = tmp;
      return *this;
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Rotation matrix multiplication.
     * \test PASSED
     */
    friend
    self operator *(const self& R1,const self& R2) {
      return self( R1.q[0] * R2.q[0] - R1.q[1] * R2.q[1],
                   R1.q[1] * R2.q[0] + R1.q[0] * R2.q[1] );
    };

    /**
     * Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
    Matrix >::type operator *(const self& R, const Matrix& M) {
      if(M.get_row_count() != 2)
	throw std::range_error("Matrix M's row count is not 2, cannot perform 2D rotation!");
      Matrix result(M);
      for(unsigned int jj=0;jj < result.get_col_count();++jj) {
	result(0,jj) = R.q[0] * M(0,jj) - R.q[1] * M(1,jj);
	result(1,jj) = R.q[1] * M(0,jj) + R.q[0] * M(1,jj);
      };
      return result;
    };

    /**
     * 2D Rotation matrix times a column 2D vector.
     * \test PASSED
     */
    friend
    vect<value_type,2> operator *(const self& R, const vect<value_type,2>& V) {
      return vect<value_type,2>(V[0] * R.q[0] - V[1] * R.q[1],V[0] * R.q[1] + V[1] * R.q[0]);
    };

    /**
     * Row 2D vector times a rotation matrix.
     * \test PASSED
     */
    friend vect<value_type,2> operator *(const vect<value_type,2>& V,const self& R) {
      return vect<value_type,2>(V[0] * R.q[0] + V[1] * R.q[1], V[1] * R.q[0] - V[0] * R.q[1]);
    };
    
    /**
     * Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix> 
    friend 
    typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
    Matrix >::type operator *(const Matrix& M,const self& R) {
      if(M.get_col_count() != 2)
        throw std::range_error("Matrix M's column count is not 2, cannot perform 2D rotation!");
      Matrix result(M);
      for(unsigned int i=0;i<result.get_row_count();++i) {
        result(i,0) = R.q[0] * M(i,0) + R.q[1] * M(i,1);
        result(i,1) = R.q[0] * M(i,1) - R.q[1] * M(i,0);
      };
      return result;
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality Comparison operator.
     * \test PASSED
     */
    friend
    bool operator ==(const self& R1, const self& R2) {
      return ((R1.q[0] == R2.q[0]) && (R1.q[1] == R2.q[1]));
    };

    /**
     * Inequality Comparison operator.
     * \test PASSED
     */
    friend
    bool operator !=(const self& R1, const self& R2) {
      return ((R1.q[0] != R2.q[0]) || (R1.q[1] != R2.q[1]));
    };

/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/

    /**
     * Creates the transpose matrix.
     * \test PASSED
     */
    friend 
    self transpose(const self& R) {
      return self(R.q[0],-R.q[1]);
    };
    
    /**
     * Creates the transpose matrix.
     * \test PASSED
     */
    friend 
    self transpose_move(const self& R) {
      return self(R.q[0],-R.q[1]);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& R) {
      return value_type(2.0)*R.q[0];
    };

    /**
     * Gets the determinant of the matrix.
     * \test PASSED
     */
    friend
    value_type determinant(const self&) {
      return value_type(1.0);
    };

    /**
     * Gets the inverse of the matrix.
     * \test PASSED
     */
    friend
    self invert(const self& R) {
      return transpose(R);
    };

    /**
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      return mat<value_type,mat_structure::symmetric>(q[0],value_type(0.0),q[0]);
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      return mat<value_type,mat_structure::skew_symmetric>(-q[1]);
    };


/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_ALIAS("cos",q[0])
        & RK_SERIAL_SAVE_WITH_ALIAS("sin",q[1]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("cos",q[0])
        & RK_SERIAL_LOAD_WITH_ALIAS("sin",q[1]);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x00000016,1,"rot_mat_2D",serialization::serializable)

};

/**
 * Prints a 2D rotation matrix to a standard output stream (<<) as "(angle = a)".
 * \test PASSED
 */
template <typename T>
std::ostream& operator <<(std::ostream& out_stream,const rot_mat_2D<T>& R) {
  out_stream << "(angle = " << R.getAngle() << ")";
  return out_stream;
};


template <typename T>
struct is_readable_matrix< rot_mat_2D<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< rot_mat_2D<T> > type;
};







/**
 * This class is a transformation matrix 3 by 3, i.e. to rotate and translate a 2D vector.
 * \test All unit tests for this class have been passed!
 */
template <typename T>
class trans_mat_2D : public serialization::serializable {
  public:
    typedef trans_mat_2D<T> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef void col_iterator;
    typedef void const_col_iterator;
    typedef void row_iterator;
    typedef void const_row_iterator;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 3);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 3);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::square);
    
    
  private:
    value_type q[9];

    explicit trans_mat_2D(const_reference a11, const_reference a12, const_reference a13, 
			  const_reference a21, const_reference a22, const_reference a23) {
      q[0] = a11; q[3] = a12; q[6] = a13;
      q[1] = a21; q[4] = a22; q[7] = a23;
      q[2] = 0.0; q[5] = 0.0; q[8] = 1.0;
    };
  public:

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor.
     * \test PASSED
     */
    trans_mat_2D() {
      q[0] = 1.0; q[3] = 0.0; q[6] = 0.0;
      q[1] = 0.0; q[4] = 1.0; q[7] = 0.0;
      q[2] = 0.0; q[5] = 0.0; q[8] = 1.0;
    };

    /**
     * Constructor from a rotation angle and a translation vector.
     * \test PASSED
     */
    explicit trans_mat_2D(const_reference Angle, vect<value_type,2> Translation = vect<value_type,2>()) {
      q[4] =  (q[0] = cos(Angle));
      q[3] = -(q[1] = sin(Angle));
      q[6] = Translation[0];
      q[7] = Translation[1];
      q[2] = 0.0;
      q[5] = 0.0;
      q[8] = 1.0;
    };
    
    explicit trans_mat_2D(const rot_mat_2D<value_type>& R, vect<value_type,2> Translation = vect<value_type,2>()) {
      q[4] =  (q[0] = R(0,0));
      q[3] = -(q[1] = R(1,0));
      q[6] = Translation[0];
      q[7] = Translation[1];
      q[2] = 0.0;
      q[5] = 0.0;
      q[8] = 1.0;
    };
    
    template <typename Matrix>
    explicit trans_mat_2D(const Matrix& M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value && 
                                                                        !boost::is_same<self, Matrix>::value, void* >::type dummy = NULL ) {
      if((M.get_col_count() != 3) || (M.get_row_count() != 3))
	throw std::range_error("Right-hand-side of 2D transformation matrix assignment is not a 3x3 matrix!");
      vect<value_type,2> v = unit(vect<value_type,2>(M(0,0),M(1,0)));
      q[0] = v[0]; q[1] = v[1]; q[2] = 0.0;
      q[3] = -q[1]; q[4] = q[0]; q[5] = 0.0;
      q[6] = M(0,2);
      q[7] = M(1,2);
      q[8] = 1.0;
    };

    // Copy Constructor. Default is fine. \test PASSED

    /**
     * Destructor.
     * \test PASSED
     */
    ~trans_mat_2D() { };


/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides a copy of the transformation matrix as an ordinary 3x3 matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      return mat<value_type,mat_structure::square>(q[0],q[3],q[6],q[1],q[4],q[7],0.0,0.0,1.0);
    };

    /**
     * Provides a copy of the rotation matrix part of the transformation matrix.
     * \test PASSED
     */
    rot_mat_2D<value_type> getRotMat() const {
      return rot_mat_2D<value_type>(q[0],q[1]);
    };

    /**
     * Sets the rotation part of the transformation matrix.
     * \test PASSED
     */
    void setRotMat(const rot_mat_2D<value_type>& R) {
      q[4] = (q[0] = R.q[0]);
      q[3] = -(q[1] = R.q[1]);
    };

    /**
     * Returns the angle of the rotation matrix.
     * \test PASSED
     */
    value_type getAngle() const {
      return atan2(q[1],q[0]);
    };

    /**
     * Returns the angle of the rotation matrix.
     * \test PASSED
     */
    void setAngle(const_reference Angle) {
      q[4] = (q[0] = cos(Angle));
      q[3] = -(q[1] = sin(Angle));
    };

    /**
     * Provides a copy of the translation part of the transformation matrix.
     * \test PASSED
     */
    vect<value_type,2> getTranslation() const {
      return vect<value_type,2>(q[6],q[7]);
    };

    /**
     * Sets the translation part of the transformation matrix.
     * \test PASSED
     */
    void setTranslation(const vect<value_type,2>& Translation) {
      q[6] = Translation.q[0];
      q[7] = Translation.q[1];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      if(i >= 9)
	throw std::range_error("Matrix index out of range.");
      return q[i];
    };

    /**
     * Array double-indexing operator, ith row and jth column, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i,size_type j) const {
      if((i >= 3) || (j >= 3))
	throw std::range_error("Matrix index out of range.");
      return q[j*3+i];
    };
    
    size_type get_row_count() const { return 3; };
    size_type get_col_count() const { return 3; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard Assignment operator.
     * \test PASSED
     */
    self& operator =(const self& M) {
      q[0] = M.q[0];
      q[1] = M.q[1];
      q[2] = 0.0;
      q[3] = M.q[3];
      q[4] = M.q[4];
      q[5] = 0.0;
      q[6] = M.q[6];
      q[7] = M.q[7];
      q[8] = 1.0;
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value, 
    self& >::type operator =(const Matrix& M) {
      if((M.get_col_count() != 3) || (M.get_row_count() != 3))
	throw std::range_error("Right-hand-side of 2D transformation matrix assignment is not a 3x3 matrix!");
      vect<value_type,2> v = unit(vect<value_type,2>(M(0,0),M(1,0)));
      q[0] = v[0]; q[1] = v[1]; q[2] = 0.0;
      q[3] = -q[1]; q[4] = q[0]; q[5] = 0.0;
      q[6] = M(0,2);
      q[7] = M(1,2);
      q[8] = 1.0;
    };

    /**
     * Multiply-and-store operator.
     * \test PASSED
     */
    self& operator *=(const self& M) {
      return (*this = (*this) * M);
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value, 
    self& >::type operator *=(const Matrix& M) {
      return (*this = (*this) * M);
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/


    /**
     * Multiplication with rotation matrix.
     */
    friend 
    self operator *(const rot_mat_2D<value_type>& R, const self& M) {
      return self(R) * M;
    };
    
    /**
     * Multiplication with rotation matrix.
     */
    friend 
    self operator *(const self& M, const rot_mat_2D<value_type>& R) {
      return M * self(R);
    };
    
    /**
     * Matrix multiplication.
     * \test PASSED
     */
    friend
    self operator *(const self& M1,const self& M2) {
      return trans_mat_2D<T>(M1.q[0]*M2.q[0] + M1.q[3]*M2.q[1],
                             M1.q[0]*M2.q[3] + M1.q[3]*M2.q[4],
                             M1.q[0]*M2.q[6] + M1.q[3]*M2.q[7] + M1.q[6],
                             M1.q[1]*M2.q[0] + M1.q[4]*M2.q[1],
                             M1.q[1]*M2.q[3] + M1.q[4]*M2.q[4],
                             M1.q[1]*M2.q[6] + M1.q[4]*M2.q[7] + M1.q[7]);
    };

    /**
     * Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
    Matrix >::type operator *(const self& M1, const Matrix& M2) {
      if(M2.get_row_count() != 3)
	throw std::range_error("Matrix M's row count is not 3, 2D transformation impossible!");
      Matrix result(M2);
      for(size_type i=0;i<3;++i)
	for(size_type jj=0;jj<result.get_col_count();++jj) {
	  result(i,jj) = 0;
	  for(size_type j=0;j<3;++j)
	    result(i,jj) += M1.q[j*3+i] * M2(j,jj);
	};
      return result;
    };
    
    /**
     * Matrix multiplication.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value, 
    Matrix >::type operator *(const Matrix& M1,const self& M2) {
      if(M1.get_col_count() != 3)
        throw std::range_error("Matrix M1's column count is not 3, 2D transformation impossible!");
      Matrix result(M1.get_row_count(),3);
      for(size_type i=0;i<result.get_row_count();++i)
        for(size_type jj=0;jj<3;++jj) {
	  result(i,jj) = 0;
	  for(size_type j=0;j<3;++j)
            result(i,jj) += M1(i,j) * M2.q[jj*3+j];
	};
      return result;
    };
    

    /**
     * 2D Transformation matrix times a column 2D vector.
     * \test PASSED
     */
    friend
    vect<value_type,2> operator *(const self& M, const vect<value_type,2>& V) {
      return vect<value_type,2>(V[0] * M.q[0] + V[1] * M.q[3] + M.q[6],V[0] * M.q[1] + V[1] * M.q[4] + M.q[7]);
    };

    /**
     * 2D Transformation matrix times a column 2D augmented vector.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& M, const vect<value_type,3>& V) {
      return vect<value_type,3>(V[0] * M.q[0] + V[1] * M.q[3] + V[2] * M.q[6],V[0] * M.q[1] + V[1] * M.q[4] + V[2] * M.q[7], V[2] * M.q[8]);
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Standard equality operator.
     * \test PASSED
     */
    friend
    bool operator ==(const self& M1, const self& M2) {
      return ((M1.q[0] == M2.q[0]) && (M1.q[1] == M2.q[1]) && (M1.q[6] == M2.q[6]) && (M1.q[7] == M2.q[7]));
    };

    /**
     * Standard inequality operator.
     * \test PASSED
     */
    friend
    bool operator !=(const self& M1, const self& M2) {
      return ((M1.q[0] != M2.q[0]) || (M1.q[1] != M2.q[1]) || (M1.q[6] != M2.q[6]) || (M1.q[7] != M2.q[7]));
    };

/*******************************************************************************
                         Special Methods
*******************************************************************************/

    /**
     * Rotate-only a 2D vector.
     * \test PASSED
     */
    vect<value_type,2> rotate(const vect<value_type,2>& V) const {
      return vect<value_type,2>(V[0] * q[0] + V[1] * q[3],V[0] * q[1] + V[1] * q[4]);
    };


/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/

    /**
     * Creates the transpose matrix.
     * \note the matrix is no longer a transformation matrix.
     * \test PASSED
     */
    friend
    mat<value_type,mat_structure::square> transpose(const self& M) {
      return mat<value_type,mat_structure::square>(M.q[0],M.q[1],0.0,M.q[3],M.q[4],0.0,M.q[6],M.q[7],1.0);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& M) {
      return M.q[0] + M.q[4] + value_type(1.0);
    };

    /**
     * Gets the determinant of the matrix.
     * \test PASSED
     */
    friend
    value_type determinant(const self&) {
      return value_type(1.0);
    };

    /**
     * Invert the transformation.
     * \test PASSED
     */
    friend
    self invert(const self& M) {
      return self(M.q[0],M.q[1],-M.q[0]*M.q[6]-M.q[1]*M.q[7],M.q[3],M.q[4],-M.q[3]*M.q[6]-M.q[4]*M.q[7]);
    };

    /**
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      return mat<value_type,mat_structure::symmetric>(q[0],value_type(0.0),value_type(0.5)*q[6],q[0],value_type(0.5)*q[7],value_type(1.0));
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      return mat<value_type,mat_structure::skew_symmetric>(-q[1],value_type(0.5)*q[6],value_type(0.5)*q[7]);
    };


/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_ALIAS("cos",q[0])
        & RK_SERIAL_SAVE_WITH_ALIAS("sin",q[1])
	& RK_SERIAL_SAVE_WITH_ALIAS("t_x",q[6])
	& RK_SERIAL_SAVE_WITH_ALIAS("t_y",q[7]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("cos",q[0])
        & RK_SERIAL_LOAD_WITH_ALIAS("sin",q[1])
	& RK_SERIAL_LOAD_WITH_ALIAS("t_x",q[6])
	& RK_SERIAL_LOAD_WITH_ALIAS("t_y",q[7]);
      q[3] = -q[1];
      q[4] = q[0];
      q[2] = 0.0;
      q[5] = 0.0;
      q[8] = 1.0;
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x00000017,1,"trans_mat_2D",serialization::serializable)

};


/**
 * Prints a 2D rotation matrix to a standard output stream (<<) as "(angle = a; translation = (tx; ty))".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const trans_mat_2D<T>& M) {
  out_stream << "(angle = " << M.getAngle() << "; translation = " << M.getTranslation() << ")";
  return out_stream;
};



template <typename T>
struct is_readable_matrix< trans_mat_2D<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< trans_mat_2D<T> > type;
};









};


#endif













