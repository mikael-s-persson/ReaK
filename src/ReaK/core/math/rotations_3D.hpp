/**
 * \file rotations_3D.hpp
 *
 * This library declares all geometric 3D rotation classes for fixed (2,3) and variable dimensions.
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

#ifndef ROTATIONS_3D_HPP
#define ROTATIONS_3D_HPP

#include "mat_alg.hpp"

#include "vect_concepts.hpp"


namespace ReaK {



//Forward declaration
template <class T>
class quaternion;

template <class T>
class euler_angles_TB;

template <class T>
class axis_angle;

template <class T>
class trans_mat_3D;



/**
 * This class is a rotation matrix 3 by 3.
 * \test All tests for this class have been passed!
 */
template <typename T>
class rot_mat_3D : public serialization::serializable {
  public:
    typedef rot_mat_3D<T> self;
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
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::orthogonal);
    
  private:
    value_type q[9];

    /**
     * Constructor from the components of the rotation matrix.
     * \test PASSED
     */
    explicit rot_mat_3D(const_reference a11, const_reference a12, const_reference a13, 
			const_reference a21, const_reference a22, const_reference a23, 
			const_reference a31, const_reference a32, const_reference a33) {
      q[0] = a11; q[1] = a21; q[2] = a31;
      q[3] = a12; q[4] = a22; q[5] = a32;
      q[6] = a13; q[7] = a23; q[8] = a33;
    };

    /**
     * Assignment to a matrix.
     * \test PASSED
     */
//     self& operator =(const mat_cm<T>& M) {
//       q[0] = M.q[0];
//       q[1] = M.q[1];
//       q[2] = M.q[2];
//       q[3] = M.q[3];
//       q[4] = M.q[4];
//       q[5] = M.q[5];
//       q[6] = M.q[6];
//       q[7] = M.q[7];
//       q[8] = M.q[8];
//       return *this;
//     };

  public:

    friend class quaternion<value_type>;
    friend class euler_angles_TB<value_type>;
    friend class axis_angle<value_type>;
    friend class trans_mat_3D<value_type>;

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor, sets the matrix to identity (no rotation).
     * \test PASSED
     */
    rot_mat_3D() {
      q[0] = 1.0; q[1] = 0.0; q[2] = 0.0;
      q[3] = 0.0; q[4] = 1.0; q[5] = 0.0;
      q[6] = 0.0; q[7] = 0.0; q[8] = 1.0;
    };

    /**
     * Constructor from an array of components.
     * \test PASSED
     */
    rot_mat_3D(const_pointer M) {
      vect<value_type,3> v1 = unit(vect<value_type,3>(M));
      q[0] = v1[0];
      q[1] = v1[1];
      q[2] = v1[2];
      vect<value_type,3> v2(&M[3]);
      v2 = unit(v2 - (v2 * v1) * v1);
      q[3] = v2[0];
      q[4] = v2[1];
      q[5] = v2[2];
      v2 = v1 % v2;
      q[6] = v2[0];
      q[7] = v2[1];
      q[8] = v2[2];
    };
    
    rot_mat_3D(const self& R) {
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = R.q[3];
      q[4] = R.q[4];
      q[5] = R.q[5];
      q[6] = R.q[6];
      q[7] = R.q[7];
      q[8] = R.q[8];
    };
    
    template <typename Matrix>
    explicit rot_mat_3D(const Matrix& M, typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                                                      !boost::is_same<Matrix,self>::value, void* >::type dummy = NULL) {
      if((M.get_col_count() != 3) || (M.get_row_count() != 3))
	throw std::range_error("Right-hand-side of assignment to a 3D rotation matrix is not of dimension 3x3!");
      vect<value_type,3> v1(M(0,0),M(1,0),M(2,0));
      q[0] = v1[0];
      q[1] = v1[1];
      q[2] = v1[2];
      vect<value_type,3> v2(M(0,1),M(1,1),M(2,1));
      v2 = unit(v2 - (v2 * v1) * v1);
      q[3] = v2[0];
      q[4] = v2[1];
      q[5] = v2[2];
      v2 = v1 % v2;
      q[6] = v2[0];
      q[7] = v2[1];
      q[8] = v2[2];
    };

    // Copy constructor. Default is good. \test PASSED
    /**
     * Destructor.
     * \test PASSED
     */
    ~rot_mat_3D() { };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides a copy of the rotation matrix as an ordinary 3x3 matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      return mat<value_type,mat_structure::square>(q[0],q[3],q[6],q[1],q[4],q[7],q[2],q[5],q[8]);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    value_type operator [](size_type i) const {
      if(i >= 9)
	throw std::range_error("Matrix index out of range.");
      return q[i];
    };

    /**
     * Array double-indexing operator, ith row and jth column, accessor for read only.
     * \test PASSED
     */
    value_type operator ()(size_type i,size_type j) const {
      if((i >= 3) || (j >= 3))
	throw std::range_error("Matrix index out of range.");
      return q[j*3+i];
    };
    
    size_type get_row_count() const { return 3; };
    size_type get_col_count() const { return 3; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    self& operator =(const self& M) {
      q[0] = M.q[0];
      q[1] = M.q[1];
      q[2] = M.q[2];
      q[3] = M.q[3];
      q[4] = M.q[4];
      q[5] = M.q[5];
      q[6] = M.q[6];
      q[7] = M.q[7];
      q[8] = M.q[8];
      return *this;
    }; 
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value,
    self& >::type operator =(const Matrix& M) {
      if((M.get_col_count() != 3) || (M.get_row_count() != 3))
	throw std::range_error("Right-hand-side of assignment to a 3D rotation matrix is not of dimension 3x3!");
      vect<value_type,3> v1(M(0,0),M(1,0),M(2,0));
      q[0] = v1[0];
      q[1] = v1[1];
      q[2] = v1[2];
      vect<value_type,3> v2(M(0,1),M(1,1),M(2,1));
      v2 = unit(v2 - (v2 * v1) * v1);
      q[3] = v2[0];
      q[4] = v2[1];
      q[5] = v2[2];
      v2 = v1 % v2;
      q[6] = v2[0];
      q[7] = v2[1];
      q[8] = v2[2];
    };

    /**
     * Assignment operator from a quaternion representation.
     */
    self& operator =(const quaternion<value_type>& Q) {
      return (*this = Q.getRotMat());
    };

    /**
     * Assignment operator from a euler angles TB representation.
     */
    self& operator =(const euler_angles_TB<value_type>& E) {
      return (*this = E.getRotMat());
    };

    /**
     * Assignment operator from an axis / angle representation.
     */
    self& operator =(const axis_angle<value_type>& A) {
      return (*this = A.getRotMat());
    };

    /**
     * Multiplication by a rotation matrix and store.
     * \test PASSED
     */
    self& operator *=(const self& M) {
      *this = self(q[0]*M.q[0] + q[3]*M.q[1] + q[6]*M.q[2],q[0]*M.q[3] + q[3]*M.q[4] + q[6]*M.q[5],q[0]*M.q[6] + q[3]*M.q[7] + q[6]*M.q[8],
                   q[1]*M.q[0] + q[4]*M.q[1] + q[7]*M.q[2],q[1]*M.q[3] + q[4]*M.q[4] + q[7]*M.q[5],q[1]*M.q[6] + q[4]*M.q[7] + q[7]*M.q[8],
                   q[2]*M.q[0] + q[5]*M.q[1] + q[8]*M.q[2],q[2]*M.q[3] + q[5]*M.q[4] + q[8]*M.q[5],q[2]*M.q[6] + q[5]*M.q[7] + q[8]*M.q[8]);
      return *this;
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Multiplication by a rotation matrix.
     * \test PASSED
     */
    friend
    self operator *(const self& M1, const self& M2) {
      return self(M1.q[0]*M2.q[0] + M1.q[3]*M2.q[1] + M1.q[6]*M2.q[2], M1.q[0]*M2.q[3] + M1.q[3]*M2.q[4] + M1.q[6]*M2.q[5], M1.q[0]*M2.q[6] + M1.q[3]*M2.q[7] + M1.q[6]*M2.q[8],
                  M1.q[1]*M2.q[0] + M1.q[4]*M2.q[1] + M1.q[7]*M2.q[2], M1.q[1]*M2.q[3] + M1.q[4]*M2.q[4] + M1.q[7]*M2.q[5], M1.q[1]*M2.q[6] + M1.q[4]*M2.q[7] + M1.q[7]*M2.q[8],
                  M1.q[2]*M2.q[0] + M1.q[5]*M2.q[1] + M1.q[8]*M2.q[2], M1.q[2]*M2.q[3] + M1.q[5]*M2.q[4] + M1.q[8]*M2.q[5], M1.q[2]*M2.q[6] + M1.q[5]*M2.q[7] + M1.q[8]*M2.q[8]);
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
	throw std::range_error("Matrix M's row count is not 3, 3D rotation impossible!");
      Matrix result(M2);
      for(size_type i=0;i<3;++i)
	for(size_type jj=0;jj<result.get_col_count();++jj) {
	  result(i,jj) = 0;
	  for(size_type j=0;j<3;++j)
	    result(i,jj) += M1.q[j*3+i] * M2(j,jj);
	};
      return result;
    };

    template <typename Matrix> 
    friend 
    typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
    Matrix >::type operator *(const Matrix& M1,const self& M2)  {
      if(M1.get_col_count() != 3)
        throw std::range_error("Matrix M1's column count is not 3, 3D rotation impossible!");
      Matrix result(M1);
      for(size_type i=0;i<result.get_row_count();++i)
        for(size_type jj=0;jj<3;++jj) {
          result(i,jj) = 0;
	  for(size_type j=0;j<3;++j)
            result(i,jj) += M1(i,j) * M2.q[jj*3+j];
        };
      return result;
    };

    /**
     * Multiplication with a column vector.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& R, const vect<value_type,3>& V) {
      return vect<value_type,3>( R.q[0]*V[0] + R.q[3]*V[1] + R.q[6]*V[2],
                                 R.q[1]*V[0] + R.q[4]*V[1] + R.q[7]*V[2],
                                 R.q[2]*V[0] + R.q[5]*V[1] + R.q[8]*V[2]);
    };

    friend 
    vect<value_type,3> operator *(const vect<value_type,3>& V, const self& R) {
      return vect<value_type,3>( R.q[0]*V[0] + R.q[1]*V[1] + R.q[2]*V[2],
                                 R.q[3]*V[0] + R.q[4]*V[1] + R.q[5]*V[2],
                                 R.q[6]*V[0] + R.q[7]*V[1] + R.q[8]*V[2]);
    };


/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality operator for a rotation matrix.
     * \test PASSED
     */
    friend
    bool operator ==(const self& M1, const self& M2) {
      return ((M1.q[0] == M2.q[0]) &&
              (M1.q[1] == M2.q[1]) &&
              (M1.q[2] == M2.q[2]) &&
              (M1.q[3] == M2.q[3]) &&
              (M1.q[4] == M2.q[4]) &&
              (M1.q[5] == M2.q[5]));
    };

    /**
     * Inequality operator for a rotation matrix.
     * \test PASSED
     */
    friend
    bool operator !=(const self& M1,const self& M2) {
      return ((M1.q[0] != M2.q[0]) ||
              (M1.q[1] != M2.q[1]) ||
              (M1.q[2] != M2.q[2]) ||
              (M1.q[3] != M2.q[3]) ||
              (M1.q[4] != M2.q[4]) ||
              (M1.q[5] != M2.q[5]));
    };


/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/

    /**
     * Produces a transpose matrix which is the inverse rotation.
     * \test PASSED
     */
    friend
    self transpose(const self& R) {
      return self(R.q[0],R.q[1],R.q[2],R.q[3],R.q[4],R.q[5],R.q[6],R.q[7],R.q[8]);
    };
    
    friend
    self transpose_move(const self& R) {
      return self(R.q[0],R.q[1],R.q[2],R.q[3],R.q[4],R.q[5],R.q[6],R.q[7],R.q[8]);
    };

    /**
     * Produces a cofactor matrix which is the same as the rotation matrix itself.
     * \test PASSED
     */
    friend
    self cofactor(const self& R) {
      return R;
    };

    /**
     * Invert the transformation.
     * \test PASSED
     */
    friend
    self invert(const self& R) {
      return self(R.q[0],R.q[1],R.q[2],R.q[3],R.q[4],R.q[5],R.q[6],R.q[7],R.q[8]);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& R) {
      return R.q[0] + R.q[4] + R.q[8];
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
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      return mat<value_type,mat_structure::symmetric>(*this);
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      return mat<value_type,mat_structure::skew_symmetric>(*this);
    };


/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_ALIAS("r11",q[0])
        & RK_SERIAL_SAVE_WITH_ALIAS("r21",q[1])
	& RK_SERIAL_SAVE_WITH_ALIAS("r31",q[2])
	& RK_SERIAL_SAVE_WITH_ALIAS("r12",q[3])
	& RK_SERIAL_SAVE_WITH_ALIAS("r22",q[4])
	& RK_SERIAL_SAVE_WITH_ALIAS("r32",q[5])
	& RK_SERIAL_SAVE_WITH_ALIAS("r13",q[6])
	& RK_SERIAL_SAVE_WITH_ALIAS("r23",q[7])
	& RK_SERIAL_SAVE_WITH_ALIAS("r33",q[8]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("r11",q[0])
        & RK_SERIAL_LOAD_WITH_ALIAS("r21",q[1])
	& RK_SERIAL_LOAD_WITH_ALIAS("r31",q[2])
	& RK_SERIAL_LOAD_WITH_ALIAS("r12",q[3])
	& RK_SERIAL_LOAD_WITH_ALIAS("r22",q[4])
	& RK_SERIAL_LOAD_WITH_ALIAS("r32",q[5])
	& RK_SERIAL_LOAD_WITH_ALIAS("r13",q[6])
	& RK_SERIAL_LOAD_WITH_ALIAS("r23",q[7])
	& RK_SERIAL_LOAD_WITH_ALIAS("r33",q[8]);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x00000018,1,"rot_mat_3D",serialization::serializable)

};

/**
 * Prints a rotation matrix to a standard output stream (<<) as
 * "((a11; a12; a13); (a21; a22; a23); (a31; a32; a33))".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const rot_mat_3D<T>& R) {
  return out_stream << "((" << R(0,0) << "; " << R(0,1) << "; " << R(0,2) << "); ("
                            << R(1,0) << "; " << R(1,1) << "; " << R(1,2) << "); ("
                            << R(2,0) << "; " << R(2,1) << "; " << R(2,2) << "))";
};


template <typename T>
struct is_readable_matrix< rot_mat_3D<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< rot_mat_3D<T> > type;
};







/**
 * This class represents a rotation using quaternions (or Euler-Rodriguez parameters).
 * The convention used is with the leading scalar.
 * \test All tests for this class have been passed!
 */
template <typename T>
class quaternion : public serialization::serializable {
  public:
    typedef quaternion<T> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
  private:
    value_type q[4];

    quaternion(const_reference q0, const_reference q1, const_reference q2, const_reference q3) { q[0] = q0; q[1] = q1; q[2] = q2; q[3] = q3; return; };
  public:

    friend class euler_angles_TB<value_type>;
    friend class axis_angle<value_type>;
    friend class trans_mat_3D<value_type>;

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default Constructor.
     * \test PASSED
     */
    quaternion() { q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0; };

    /**
     * Copy-constructor.
     * \test PASSED
     */
    quaternion(const self& Q) { q[0] = Q.q[0]; q[1] = Q.q[1]; q[2] = Q.q[2]; q[3] = Q.q[3]; };
    
    template <typename Vector>
    explicit quaternion(const Vector& aV, typename boost::enable_if_c< is_readable_vector<Vector>::value, void* >::type dummy = NULL) { RK_UNUSED(dummy); 
      vect<value_type,4> v = unit(vect<value_type,4>(aV[0],aV[1],aV[2],aV[3])); 
      q[0] = v[0]; q[1] = v[1]; q[2] = v[2]; q[3] = v[3]; 
    };
    
    /**
     * Constructor from a rotation matrix.
     * \test PASSED
     */
    explicit quaternion(const rot_mat_3D<value_type>& R) {
      using std::sqrt;
      value_type tra = R.q[0] + R.q[4] + R.q[8];
      if (tra > 0.01) {
        q[0] = value_type(0.5)*sqrt(value_type(1.0) + tra);
        q[1] = value_type(0.25)*(R.q[5] - R.q[7])/q[0];
        q[2] = value_type(0.25)*(R.q[6] - R.q[2])/q[0];
        q[3] = value_type(0.25)*(R.q[1] - R.q[3])/q[0];
      } else if ((R.q[0] > R.q[4]) && (R.q[0] > R.q[8])) {
        q[1] = value_type(0.5)*sqrt(value_type(1.0) + R.q[0] - R.q[4] - R.q[8]);
        q[0] = value_type(0.25)*(R.q[7] - R.q[5])/q[1];
        q[2] = value_type(0.25)*(R.q[3] + R.q[1])/q[1];
        q[3] = value_type(0.25)*(R.q[6] + R.q[2])/q[1];
      } else if (R.q[4] > R.q[8]) {
        q[2] = value_type(0.5)*sqrt(value_type(1.0) + R.q[4] - R.q[0] - R.q[8]);
        q[0] = value_type(0.25)*(R.q[6] - R.q[2])/q[2];
        q[1] = value_type(0.25)*(R.q[3] + R.q[1])/q[2];
        q[3] = value_type(0.25)*(R.q[7] + R.q[5])/q[2];
      } else {
        q[3] = value_type(0.5)*sqrt(value_type(1.0) + R.q[8] - R.q[0] - R.q[4]);
        q[0] = value_type(0.25)*(R.q[3] - R.q[1])/q[3];
        q[1] = value_type(0.25)*(R.q[6] + R.q[2])/q[3];
        q[2] = value_type(0.25)*(R.q[7] + R.q[5])/q[3];
      };
      return;
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~quaternion() { };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides the rotation matrix as an ordinary 3x3 matrix.
     * \test PASSED
     */
    mat<value_type, mat_structure::square> getMat() const {
      value_type t01(value_type(2.0)*q[0]*q[1]);
      value_type t02(value_type(2.0)*q[0]*q[2]);
      value_type t03(value_type(2.0)*q[0]*q[3]);
      value_type t11(value_type(2.0)*q[1]*q[1]);
      value_type t12(value_type(2.0)*q[1]*q[2]);
      value_type t13(value_type(2.0)*q[1]*q[3]);
      value_type t22(value_type(2.0)*q[2]*q[2]);
      value_type t23(value_type(2.0)*q[2]*q[3]);
      value_type t33(value_type(2.0)*q[3]*q[3]);
      return mat<value_type, mat_structure::square>(value_type(1.0) - t22 - t33,       t12 - t03,       t02 + t13,
						          t12 + t03, value_type(1.0) - t11 - t33,       t23 - t01,
						          t13 - t02,       t01 + t23, value_type(1.0) - t11 - t22);
    };

    /**
     * Provides the rotation matrix corresponding to the quaternion.
     * \test PASSED
     */
    rot_mat_3D<value_type> getRotMat() const {
      value_type t01(value_type(2.0)*q[0]*q[1]);
      value_type t02(value_type(2.0)*q[0]*q[2]);
      value_type t03(value_type(2.0)*q[0]*q[3]);
      value_type t11(value_type(2.0)*q[1]*q[1]);
      value_type t12(value_type(2.0)*q[1]*q[2]);
      value_type t13(value_type(2.0)*q[1]*q[3]);
      value_type t22(value_type(2.0)*q[2]*q[2]);
      value_type t23(value_type(2.0)*q[2]*q[3]);
      value_type t33(value_type(2.0)*q[3]*q[3]);
      return rot_mat_3D<value_type>(value_type(1.0) - t22 - t33,       t12 - t03,       t02 + t13,
			                  t12 + t03, value_type(1.0) - t11 - t33,       t23 - t01,
			                  t13 - t02,       t01 + t23, value_type(1.0) - t11 - t22);
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      if(i >= 4)
	throw std::range_error("Matrix index out of range.");
      return q[i];
    };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Assignment operator from a quaternion.
     * \test PASSED
     */
    self& operator =(const self& Q) {
      q[0] = Q.q[0];
      q[1] = Q.q[1];
      q[2] = Q.q[2];
      q[3] = Q.q[3];
      return *this;
    };

    /**
     * Assignment operator from a rotation matrix.
     * \test PASSED
     */
    self& operator =(const rot_mat_3D<value_type>& R) {
      return (*this = self(R));
    };

    /**
     * Assignment operator from a euler angles TB representation.
     * \test PASSED
     */
    self& operator =(const euler_angles_TB<value_type>& E) {
      return (*this = E.getQuaternion());
    };

    /**
     * Assignment operator from an axis / angle representation.
     * \test PASSED
     */
    self& operator =(const axis_angle<value_type>& A) {
      return (*this = A.getQuaternion());
    };

    /**
     * Multiply-and-store operator from a quaternion.
     * \test PASSED
     */
    self& operator *=(const self& Q) {
      return (*this = *this * Q);
    };
    
    /**
     * Multiply-and-store operator from a rotation matrix.
     * \test PASSED
     */
    self& operator *=(const rot_mat_3D<value_type>& R) {
      return (*this *= self(R));
    };

    /**
     * Multiply-and-store operator from a euler angles TB representation.
     * \test PASSED
     */
    self& operator *=(const euler_angles_TB<value_type>& E) {
      return (*this *= E.getQuaternion());
    };

    /**
     * Multiply-and-store operator from an axis / angle representation.
     * \test PASSED
     */
    self& operator *=(const axis_angle<value_type>& A) {
      return (*this *= A.getQuaternion());
    };


/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Multiplication by a quaternion.
     * \test PASSED
     */
    friend
    self operator *(const self& Q1, const self& Q2) {
      return self(Q2.q[0]*Q1.q[0]-Q2.q[1]*Q1.q[1]-Q2.q[2]*Q1.q[2]-Q2.q[3]*Q1.q[3],
                  Q2.q[0]*Q1.q[1]+Q2.q[3]*Q1.q[2]-Q2.q[2]*Q1.q[3]+Q2.q[1]*Q1.q[0],
                  Q2.q[0]*Q1.q[2]-Q2.q[3]*Q1.q[1]+Q2.q[1]*Q1.q[3]+Q2.q[2]*Q1.q[0],
                  Q2.q[0]*Q1.q[3]+Q2.q[2]*Q1.q[1]-Q2.q[1]*Q1.q[2]+Q2.q[3]*Q1.q[0]);
    };

    /**
     * Multiplication by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const self& Q, const Matrix& M) {
      return Q.getRotMat() * M;
    };
    
    /**
     * Multiplication by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const Matrix& M, const self& Q) {
      return M * Q.getRotMat();
    };

    /**
     * Multiplication by a column vector.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& Q, const vect<value_type,3>& V) {
      value_type t[9];
      t[0] =   Q.q[0]*Q.q[1];
      t[1] =   Q.q[0]*Q.q[2];
      t[2] =   Q.q[0]*Q.q[3];
      t[3] =  -Q.q[1]*Q.q[1];
      t[4] =   Q.q[1]*Q.q[2];
      t[5] =   Q.q[1]*Q.q[3];
      t[6] =  -Q.q[2]*Q.q[2];
      t[7] =   Q.q[2]*Q.q[3];
      t[8] =  -Q.q[3]*Q.q[3];
      return vect<T,3>( value_type(2.0)*( (t[6] + t[8])*V[0] + (t[4] - t[2])*V[1] + (t[1] + t[5])*V[2] ) + V[0],
                        value_type(2.0)*( (t[2] + t[4])*V[0] + (t[3] + t[8])*V[1] + (t[7] - t[0])*V[2] ) + V[1],
                        value_type(2.0)*( (t[5] - t[1])*V[0] + (t[0] + t[7])*V[1] + (t[3] + t[6])*V[2] ) + V[2]);
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality operator with a quaternion representation.
     * \test PASSED
     */
    friend
    bool operator ==(const self& Q1, const self& Q2) {
      return ((Q1.q[0] == Q2.q[0]) && (Q1.q[1] == Q2.q[1]) && (Q1.q[2] == Q2.q[2]) && (Q1.q[3] == Q2.q[3]));
    };

    /**
     * Inequality operator with a quaternion representation.
     * \test PASSED
     */
    friend
    bool operator !=(const self& Q1, const self& Q2) {
      return ((Q1.q[0] != Q2.q[0]) || (Q1.q[1] != Q2.q[1]) || (Q1.q[2] != Q2.q[2]) || (Q1.q[3] != Q2.q[3]));
    };



/*******************************************************************************
                         Special Methods
*******************************************************************************/

    /**
     * Gets the time-derivative of the quaternion that corresponds to the angular velocity Omega.
     * \test PASSED
     */
    vect<value_type,4> getQuaternionDot(const vect<value_type,3>& Omega) const {
      return vect<value_type,4>(-value_type(0.5)*(q[1]*Omega.q[0] + q[2]*Omega.q[1] + q[3]*Omega.q[2]),
				 value_type(0.5)*(q[0]*Omega.q[0] - q[3]*Omega.q[1] + q[2]*Omega.q[2]),
				 value_type(0.5)*(q[0]*Omega.q[1] + q[3]*Omega.q[0] - q[1]*Omega.q[2]),
				 value_type(0.5)*(q[0]*Omega.q[2] - q[2]*Omega.q[0] + q[1]*Omega.q[1]));
    };

    /**
     * Gets the angular velocity that corresponds to the time-derivative of the quaternion.
     * \test PASSED
     */
    vect<value_type,3> getOmega(const vect<value_type,4>& QuaternionDot) const {
      return vect<value_type,3>(value_type(2.0)*(-q[1]*QuaternionDot.q[0] + q[0]*QuaternionDot.q[1] + q[3]*QuaternionDot.q[2] - q[2]*QuaternionDot.q[3]),
				value_type(2.0)*(-q[2]*QuaternionDot.q[0] - q[3]*QuaternionDot.q[1] + q[0]*QuaternionDot.q[2] + q[1]*QuaternionDot.q[3]),
				value_type(2.0)*(-q[3]*QuaternionDot.q[0] + q[2]*QuaternionDot.q[1] - q[1]*QuaternionDot.q[2] + q[0]*QuaternionDot.q[3]));
    };

    /**
     * Gets the 2-time-derivative of the quaternion that corresponds to the angular velocity Omega.
     * \test PASSED
     */
    vect<value_type,4> getQuaternionDotDot(const vect<value_type,4>& QD, const vect<value_type,3>& W, const vect<value_type,3>& WD) const {
      return vect<value_type,4>(-value_type(0.5)*( q[1]*WD.q[0] + q[2]*WD.q[1] + q[3]*WD.q[2] + QD.q[1]*W.q[0] + QD.q[2]*W.q[1] + QD.q[3]*W.q[2]),
				 value_type(0.5)*( q[0]*WD.q[0] - q[3]*WD.q[1] + q[2]*WD.q[2] + QD.q[0]*W.q[0] - QD.q[3]*W.q[1] + QD.q[2]*W.q[2]),
				 value_type(0.5)*( q[3]*WD.q[0] + q[0]*WD.q[1] - q[1]*WD.q[2] + QD.q[3]*W.q[0] + QD.q[0]*W.q[1] - QD.q[1]*W.q[2]),
				 value_type(0.5)*(-q[2]*WD.q[0] + q[1]*WD.q[1] + q[0]*WD.q[2] - QD.q[2]*W.q[0] + QD.q[1]*W.q[1] + QD.q[0]*W.q[2]));
    };

    /**
     * Gets the angular acceleration that corresponds to the 2-time-derivative of the quaternion.
     * \test PASSED
     */
    vect<value_type,3> getOmegaDot(const vect<value_type,4>& QD, const vect<value_type,4>& QDD) const {
      return vect<value_type,3>(value_type(2.0)*(-q[1]*QDD.q[0] + q[0]*QDD.q[1] + q[3]*QDD.q[2] - q[2]*QDD.q[3] - QD.q[1]*QD.q[0] + QD.q[0]*QD.q[1] + QD.q[3]*QD.q[2] - QD.q[2]*QD.q[3]),
				value_type(2.0)*(-q[2]*QDD.q[0] - q[3]*QDD.q[1] + q[0]*QDD.q[2] + q[1]*QDD.q[3] - QD.q[2]*QD.q[0] - QD.q[3]*QD.q[1] + QD.q[0]*QD.q[2] + QD.q[1]*QD.q[3]),
				value_type(2.0)*(-q[3]*QDD.q[0] + q[2]*QDD.q[1] - q[1]*QDD.q[2] + q[0]*QDD.q[3] - QD.q[3]*QD.q[0] + QD.q[2]*QD.q[1] - QD.q[1]*QD.q[2] + QD.q[0]*QD.q[3]));
    };

/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/

    /**
     * Produces a transpose quaternion which is the inverse rotation.
     * \test PASSED
     */
    friend
    self transpose(const self& Q) {
      return self(Q.q[0],-Q.q[1],-Q.q[2],-Q.q[3]);
    };
    
    /**
     * Produces a transpose quaternion which is the inverse rotation.
     * \test PASSED
     */
    friend
    self transpose_move(const self& Q) {
      return self(Q.q[0],-Q.q[1],-Q.q[2],-Q.q[3]);
    };

    /**
     * Produces a cofactor matrix which is the same as the rotation matrix itself.
     * \test PASSED
     */
    friend
    self cofactor(const self& Q) {
      return Q;
    };

    /**
     * Invert the rotation.
     * \test PASSED
     */
    friend
    self invert(const self& Q) {
      return self(Q.q[0],-Q.q[1],-Q.q[2],-Q.q[3]);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& Q) {
      return value_type(4.0)*Q.q[0]*Q.q[0] - value_type(1.0);
    };

    /**
     * Gets the determinant of the matrix.
     * \test PASSED
     */
    friend
    value_type determinant(const self& Q) {
      return value_type(1.0);
    };

    /**
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      value_type t11(value_type(2.0)*q[1]*q[1]);
      value_type t12(value_type(2.0)*q[1]*q[2]);
      value_type t13(value_type(2.0)*q[1]*q[3]);
      value_type t22(value_type(2.0)*q[2]*q[2]);
      value_type t23(value_type(2.0)*q[2]*q[3]);
      value_type t33(value_type(2.0)*q[3]*q[3]);
      return mat<value_type,mat_structure::symmetric>(value_type(1.0) - t22 - t33,             t12,             t13,
			                                               value_type(1.0) - t11 - t33,             t23,
			                                                                value_type(1.0) - t11 - t22);
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      value_type t01(value_type(2.0)*q[0]*q[1]);
      value_type t02(value_type(2.0)*q[0]*q[2]);
      value_type t03(value_type(2.0)*q[0]*q[3]);
      return mat<value_type,mat_structure::skew_symmetric>( -t03, t02, -t01);
    };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(q[0])
        & RK_SERIAL_SAVE_WITH_NAME(q[1])
	& RK_SERIAL_SAVE_WITH_NAME(q[2])
	& RK_SERIAL_SAVE_WITH_NAME(q[3]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(q[0])
        & RK_SERIAL_LOAD_WITH_NAME(q[1])
	& RK_SERIAL_LOAD_WITH_NAME(q[2])
	& RK_SERIAL_LOAD_WITH_NAME(q[3]);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x0000001A,1,"quaternion",serialization::serializable)

};


/**
 * Prints a quaternion to a standard output stream (<<) as "(q0; q1; q2; q3)".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const quaternion<T>& Q) {
  return (out_stream << "(" << Q[0] << "; " << Q[1] << "; " << Q[2] << "; " << Q[3] << ")");
};








/**
 * This class repressents a rotation using Euler angles (Tait-Bryan), 321-body-fixed, in body frame.
 */
template <class T>
class euler_angles_TB : public serialization::serializable {
  public:
    typedef euler_angles_TB<T> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  private:
    value_type q[3];

  public:
    friend class axis_angle<value_type>;
    friend class trans_mat_3D<value_type>;

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default Constructor.
     * \test PASSED
     */
    euler_angles_TB() {
      q[0] = 0.0;
      q[1] = 0.0;
      q[2] = 0.0;
    };

    /**
     * Constructor from three euler angles.
     * \test PASSED
     */
    euler_angles_TB(const_reference Yaw_,const_reference Pitch_,const_reference Roll_) {
      q[0] = Yaw_;
      q[1] = Pitch_;
      q[2] = Roll_;
    };

    /**
     * Copy-constructor.
     * \test PASSED
     */
    euler_angles_TB(const self& E) {
      q[0] = E.q[0];
      q[1] = E.q[1];
      q[2] = E.q[2];
    };

    /**
     * Constructor from a quaternion.
     * \test PASSED
     */
    explicit euler_angles_TB(const quaternion<value_type>& Q) {
      using std::asin;
      using std::atan2;
      using std::cos;
      q[1] = value_type(2.0)*(Q.q[0]*Q.q[2] - Q.q[1]*Q.q[3]);
      if((q[1] != value_type(1.0)) && (q[1] != value_type(-1.0))) {
        q[1] = asin(q[1]);
        value_type cp = value_type(1.0)/cos(q[1]);
        q[2] = atan2(value_type(2.0)*cp*(Q.q[2]*Q.q[3] + Q.q[0]*Q.q[1]),cp*(value_type(1.0) - value_type(2.0)*(Q.q[1]*Q.q[1] + Q.q[2]*Q.q[2])));
        q[0] = atan2(value_type(2.0)*cp*(Q.q[1]*Q.q[2] + Q.q[0]*Q.q[3]),cp*(value_type(1.0) - value_type(2.0)*(Q.q[2]*Q.q[2] + Q.q[3]*Q.q[3])));
      } else {
        q[0] = value_type(0.0);
        q[2] = atan2(q[1]*value_type(2.0)*(Q.q[1]*Q.q[2] - Q.q[0]*Q.q[3]),q[1]*value_type(2.0)*(Q.q[1]*Q.q[3] + Q.q[0]*Q.q[2]));
        q[1] *= value_type(1.57079632679489662);
      };
      return;
    };

    /**
     * Constructor from a rotation matrix.
     * \test PASSED
     */
    explicit euler_angles_TB(const rot_mat_3D<value_type>& R) {
      using std::asin;
      using std::atan2;
      using std::cos;
      if((R.q[2] != value_type(1.0)) && (R.q[2] != value_type(-1.0))) {
        q[1] = asin(-R.q[2]);
        value_type cp = value_type(1.0)/cos(q[1]);
        q[2] = atan2(cp*R.q[5],cp*R.q[8]);
        q[0] = atan2(cp*R.q[1],cp*R.q[0]);
      } else {
        q[0] = value_type(0.0);
        q[2] = atan2(-R.q[2]*R.q[3],-R.q[2]*R.q[6]);
        q[1] = -R.q[2]*value_type(1.57079632679489662);
      };
      return;
    };

    euler_angles_TB(const rot_mat_3D<value_type>& R,const euler_angles_TB<value_type>& Predicted) {
      using std::asin;
      using std::atan2;
      using std::cos;
      if((R.q[2] != value_type(1.0)) && (R.q[2] != value_type(-1.0))) {
        q[1] = asin(-R.q[2]);
        value_type cp = value_type(1.0)/cos(q[1]);
        q[2] = atan2(cp*R.q[5],cp*R.q[8]);
        q[0] = atan2(cp*R.q[1],cp*R.q[0]);
      } else {
        q[0] = value_type(0.0);
        q[2] = atan2(-R.q[2]*R.q[3],-R.q[2]*R.q[6]);
        q[1] = -R.q[2]*value_type(1.57079632679489662);
      };
      return;

      //UNTESTED CODE: ...
      /*T s, c;
      T d1, d2, d3, d4, d5, m1;

      vect<T,3> Turns(floor((Predicted.q[0] + M_PI) / T(2.0*M_PI)), floor((Predicted.q[1] + M_PI) / T(2.0*M_PI)), floor((Predicted.q[2] + M_PI) / T(2.0*M_PI)));
      vect<T,3> Pred = vect<T,3>( Predicted.q[0] - Turns.q[0]*T(2.0*M_PI), Predicted.q[1] - Turns.q[1]*T(2.0*M_PI), Predicted.q[2] - Turns.q[2]*T(2.0*M_PI));

      vect<T,3> Result(0.0,asin(-m[6]),0.0);
      c = cos(Result.q[1]);

      if(c != 0.0) {
        Result.q[0] = atan2(m[3],m[0]);
        Result.q[2] = atan2(m[7],m[8]);
        d1 = fabs(Result.q[2] - Predicted.q[2]) + fabs(Result.q[0] - Predicted.q[0]);
        d2 = fabs(Result.q[2] + T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] + T(M_PI) - Predicted.q[0]);
        d3 = fabs(Result.q[2] - T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] - T(M_PI) - Predicted.q[0]);
        d4 = fabs(Result.q[2] - T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] + T(M_PI) - Predicted.q[0]);
        d5 = fabs(Result.q[2] + T(M_PI) - Predicted.q[2]) + fabs(Result.q[0] - T(M_PI) - Predicted.q[0]);
        m1 = MIN(MIN(d1,d2),MIN(d3,MIN(d4,d5)));
        if(m1 == d5)
          Result = TVect3<T>(Result.q[0] - T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) + Turns.q[1]*T(2.0*M_PI), Result.q[2] + T(M_PI) + Turns.q[2]*T(2.0*M_PI));
        else if(m1 == d4)
          Result = TVect3<T>(Result.q[0] + T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) + Turns.q[1]*T(2.0*M_PI), Result.q[2] - T(M_PI) + Turns.q[2]*T(2.0*M_PI));
        else if(m1 == d3)
          Result = TVect3<T>(Result.q[0] - T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) + Turns.q[1]*T(2.0*M_PI), Result.q[2] - T(M_PI) + Turns.q[2]*T(2.0*M_PI));
        else if(m1 == d2)
          Result = TVect3<T>(Result.q[0] + T(M_PI) + Turns.q[0]*T(2.0*M_PI), -Result.q[1] + T(M_PI) + Turns.q[1]*T(2.0*M_PI), Result.q[2] + T(M_PI) + Turns.q[2]*T(2.0*M_PI));
        else
          Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] + Turns.q[2]*T(2.0*M_PI));
      } else {
        //Singularity Arises here because Roll becomes equivalent to yaw.
        //Here, one of them (Yaw) is set to the prediction.
        Result.q[0] = Predicted.q[0];
        s = sin(Result.q[0]);
        c = cos(Result.q[0]);
        Result.q[2] = atan2(m[2]*s-m[5]*c,m[2]*c+m[5]*s);
        d1 = fabs(Result.q[2] - Predicted.q[2]);
        d2 = fabs(Result.q[2] + T(2.0*M_PI) - Predicted.q[2]);
        d3 = fabs(Result.q[2] - T(2.0*M_PI) - Predicted.q[2]);
        m1 = MIN(d1,MIN(d2,d3));
        if(m1 == d3)
          Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] - T(2.0*M_PI) + Turns.q[2]*T(2.0*M_PI));
        else if(m1 == d2)
          Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] + T(2.0*M_PI) + Turns.q[2]*T(2.0*M_PI));
        else
          Result = TVect3<T>(Result.q[0] + Turns.q[0]*T(2.0*M_PI), Result.q[1] + Turns.q[1]*T(2.0*M_PI), Result.q[2] + Turns.q[2]*T(2.0*M_PI));
      };
      return Result;*/
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~euler_angles_TB() { };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Get yaw, read-write.
     * \test PASSED
     */
    reference yaw() { return q[0]; };
    /**
     * Get pitch, read-write.
     * \test PASSED
     */
    reference pitch() { return q[1]; };
    /**
     * Get roll, read-write.
     * \test PASSED
     */
    reference roll() { return q[2]; };

    /**
     * Get yaw, read-only.
     * \test PASSED
     */
    const_reference yaw() const { return q[0]; };
    /**
     * Get pitch, read-only.
     * \test PASSED
     */
    const_reference pitch() const { return q[1]; };
    /**
     * Get roll, read-only.
     * \test PASSED
     */
    const_reference roll() const { return q[2]; };

    /**
     * Provides a quaternion corresponding to this rotation.
     * \test PASSED
     */
    quaternion<value_type> getQuaternion() const {
      using std::cos;
      using std::sin;
      value_type cpsi = cos(value_type(0.5)*q[0]);
      value_type spsi = sin(value_type(0.5)*q[0]);
      value_type ctheta = cos(value_type(0.5)*q[1]);

      value_type stheta = sin(value_type(0.5)*q[1]);
      value_type cphi = cos(value_type(0.5)*q[2]);
      value_type sphi = sin(value_type(0.5)*q[2]);

      return quaternion<value_type>(cphi*ctheta*cpsi + sphi*stheta*spsi,
				    sphi*ctheta*cpsi - cphi*stheta*spsi,
				    cphi*stheta*cpsi + sphi*ctheta*spsi,
				    cphi*ctheta*spsi - sphi*stheta*cpsi);
    };

    /**
     * Provides a rotation matrix corresponding to this rotation.
     * \test PASSED
     */
    rot_mat_3D<value_type> getRotMat() const {
      using std::sin;
      using std::cos;
      value_type s1(sin(q[0]));
      value_type c1(cos(q[0]));
      value_type s2(sin(q[1]));
      value_type c2(cos(q[1]));
      value_type s3(sin(q[2]));
      value_type c3(cos(q[2]));

      return rot_mat_3D<value_type>(c1 * c2,-(s1 * c3) + (c1 * s2 * s3), (s1 * s3) + (c1 * s2 * c3),
				    s1 * c2, (c1 * c3) + (s1 * s2 * s3),-(c1 * s3) + (s1 * s2 * c3),
				    -s2    ,                    c2 * s3,                    c2 * c3);
    };

    /**
     * Provides a rotation matrix as a regular 3x3 matris corresponding to this rotation.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      using std::sin;
      using std::cos;
      value_type s1(sin(q[0]));
      value_type c1(cos(q[0]));
      value_type s2(sin(q[1]));
      value_type c2(cos(q[1]));
      value_type s3(sin(q[2]));
      value_type c3(cos(q[2]));

      return mat<value_type,mat_structure::square>( c1 * c2,-(s1 * c3) + (c1 * s2 * s3), (s1 * s3) + (c1 * s2 * c3),
						    s1 * c2, (c1 * c3) + (s1 * s2 * s3),-(c1 * s3) + (s1 * s2 * c3),
						    -s2    ,                    c2 * s3,                    c2 * c3);
    };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard assignment operator.
     * \test PASSED
     */
    self& operator =(const self& E) {
      q[0] = E.q[0];
      q[1] = E.q[1];
      q[2] = E.q[2];
      return *this;
    };

    /**
     * Assignment from a quaternion.
     * \test PASSED
     */
    self& operator =(const quaternion<value_type>& Q) {
      return *this = self(Q);
    };

    /**
     * Assignment from a rotation matrix.
     * \test PASSED
     */
    self& operator =(const rot_mat_3D<value_type>& R) {
      return *this = self(R);
    };

    /**
     * Assignment from an axis / angle representation.
     * \test PASSED
     */
    self& operator =(const axis_angle<T>& A) {
      return (*this = A.getEulerAnglesTB());
    };

    /**
     * Multiply-and-store from a euler angles.
     * \test PASSED
     */
    self& operator *=(const self& E) {
      return (*this = (this->getRotMat() * E.getRotMat()));
    };

    /**
     * Multiply-and-store from a rotation matrix.
     * \test PASSED
     */
    self& operator *=(const rot_mat_3D<value_type>& R) {
      return (*this = (this->getRotMat() * R));
    };

    /**
     * Multiply-and-store from a quaternion.
     * \test PASSED
     */
    self& operator *=(const quaternion<value_type>& Q) {
      return (*this = (this->getQuaternion() * Q));
    };
    
    /**
     * Assignment from an axis / angle representation.
     * \test PASSED
     */
    self& operator *=(const axis_angle<T>& A) {
      return (*this = this->getRotMat() * A.getRotMat());
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Multiply by a euler angle representation.
     * \test PASSED
     */
    friend
    rot_mat_3D<value_type> operator *(const self& E1, const self& E2) {
      return E1.getRotMat() * E2.getRotMat();
    };

    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const self& E, const Matrix& M) {
      return E.getRotMat() * M;
    };
    
    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const Matrix& M, const self& E) {
      return M * E.getRotMat();
    };

    /**
     * Multiply by a vector, rotating it.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& E, const vect<value_type,3>& V) {
      return E.getRotMat() * V;
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality comparison operator with a euler angle representation.
     * \test PASSED
     */
    friend
    bool operator ==(const self& E1, const self& E2) {
      return ((E1.q[0] == E2.q[0]) && (E1.q[1] == E2.q[1]) && (E1.q[2] == E2.q[2]));
    };

    /**
     * Inequality comparison operator with a euler angle representation.
     * \test PASSED
     */
    friend
    bool operator !=(const self& E1, const self& E2) {
      return ((E1.q[0] != E2.q[0]) || (E1.q[1] != E2.q[1]) || (E1.q[2] != E2.q[2]));
    };



/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/
    /**
     * Produces a transpose quaternion which is the inverse rotation.
     * \test PASSED
     */
    friend 
    self transpose(const self& E) {
      using std::sin;
      using std::cos;
      using std::asin;
      using std::atan2;
      
      value_type s1(sin(E.q[0]));
      value_type c1(cos(E.q[0]));
      value_type s2(sin(E.q[1]));
      value_type c2(cos(E.q[1]));
      value_type s3(sin(E.q[2]));
      value_type c3(cos(E.q[2]));

      value_type R2((s1 * s3) + (c1 * s2 * c3));

      self result;

      if((R2 != value_type(1.0)) && (R2 != value_type(-1.0))) {
	value_type R0(c1 * c2);
	value_type R1(-(s1 * c3) + (c1 * s2 * s3));
	value_type R5(-(c1 * s3) + (s1 * s2 * c3));
	value_type R8(c2 * c3);
        result.q[1] = asin(-R2);
        value_type cp = value_type(1.0)/cos(result.q[1]);
        result.q[2] = atan2(cp*R5,cp*R8);
        result.q[0] = atan2(cp*R1,cp*R0);
      } else {
	value_type R3(s1 * c2);
	result.q[0] = value_type(0.0);
        result.q[2] = atan2(-R2*R3,R2*s2);
        result.q[1] = -R2*value_type(1.57079632679489662);
      };
      return result;
    };

    /**
     * Produces a cofactor matrix which is the same as the rotation matrix itself.
     * \test PASSED
     */
    friend
    self cofactor(const self& E) {
      return E;
    };

    /**
     * Invert the rotation.
     * \test PASSED
     */
    friend
    self invert(const self& E) {
      return transpose(E);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& E) {
      using std::sin;
      using std::cos;
      value_type t = cos(value_type(0.5)*E.q[2])*cos(value_type(0.5)*E.q[1])*cos(value_type(0.5)*E.q[0]) + sin(value_type(0.5)*E.q[2])*sin(value_type(0.5)*E.q[1])*sin(value_type(0.5)*E.q[0]);
      return value_type(4.0)*t*t - value_type(1.0);
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
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      return mat<value_type,mat_structure::symmetric>(this->getMat());
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      return mat<value_type,mat_structure::skew_symmetric>(this->getMat());
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(q[0])
        & RK_SERIAL_SAVE_WITH_NAME(q[1])
	& RK_SERIAL_SAVE_WITH_NAME(q[2]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(q[0])
        & RK_SERIAL_LOAD_WITH_NAME(q[1])
	& RK_SERIAL_LOAD_WITH_NAME(q[2]);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x0000001B,1,"euler_angles_TB",serialization::serializable)

};


/**
 * Prints a euler angles to a standard output stream (<<) as "(Yaw = value; Pitch = value; Roll = value)".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const euler_angles_TB<T>& E) {
  return out_stream << "(Yaw = " << E.yaw() << "; Pitch = " << E.pitch() << "; Roll = " << E.roll() << ")";
};







/**
 * This class is a 3D rotation represented by an axis and angle.
 * \test All tests for this class have been passed!
 */
template <class T>
class axis_angle : public serialization::serializable {
  public:
    typedef axis_angle<T> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  private:
    value_type mAngle;
    vect<value_type,3> mAxis;

  public:

    friend class trans_mat_3D<value_type>;

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default Constructor.
     * \test PASSED
     */
    axis_angle() : mAngle(0.0), mAxis(1.0,0.0,0.0) { };

    /**
     * Constructor from angle and axis.
     * \test PASSED
     */
    axis_angle(const value_type& aAngle,const vect<value_type,3>& aAxis) : mAngle(aAngle), mAxis(unit(aAxis)) { };

    /**
     * Copy-constructor.
     * \test PASSED
     */
    axis_angle(const self& A) : mAngle(A.mAngle), mAxis(A.mAxis) { };
    
    /**
     * Constructor from a quaternion.
     * \test PASSED
     */
    explicit axis_angle(const quaternion<value_type>& Q) {
      using std::sqrt;
      using std::acos;
      value_type tmp(sqrt(Q.q[1]*Q.q[1] + Q.q[2]*Q.q[2] + Q.q[3]*Q.q[3]));
      if(tmp > value_type(0.0000001)) {
	mAxis.q[0] = Q.q[1] / tmp;
	mAxis.q[1] = Q.q[2] / tmp;
	mAxis.q[2] = Q.q[3] / tmp;
	mAngle = value_type(2.0) * acos(Q.q[0]);
      } else {
	mAxis.q[0] = value_type(1.0);
	mAxis.q[1] = value_type(0.0);
	mAxis.q[2] = value_type(0.0);
	mAngle = value_type(0.0);
      };
    };

    /**
     * Constructor from a rotation matrix.
     * \test PASSED
     */
    explicit axis_angle(const rot_mat_3D<value_type>& R) {
      using std::sin;
      using std::acos;
      value_type tmp(value_type(0.5)*(trace(R) - value_type(1.0)));
      if(tmp > value_type(0.0000001)) {
	mAngle = acos(tmp);
	value_type cosec_a = value_type(0.5) / sin(mAngle);
	mAxis.q[0] = (R.q[5] - R.q[7]) * cosec_a;
	mAxis.q[1] = (R.q[6] - R.q[2]) * cosec_a;
	mAxis.q[2] = (R.q[1] - R.q[3]) * cosec_a;
      } else {
	mAxis.q[0] = value_type(1.0);
	mAxis.q[1] = value_type(0.0);
	mAxis.q[2] = value_type(0.0);
	mAngle = value_type(0.0);
      };
    };

    /**
     * Constructor from euler angles.
     * \test PASSED
     */
    explicit axis_angle(const euler_angles_TB<value_type>& E) {
      using std::sin;
      using std::cos;
      using std::sqrt;
      using std::acos;
      value_type cpsi = cos(value_type(0.5)*E.q[0]);
      value_type spsi = sin(value_type(0.5)*E.q[0]);
      value_type ctheta = cos(value_type(0.5)*E.q[1]);
      value_type stheta = sin(value_type(0.5)*E.q[1]);
      value_type cphi = cos(value_type(0.5)*E.q[2]);
      value_type sphi = sin(value_type(0.5)*E.q[2]);

      value_type q[4];
      q[0] = cphi*ctheta*cpsi + sphi*stheta*spsi;
      q[1] = sphi*ctheta*cpsi - cphi*stheta*spsi;
      q[2] = cphi*stheta*cpsi + sphi*ctheta*spsi;
      q[3] = cphi*ctheta*spsi - sphi*stheta*cpsi;

      value_type tmp(sqrt(value_type(1.0)-q[0]*q[0]));
      if(tmp > value_type(0.0000001)) {
	mAxis.q[0] = q[1] / tmp;
	mAxis.q[1] = q[2] / tmp;
	mAxis.q[2] = q[3] / tmp;
	mAngle = value_type(2.0) * acos(q[0]);
      } else {
	mAxis.q[0] = value_type(1.0);
	mAxis.q[1] = value_type(0.0);
	mAxis.q[2] = value_type(0.0);
	mAngle = value_type(0.0);
      };
    };


/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides the angle, read-write.
     * \test PASSED
     */
    reference angle() {
      return mAngle;
    };

    /**
     * Provides the axis, read-write.
     * \test PASSED
     */
    vect<value_type,3>& axis() {
      return mAxis;
    };

    /**
     * Provides the angle, read-only.
     * \test PASSED
     */
    const_reference angle() const {
      return mAngle;
    };

    /**
     * Provides the axis, read-only.
     * \test PASSED
     */
    const vect<value_type,3>& axis() const {
      return mAxis;
    };

    /**
     * Provides a quaternion representation of this rotation.
     * \test PASSED
     */
    quaternion<value_type> getQuaternion() const {
      using std::cos;
      using std::sin;
      value_type t = norm(mAxis);
      if(t == value_type(0.0))
        return quaternion<value_type>(value_type(1.0),value_type(0.0),value_type(0.0),value_type(0.0));
      t = sin(value_type(0.5)*mAngle);
      return quaternion<value_type>(cos(value_type(0.5)*mAngle),mAxis.q[0] * t,mAxis.q[1] * t,mAxis.q[2] * t);
    };

    /**
     * Provides a euler angles representation of this rotation.
     * \test PASSED
     */
    euler_angles_TB<value_type> getEulerAnglesTB() const {
      using std::sin;
      using std::cos;
      using std::asin;
      using std::atan2;
      value_type quat[4];
      value_type t = norm(mAxis);
      if(t == value_type(0.0)) {
        quat[0] = value_type(1.0);
        quat[1] = value_type(0.0);
        quat[2] = value_type(0.0);
        quat[3] = value_type(0.0);
      } else {
        quat[0] = sin(mAngle / value_type(2.0));
        quat[1] = mAxis.q[0] * quat[0];
        quat[2] = mAxis.q[1] * quat[0];
        quat[3] = mAxis.q[2] * quat[0];
        quat[0] = cos(mAngle / value_type(2.0));
      };
      euler_angles_TB<value_type> result;
      result.q[1] = value_type(2.0)*(quat[0]*quat[2] - quat[1]*quat[3]);
      if((result.q[1] != value_type(1.0)) && (result.q[1] != value_type(-1.0))) {
        result.q[1] = asin(result.q[1]);
        value_type cp = value_type(1.0)/cos(result.q[1]);
        result.q[2] = atan2(value_type(2.0)*cp*(quat[2]*quat[3] + quat[0]*quat[1]),cp*(value_type(1.0) - value_type(2.0)*(quat[1]*quat[1] + quat[2]*quat[2])));
        result.q[0] = atan2(value_type(2.0)*cp*(quat[1]*quat[2] + quat[0]*quat[3]),cp*(value_type(1.0) - value_type(2.0)*(quat[2]*quat[2] + quat[3]*quat[3])));
      } else {
        result.q[0] = value_type(0.0);
        result.q[2] = atan2(result.q[1]*value_type(2.0)*(quat[1]*quat[2] - quat[0]*quat[3]),result.q[1]*value_type(2.0)*(quat[1]*quat[3] + quat[0]*quat[2]));
        result.q[1] *= value_type(1.57079632679489662);
      };
      return result;
    };

    /**
     * Provides a rotation matrix representation of this rotation.
     * \test PASSED
     */
    rot_mat_3D<value_type> getRotMat() const {
      using std::cos;
      using std::sin;
      value_type ca(cos(mAngle));
      value_type one_minus_ca(value_type(1.0) - ca);
      value_type t11(ca + one_minus_ca*mAxis.q[0]*mAxis.q[0]);
      value_type t22(ca + one_minus_ca*mAxis.q[1]*mAxis.q[1]);
      value_type t33(ca + one_minus_ca*mAxis.q[2]*mAxis.q[2]);
      value_type t12(one_minus_ca*mAxis.q[0]*mAxis.q[1]);
      value_type t13(one_minus_ca*mAxis.q[0]*mAxis.q[2]);
      value_type t23(one_minus_ca*mAxis.q[1]*mAxis.q[2]);
      value_type sin_a(sin(mAngle));
      value_type t01(sin_a*mAxis.q[0]);
      value_type t02(sin_a*mAxis.q[1]);
      value_type t03(sin_a*mAxis.q[2]);

      return rot_mat_3D<value_type>(t11,t12-t03,t13+t02,
				    t12+t03,t22,t23-t01,
				    t13-t02,t23+t01,t33);
    };

    /**
     * Provides a 3x3 matrix representation of this rotation.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      using std::cos;
      using std::sin;
      value_type ca(cos(mAngle));
      value_type one_minus_ca(value_type(1.0) - ca);
      value_type t11(ca + one_minus_ca*mAxis.q[0]*mAxis.q[0]);
      value_type t22(ca + one_minus_ca*mAxis.q[1]*mAxis.q[1]);
      value_type t33(ca + one_minus_ca*mAxis.q[2]*mAxis.q[2]);
      value_type t12(one_minus_ca*mAxis.q[0]*mAxis.q[1]);
      value_type t13(one_minus_ca*mAxis.q[0]*mAxis.q[2]);
      value_type t23(one_minus_ca*mAxis.q[1]*mAxis.q[2]);
      value_type sin_a(sin(mAngle));
      value_type t01(sin_a*mAxis.q[0]);
      value_type t02(sin_a*mAxis.q[1]);
      value_type t03(sin_a*mAxis.q[2]);

      return mat<value_type,mat_structure::square>(t11,t12-t03,t13+t02,
						   t12+t03,t22,t23-t01,
						   t13-t02,t23+t01,t33);
    };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Assignment from an axis / angle representation.
     * \test PASSED
     */
    self& operator =(const self& A) {
      mAngle = A.mAngle;
      mAxis = A.mAxis;
      return *this;
    };

    /**
     * Standard assignment operator.
     * \test PASSED
     */
    self& operator =(const euler_angles_TB<value_type>& E) {
      return (*this = self(E));
    };

    /**
     * Assignment from a quaternion.
     * \test PASSED
     */
    self& operator =(const quaternion<value_type>& Q) {
      return (*this = self(Q));
    };

    /**
     * Assignment from a rotation matrix.
     * \test PASSED
     */
    self& operator =(const rot_mat_3D<value_type>& R) {
      return (*this = self(R));
    };

    /**
     * Multiply-and-store from a axis / angle.
     * \test PASSED
     */
    self& operator *=(const self& A) {
      return (*this = this->getQuaternion() * A.getQuaternion());
    };

    /**
     * Multiply-and-store from a euler angles.
     * \test PASSED
     */
    self& operator *=(const euler_angles_TB<value_type>& E) {
      return (*this = (this->getRotMat() * E.getRotMat()));
    };

    /**
     * Multiply-and-store from a rotation matrix.
     * \test PASSED
     */
    self& operator *=(const rot_mat_3D<value_type>& R) {
      return (*this = (this->getRotMat() * R));
    };

    /**
     * Multiply-and-store from a quaternion.
     * \test PASSED
     */
    self& operator *=(const quaternion<value_type>& Q) {
      return (*this = (this->getQuaternion() * Q));
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Multiplication with an axis / angle representation.
     * \test PASSED
     */
    friend
    quaternion<value_type> operator *(const self& A1, const self& A2) {
      return A1.getQuaternion() * A2.getQuaternion();
    };

    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const self& A, const Matrix& M) {
      return A.getRotMat() * M;
    };
    
    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
    Matrix >::type operator *(const Matrix& M, const self& A) {
      return M * A.getRotMat();
    };

    /**
     * Multiplication with a column vector.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& A, const vect<value_type,3>& V) {
      return A.getRotMat() * V;
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality comparison operator with a axis / angle representation.
     * \test PASSED
     */
    friend
    bool operator ==(const self& A1, const self& A2) {
      return ((A1.mAngle == A2.mAngle) && (A1.mAxis == A2.mAxis));
    };

    /**
     * Inequality comparison operator with a axis / angle representation.
     * \test PASSED
     */
    friend
    bool operator !=(const self& A1, const self& A2) {
      return ((A1.mAngle != A2.mAngle) || (A1.mAxis != A2.mAxis));
    };

/*******************************************************************************
                         Standard Matrix Methods
*******************************************************************************/

    /**
     * Produces a transpose axis/angle which is the inverse rotation.
     * \test PASSED
     */
    friend 
    self transpose(const self& A) {
      return self(-A.mAngle,A.mAxis);
    };

    /**
     * Produces a cofactor matrix which is the same as the rotation itself.
     * \test PASSED
     */
    friend
    self cofactor(const self& A) {
      return A;
    };

    /**
     * Invert the rotation.
     * \test PASSED
     */
    friend 
    self invert(const self& A) {
      return self(-A.mAngle,A.mAxis);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& A) {
      using std::cos;
      return value_type(2.0)*cos(A.mAngle) + value_type(1.0);
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
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      using std::cos;
      value_type ca(cos(mAngle));
      value_type one_minus_ca(value_type(1.0) - ca);
      value_type t11(ca + one_minus_ca*mAxis.q[0]*mAxis.q[0]);
      value_type t22(ca + one_minus_ca*mAxis.q[1]*mAxis.q[1]);
      value_type t33(ca + one_minus_ca*mAxis.q[2]*mAxis.q[2]);
      value_type t12(one_minus_ca*mAxis.q[0]*mAxis.q[1]);
      value_type t13(one_minus_ca*mAxis.q[0]*mAxis.q[2]);
      value_type t23(one_minus_ca*mAxis.q[1]*mAxis.q[2]);
      return mat<value_type,mat_structure::symmetric>(t11,t12,t13,
				      	                  t22,t23,  
					                      t33);
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      using std::sin;
      value_type sin_a(sin(mAngle));
      value_type t01(sin_a*mAxis.q[0]);
      value_type t02(sin_a*mAxis.q[1]);
      value_type t03(sin_a*mAxis.q[2]);
      return mat<value_type,mat_structure::skew_symmetric>(-t03, t02,
		                                                -t01);
    };
    
        
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(mAngle)
        & RK_SERIAL_SAVE_WITH_NAME(mAxis);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(mAngle)
        & RK_SERIAL_LOAD_WITH_NAME(mAxis);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x0000001C,1,"axis_angle",serialization::serializable)
    
};


/**
 * Prints a axis / angle to a standard output stream (<<) as "(Angle = value; Axis = (a1; a2; a3))".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const axis_angle<T>& A) {
  return out_stream << "(Angle = " << A.angle() << "; Axis = " << A.axis() << ")";
};







/**
 * This class is a transformation matrix 4 by 4, i.e. to rotate and translate a 3D vector.
 * \test All tests for this class have been passed!
 */
template <class T>
class trans_mat_3D : public serialization::serializable {
  public:
    typedef trans_mat_3D<T> self;
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
  
    BOOST_STATIC_CONSTANT(std::size_t, static_row_count = 4);
    BOOST_STATIC_CONSTANT(std::size_t, static_col_count = 4);
    BOOST_STATIC_CONSTANT(mat_alignment::tag, alignment = mat_alignment::column_major);
    BOOST_STATIC_CONSTANT(mat_structure::tag, structure = mat_structure::orthogonal);
    
  private:
    value_type q[16];

    trans_mat_3D(const_reference a11, const_reference a12, const_reference a13, const_reference a14, 
		 const_reference a21, const_reference a22, const_reference a23, const_reference a24, 
		 const_reference a31, const_reference a32, const_reference a33, const_reference a34) {
      q[0] = a11; q[1] = a21; q[2] = a31; q[3] = value_type(0.0);
      q[4] = a12; q[5] = a22; q[6] = a32; q[7] = value_type(0.0);
      q[8] = a13; q[9] = a23; q[10] = a33; q[11] = value_type(0.0);
      q[12] = a14; q[13] = a24; q[14] = a34; q[15] = value_type(1.0);
    };

//     trans_mat_3D operator =(const mat<T>& M) {
//       q[0] = M[0];
//       q[1] = M[1];
//       q[2] = M[2];
//       q[3] = 0.0;
//       q[4] = M[4];
//       q[5] = M[5];
//       q[6] = M[6];
//       q[7] = 0.0;
//       q[8] = M[8];
//       q[9] = M[9];
//       q[10] = M[10];
//       q[11] = 0.0;
//       q[12] = M[12];
//       q[13] = M[13];
//       q[14] = M[14];
//       q[15] = 1.0;
//       return *this;
//     };

  public:

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default Constructor.
     * \test PASSED
     */
    trans_mat_3D() {
      q[0] = value_type(1.0); q[1] = value_type(0.0); q[2] = value_type(0.0); q[3] = value_type(0.0);
      q[4] = value_type(0.0); q[5] = value_type(1.0); q[6] = value_type(0.0); q[7] = value_type(0.0);
      q[8] = value_type(0.0); q[9] = value_type(0.0); q[10] = value_type(1.0); q[11] = value_type(0.0);
      q[12] = value_type(0.0); q[13] = value_type(0.0); q[14] = value_type(0.0); q[15] = value_type(1.0);
    };

    /**
     * Constructor from a 4x4 array (16 values).
     * \test PASSED
     */
    explicit trans_mat_3D(const_pointer M) {
      vect<value_type,3> v1 = unit(vect<value_type,3>(M));
      q[0] = v1[0];
      q[1] = v1[1];
      q[2] = v1[2];
      q[3] = value_type(0.0);
      vect<value_type,3> v2 = unit(v2 - (v2 * v1) * v1);
      q[4] = v2[0];
      q[5] = v2[1];
      q[6] = v2[2];
      q[7] = value_type(0.0);
      v2 = v1 % v2;
      q[8] = v2[0];
      q[9] = v2[1];
      q[10] = v2[2];
      q[11] = value_type(0.0);
      q[12]= M[12];
      q[13]= M[13];
      q[14]= M[14];
      q[15]= value_type(1.0);
      return;
    };

    /**
     * Constructor from a quaternion representation and an optional translation vector V.
     * \test PASSED
     */
    explicit trans_mat_3D(const quaternion<value_type>& Q, 
			  const vect<value_type,3>& V = vect<value_type,3>(value_type(0.0),value_type(0.0),value_type(0.0))) {
      rot_mat_3D<value_type> R(Q.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = V[0];
      q[13] = V[1];
      q[14] = V[2];
      q[15] = value_type(1.0);
    };

    /**
     * Constructor from a euler angles TB representation and an optional translation vector V.
     * \test PASSED
     */
    explicit trans_mat_3D(const euler_angles_TB<value_type>& E, 
			  const vect<value_type,3>& V = vect<value_type,3>(value_type(0.0),value_type(0.0),value_type(0.0))) {
      rot_mat_3D<value_type> R(E.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = V[0];
      q[13] = V[1];
      q[14] = V[2];
      q[15] = value_type(1.0);
    };

    /**
     * Constructor from an axis / angle representation and an optional translation vector V.
     * \test PASSED
     */
    explicit trans_mat_3D(const axis_angle<value_type>& A, 
			  const vect<value_type,3>& V = vect<value_type,3>(value_type(0.0),value_type(0.0),value_type(0.0))) {
      rot_mat_3D<value_type> R(A.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = V[0];
      q[13] = V[1];
      q[14] = V[2];
      q[15] = value_type(1.0);
    };

    /**
     * Constructor from a rotation matrix and an optional translation vector V.
     * \test PASSED
     */
    explicit trans_mat_3D(const rot_mat_3D<value_type>& R, 
			  const vect<value_type,3>& V = vect<value_type,3>(value_type(0.0),value_type(0.0),value_type(0.0))) {
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = V[0];
      q[13] = V[1];
      q[14] = V[2];
      q[15] = value_type(1.0);
    };

    // Copy-constructor. Default is good. \test PASSED

    /**
     * Destructor.
     * \test PASSED
     */
    ~trans_mat_3D() { };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Provides a copy of the transformation matrix as an ordinary 4x4 matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::square> getMat() const {
      return mat<value_type,mat_structure::square>(q[0],q[4],q[8],q[12],q[1],q[5],q[9],q[13],q[2],q[6],q[10],q[14],value_type(0.0),value_type(0.0),value_type(0.0),value_type(1.0));
    };

    /**
     * Provides the rotation part of the transformation as a rotation matrix.
     * \test PASSED
     */
    rot_mat_3D<value_type> getRotMat() const {
      return rot_mat_3D<value_type>(q[0],q[4],q[8],
			            q[1],q[5],q[9],
			            q[2],q[6],q[10]);
    };

    /**
     * Sets the rotation part of the transformation from a rotation matrix.
     * \test PASSED
     */
    void setRotMat(const rot_mat_3D<value_type>& R) {
      q[0]  = R.q[0];
      q[1]  = R.q[1];
      q[2]  = R.q[2];
      q[4]  = R.q[3];
      q[5]  = R.q[4];
      q[6]  = R.q[5];
      q[8]  = R.q[6];
      q[9]  = R.q[7];
      q[10] = R.q[8];
    };

    /**
     * Returns the quaternion of the rotation matrix.
     * \test PASSED
     */
    quaternion<value_type> getQuaternion() const {
      return quaternion<value_type>(getRotMat());
    };

    /**
     * Sets the quaternion of the rotation matrix.
     * \test PASSED
     */
    void setQuaternion(const quaternion<value_type>& Q) {
      setRotMat(Q.getRotMat());
    };

    /**
     * Returns the euler angles TB of the rotation matrix.
     * \test PASSED
     */
    euler_angles_TB<value_type> getEulerAnglesTB() const {
      return euler_angles_TB<value_type>(getRotMat());
    };

    /**
     * Sets the euler angles TB of the rotation matrix.
     * \test PASSED
     */
    void setEulerAnglesTB(const euler_angles_TB<value_type>& E) {
      setRotMat(E.getRotMat());
    };

    /**
     * Returns the axis / angle of the rotation matrix.
     * \test PASSED
     */
    axis_angle<value_type> getAxisAngle() const {
      return axis_angle<value_type>(getRotMat());
    };

    /**
     * Sets the axis / angle of the rotation matrix.
     * \test PASSED
     */
    void setAxisAngle(const axis_angle<value_type>& A) {
      setRotMat(A.getRotMat());
    };

    /**
     * Provides the translation part of the transformation matrix as a vector.
     * \test PASSED
     */
    vect<value_type,3> getTranslation() const {
      return vect<value_type,3>(q[12],q[13],q[14]);
    };

    /**
     * Sets the translation part of the transformation matrix to a vector.
     * \test PASSED
     */
    void setTranslation(const vect<value_type,3>& Translation) {
      q[12] = Translation[0];
      q[13] = Translation[1];
      q[14] = Translation[2];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      if(i >= 16)
	throw std::range_error("Matrix index out of range.");
      return q[i];
    };

    /**
     * Array double-indexing operator, ith row and jth column, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i,size_type j) const {
      if((i >= 4) || (j >= 4))
	throw std::range_error("Matrix index out of range.");
      return q[j*4+i];
    };
    
    size_type get_row_count() const { return 4; };
    size_type get_col_count() const { return 4; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    self& operator =(const self& M) {
      q[0] = M.q[0];
      q[1] = M.q[1];
      q[2] = M.q[2];
      q[3] = value_type(0.0);
      q[4] = M.q[4];
      q[5] = M.q[5];
      q[6] = M.q[6];
      q[7] = value_type(0.0);
      q[8] = M.q[8];
      q[9] = M.q[9];
      q[10]= M.q[10];
      q[11]= value_type(0.0);
      q[12]= M.q[12];
      q[13]= M.q[13];
      q[14]= M.q[14];
      q[15]= value_type(1.0);
      return *this;
    }; 

    /**
     * Assignment operator with rotation matrix.
     * \test PASSED
     */
    self& operator =(const rot_mat_3D<value_type>& R)  {
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = value_type(0.0);
      q[13] = value_type(0.0);
      q[14] = value_type(0.0);
      q[15] = value_type(1.0);
      return *this;
    };

    /**
     * Assignment operator with a quaternion representation.
     * \test PASSED
     */
    self& operator =(const quaternion<value_type>& Q) {
      rot_mat_3D<value_type> R(Q.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = value_type(0.0);
      q[13] = value_type(0.0);
      q[14] = value_type(0.0);
      q[15] = value_type(1.0);
      return *this;
    };

    /**
     * Assignment operator with euler angles TB representation.
     * \test PASSED
     */
    self& operator =(const euler_angles_TB<value_type>& E) {
      rot_mat_3D<value_type> R(E.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = value_type(0.0);
      q[13] = value_type(0.0);
      q[14] = value_type(0.0);
      q[15] = value_type(1.0);
      return *this;
    };

    /**
     * Assignment operator with an axis / angle representation.
     * \test PASSED
     */
    self& operator =(const axis_angle<value_type>& A) {
      rot_mat_3D<value_type> R(A.getRotMat());
      q[0] = R.q[0];
      q[1] = R.q[1];
      q[2] = R.q[2];
      q[3] = value_type(0.0);
      q[4] = R.q[3];
      q[5] = R.q[4];
      q[6] = R.q[5];
      q[7] = value_type(0.0);
      q[8] = R.q[6];
      q[9] = R.q[7];
      q[10] = R.q[8];
      q[11] = value_type(0.0);
      q[12] = value_type(0.0);
      q[13] = value_type(0.0);
      q[14] = value_type(0.0);
      q[15] = value_type(1.0);
      return *this;
    };

    /**
     * Multiply-and-store with a transformation matrix.
     * \test PASSED
     */
    self& operator *=(const self& M) {
      (*this) = self(q[0]*M.q[0] + q[4]*M.q[1] + q[8] *M.q[2], q[0]*M.q[4] + q[4]*M.q[5] + q[8] *M.q[6], q[0]*M.q[8] + q[4]*M.q[9] + q[8] *M.q[10], q[0]*M.q[12] + q[4]*M.q[13] + q[8]*M.q[14]  + q[12],
                     q[1]*M.q[0] + q[5]*M.q[1] + q[9] *M.q[2], q[1]*M.q[4] + q[5]*M.q[5] + q[9] *M.q[6], q[1]*M.q[8] + q[5]*M.q[9] + q[9] *M.q[10], q[1]*M.q[12] + q[5]*M.q[13] + q[9]*M.q[14]  + q[13],
                     q[2]*M.q[0] + q[6]*M.q[1] + q[10]*M.q[2], q[2]*M.q[4] + q[6]*M.q[5] + q[10]*M.q[6], q[2]*M.q[8] + q[6]*M.q[9] + q[10]*M.q[10], q[2]*M.q[12] + q[6]*M.q[13] + q[10]*M.q[14] + q[14]);
      return *this;
    };

    /**
     * Multiply-and-store with a rotation matrix.
     * \test PASSED
     */
    self& operator *=(const rot_mat_3D<value_type>& R) {
      (*this) = self(q[0]*R.q[0] + q[4]*R.q[1] + q[8] *R.q[2], q[0]*R.q[3] + q[4]*R.q[4] + q[8] *R.q[5], q[0]*R.q[6] + q[4]*R.q[7] + q[8] *R.q[8], q[12],
                     q[1]*R.q[0] + q[5]*R.q[1] + q[9] *R.q[2], q[1]*R.q[3] + q[5]*R.q[4] + q[9] *R.q[5], q[1]*R.q[6] + q[5]*R.q[7] + q[9] *R.q[8], q[13],
                     q[2]*R.q[0] + q[6]*R.q[1] + q[10]*R.q[2], q[2]*R.q[3] + q[6]*R.q[4] + q[10]*R.q[5], q[2]*R.q[6] + q[6]*R.q[7] + q[10]*R.q[8], q[14]);
      return *this;
    };

    /**
     * Multiply-and-store with a quaternion representation.
     * \test PASSED
     */
    self& operator *=(const quaternion<value_type>& Q) {
      return (*this *= Q.getRotMat());
    };

    /**
     * Multiply-and-store with a euler angles TB representation.
     * \test PASSED
     */
    self& operator *=(const euler_angles_TB<value_type>& E) {
      return (*this *= E.getRotMat());
    };

    /**
     * Multiply-and-store with an axis / angle representation.
     * \test PASSED
     */
    self& operator *=(const axis_angle<value_type>& A) {
      return (*this *= A.getRotMat());
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Multiplication with a transformation matrix.
     * \test PASSED
     */
    friend
    self operator *(const self& M1, const self& M2) {
      return self(M1.q[0]*M2.q[0] + M1.q[4]*M2.q[1] + M1.q[8] *M2.q[2], M1.q[0]*M2.q[4] + M1.q[4]*M2.q[5] + M1.q[8] *M2.q[6], M1.q[0]*M2.q[8] + M1.q[4]*M2.q[9] + M1.q[8] *M2.q[10], M1.q[0]*M2.q[12] + M1.q[4]*M2.q[13] + M1.q[8]*M2.q[14]  + M1.q[12],
                  M1.q[1]*M2.q[0] + M1.q[5]*M2.q[1] + M1.q[9] *M2.q[2], M1.q[1]*M2.q[4] + M1.q[5]*M2.q[5] + M1.q[9] *M2.q[6], M1.q[1]*M2.q[8] + M1.q[5]*M2.q[9] + M1.q[9] *M2.q[10], M1.q[1]*M2.q[12] + M1.q[5]*M2.q[13] + M1.q[9]*M2.q[14]  + M1.q[13],
                  M1.q[2]*M2.q[0] + M1.q[6]*M2.q[1] + M1.q[10]*M2.q[2], M1.q[2]*M2.q[4] + M1.q[6]*M2.q[5] + M1.q[10]*M2.q[6], M1.q[2]*M2.q[8] + M1.q[6]*M2.q[9] + M1.q[10]*M2.q[10], M1.q[2]*M2.q[12] + M1.q[6]*M2.q[13] + M1.q[10]*M2.q[14] + M1.q[14]);
    };

    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value &&
                                 !boost::is_same<Matrix,rot_mat_3D<value_type> >::value,
    Matrix >::type operator *(const self& M1, const Matrix& M2) {
      return Matrix(M1.getMat() * M2);
    };
    
    /**
     * Multiply by a matrix.
     * \test PASSED
     */
    template <typename Matrix>
    friend
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value &&
                                 !boost::is_same<Matrix,self>::value &&
                                 !boost::is_same<Matrix,rot_mat_3D<value_type> >::value,
    Matrix >::type operator *(const Matrix& M1, const self& M2) {
      return Matrix(M1 * M2.getMat());
    };

    /**
     * Multiplication with a rotation matrix.
     * \test PASSED
     */
    friend
    self operator *(const self& M, const rot_mat_3D<T>& R) {
      return self(M.q[0]*R[0] + M.q[4]*R[1] + M.q[8] *R[2], M.q[0]*R[3] + M.q[4]*R[4] + M.q[8] *R[5], M.q[0]*R[6] + M.q[4]*R[7] + M.q[8] *R[8], M.q[12],
                  M.q[1]*R[0] + M.q[5]*R[1] + M.q[9] *R[2], M.q[1]*R[3] + M.q[5]*R[4] + M.q[9] *R[5], M.q[1]*R[6] + M.q[5]*R[7] + M.q[9] *R[8], M.q[13],
                  M.q[2]*R[0] + M.q[6]*R[1] + M.q[10]*R[2], M.q[2]*R[3] + M.q[6]*R[4] + M.q[10]*R[5], M.q[2]*R[6] + M.q[6]*R[7] + M.q[10]*R[8], M.q[14]);
    };

    /**
     * Multiplication with a 3D column vector.
     * \test PASSED
     */
    friend
    vect<value_type,3> operator *(const self& M, const vect<value_type,3>& V) {
      return vect<value_type,3>( M.q[0]*V[0] + M.q[4]*V[1] + M.q[8] *V[2] + M.q[12],
                                 M.q[1]*V[0] + M.q[5]*V[1] + M.q[9] *V[2] + M.q[13],
                                 M.q[2]*V[0] + M.q[6]*V[1] + M.q[10]*V[2] + M.q[14]);
    };

    /**
     * Multiplication with a 4D column vector.
     * \test PASSED
     */
    friend
    vect<value_type,4> operator *(const self& M, const vect<value_type,4>& V) {
      return vect<value_type,4>( M.q[0]*V[0] + M.q[4]*V[1] + M.q[8] *V[2] + M.q[12]*V[3],
				 M.q[1]*V[0] + M.q[5]*V[1] + M.q[9] *V[2] + M.q[13]*V[3],
				 M.q[2]*V[0] + M.q[6]*V[1] + M.q[10]*V[2] + M.q[14]*V[3],
				 V[3]);
    };

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality operator with a transformation matrix.
     * \test PASSED
     */
    friend
    bool operator ==(const self& M1, const self& M2) {
      return ((M1.q[0] == M2.q[0]) && (M1.q[1] == M2.q[1]) && (M1.q[2] == M2.q[2]) &&
              (M1.q[4] == M2.q[4]) && (M1.q[5] == M2.q[5]) && (M1.q[6] == M2.q[6]) &&
              (M1.q[12] == M2.q[12]) && (M1.q[13] == M2.q[13]) && (M1.q[14] == M2.q[14]));
    };

    /**
     * Inequality operator with a transformation matrix.
     * \test PASSED
     */
    friend
    bool operator !=(const self& M1, const self& M2) {
      return ((M1.q[0] != M2.q[0]) || (M1.q[1] != M2.q[1]) || (M1.q[2] != M2.q[2]) ||
              (M1.q[4] != M2.q[4]) || (M1.q[5] != M2.q[5]) || (M1.q[6] != M2.q[6]) ||
              (M1.q[12] != M2.q[12]) || (M1.q[13] != M2.q[13]) || (M1.q[14] != M2.q[14]));
    };

/*******************************************************************************
                         Special Methods
*******************************************************************************/

    /**
     * Rotate-only a 3D column vector.
     * \test PASSED
     */
    vect<value_type,3> rotate(const vect<value_type,3>& V) const {
      return vect<value_type,3>(q[0]*V[0] + q[4]*V[1] + q[8]*V[2],
				q[1]*V[0] + q[5]*V[1] + q[9]*V[2],
				q[2]*V[0] + q[6]*V[1] + q[10]*V[2]);
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
      return mat<value_type,mat_structure::square>(M.q[0],M.q[1],M.q[2],0.0,M.q[4],M.q[5],M.q[6],0.0,M.q[8],M.q[9],M.q[10],0.0,M.q[11],M.q[12],M.q[13],1.0);
    };
    
    /**
     * Creates the transpose matrix.
     * \note the matrix is no longer a transformation matrix.
     * \test PASSED
     */
    friend
    mat<value_type,mat_structure::square> transpose_move(const self& M) {
      return mat<value_type,mat_structure::square>(M.q[0],M.q[1],M.q[2],0.0,M.q[4],M.q[5],M.q[6],0.0,M.q[8],M.q[9],M.q[10],0.0,M.q[11],M.q[12],M.q[13],1.0);
    };

    /**
     * Gets the trace of the matrix.
     * \test PASSED
     */
    friend
    value_type trace(const self& M) {
      return M.q[0] + M.q[5] + M.q[10] + value_type(1.0);
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
      return self(M.q[0],M.q[1],M.q[2],-M.q[0]*M.q[12]-M.q[1]*M.q[13]-M.q[2]*M.q[14],
		  M.q[4],M.q[5],M.q[6],-M.q[4]*M.q[12]-M.q[5]*M.q[13]-M.q[6]*M.q[14],
		  M.q[8],M.q[9],M.q[10],-M.q[8]*M.q[12]-M.q[9]*M.q[13]-M.q[10]*M.q[14]);
    };

    /**
     * Gets the symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::symmetric> getSymPart() const {
      return mat<value_type,mat_structure::symmetric>(getMat());
    };

    /**
     * Gets the skew-symmetric part of the matrix.
     * \test PASSED
     */
    mat<value_type,mat_structure::skew_symmetric> getSkewSymPart() const {
      return mat<value_type,mat_structure::skew_symmetric>(getMat());
    };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_ALIAS("r11",q[0])
        & RK_SERIAL_SAVE_WITH_ALIAS("r21",q[1])
	& RK_SERIAL_SAVE_WITH_ALIAS("r31",q[2])
	& RK_SERIAL_SAVE_WITH_ALIAS("r12",q[4])
	& RK_SERIAL_SAVE_WITH_ALIAS("r22",q[5])
	& RK_SERIAL_SAVE_WITH_ALIAS("r32",q[6])
	& RK_SERIAL_SAVE_WITH_ALIAS("r13",q[8])
	& RK_SERIAL_SAVE_WITH_ALIAS("r23",q[9])
	& RK_SERIAL_SAVE_WITH_ALIAS("r33",q[10])
	& RK_SERIAL_SAVE_WITH_ALIAS("t_x",q[12])
	& RK_SERIAL_SAVE_WITH_ALIAS("t_y",q[13])
	& RK_SERIAL_SAVE_WITH_ALIAS("t_z",q[14]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("r11",q[0])
        & RK_SERIAL_LOAD_WITH_ALIAS("r21",q[1])
	& RK_SERIAL_LOAD_WITH_ALIAS("r31",q[2])
	& RK_SERIAL_LOAD_WITH_ALIAS("r12",q[4])
	& RK_SERIAL_LOAD_WITH_ALIAS("r22",q[5])
	& RK_SERIAL_LOAD_WITH_ALIAS("r32",q[6])
	& RK_SERIAL_LOAD_WITH_ALIAS("r13",q[8])
	& RK_SERIAL_LOAD_WITH_ALIAS("r23",q[9])
	& RK_SERIAL_LOAD_WITH_ALIAS("r33",q[10])
	& RK_SERIAL_LOAD_WITH_ALIAS("t_x",q[12])
	& RK_SERIAL_LOAD_WITH_ALIAS("t_y",q[13])
	& RK_SERIAL_LOAD_WITH_ALIAS("t_z",q[14]);
      q[3] = 0.0; q[7] = 0.0; q[11] = 0.0; q[15] = 1.0;
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0x00000019,1,"trans_mat_3D",serialization::serializable)

};

/**
 * Prints a 3D transformation matrix to a standard output stream (<<)
 * as "(quaternion = (q0; q1; q2; q3); translation = (tx; ty; tz))".
 * \test PASSED
 */
template <class T>
std::ostream& operator <<(std::ostream& out_stream,const trans_mat_3D<T>& M) {
  out_stream << "(quaternion = " << M.getQuaternion() << "; translation = " << M.getTranslation() << ")";
  return out_stream;
};



template <typename T>
struct is_readable_matrix< trans_mat_3D<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_matrix< trans_mat_3D<T> > type;
};















/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const rot_mat_3D<T>& R, const quaternion<T>& Q) {
//   return R * Q.getRotMat();
// };

/**
 * Multiplication by a rotation matrix.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const quaternion<T>& Q, const rot_mat_3D<T>& R) {
//   return Q.getRotMat() * R;
// };

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const rot_mat_3D<T>& R, const euler_angles_TB<T>& E) {
//   return R * E.getRotMat();
// };

/**
 * Multiply by a rotation matrix representation.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const euler_angles_TB<T>& E, const rot_mat_3D<T>& R) {
//   return E.getRotMat() * R;
// };

/**
 * Multiplication by a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
quaternion<T> operator *(const quaternion<T>& Q, const euler_angles_TB<T>& E) {
  return Q * E.getQuaternion();
};

/**
 * Multiply by a quaternion representation.
 * \test PASSED
 */
template <typename T>
quaternion<T> operator *(const euler_angles_TB<T>& E, const quaternion<T>& Q) {
  return E.getQuaternion() * Q;
};

/**
 * Multiplication with an axis / angle representation.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const rot_mat_3D<T>& R, const axis_angle<T>& A) {
//   return R * A.getRotMat();
// };

/**
 * Multiplication with a rotation matrix.
 * \test PASSED
 */
// template <typename T>
// rot_mat_3D<T> operator *(const axis_angle<T>& A, const rot_mat_3D<T>& R) {
//   return A.getRotMat() * R;
// };

/**
 * Multiplication by an axis / angle representation.
 * \test PASSED
 */
template <typename T>
quaternion<T> operator *(const quaternion<T>& Q, const axis_angle<T>& A) {
  return Q * A.getQuaternion();
};

/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
template <typename T>
quaternion<T> operator *(const axis_angle<T>& A, const quaternion<T>& Q) {
  return A.getQuaternion() * Q;
};
    
/**
 * Multiply by an axis / angle representation.
 * \test PASSED
 */
template <typename T>
rot_mat_3D<T> operator *(const euler_angles_TB<T>& E, const axis_angle<T>& A) {
  return E.getRotMat() * A.getRotMat();
};

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
rot_mat_3D<T> operator *(const axis_angle<T>& A, const euler_angles_TB<T>& E) {
  return A.getRotMat() * E.getRotMat();
};

/**
 * Multiplication with a transformation matrix.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const rot_mat_3D<T>& R, const trans_mat_3D<T>& M) {
  return trans_mat_3D<T>(R) * M;
};

/**
 * Multiplication by a transformation matrix.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const quaternion<T>& Q, const trans_mat_3D<T>& M) {
  return trans_mat_3D<T>(Q.getRotMat()) * M;
};

/**
 * Multiplication with a quaternion representation.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const trans_mat_3D<T>& M, const quaternion<T>& Q) {
  return M * Q.getRotMat();
};

/**
 * Multiply by a transformation matrix.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const euler_angles_TB<T>& E, const trans_mat_3D<T>& M) {
  return trans_mat_3D<T>(E) * M;
};

/**
 * Multiplication with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const trans_mat_3D<T>& M, const euler_angles_TB<T>& E) {
  return M * E.getRotMat();
};

/**
 * Multiplication with a transformation matrix.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const axis_angle<T>& A, const trans_mat_3D<T>& M) {
  return trans_mat_3D<T>(A) * M;
};

/**
 * Multiplication with an axis / angle representation.
 * \test PASSED
 */
template <typename T>
trans_mat_3D<T> operator *(const trans_mat_3D<T>& M, const axis_angle<T>& A) {
  return M * A.getRotMat();
};
    








/**
 * Equality operator for a quaternion representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const rot_mat_3D<T>& R, const quaternion<T>& Q) {
  return Q.getRotMat() == R;
};

/**
 * Inequality operator for a quaternion representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const rot_mat_3D<T>& R, const quaternion<T>& Q) {
  return Q.getRotMat() != R;
};

/**
 * Equality operator with a rotation matrix.
 * \test PASSED
 */
template <typename T>
bool operator ==(const quaternion<T>& Q, const rot_mat_3D<T>& R) {
  return Q.getRotMat() == R;
};

/**
 * Inequality operator with a rotation matrix.
 * \test PASSED
 */
template <typename T>
bool operator !=(const quaternion<T>& Q, const rot_mat_3D<T>& R) {
  return Q.getRotMat() != R;
};

/**
 * Equality operator with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const quaternion<T>& Q, const euler_angles_TB<T>& E) {
  return Q == E.getQuaternion();
};

/**
 * Inequality operator with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const quaternion<T>& Q, const euler_angles_TB<T>& E) {
  return Q != E.getQuaternion();
};

/**
 * Equality operator with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const euler_angles_TB<T>& E, const quaternion<T>& Q) {
  return Q == E.getQuaternion();
};

/**
 * Inequality operator with a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const euler_angles_TB<T>& E, const quaternion<T>& Q) {
  return Q != E.getQuaternion();
};

/**
 * Equality operator for a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const rot_mat_3D<T>& R, const euler_angles_TB<T>& E) {
  return R == E.getRotMat();
};

/**
 * Inequality operator for a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const rot_mat_3D<T>& R, const euler_angles_TB<T>& E) {
  return R != E.getRotMat();
};

/**
 * Equality operator for a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const euler_angles_TB<T>& E, const rot_mat_3D<T>& R) {
  return R == E.getRotMat();
};

/**
 * Inequality operator for a euler angles TB representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const euler_angles_TB<T>& E, const rot_mat_3D<T>& R) {
  return R != E.getRotMat();
};

/**
 * Equality operator for an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const rot_mat_3D<T>& R, const axis_angle<T>& A) {
  return R == A.getRotMat();
};

/**
 * Inequality operator for an axis /angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const rot_mat_3D<T>& R, const axis_angle<T>& A) {
  return R != A.getRotMat();
};

/**
 * Equality operator for an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const axis_angle<T>& A, const rot_mat_3D<T>& R) {
  return R == A.getRotMat();
};

/**
 * Inequality operator for an axis /angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const axis_angle<T>& A, const rot_mat_3D<T>& R) {
  return R != A.getRotMat();
};

/**
 * Equality operator with an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const quaternion<T>& Q, const axis_angle<T>& A) {
  return Q == A.getQuaternion();
};

/**
 * Inequality operator with an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const quaternion<T>& Q, const axis_angle<T>& A) {
  return Q != A.getQuaternion();
};

/**
 * Equality operator with an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const axis_angle<T>& A, const quaternion<T>& Q) {
  return Q == A.getQuaternion();
};

/**
 * Inequality operator with an axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const axis_angle<T>& A, const quaternion<T>& Q) {
  return Q != A.getQuaternion();
};

/**
 * Equality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const euler_angles_TB<T>& E, const axis_angle<T>& A) {
  return E == A.getEulerAnglesTB();
};

/**
 * Inequality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const euler_angles_TB<T>& E, const axis_angle<T>& A) {
  return E != A.getEulerAnglesTB();
};

/**
 * Equality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator ==(const axis_angle<T>& A, const euler_angles_TB<T>& E) {
  return E == A.getEulerAnglesTB();
};

/**
 * Inequality comparison operator with a axis / angle representation.
 * \test PASSED
 */
template <typename T>
bool operator !=(const axis_angle<T>& A, const euler_angles_TB<T>& E) {
  return E != A.getEulerAnglesTB();
};







};


#endif












