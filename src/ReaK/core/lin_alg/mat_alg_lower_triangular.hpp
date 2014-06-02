/**
 * \file mat_alg_lower_triangular.hpp
 * 
 * This library contains old code which is specific to handing, storing and operating 
 * with lower-triangular matrices. It stores only the lower-triangular part and operate
 * more efficiently than if it was a dense matrix.
 * 
 * \todo Port all this code to the newer framework for matrices.
 * 
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

#ifndef REAK_MAT_ALG_LOWER_TRIANGULAR_HPP
#define REAK_MAT_ALG_LOWER_TRIANGULAR_HPP

#include "mat_alg_general.hpp"


namespace ReaK {



#if 0



/**
 * This class holds a lower-triangular matrix. This class will hold only the lower-triangular part.
 */
template <class T>
class mat_lo_tri : public mat<T> {
  protected:
    std::vector<T> q; ///< Holds the array of scalar entries.
    unsigned int size; ///< Holds the dimension, both row and column count are equal to size.

  public:

    /**
     * This embedded class allows for direct access to the array q, for faster access than general channels.
     */
    class access {
      public:
        std::vector<T>& q; ///< Public reference to the q array.
        unsigned int& size; ///< Public reference to the size.

        /**
         * Explicit RAII association of this object with a mat_lo_tri object.
         */
        explicit access(mat_lo_tri<T>& M) : q(M.q), size(M.size) { };

        /**
         * Direct array access with flat index.
         * \param i flat index to the array q.
         * \return reference to the entry.
         */
        T& operator [](unsigned int i) { return q[i]; };
        /**
         * Direct array access with flat index. Read-only version.
         * \param i flat index to the array q.
         * \return read-only reference to the entry.
         */
        const T& operator [](unsigned int i) const { return q[i]; };

        /**
         * Direct array access with row-column indices.
         * \param i row index.
         * \param j column index.
         * \return reference to the entry.
         */
        T& operator ()(unsigned int i,unsigned int j) { return q[mat_triangular_size(i) + j]; };
        /**
         * Direct array access with row-column indices. Read-only version.
         * \param i row index.
         * \param j column index.
         * \return read-only reference to the entry.
         */
        const T& operator ()(unsigned int i,unsigned int j) const { return q[mat_triangular_size(i) + j]; };
    };

    /**
     * This embedded class allows for direct read-only access to the array q, for faster access than general channels.
     */
    class const_access {
      public:
        const std::vector<T>& q; ///< Public read-only reference to the q array.
        const unsigned int& size; ///< Public read-only reference to the size.

        /**
         * Explicit RAII association of this object with a mat_lo_tri object.
         */
        explicit const_access(const mat_lo_tri<T>& M) : q(M.q), size(M.size) { };

        /**
         * Direct array access with flat index. Read-only version.
         * \param i flat index to the array q.
         * \return read-only reference to the entry.
         */
        const T& operator [](unsigned int i) const { return q[i]; };

        /**
         * Direct array access with row-column indices. Read-only version.
         * \param i row index.
         * \param j column index.
         * \return read-only reference to the entry.
         */
        const T& operator ()(unsigned int i,unsigned int j) const { return q[mat_triangular_size(i) + j]; };
    };

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
    /**
     * Default constructor: sets all to zero.
     */
    mat_lo_tri() :
             q(0),
             size(0) { };

    /**
     * Constructor for a sized matrix.
     */
    mat_lo_tri(unsigned int Size, T aFill = 0) :
             q(mat_triangular_size(Size),aFill),
             size(Size) { };

    /**
     * Constructor for an identity matrix.
     */
    mat_lo_tri(unsigned int Size, bool aIdentity) :
             q(mat_triangular_size(Size),T(0.0)),
             size(Size) {
      if(aIdentity) {
        unsigned int k=0;
        for(unsigned int i=0; i < size; k += ++i)
          q[k+i] = 1.0;
      };
    };

    /**
     * Standard Copy Constructor with standard semantics.
     */
    mat_lo_tri(const mat_lo_tri<T>& M) :
             q(M.q.begin(),M.q.end()),
             size(M.size) { };

    /**
     * Constructor from a general matrix, copying only the lower triangular part.
     */
    explicit mat_lo_tri(const mat<T>& M) :
             q(mat_triangular_size((M.getRowCount() > M.getColCount() ? M.getRowCount() : M.getColCount())),T(0.0)),
             size((M.getRowCount() > M.getColCount() ? M.getRowCount() : M.getColCount())) {
      unsigned int k=0;
      unsigned int i=0;
      unsigned int min_size = (M.getRowCount() > M.getColCount() ? M.getColCount() : M.getRowCount());
      for(;i<min_size;k += ++i) {
        for(unsigned int j=0;j<i;++j) {
          q[k+j] = M(i,j);
        };
        q[k+i] = M(i,i);
      };
      if(M.getRowCount() > M.getColCount()) {
        for(;i<size;k += ++i) {
             for(unsigned int j=0;j<min_size;++j)
            q[k+j] = M(i,j);
        };
      };
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~mat_lo_tri() { };

    /**
     * Constructs a 2x2 lower-triangular matrix from three elements.
     * \test PASSED
     */
    mat_lo_tri(T a11,T a21,T a22) : q(3), size(2) {
      q[0] = a11;
      q[1] = a21;
      q[2] = a22;
    };

    /**
     * Constructs a 3x3 lower-triangular matrix from six elements.
     * \test PASSED
     */
    mat_lo_tri(T a11,T a21,T a31,T a22,T a32,T a33) : q(6), size(3) {
      q[0] = a11;
      q[1] = a21;
      q[2] = a22;
      q[3] = a31;
      q[4] = a32;
      q[5] = a33;
    };

    /**
     * Constructs a 4x4 lower-triangular matrix from ten elements.
     * \test PASSED
     */
    mat_lo_tri(T a11,T a21,T a31,T a41,T a22,T a32,T a42,T a33,T a43,T a44) : q(10), size(4) {
      q[0] = a11;
      q[1] = a21;
      q[2] = a22;
      q[3] = a31;
      q[4] = a32;
      q[5] = a33;
      q[6] = a41;
      q[7] = a42;
      q[8] = a43;
      q[9] = a44;
    };


/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    const T& set(unsigned int i,unsigned int j,const T& Value) {
      if(i >= j)
        q[mat_triangular_size(i) + j] = Value;
      return Value;
    };

    T operator()(unsigned int i,unsigned int j) const {
      if(i >= j)
        return q[mat_triangular_size(i) + j];
      else
        return T(0.0);
    };

    unsigned int getRowCount() const {
      return size;
    };

    void setRowCount(unsigned int aRowCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(mat_triangular_size(aRowCount),T(0.0));
      size = aRowCount;
    };

    unsigned int getColCount() const {
      return size;
    };

    void setColCount(unsigned int aColCount,bool aPreserveData = false) { RK_UNUSED(aPreserveData);
      q.resize(mat_triangular_size(aColCount),T(0.0));
      size = aColCount;
    };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard Assignment operator with a lower-triangular matrix.
     */
    mat_lo_tri<T>& operator =(const mat_lo_tri<T>& M) {
      size = M.size;
      q.assign(M.q.begin(),M.q.end());
      return *this;
    };

    /**
     * Standard Assignment operator with a general matrix. Copying only the lower-triangular part of M.
     */
    mat_lo_tri<T>& operator =(const mat<T>& M) {
      q.resize(mat_triangular_size((M.getRowCount() > M.getColCount() ? M.getRowCount() : M.getColCount())),T(0.0));
      size = (M.getRowCount() > M.getColCount() ? M.getRowCount() : M.getColCount());
      unsigned int k=0;
      unsigned int i=0;
      unsigned int min_size = (M.getRowCount() > M.getColCount() ? M.getColCount() : M.getRowCount());
      for(;i<min_size;k += ++i) {
        for(unsigned int j=0;j<i;++j)
          q[k+j] = M(i,j);
        q[k+i] = M(i,i);
      };
      if(M.getRowCount() > M.getColCount()) {
        for(;i<size;k += ++i) {
             for(unsigned int j=0;j<min_size;++j)
            q[k+j] = M(i,j);
          for(unsigned int j=min_size;j<=i;++j)
            q[k+j] = 0.0;
        };
      } else
        for(;i<size;k += ++i)
          for(unsigned int j=0;j<=i;++j)
            q[k+j] = 0.0;
      return *this;
    };

    /**
     * Add-and-store operator with standard semantics.
     * \param M the other matrix to be added to this.
     * \return this matrix by reference.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T>& operator +=(const mat_lo_tri<T>& M) throw(std::range_error) {
      if(M.size != size)
        throw std::range_error("Matrix dimension mismatch.");
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        q[i] += M.q[i];
      return *this;
    };

    /**
     * Sub-and-store operator with standard semantics.
     * \param M the other matrix to be substracted from this.
     * \return this matrix by reference.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T>& operator -=(const mat_lo_tri<T>& M) throw(std::range_error) {
      if(M.size != size)
        throw std::range_error("Matrix dimension mismatch.");
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        q[i] -= M.q[i];
      return *this;
    };

    /**
     * Scalar-multiply-and-store operator with standard semantics.
     * \param S the scalar to be multiplied to this.
     * \return this matrix by reference.
     */
    mat_lo_tri<T>& operator *=(const T& S) {
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        q[i] *= S;
      return *this;
    };

    /**
     * Matrix-multiply-and-store operator with an lower-triangular.
     * \param M the other matrix to be multiplied with this.
     * \return this matrix by reference.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T>& operator *=(const mat_lo_tri<T>& M) throw(std::range_error) {
      if(size != M.size)
        throw std::range_error("Matrix dimension mismatch.");
      mat_lo_tri<T> result(size);
      unsigned int k=0;unsigned int i=0;
      for(;i<size;k += ++i) {
        unsigned int h=mat_triangular_size(i);
        for(unsigned int l=i;l<size;h += ++l) {
          for(unsigned int j=0;j<=i;++j)
            result.q[h+j] += M.q[k+j] * q[h+i];
        };
      };
      q.swap(result.q);
      return *this;
    };

/*******************************************************************************
                         Basic Operators
*******************************************************************************/

    /**
     * Addition operator with standard semantics.
     * \param M the other matrix to be added to this.
     * \return the matrix sum of this and M.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T> operator +(const mat_lo_tri<T>& M) const throw(std::range_error) {
      if(M.size != size)
        throw std::range_error("Matrix dimension mismatch.");
      mat_lo_tri<T> result(size);
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        result.q[i] = q[i] + M.q[i];
      return result;
    };

    /**
     * Negation operator with standard semantics.
     * \return the negative of this matrix sum.
     */
    mat_lo_tri<T> operator -() const {
      mat_lo_tri<T> result(size);
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        result.q[i] = -q[i];
      return result;
    };

    /**
     * Substraction operator with standard semantics.
     * \param M the other matrix to be substracted from this.
     * \return the matrix difference of this and M.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T> operator -(const mat_lo_tri<T>& M) const throw(std::range_error) {
      if(M.size != size)
        throw std::range_error("Matrix dimension mismatch.");
      mat_lo_tri<T> result(size);
      for(unsigned int i=0;i<mat_triangular_size(size);++i)
        result.q[i] = q[i] - M.q[i];
      return result;
    };

    /**
     * Matrix multiplication operator with standard semantics.
     * \param M the other matrix to be multiplied with this.
     * \return the matrix product of this and M.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_cm<T> operator *(const mat<T>& M) const throw(std::range_error) {
      if(size != M.getRowCount())
        throw std::range_error("Matrix dimension mismatch.");
      mat_cm<T> result(size,M.getColCount());
      unsigned int k=0;unsigned int i=0;
      for(;i<size;k += ++i)
        for(unsigned int l=0;l<M.getColCount();++l)
          for(unsigned int j=0;j<=i;++j)
            result.q[l*size+i] += q[k+j] * M(j,l);
      return result;
    };

    /**
     * Matrix multiplication operator with a lower-triangular.
     * \param M the other matrix to be multiplied with this.
     * \return the matrix product of this and M.
     * \throw std::range_error if the matrix dimensions don't match.
     */
    mat_lo_tri<T> operator *(const mat_lo_tri<T>& M) const throw(std::range_error) {
      if(size != M.size)
        throw std::range_error("Matrix dimension mismatch.");
      mat_lo_tri<T> result(size);
      unsigned int k=0;unsigned int i=0;
      for(;i<size;k += ++i) {
        unsigned int h=mat_triangular_size(i);
        for(unsigned int l=i;l<size;h += ++l) {
          for(unsigned int j=0;j<=i;++j)
            result.q[h+j] += M.q[k+j] * q[h+i];
        };
      };
      return result;
    };


/*******************************************************************************
                         Special Methods
*******************************************************************************/

    /**
     * Extracts a sub-matrix from this matrix.
     * \param aRowOffset Number of rows before the start of the sub-matrix rows.
     * \param aColOffset Number of columns before the start of the sub-matrix columns.
     * \param aRowCountOut Number of rows of the sub-matrix.
     * \param aColCountOut Number of columns of the sub-matrix.
     * \return The sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    mat_cm<T> getSubMat(unsigned int aRowOffset,unsigned int aColOffset,unsigned int aRowCountOut,unsigned int aColCountOut) const throw(std::range_error) {
      if((aRowOffset + aRowCountOut > size) || (aColOffset + aColCountOut > size))
        throw std::range_error("Matrix dimension mismatch.");
      mat_cm<T> result(aRowCountOut,aColCountOut);
      unsigned int k=mat_triangular_size(aColOffset);
      for(unsigned int j=0;j<aColCountOut;k += (++j + aColOffset)) {
        unsigned int h=mat_triangular_size(aRowOffset);unsigned int i=0;
        for(;((i<aRowCountOut) && (i+aRowOffset < j+aColOffset));h += (++i + aRowOffset))
          result.q[j*aRowCountOut+i] = 0.0;
        for(;i<aRowCountOut;h += (++i + aRowOffset)) {
          result.q[j*aRowCountOut+i] = q[h+j+aColOffset];
        };
      };
      return result;
    };

    /**
     * Extracts a lower-triangular sub-matrix from this matrix.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \param aSizeOut Number of rows/columns of the sub-matrix.
     * \return The lower-triangular sub-matrix contained in this matrix.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    mat_lo_tri<T> getSubMat(unsigned int aDiagOffset,unsigned int aSizeOut) const throw(std::range_error) {
      if(aDiagOffset + aSizeOut > size)
        throw std::range_error("Matrix dimension mismatch.");
      mat_lo_tri<T> result(aSizeOut);
      unsigned int k=mat_triangular_size(aDiagOffset);
      unsigned int k_out=0;
      for(unsigned int i=0;i<aSizeOut;k_out += ++i, k += (i + aDiagOffset)) {
        for(unsigned int j=0;j<=i;++j) {
          result.q[k_out+j] = q[k+j+aDiagOffset];
        };
      };
      return result;
    };


    /** Sets the sub-part of this matrix to a lower-triangular sub-matrix M.
     * \param M An lower-triangular sub-matrix that will be written in the sub-part of this matrix.
     * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
     * \return This matrix, by reference.
     * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
     */
    mat_lo_tri<T>& setSubMat(const mat_lo_tri<T>& M,unsigned int aDiagOffset) throw(std::range_error) {
      if(aDiagOffset + M.size > size)
        throw std::range_error("Matrix dimension mismatch.");
      unsigned int k=mat_triangular_size(aDiagOffset);
      unsigned int k_in=0;
      for(unsigned int i=0;i<M.size;k_in += ++i, k += (i + aDiagOffset)) {
        for(unsigned int j=0;j<=i;++j) {
           q[k+j+aDiagOffset] = M.q[k_in+j];
        };
      };
      return *this;
    };

};




#endif


};


#endif











