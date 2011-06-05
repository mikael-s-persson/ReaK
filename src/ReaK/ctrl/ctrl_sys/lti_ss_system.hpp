
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

#ifndef LTI_SS_SYSTEM_HPP
#define LTI_SS_SYSTEM_HPP

#include "linear_ss_system_concept.hpp"

#include "base/named_object.hpp"

#include "math/mat_alg.hpp"
#include "math/vect_alg.hpp"

namespace ReaK {

namespace ctrl {


template <typename T>
class lti_system_ss : public named_object {
  private:
    mat<T,mat_structure::square> A;
    mat<T,mat_structure::rectangular> B;
    mat<T,mat_structure::rectangular> C;
    mat<T,mat_structure::rectangular> D;
    
  public:
    typedef lti_system_ss<T> self;
    typedef T value_type;
    typedef std::size_t size_type;
    
    typedef vect_n<T> point_type;
    typedef vect_n<T> point_difference_type;
    typedef vect_n<T> point_derivative_type;
    typedef self topology;
  
    typedef T time_type;
    typedef T time_difference_type;
  
    typedef vect_n<T> input_type;
    typedef vect_n<T> output_type;
    
    typedef mat<T,mat_structure::square> matrixA_type;
    typedef mat<T,mat_structure::rectangular> matrixB_type;
    typedef mat<T,mat_structure::rectangular> matrixC_type;
    typedef mat<T,mat_structure::rectangular> matrixD_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 0);
        
    lti_system_ss(size_type aX_size = 0, size_type aU_size = 0, size_type aY_size = 0, const std::string& aName = "") :
                  A(aX_size), B(aX_size,aU_size), C(aY_size,aX_size), D(aY_size,aU_size) {
      setName(aName);
    };

    lti_system_ss(const self& rhs) : A(rhs.A), B(rhs.B), C(rhs.C), D(rhs.D) {
      setName(rhs.getName());
    };
    
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    lti_system_ss(const MatrixA& aA, const MatrixB& aB, const MatrixC& aC, const MatrixD& aD, const std::string& aName = "",
                  typename boost::enable_if_c< is_readable_matrix<MatrixA>::value &&
                                               is_readable_matrix<MatrixB>::value &&
                                               is_readable_matrix<MatrixC>::value &&
                                               is_readable_matrix<MatrixD>::value , void* >::type dummy = NULL) : 
                  A(aA), B(aB), C(aC), D(aD) { 
      setName(aName);
    };
  
  
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.A,rhs.A);
      swap(lhs.B,rhs.B);
      swap(lhs.C,rhs.C);
      swap(lhs.D,rhs.D);
    };
		  
    self& operator =(self rhs) {
      swap(*this,rhs);
      return *this;
    };
		  
    template <typename MatrixA>
    typename boost::enable_if_c< is_readable_matrix<MatrixA>::value,
    void >::type setA(const MatrixA& aA) {
      A = aA;
    };
    
    template <typename MatrixB>
    typename boost::enable_if_c< is_readable_matrix<MatrixB>::value,
    void >::type setB(const MatrixB& aB) {
      B = aB;
    };
    
    template <typename MatrixC>
    typename boost::enable_if_c< is_readable_matrix<MatrixC>::value,
    void >::type setC(const MatrixC& aC) {
      C = aC;
    };
    
    template <typename MatrixD>
    typename boost::enable_if_c< is_readable_matrix<MatrixD>::value,
    void >::type setD(const MatrixD& aD) {
      D = aD;
    };
    
    const matrixA_type& getA() const { return A; };
    const matrixB_type& getB() const { return B; };
    const matrixC_type& getC() const { return C; };
    const matrixD_type& getD() const { return D; };
    
    size_type get_state_count() const { return A.get_col_count(); };
    size_type get_input_count() const { return B.get_col_count(); };
    size_type get_output_count() const { return C.get_row_count(); };

    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD) const {
      aA = A;
      aB = B;
      aC = C;
      aD = D;
    };
    
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(const time_type&, MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD) const {
      aA = A;
      aB = B;
      aC = C;
      aD = D;
    };
    
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(const point_type&, const input_type&, const time_type&, MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD) const {
      aA = A;
      aB = B;
      aC = C;
      aD = D;
    };
    
    
    point_derivative_type get_state_derivative(const point_type& p, const input_type& u, const time_type& t = 0) const {
      return A * p + B * u;
    };
    
    output_type get_output(const point_type& p, const input_type& u, const time_type& t = 0) const {
      return C * p + D * u;
    };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return p + dp;
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(A)
         & RK_SERIAL_SAVE_WITH_NAME(B)
         & RK_SERIAL_SAVE_WITH_NAME(C)
         & RK_SERIAL_SAVE_WITH_NAME(D);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(A)
         & RK_SERIAL_LOAD_WITH_NAME(B)
         & RK_SERIAL_LOAD_WITH_NAME(C)
         & RK_SERIAL_LOAD_WITH_NAME(D);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300001,1,"lti_system_ss",named_object)
  
};
  
  



};

};

#endif






