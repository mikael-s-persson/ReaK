
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

#ifndef DECOMP_COVARIANCE_MATRIX_HPP
#define DECOMP_COVARIANCE_MATRIX_HPP

#include "math/mat_alg.hpp"
#include "base/named_object.hpp"

#include "covariance_concept.hpp"
#include "math/mat_gaussian_elim.hpp"
#include "math/mat_qr_decomp.hpp"


namespace ReaK {

namespace ctrl {



template <typename T>
class decomp_covariance_matrix : public named_object {
  public:
    typedef decomp_covariance_matrix<T> self;
    
    typedef T value_type;
    typedef std::size_t size_type;
    
    typedef vect_n<T> point_type;
    typedef vect_n<T> point_difference_type;
    
    typedef mat<T, mat_structure::symmetric> matrix_type;
    typedef mat<T, mat_structure::square> matrix_block_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(covariance_storage::tag, storage = covariance_storage::other);
    
  private:
    matrix_block_type mat_X;
    matrix_block_type mat_Y;
    
  public:
    
    explicit decomp_covariance_matrix(const matrix_block_type& aMatX = matrix_block_type(), 
				      const matrix_block_type& aMatY = matrix_block_type(), 
				      const std::string& aName = "") : 
				      mat_X(aMatX), 
				      mat_Y(aMatY) { 
      setName(aName); 
    };
    
    explicit decomp_covariance_matrix(size_type aSize, 
				      covariance_initial_level aLevel = covariance_initial_level::full_info, 
				      const std::string& aName = "") : 
				      mat_X(aSize, value_type(0)), 
				      mat_Y(aSize, value_type(0)) { 
      setName(aName); 
      if(aLevel == covariance_initial_level::full_info) {
	mat_Y = mat<value_type, mat_structure::identity>(aSize);
      } else {
	mat_X = mat<value_type, mat_structure::identity>(aSize);
      };
    };
    
    matrix_type get_matrix() const { 
      try {
	mat<value_type, mat_structure::square> mY_inv;
	invert_PLU(mat_Y,mY_inv);
	return matrix_type(mat_X * mY_inv);
      } catch(singularity_error& e) {
	mat<value_type, mat_structure::square> mY_inv;
	pseudoinvert_QR(mat_Y,mY_inv);
	return matrix_type(mat_X * mY_inv);
      }; 
    };
    matrix_type get_inverse_matrix() const { 
      try {
	mat<value_type, mat_structure::square> mX_inv;
	invert_PLU(mat_X,mX_inv);
	return matrix_type(mat_Y * mX_inv);
      } catch(singularity_error& e) {
	mat<value_type, mat_structure::square> mX_inv;
	pseudoinvert_QR(mat_X,mX_inv);
	return matrix_type(mat_Y * mX_inv);
      }; 
    };
    
    const matrix_block_type& get_covarying_block() const { return mat_X; };
    const matrix_block_type& get_informing_inv_block() const { return mat_Y; };
    
    friend 
    void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mat_X,rhs.mat_X);
      swap(lhs.mat_Y,rhs.mat_Y);
    };
    
    self& operator =(self rhs) {
      swap(rhs,*this);
      return *this;
    };
        
    operator matrix_type() const { return get_matrix(); };
    
    friend 
    matrix_type invert(const self& aObj) {
      return aObj.get_inverse_matrix();
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(mat_X)
         & RK_SERIAL_SAVE_WITH_NAME(mat_Y);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(mat_X)
         & RK_SERIAL_LOAD_WITH_NAME(mat_Y);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC230000A,1,"decomp_covariance_matrix",named_object)
    
};




  

};

};


#endif











