
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

#ifndef COVARIANCE_MATRIX_HPP
#define COVARIANCE_MATRIX_HPP

#include "math/mat_alg.hpp"
#include "math/mat_cholesky.hpp"
#include "base/named_object.hpp"

#include "covariance_concept.hpp"


namespace ReaK {

namespace ctrl {



template <typename T>
class covariance_matrix : public named_object {
  public:
    typedef covariance_matrix<T> self;
    
    typedef T value_type;
    typedef std::size_t size_type;
    
    typedef vect_n<T> point_type;
    typedef vect_n<T> point_difference_type;
    
    typedef mat<T, mat_structure::symmetric> matrix_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(covariance_storage::tag, storage = covariance_storage::covariance);
    
  private:
    matrix_type mat_cov;
    
  public:
    
    explicit covariance_matrix(const matrix_type& aMat = matrix_type(), const std::string& aName = "") : mat_cov(aMat) { setName(aName); };
    
    explicit covariance_matrix(size_type aSize, 
			       covariance_initial_level::tag aLevel = covariance_initial_level::full_info, 
			       const std::string& aName = "") : 
			       mat_cov(aSize, value_type( ( aLevel == covariance_initial_level::full_info ? 0 : std::numeric_limits< value_type >::infinity() ) )) { 
      setName(aName); 
    };
    
    const matrix_type& get_matrix() const { return mat_cov; };
    matrix_type get_inverse_matrix() const { 
      matrix_type m_inv; 
      invert_Cholesky(mat_cov, m_inv, std::numeric_limits< value_type >::epsilon()); 
      return m_inv;
    };
    
    friend 
    void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mat_cov,rhs.mat_cov);
    };
    
    self& operator =(self rhs) {
      swap(rhs,*this);
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
    self& >::type operator =(const Matrix& rhs) {
      mat_cov = rhs;
      return *this;
    };
    
    operator const matrix_type&() const { return mat_cov; };
    
    friend 
    matrix_type invert(const self& aObj) {
      return aObj.get_inverse_matrix();
    };
    
    size_type size() const { return mat_cov.get_row_count(); };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(mat_cov);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(mat_cov);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300008,1,"covariance_matrix",named_object)
    
};





};

};

#endif









