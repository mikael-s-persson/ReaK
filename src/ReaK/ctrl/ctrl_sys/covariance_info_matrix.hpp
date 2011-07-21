
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

#ifndef INFORMATION_MATRIX_HPP
#define INFORMATION_MATRIX_HPP

#include "math/mat_alg.hpp"
#include "math/mat_cholesky.hpp"
#include "base/named_object.hpp"

#include "covariance_concept.hpp"


namespace ReaK {

namespace ctrl {



template <typename StateType>
class covariance_info_matrix : public named_object {
  public:
    typedef covariance_info_matrix<T> self;
    
    typedef StateType point_type;
    typedef typename state_vector_traits<StateType>::state_difference_type point_difference_type;
    
    typedef typename state_vector_traits<StateType>::value_type value_type;
    typedef mat<value_type, mat_structure::symmetric> matrix_type;
    typedef typename matrix_type::size_type size_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = vect_traits<point_difference_type>::dimensions);
    BOOST_STATIC_CONSTANT(covariance_storage::tag, storage = covariance_storage::information);
    
  private:
    matrix_type mat_info;
    
  public:
    
    explicit covariance_info_matrix(const matrix_type& aMat = matrix_type(), 
				    const std::string& aName = "") : mat_info(aMat) { 
      setName(aName); 
      invert_Cholesky(aMat, mat_info, std::numeric_limits< value_type >::epsilon());
    };
    
    explicit covariance_info_matrix(size_type aSize, 
			        covariance_initial_level aLevel = covariance_initial_level::full_info, 
			        const std::string& aName = "") : 
			        mat_info(aSize, value_type( ( aLevel == covariance_initial_level::no_info ? 0 : std::numeric_limits< value_type >::infinity() ) )) { 
      setName(aName); 
    };
    
    matrix_type get_matrix() const { 
      matrix_type m_inv; 
      invert_Cholesky(mat_info, m_inv, std::numeric_limits< value_type >::epsilon()); 
      return m_inv;
    };
    const matrix_type& get_inverse_matrix() const { 
      return mat_info;
    };
      
    friend 
    void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mat_info,rhs.mat_info);
    };
    
    self& operator =(self rhs) {
      swap(rhs,*this);
      return *this;
    };
    
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
    self& >::type operator =(const Matrix& rhs) {
      invert_Cholesky(rhs, mat_info, std::numeric_limits< value_type >::epsilon());
      return *this;
    };
    
    operator matrix_type() const { return get_matrix(); };
    
    friend 
    const matrix_type& invert(const self& aObj) {
      return aObj.get_inverse_matrix();
    };
    
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(mat_info);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(mat_info);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300009,1,"covariance_info_matrix",named_object)
    
};





};

};

#endif

















