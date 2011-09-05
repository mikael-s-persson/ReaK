/**
 * \file covariance_info_matrix.hpp
 * 
 * This library provides a class template to represent a covariance matrix as an 
 * information matrix (inverse of the covariance matrix).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_INFORMATION_MATRIX_HPP
#define REAK_INFORMATION_MATRIX_HPP

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_cholesky.hpp"
#include "base/named_object.hpp"

#include "covariance_concept.hpp"


namespace ReaK {

namespace ctrl {


/**
 * This class template can represent a covariance matrix (as of the CovarianceMatrixConcept) 
 * by containing the information matrix. For some algorithms, only the inverse of the covariance
 * matrix is needed and it is more efficient to store the covariance matrix by its inverse, 
 * the information matrix. The information matrix also has the advantage of being able to 
 * represent the no-information case (corresponding to an infinite covariance matrix). Outputing 
 * the information matrix induces no computation while outputing the covariance matrix induces
 * a matrix inversion.
 * 
 * Models: CovarianceMatrixConcept
 * 
 * \tparam StateType The state-vector type which the covariance matrix is the covariance of, should model StateVectorConcept.
 */
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
    
    /**
     * Parametrized constructor.
     * \param aMat The covariance matrix to which to initialize this object.
     */
    explicit covariance_info_matrix(const matrix_type& aMat = matrix_type(), 
				    const std::string& aName = "") : mat_info(aMat) { 
      setName(aName); 
      invert_Cholesky(aMat, mat_info, std::numeric_limits< value_type >::epsilon());
    };
    
    /**
     * Parametrized constructor.
     * \param aSize The size of the covariance matrix.
     * \param aLevel The information level to initialize this object with.
     */
    explicit covariance_info_matrix(size_type aSize, 
			        covariance_initial_level aLevel = covariance_initial_level::full_info, 
			        const std::string& aName = "") : 
			        mat_info(aSize, value_type( ( aLevel == covariance_initial_level::no_info ? 0 : std::numeric_limits< value_type >::infinity() ) )) { 
      setName(aName); 
    };
    
    /**
     * Returns the covariance matrix (as a matrix object).
     * \return The covariance matrix (as a matrix object).
     */
    matrix_type get_matrix() const { 
      matrix_type m_inv; 
      invert_Cholesky(mat_info, m_inv, std::numeric_limits< value_type >::epsilon()); 
      return m_inv;
    };
    /**
     * Returns the inverse covariance matrix (information matrix) (as a matrix object).
     * \return The inverse covariance matrix (information matrix) (as a matrix object).
     */
    const matrix_type& get_inverse_matrix() const { 
      return mat_info;
    };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.mat_info,rhs.mat_info);
    };
    
    /**
     * Standard assignment operator.
     */
    self& operator =(self rhs) {
      swap(rhs,*this);
      return *this;
    };
    
    /**
     * Assignment to a readable matrix (covariance matrix).
     */
    template <typename Matrix>
    typename boost::enable_if_c< is_readable_matrix<Matrix>::value, 
    self& >::type operator =(const Matrix& rhs) {
      invert_Cholesky(rhs, mat_info, std::numeric_limits< value_type >::epsilon());
      return *this;
    };
    
    /**
     * Implicit conversion to a covariance matrix type.
     */
    operator matrix_type() const { return get_matrix(); };
    
    /**
     * Conversion to an information matrix type.
     */
    friend const matrix_type& invert(const self& aObj) {
      return aObj.get_inverse_matrix();
    };
    
    /**
     * Returns the size of the covariance matrix.
     * \return The size of the covariance matrix.
     */
    size_type size() const { return mat_info.get_row_count(); };
    
    
    
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

















