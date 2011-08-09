/**
 * \file covar_topology.hpp
 * 
 * This library provides a class template which creates a topology for a covariance 
 * matrix type. Many algorithms in ReaK::ctrl and ReaK::pp require topologies to 
 * represent the space in which a point-type can exist, this library simply implements
 * a topology for a covariance matrix which can be used in a Gaussian belief-state 
 * topology (see gaussian_belief_space).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef COVAR_TOPOLOGY_HPP
#define COVAR_TOPOLOGY_HPP

#include "covariance_concept.hpp"
#include "math/mat_norms.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT
#include <boost/shared_ptr.hpp>

#include "math/mat_qr_decomp.hpp"
#include "math/mat_exp_methods.hpp"

namespace ReaK {

namespace ctrl {


/**
 * This class template creates a topology for a covariance 
 * matrix type. Many algorithms in ReaK::ctrl and ReaK::pp require topologies to 
 * represent the space in which a point-type can exist, this class simply implements
 * a topology for a covariance matrix which can be used in a Gaussian belief-state 
 * topology (see gaussian_belief_space) or elsewhere.
 * 
 * Models: MetricSpaceConcept.
 * 
 * \tparam Covariance A covariance matrix type, modeling CovarianceMatrixConcept.
 * \tparam RandomNumberGenerator A standard random number generator functor type (e.g. boost::minstd_rand).
 */
template <typename Covariance, typename RandomNumberGenerator = boost::minstd_rand>
class covar_topology {
  typedef boost::uniform_01<RandomNumberGenerator, double> rand_t;

  public:
    typedef covar_topology<Covariance,RandomNumberGenerator> self;
    typedef Covariance point_type; 
    typedef typename covariance_mat_traits<Covariance>::matrix_type matrix_type;
    typedef typename mat_traits<matrix_type>::size_type size_type;
    typedef typename covariance_mat_traits<Covariance>::value_type value_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = point_type::dimensions);
    
    /**
     * This nested class implements the point-difference type for the covariance topology.
     * It implements a simple linear interpolation of the covariance matrices and uses 
     * the Frobenius norm as a metric. It also implements all required arithmetic operators
     * with standard semantics.
     */
    struct point_difference_type {
      matrix_type M;
      
      point_difference_type() { };
      
      point_difference_type(const point_type& aSrc, 
			    const point_type& aDst) :
                            M(aSrc.get_matrix() - aDst.get_matrix()) { };
      
      friend value_type norm(const self& dp) {
	return frobenius_norm(dp.M);
      };
      
      friend point_difference_type operator+(const point_difference_type& a, const point_difference_type& b) {
        point_difference_type result;
	result.M = a.M + b.M;
        return result;
      };

      point_difference_type& operator+=(const point_difference_type& b) {
        M += b.M;
        return *this;
      };

      point_difference_type operator-() const {
        point_difference_type result;
        result.M = -M;
        return result;
      };

      friend point_difference_type operator-(const point_difference_type& a, const point_difference_type& b) {
        point_difference_type result;
	result.M = a.M - b.M;
	return result;
      };

      point_difference_type& operator-=(const point_difference_type& b) {
        M -= b.M;
        return *this;
      };

      friend point_difference_type operator*(point_difference_type a, double b) {
        a.M *= b;
	return a;
      };

      friend point_difference_type operator*(double a, point_difference_type b) {
        b.M *= a;
	return b;
      };
      
      
    };
    
    
  private:
    boost::shared_ptr<RandomNumberGenerator> gen_ptr;
    boost::shared_ptr<rand_t> rand;
    value_type max_eigenvalue;
    size_type mat_size;
    
  public:
    
    /**
     * Parametrized constructor with default random-number generator.
     * \param aSize The size of the covariance matrix.
     * \param aMaxEigenValue The maximum eigen-value of the random covariance matrices.
     */
    explicit covar_topology(size_type aSize = 0, const value_type& aMaxEigenValue = value_type(1.0)) 
      : gen_ptr(new RandomNumberGenerator), rand(new rand_t(*gen_ptr)), 
        max_eigenvalue(aMaxEigenValue), mat_size(aSize) { };

    /**
     * Parametrized constructor with default random-number generator.
     * \param aGen The random-number generator used to generate random covariance matrices.
     * \param aSize The size of the covariance matrix.
     * \param aMaxEigenValue The maximum eigen-value of the random covariance matrices.
     */
    explicit covar_topology(RandomNumberGenerator& aGen, size_type aSize = 0, const value_type& aMaxEigenValue = value_type(1.0)) 
      : gen_ptr(), rand(new rand_t(aGen)), max_eigenvalue(aMaxEigenValue), mat_size(aSize) { };
                     
    /**
     * Generates a random covariance matrix.
     * \return A random covariance matrix.
     */
    point_type random_point() const 
    {
      mat<value_type,mat_structure::diagonal> D(mat_size);
      for(size_type i = 0; i < mat_size; ++i)
	D(i,i) = (*rand)() * max_eigenvalue;
      
      mat<value_type,mat_structure::skew_symmetric> S(mat_size);
      for(size_type i = 1; i < mat_size; ++i)
	for(size_type j = 0; j < i; ++j)
	  S(j,i) = (*rand)() * value_type(10.0);
      mat<value_type,mat_structure::square> Q(S);
      exp_PadeSAS(S,Q,QR_linlsqsolver());
      
      return point_type(matrix_type(transpose(Q) * (D * Q)));
    }
    
    /**
     * Computes the distance between to covariance matrices, implemented as the Forbenius norm of the 
     * difference between the two matrices.
     * \param a The first covariance matrix.
     * \param b The second covariance matrix.
     * \return The distance between the covariance matrices.
     */
    double distance(const point_type& a, const point_type& b) const {
      return norm(point_difference_type(a,b));
    };

    /**
     * Computes the covariance matrix which is the linear interpolation between two matrices, by a fraction.
     * \param a The first covariance matrix.
     * \param fraction The scalar fraction at which the evaluate the linear interpolation.
     * \param b The second covariance matrix.
     * \return The interpolated covariance matrix.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return point_type(matrix_type(value_type(1.0 - fraction) * a.get_matrix() + value_type(fraction) * b.get_matrix()));
    };

    /**
     * Computes the difference between two covariance matrices.
     * \param a The first covariance matrix.
     * \param b The second covariance matrix.
     * \return The difference between two covariance matrices
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      return point_difference_type(a,b);
    };

    /**
     * Adds a given covariance matrix difference to a given covariance matrix.
     * \param a A covariance matrix.
     * \param delta A covariance matrix difference to add to a.
     * \return The adjusted covariance matrix.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return point_type(matrix_type( a.get_matrix() + delta.M ));
    };
  
    /**
     * Returns the origin of the topology (a nil covariance matrix).
     * \return The origin of the topology (a nil covariance matrix).
     */
    point_type origin() const {
      return point_type(matrix_type( mat<value_type,mat_structure::nil>(mat_size) ));
    };

    /**
     * Returns the norm of a covariance matrix difference.
     * \param a The covariance matrix difference.
     * \return The norm of a covariance matrix difference (Frobenius norm).
     */
    double norm(const point_difference_type& a) const {
      return ::ReaK::ctrl::norm(a);
    };
    
    
    
    
    
};


};

};

#endif


















