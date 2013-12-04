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

#ifndef REAK_COVAR_TOPOLOGY_HPP
#define REAK_COVAR_TOPOLOGY_HPP

#include "base/defs.hpp"

#include "covariance_concept.hpp"
#include "lin_alg/mat_norms.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_exp_methods.hpp"
#include "path_planning/global_rng.hpp"
#include "path_planning/metric_space_concept.hpp"
#include "topologies/default_random_sampler.hpp"

namespace ReaK {

namespace ctrl {


/**
 * This class template creates a topology for a covariance 
 * matrix type. Many algorithms in ReaK::ctrl and ReaK::pp require topologies to 
 * represent the space in which a point-type can exist, this class simply implements
 * a topology for a covariance matrix which can be used in a Gaussian belief-state 
 * topology (see gaussian_belief_space) or elsewhere.
 * 
 * Models: TopologyConcept, MetricSpaceConcept, LieGroupConcept, and PointDistributionConcept.
 * 
 * \tparam Covariance A covariance matrix type, modeling CovarianceMatrixConcept.
 */
template <typename Covariance>
class covar_topology {
  public:
    typedef covar_topology<Covariance> self;
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
    
    
    struct distance_metric_type {
      
      double operator()(const point_difference_type& dp, const self&) const {
        return norm(dp);
      };
      
      double operator()(const point_type& a, const point_type& b, const self& s) const {
        return this->operator()(point_difference_type(a,b),s);
      };
      
    };
    
    typedef ReaK::pp::default_random_sampler random_sampler_type;
    
    
  private:
    value_type max_eigenvalue;
    size_type mat_size;
    
  public:
    
    /**
     * Parametrized constructor with default random-number generator.
     * \param aSize The size of the covariance matrix.
     * \param aMaxEigenValue The maximum eigen-value of the random covariance matrices.
     */
    explicit covar_topology(size_type aSize = 0, const value_type& aMaxEigenValue = value_type(1.0)) 
      : max_eigenvalue(aMaxEigenValue), mat_size(aSize) { };
                     
    /**
     * Generates a random covariance matrix.
     * \return A random covariance matrix.
     */
    point_type random_point() const 
    {
      boost::uniform_01<pp::global_rng_type&, double> uniform_rng(pp::get_global_rng());
      mat<value_type,mat_structure::diagonal> D(mat_size);
      for(size_type i = 0; i < mat_size; ++i)
        D(i,i) = uniform_rng() * max_eigenvalue;
      
      mat<value_type,mat_structure::skew_symmetric> S(mat_size);
      for(size_type i = 1; i < mat_size; ++i)
        for(size_type j = 0; j < i; ++j)
          S(j,i) = uniform_rng() * value_type(10.0);
      mat<value_type,mat_structure::square> Q(S);
      exp_PadeSAS(S,Q,QR_linlsqsolver());
      
      return point_type(matrix_type(transpose_view(Q) * (D * Q)));
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
     * Computes the covariance matrix which is the linear interpolation between two matrices, by a fraction.
     * \param a The first covariance matrix.
     * \param fraction The scalar fraction at which the evaluate the linear interpolation.
     * \param b The second covariance matrix.
     * \return The interpolated covariance matrix.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return point_type(matrix_type(value_type(1.0 - fraction) * a.get_matrix() + value_type(fraction) * b.get_matrix()));
    };
    
    
    
};


template <typename Covariance>
typename covar_topology<Covariance>::distance_metric_type get(ReaK::pp::distance_metric_t, const covar_topology<Covariance>&) {
  return typename covar_topology<Covariance>::distance_metric_type();
};


};

};

#endif


















