/**
 * \file gaussian_belief_space.hpp
 * 
 * This library provides a class template which can join together a state-vector topology
 * and a covariance matrix topology (see covar_topology.hpp) to create a Gaussian belief-state
 * topology. The distance pseudo-metric used is the symmetric KL-divergence between two 
 * belief-states.
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

#ifndef REAK_GAUSSIAN_BELIEF_SPACE_HPP
#define REAK_GAUSSIAN_BELIEF_SPACE_HPP

#include "gaussian_belief_state.hpp"
#include "path_planning/metric_space_concept.hpp"

namespace ReaK {

namespace ctrl {



/**
 * This class template can join together a state-vector topology
 * and a covariance matrix topology (see covar_topology.hpp) to create a Gaussian belief-state
 * topology. The distance pseudo-metric used is the symmetric KL-divergence between two 
 * belief-states.
 * 
 * Models: ReaK::pp::MetricSpaceConcept.
 * 
 * \tparam StateTopology The topology to which the state-vector belongs, should model the ReaK::pp::MetricSpaceConcept.
 * \tparam CovarianceTopology The topology to which the covariance matrix belongs, should model the ReaK::pp::MetricSpaceConcept.
 */
template <typename StateTopology, typename CovarianceTopology>
class gaussian_belief_space {
  public:
    typedef gaussian_belief_space< StateTopology, CovarianceTopology > self;
    typedef StateTopology mean_state_topology;
    typedef CovarianceTopology covariance_topology;
    
    typedef typename pp::metric_topology_traits<CovarianceTopology>::point_type covariance_type;
    typedef typename pp::metric_topology_traits<CovarianceTopology>::point_difference_type covariance_diff_type;
    typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
    typedef typename covariance_mat_traits<covariance_type>::value_type value_type;
    
    typedef typename pp::metric_topology_traits<StateTopology>::point_type mean_state_type;
    typedef typename pp::metric_topology_traits<StateTopology>::point_difference_type mean_state_diff_type;
    
    typedef typename gaussian_belief_state< covariance_type, mean_state_type > point_type;
    
    BOOST_CONCEPT_ASSERT((pp::MetricSpaceConcept<StateTopology>));
    BOOST_CONCEPT_ASSERT((pp::MetricSpaceConcept<CovarianceTopology>));

    /**
     * This nested class represents the difference between two belief-states (with some fraction).
     */
    struct point_difference_type {
      point_type b0;
      point_type b1;
      value_type fraction;
            
      point_difference_type(const point_type& aB0, 
			    const point_type& aB1,
			    const value_type& aFraction = value_type(1.0)) :
                            b0(aB0), b1(aB1), fraction(aFraction) { };
      
      point_difference_type operator-() const {
        return point_difference_type( b1, b0, fraction );
      };

      friend point_difference_type operator*(const point_difference_type& a, double b) {
        return point_difference_type(a.b0, a.b1, a.fraction * b);
      };

      friend point_difference_type operator*(double a, const point_difference_type& b) {
        return point_difference_type(b.b0, b.b1, b.fraction * a);
      };
      
    };
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
  private:
    mean_state_topology mean_state_space;
    covariance_topology covariance_space;
    
  public:  
    
    /**
     * Parametric and default constructor.
     * \param aMeanStateSpace The topology used for the mean-state.
     * \param aCovarianceSpace The topology used for the covariance matrix.
     */
    gaussian_belief_space(const mean_state_topology& aMeanStateSpace = mean_state_topology(),
                          const covariance_topology& aCovarianceSpace = covariance_topology()) :
			  mean_state_space(aMeanStateSpace),
			  covariance_space(aCovarianceSpace) { };
			  
    
    
    /**
     * Computes the distance between two belief-states, using symmetric KL-divergence.
     * \param p1 The first belief-state.
     * \param p2 The second belief-state.
     * \return The symmetric KL-divergence between the two belief-states.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      return double( symKL_divergence(p1,p2) );
    };
    
    /**
     * Computes the norm of a belief-state difference, using symmetric KL-divergence.
     * \param dp The belief-state difference.
     * \return The symmetric KL-divergence between the two end belief-states.
     */
    double norm(const point_difference_type& dp) const {
      return double( symKL_divergence(dp.b0,dp.b1) * dp.fraction );
    };
    
    /**
     * Computes a random belief-state from the underlying state and covariance topologies.
     * \return a random belief-state from the underlying state and covariance topologies.
     */
    point_type random_point() const {
      return point_type( mean_state_space.random_point(),
			 covariance_space.random_point());
    };
    
    /**
     * Computes the difference between two belief-states.
     * \param p1 The first belief-state.
     * \param p2 The second belief-state.
     * \return The difference between the two belief-states.
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return point_difference_type(p1,p2);
    };
    
    /**
     * Computes the belief-state which lies a fraction between two belief-states (also works for 
     * a fraction that is not between 0 and 1).
     * \param p1 The starting belief-state (d = 0).
     * \param d The fraction between the two belief-states.
     * \param p2 The final belief-state (d = 1).
     * \return The belief-state which lies a fraction between two belief-states.
     */
    point_type move_position_toward(const point_type& p1, double d, const point_type& p2) const {
      return point_type( mean_state_space.move_position_toward(p1.get_mean_state(), d, p2.get_mean_state()),
			 covariance_space.move_position_toward(p1.get_covariance(), d, p2.get_covariance()) );
    };
    
    /**
     * Computes the origin of the belief-space.
     * \return the origin of the belief-space, as a belief-state with the mean-state at the origin of the mean-state topology and the covariance at the origin of the covariance topology.
     */
    point_type origin() const {
      return point_type( mean_state_space.origin(), covariance_space.origin() );
    };
    
    /**
     * Adjusts a belief-state with a belief-state difference.
     * \param p1 The starting belief-state.
     * \param dp The belief-state difference to add to the starting belief-state.
     * \return The adjusted belief-state (semantically p1 + dp).
     */
    point_type adjust(const point_type& p1, const point_difference_type& dp) const {
      return point_type( mean_state_space.adjust( p1.get_mean_state(),
						  dp.fraction * mean_state_space.difference( dp.b0.get_mean_state(),
											     dp.b1.get_mean_state() ) ),
			 covariance_space.adjust( p1.get_covariance(),
						  dp.fraction * covariance_space.difference( dp.b0.get_covariance(),
											     dp.b1.get_covariance() ) ) );
    };
    
};




};

};

#endif











