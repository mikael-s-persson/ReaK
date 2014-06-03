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

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/default_random_sampler.hpp>

#include "gaussian_belief_state.hpp"

#include <boost/mpl/and.hpp>

namespace ReaK {

namespace ctrl {



/**
 * This class template can join together a state-vector topology
 * and a covariance matrix topology (see covar_topology.hpp) to create a Gaussian belief-state
 * topology. The distance pseudo-metric used is the symmetric KL-divergence between two 
 * belief-states.
 * 
 * Models: TopologyConcept, MetricSpaceConcept, PointDistributionConcept, and ContinuousBeliefSpaceConcept.
 * 
 * \tparam StateTopology The topology to which the state-vector belongs, should model the ReaK::pp::TopologyConcept.
 * \tparam CovarianceTopology The topology to which the covariance matrix belongs, should model the ReaK::pp::TopologyConcept.
 */
template <typename StateTopology, typename CovarianceTopology>
class gaussian_belief_space : public named_object {
  public:
    typedef gaussian_belief_space< StateTopology, CovarianceTopology > self;
    typedef StateTopology state_topology;
    typedef CovarianceTopology covariance_topology;
    
    typedef typename pp::topology_traits<CovarianceTopology>::point_type covariance_type;
    typedef typename pp::topology_traits<CovarianceTopology>::point_difference_type covariance_diff_type;
    typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
    typedef typename covariance_mat_traits<covariance_type>::value_type value_type;
    
    typedef typename pp::topology_traits<StateTopology>::point_type mean_state_type;
    typedef typename pp::topology_traits<StateTopology>::point_difference_type mean_state_diff_type;
    
    typedef gaussian_belief_state< mean_state_type, covariance_type > point_type;
    typedef std::pair< point_type, point_type > point_difference_type;
    
    BOOST_CONCEPT_ASSERT((pp::TopologyConcept<StateTopology>));
    BOOST_CONCEPT_ASSERT((pp::TopologyConcept<CovarianceTopology>));

    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
    typedef pp::default_distance_metric distance_metric_type;
    typedef pp::default_random_sampler random_sampler_type;
    
    typedef ReaK::shared_ptr<state_topology> state_topology_ptr;
    typedef ReaK::shared_ptr<covariance_topology> covariance_topology_ptr;
    
  private:
    state_topology_ptr mean_state_space;
    covariance_topology_ptr covariance_space;
    
  public:  
    
    /**
     * Parametric and default constructor.
     * \param aMeanStateSpace The topology used for the mean-state.
     * \param aCovarianceSpace The topology used for the covariance matrix.
     */
    gaussian_belief_space(const state_topology_ptr& aMeanStateSpace = state_topology_ptr(new state_topology()),
                          const covariance_topology_ptr& aCovarianceSpace = covariance_topology_ptr(new covariance_topology()),
                          const std::string& aName = "") :
                          mean_state_space(aMeanStateSpace),
                          covariance_space(aCovarianceSpace) {
      setName(aName);
    };
    
    /**
     * Computes the distance between two belief-states, using symmetric KL-divergence.
     * \param p1 The first belief-state.
     * \param p2 The second belief-state.
     * \return The symmetric KL-divergence between the two belief-states.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      return double( symKL_divergence(p1,p2,*this) );
    };
    
    /**
     * Computes the norm of a belief-state difference, using symmetric KL-divergence.
     * \param dp The belief-state difference.
     * \return The symmetric KL-divergence between the two end belief-states.
     */
    double norm(const point_difference_type& dp) const {
      return double( symKL_divergence(dp.first, dp.second, *this) );
    };
    
    /**
     * Computes a random belief-state from the underlying state and covariance topologies.
     * \return a random belief-state from the underlying state and covariance topologies.
     */
    point_type random_point() const {
      return point_type( get(pp::random_sampler, *mean_state_space)(*mean_state_space),
                         get(pp::random_sampler, *covariance_space)(*covariance_space));
    };
    
    /**
     * Computes the difference between two belief-states.
     * \param p1 The first belief-state.
     * \param p2 The second belief-state.
     * \return The difference between the two belief-states.
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return point_difference_type(p1, p2);
    };
    
    /**
     * Computes the origin of the belief-space.
     * \return the origin of the belief-space, as a belief-state with the mean-state at the origin of the mean-state topology and the covariance at the origin of the covariance topology.
     */
    point_type origin() const {
      return point_type( mean_state_space->origin(), covariance_space->origin() );
    };
    
    /**
     * Adjusts a belief-state with a belief-state difference.
     * \param p1 The starting belief-state.
     * \param dp The belief-state difference to add to the starting belief-state.
     * \return The adjusted belief-state (semantically p1 + dp).
     */
    point_type adjust(const point_type& p1, const point_difference_type& dp) const {
      return point_type( mean_state_space->adjust( p1.get_mean_state(),
                                                   mean_state_space->difference( dp.first.get_mean_state(),
                                                                                 dp.second.get_mean_state() ) ),
                         covariance_space->adjust( p1.get_covariance(),
                                                   covariance_space->difference( dp.first.get_covariance(),
                                                                                 dp.second.get_covariance() ) ) );
    };
    
    
    /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& p1, double d, const point_type& p2) const {
      return point_type( 
        mean_state_space->move_position_toward( p1.get_mean_state(), d, p2.get_mean_state()),
        covariance_space->move_position_toward( p1.get_covariance(), d, p2.get_covariance()) 
      );
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_back_to(const point_type& p1, double d, const point_type& p2) const {
      return point_type( 
        mean_state_space->move_position_back_to( p1.get_mean_state(), d, p2.get_mean_state()),
        covariance_space->move_position_back_to( p1.get_covariance(), d, p2.get_covariance()) 
      );
    };
    
    
    /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Brings a given point back with the bounds of the space.
     */
    void bring_point_in_bounds(point_type& p1) const {
      mean_state_type m = p1.get_mean_state();
      mean_state_space->bring_point_in_bounds(m);
      p1.set_mean_state(m);
      covariance_type C = p1.get_covariance();
      covariance_space->bring_point_in_bounds(C);
      p1.set_covariance(C);
    };
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    bool is_in_bounds(const point_type& p1) const {
      return (mean_state_space->is_in_bounds(p1.get_mean_state()) &&
              covariance_space->is_in_bounds(p1.get_covariance()));
    };
    
    
    /**
     * Returns the state-topology on which the mean-states lie.
     * \return The state-topology on which the mean-states lie.
     */
    const state_topology& get_state_topology() const {
      return *mean_state_space;
    };

    /**
     * Returns the covariance-topology on which the covariances lie.
     * \return The covariance-topology on which the covariances lie.
     */
    const covariance_topology& get_covariance_topology() const {
      return *covariance_space;
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(mean_state_space)
         & RK_SERIAL_SAVE_WITH_NAME(covariance_space);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(mean_state_space)
         & RK_SERIAL_LOAD_WITH_NAME(covariance_space);
      
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300011,1,"gaussian_belief_space",named_object)

};


/**
 * This class is a bijection mapping from a Gaussian belief-space to a state-space (topology) 
 * by making the assumption of maximum likelihood. In other words, this mapping reduces Gaussian
 * belief-states into their mean value only.
 * 
 * Models: BijectionConcept between a gaussian_belief_space class template and a compatible state topology.
 */
struct gaussian_ML_reduction {
  
  /**
   * This function extracts the mean-value (most likely value) from a Gaussian belief-state.
   * \tparam BeliefPoint The type of the Gaussian belief-state.
   * \tparam StateSpace The original state-space from which a gaussian_belief_space was constructed.
   * \tparam CovarSpace The original covariance-space from which a gaussian_belief_space was constructed.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with (constructable from) the mean-state type that the BeliefPoint type would produce.
   * \param b The belief-state from which the maximum likelihood value is sought.
   * \return The maximum likelihood value of the belief-state.
   */
  template <typename BeliefPoint, typename StateSpace, typename CovarSpace, typename StateSpaceOut>
  typename pp::topology_traits<StateSpaceOut>::point_type map_to_space(const BeliefPoint& b,
                                                                       const gaussian_belief_space<StateSpace, CovarSpace>&,
                                                                       const StateSpaceOut&) const {
    return typename pp::topology_traits<StateSpaceOut>::point_type(b.get_mean_state());
  };
  
};



};


namespace pp {

template <typename StateTopology, typename CovarianceTopology>
struct is_metric_space< ctrl::gaussian_belief_space<StateTopology, CovarianceTopology> > : 
  boost::mpl::and_< is_metric_space<StateTopology>, is_metric_space<CovarianceTopology> > { };

template <typename StateTopology, typename CovarianceTopology>
struct is_reversible_space< ctrl::gaussian_belief_space<StateTopology, CovarianceTopology> > : 
  boost::mpl::and_< is_reversible_space<StateTopology>, is_reversible_space<CovarianceTopology> > { };

template <typename StateTopology, typename CovarianceTopology>
struct is_point_distribution< ctrl::gaussian_belief_space<StateTopology, CovarianceTopology> > : 
  boost::mpl::and_< is_point_distribution<StateTopology>, is_point_distribution<CovarianceTopology> > { };

template <typename StateTopology, typename CovarianceTopology>
struct is_metric_symmetric< ctrl::gaussian_belief_space<StateTopology, CovarianceTopology> > : 
  boost::mpl::and_< is_metric_symmetric<StateTopology>, is_metric_symmetric<CovarianceTopology> > { };

};


};

#endif











