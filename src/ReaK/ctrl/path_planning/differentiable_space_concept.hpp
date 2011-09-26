/**
 * \file differentiable_space_concept.hpp
 * 
 * This library defines the meta-functions, traits and concepts related to the construction of a 
 * differentiable space, i.e. a topology which is a metric-space from which a metric-space 
 * which represents the derivative of that topology can be obtained given a topology against
 * which the derivative is taken (e.g. a time-topology).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_DIFFERENTIABLE_SPACE_CONCEPT_HPP
#define REAK_DIFFERENTIABLE_SPACE_CONCEPT_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

#include <boost/mpl/bool.hpp>

namespace ReaK {

namespace pp {


/**
 * This meta-function computes the metric-space type of the derivative of a metric-space 
 * with respect to another space.
 * \tparam DependentSpace The space type for which the derivative space type is sought.
 * \tparam IndependentSpace The space type against which the derivative is taken.
 */
template <typename DependentSpace, typename IndependentSpace>
struct derivative_space {
  /** The type that describes the derivative space of the DependentSpace against the IndependentSpace. */
  typedef void type;
  /** The bool-type that tells if this derivative space is available. */
  typedef boost::mpl::false_ is_available;
};
  
/**
 * This meta-function computes the metric-space type of the integral of a metric-space 
 * with respect to another space.
 * \tparam DependentSpace The space type for which the integral space type is sought.
 * \tparam IndependentSpace The space type against which the integral is taken.
 */
template <typename DependentSpace, typename IndependentSpace>
struct integral_space {
  /** The type that describes the integral space of the DependentSpace against the IndependentSpace. */
  typedef void type;
  /** The bool-type that tells if this integral space is available. */
  typedef boost::mpl::false_ is_available;
};


/**
 * This concept defines the requirements to fulfill in order to model a differential relation 
 * as used in ReaK::pp. A differentiable relation serves to map a spatial metric-space 
 * its derivative metric-space with respect to a given independent space (e.g. time). The mapping 
 * also serves to map elements between to two spaces.
 * 
 * Required concepts:
 * 
 * All topologies involved should model the MetricSpaceConcept.
 * 
 * Valid expressions:
 * 
 * v_space = diff_relation.get_derivative_space(space, t_space);  The derivative space (v_space) of the metric space (space) can be obtained given an independent space (e.g. time topology, t_space).
 * 
 * space = diff_relation.get_integral_space(v_space, t_space);  The metric space (space) that is an integral of the derivative space (v_space) can be obtained given an independent space (e.g. time topology, t_space).
 * 
 * v = diff_relation.lift_to_derivative(dp,dt);  A derivative-point (v) can be obtained from lifting a point-difference (dp) from the space via a difference-point on the independent space (dt). This expression is analogous to v = dp / dt.
 * 
 * dp = diff_relation.descend_from_derivative(v,dt);  A point-difference (dp) can be obtained from descending a derivative-point (v) to the space via a difference-point on the independent space (dt). This expression is analogous to dp = v * dt.
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename DifferentialRelation, typename DependentTopology, typename IndependentTopology>
struct DifferentialRelationConcept {
  DifferentialRelation diff_relation;
  
  DependentTopology space;
  typename metric_topology_traits<DependentTopology>::point_difference_type dp;
  
  typedef typename derivative_space<DependentTopology,IndependentTopology>::type DerivativeSpace;
  DerivativeSpace v_space;
  typename metric_topology_traits<DerivativeSpace>::point_type v;
  
  IndependentTopology t_space;
  typename metric_topology_traits<IndependentTopology>::point_difference_type dt;
  
  void constraints() {
    boost::function_requires< MetricSpaceConcept< DependentTopology > >();
    boost::function_requires< MetricSpaceConcept< IndependentTopology > >();
    boost::function_requires< MetricSpaceConcept< DerivativeSpace > >();
    const DerivativeSpace& v_space = diff_relation.get_derivative_space(space,t_space);
    const DependentTopology& space = diff_relation.get_derivative_space(v_space,t_space);
    v = diff_relation.lift_to_derivative(dp,dt);
    dp = diff_relation.descend_from_derivative(v,dt);
  };
  
};



};

};

#endif














