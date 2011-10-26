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

#include "base/defs.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

#include <boost/mpl/bool.hpp>

namespace ReaK {

namespace pp {


/**
 * This traits class defines the traits that characterize a differentiable space. The only 
 * property, currently, is the "order" constant that gives the order of differentiation possible
 * on the given differentiable space (see DifferentiableSpaceConcept).
 * \tparam DifferentiableSpace The differentiable space for which the traits are sought.
 */
template <typename DifferentiableSpace>
struct differentiable_space_traits {
  
  BOOST_STATIC_CONSTANT(std::size_t, order = DifferentiableSpace::differential_order);
};

/**
 * This meta-function provides the type of N-order differential space with a differentiable space bundle.
 * \tparam DifferentiableSpace The differentiable space for which the N-order differential space is sought.
 * \tparam IndependentTopology The independent space against which the differentiation is applied.
 * \tparam Order The order of differentiation of the differential space type that is sought.
 */
template <typename DifferentiableSpace, typename IndependentTopology, unsigned int Order>
struct derived_N_order_space {
  typedef typename DifferentiableSpace::template space<Order,IndependentTopology>::type type;
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
 * space = get_space<0..N>(diff_space,t_space);  The metric space (space) corresponding to the 0 to Nth order derivative space can be obtained given an independent space (e.g. time topology, t_space).
 * 
 * v = lift_to_space<1..N>(dp,dt,diff_space,t_space);  A derivative-point (v) can be obtained from lifting a point-difference (dp) from the space via a difference-point on the independent space (dt). This expression is analogous to v = dp / dt.
 * 
 * dp = descend_to_space<0..N-1>(v,dt,diff_space,t_space);  A point-difference (dp) can be obtained from descending a derivative-point (v) to the space via a difference-point on the independent space (dt). This expression is analogous to dp = v * dt.
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename DifferentiableSpace, unsigned int Order, typename IndependentTopology>
struct DifferentiableSpaceConcept : DifferentiableSpaceConcept<DifferentiableSpace, Order-1, IndependentTopology> {
  
  BOOST_STATIC_ASSERT(differentiable_space_traits<DifferentiableSpace>::order >= Order);
  
  typedef typename derived_N_order_space<DifferentiableSpace,IndependentTopology,Order-1>::type base_space_type;
  typedef typename derived_N_order_space<DifferentiableSpace,IndependentTopology,Order>::type derived_space_type;
  
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept< derived_space_type >));
  
  typename metric_topology_traits<base_space_type>::point_difference_type dp;
  
  typename metric_topology_traits<derived_space_type>::point_type v;
  
  typename metric_topology_traits<IndependentTopology>::point_difference_type dt;
  
  BOOST_CONCEPT_USAGE(DifferentiableSpaceConcept) 
  {
    const derived_space_type& space = get_space<Order>(this->diff_space,this->t_space); RK_UNUSED(space);
    v = lift_to_space<Order>(dp,dt,this->diff_space,this->t_space);
    dp = descend_to_space<Order-1>(v,dt,this->diff_space,this->t_space);
  };
  
};

template <typename DifferentiableSpace, typename IndependentTopology>
struct DifferentiableSpaceConcept<DifferentiableSpace, 0, IndependentTopology> {
  
  typedef typename derived_N_order_space<DifferentiableSpace,IndependentTopology,0>::type base_space_type;
  
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept< base_space_type >));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept< IndependentTopology >));
  
  DifferentiableSpace diff_space;
  IndependentTopology t_space;
  
  BOOST_CONCEPT_USAGE(DifferentiableSpaceConcept) 
  { 
    const base_space_type& space = get_space<0>(this->diff_space,this->t_space); RK_UNUSED(space);
  };
};


};

};

#endif














