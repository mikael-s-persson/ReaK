/**
 * \file bounded_space_concept.hpp
 * 
 * This library defines the concepts that pertain to what can be considered 
 * a bounded metric-space, as used in ReaK::pp. A bounded metric space can be 
 * bounded 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_BOUNDED_SPACE_CONCEPT_HPP
#define REAK_BOUNDED_SPACE_CONCEPT_HPP

#include "metric_space_concept.hpp"

#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  

/**
 * This concept defines the requirements to fulfill in order to model a bounded metric-space 
 * as used in ReaK::pp.
 * 
 * Required Models:
 * 
 * The topology should model the MetricSpaceConcept.
 * 
 * Valid expressions:
 * 
 * pd = space.get_diff_to_boundary(p);  The difference (pd) between a point (p) and the closest boundary-point of the metric-space can be obtained.
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct BoundedSpaceConcept {
  typename metric_topology_traits<Topology>::point_type p;
  typename metric_topology_traits<Topology>::point_difference_type pd;
  Topology space;
  
  BOOST_CONCEPT_USAGE(BoundedSpaceConcept) 
  {
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    pd = space.get_diff_to_boundary(p);
  };
  
};



/**
 * This concept defines the requirements to fulfill in order to model a sphere-bounded metric-space 
 * as used in ReaK::pp.
 * 
 * Required Models:
 * 
 * The topology should model the BoundedSpaceConcept.
 * 
 * Valid expressions:
 * 
 * d = space.get_radius();  The radius (d) can be obtained (maximum distance between a point and the origin of the metric-space).
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct SphereBoundedSpaceConcept {
  Topology space;
  double d;
  
  BOOST_CONCEPT_USAGE(SphereBoundedSpaceConcept) 
  {
    BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<Topology>));
    d = space.get_radius();
  };
  
};




};

};


#endif



