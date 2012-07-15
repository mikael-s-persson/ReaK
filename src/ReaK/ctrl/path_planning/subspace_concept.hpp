/**
 * \file subspace_concept.hpp
 * 
 * This library defines the traits and concepts that pertain to what can be considered 
 * a sub-space, as used in ReaK::pp. Sub-spaces are topologies which are embedded in a 
 * larger topology (or super-space). A typical example is the free configuration space (e.g. C-free)
 * embedded in the overall configuration space. It is useful to be able to relate a sub-space
 * to its super-space, for example, if a motion planner plans a trajectory through C-free, it 
 * is useful to cast that trajectory onto the configuration space such that collision detections
 * are avoided when executing the trajectory.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_SUBSPACE_CONCEPT_HPP
#define REAK_SUBSPACE_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {
  
  
  
/**
 * This traits class defines the types and constants associated to a sub-space.
 * \tparam Topology The topology type for which the topology traits are sought.
 */
template <typename Topology>
struct subspace_traits {
  /** The type of the super-space in which this sub-space is embedded. */
  typedef typename Topology::super_space_type super_space_type;
  
};


/**
 * This concept defines the requirements to fulfill in order to model a sub-space 
 * as used in ReaK::pp.
 * 
 * Valid expressions:
 * 
 * super_space = space.get_super_space();  The super-space can be obtained (at least, by const-reference).
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct SubSpaceConcept {
  Topology space;
  
  BOOST_CONCEPT_USAGE(SubSpaceConcept) 
  {
    const typename subspace_traits<Topology>::super_space_type& super_space = space.get_super_space();
  };
  
};



};

};


#endif


