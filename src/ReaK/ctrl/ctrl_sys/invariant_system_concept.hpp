
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

#ifndef INVARIANT_SYSTEM_CONCEPT_HPP
#define INVARIANT_SYSTEM_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace ctrl {


template <typename InvariantSystem>
struct invariant_system_traits {
  typedef typename InvariantSystem::point_type point_type;
  
  typedef typename InvariantSystem::time_type time_type;
  
  typedef typename InvariantSystem::input_type input_type;
  typedef typename InvariantSystem::output_type output_type;
  
  typedef typename InvariantSystem::invariant_type invariant_type;
  typedef typename InvariantSystem::output_error_type output_error_type;
  typedef typename InvariantSystem::invariant_frame_type invariant_frame_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = InvariantSystem::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = InvariantSystem::input_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = InvariantSystem::output_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, invariant_dimensions = InvariantSystem::invariant_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, output_error_dimensions = InvariantSystem::output_error_dimensions);
  
};
  
  
template <typename InvariantSystem, typename Topology>
struct InvariantSystemConcept {
  InvariantSystem sys;
  typename invariant_system_traits<InvariantSystem>::point_type p;
  typename invariant_system_traits<InvariantSystem>::time_type t;
  typename invariant_system_traits<InvariantSystem>::input_type u;
  typename invariant_system_traits<InvariantSystem>::output_type y;
  typename invariant_system_traits<InvariantSystem>::invariant_type i;
  typename invariant_system_traits<InvariantSystem>::output_error_type e;
  typename invariant_system_traits<InvariantSystem>::invariant_frame_type W;
  void constraints() {
    sys.get_invariant(i,t,p,u);
    sys.get_output_error(e,t,p,u,y);
    sys.get_invariant_frame(W,t,p);
  };
  
};




};

};

#endif






