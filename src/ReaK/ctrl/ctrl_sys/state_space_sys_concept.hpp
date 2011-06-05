
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

#ifndef STATE_SPACE_SYS_CONCEPT_HPP
#define STATE_SPACE_SYS_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace ctrl {


template <typename SSSystem>
struct ss_system_traits {
  typedef typename SSSystem::point_type point_type;
  typedef typename SSSystem::point_difference_type point_difference_type;
  typedef typename SSSystem::point_derivative_type point_derivative_type;
  
  typedef typename SSSystem::time_type time_type;
  typedef typename SSSystem::time_difference_type time_difference_type;
  
  typedef typename SSSystem::input_type input_type;
  typedef typename SSSystem::output_type output_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = SSSystem::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = SSSystem::input_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = SSSystem::output_dimensions);
  
};
  
  
template <typename SSSystem, typename Topology>
struct SSSystemConcept {
  SSSystem sys;
  typename ss_system_traits<SSSystem>::point_type p;
  typename ss_system_traits<SSSystem>::point_difference_type dp;
  typename ss_system_traits<SSSystem>::point_derivative_type dp_dt;
  typename ss_system_traits<SSSystem>::time_type t;
  typename ss_system_traits<SSSystem>::time_difference_type dt;
  typename ss_system_traits<SSSystem>::input_type u;
  typename ss_system_traits<SSSystem>::output_type y;
  void constraints() {
    dp = -dp;
    p = p + dp;
    
    dp = dp_dt * dt;  //state-space system requirements
    dp_dt = sys.get_state_derivative(p,u,t);
    y = sys.get_output(p,u,t);
  };
  
};




};

};

#endif





