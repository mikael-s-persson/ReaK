
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

#ifndef DISCRETE_SSS_CONCEPT_HPP
#define DISCRETE_SSS_CONCEPT_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace ctrl {



template <typename DiscreteSystem>
struct discrete_sss_traits {
  typedef typename DiscreteSystem::point_type point_type;
  typedef typename DiscreteSystem::point_difference_type point_difference_type;
  
  typedef typename DiscreteSystem::time_type time_type;
  typedef typename DiscreteSystem::time_difference_type time_difference_type;
  
  typedef typename DiscreteSystem::input_type input_type;
  typedef typename DiscreteSystem::output_type output_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = DiscreteSystem::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = DiscreteSystem::input_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = DiscreteSystem::output_dimensions);
  
};
  
  
template <typename DiscreteSystem>
struct DiscreteSSSConcept {
  DiscreteSystem sys;
  typename discrete_sss_traits<DiscreteSystem>::point_type p;
  typename discrete_sss_traits<DiscreteSystem>::point_difference_type dp;
  typename discrete_sss_traits<DiscreteSystem>::time_type t;
  typename discrete_sss_traits<DiscreteSystem>::time_difference_type dt;
  typename discrete_sss_traits<DiscreteSystem>::input_type u;
  typename discrete_sss_traits<DiscreteSystem>::output_type y;
  void constraints() {
    dp = -dp;
    dp = p - p;
    p = p + dp;
    
    dt = sys.get_time_step();
    p = sys.get_next_state(p,u,t);
    y = sys.get_output(p,u,t);
  };
  
};





};

};

#endif










