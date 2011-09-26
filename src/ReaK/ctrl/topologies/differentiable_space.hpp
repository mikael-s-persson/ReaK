/**
 * \file differentiable_space.hpp
 * 
 * This library provides classes that define a differentiable-space. A differentiable-space is 
 * a simple association of two topologies to relate them by a derivative-integral relationship.
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

#ifndef REAK_DIFFERENTIABLE_SPACE_HPP
#define REAK_DIFFERENTIABLE_SPACE_HPP

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include <cmath>
#include "time_topology.hpp"
#include "path_planning/metric_space_concept.hpp"

namespace ReaK {

namespace pp {
  

struct default_differentiation_rule {
  template <typename T, typename U, typename V>
  static void lift(T& v, const U& dp, const V& dt) const {
    v = dp / dt;
  };
  template <typename T, typename U, typename V>
  static void descend(T& dp, const U& v, const V& dt) const {
    dp = v * dt;
  };
};



template <typename FirstOrderSpace, typename SecondOrderSpace, typename IndependentSpace = time_topology, typename DifferentiationRule = default_differentiation_rule>
class differentiable_space : public FirstOrderSpace {
  public:
    typedef differentiable_space<FirstOrderSpace,SecondOrderSpace,IndependentSpace> self;
    typedef SecondOrderSpace derivative_space;
    typedef DifferentiationRule diff_rule;
    
    typedef typename metric_topology_traits<derivative_space>::point_type point_derivative_type;
    
    
    
  protected:
    derivative_space tangent_space;
};




};

};

#endif








