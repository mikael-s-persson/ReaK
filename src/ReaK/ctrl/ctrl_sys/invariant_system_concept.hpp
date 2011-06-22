
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

#include "discrete_linear_sss_concept.hpp"
#include "linear_ss_system_concept.hpp"

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
  
  typedef typename InvariantSystem::invariant_error_type invariant_error_type;
  typedef typename InvariantSystem::invariant_correction_type invariant_correction_type;
  
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = InvariantSystem::dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = InvariantSystem::input_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = InvariantSystem::output_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = InvariantSystem::invariant_error_dimensions);
  BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = InvariantSystem::invariant_correction_dimensions);
  
};
  
  
template <typename InvariantDiscreteSystem>
struct InvariantDiscreteSystemConcept {
  InvariantDiscreteSystem sys;
  typename invariant_system_traits<InvariantDiscreteSystem>::point_type p;
  typename invariant_system_traits<InvariantDiscreteSystem>::time_type t;
  typename invariant_system_traits<InvariantDiscreteSystem>::input_type u;
  typename invariant_system_traits<InvariantDiscreteSystem>::output_type y;
  typename invariant_system_traits<InvariantDiscreteSystem>::invariant_error_type e;
  typename invariant_system_traits<InvariantDiscreteSystem>::invariant_correction_type c;
    
  typename discrete_linear_sss_traits<InvariantDiscreteSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantDiscreteSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantDiscreteSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantDiscreteSystem>::matrixD_type D;
  
  void constraints() {
    boost::function_requires< DiscreteSSSConcept<InvariantDiscreteSystem> >();
    
    DiscreteLinearizedSystemType().constraints(sys, p, u, t, A, B, C, D);
    
    e = sys.get_invariant_error(p,u,y,t);
    c = transpose(C) * e;
    p = sys.apply_correction(p,c,u,t);
  };
  
};


  
template <typename InvariantContinuousSystem>
struct InvariantContinuousSystemConcept {
  InvariantContinuousSystem sys;
  typename invariant_system_traits<InvariantContinuousSystem>::point_type p;
  typename ss_system_traits<InvariantContinuousSystem>::point_derivative_type dp_dt;
  typename invariant_system_traits<InvariantContinuousSystem>::time_type t;
  typename invariant_system_traits<InvariantContinuousSystem>::input_type u;
  typename invariant_system_traits<InvariantContinuousSystem>::output_type y;
  typename invariant_system_traits<InvariantContinuousSystem>::invariant_error_type e;
  typename invariant_system_traits<InvariantContinuousSystem>::invariant_correction_type c;
    
  typename linear_ss_system_traits<InvariantContinuousSystem>::matrixA_type A;
  typename linear_ss_system_traits<InvariantContinuousSystem>::matrixB_type B;
  typename linear_ss_system_traits<InvariantContinuousSystem>::matrixC_type C;
  typename linear_ss_system_traits<InvariantContinuousSystem>::matrixD_type D;
  
  void constraints() {
    boost::function_requires< SSSystemConcept<InvariantContinuousSystem> >();
    
    LinearizedSystemType().constraints(sys, p, u, t, A, B, C, D);
    
    e     = sys.get_invariant_error(p,u,y,t);
    c     = transpose(C) * e;
    dp_dt = sys.apply_correction(p,dp_dt,c,u,t);
  };
  
};


};

};

#endif






