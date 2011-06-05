
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

#ifndef LINEAR_SS_SYSTEM_CONCEPT_HPP
#define LINEAR_SS_SYSTEM_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "state_space_sys_concept.hpp"

namespace ReaK {

namespace ctrl {


template <typename SSSystem>
struct linear_ss_system_traits {
  typedef typename SSSystem::matrixA_type matrixA_type;
  typedef typename SSSystem::matrixB_type matrixB_type;
  typedef typename SSSystem::matrixC_type matrixC_type;
  typedef typename SSSystem::matrixD_type matrixD_type;
  
};
  

struct LTISystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point&, const Input&, const Time&, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D);
  };
};

struct LTVSystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point&, const Input&, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D,t);
  };
};

struct LinearizedSystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point& p, const Input& u, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D,t,p,u);
  };
};

  
template <typename LinearSSSystem, typename Topology, typename SystemType = LTISystemType >
struct LinearSSSystemConcept {
  LinearSSSystem sys;
  SystemType sys_type;
  typename ss_system_traits<LinearSSSystem>::point_type p;
  typename ss_system_traits<LinearSSSystem>::point_difference_type dp; 
  typename ss_system_traits<LinearSSSystem>::point_derivative_type dp_dt;
  typename ss_system_traits<LinearSSSystem>::time_type t;
  typename ss_system_traits<LinearSSSystem>::time_difference_type dt;
  typename ss_system_traits<LinearSSSystem>::input_type u;
  typename ss_system_traits<LinearSSSystem>::output_type y;
  
  typename linear_ss_system_traits<LinearSSSystem>::matrixA_type A;
  typename linear_ss_system_traits<LinearSSSystem>::matrixB_type B;
  typename linear_ss_system_traits<LinearSSSystem>::matrixC_type C;
  typename linear_ss_system_traits<LinearSSSystem>::matrixD_type D;
  
  void constraints() {
    boost::function_requires< SSSystemConcept<LinearSSSystem,Topology> >();
    
    sys_type.constraints(sys,p,u,t, A, B, C, D);
    dp_dt = A * p + B * u;
    y     = C * p + D * u;
  };
  
};



};

};

#endif




