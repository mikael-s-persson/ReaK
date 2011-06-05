
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

#ifndef DISCRETE_LINEAR_SSS_CONCEPT_HPP
#define DISCRETE_LINEAR_SSS_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "discrete_sss_concept.hpp"
#include "linear_ss_system_concept.hpp"

namespace ReaK {

namespace ctrl {


template <typename DiscreteSystem>
struct discrete_linear_sss_traits {
  typedef DiscreteSystem::matrixA_type matrixA_type;
  typedef DiscreteSystem::matrixB_type matrixB_type;
  typedef DiscreteSystem::matrixC_type matrixC_type;
  typedef DiscreteSystem::matrixD_type matrixD_type;
  
};

struct DiscreteLTISystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point&, const Input&, const Time&,  
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A, B, C, D);
  };
};

struct DiscreteLTVSystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point&, const Input&, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A, B, C, D, t);
  };
};

struct DiscreteLinearizedSystemType {
  template <typename System, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const Point& p, const Input& u, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A, B, C, D, t, p, u);
  };
};
  
  
template <typename DiscreteSystem, typename SystemType = DiscreteLTISystemType >
struct DiscreteLinearSSSConcept {
  DiscreteSystem sys;
  SystemType sys_type;
  typename discrete_sss_traits<DiscreteSystem>::point_type p;
  typename discrete_sss_traits<DiscreteSystem>::time_type t;
  typename discrete_sss_traits<DiscreteSystem>::input_type u;
  typename discrete_sss_traits<DiscreteSystem>::output_type y;
  
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixD_type D;
  
  void constraints() {
    boost::function_requires< DiscreteSSSConcept<DiscreteSystem> >();
    
    sys_type.constraints(sys, p, u, t, A, B, C, D);
    p = A * p + B * u;
    y = C * p + D * u;
  };
  
};




};

};

#endif








