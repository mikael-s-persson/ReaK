/**
 * \file linear_ss_system_concept.hpp
 * 
 * This library defines a number of traits and concept classes related to the definition 
 * of a continuous-time linear state-space system, with different system types such as a 
 * linear-time-invariant (LTI), linear-time-varying (LTV), and linearized. The main 
 * characteristic of such systems is that they can provide system matrices (A,B,C,D), 
 * whether they are independent of state and time (LTI), independent of state (LTV), 
 * or a linear approximation based on the state and time (Linearized).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_LINEAR_SS_SYSTEM_CONCEPT_HPP
#define REAK_LINEAR_SS_SYSTEM_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "state_space_sys_concept.hpp"
#include <lin_alg/arithmetic_tuple.hpp>

namespace ReaK {

namespace ctrl {


/**
 * This traits class defines the traits of a continuous-time linear state-space system.
 * This traits class only includes the traits that are not already included in the 
 * ss_system_traits class.
 * \tparam SSSystem The continuous-time state-space system whose traits are sought.
 */
template <typename SSSystem>
struct linear_ss_system_traits {
  typedef typename SSSystem::matrixA_type matrixA_type;
  typedef typename SSSystem::matrixB_type matrixB_type;
  typedef typename SSSystem::matrixC_type matrixC_type;
  typedef typename SSSystem::matrixD_type matrixD_type;
  
};


/**
 * This concept class defines no requirements for a system to be able to provide 
 * linear(ized) system matrices.
 */
struct NonLinearSystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System&, const StateSpaceType&, const Point&, const Input&, const Time&, 
		   A_t&, B_t&, C_t&, D_t&) { };
};
  
/**
 * This concept class defines the requirement for a system to be able to provide 
 * linear system matrices for an LTI system, that is, the matrices are independent 
 * of time or state or input.
 * 
 * Valid expression:
 * 
 * sys.get_linear_blocks(A,B,C,D);  The system matrices can be obtained without providing a time or state.
 */
struct LTISystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point&, const Input&, const Time&, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D,state_space);
  };
};

/**
 * This concept class defines the requirement for a system to be able to provide 
 * linear system matrices for an LTV system, that is, the matrices are independent 
 * of state or input, but dependent on time. Note that an LTI system is a subset of 
 * an LTV system.
 * 
 * Valid expression:
 * 
 * sys.get_linear_blocks(A,B,C,D,t);  The system matrices can be obtained without providing a state.
 */
struct LTVSystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point&, const Input&, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D,state_space,t);
  };
};

/**
 * This concept class defines the requirement for a system to be able to provide 
 * linear system matrices for a linearized system, that is, the matrices are dependent 
 * on time, state and-or input, but dependent on time. Note that an LTI system and LTV 
 * system are subsets of a linearized system.
 * 
 * Valid expression:
 * 
 * sys.get_linear_blocks(A,B,C,D,t,p,u);  The system matrices can be obtained by providing a time, state and input.
 */
struct LinearizedSystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point& p, const Input& u, const Time& t, 
		   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_linear_blocks(A,B,C,D,state_space,t,p,u);
  };
};

/**
 * This concept class defines the requirements for a continuous-time state-space system to be 
 * considered a linear system, as used in ReaK::ctrl. This concept class depends on a helping
 * concept class which is given as SystemType which defines what kind of linear system it is 
 * (LTISystemType, LTVSystemType or LinearizedSystemType).
 * 
 * Required concepts:
 * 
 * The state-space system should model SSSystemConcept.
 * 
 * Valid expressions:
 * 
 * sys_type.constraints(sys, p, u, t, A, B, C, D);  The system should comply to the constraints of the SystemType concept class.
 * 
 * dp_dt = A * p + B * u;  The state-derivative can be obtained by linear transformation of the state and input using the system matrices.
 * 
 * y = C * p + D * u;  The output can be obtained by linear transformation of the state and input using the system matrices.
 * 
 * \tparam LinearSSSystem The continuous-time state-space system to be tested for linearity.
 * \tparam SystemType The concept class that tests the system-type, can be either LTISystemType, LTVSystemType or LinearizedSystemType.
 */
template <typename LinearSSSystem, typename StateSpaceType, typename SystemType = LTISystemType >
struct LinearSSSystemConcept : SSSystemConcept<LinearSSSystem,StateSpaceType> {
  SystemType sys_type;
  
  typename linear_ss_system_traits<LinearSSSystem>::matrixA_type A;
  typename linear_ss_system_traits<LinearSSSystem>::matrixB_type B;
  typename linear_ss_system_traits<LinearSSSystem>::matrixC_type C;
  typename linear_ss_system_traits<LinearSSSystem>::matrixD_type D;
  
  BOOST_CONCEPT_USAGE(LinearSSSystemConcept)
  {
    using ReaK::to_vect; using ReaK::from_vect;
    typedef typename ss_system_traits<LinearSSSystem>::point_type StateType;
    typedef typename ss_system_traits<LinearSSSystem>::point_derivative_type StateDerivType;
    typedef typename ss_system_traits<LinearSSSystem>::output_type OutputType;
    typedef typename mat_traits<typename linear_ss_system_traits<LinearSSSystem>::matrixA_type>::value_type ValueType;
    
    sys_type.constraints(this->sys, this->state_space, this->p, this->u, this->t, A, B, C, D);
    this->dp_dt = from_vect<StateDerivType>(A * to_vect<ValueType>(this->p) + B * to_vect<ValueType>(this->u));
    this->y     = from_vect<OutputType>(C * to_vect<ValueType>(this->p) + D * to_vect<ValueType>(this->u));
  };
  
};



};

};

#endif




