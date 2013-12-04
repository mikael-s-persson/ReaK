/**
 * \file discrete_linear_sss_concept.hpp
 * 
 * This library defines a number of traits and concept classes related to the definition 
 * of a discrete-time linear state-space system, with different system types such as a 
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

#ifndef REAK_DISCRETE_LINEAR_SSS_CONCEPT_HPP
#define REAK_DISCRETE_LINEAR_SSS_CONCEPT_HPP


#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "discrete_sss_concept.hpp"
#include "linear_ss_system_concept.hpp"
#include <lin_alg/arithmetic_tuple.hpp>

namespace ReaK {

namespace ctrl {

/**
 * This traits class defines the traits of a discrete-time linear state-space system.
 * This traits class only includes the traits that are not already included in the 
 * discrete_sss_traits class.
 * \tparam DiscreteSystem The discrete-time state-space system whose traits are sought.
 */
template <typename DiscreteSystem>
struct discrete_linear_sss_traits {
  typedef typename DiscreteSystem::matrixA_type matrixA_type;
  typedef typename DiscreteSystem::matrixB_type matrixB_type;
  typedef typename DiscreteSystem::matrixC_type matrixC_type;
  typedef typename DiscreteSystem::matrixD_type matrixD_type;
  
};


/**
 * This concept class defines no requirements for a system to be able to provide 
 * linear(ized) system matrices.
 */
struct DiscreteNonLinearSystemType {
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
 * sys.get_state_transition_blocks(A,B);  The system's state transition matrices can be obtained without providing a time or state.
 * 
 * sys.get_output_function_blocks(C,D);  The system's output function matrices can be obtained without providing a time or state.
 */
struct DiscreteLTISystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point&, const Input&, const Time&,  
                   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space);
    sys.get_output_function_blocks(C, D, state_space);
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
 * sys.get_state_transition_blocks(A,B,t_0,t_1);  The system's state transition matrices can be obtained without providing states, but providing the time before and after the step.
 * 
 * sys.get_output_function_blocks(C,D,t_0);  The system's output function matrices can be obtained without providing a state.
 */
struct DiscreteLTVSystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point&, const Input&, const Time& t, 
                   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space, t, t);
    sys.get_output_function_blocks(C, D, state_space, t);
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
 * sys.get_state_transition_blocks(A,B,t_0,t_1,p_0,p_1,u_0,u_1);  The system's state transition matrices can be obtained by providing the time, state and input, before and after the step is taken.
 * 
 * sys.get_output_function_blocks(C,D,t_0,p_0,u_0);  The system's output function matrices can be obtained from the current time, state and input.
 */
struct DiscreteLinearizedSystemType {
  template <typename System, typename StateSpaceType, typename Point, typename Input, typename Time, 
            typename A_t, typename B_t, typename C_t, typename D_t>
  void constraints(const System& sys, const StateSpaceType& state_space, const Point& p, const Input& u, const Time& t, 
                   A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space, t, t, p, p, u, u);
    sys.get_output_function_blocks(C, D, state_space, t, p, u);
  };
};
  
  
/**
 * This concept class defines the requirements for a discrete-time state-space system to be 
 * considered a linear system, as used in ReaK::ctrl. This concept class depends on a helping
 * concept class which is given as SystemType which defines what kind of linear system it is 
 * (DiscreteLTISystemType, DiscreteLTVSystemType or DiscreteLinearizedSystemType).
 * 
 * Required concepts:
 * 
 * The state-space system should model DiscreteSSSConcept.
 * 
 * Valid expressions:
 * 
 * sys_type.constraints(sys, state_space, p, u, t, A, B, C, D);  The system should comply to the constraints of the SystemType concept class.
 * 
 * p = A * p + B * u;  The next state can be obtained by linear transformation of the state and input using the system matrices.
 * 
 * y = C * p + D * u;  The output can be obtained by linear transformation of the state and input using the system matrices.
 * 
 * \tparam DiscreteSystem The discrete-time state-space system to be tested for linearity.
 * \tparam StateSpaceType The type of the state-space topology on which the state-space system should be able to act.
 * \tparam SystemType The concept class that tests the system-type, can be either DiscreteLTISystemType, DiscreteLTVSystemType or DiscreteLinearizedSystemType.
 */
template <typename DiscreteSystem, typename StateSpaceType, typename SystemType = DiscreteLTISystemType >
struct DiscreteLinearSSSConcept : DiscreteSSSConcept<DiscreteSystem,StateSpaceType> {
  SystemType sys_type;
  
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<DiscreteSystem>::matrixD_type D;
  
  BOOST_CONCEPT_USAGE(DiscreteLinearSSSConcept)
  {
    using ReaK::to_vect;
    using ReaK::from_vect;
    typedef typename discrete_sss_traits<DiscreteSystem>::point_type StateType;
    typedef typename discrete_sss_traits<DiscreteSystem>::output_type OutputType;
    typedef typename mat_traits<typename discrete_linear_sss_traits<DiscreteSystem>::matrixA_type>::value_type ValueType;
    
    sys_type.constraints(this->sys, this->state_space, this->p, this->u, this->t, A, B, C, D);
    this->p = from_vect<StateType>(A * to_vect<ValueType>(this->p) + B * to_vect<ValueType>(this->u));
    this->y = from_vect<OutputType>(C * to_vect<ValueType>(this->p) + D * to_vect<ValueType>(this->u));
  };
  
};




};

};

#endif








