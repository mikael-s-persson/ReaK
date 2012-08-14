/**
 * \file invariant_system_concept.hpp
 * 
 * This library defines the traits and concepts related to the definition of an invariant state-space 
 * system. An invariant system is basically one which has an invariant frame, output error and correction 
 * function which allows for the system matrices to be given certain characteristics such as being 
 * symplectic, symmetry-preserving and/or momentum/energy conserving. For certain systems, an invariant
 * formulation can allow for certain specialized algorithms (e.g. the Invariant Kalman Filter) to perform 
 * much better and even be provably optimal for even a non-linear system (given that invariants are properly
 * chosen).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_INVARIANT_SYSTEM_CONCEPT_HPP
#define REAK_INVARIANT_SYSTEM_CONCEPT_HPP

#include "discrete_linear_sss_concept.hpp"
#include "linear_ss_system_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace ctrl {

/**
 * This traits class defines the traits that an invariant state-space system should have (continuous or discrete),
 * in addition to those of a normal state-space system (ss_system_traits and discrete_sss_traits).
 * \tparam InvariantSystem The state-space system type for which the traits are sought.
 */
template <typename InvariantSystem>
struct invariant_system_traits {
  /** This type defines the state-vector type for the system. */
  typedef typename InvariantSystem::point_type point_type;
  
  /** This type describes time. */
  typedef typename InvariantSystem::time_type time_type;
  
  /** This type describes the input vector. */
  typedef typename InvariantSystem::input_type input_type;
  /** This type describes the output vector. */
  typedef typename InvariantSystem::output_type output_type;
  
  /** This type describes the invariant output error vector. */
  typedef typename InvariantSystem::invariant_error_type invariant_error_type;
  /** This type describes the invariant state-correction vector. */
  typedef typename InvariantSystem::invariant_correction_type invariant_correction_type;
  /** This type describes the invariant frame (matrix, usually orthogonal, but not required to be). */
  typedef typename InvariantSystem::invariant_frame_type invariant_frame_type;
  
  /** This constant defines the dimensions of the state vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = InvariantSystem::dimensions);
  /** This constant defines the dimensions of the input vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = InvariantSystem::input_dimensions);
  /** This constant defines the dimensions of the output vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = InvariantSystem::output_dimensions);
  /** This constant defines the dimensions of the output error vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = InvariantSystem::invariant_error_dimensions);
  /** This constant defines the dimensions of the state correction vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = InvariantSystem::invariant_correction_dimensions);
  
};
  
  

/**
 * This concept class template defines the requirements for a discrete-time state-space system to 
 * be considered a valid invariant system.
 * 
 * Required concepts:
 * 
 * the state-space system should model the DiscreteSSSConcept.
 * 
 * the state-space system should be a DiscreteLinearizedSystemType.
 * 
 * Valid expressions:
 * 
 * e = sys.get_invariant_error(p,u,y,t);  The invariant output error vector can be obtained from the current state (p), input (u), output (y) and time (t).
 * 
 * c = transpose(C) * e;  The state-correction vector can be obtained by multiplying the transpose of the state-to-output system matrix (C) with the invariant output error vector.
 * 
 * p = sys.apply_correction(p,c,u,t);  The state-vector (p) can be corrected by a state-correction vector (c), given the input (u) and time (t).
 *
 * W = sys.get_invariant_prior_frame(prev_p,p,u,t);  The invariant frame between a state (prev_p) and its prior prediction (p) can be obtained, given the input (u) and time (t). 
 *
 * W = sys.get_invariant_posterior_frame(prev_p,p,u,t);  The invariant frame between a state (prev_p) and its posterior corrected value (p) can be obtained, given the input (u) and time (t).
 * 
 * \tparam InvariantDiscreteSystem The state-space system to be checked for compliance to this concept.
 * \tparam StateSpaceType The state-space topology type with which the invariant system must be compatible.
 */
template <typename InvariantDiscreteSystem, typename StateSpaceType>
struct InvariantDiscreteSystemConcept : DiscreteLinearSSSConcept<InvariantDiscreteSystem,StateSpaceType,DiscreteLinearizedSystemType> {
  typename invariant_system_traits<InvariantDiscreteSystem>::invariant_error_type e;
  typename invariant_system_traits<InvariantDiscreteSystem>::invariant_correction_type c;
  typename invariant_system_traits<InvariantDiscreteSystem>::invariant_frame_type W;
  
  BOOST_CONCEPT_USAGE(InvariantDiscreteSystemConcept)
  {
    using ReaK::to_vect; using ReaK::from_vect;
    typedef typename invariant_system_traits<InvariantDiscreteSystem>::invariant_correction_type InvCorr;
    typedef typename mat_traits<typename invariant_system_traits<InvariantDiscreteSystem>::invariant_frame_type>::value_type ValueType;
    
    this->e = this->sys.get_invariant_error(this->state_space,this->p,this->u,this->y,this->t);
    this->c = from_vect<InvCorr>(transpose_view(this->C) * to_vect<ValueType>(this->e));
    this->p = this->sys.apply_correction(this->state_space,this->p,this->c,this->u,this->t);
    this->W = this->sys.get_invariant_prior_frame(this->state_space,this->p,this->p,this->u,this->t);
    this->W = this->sys.get_invariant_posterior_frame(this->state_space,this->p,this->p,this->u,this->t);
  };
  
};


/**
 * This concept class template defines the requirements for a continuous-time state-space system to 
 * be considered a valid invariant system.
 * 
 * Required concepts:
 * 
 * the state-space system should model the SSSystemConcept.
 * 
 * the state-space system should be a LinearizedSystemType.
 * 
 * Valid expressions:
 * 
 * e = sys.get_invariant_error(p,u,y,t);  The invariant output error vector can be obtained from the current state (p), input (u), output (y) and time (t).
 * 
 * c = transpose(C) * e;  The state-correction vector can be obtained by multiplying the transpose of the state-to-output system matrix (C) with the invariant output error vector.
 * 
 * p = sys.apply_correction(p,c,u,t);  The state-vector (p) can be corrected by a state-correction vector (c), given the input (u) and time (t).
 * 
 * \tparam InvariantContinuousSystem The state-space system to be checked for compliance to this concept.
 * \tparam StateSpaceType The state-space topology type with which the invariant system must be compatible.
 */
template <typename InvariantContinuousSystem, typename StateSpaceType>
struct InvariantContinuousSystemConcept : LinearSSSystemConcept<InvariantContinuousSystem,StateSpaceType,LinearizedSystemType> {
  typename invariant_system_traits<InvariantContinuousSystem>::invariant_error_type e;
  typename invariant_system_traits<InvariantContinuousSystem>::invariant_correction_type c;
  
  BOOST_CONCEPT_USAGE(InvariantContinuousSystemConcept)
  {
    using ReaK::to_vect; using ReaK::from_vect;
    typedef typename invariant_system_traits<InvariantContinuousSystem>::invariant_correction_type InvCorr;
    typedef typename mat_traits<typename invariant_system_traits<InvariantContinuousSystem>::invariant_frame_type>::value_type ValueType;
    
    this->e     = this->sys.get_invariant_error(this->state_space,this->p,this->u,this->y,this->t);
    this->c     = from_vect<InvCorr>(transpose_view(this->C) * to_vect<ValueType>(this->e));
    this->dp_dt = this->sys.apply_correction(this->state_space,this->p,this->dp_dt,this->c,this->u,this->t);
  };
  
};


};

};

#endif






