/**
 * \file airship_mixed_models.hpp
 * 
 * This library contains a number of building blocks for invariantized discrete-time state-space systems 
 * to describe the dynamics of an airship (UAV). These are simplified models, 
 * with forces of limited complexity applied, just free-floating dynamics with 6 dof actuation forces 
 * and some simple augmented states for drag and imbalances. These systems
 * benefit from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum 
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_AIRSHIP_MIXED_MODELS_HPP
#define REAK_AIRSHIP_MIXED_MODELS_HPP

#include "base/named_object.hpp"

#include "ctrl_sys/invariant_system_concept.hpp"
#include "ctrl_sys/augmented_sss_concept.hpp"

#include "topologies/se3_topologies.hpp"
#include "topologies/hyperball_topology.hpp"
#include "topologies/line_topology.hpp"
#include "topologies/temporal_space.hpp"
#include "topologies/time_poisson_topology.hpp"

#include "lin_alg/mat_alg.hpp"
#include "ctrl_sys/gaussian_belief_space.hpp"
#include "ctrl_sys/covariance_matrix.hpp"
#include "ctrl_sys/covar_topology.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/static_assert.hpp>


namespace ReaK {

namespace ctrl {



template <typename InputTuple>
class ss_system_input_tuple : public named_object {
  public:
    typedef ss_system_input_tuple<InputTuple> self;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    typedef vect_n<double> input_type;
    typedef covariance_matrix< vect_n<double> > covar_type;
    typedef gaussian_belief_state< input_type,  covar_type > input_belief_type;
    
    
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 0);
    
  private:
    InputTuple data;
    std::size_t total_input_dim;
    
    BOOST_STATIC_CONSTANT(unsigned int, input_tuple_size = arithmetic_tuple_size<InputTuple>::value);
    
    
    template <unsigned int I>
    typename boost::enable_if_c< (I == 0),
    void >::type construct_input_dimensions_impl(std::size_t& dim) {
      using ReaK::get;
      get<0>(data).construct_input_dimensions(dim);
    };
    
    template <unsigned int I>
    typename boost::enable_if_c< (I != 0),
    void >::type construct_input_dimensions_impl(std::size_t& dim) const {
      using ReaK::get;
      construct_input_dimensions_impl<I-1>(dim);
      get<I>(data).construct_input_dimensions(dim);
    };
    
    void construct_input_dimensions() {
      total_input_dim = 0;
      construct_input_dimensions_impl< input_tuple_size - 1 >(total_input_dim);
    };
    
    
    template <unsigned int I, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I == 0),
    void >::type add_state_difference_impl(const FlyWeight& params,
                                           const StateSpaceType& space, 
                                           const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                           typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                                           const input_type& u, time_difference_type dt, time_type t) const {
      using ReaK::get;
      get<0>(data).add_state_difference(params, space, x, dx, u, dt, t);
    };
    
    template <unsigned int I, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I != 0),
    void >::type add_state_difference_impl(const FlyWeight& params,
                                           const StateSpaceType& space, 
                                           const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                           typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                                           const input_type& u, time_difference_type dt, time_type t) const {
      using ReaK::get;
      add_state_difference_impl<I-1>(params, space, x, dx, u, dt, t);
      get<I>(data).add_state_difference(params, space, x, dx, u, dt, t);
    };
    
    
    template <unsigned int I, typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I == 0),
    void >::type add_state_transition_blocks_impl(MatrixA& A, MatrixB& B, 
                                                  const FlyWeight& params, const StateSpaceType& space, 
                                                  time_type t_0, time_type t_1,
                                                  const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                                  const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                                  const input_type& u_0, const input_type& u_1) const {
      using ReaK::get;
      get<0>(data).add_state_transition_blocks(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
    template <unsigned int I, typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I != 0),
    void >::type add_state_transition_blocks_impl(MatrixA& A, MatrixB& B, 
                                                  const FlyWeight& params, const StateSpaceType& space, 
                                                  time_type t_0, time_type t_1,
                                                  const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                                  const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                                  const input_type& u_0, const input_type& u_1) const {
      using ReaK::get;
      add_state_transition_blocks_impl<I-1>(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
      get<I>(data).add_state_transition_blocks(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
  public:
    
    /**
     * Returns a belief point for zero input values and a given uniform covariance value.
     * \param aCovValue A uniform covariance value to give to all the input component beliefs.
     * \return A belief point for zero input values and a given uniform covariance value.
     */
    input_belief_type get_zero_input_belief(double aCovValue = 1.0) const {
      return input_belief_type(input_type(vect_n<double>(total_input_dim, 0.0)), 
                               covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(total_input_dim,aCovValue))));
    };
    
    /**
     * Returns the dimensions of the input of the system.
     * \return The dimensions of the input of the system.
     */
    std::size_t get_input_dimensions() const { return total_input_dim; };
    
    
    ss_system_input_tuple(const InputTuple& aData = InputTuple()) : data(aData) {
      construct_input_dimensions();
    };
    
#if (!defined(BOOST_NO_RVALUE_REFERENCES) && !defined(BOOST_NO_VARIADIC_TEMPLATES))
    
    template <typename... Args>
    ss_system_input_tuple(Args&&... args) : data(std::forward<Args>(args)...) {
      construct_input_dimensions();
    };
    
#endif
    
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param x The current state object for the system.
     * \param dx The state-difference accumulated (as output) for the different effects.
     * \param u The current input vector for the system.
     * \param dt The time period over which to compute the effect of the input vector.
     * \param t The current time corresponding to the current state of the system.
     */
    template <typename FlyWeight, typename StateSpaceType>
    void add_state_difference(const FlyWeight& params,
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const input_type& u, time_difference_type dt, 
                              time_type t = 0.0) const {
      add_state_difference_impl< input_tuple_size - 1 >(params, space, x, dx, u, dt, t);
    };
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
     * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param t_0 The previous time corresponding to the previous state of the system.
     * \param t_1 The next time corresponding to the next state of the system.
     * \param p_0 The previous state object for the system.
     * \param p_1 The next state object for the system.
     * \param u_0 The previous input vector for the system.
     * \param u_1 The next input vector for the system.
     */
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params,
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const input_type& u_0, const input_type& u_1) const {
      add_state_transition_blocks_impl< input_tuple_size - 1 >(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
    
    
};




template <typename OutputTuple>
class ss_system_output_tuple : public named_object {
  public:
    typedef ss_system_output_tuple<OutputTuple> self;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    typedef vect_n<double> output_type;
    typedef covariance_matrix< vect_n<double> > covar_type;
    typedef gaussian_belief_state< output_type,  covar_type > output_belief_type;
    
    
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 0);
    
  private:
    OutputTuple data;
    std::size_t total_output_dim;
    
    BOOST_STATIC_CONSTANT(unsigned int, output_tuple_size = arithmetic_tuple_size<OutputTuple>::value);
    
    
    template <unsigned int I>
    typename boost::enable_if_c< (I == 0),
    void >::type construct_output_dimensions_impl(std::size_t& dim) {
      using ReaK::get;
      get<0>(data).construct_output_dimensions(dim);
    };
    
    template <unsigned int I>
    typename boost::enable_if_c< (I != 0),
    void >::type construct_output_dimensions_impl(std::size_t& dim) const {
      using ReaK::get;
      construct_output_dimensions_impl<I-1>(dim);
      get<I>(data).construct_output_dimensions(dim);
    };
    
    void construct_output_dimensions() {
      total_output_dim = 0;
      construct_output_dimensions_impl< output_tuple_size - 1 >(total_output_dim);
    };
    
    
    template <unsigned int I, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I == 0),
    void >::type set_output_from_state_impl(const FlyWeight& params,
                                            const StateSpaceType& space, 
                                            const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                            output_type& y, time_type t) const {
      using ReaK::get;
      get<0>(data).set_output_from_state(params, space, x, y, t);
    };
    
    template <unsigned int I, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I != 0),
    void >::type set_output_from_state_impl(const FlyWeight& params,
                                           const StateSpaceType& space, 
                                           const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                           typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                                           const input_type& u, time_difference_type dt, time_type t) const {
      using ReaK::get;
      set_output_from_state_impl<I-1>(params, space, x, y, t);
      get<I>(data).set_output_from_state(params, space, x, y, t);
    };
    
    
    template <unsigned int I, typename MatrixC, typename MatrixD, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I == 0),
    void >::type add_output_function_blocks_impl(MatrixC& C, MatrixD& D, 
                                                 const FlyWeight& params, const StateSpaceType& space,
                                                 time_type t,
                                                 const typename pp::topology_traits<StateSpaceType>::point_type& p, 
                                                 const input_type& u) const {
      using ReaK::get;
      get<0>(data).add_output_function_blocks(A, B, params, space, t, p, u);
    };
    
    template <unsigned int I, typename MatrixC, typename MatrixD, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if_c< (I != 0),
    void >::type add_output_function_blocks_impl(MatrixC& C, MatrixD& D, 
                                                 const FlyWeight& params, const StateSpaceType& space, 
                                                 time_type t,
                                                 const typename pp::topology_traits<StateSpaceType>::point_type& p, 
                                                 const input_type& u) const {
      using ReaK::get;
      add_output_function_blocks_impl<I-1>(C, D, params, space, t, p, u);
      get<I>(data).add_output_function_blocks(C, D, params, space, t, p, u);
    };
    
  public:
    
    /**
     * Returns a belief point for zero input values and a given uniform covariance value.
     * \param aCovValue A uniform covariance value to give to all the input component beliefs.
     * \return A belief point for zero input values and a given uniform covariance value.
     */
    output_belief_type get_zero_output_belief(double aCovValue = 1.0) const {
      return output_belief_type(output_type(vect_n<double>(total_output_dim, 0.0)), 
                                covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(total_output_dim,aCovValue))));
    };
    
    /**
     * Returns the dimensions of the input of the system.
     * \return The dimensions of the input of the system.
     */
    std::size_t get_output_dimensions() const { return total_output_dim; };
    
    
    ss_system_output_tuple(const OutputTuple& aData = OutputTuple()) : data(aData) {
      construct_output_dimensions();
    };
    
#if (!defined(BOOST_NO_RVALUE_REFERENCES) && !defined(BOOST_NO_VARIADIC_TEMPLATES))
    
    template <typename... Args>
    ss_system_output_tuple(Args&&... args) : data(std::forward<Args>(args)...) {
      construct_output_dimensions();
    };
    
#endif
    
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param x The current state object for the system.
     * \param y The current output vector for the system (accumulated as output).
     * \param t The current time corresponding to the current state of the system.
     */
    template <typename FlyWeight, typename StateSpaceType>
    void set_output_from_state(const FlyWeight& params,
                               const StateSpaceType& space, 
                               const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                               output_type& y, time_type t) const {
      set_output_from_state_impl< output_tuple_size - 1 >(params, space, x, y, t);
    };
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param C Holds, as output, the state-to-output jacobian matrix of the output-function of the system.
     * \param D Holds, as output, the input-to-input jacobian matrix of the output-function of the system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param t The current time corresponding to the current state of the system.
     * \param p The current state object for the system.
     * \param u The current input vector for the system.
     */
    template <typename MatrixC, typename MatrixD, typename FlyWeight, typename StateSpaceType>
    void add_output_function_blocks(MatrixC& C, MatrixD& D, 
                                    const FlyWeight& params, const StateSpaceType& space, 
                                    time_type t,
                                    const typename pp::topology_traits<StateSpaceType>::point_type& p, 
                                    const input_type& u) const {
      add_output_function_blocks_impl< output_tuple_size - 1 >(C, D, params, space, t, p, u);
    };
    
    
    
};




template <typename StateSysTuple, typename StateSpaceTuple>
class ss_system_state_tuple : public named_object {
  public:
    typedef ss_system_state_tuple<StateSysTuple, StateSpaceTuple> self;
    
    typedef pp::metric_space_tuple< 
      StateSpaceTuple,
      pp::manhattan_tuple_distance > state_space_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef covariance_matrix< vect_n<double> > covar_type;
    typedef covar_topology< covar_type > covar_space_type;
    typedef pp::temporal_space<state_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_state_space_type;
    typedef gaussian_belief_space<state_space_type, covar_space_type> belief_space_type;
    typedef pp::temporal_space<belief_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_belief_space_type;
    typedef gaussian_belief_state< point_type,  covar_type > state_belief_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, actual_state_dimensions = 0);
    
  private:
    StateSysTuple systems;
    std::size_t total_state_dim;
    std::size_t total_inv_corr_dim;
    std::size_t total_actual_state_dim;
    
    BOOST_STATIC_CONSTANT(unsigned int, system_tuple_size = arithmetic_tuple_size<StateSysTuple>::value);
    
    template <unsigned int I>
    typename boost::enable_if_c< (I == 0),
    void >::type construct_all_dimensions_impl(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) {
      using ReaK::get;
      get<0>(systems).construct_all_dimensions(state_dim, inv_corr_dim, actual_dim);
    };
    
    template <unsigned int I>
    typename boost::enable_if_c< (I != 0),
    void >::type construct_all_dimensions_impl(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) const {
      using ReaK::get;
      construct_all_dimensions_impl<I-1>(state_dim, inv_corr_dim, actual_dim);
      get<I>(systems).construct_all_dimensions(state_dim, inv_corr_dim, actual_dim);
    };
    
    void construct_all_dimensions() {
      total_state_dim = 0;
      total_inv_corr_dim = 0;
      total_actual_state_dim = 0;
      construct_all_dimensions_impl< system_tuple_size - 1 >(total_state_dim, total_inv_corr_dim, total_actual_state_dim);
    };
    
    
    template <unsigned int I>
    typename boost::enable_if_c< (I == 0),
    void >::type get_zero_state_impl(point_type& x) {
      using ReaK::get;
      get<0>(systems).get_zero_state(get<0>(x));
    };
    
    template <unsigned int I>
    typename boost::enable_if_c< (I != 0),
    void >::type get_zero_state_impl(point_type& x) const {
      using ReaK::get;
      get_zero_state_impl<I-1>(x);
      get<I>(systems).get_zero_state(get<I>(x));
    };
    
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    const typename System::point_type& >::type get_state_for_system_impl(const point_type& x) {
      using ReaK::get;
      return get<I>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const typename System::point_type& >::type get_state_for_system_impl(const point_type& x) {
      return get_state_for_system_impl<I-1, System>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const typename System::point_type& >::type get_state_for_system_impl(const point_type&) {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    typename System::point_type& >::type get_state_for_system_impl(point_type& x) {
      using ReaK::get;
      return get<I>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    typename System::point_type& >::type get_state_for_system_impl(point_type& x) {
      return get_state_for_system_impl<I-1, System>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    typename System::point_type& >::type get_state_for_system_impl(point_type&) {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    const typename System::point_difference_type& >::type get_state_diff_for_system_impl(const point_difference_type& x) {
      using ReaK::get;
      return get<I>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const typename System::point_difference_type& >::type get_state_diff_for_system_impl(const point_difference_type& x) {
      return get_state_diff_for_system_impl<I-1, System>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const typename System::point_difference_type& >::type get_state_diff_for_system_impl(const point_difference_type&) {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    typename System::point_difference_type& >::type get_state_diff_for_system_impl(point_difference_type& x) {
      using ReaK::get;
      return get<I>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    typename System::point_difference_type& >::type get_state_diff_for_system_impl(point_difference_type& x) {
      return get_state_diff_for_system_impl<I-1, System>(x);
    };
    
    template <unsigned int I, typename System>
    static typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    typename System::point_difference_type& >::type get_state_diff_for_system_impl(point_difference_type&) {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    const System& >::type get_system_impl() const {
      using ReaK::get;
      return get<I>(systems);
    };
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const System& >::type get_system_impl() const {
      return get_system_impl<I-1, System>();
    };
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    const System& >::type get_system_impl() const {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type >,
    System& >::type get_system_impl() {
      using ReaK::get;
      return get<I>(systems);
    };
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I != 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    System& >::type get_system_impl() {
      return get_system_impl<I-1, System>(x);
    };
    
    template <unsigned int I, typename System>
    typename boost::enable_if<
      boost::mpl::and_< 
        boost::mpl::bool_< I == 0 >, 
        boost::mpl::not_< boost::is_same< System, typename arithmetic_tuple_element<I, StateSysTuple>::type > > >,
    System& >::type get_system_impl() {
      BOOST_STATIC_ASSERT_MSG(false, "The given system type was not found!");
    };
    
    
    
    template <unsigned int I, typename FlyWeight, typename InputType>
    typename boost::enable_if_c< (I == 0),
    void >::type add_state_difference_impl(const FlyWeight& params,
                                           const state_space_type& space, 
                                           const point_type& x, 
                                           point_difference_type& dx,
                                           const InputType& u, time_difference_type dt, time_type t) const {
      using ReaK::get;
      get<0>(systems).add_state_difference(params, space, x, dx, u, dt, t);
    };
    
    template <unsigned int I, typename FlyWeight, typename InputType>
    typename boost::enable_if_c< (I != 0),
    void >::type add_state_difference_impl(const FlyWeight& params,
                                           const state_space_type& space, 
                                           const point_type& x, 
                                           point_difference_type& dx,
                                           const InputType& u, time_difference_type dt, time_type t) const {
      using ReaK::get;
      add_state_difference_impl<I-1>(params, space, x, dx, u, dt, t);
      get<I>(systems).add_state_difference(params, space, x, dx, u, dt, t);
    };
    
    
    template <unsigned int I, typename MatrixA, typename MatrixB, typename FlyWeight, typename InputType>
    typename boost::enable_if_c< (I == 0),
    void >::type add_state_transition_blocks_impl(MatrixA& A, MatrixB& B, 
                                                  const FlyWeight& params, 
                                                  const state_space_type& space, 
                                                  time_type t_0, time_type t_1,
                                                  const point_type& p_0,
                                                  const point_type& p_1, 
                                                  const InputType& u_0, const InputType& u_1) const {
      using ReaK::get;
      get<0>(systems).add_state_transition_blocks(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
    template <unsigned int I, typename MatrixA, typename MatrixB, typename FlyWeight, typename InputType>
    typename boost::enable_if_c< (I != 0),
    void >::type add_state_transition_blocks_impl(MatrixA& A, MatrixB& B, 
                                                  const FlyWeight& params, 
                                                  const state_space_type& space, 
                                                  time_type t_0, time_type t_1,
                                                  const point_type& p_0,
                                                  const point_type& p_1, 
                                                  const InputType& u_0, const InputType& u_1) const {
      using ReaK::get;
      add_state_transition_blocks_impl<I-1>(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
      get<I>(systems).add_state_transition_blocks(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
    
  public:
    
    template <typename System>
    static const typename System::point_type& get_state_for_system(const point_type& x) const {
      return get_state_for_system_impl< system_tuple_size - 1, System >(x);
    };
    
    template <typename System>
    static typename System::point_type& get_state_for_system(point_type& x) const {
      return get_state_for_system_impl< system_tuple_size - 1, System >(x);
    };
    
    
    template <typename System>
    static const typename System::point_difference_type& get_state_diff_for_system(const point_difference_type& x) const {
      return get_state_diff_for_system_impl< system_tuple_size - 1, System >(x);
    };
    
    template <typename System>
    static typename System::point_difference_type& get_state_diff_for_system(point_difference_type& x) const {
      return get_state_diff_for_system_impl< system_tuple_size - 1, System >(x);
    };
    
    
    template <typename System>
    const System& get_system() const {
      return get_system_impl< system_tuple_size - 1, System >(x);
    };
    
    template <typename System>
    System& get_system() {
      return get_system_impl< system_tuple_size - 1, System >(x);
    };
    
    
    /**
     * Returns a belief point for zero input values and a given uniform covariance value.
     * \param aCovValue A uniform covariance value to give to all the input component beliefs.
     * \return A belief point for zero input values and a given uniform covariance value.
     */
    state_belief_type get_zero_state_belief(double aCovValue = 1.0) const {
      point_type x_init;
      get_zero_state_impl< system_tuple_size - 1 >(x_init);
      return state_belief_type(x_init, 
                               covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(total_inv_corr_dim, aCovValue))));
    };
    
    /**
     * Returns the dimensions of the state of the system.
     * \return The dimensions of the state of the system.
     */
    std::size_t get_state_dimensions() const { return total_state_dim; };
    
    /**
     * Returns the dimensions of the invariant correction of the system.
     * \return The dimensions of the invariant correction of the system.
     */
    std::size_t get_correction_dimensions() const { return total_inv_corr_dim; };
    
    /**
     * Returns the dimensions of the actual state of the system.
     * \return The dimensions of the actual state of the system.
     */
    std::size_t get_actual_state_dimensions() const { return total_actual_state_dim; };
    
    
    ss_system_state_tuple(const StateSysTuple& aSystems = StateSysTuple()) : systems(aSystems) {
      construct_all_dimensions();
    };
    
#if (!defined(BOOST_NO_RVALUE_REFERENCES) && !defined(BOOST_NO_VARIADIC_TEMPLATES))
    
    template <typename... Args>
    ss_system_state_tuple(Args&&... args) : systems(std::forward<Args>(args)...) {
      construct_all_dimensions();
    };
    
#endif
    
    
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param x The current state object for the system.
     * \param dx The state-difference accumulated (as output) for the different effects.
     * \param u The current input vector for the system.
     * \param dt The time period over which to compute the effect of the input vector.
     * \param t The current time corresponding to the current state of the system.
     */
    template <typename FlyWeight, typename StateSpaceType>
    void add_state_difference(const FlyWeight& params,
                              const state_space_type& space, 
                              const point_type& x, 
                              point_difference_type& dx,
                              const input_type& u, time_difference_type dt, 
                              time_type t = 0.0) const {
      add_state_difference_impl< system_tuple_size - 1 >(params, space, x, dx, u, dt, t);
    };
    
    /**
     * Fills in the state-difference object with the effects of the given input on the state 
     * over the given time difference.
     * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
     * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
     * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
     * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
     * \param params The fly-weight parameters that describe the dynamics of the system.
     * \param space The space object in which to operate.
     * \param t_0 The previous time corresponding to the previous state of the system.
     * \param t_1 The next time corresponding to the next state of the system.
     * \param p_0 The previous state object for the system.
     * \param p_1 The next state object for the system.
     * \param u_0 The previous input vector for the system.
     * \param u_1 The next input vector for the system.
     */
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params,
                                     const state_space_type& space, 
                                     time_type t_0, time_type t_1,
                                     const point_type& p_0,
                                     const point_type& p_1, 
                                     const input_type& u_0, const input_type& u_1) const {
      add_state_transition_blocks_impl< system_tuple_size - 1 >(A, B, params, space, t_0, t_1, p_0, p_1, u_0, u_1);
    };
    
    
    
};




namespace detail {
  
  typedef arithmetic_tuple< 
            arithmetic_tuple< vect<double,3>, vect<double,3> >, 
            arithmetic_tuple< unit_quat<double>, vect<double,3> > > sat3D_state_type;
  
  inline const sat3D_state_type& get_sat3D_state(const sat3D_state_type& x) { return x; };
  
  template <typename StateTuple>
  const sat3D_state_type& get_sat3D_state(const StateTuple& x) { using ReaK::get; return get<0>(x); };
  
  
  inline sat3D_state_type& get_sat3D_state(sat3D_state_type& x) { return x; };
  
  template <typename StateTuple>
  sat3D_state_type& get_sat3D_state(StateTuple& x) { using ReaK::get; return get<0>(x); };
  
  
  typedef arithmetic_tuple< 
            arithmetic_tuple< vect<double,3>, vect<double,3> >, 
            arithmetic_tuple< vect<double,3>, vect<double,3> > > sat3D_state_diff_type;
  
  inline const sat3D_state_diff_type& get_sat3D_state_diff(const sat3D_state_diff_type& x) { return x; };
  
  template <typename StateTuple>
  const sat3D_state_diff_type& get_sat3D_state_diff(const StateTuple& x) { using ReaK::get; return get<0>(x); };
  
  
  inline sat3D_state_diff_type& get_sat3D_state_diff(sat3D_state_diff_type& x) { return x; };
  
  template <typename StateTuple>
  sat3D_state_diff_type& get_sat3D_state_diff(StateTuple& x) { using ReaK::get; return get<0>(x); };
  
};



class airship3D_6dof_thrusters : public named_object {
  public:
    
    
  private:
    std::size_t start_index;
    
  public:
    
    airship3D_6dof_thrusters() : start_index(0) { };
    
    void construct_input_dimensions(std::size_t& cur_dim) {
      start_index = cur_dim;
      cur_dim += 6;
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const input_type& u, time_difference_type dt, time_type t) const {
      using ReaK::get; // for ADL
      
      vect<double,3> f(u[start_index], u[start_index+1], u[start_index+2]);
      vect<double,3> tau(u[start_index+3], u[start_index+4], u[start_index+5]);
      
      detail::sat3D_state_diff_type& ds = detail::get_sat3D_state_diff(dx);
      
      // velocity:
      f = get_quaternion(detail::get_sat3D_state(x)).as_rotation() * f;
      f *= (dt / params.effective_mass);
      get<1>(get<0>(ds)) += f;
      // position:
      get<0>(get<0>(ds)) += (0.5 * dt) * f;
      
      // ang-velocity
      tau = params.effective_J_inv * tau;
      tau *= dt;
      get<1>(get<1>(ds)) += tau;
      // quaternion-diff (Lie alg.):
      get<0>(get<1>(ds)) += (0.5 * dt) * tau;
      
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const input_type& u_0, const input_type& u_1) const {
      double dt = t_1 - t_0;
      
      mat<double,mat_structure::square> R_0 = get_quaternion(detail::get_sat3D_state(p_0)).as_rotation().getMat();
      
#ifdef USE_HOT_DEL_Q_TERMS
      // TODO Add the d(p + f)/dq and d(v + f)/dq
#endif
      
      // (p,v)-f block:
      sub(B)(range(0,2), range(start_index, start_index+2)) += (0.5 * dt * dt / params.effective_mass) * R_0;
      sub(B)(range(3,5), range(start_index, start_index+2)) += (dt / params.effective_mass) * R_0;
      
      // (q,w)-t block:
      sub(B)(range(6,8), range(start_index+3, start_index+5))  += (0.5 * dt * dt) * params.effective_J_inv;
      sub(B)(range(9,11), range(start_index+3, start_index+5)) += dt * params.effective_J_inv;
    };
    
    
    
};



class tryphon_8_thrusters : public named_object {
  public:
    
    
  private:
    std::size_t start_index;
    
  public:
    
    double side_length;
    
    tryphon_8_thrusters(double aSideLength = 2.0) : start_index(0), side_length(aSideLength) { };
    
    void construct_input_dimensions(std::size_t& cur_dim) {
      start_index = cur_dim;
      cur_dim += 8;
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const input_type& u, time_difference_type dt, time_type t) const {
      vect<double,3> f(
        u[start_index], 
        u[start_index+1], 
        u[start_index+2]);
      vect<double,3> tau(
        u[start_index+3], 
        u[start_index+4], 
        u[start_index+5]);
      
      detail::sat3D_state_diff_type& ds = detail::get_sat3D_state_diff(dx);
      
      // velocity:
      f = get_quaternion(detail::get_sat3D_state(x)).as_rotation() * f;
      f *= (dt / params.effective_mass);
      get<1>(get<0>(ds)) += f;
      // position:
      get<0>(get<0>(ds)) += (0.5 * dt) * f;
      
      // ang-velocity
      tau = params.effective_J_inv * tau;
      tau *= dt;
      get<1>(get<1>(ds)) += tau;
      // quaternion-diff (Lie alg.):
      get<0>(get<1>(ds)) += (0.5 * dt) * tau;
      
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const input_type& u_0, const input_type& u_1) const {
      double dt = t_1 - t_0;
      
      mat<double,mat_structure::square> R_0 = get_quaternion(detail::get_sat3D_state(p_0)).as_rotation().getMat();
      
#ifdef USE_HOT_DEL_Q_TERMS
      // TODO Add the d(p + f)/dq and d(v + f)/dq
#endif
      
      // (p,v)-f block:
      sub(B)(range(0,2), range(start_index, start_index+2)) += (0.5 * dt * dt / params.effective_mass) * R_0;
      sub(B)(range(3,5), range(start_index, start_index+2)) += (dt / params.effective_mass) * R_0;
      
      // (q,w)-t block:
      sub(B)(range(6,8), range(start_index+3, start_index+5))  += (0.5 * dt * dt) * params.effective_J_inv;
      sub(B)(range(9,11), range(start_index+3, start_index+5)) += dt * params.effective_J_inv;
    };
    
    
    
};


class tryphon_12_thrusters : public named_object {
  public:
    
    
  private:
    std::size_t start_index;
    
  public:
    
    double side_length;
    
    tryphon_12_thrusters(double aSideLength = 2.0) : start_index(0), side_length(aSideLength) { };
    
    void construct_input_dimensions(std::size_t& cur_dim) {
      start_index = cur_dim;
      cur_dim += 8;
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const input_type& u, time_difference_type dt, time_type t) const {
      vect<double,3> f(
        u[start_index], 
        u[start_index+1], 
        u[start_index+2]);
      vect<double,3> tau(
        u[start_index+3], 
        u[start_index+4], 
        u[start_index+5]);
      
      detail::sat3D_state_diff_type& ds = detail::get_sat3D_state_diff(dx);
      
      // velocity:
      f = get_quaternion(detail::get_sat3D_state(x)).as_rotation() * f;
      f *= (dt / params.effective_mass);
      get<1>(get<0>(ds)) += f;
      // position:
      get<0>(get<0>(ds)) += (0.5 * dt) * f;
      
      // ang-velocity
      tau = params.effective_J_inv * tau;
      tau *= dt;
      get<1>(get<1>(ds)) += tau;
      // quaternion-diff (Lie alg.):
      get<0>(get<1>(ds)) += (0.5 * dt) * tau;
      
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const input_type& u_0, const input_type& u_1) const {
      double dt = t_1 - t_0;
      
      mat<double,mat_structure::square> R_0 = get_quaternion(detail::get_sat3D_state(p_0)).as_rotation().getMat();
      
#ifdef USE_HOT_DEL_Q_TERMS
      // TODO Add the d(p + f)/dq and d(v + f)/dq
#endif
      
      // (p,v)-f block:
      sub(B)(range(0,2), range(start_index, start_index+2)) += (0.5 * dt * dt / params.effective_mass) * R_0;
      sub(B)(range(3,5), range(start_index, start_index+2)) += (dt / params.effective_mass) * R_0;
      
      // (q,w)-t block:
      sub(B)(range(6,8), range(start_index+3, start_index+5))  += (0.5 * dt * dt) * params.effective_J_inv;
      sub(B)(range(9,11), range(start_index+3, start_index+5)) += dt * params.effective_J_inv;
    };
    
};







/**
 * This class implements an invariantized momentum-tracking discrete-time state-space system describe the 
 * dynamics of a free-floating, near-buoyant, near-balanced 6-dof airship. This is a simplified model, 
 * with just free-floating dynamics with 6 dof actuation forces and some simple augmented states for imbalances. 
 * This system benefits from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum 
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 * This system incorporates augmented states for the mass-imbalance and the eccentricity vector.
 * Also, this system operates within a first-order (once-differentiable) SE(3) topology, augmented by 
 * the parameters.
 */
class airship3D_imdt_em_sys : public named_object {
  public:
    
    typedef pp::metric_space_tuple< 
      arithmetic_tuple< 
        pp::se3_1st_order_topology<double>::type,
        pp::line_segment_topology<double>,
        pp::hyperball_topology< vect<double,3> > >,
      pp::manhattan_tuple_distance > state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
    
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::square> invariant_frame_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 17);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 16);
    BOOST_STATIC_CONSTANT(std::size_t, actual_state_dimensions = 12);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::rectangular> matrixD_type;
    
    struct zero_input_trajectory {
      input_type get_point(time_type) const {
        return input_type(0.0,0.0,0.0,0.0,0.0,0.0);
      };
    };
    
    typedef covariance_matrix< vect_n<double> > covar_type;
    typedef covar_topology< covar_type > covar_space_type;
    typedef pp::temporal_space<state_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_state_space_type;
    typedef gaussian_belief_space<state_space_type, covar_space_type> belief_space_type;
    typedef pp::temporal_space<belief_space_type, pp::time_poisson_topology, pp::time_distance_only> temporal_belief_space_type;
    typedef gaussian_belief_state< point_type,  covar_type > state_belief_type;
    typedef gaussian_belief_state< input_type,  covar_type > input_belief_type;
    typedef gaussian_belief_state< output_type, covar_type > output_belief_type;
    
    virtual shared_ptr< temporal_state_space_type > get_temporal_state_space(double aStartTime = 0.0, double aEndTime = 1.0) const;
    virtual shared_ptr< state_space_type > get_state_space() const;
    
    virtual shared_ptr< temporal_belief_space_type > get_temporal_belief_space(double aStartTime = 0.0, double aEndTime = 1.0) const;
    virtual shared_ptr< belief_space_type > get_belief_space() const;
    
    virtual state_belief_type get_zero_state_belief(double aCovValue = 10.0) const;
    virtual input_belief_type get_zero_input_belief(double aCovValue = 1.0) const;
    virtual output_belief_type get_zero_output_belief(double aCovValue = 1.0) const;
    
  protected:
    double mMass;
    mat<double,mat_structure::symmetric> mInertiaMoment;
    mat<double,mat_structure::symmetric> mInertiaMomentInv;
    time_difference_type mDt;
    vect<double,3> mGravityAcc;
    
  public:  
    
    /**
     * Returns the dimensions of the states of the system.
     * \return The dimensions of the states of the system.
     */
    virtual std::size_t get_state_dimensions() const { return 17; };
    
    /**
     * Returns the dimensions of the input of the system.
     * \return The dimensions of the input of the system.
     */
    virtual std::size_t get_input_dimensions() const { return 6; };
    
    /**
     * Returns the dimensions of the output of the system.
     * \return The dimensions of the output of the system.
     */
    virtual std::size_t get_output_dimensions() const { return 7; };
    
    /**
     * Returns the dimensions of the invariant errors of the system.
     * \return The dimensions of the invariant errors of the system.
     */
    virtual std::size_t get_invariant_error_dimensions() const { return 6; };
    
    /**
     * Returns the dimensions of the corrections to the states of the system.
     * \return The dimensions of the corrections to the states of the system.
     */
    virtual std::size_t get_correction_dimensions() const { return 16; };
    
    /**
     * Returns the dimensions of the actual states of the system.
     * \return The dimensions of the actual states of the system.
     */
    virtual std::size_t get_actual_state_dimensions() const { return 12; };
    
    /**
     * Constructor.
     * \param aName The name for this object.
     * \param aMass The mass of the airship.
     * \param aInertiaMoment The inertia tensor of the airship.
     * \param aDt The time-step for this discrete-time system.
     * \param aGravityAcc The gravitational acceleration vector (usually (0,0,-9.8) for a z-up world).
     */
    airship3D_imdt_em_sys(const std::string& aName = "", 
                          double aMass = 1.0, 
                          const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
                          double aDt = 0.01,
                          const vect<double,3>& aGravityAcc = (vect<double,3>(0.0,0.0,-9.81))); 
    
    /**
     * This function returns the time-step for this discrete-time system.
     * \return The time-step for this discrete-time system.
     */
    time_difference_type get_time_step() const { return mDt; };
    
    /**
     * This function sets the time-step for this discrete-time system.
     * \param aDt The new time-step for this discrete-time system.
     */
    virtual void set_time_step(time_difference_type aDt) { mDt = aDt; };
    
    /**
     * This function computes the next state of the system, i.e., the state at one time-step after the current time.
     * \param space The state-space within which the states reside.
     * \param x The current state of the system.
     * \param u The current input being applied to the system.
     * \param t The current time.
     * \return The state after one time-step beyond the given current state of the system.
     */
    virtual point_type get_next_state(const state_space_type& space, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    /**
     * This function computes the linearization of the state-transitions of the system.
     * In other words, it populates the system matrices with the values appropriate for 
     * the given state-transition.
     * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
     * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
     * \param space The state-space within which the states reside.
     * \param t_0 The time before the state-transition occurred.
     * \param t_1 The time after the state-transition occurred.
     * \param p_0 The state before the state-transition occurred.
     * \param p_1 The state after the state-transition occurred.
     * \param u_0 The input before the state-transition occurred.
     * \param u_1 The input after the state-transition occurred.
     */
    virtual void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, 
                                             const state_space_type& space, 
                                             const time_type& t_0, const time_type& t_1,
                                             const point_type& p_0, const point_type& p_1,
                                             const input_type& u_0, const input_type& u_1) const;
    
    /**
     * This function computes the output of the system corresponding to the current state.
     * \param space The state-space within which the states reside.
     * \param x The current state of the system.
     * \param u The current input being applied to the system.
     * \param t The current time.
     * \return The output for the given current state of the system.
     */
    virtual output_type get_output(const state_space_type& space, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    /**
     * This function computes the linearization of the output-function of the system.
     * In other words, it populates the system matrices with the values appropriate at 
     * the given state.
     * \param C Holds, as output, the state-to-output jacobian matrix of the output-function of the system.
     * \param D Holds, as output, the input-to-output jacobian matrix of the output-function of the system.
     * \param space The state-space within which the states reside.
     * \param t The current time.
     * \param p The current state of the system.
     * \param u The input at the current time.
     */
    virtual void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type& space, 
                                            const time_type& t, const point_type& p, const input_type& u) const;
    
    /**
     * This function computes the invariant output-error of the system corresponding to the current state and the given output.
     * \param space The state-space within which the states reside.
     * \param x The current state of the system.
     * \param u The current input being applied to the system.
     * \param y The output against which to compute the invariant error.
     * \param t The current time.
     * \return The invariant output-error for the given state and output.
     */
    virtual invariant_error_type get_invariant_error(const state_space_type& space, 
                                                     const point_type& x, const input_type& u, 
                                                     const output_type& y, const time_type& t) const;
    
    /**
     * This function computes a state corresponding to the given state corrected by a given invariant term.
     * \param space The state-space within which the states reside.
     * \param x The current state of the system.
     * \param c The invariant correction term to apply to the state.
     * \param u The current input being applied to the system.
     * \param t The current time.
     * \return The corrected state of the system.
     */
    virtual point_type apply_correction(const state_space_type& space, const point_type& x, const invariant_correction_type& c, 
                                        const input_type& u, const time_type& t) const;
    
    /**
     * This function computes the invariant frame transition matrix for the prior stage, 
     * i.e., during a state transition from x_0 to x_1, what invariant frame transition matrix
     * describes the shift from one frame to the other.
     * \param space The state-space within which the states reside.
     * \param x_0 The state of the system before the state-transition.
     * \param x_1 The state of the system after the state-transition.
     * \param u The input being applied to the system before the state-transition.
     * \param t The time before the state-transition.
     * \return The invariant frame transition matrix for the prior stage.
     */
    virtual invariant_frame_type get_invariant_prior_frame(const state_space_type& space, const point_type& x_0, const point_type& x_1, const input_type& u, const time_type& t) const;
    
    /**
     * This function computes the invariant frame transition matrix for the posterior stage, 
     * i.e., during a state correction from x_0 to x_1, what invariant frame transition matrix
     * describes the shift from one frame to the other.
     * \param space The state-space within which the states reside.
     * \param x_0 The state of the system before the correction.
     * \param x_1 The state of the system after the correction.
     * \param u The current input being applied to the system.
     * \param t The current time.
     * \return The invariant frame transition matrix for the posterior stage.
     */
    virtual invariant_frame_type get_invariant_posterior_frame(const state_space_type& space, const point_type& x_0, const point_type& x_1, const input_type& u, const time_type& t) const {
      return get_invariant_prior_frame(space,x_0,x_1,u,t);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_imdt_em_sys,0xC231001A,1,"airship3D_imdt_em_sys",named_object)
    
};

template <>
struct is_invariant_system< airship3D_imdt_em_sys > : boost::mpl::true_ { };

template <>
struct is_augmented_ss_system< airship3D_imdt_em_sys > : boost::mpl::true_ { };



};

};

#endif




