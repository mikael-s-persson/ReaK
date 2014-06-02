/**
 * \file airship_basic_mixins.hpp
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

#ifndef REAK_AIRSHIP_BASIC_MIXINS_HPP
#define REAK_AIRSHIP_BASIC_MIXINS_HPP

#include <ReaK/core/base/named_object.hpp>

#include "state_space_system_tuple.hpp"

#include <ReaK/ctrl/topologies/se3_topologies.hpp>
#include <ReaK/ctrl/topologies/hyperball_topology.hpp>
#include <ReaK/ctrl/topologies/line_topology.hpp>

namespace ReaK {

namespace ctrl {


class airship_parameter_pack : public named_object {
  public:
    
    double mass;
    double added_mass_factor;
    mat<double, mat_structure::symmetric> J;
    
    double effective_mass;
    double added_mass;
    mat<double, mat_structure::symmetric> effective_J;
    mat<double, mat_structure::symmetric> effective_J_inv;
    
    bool use_hot_del_q_terms;
    bool use_momentum_transfer_terms;
    
    vect<double,3> gravity_acc_vect;
    vect<double,3> magnetic_field_vect;
    
    vect<double,3> IMU_position;
    quaternion<double> IMU_orientation;
    
    double compass_offset;
    
    airship_parameter_pack() : named_object(), mass(1.0), added_mass_factor(0.5), J(1.0, 0.0, 0.0, 1.0, 0.0, 1.0), 
                               effective_mass(mass), effective_J(J), effective_J_inv(J),
                               use_hot_del_q_terms(false), use_momentum_transfer_terms(true),
                               gravity_acc_vect(0.0, 0.0, -9.81), magnetic_field_vect(0.0, 0.0, 0.0),
                               IMU_position(0.0,0.0,0.0), IMU_orientation(), compass_offset(0.0) { setName("airship_parameter_pack"); };
    
    void reset_parameters() {
      effective_mass = mass;
      effective_J = J;
      added_mass = added_mass_factor * mass;
    };
    
    void finalize_parameters() {
      if((effective_J.get_row_count() != 3) || (effective_mass < std::numeric_limits< double >::epsilon()))
        throw system_incoherency("Inertial information is improper in the airship_parameter_pack!");
      
      try {
        invert_Cholesky(effective_J, effective_J_inv);
      } catch(singularity_error&) {
        throw system_incoherency("Inertial tensor is singular in the airship_parameter_pack!");
      };
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mass)
        & RK_SERIAL_SAVE_WITH_NAME(added_mass_factor)
        & RK_SERIAL_SAVE_WITH_NAME(J)
        & RK_SERIAL_SAVE_WITH_NAME(use_hot_del_q_terms)
        & RK_SERIAL_SAVE_WITH_NAME(use_momentum_transfer_terms)
        & RK_SERIAL_SAVE_WITH_NAME(gravity_acc_vect)
        & RK_SERIAL_SAVE_WITH_NAME(magnetic_field_vect)
        & RK_SERIAL_SAVE_WITH_NAME(IMU_position)
        & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
        & RK_SERIAL_SAVE_WITH_NAME(compass_offset);
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mass)
        & RK_SERIAL_LOAD_WITH_NAME(added_mass_factor)
        & RK_SERIAL_LOAD_WITH_NAME(J)
        & RK_SERIAL_LOAD_WITH_NAME(use_hot_del_q_terms)
        & RK_SERIAL_LOAD_WITH_NAME(use_momentum_transfer_terms)
        & RK_SERIAL_LOAD_WITH_NAME(gravity_acc_vect)
        & RK_SERIAL_LOAD_WITH_NAME(magnetic_field_vect)
        & RK_SERIAL_LOAD_WITH_NAME(IMU_position)
        & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
        & RK_SERIAL_LOAD_WITH_NAME(compass_offset);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(airship_parameter_pack,0xC2310020,1,"airship_parameter_pack",named_object)
    
};




class satellite_state_model : public named_object {
  public:
    typedef pp::se3_1st_order_topology<double>::type state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
  private:
    std::size_t state_start_index;
    std::size_t inv_corr_start_index;
    std::size_t actual_state_start_index;
    
  public:
    
    state_space_type create_state_space() const {
      return pp::make_se3_space(
        "satellite_state_space",
        vect<double,3>(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()),
        vect<double,3>( std::numeric_limits<double>::infinity(),  std::numeric_limits<double>::infinity(),  std::numeric_limits<double>::infinity()),
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    };
    
    void get_zero_state(point_type& x) const {
      x = point_type(
        make_arithmetic_tuple(
          vect<double,3>(0.0,0.0,0.0),
          vect<double,3>(0.0,0.0,0.0)
        ),
        make_arithmetic_tuple(
          unit_quat<double>(),
          vect<double,3>(0.0,0.0,0.0)
        )
      );
    };
    
    std::size_t get_state_start_index() const { return state_start_index; };
    std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; };
    std::size_t get_actual_state_start_index() const { return actual_state_start_index; };
    
    satellite_state_model() : named_object() {
      setName("satellite_state_model");
    };
    
    void construct_all_dimensions(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) {
      state_start_index = state_dim;
      state_dim += 13;
      inv_corr_start_index = inv_corr_dim;
      inv_corr_dim += 12;
      actual_state_start_index = actual_dim;
      actual_dim += 12;
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_to_fly_weight_params(const FlyWeight& params, 
                                  const StateSpaceType& space, const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                  const InputType&, time_difference_type dt, time_type t) const { };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const InputType&, time_difference_type dt, time_type t) const {
      
      const point_type& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      point_difference_type& dx_se3 = params.get_state_models().template get_state_diff_for_system<satellite_state_model>(dx);
      
      // position:
      get<0>(get<0>(dx_se3)) += dt * get_velocity(x_se3);
      
      // quaternion-diff (Lie alg.):
      get<0>(get<1>(dx_se3)) += dt * get_ang_velocity(x_se3);
      
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const InputType& u_0, const InputType& u_1) const {
      const point_type& x0_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_0);
      const point_type& x1_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_1);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      double dt = t_1 - t_0;
      
      const std::pair<std::size_t, std::size_t> p_r(inv_corr_start_index, inv_corr_start_index+3);
      const std::pair<std::size_t, std::size_t> v_r(inv_corr_start_index+3, inv_corr_start_index+6);
      const std::pair<std::size_t, std::size_t> q_r(inv_corr_start_index+6, inv_corr_start_index+9);
      const std::pair<std::size_t, std::size_t> w_r(inv_corr_start_index+9, inv_corr_start_index+12);
      
      // Position row:
      // p-p block:
      sub(A)(p_r, p_r) += mat_ident<double>(3);
      // p-v block:
      sub(A)(p_r, v_r) += dt * mat_ident<double>(3);
      // v-v block:
      sub(A)(v_r, v_r) += mat_ident<double>(3);
      
      
      mat<double,mat_structure::square> R_0_1((invert(get_quaternion(x1_se3).as_rotation()) * get_quaternion(x0_se3).as_rotation()).getMat());
      mat<double,mat_structure::square> JRJ(sys_params.effective_J_inv * R_0_1 * sys_params.effective_J);
      
      sub(A)(q_r, q_r) += R_0_1;
      sub(A)(q_r, w_r) += dt * R_0_1;
      sub(A)(w_r, w_r) += JRJ;
      
      if( sys_params.use_hot_del_q_terms ) {
        vect<double,3> l_net_0 = sys_params.effective_J * get_ang_velocity(x0_se3);
        vect<double,3> l_net_1 = sys_params.effective_J * get_ang_velocity(x1_se3);
        
        sub(A)(q_r, q_r) -= (0.5 * dt) * sys_params.effective_J_inv * R_0_1 * mat<double,mat_structure::skew_symmetric>(l_net_0);
        
        sub(A)(q_r, w_r) += (0.5 * dt) * (JRJ - R_0_1);
        
        sub(A)(w_r, q_r) += sys_params.effective_J_inv * (mat<double,mat_structure::skew_symmetric>(l_net_1) * R_0_1 
                                                      - R_0_1 * mat<double,mat_structure::skew_symmetric>(l_net_0));
        
        sub(A)(w_r, w_r) += (0.5 * dt) * sys_params.effective_J_inv * mat<double,mat_structure::skew_symmetric>(l_net_1) * R_0_1;
      };
      
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InvCorrType, typename InputType>
    void apply_correction_to_state(const FlyWeight& params, const StateSpaceType& space, 
                                   const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                   typename pp::topology_traits<StateSpaceType>::point_type& x_c, 
                                   const InvCorrType& c, const InputType& u, const time_type& t) const {
      const point_type& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      point_type& x_c_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x_c);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[inv_corr_start_index+6], c[inv_corr_start_index+7], c[inv_corr_start_index+8]) );
      unit_quat<double> q_new = get_quaternion(x_se3) * q_diff;
      
      vect<double,3> w_new = sys_params.effective_J_inv * (invert(q_diff).as_rotation() * (sys_params.effective_J * get_ang_velocity(x_se3) + vect<double,3>(c[inv_corr_start_index+9],c[inv_corr_start_index+10],c[inv_corr_start_index+11])));
      
      x_c_se3 = point_type(
        make_arithmetic_tuple(
          get_position(x_se3) + vect<double,3>(c[inv_corr_start_index],c[inv_corr_start_index+1],c[inv_corr_start_index+2]),
          get_velocity(x_se3) + vect<double,3>(c[inv_corr_start_index+3],c[inv_corr_start_index+4],c[inv_corr_start_index+5])
        ),
        make_arithmetic_tuple(
          q_new, 
          w_new
        )
      );
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType, typename InvarFrameType>
    void set_invariant_frame_blocks(const FlyWeight& params, const StateSpaceType& space, 
                                    InvarFrameType& invar_frame, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_0, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_1, 
                                    const InputType& u, const time_type& t) const {
      /* identity is OK */
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(satellite_state_model,0xC2310021,1,"satellite_state_model",named_object)
    
};


class near_buoyancy_state_model : public named_object {
  public:
    typedef pp::line_segment_topology<double> state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
  private:
    
    std::size_t state_start_index;
    std::size_t inv_corr_start_index;
    
  public:
    
    state_space_type create_state_space() const {
      return pp::line_segment_topology<double>("mass_imbalance_param_space", 0.0, std::numeric_limits<double>::infinity());
    };
    
    void get_zero_state(point_type& x) const {
      x = 0.0;
    };
    
    std::size_t get_state_start_index() const { return state_start_index; };
    std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; };
    
    near_buoyancy_state_model() { };
    
    void construct_all_dimensions(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) {
      state_start_index = state_dim;
      state_dim += 1;
      inv_corr_start_index = inv_corr_dim;
      inv_corr_dim += 1;
      RK_UNUSED(actual_dim);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_to_fly_weight_params(const FlyWeight& params, 
                                  const StateSpaceType& space, const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                  const InputType&, time_difference_type dt, time_type t) const {
      const point_type dm = params.get_state_models().template get_state_for_system<near_buoyancy_state_model>(x);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      sys_params.effective_mass += dm;
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const InputType&, time_difference_type dt, time_type t) const {
      
      typedef satellite_state_model::point_difference_type SE3StateDiff;
      
      SE3StateDiff& dx_se3 = params.get_state_models().template get_state_diff_for_system<satellite_state_model>(dx);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const point_type dm = params.get_state_models().template get_state_for_system<near_buoyancy_state_model>(x);
      
      vect<double,3> gf = dm * sys_params.gravity_acc_vect;
      
      // velocity:
      gf *= (dt / (sys_params.effective_mass + sys_params.added_mass));
      get<1>(get<0>(dx_se3)) += gf;
      // position:
      get<0>(get<0>(dx_se3)) += (0.5 * dt) * gf;
      
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const InputType& u_0, const InputType& u_1) const {
      
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      
      const double dt = t_1 - t_0;
      
      const std::pair<std::size_t, std::size_t> p_r(sat3d_state_index, sat3d_state_index+3);
      const std::pair<std::size_t, std::size_t> v_r(sat3d_state_index+3, sat3d_state_index+6);
      
      const std::size_t m_r = inv_corr_start_index;
      
      // p-m block:
      slice(A)(p_r, m_r) += (0.5 * dt * dt) * sys_params.gravity_acc_vect;
      // v-m block:
      slice(A)(v_r, m_r) += dt * sys_params.gravity_acc_vect;
      
      A(m_r, m_r) += 1.0;
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InvCorrType, typename InputType>
    void apply_correction_to_state(const FlyWeight& params, const StateSpaceType& space, 
                                   const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                   typename pp::topology_traits<StateSpaceType>::point_type& x_c, 
                                   const InvCorrType& c, const InputType& u, const time_type& t) const {
      const point_type& dm = params.get_state_models().template get_state_for_system<near_buoyancy_state_model>(x);
      point_type& dm_c = params.get_state_models().template get_state_for_system<near_buoyancy_state_model>(x_c);
      
      dm_c = dm + c[inv_corr_start_index];
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType, typename InvarFrameType>
    void set_invariant_frame_blocks(const FlyWeight& params, const StateSpaceType& space, 
                                    InvarFrameType& invar_frame, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_0, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_1, 
                                    const InputType& u, const time_type& t) const {
      /* identity is OK */
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(near_buoyancy_state_model,0xC2310022,1,"near_buoyancy_state_model",named_object)
    
};


class eccentricity_state_model : public named_object {
  public:
    typedef pp::hyperball_topology< vect<double,3> > state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
  private:
    
    std::size_t state_start_index;
    std::size_t inv_corr_start_index;
    
  public:
    
    state_space_type create_state_space() const {
      return pp::hyperball_topology< vect<double,3> >("eccentricity_param_space", vect<double,3>(0.0,0.0,0.0), std::numeric_limits<double>::infinity());
    };
    
    void get_zero_state(point_type& x) const {
      x = vect<double,3>(0.0,0.0,0.0);
    };
    
    std::size_t get_state_start_index() const { return state_start_index; };
    std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; };
    
    eccentricity_state_model() { };
    
    void construct_all_dimensions(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) {
      state_start_index = state_dim;
      state_dim += 3;
      inv_corr_start_index = inv_corr_dim;
      inv_corr_dim += 3;
      RK_UNUSED(actual_dim);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_to_fly_weight_params(const FlyWeight& params, 
                                  const StateSpaceType& space, const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                  const InputType&, time_difference_type dt, time_type t) const {
      const point_type& r = params.get_state_models().template get_state_for_system<eccentricity_state_model>(x);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      mat<double, mat_structure::skew_symmetric> r_cross(r);
      sys_params.effective_J -= mat<double, mat_structure::symmetric>(sys_params.effective_mass * r_cross * r_cross);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const InputType&, time_difference_type dt, time_type t) const {
      
      typedef satellite_state_model::point_type SE3State;
      typedef satellite_state_model::point_difference_type SE3StateDiff;
      
      const SE3State& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      SE3StateDiff& dx_se3 = params.get_state_models().template get_state_diff_for_system<satellite_state_model>(dx);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const vect<double,3>& r = params.get_state_models().template get_state_for_system<eccentricity_state_model>(x);
      mat<double,mat_structure::skew_symmetric> r_cross(r);
      
      quaternion<double> q = get_quaternion(x_se3).as_rotation();
      vect<double,3> gt_impulse = (dt * sys_params.effective_mass) * (sys_params.effective_J_inv * (r % (invert(q) * sys_params.gravity_acc_vect)));
      
      // ang-velocity
      get<1>(get<1>(dx_se3)) += gt_impulse;
      // quat-diff:
      get<0>(get<1>(dx_se3)) += (0.5 * dt) * gt_impulse;
      
      if( sys_params.use_momentum_transfer_terms ) {
        vect<double,3> l_transfer = r % gt_impulse;
        // neglects HOT in (I + m [rx] J^-1 [rx])^-1 (r % gt_impulse) ... by approx I + mrJr ~= I
        //  here is the adjustment for the HOTs:
        try {
          mat<double, mat_structure::symmetric> X(mat_ident<double>(3) + sys_params.effective_mass * (r_cross * sys_params.effective_J_inv * r_cross));
          mat<double, mat_structure::rectangular> b(3,1);
          slice(b)(range(0,3),0) = l_transfer;
          linsolve_Cholesky(X, b, 1e-6);
          l_transfer = slice(b)(range(0,3),0); // <-- commit change if cholesky succeeded.
        } catch(...) { /* if cholesky failed, no HOT adjustment is applied to l_transfer */ };
        
        // The transfer fraction reflects the fact that you can't have the linear momentum of the airflow to transfer as angular momentum of the airship. Or can you???
        double transfer_frac = sys_params.effective_mass / (sys_params.effective_mass + sys_params.added_mass);
        
        vect<double,3> p_transfer = (-sys_params.effective_mass * transfer_frac) * (sys_params.effective_J_inv * (r % l_transfer));
        
        l_transfer = transfer_frac * (q * l_transfer);
        
        get<1>(get<0>(dx_se3)) += l_transfer;
        get<0>(get<0>(dx_se3)) += (0.5 * dt) * l_transfer;
        
        get<1>(get<1>(dx_se3)) += p_transfer;
        get<0>(get<1>(dx_se3)) += (0.5 * dt) * p_transfer;
      };
      
    };
    
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    typename boost::enable_if< 
      typename FlyWeight::state_models_type::template has_system<near_buoyancy_state_model>,
    void >::type add_state_transition_blocks_for_dm(
        MatrixA& A, MatrixB& B,
        const FlyWeight& params, 
        const StateSpaceType& space, 
        double dt,
        const mat<double,mat_structure::square>& R_0,
        const mat<double,mat_structure::square>& R_0_1,
        const vect<double,3>& r,
        const mat<double,mat_structure::skew_symmetric>& r_cross,
        const vect<double,3>& local_g) const {
      
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const std::pair<std::size_t, std::size_t> p_r(sat3d_state_index, sat3d_state_index+3);
      const std::pair<std::size_t, std::size_t> v_r(sat3d_state_index+3, sat3d_state_index+6);
      const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index+6, sat3d_state_index+9);
      const std::pair<std::size_t, std::size_t> w_r(sat3d_state_index+9, sat3d_state_index+12);
      
      const std::size_t m_r = params.get_state_models().template get_system<near_buoyancy_state_model>().get_inv_corr_start_index();
      
      // The transfer fraction reflects the fact that you can't have the linear momentum of the airflow to transfer as angular momentum of the airship. Or can you???
      double transfer_frac = sys_params.effective_mass / (sys_params.effective_mass + sys_params.added_mass);
      
      if( sys_params.use_momentum_transfer_terms ) {
        // v-m :
        vect<double,3> rJrg = r % (sys_params.effective_J_inv * (r % local_g));
        try {
          mat<double, mat_structure::symmetric> X(mat_ident<double>(3) + sys_params.effective_mass * (r_cross * sys_params.effective_J_inv * r_cross));
          mat<double, mat_structure::rectangular> b(3,1);
          slice(b)(range(0,3),0) = rJrg;
          linsolve_Cholesky(X, b, 1e-6);
          rJrg = slice(b)(range(0,3),0); // <-- commit change if cholesky succeeded.
        } catch(...) { /* if cholesky failed, no HOT adjustment is applied to rJrg */ };
        rJrg = R_0 * rJrg;
        rJrg *= transfer_frac;
        slice(A)(v_r, m_r) += rJrg;
        slice(A)(p_r, m_r) += (0.5 * dt) * rJrg;
      };
      
      // w-m :
      vect<double,3> Jrg = sys_params.effective_J_inv * (r % local_g);
      if( sys_params.use_momentum_transfer_terms ) {
        //  here is the adjustment for the HOTs:
        try {
          mat<double, mat_structure::symmetric> Y(mat_ident<double>(3) + (sys_params.effective_mass * transfer_frac) * (sys_params.effective_J_inv * r_cross * r_cross));
          mat<double, mat_structure::rectangular> b(3,1);
          slice(b)(range(0,3),0) = Jrg;
          linsolve_Cholesky(Y, b, 1e-6);
          Jrg = slice(b)(range(0,3),0); // <-- commit change if cholesky succeeded.
        } catch(...) { /* if cholesky failed, no HOT adjustment is applied to delw */ };
      };
      Jrg = R_0_1 * Jrg;
      slice(A)(w_r, m_r) += Jrg;
      slice(A)(q_r, m_r) += (0.5 * dt) * Jrg;
      
    };
    
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType>
    typename boost::disable_if< 
      typename FlyWeight::state_models_type::template has_system<near_buoyancy_state_model>,
    void >::type add_state_transition_blocks_for_dm(
        MatrixA&, MatrixB&, const FlyWeight&, const StateSpaceType&, 
        double, const mat<double,mat_structure::square>&, const mat<double,mat_structure::square>&,
        const vect<double,3>&, const mat<double,mat_structure::skew_symmetric>&, const vect<double,3>&) const { };
    
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const InputType& u_0, const InputType& u_1) const {
      typedef satellite_state_model::point_type SE3State;
      
      const SE3State& x0_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_0);
      const SE3State& x1_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_1);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const vect<double,3>& r = params.get_state_models().template get_state_for_system<eccentricity_state_model>(p_0);
      mat<double,mat_structure::skew_symmetric> r_cross(r);
      
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      
      const double dt = t_1 - t_0;
      
      const std::pair<std::size_t, std::size_t> p_r(sat3d_state_index, sat3d_state_index+3);
      const std::pair<std::size_t, std::size_t> v_r(sat3d_state_index+3, sat3d_state_index+6);
      const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index+6, sat3d_state_index+9);
      const std::pair<std::size_t, std::size_t> w_r(sat3d_state_index+9, sat3d_state_index+12);
      
      const std::pair<std::size_t, std::size_t> r_r(inv_corr_start_index, inv_corr_start_index+3);
      
      mat<double,mat_structure::square> R_0_1((invert(get_quaternion(x1_se3).as_rotation()) * get_quaternion(x0_se3).as_rotation()).getMat());
      
      mat<double,mat_structure::square> R_0(get_quaternion(x0_se3).as_rotation().getMat());
      vect<double,3> local_g = dt * (transpose_view(R_0) * sys_params.gravity_acc_vect);
      mat<double,mat_structure::skew_symmetric> local_g_cross(local_g);
      
      // The transfer fraction reflects the fact that you can't have the linear momentum of the airflow to transfer as angular momentum of the airship. Or can you???
      double transfer_frac = sys_params.effective_mass / (sys_params.effective_mass + sys_params.added_mass);
      
      if( sys_params.use_momentum_transfer_terms ) {
        // p/v-r
        mat<double,mat_structure::skew_symmetric> local_Jrg_cross(sys_params.effective_J_inv * (r % local_g));
        mat<double,mat_structure::square> delv(sys_params.effective_mass * (local_Jrg_cross + r_cross * sys_params.effective_J_inv * local_g_cross));
        // neglects the HOT in : ... + dt m R_0 d/dr((I + m [rx] J^-1 [rx])^-1) [rx] J^-1 [rx] R_0^T g
        //   by assuming that (I + m [rx] J^-1 [rx]) ~= I
        //  here is the adjustment for the HOTs:
        try {
          mat<double, mat_structure::symmetric> X(mat_ident<double>(3) + sys_params.effective_mass * (r_cross * sys_params.effective_J_inv * r_cross));
          mat<double, mat_structure::square> b(delv);
          linsolve_Cholesky(X, b, 1e-6);
          delv = b; // <-- commit change if cholesky succeeded.
        } catch(...) { /* if cholesky failed, no HOT adjustment is applied to delv */ };
        delv = R_0 * delv;
        delv *= transfer_frac;
        sub(A)(v_r, r_r) -= delv;
        sub(A)(p_r, r_r) -= (0.5 * dt) * delv;
          
      };
      
      // q/w-r
      mat<double,mat_structure::square> delw(sys_params.effective_mass * (sys_params.effective_J_inv * local_g_cross));
      // neglects the HOT in : ... + dt m d/dr((I + m J^-1 [rx] [rx])^-1) J^-1 [rx] R_0^T g
      //   by assuming that (I + m J^-1 [rx] [rx]) ~= I
      if( sys_params.use_momentum_transfer_terms ) {
        //  here is the adjustment for the HOTs:
        try {
          mat<double, mat_structure::symmetric> Y(mat_ident<double>(3) + (sys_params.effective_mass * transfer_frac) * (sys_params.effective_J_inv * r_cross * r_cross));
          mat<double, mat_structure::square> b(delw);
          linsolve_Cholesky(Y, b, 1e-6);
          delw = b; // <-- commit change if cholesky succeeded.
        } catch(...) { /* if cholesky failed, no HOT adjustment is applied to delw */ };
      };
      delw = R_0_1 * delw;
      
      sub(A)(w_r, r_r) -= delw;
      sub(A)(q_r, r_r) -= (0.5 * dt) * delw;
      
      add_state_transition_blocks_for_dm(A, B, params, space, dt, R_0, R_0_1, r, r_cross, local_g);
      
      sub(A)(r_r, r_r) += mat_ident<double>(3);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InvCorrType, typename InputType>
    void apply_correction_to_state(const FlyWeight& params, const StateSpaceType& space, 
                                   const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                   typename pp::topology_traits<StateSpaceType>::point_type& x_c, 
                                   const InvCorrType& c, const InputType& u, const time_type& t) const {
      const point_type& r = params.get_state_models().template get_state_for_system<eccentricity_state_model>(x);
      point_type& r_c = params.get_state_models().template get_state_for_system<eccentricity_state_model>(x_c);
      
      r_c = r + vect<double,3>(c[inv_corr_start_index], c[inv_corr_start_index+1], c[inv_corr_start_index+2]);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType, typename InvarFrameType>
    void set_invariant_frame_blocks(const FlyWeight& params, const StateSpaceType& space, 
                                    InvarFrameType& invar_frame, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_0, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_1, 
                                    const InputType& u, const time_type& t) const {
      /* identity is OK */
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(eccentricity_state_model,0xC2310023,1,"eccentricity_state_model",named_object)
    
};


class cross_inertia_state_model : public named_object {
  public:
    typedef pp::hyperball_topology< vect<double,3> > state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
  private:
    
    std::size_t state_start_index;
    std::size_t inv_corr_start_index;
    
  public:
    
    state_space_type create_state_space() const {
      return pp::hyperball_topology< vect<double,3> >("cross_inertia_param_space", vect<double,3>(0.0,0.0,0.0), std::numeric_limits<double>::infinity());
    };
    
    void get_zero_state(point_type& x) const {
      x = vect<double,3>(0.0,0.0,0.0);
    };
    
    std::size_t get_state_start_index() const { return state_start_index; };
    std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; };
    
    cross_inertia_state_model() { };
    
    void construct_all_dimensions(std::size_t& state_dim, std::size_t& inv_corr_dim, std::size_t& actual_dim) {
      state_start_index = state_dim;
      state_dim += 3;
      inv_corr_start_index = inv_corr_dim;
      inv_corr_dim += 3;
      RK_UNUSED(actual_dim);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_to_fly_weight_params(const FlyWeight& params, 
                                  const StateSpaceType& space, const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                  const InputType&, time_difference_type dt, time_type t) const {
      const point_type& s = params.get_state_models().template get_state_for_system<cross_inertia_state_model>(x);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      sys_params.effective_J += mat<double, mat_structure::symmetric>(0.0, s[0], s[1], 0.0, s[2], 0.0);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_difference(const FlyWeight& params, 
                              const StateSpaceType& space, 
                              const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                              typename pp::topology_traits<StateSpaceType>::point_difference_type& dx,
                              const InputType&, time_difference_type dt, time_type t) const {
      /* nothing to do, the main effect is the addition to the effective inertia */
    };
    
    template <typename MatrixA, typename MatrixB, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                     const FlyWeight& params, 
                                     const StateSpaceType& space, 
                                     time_type t_0, time_type t_1,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_0,
                                     const typename pp::topology_traits<StateSpaceType>::point_type& p_1, 
                                     const InputType& u_0, const InputType& u_1) const {
      typedef satellite_state_model::point_type SE3State;
      
      const SE3State& x0_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_0);
      const SE3State& x1_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(p_1);
      typename FlyWeight::system_param_type& sys_params = params.get_system_parameters();
      
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      
      const double dt = t_1 - t_0;
      
      const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index+6, sat3d_state_index+9);
      const std::pair<std::size_t, std::size_t> w_r(sat3d_state_index+9, sat3d_state_index+12);
      
      const std::pair<std::size_t, std::size_t> s_r(inv_corr_start_index, inv_corr_start_index+3);
      
      mat<double,mat_structure::square> R_0_1((invert(get_quaternion(x1_se3).as_rotation()) * get_quaternion(x0_se3).as_rotation()).getMat());
      
      // w-sigma block:
      const vect<double,3>& w_0 = get_ang_velocity(x0_se3);
      const vect<double,3>& w_1 = get_ang_velocity(x1_se3);
      mat<double,mat_structure::square> del_sig_0(w_0[1], w_0[2], 0.0, 
                                                  w_0[0], 0.0, w_0[2], 
                                                  0.0, w_0[0], w_0[1]);
      mat<double,mat_structure::square> del_sig_1(w_1[1], w_1[2], 0.0, 
                                                  w_1[0], 0.0, w_1[2], 
                                                  0.0, w_1[0], w_1[1]);
      mat<double,mat_structure::square> delw(sys_params.effective_J_inv * (R_0_1 * del_sig_0 - del_sig_1));
      sub(A)(w_r, s_r) += delw;
      sub(A)(q_r, s_r) += (0.5 * dt) * delw;
      // neglects the cross terms and any other non-trivial place where J appears.
      
      sub(A)(s_r, s_r) += mat_ident<double>(3);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InvCorrType, typename InputType>
    void apply_correction_to_state(const FlyWeight& params, const StateSpaceType& space, 
                                   const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                   typename pp::topology_traits<StateSpaceType>::point_type& x_c, 
                                   const InvCorrType& c, const InputType& u, const time_type& t) const {
      const point_type& s = params.get_state_models().template get_state_for_system<cross_inertia_state_model>(x);
      point_type& s_c = params.get_state_models().template get_state_for_system<cross_inertia_state_model>(x_c);
      
      s_c = s + vect<double,3>(c[inv_corr_start_index], c[inv_corr_start_index+1], c[inv_corr_start_index+2]);
    };
    
    template <typename FlyWeight, typename StateSpaceType, typename InputType, typename InvarFrameType>
    void set_invariant_frame_blocks(const FlyWeight& params, const StateSpaceType& space, 
                                    InvarFrameType& invar_frame, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_0, 
                                    const typename pp::topology_traits<StateSpaceType>::point_type& x_1, 
                                    const InputType& u, const time_type& t) const {
      /* identity is OK */
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(cross_inertia_state_model,0xC2310026,1,"cross_inertia_state_model",named_object)
    
};


class sat_position_output_model : public named_object {
  public:
    
    typedef vect_n<double> output_type;
    typedef vect_n<double> invariant_error_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    
  private:
    std::size_t start_index;
    std::size_t inv_start_index;
    
  public:
    
    sat_position_output_model() : start_index(0), inv_start_index(0) { };
    
    void construct_output_dimensions(std::size_t& cur_dim, std::size_t& cur_inv_dim) {
      start_index = cur_dim;
      cur_dim += 3;
      inv_start_index = cur_inv_dim;
      cur_inv_dim += 3;
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void set_output_from_state(const FlyWeight& params,
                               const StateSpaceType& space, 
                               const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                               output_type& y, time_type t) const {
      typedef satellite_state_model::point_type SE3State;
      const SE3State& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      
      y[range(start_index, start_index+3)] = get_position(x_se3);
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void set_inv_err_from_output(const FlyWeight& params,
                                 const StateSpaceType& space, 
                                 const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                 const output_type& y, invariant_error_type& e, time_type t) const {
      typedef satellite_state_model::point_type SE3State;
      const SE3State& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      
      e[range(inv_start_index, inv_start_index+3)] = y[range(start_index, start_index+3)] - get_position(x_se3);
    };
    
    template <typename MatrixC, typename MatrixD, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_output_function_blocks(MatrixC& C, MatrixD& D, 
                                    const FlyWeight& params, const StateSpaceType& space, 
                                    time_type t,
                                    const typename pp::topology_traits<StateSpaceType>::point_type& p, 
                                    const InputType& u) const {
      
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      const std::pair<std::size_t, std::size_t> p_r(sat3d_state_index, sat3d_state_index+3);
      const std::pair<std::size_t, std::size_t> pm_r(inv_start_index, inv_start_index+3);
      
      sub(C)(pm_r, p_r) += mat_ident<double>(3);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(sat_position_output_model,0xC2310029,1,"sat_position_output_model",named_object)
    
};


class sat_quaternion_output_model : public named_object {
  public:
    
    typedef vect_n<double> output_type;
    typedef vect_n<double> invariant_error_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    
  private:
    std::size_t start_index;
    std::size_t inv_start_index;
    
  public:
    
    sat_quaternion_output_model() : start_index(0), inv_start_index(0) { };
    
    void construct_output_dimensions(std::size_t& cur_dim, std::size_t& cur_inv_dim) {
      start_index = cur_dim;
      cur_dim += 4;
      inv_start_index = cur_inv_dim;
      cur_inv_dim += 3;
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void set_output_from_state(const FlyWeight& params,
                               const StateSpaceType& space, 
                               const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                               output_type& y, time_type t) const {
      typedef satellite_state_model::point_type SE3State;
      const SE3State& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      
      y[range(start_index, start_index+4)] = get_quaternion(x_se3);
    };
    
    template <typename FlyWeight, typename StateSpaceType>
    void set_inv_err_from_output(const FlyWeight& params,
                                 const StateSpaceType& space, 
                                 const typename pp::topology_traits<StateSpaceType>::point_type& x, 
                                 const output_type& y, invariant_error_type& e, time_type t) const {
      typedef satellite_state_model::point_type SE3State;
      const SE3State& x_se3 = params.get_state_models().template get_state_for_system<satellite_state_model>(x);
      
      unit_quat<double> q_diff = invert(get_quaternion(x_se3))
                               * unit_quat<double>(y[start_index],y[start_index+1],y[start_index+2],y[start_index+3]);
      vect<double,3> a = 2.0 * log(q_diff);
      
      e[range(inv_start_index, inv_start_index+3)] = a;
    };
    
    template <typename MatrixC, typename MatrixD, typename FlyWeight, typename StateSpaceType, typename InputType>
    void add_output_function_blocks(MatrixC& C, MatrixD& D, 
                                    const FlyWeight& params, const StateSpaceType& space, 
                                    time_type t,
                                    const typename pp::topology_traits<StateSpaceType>::point_type& p, 
                                    const InputType& u) const {
      const std::size_t sat3d_state_index = params.get_state_models().template get_system<satellite_state_model>().get_inv_corr_start_index();
      const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index+6, sat3d_state_index+9);
      const std::pair<std::size_t, std::size_t> qm_r(inv_start_index, inv_start_index+3);
      
      sub(C)(qm_r, q_r) += mat_ident<double>(3);  // TODO Add a frame transition ? (in invariant posterior frame)
    };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(sat_quaternion_output_model,0xC231002A,1,"sat_quaternion_output_model",named_object)
        
};



};

};

#endif




