/**
 * \file MEAQR_topology.hpp
 * 
 * This library provides classes
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_MEAQR_TOPOLOGY_HPP
#define REAK_MEAQR_TOPOLOGY_HPP


#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT


#include "base/named_object.hpp"

#include "topologies/basic_distance_metrics.hpp"
#include "topologies/tuple_distance_metrics.hpp"
#include "topologies/hyperbox_topology.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "ctrl_sys/linear_ss_system_concept.hpp"

#include "IHAQR_topology.hpp"

#include "sys_integrators/runge_kutta4_integrator_sys.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_operators.hpp"
#include "lin_alg/mat_are_solver.hpp"
#include "lin_alg/mat_exp_methods.hpp"

namespace ReaK {

namespace pp {
  
namespace detail {
  
  
  template <typename StateSpace, typename StateSpaceSystem>
  class MEAQR_ZIR_system {
    public:
      
      typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
      typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
      typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
      typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;
      
      typedef arithmetic_tuple< 
        mat<double,mat_structure::square>,
        mat<double,mat_structure::square>,
        vect_n<double>,
        vect_n<double> > point_type;
      typedef point_type point_difference_type;
      typedef point_type point_derivative_type;
      
      typedef vector_topology< point_type > state_space_type;
      
      typedef double time_type;
      typedef double time_difference_type;
      
      typedef vect_n<double> input_type;
      typedef vect_n<double> output_type;
      
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
      BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 0);
      BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 0);
      
      typedef typename IHAQR_topology< StateSpace, StateSpaceSystem >::linearization_payload lin_payload;
      typedef typename IHAQR_topology< StateSpace, StateSpaceSystem >::IHAQR_payload IHAQR_payload;
      
    public:
      lin_payload* lin_data;
      IHAQR_payload* IHAQR_data;
      mat<double,mat_structure::square> BRBmatrix;
      vect_n<double> stat_drift;
      
    public:
      
      MEAQR_ZIR_system(lin_payload* aLinData = NULL, 
                       IHAQR_payload* aIHAQRData = NULL,
                       const mat<double,mat_structure::diagonal>& aR = (mat<double,mat_structure::diagonal>())) : lin_data(aLinData), IHAQR_data(aIHAQRData) {
        if(!IHAQR_data || !lin_data)
          return;
        
        BRBmatrix = lin_data->B * invert(aR) * transpose_view(lin_data->B);
        stat_drift = to_vect<double>(IHAQR_data->c);
      };
      
      
      point_derivative_type get_state_derivative(const state_space_type&, const point_type& p, input_type, double) const {
        if(!IHAQR_data || !lin_data)
          return point_derivative_type();
        
        const mat<double,mat_structure::square>& L = get<0>(p);
        const mat<double,mat_structure::square>& H = get<1>(p);
        const vect_n<double>& eta = get<2>(p);
        const vect_n<double>& zir = get<3>(p);
        
        // Compute derivative 
        vect_n<double> T_rhs = BRBmatrix * eta - stat_drift;
        mat_vect_adaptor< vect_n<double> > T_rhs_m(T_rhs);
        ReaK::detail::backsub_Cholesky_impl(H, T_rhs_m);
        
        mat<double,mat_structure::square> Temp_L = lin_data->A * L;
        ReaK::detail::forwardsub_L_impl(L, Temp_L, 1e-6);
        mat<double,mat_structure::square> F = BRBmatrix;
        ReaK::detail::forwardsub_L_impl(L, F, 1e-6);
        F = transpose(F);
        ReaK::detail::forwardsub_L_impl(L, F, 1e-6);
        F += Temp_L;
        F += transpose_view(Temp_L);
        for(std::size_t i = 0; i < F.get_row_count(); ++i)
          F(i,i) *= 0.5;
        //L_dot = L * tril(F);
        ReaK::detail::inplace_lower_multiply_with_fill_impl(L, F);
        
        Temp_L = lin_data->A * H;
        ReaK::detail::forwardsub_L_impl(H, Temp_L, 1e-6);
        mat<double,mat_structure::square> G = BRBmatrix;
        ReaK::detail::forwardsub_L_impl(H, G, 1e-6);
        G = transpose(G);
        ReaK::detail::forwardsub_L_impl(H, G, 1e-6);
        G -= Temp_L;
        G -= transpose_view(Temp_L);
        for(std::size_t i = 0; i < G.get_row_count(); ++i)
          G(i,i) *= 0.5;
        //H_dot     = H * tril(G);
        ReaK::detail::inplace_lower_multiply_with_fill_impl(H, G);
        
        return point_derivative_type(
#ifdef RK_ENABLE_CXX11_FEATURES
          std::move(F),                              // L_dot
          std::move(G),                              // H_dot
#else
          F, G,
#endif
          transpose_view(lin_data->A) * eta - T_rhs, // eta_dot
          lin_data->A * zir + stat_drift);           // zir_dot
      };
      
      output_type get_output(const state_space_type&, const point_type&, input_type, double) const {
        return output_type();
      };
      
  };
  
  
};



template <typename StateSpace, typename StateSpaceSystem>
class MEAQR_point_type : public IHAQR_point_type<StateSpace, StateSpaceSystem> {
  public:
    typedef IHAQR_point_type<StateSpace, StateSpaceSystem> base_type;
    typedef MEAQR_point_type<StateSpace, StateSpaceSystem> self;
    
    typedef typename detail::MEAQR_ZIR_system<StateSpace, StateSpaceSystem>::point_type MEAQR_bundle_type;
    
    typedef typename base_type::matrixA_type matrixA_type;
    typedef typename base_type::matrixB_type matrixB_type;
    
    typedef typename base_type::state_type state_type;
    typedef typename base_type::state_difference_type state_difference_type;
    typedef typename base_type::state_derivative_type state_derivative_type;
    typedef typename base_type::system_input_type system_input_type;
    
    typedef typename base_type::linearization_payload linearization_payload;
    typedef typename base_type::IHAQR_payload IHAQR_payload;
    
    struct MEAQR_payload {
      std::vector< std::pair<double, MEAQR_bundle_type> > pts;
    };
    
    mutable shared_ptr< MEAQR_payload > MEAQR_data;
    
    MEAQR_point_type(const base_type& rhs) : base_type(rhs) { };
    explicit MEAQR_point_type(const state_type& aX = state_type()) : base_type(aX) { };
    
    MEAQR_point_type& operator()(const base_type& rhs) {
      this->x = rhs.x;
      this->lin_data = rhs.lin_data;
      this->IHAQR_data = rhs.IHAQR_data;
      this->MEAQR_data = shared_ptr< MEAQR_payload >();
      return *this;
    };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    MEAQR_point_type(base_type&& rhs) : base_type(std::move(rhs)) { };
    explicit MEAQR_point_type(state_type&& aX) : base_type(std::move(aX)) { };
    
    MEAQR_point_type& operator()(base_type&& rhs) {
      this->x = std::move(rhs.x);
      this->lin_data = std::move(rhs.lin_data);
      this->IHAQR_data = std::move(rhs.IHAQR_data);
      this->MEAQR_data = shared_ptr< MEAQR_payload >();
      return *this;
    };
#endif
    
    
    virtual ~MEAQR_point_type() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400034,1,"MEAQR_point_type",base_type)

  
};





/**
 * This class 
 */
template <typename StateSpace, typename StateSpaceSystem>
class MEAQR_topology : public named_object
{
  public:
    typedef MEAQR_topology<StateSpace, StateSpaceSystem> self;
    
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixA_type matrixA_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixB_type matrixB_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixC_type matrixC_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixD_type matrixD_type;
    
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;
    
    typedef IHAQR_topology< StateSpace, StateSpaceSystem > IHAQR_space_type;
    typedef typename IHAQR_space_type::point_type IHAQR_point_type;
    typedef typename IHAQR_space_type::point_difference_type IHAQR_point_difference_type;
    
    typedef typename detail::MEAQR_ZIR_system<StateSpace, StateSpaceSystem>::point_type MEAQR_bundle_type;
    
    typedef MEAQR_point_type<StateSpace, StateSpaceSystem> point_type;
    typedef typename point_type::MEAQR_payload MEAQR_payload;
    
    typedef IHAQR_point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  protected:
    
    shared_ptr< IHAQR_space_type > m_IHAQR_space;
    double m_MEAQR_data_step_size;
    double m_idle_power_cost;
    
    void compute_MEAQR_data(const point_type& p) const {
      m_IHAQR_space->compute_IHAQR_data(p);
      if(p.MEAQR_data)
        return;
      p.MEAQR_data = shared_ptr< MEAQR_payload >(new MEAQR_payload());
      
      detail::MEAQR_ZIR_system<StateSpace, StateSpaceSystem> MEAQR_sys(
        p.lin_data.get(), 
        p.IHAQR_data.get(),
        m_IHAQR_space->m_R);
      
      vector_topology< MEAQR_bundle_type > MEAQR_sys_space;
      
      double init_damping_value = 1e-5 * norm_inf(MEAQR_sys.BRBmatrix);
      std::size_t N = p.lin_data->A.get_row_count();
      MEAQR_bundle_type start_point(
        mat<double,mat_structure::square>( mat<double,mat_structure::scalar>(N, init_damping_value) ),
        mat<double,mat_structure::square>( mat<double,mat_structure::scalar>(N, 100.0 * init_damping_value) ),
        vect_n<double>(N, 0.0),
        vect_n<double>(N, 0.0)
      );
      
      ReaK::detail::decompose_Cholesky_impl(MEAQR_sys.BRBmatrix * m_MEAQR_data_step_size + get<0>(start_point), get<0>(start_point), 1e-8);
      get<1>(start_point) = get<0>(start_point);
      
      MEAQR_bundle_type current_point = start_point;
      
      typedef typename detail::MEAQR_ZIR_system<StateSpace, StateSpaceSystem>::input_type InputType;
      constant_trajectory< vector_topology< InputType > > input_traj = constant_trajectory< vector_topology< InputType > >(InputType());
      
      // while not reached (fly-by) the goal yet:
      double current_time = 0.0;
      p.MEAQR_data->pts.push_back(std::make_pair(current_time, current_point));
      while( current_time < m_IHAQR_space->m_max_time_horizon )  {
        // integrate for one time-step.
//         ReaK::ctrl::detail::dormand_prince45_integrate_impl(
//           MEAQR_sys_space,
//           MEAQR_sys,
//           start_point,
//           current_point,
//           input_traj,
//           current_time,
//           current_time + m_MEAQR_data_step_size,
//           m_MEAQR_data_step_size * 1e-2,
//           1e-3,
//           m_MEAQR_data_step_size * 1e-6,
//           m_MEAQR_data_step_size * 0.1);
        ReaK::ctrl::detail::runge_kutta4_integrate_impl(
          MEAQR_sys_space,
          MEAQR_sys,
          start_point,
          current_point,
          input_traj,
          current_time,
          current_time + m_MEAQR_data_step_size,
          m_MEAQR_data_step_size * 1e-2);
        start_point = current_point;
        current_time += m_MEAQR_data_step_size;
        p.MEAQR_data->pts.push_back(std::make_pair(current_time, current_point));
      };
      
      return;
    };
    
    virtual bool is_free_impl(const state_type& a) const {
      return true;
    };
    
    
    
    point_type move_position_toward_impl(const point_type& a, double fraction, const point_type& b, bool with_collision_check) const 
    {
      if(!a.lin_data)
        m_IHAQR_space->compute_linearization_data(a);
      if(!b.MEAQR_data)
        compute_MEAQR_data(b);
      // first, find the T-optimal:
      
      typedef typename std::vector< std::pair<double, MEAQR_bundle_type> >::const_iterator Iter;
      
      double min_J = std::numeric_limits<double>::infinity();
      Iter min_it = b.MEAQR_data->pts.end();
      
      mat<double,mat_structure::square> eAdt(b.lin_data->A.get_row_count());
      ReaK::exp_PadeSAS(b.lin_data->A * m_MEAQR_data_step_size, eAdt, QR_linlsqsolver(), 1e-6);
      mat<double,mat_structure::square> eAt = eAdt;
      
      vect_n<double> a_rel = to_vect<double>(m_IHAQR_space->m_space.difference(a.x, b.x));
      
      bool found_T_optimal = false;
      for(Iter it = b.MEAQR_data->pts.begin(); it != b.MEAQR_data->pts.end(); ++it) {
        
        vect_n<double> s = eAt * a_rel + get<3>(it->second);  // ZIR(a_rel, t);
        
        // s = invert(L) * s
        mat_vect_adaptor< vect_n<double> > s_m(s);
//         std::cout << " t = \t" << std::setw(10) << it->first << std::endl;
//         std::cout << " L = \t" << get<0>(it->second) << std::endl;
        ReaK::detail::forwardsub_L_impl(get<0>(it->second), s_m, 1e-6);
        
        double J = m_idle_power_cost * it->first + 0.5 * (s * s); 
//         std::cout << " J = \t" << std::setw(10) << J << std::endl;
        if(J < min_J) {
          min_J = J;
          min_it = it;
        };
        if(min_J < m_idle_power_cost * it->first) {
          found_T_optimal = true;
          break;
        };
        
        eAt = eAdt * eAt;
      };
      if(!found_T_optimal)
        return point_type( m_IHAQR_space->move_position_toward(a, fraction, b).x );
      ++min_it;
      double T_optimal = min_it->first;
//       std::cout << " t_optimal = " << T_optimal << std::endl;
      
//       RK_NOTICE(1," reached!");
      state_type goal_point = m_IHAQR_space->m_space.move_position_toward(a.x, fraction, b.x);
//       RK_NOTICE(1," reached!");
      mat<double,mat_structure::rectangular> K = invert(m_IHAQR_space->m_R) * transpose_view(b.lin_data->B);
//       RK_NOTICE(1," reached!");
      state_type x_current = a.x;
      state_type x_next = x_current;
//       RK_NOTICE(1," reached!");
      system_input_type u_prev = a.lin_data->u;
//       RK_NOTICE(1," reached!");
      m_IHAQR_space->m_input_space.bring_point_in_bounds(u_prev);
//       RK_NOTICE(1," reached!");
      
      // Iterate through the MEAQR sequence points.
      double current_time = 0.0;
      double accum_steer_cost = 0.0;
      while(min_it != b.MEAQR_data->pts.begin()) {
        Iter prev_it = min_it;
        --min_it;
        
//       RK_NOTICE(1," reached!");
        // integrate for some time steps.
        
        bool ended_in_collision = false;
        // while not reached (fly-by) the goal yet:
        while( ( current_time < T_optimal - min_it->first ) &&
               ( m_IHAQR_space->m_space.distance(x_current, goal_point) > m_IHAQR_space->m_goal_proximity_threshold ) ) {
          // compute the current MEAQR input
          double prev_fraction = (current_time - T_optimal + prev_it->first) / m_MEAQR_data_step_size;
          double next_fraction = (T_optimal - min_it->first - current_time) / m_MEAQR_data_step_size;
          mat<double,mat_structure::square> H = prev_fraction * get<1>(prev_it->second) + next_fraction * get<1>(min_it->second);
          vect_n<double> eta = prev_fraction * get<2>(prev_it->second) + next_fraction * get<2>(min_it->second);
          std::cout << "eta = " << eta << std::endl;
//       RK_NOTICE(1," reached!");
          
          vect_n<double> HHx = to_vect<double>(m_IHAQR_space->m_space.difference(x_current, b.x));
          std::cout << "dx = " << HHx << std::endl;
          std::cout << "H = " << H << std::endl;
          mat_vect_adaptor< vect_n<double> > HHx_m(HHx);
          ReaK::detail::backsub_Cholesky_impl(H, HHx_m);
          std::cout << "M * dx = " << HHx << std::endl;
//       RK_NOTICE(1," reached!");
          
          system_input_type u_current = b.lin_data->u - K * (HHx + eta);
          std::cout << " u_current (before bounding) = " << u_current << std::endl;
          m_IHAQR_space->m_input_space.bring_point_in_bounds(u_current);
//       RK_NOTICE(1," reached!");
          
          system_input_type du_dt = (u_current - u_prev) * (1.0 / m_IHAQR_space->m_time_step);
          m_IHAQR_space->m_input_rate_space.bring_point_in_bounds(du_dt);
          u_current = u_prev + m_IHAQR_space->m_time_step * du_dt;
//       RK_NOTICE(1," reached!");
          std::cout << " u_current (after bounding) = " << u_current << std::endl;
          
          accum_steer_cost += to_vect<double>(u_current) * (m_IHAQR_space->m_R * to_vect<double>(u_current)) * m_IHAQR_space->m_time_step;
          
//       RK_NOTICE(1," reached!");
          constant_trajectory< vector_topology< system_input_type > > input_traj(u_current);
          
          // integrate for one time-step.
//           ReaK::ctrl::detail::dormand_prince45_integrate_impl(
//             m_IHAQR_space->m_space,
//             *(m_IHAQR_space->m_system),
//             x_current,
//             x_next,
//             input_traj,
//             current_time,
//             current_time + m_IHAQR_space->m_time_step,
//             m_IHAQR_space->m_time_step * 1e-2,
//             1e-3,
//             m_IHAQR_space->m_time_step * 1e-6,
//             m_IHAQR_space->m_time_step * 0.1);
          ReaK::ctrl::detail::runge_kutta4_integrate_impl(
            m_IHAQR_space->m_space,
            *(m_IHAQR_space->m_system),
            x_current,
            x_next,
            input_traj,
            current_time,
            current_time + m_IHAQR_space->m_time_step,
            m_IHAQR_space->m_time_step * 1e-2);
//       RK_NOTICE(1," reached!");
          
          if((!with_collision_check) || is_free_impl(x_next)) {
            x_current = x_next;
            current_time += m_IHAQR_space->m_time_step;
            std::cout << " time = " << current_time << "  state = " << x_current << std::endl;
            u_prev = u_current;
          } else {
            ended_in_collision = true;
            break;
          };
        };
        if(ended_in_collision)
          break;
//       RK_NOTICE(1," reached!");
      
      };
//       RK_NOTICE(1," reached!");
      
      point_type result( x_current );
      return result;
    };
    
    point_type random_point_impl(bool with_collision_check) const {
      state_type result_pt = get(random_sampler, m_IHAQR_space->m_space)(m_IHAQR_space->m_space);
      while(with_collision_check && !is_free_impl(result_pt))
        result_pt = get(random_sampler, m_IHAQR_space->m_space)(m_IHAQR_space->m_space);
      point_type result( result_pt );
      m_IHAQR_space->compute_IHAQR_data(result);
      compute_MEAQR_data(result);
      return result;
    };
    
  public:
    
    StateSpace& get_state_space() { return m_IHAQR_space->get_state_space(); };
    const StateSpace& get_state_space() const { return m_IHAQR_space->get_state_space(); };
    
    IHAQR_space_type& get_IHAQR_space() { return *m_IHAQR_space; };
    const IHAQR_space_type& get_IHAQR_space() const { return *m_IHAQR_space; };
    
    
    
    /**
     * Default constructor.
     */
    MEAQR_topology(const std::string& aName = "MEAQR_topology",
                   const shared_ptr< IHAQR_space_type >& aIHAQRSpace = shared_ptr< IHAQR_space_type >(),
                   double aMEAQRDataStepSize = 0.1,
                   double aIdlePowerCost = 1.0) : 
                   named_object(),
                   m_IHAQR_space(aIHAQRSpace),
                   m_MEAQR_data_step_size(aMEAQRDataStepSize),
                   m_idle_power_cost(aIdlePowerCost) {
      setName(aName);
    };
    
    virtual ~MEAQR_topology() { };
    
    
    /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const {
      if(!b.MEAQR_data)
        compute_MEAQR_data(b);
      
      double min_J = std::numeric_limits<double>::infinity();
      
      typedef typename std::vector< std::pair<double, MEAQR_bundle_type> >::const_iterator Iter;
      
      mat<double,mat_structure::square> eAdt(b.lin_data->A.get_row_count());
      ReaK::exp_PadeSAS(b.lin_data->A * m_MEAQR_data_step_size, eAdt, QR_linlsqsolver(), 1e-6);
      mat<double,mat_structure::square> eAt = eAdt;
      
      vect_n<double> a_rel = to_vect<double>(m_IHAQR_space->m_space.difference(a.x, b.x));
      
      for(Iter it = b.MEAQR_data->pts.begin(); it != b.MEAQR_data->pts.end(); ++it) {
        
        vect_n<double> s = eAt * a_rel + get<3>(it->second);
        
        // s = invert(L) * s
        mat_vect_adaptor< vect_n<double> > s_m(s);
        ReaK::detail::forwardsub_L_impl(get<0>(it->second), s_m, 1e-6);
        
        double J = m_idle_power_cost * it->first + 0.5 * (s * s); 
        if(J < min_J)
          min_J = J;
        if(min_J < m_idle_power_cost * it->first)
          break;
        
        eAt = eAdt * eAt;
      };
      
      return min_J;
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      vect_n<double> dx = to_vect<double>(delta.dx);
      return dx * dx;
    };
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      return random_point_impl(false);
    };
    
   /*************************************************************************
    *                             TopologyConcept
    * **********************************************************************/

    /**
     * Returns the difference between two points (analogous to a - b, but implemented in SO(3) Lie algebra).
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      return point_difference_type( m_IHAQR_space->m_space.difference(b.x, a.x) );
    };

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return point_type( m_IHAQR_space->m_space.adjust(a.x, delta.dx) );
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return point_type( m_IHAQR_space->m_space.origin() );
    };
    
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      return m_IHAQR_space->m_space.is_in_bounds(a.x);
    };
    
    // NOTE: don't know if I can get rid of this. (only seems useful in bounded interpolators (and samplers)).
    point_difference_type get_diff_to_boundary(const point_type&) const {
      return point_difference_type();
    };

    /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return move_position_toward_impl(a, fraction, b, false);
    };

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_IHAQR_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_MEAQR_data_step_size)
        & RK_SERIAL_SAVE_WITH_NAME(m_idle_power_cost);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_IHAQR_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_MEAQR_data_step_size)
        & RK_SERIAL_LOAD_WITH_NAME(m_idle_power_cost);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400035,1,"MEAQR_topology",named_object)
    
};

template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_space< MEAQR_topology<StateSpace, StateSpaceSystem> > : boost::mpl::true_ { };

template <typename StateSpace, typename StateSpaceSystem>
struct is_point_distribution< MEAQR_topology<StateSpace, StateSpaceSystem> > : 
  is_point_distribution<StateSpace> { };

template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_symmetric< MEAQR_topology<StateSpace, StateSpaceSystem> > : boost::mpl::false_ { };
  
  
  
  
/**
 * This class has collision detection.
 */
template <typename StateSpace, typename StateSpaceSystem>
class MEAQR_topology_with_CD : public MEAQR_topology<StateSpace, StateSpaceSystem> {
  public:
    typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem> self;
    typedef MEAQR_topology<StateSpace, StateSpaceSystem> base_type;
    
    typedef base_type super_space_type;
    
    typedef typename base_type::state_type state_type;
    typedef typename base_type::state_difference_type state_difference_type;
    typedef typename base_type::state_derivative_type state_derivative_type;
    typedef typename base_type::system_input_type system_input_type;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef typename base_type::IHAQR_space_type IHAQR_space_type;
    
  protected:
    shared_ptr< kte::direct_kinematics_model > m_model; 
    
    virtual bool is_free_impl(const state_type& a) const {
      
      detail::write_joint_coordinates_impl(a, this->m_IHAQR_space->get_state_space(), m_model);
      // update the kinematics model with the given joint states.
      m_model->doDirectMotion();
      
      for( std::vector< shared_ptr< geom::proxy_query_pair_2D > >::const_iterator it = m_proxy_env_2D.begin(); it != m_proxy_env_2D.end(); ++it) {
        shared_ptr< geom::proximity_finder_2D > tmp = (*it)->findMinimumDistance();
        if((tmp) && (tmp->getLastResult().mDistance < 0.0))
          return false;
      };
      for( std::vector< shared_ptr< geom::proxy_query_pair_3D > >::const_iterator it = m_proxy_env_3D.begin(); it != m_proxy_env_3D.end(); ++it) {
        shared_ptr< geom::proximity_finder_3D > tmp = (*it)->findMinimumDistance();
        if((tmp) && (tmp->getLastResult().mDistance < 0.0))
          return false;
      };
      
      return true;
    };
    
  public:
    
    std::vector< shared_ptr< geom::proxy_query_pair_2D > > m_proxy_env_2D;
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > m_proxy_env_3D;
    
    /**
     * Default constructor.
     */
    MEAQR_topology_with_CD(const std::string& aName = "MEAQR_topology_with_CD",
                           const shared_ptr< IHAQR_space_type >& aIHAQRSpace = shared_ptr< IHAQR_space_type >(),
                           double aMEAQRDataStepSize = 0.1,
                           double aIdlePowerCost = 1.0,
                           const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >()) : 
                           base_type(aName, aIHAQRSpace, aMEAQRDataStepSize, aIdlePowerCost),
                           m_model(aModel),
                           m_proxy_env_2D(),
                           m_proxy_env_3D() { };
    
    virtual ~MEAQR_topology_with_CD() { };
    
    super_space_type& get_super_space() { return *this; };
    
    const super_space_type& get_super_space() const { return *this; };
    
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      return this->random_point_impl(true);
    };
    
    
    bool is_free(const point_type& a) const {
      return this->is_free_impl(a);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return this->move_position_toward_impl(a, fraction, b, true);
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_model)
        & RK_SERIAL_SAVE_WITH_NAME(m_proxy_env_2D)
        & RK_SERIAL_SAVE_WITH_NAME(m_proxy_env_3D);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_model)
        & RK_SERIAL_LOAD_WITH_NAME(m_proxy_env_2D)
        & RK_SERIAL_LOAD_WITH_NAME(m_proxy_env_3D);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400036,1,"MEAQR_topology_with_CD",base_type)
    
    
};
  
  
template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_space< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem> > : boost::mpl::true_ { };

template <typename StateSpace, typename StateSpaceSystem>
struct is_point_distribution< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem> > : 
  is_point_distribution<StateSpace> { };

template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_symmetric< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem> > : boost::mpl::false_ { };
  




};

};


#endif








