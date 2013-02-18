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

#include "tuple_distance_metrics.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "ctrl_sys/linear_ss_system_concept.hpp"

#include "IHAQR_topology.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_are_solver.hpp"

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
      
      struct input_type { };
      struct output_type { };
      
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
      BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 0);
      BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 0);
      
      typedef typename IHAQR_topology< StateSpace, StateSpaceSystem >::linearization_payload lin_payload;
      typedef typename IHAQR_topology< StateSpace, StateSpaceSystem >::IHAQR_payload IHAQR_payload;
      
    private:
      lin_payload* lin_data;
      IHAQR_payload* IHAQR_data;
      mat<double,mat_structure::square> BRBmatrix;
      vect_n<double> stat_drift;
      
    public:
      
      MEAQR_ZIR_system(lin_payload* aLinData = NULL, 
                       IHAQR_payload* aIHAQRData = NULL,
                       const mat<double,mat_structure::diagonal>& aR) : lin_data(aLinData), IHAQR_data(aIHAQRData) {
        if(!IHAQR_data || !lin_data)
          return;
        
        BRBmatrix = lin_data->B * invert(aR) * transpose_view(lin_data->B);
        stat_drift = to_vect<double>(IHAQR_data->c);
      };
      
      
      point_derivative_type get_state_derivative(const state_space_type&, const point_type& p, input_type, double) const {
        if(!IHAQR_data || !lin_data)
          return;
        
        const mat<double,mat_structure::square>& L = get<0>(p);
        const mat<double,mat_structure::square>& H = get<1>(p);
        const vect_n<double>& eta = get<2>(p);
        const vect_n<double>& zir = get<3>(p);
        
        opts_L.LT = true;
        opts_U.LT = true;
        opts_U.TRANSA = true; 
        % Compute derivative 
        vect_n<double> T_rhs = BRBmatrix * eta + stat_drift;
        mat_vect_adaptor< vect_n<double> > T_rhs_m(T_rhs);
        detail::backsub_Cholesky_impl(H, T_rhs);
        
        mat<double,mat_structure::square> Temp_L = lin_data->A * L;
        detail::forwardsub_L_impl(L, Temp_L);
        mat<double,mat_structure::square> F = BRBmatrix;
        detail::forwardsub_L_impl(L, F);
        F = transpose(F);
        detail::forwardsub_L_impl(L, F);
        F += Temp_L;
        F += transpose_view(Temp_L);
        //L_dot = L * tril(F);
        detail::inplace_lower_multiply_with_fill_impl(L, F);
        
        Temp_L = lin_data->A * M;
        detail::forwardsub_L_impl(H, Temp_L);
        mat<double,mat_structure::square> G = BRBmatrix;
        detail::forwardsub_L_impl(H, G);
        G = transpose(G);
        detail::forwardsub_L_impl(H, G);
        G -= Temp_L;
        G -= transpose_view(Temp_L);
        //H_dot     = H * tril(G);
        detail::inplace_lower_multiply_with_fill_impl(H, G);
        
        return point_derivative_type(
#ifdef RK_ENABLE_CXX11_FEATURES
          std::move(F),                              // L_dot
          std::move(G),                              // H_dot
#else
          F, G,
#endif
          T_rhs - transpose_view(lin_data->A) * eta, // eta_dot
          lin_data->A * zir + stat_drift);           // zir_dot
      };
      
      output_type get_output(const state_space_type&, const point_type&, input_type, double) const {
        return output_type();
      };
      
  };
  
  
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
    
    typedef typename detail::MEAQR_ZIR_system<StateSpace, StateSpaceSystem>::point_type MEAQR_point_type;
    
    struct MEAQR_payload {
      std::vector< std::pair<double, MEAQR_point_type> > pts;
    };
    
    struct point_type : IHAQR_point_type {
      mutable shared_ptr< MEAQR_payload > MEAQR_data;
      
      explicit point_type(const state_type& aX) : IHAQR_point_type(aX) { };
#ifdef RK_ENABLE_CXX11_FEATURES
      explicit point_type(state_type&& aX) : IHAQR_point_type(std::move(aX)) { };
#endif
    };
    
    typedef IHAQR_point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
  protected:
    
    shared_ptr< IHAQR_space_type > m_IHAQR_space;
    double m_MEAQR_data_step_size;
    
    void compute_MEAQR_data(const point_type& p) const {
      m_IHAQR_space->compute_IHAQR_data(p);
      
      
      
    };
    
    
  public:
    
    /**
     * Default constructor.
     */
    MEAQR_topology(const std::string& aName = "MEAQR_topology"
                   const shared_ptr< IHAQR_space_type >& aIHAQRSpace = shared_ptr< IHAQR_space_type >(),
                   double aMEAQRDataStepSize) : 
                   named_object(),
                   m_IHAQR_space(aIHAQRSpace),
                   m_MEAQR_data_step_size(aMEAQRDataStepSize) {
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
      vect_n<double> dx =  to_vect<double>(m_IHAQR_space->m_space.difference(b.x, a.x));
      return dx * b.M * dx;
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
      point_type result( get(random_sampler, m_IHAQR_space->m_space)(m_IHAQR_space->m_space) );
      m_IHAQR_space->compute_IHAQR_data(result);
      return result;
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
      return m_space.is_in_bounds(a.x);
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
      if(b.M.get_row_count() == 0)
        compute_IHAQR_data(b);
      state_type goal_point = m_space.move_position_toward(a.x, fraction, b.x);
      mat<double,mat_structure::rectangular> K = invert(m_R) * transpose_view(b.B);
      state_type x_current = a.x;
      state_type x_next = x_current;
      system_input_type u_prev = a.u - invert(m_R) * transpose_view(a.B) * a.eta;
      
      // while not reached (fly-by) the goal yet:
      double current_time = 0.0;
      while( ( current_time < m_max_time_horizon ) &&
             ( m_space.distance(x_current, goal_point) > m_goal_proximity_threshold ) ) {
        // compute the current IHAQR input
        system_input_type u_current = b.u - K * (b.M * to_vect<double>(m_space.difference(x_current, goal_point)) + b.eta);
        
        u_current = m_input_space.adjust(u_current, space.get_diff_to_boundary(u_current));
        
        system_input_type du_dt = (u_current - u_prev) * (1.0 / m_time_step);
        for(std::size_t i = 0; i < du_dt.size(); ++i)
          if(fabs(du_dt[i]) > m_input_bandwidth[i])
            u_current[i] = u_prev[i] + (du_dt[i] > 0.0 ? 1.0 : -1.0) * m_input_bandwidth[i] * m_time_step;
        
        constant_trajectory< vector_topology< system_input_type > > input_traj(u_current);
        
        // integrate for one time-step.
        ctrl::detail::dormand_prince45_integrate_impl(
          m_space,
          m_system,
          x_current,
          x_next,
          input_traj,
          current_time,
          current_time + m_time_step,
          m_time_step * 1e-2,
          1e-3,
          m_time_step * 1e-6,
          m_time_step * 0.1);
        x_current = x_next;
        current_time += m_time_step;
        u_prev = u_current;
      };
      
      point_type result( x_current );
      result.u = u_prev;
      return result;
    };

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_IHAQR_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_MEAQR_data_step_size);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_IHAQR_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_MEAQR_data_step_size);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240000C,1,"MEAQR_topology",named_object)
    
};

template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_space< MEAQR_topology<StateSpace, StateSpaceSystem> > : boost::mpl::true_ { };
        
template <typename StateSpace, typename StateSpaceSystem>
struct is_point_distribution< MEAQR_topology<StateSpace, StateSpaceSystem> > : 
  is_point_distribution<StateSpace> { };



};

};


#endif








