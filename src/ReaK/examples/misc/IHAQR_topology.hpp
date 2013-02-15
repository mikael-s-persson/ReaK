/**
 * \file IHAQR_topology.hpp
 * 
 * This library provides classes that define topologies on a system controlled by an infinite-horizon affine
 * quadratic regulator.
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

#ifndef REAK_IHAQR_TOPOLOGY_HPP
#define REAK_IHAQR_TOPOLOGY_HPP


#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT


#include "base/named_object.hpp"

#include "tuple_distance_metrics.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "ctrl_sys/linear_ss_system_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_are_solver.hpp"

namespace ReaK {

namespace pp {

/**
 * This class implements a quaternion-topology. Because quaternions are constrained on the unit 
 * hyper-sphere, this topology is indeed bounded (yet infinite at the same time). This class
 * models the MetricSpaceConcept, the LieGroupConcept, and the PointDistributionConcept.
 */
template <typename StateSpace, typename StateSpaceSystem>
class IHAQR_topology : public named_object
{
  public:
    typedef IHAQR_topology<StateSpace, StateSpaceSystem> self;
    
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixA_type matrixA_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixB_type matrixB_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixC_type matrixC_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixD_type matrixD_type;
    
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;
    
    struct point_type {
      state_type x;
      mutable system_input_type u;
      mutable matrixA_type A;
      mutable matrixB_type B;
      mutable mat<double,mat_structure::square> M;
      mutable state_derivative_type c;
      mutable state_derivative_type eta;
      
      explicit point_type(const state_type& aX) : x(aX) { };
#ifdef RK_ENABLE_CXX11_FEATURES
      explicit point_type(state_type&& aX) : x(std::move(aX)) { };
#endif
      
      
    };
    
    struct point_difference_type {
      state_difference_type dx;
      
      explicit point_difference_type(const state_difference_type& aDX) : dx(aDX) { };
#ifdef RK_ENABLE_CXX11_FEATURES
      explicit point_difference_type(state_difference_type&& aDX) : dx(std::move(aDX)) { };
#endif
      
    };
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    
    
  
  protected:
    
    shared_ptr< StateSpaceSystem > m_system;
    StateSpace m_space;
    hyperbox_topology< system_input_type > m_input_space;
    
    mat<double,mat_structure::diagonal> m_R;
    mat<double,mat_structure::diagonal> m_Q;
    double m_time_step;
    system_input_type m_input_bandwidth;
    double m_max_time_horizon;
    double m_goal_proximity_threshold;
    
    void compute_IHAQR_data(const point_type& a) const {
      // compute u
      matrixC_type Ctmp;
      matrixD_type Dtmp;
      a.u = system_input_type();
      vect_n<double> f = to_vect<double>(m_system->get_state_derivative(m_space, a.x, a.u, 0.0));
      mat_vect_adaptor< vect_n<double> > f_m(f);
      m_system->get_linear_blocks(a.A, a.B, Ctmp, Dtmp, m_space, 0.0, a.x, a.u);
      mat<double,mat_structure::rectangular> u_m(to_vect<double>(a.u).size(), 1);
      linlsq_QR(a.B, u_m, f_m);
      a.u = from_vect< system_input_type >(u_m);
      
      // compute c
      a.c = m_system->get_state_derivative(m_space, a.x, a.u, 0.0);
      vect_n<double> c_v = to_vect<double>(a.c);
      
      // fill in A and B
      m_system->get_linear_blocks(a.A, a.B, Ctmp, Dtmp, m_space, 0.0, a.x, a.u);
      
      // solve for M
      try {
        solve_care_problem(a.A, a.B, m_R, m_Q, a.M);
      } catch(std::exception& e) {
        std::cout << "Warning! Solution to the CARE problem could not be found for the given state point: " << a.x << std::endl
                  << "  The following exception was raised: " << e.what() << std::endl;
      };
      
      // solve for eta
      // To be neglected in distance metric
      mat<double,mat_structure::square> lhsMat = a.M * a.B * invert(m_R) * transpose_view(a.B) - transpose_view(a.A);
      mat_vect_adaptor< vect_n<double> > c_v_m(c_v);
      mat<double,mat_structure::rectangular> eta_m(c_v.size(), 1);
      linlsq_QR(lhsMat, eta_m, c_v_m);
      a.eta = from_vect<state_derivative_type>(mat_row_slice< mat<double,mat_structure::rectangular> >(eta_m));
      
    };
    
  public:
    
    /**
     * Default constructor.
     */
    IHAQR_topology(const std::string& aName = "IHAQR_topology"
                   const shared_ptr< StateSpaceSystem >& aSystem = shared_ptr< StateSpaceSystem >(),
                   const StateSpace& aSpace = StateSpace(),
                   const system_input_type& aMinInput = system_input_type(),
                   const system_input_type& aMaxInput = system_input_type(),
                   const mat<double,mat_structure::diagonal>& aR = (mat<double,mat_structure::diagonal>()),
                   const mat<double,mat_structure::diagonal>& aQ = (mat<double,mat_structure::diagonal>()),
                   double aTimeStep = 0.1,
                   const system_input_type& aInputBandwidth = system_input_type(),
                   double aMaxTimeHorizon = 10.0,
                   double aGoalProximityThreshold = 1.0) : 
                   named_object(),
                   m_system(aSystem),
                   m_space(aSpace),
                   m_input_space(aName + "_input_space", aMinInput, aMaxInput),
                   m_R(aR),
                   m_Q(aQ),
                   m_time_step(aTimeStep),
                   m_input_bandwidth(aInputBandwidth),
                   m_max_time_horizon(aMaxTimeHorizon),
                   m_goal_proximity_threshold(aGoalProximityThreshold) {
      setName(aName);
    };
    
    virtual ~IHAQR_topology() { };
    
    
    /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const {
      if(b.M.get_row_count() == 0)
        compute_IHAQR_data(b);
      vect_n<double> dx =  to_vect<double>(m_space.difference(b.x, a.x));
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
      point_type result( get(random_sampler, m_space)(m_space) );
      compute_IHAQR_data(result);
      return result;
    };
    
   /*************************************************************************
    *                             TopologyConcept
    * **********************************************************************/

    /**
     * Returns the difference between two points (analogous to a - b, but implemented in SO(3) Lie algebra).
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      point_difference_type result( m_space.difference(b.x, a.x) );
      return result;
    };

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      point_type result( m_space.adjust(a.x, delta.dx) );
      return result;
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      point_type result( m_space.origin() );
      return result;
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
      A & RK_SERIAL_SAVE_WITH_NAME(m_system)
        & RK_SERIAL_SAVE_WITH_NAME(m_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_input_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_R)
        & RK_SERIAL_SAVE_WITH_NAME(m_Q)
        & RK_SERIAL_SAVE_WITH_NAME(m_time_step)
        & RK_SERIAL_SAVE_WITH_NAME(m_input_bandwidth)
        & RK_SERIAL_SAVE_WITH_NAME(m_max_time_horizon)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_proximity_threshold);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_system)
        & RK_SERIAL_LOAD_WITH_NAME(m_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_R)
        & RK_SERIAL_LOAD_WITH_NAME(m_Q)
        & RK_SERIAL_LOAD_WITH_NAME(m_time_step)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_bandwidth)
        & RK_SERIAL_LOAD_WITH_NAME(m_max_time_horizon)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_proximity_threshold);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240000C,1,"IHAQR_topology",named_object)
    
};

template <typename StateSpace, typename StateSpaceSystem>
struct is_metric_space< IHAQR_topology<StateSpace, StateSpaceSystem> > : boost::mpl::true_ { };
        
template <typename StateSpace, typename StateSpaceSystem>
struct is_point_distribution< IHAQR_topology<StateSpace, StateSpaceSystem> > : 
  is_point_distribution<StateSpace> { };



};

};


#endif








