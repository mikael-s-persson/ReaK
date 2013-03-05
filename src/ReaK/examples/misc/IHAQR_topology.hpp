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

#include "topologies/tuple_distance_metrics.hpp"
#include "topologies/basic_distance_metrics.hpp"
#include "topologies/hyperbox_topology.hpp"

#include "interpolation/constant_trajectory.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "ctrl_sys/linear_ss_system_concept.hpp"

#include "sys_integrators/dormand_prince45_integrator_sys.hpp"
#include "sys_integrators/runge_kutta4_integrator_sys.hpp"

#include "topologies/direct_kinematics_topomap.hpp"       // for write_joint_coordinates_impl
#include "kte_models/direct_kinematics_model.hpp"
#include "proximity/proxy_query_model.hpp"  // for proxy-query class

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_are_solver.hpp"

namespace ReaK {

namespace pp {


// forward-declaration.
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
class MEAQR_topology;



template <typename StateSpace, typename StateSpaceSystem>
class IHAQR_point_type : public shared_object {
  public:
    typedef IHAQR_point_type<StateSpace, StateSpaceSystem> self;
    
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixA_type matrixA_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixB_type matrixB_type;
    
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;
    
    struct linearization_payload {
      system_input_type u;
      matrixA_type A;
      matrixB_type B;
      state_derivative_type c;
    };
    
    struct IHAQR_payload {
      mat<double,mat_structure::square> M;      // IH-LQR cost-to-go matrix.
      mat<double,mat_structure::rectangular> K; // IH-LQR optimal gain matrix.
      system_input_type u_bias;
    };
    
    state_type x;
    mutable shared_ptr< linearization_payload > lin_data;
    mutable shared_ptr< IHAQR_payload > IHAQR_data;
    
    explicit IHAQR_point_type(const state_type& aX = state_type()) : x(aX) { };
#ifdef RK_ENABLE_CXX11_FEATURES
    explicit IHAQR_point_type(state_type&& aX) : x(std::move(aX)) { };
#endif
    
    virtual ~IHAQR_point_type() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(x);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(x);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400031,1,"IHAQR_point_type",shared_object)

};



class IHAQR_to_state_mapper : public named_object {
  public:
    
    IHAQR_to_state_mapper() : named_object() { setName("IHAQR_to_state_mapper"); };
    
    template <typename SourceSpace, typename StateSpace, typename StateSpaceSystem>
    typename topology_traits< StateSpace >::point_type map_to_space(
      const IHAQR_point_type<StateSpace,StateSpaceSystem>& pt,
      const SourceSpace&, const StateSpace&) const {
      return pt.x;
    };
    
    template <typename DestSpace, typename StateSpace>
    typename topology_traits< DestSpace >::point_type map_to_space(
      const typename topology_traits< StateSpace >::point_type& pt,
      const StateSpace&, const DestSpace&) const {
      return typename topology_traits< DestSpace >::point_type(pt);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(IHAQR_to_state_mapper,0xC2400037,1,"IHAQR_to_state_mapper",named_object)
    
};




/**
 * This class implements a quaternion-topology. Because quaternions are constrained on the unit 
 * hyper-sphere, this topology is indeed bounded (yet infinite at the same time). This class
 * models the MetricSpaceConcept, the LieGroupConcept, and the PointDistributionConcept.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
class IHAQR_topology : public named_object
{
  public:
    typedef IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> self;
    
    friend class MEAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler>;
    
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixA_type matrixA_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixB_type matrixB_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixC_type matrixC_type;
    typedef typename ctrl::linear_ss_system_traits< StateSpaceSystem >::matrixD_type matrixD_type;
    
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_type state_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_difference_type state_difference_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::point_derivative_type state_derivative_type;
    typedef typename ctrl::ss_system_traits< StateSpaceSystem >::input_type system_input_type;
    
    typedef StateSpaceSampler state_sampler_type;
    
    typedef IHAQR_point_type<StateSpace, StateSpaceSystem> point_type;
    typedef typename point_type::linearization_payload linearization_payload;
    typedef typename point_type::IHAQR_payload IHAQR_payload;
    
    struct point_difference_type {
      state_difference_type dx;
      
      explicit point_difference_type(const state_difference_type& aDX = state_difference_type()) : dx(aDX) { };
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
    state_sampler_type m_get_sample;
    hyperbox_topology< system_input_type > m_input_space;
    hyperbox_topology< system_input_type > m_input_rate_space;
    
    mat<double,mat_structure::diagonal> m_R;
    mat<double,mat_structure::diagonal> m_Q;
    double m_time_step;
    double m_max_time_horizon;
    double m_goal_proximity_threshold;
    
    void compute_linearization_data(const point_type& a) const {
      a.lin_data = shared_ptr< linearization_payload >(new linearization_payload());
      
      // compute u
      matrixC_type Ctmp;
      matrixD_type Dtmp;
      vect_n<double> f = - to_vect<double>(m_system->get_state_derivative(m_space, a.x, a.lin_data->u, 0.0));
      mat_vect_adaptor< vect_n<double> > f_m(f);
      m_system->get_linear_blocks(a.lin_data->A, a.lin_data->B, Ctmp, Dtmp, m_space, 0.0, a.x, a.lin_data->u);
      vect_n<double> u = to_vect<double>(a.lin_data->u);
      mat_vect_adaptor< vect_n<double> > u_m(u);
      linlsq_QR(a.lin_data->B, u_m, f_m);
      a.lin_data->u = from_vect< system_input_type >(u);
      
      // fill in A and B
      m_system->get_linear_blocks(a.lin_data->A, a.lin_data->B, Ctmp, Dtmp, m_space, 0.0, a.x, a.lin_data->u);
      
      a.lin_data->c = m_system->get_state_derivative(m_space, a.x, a.lin_data->u, 0.0);
      
    };
    
    void compute_IHAQR_data(const point_type& a) const {
      if(!a.lin_data)
        compute_linearization_data(a);
      if(a.IHAQR_data)
        return;
      a.IHAQR_data = shared_ptr< IHAQR_payload >(new IHAQR_payload());
      
      // compute c
      vect_n<double> c_v = to_vect<double>(a.lin_data->c);
      
      // solve for M, K, and u_bias
      try {
        vect_n<double> u_bias_v = to_vect<double>(a.lin_data->u);
//         std::cout << " Asys = " << a.lin_data->A << std::endl;
//         std::cout << " Bsys = " << a.lin_data->B << std::endl;
//         std::cout << " u_sys = " << a.lin_data->u << std::endl;
//         solve_IHCT_AQR_with_reduction(a.lin_data->A, a.lin_data->B, c_v, m_Q, m_R, 
//                                       a.IHAQR_data->K, a.IHAQR_data->M, u_bias_v, 1e-3, true);
        solve_IHCT_AQR(a.lin_data->A, a.lin_data->B, c_v, m_Q, m_R, 
                       a.IHAQR_data->K, a.IHAQR_data->M, u_bias_v, 1e-4, true);
        a.IHAQR_data->u_bias = from_vect<system_input_type>(u_bias_v);
//         std::cout << " IHAQR Gain = " << a.IHAQR_data->K << std::endl;
//         std::cout << " IHAQR bias = " << a.IHAQR_data->u_bias << std::endl;
      } catch(std::exception& e) {
        std::cout << "Warning! Solution to the CARE problem could not be found for the given state point: " << a.x << std::endl
                  << "  The following exception was raised: " << e.what() << std::endl;
      };
    };
    
    virtual bool is_free_impl(const state_type& a) const {
      return true;
    };
    
    system_input_type get_bounded_input(const system_input_type& u_prev, system_input_type u_bias, system_input_type u_correction) const {
      
      m_input_space.bring_point_in_bounds(u_bias);
      
      system_input_type u_current = u_bias + u_correction;
      if(m_input_space.is_in_bounds(u_current)) {
        system_input_type du_dt = (u_current - u_prev) * (1.0 / m_time_step);
        m_input_rate_space.bring_point_in_bounds(du_dt);
        return u_prev + m_time_step * du_dt;
      };
      
      for(std::size_t j = 0; j < 10; ++j) {
        u_correction *= 0.5;
        u_current -= u_correction;
        if(m_input_space.is_in_bounds(u_current)) {
          u_bias = u_current;
          u_current += u_correction;
        };
      };
      
      system_input_type du_dt = (u_bias - u_prev) * (1.0 / m_time_step);
      m_input_rate_space.bring_point_in_bounds(du_dt);
      return u_prev + m_time_step * du_dt;
    };
    
    point_type move_position_toward_impl(const point_type& a, double fraction, const point_type& b, bool with_collision_check) const 
    {
      if(!a.IHAQR_data)
        compute_IHAQR_data(a);
      if(!b.IHAQR_data)
        compute_IHAQR_data(b);
      state_type goal_point = m_space.move_position_toward(a.x, fraction, b.x);
      state_type x_current = a.x;
      state_type x_next = x_current;
      system_input_type u_prev = a.lin_data->u - a.IHAQR_data->u_bias;
      m_input_space.bring_point_in_bounds(u_prev);
      
      // while not reached (fly-by) the goal yet:
      double current_time = 0.0;
      while( ( current_time < m_max_time_horizon ) &&
             ( m_space.distance(x_current, goal_point) > m_goal_proximity_threshold ) ) {
        // compute the current IHAQR input
//         std::cout << " t = " << current_time << std::endl;
//         std::cout << " x_difference = " << m_space.difference(x_current, goal_point) << std::endl;
        system_input_type u_current = get_bounded_input(u_prev, 
                                                        b.lin_data->u - b.IHAQR_data->u_bias, 
                                                        - from_vect< system_input_type >( b.IHAQR_data->K * to_vect<double>(m_space.difference(x_current, goal_point)) )
                                                       );
        
        constant_trajectory< vector_topology< system_input_type > > input_traj(u_current);
        
        // integrate for one time-step.
        ctrl::detail::runge_kutta4_integrate_impl(
          m_space,
          *m_system,
          x_current,
          x_next,
          input_traj,
          current_time,
          current_time + m_time_step,
          m_time_step * 1e-2);
//         std::cout << " x_next = " << x_next << std::endl;
        if(!with_collision_check || is_free_impl(x_next)) {
          x_current = x_next;
          current_time += m_time_step;
          u_prev = u_current;
        } else
          break;
      };
      
      point_type result( x_current );
      return result;
    };
    
    point_type random_point_impl(bool with_collision_check) const {
      state_type result_pt = m_get_sample(m_space);
      while(with_collision_check && !is_free_impl(result_pt))
        result_pt = m_get_sample(m_space);
      point_type result( result_pt );
      compute_IHAQR_data(result);
      return result;
    };
    
  public:
    
    StateSpace& get_state_space() { return m_space; };
    const StateSpace& get_state_space() const { return m_space; };
    
    state_sampler_type get_state_sampler() const { return m_get_sample; };
    
    const mat<double,mat_structure::diagonal>& get_input_cost_matrix() const { return m_R; };
    const mat<double,mat_structure::diagonal>& get_state_cost_matrix() const { return m_Q; };
    
    /**
     * Default constructor.
     */
    IHAQR_topology(const std::string& aName = "IHAQR_topology",
                   const shared_ptr< StateSpaceSystem >& aSystem = shared_ptr< StateSpaceSystem >(),
                   const StateSpace& aSpace = StateSpace(),
                   const system_input_type& aMinInput = system_input_type(),
                   const system_input_type& aMaxInput = system_input_type(),
                   const system_input_type& aInputBandwidth = system_input_type(),
                   const mat<double,mat_structure::diagonal>& aR = (mat<double,mat_structure::diagonal>()),
                   const mat<double,mat_structure::diagonal>& aQ = (mat<double,mat_structure::diagonal>()),
                   double aTimeStep = 0.1,
                   double aMaxTimeHorizon = 10.0,
                   double aGoalProximityThreshold = 1.0,
                   state_sampler_type aGetSample = state_sampler_type()) : 
                   named_object(),
                   m_system(aSystem),
                   m_space(aSpace),
                   m_get_sample(aGetSample),
                   m_input_space(aName + "_input_space", aMinInput, aMaxInput),
                   m_input_rate_space(aName + "_input_rate_space", -aInputBandwidth, aInputBandwidth),
                   m_R(aR),
                   m_Q(aQ),
                   m_time_step(aTimeStep),
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
      if(!b.IHAQR_data)
        compute_IHAQR_data(b);
      vect_n<double> dx =  to_vect<double>(m_space.difference(b.x, a.x));
      return dx * b.IHAQR_data->M * dx;
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
      return point_difference_type( m_space.difference(b.x, a.x) );
    };

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return point_type( m_space.adjust(a.x, delta.dx) );
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return point_type( m_space.origin() );
    };
      
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      return m_space.is_in_bounds(a.x);
    };
    
    // NOTE: don't know if I can get rid of this. (only seems useful in bounded interpolators (and samplers)).
    void bring_point_in_bounds(point_type& p) const {
      if(!m_space.is_in_bounds(p.x)) {
        m_space.bring_point_in_bounds(p.x);
        p.lin_data = shared_ptr< linearization_payload >();
        p.IHAQR_data = shared_ptr< IHAQR_payload >();
      };
    };
    
    // NOTE: don't know if I can get rid of this. (only seems useful in bounded interpolators (and samplers)).
    point_difference_type get_diff_to_boundary(const point_type& p) const {
      return point_difference_type( m_space.get_diff_to_boundary(p.x) );
    };

    /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return move_position_toward_impl(a,fraction,b,false);
    };

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_system)
        & RK_SERIAL_SAVE_WITH_NAME(m_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_get_sample)
        & RK_SERIAL_SAVE_WITH_NAME(m_input_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_input_rate_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_R)
        & RK_SERIAL_SAVE_WITH_NAME(m_Q)
        & RK_SERIAL_SAVE_WITH_NAME(m_time_step)
        & RK_SERIAL_SAVE_WITH_NAME(m_max_time_horizon)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_proximity_threshold);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_system)
        & RK_SERIAL_LOAD_WITH_NAME(m_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_get_sample)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_rate_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_R)
        & RK_SERIAL_LOAD_WITH_NAME(m_Q)
        & RK_SERIAL_LOAD_WITH_NAME(m_time_step)
        & RK_SERIAL_LOAD_WITH_NAME(m_max_time_horizon)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_proximity_threshold);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400032,1,"IHAQR_topology",named_object)
    
};

template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_metric_space< IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> > : boost::mpl::true_ { };
        
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_point_distribution< IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> > : 
  is_point_distribution<StateSpace> { };

template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_metric_symmetric< IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> > : boost::mpl::false_ { };
  
  

  
/**
 * This class has collision detection.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
class IHAQR_topology_with_CD : public IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> {
  public:
    typedef IHAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> self;
    typedef IHAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> base_type;
    
    typedef base_type super_space_type;
    
    typedef typename base_type::state_type state_type;
    typedef typename base_type::state_difference_type state_difference_type;
    typedef typename base_type::state_derivative_type state_derivative_type;
    typedef typename base_type::system_input_type system_input_type;
    
    typedef typename base_type::linearization_payload linearization_payload;
    typedef typename base_type::IHAQR_payload IHAQR_payload;
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  protected:
    shared_ptr< kte::direct_kinematics_model > m_model; 
    
    virtual bool is_free_impl(const state_type& a) const {
      
      detail::write_joint_coordinates_impl(a, this->m_space, m_model);
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
    IHAQR_topology_with_CD(const std::string& aName = "IHAQR_topology_with_CD",
                           const shared_ptr< StateSpaceSystem >& aSystem = shared_ptr< StateSpaceSystem >(),
                           const StateSpace& aSpace = StateSpace(),
                           const system_input_type& aMinInput = system_input_type(),
                           const system_input_type& aMaxInput = system_input_type(),
                           const system_input_type& aInputBandwidth = system_input_type(),
                           const mat<double,mat_structure::diagonal>& aR = (mat<double,mat_structure::diagonal>()),
                           const mat<double,mat_structure::diagonal>& aQ = (mat<double,mat_structure::diagonal>()),
                           double aTimeStep = 0.1,
                           double aMaxTimeHorizon = 5.0,
                           double aGoalProximityThreshold = 1.0,
                           const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >()) :
                           base_type(aName, aSystem, aSpace, aMinInput, aMaxInput, aInputBandwidth, aR, aQ, aTimeStep,aMaxTimeHorizon, aGoalProximityThreshold) { };
    
    /**
     * Default constructor.
     */
    IHAQR_topology_with_CD(const base_type& aBaseSpace,
                           const shared_ptr< kte::direct_kinematics_model >& aModel) : 
                           base_type(aBaseSpace),
                           m_model(aModel),
                           m_proxy_env_2D(),
                           m_proxy_env_3D() { };
    
    virtual ~IHAQR_topology_with_CD() { };
    
    super_space_type& get_super_space() { return *this; };
    
    const super_space_type& get_super_space() const { return *this; };
    
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const {
      point_type result = this->move_position_toward_impl(a, 1.0, b, true);
      if(get(distance_metric,this->get_state_space())(a.x, b.x, this->get_state_space()) * 0.1 > get(distance_metric,this->get_state_space())(result.x, b.x, this->get_state_space()) ) {
        return base_type::distance(a,b);
      } else {
        return std::numeric_limits<double>::infinity();
      };
    };
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      return this->random_point_impl(true);
    };
    
    
    bool is_free(const point_type& a) const {
      return this->is_free_impl(a.x);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const {
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400033,1,"IHAQR_topology_with_CD",base_type)
    
    
};
  
  
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_metric_space< IHAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > : boost::mpl::true_ { };
        
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_point_distribution< IHAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > : 
  is_point_distribution<StateSpace> { };
  
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct is_metric_symmetric< IHAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > : boost::mpl::false_ { };
  
  



};

};


#endif








