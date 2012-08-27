/**
 * \file euler_integrator_sys.hpp
 * 
 * This library implements an integrator that uses the Forward-Euler Method. This
 * method is a single-step fixed-step algorithm for numerical integration based on a 1st order
 * finite differencing scheme (forward-Euler rule). This is the simplest of all integrators, 
 * it is extremely cheap computationally, but it is also very inaccurate and guaranteed to rapidly 
 * divergence for all but the most trivial systems (which is the case for most explicit single-step methods).
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date August 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_EULER_INTEGRATOR_SYS_HPP
#define REAK_EULER_INTEGRATOR_SYS_HPP

#include "ctrl_sys/state_space_sys_concept.hpp"
#include "path_planning/metric_space_concept.hpp"
#include "path_planning/temporal_space_concept.hpp"
#include "path_planning/spatial_trajectory_concept.hpp"
#include "base/named_object.hpp"

namespace ReaK {

namespace ctrl {


namespace detail {
  
  template <typename StateSpace, 
            typename StateSpaceSystem,
            typename InputTrajectory> 
  void euler_integrate_impl(
      const StateSpace& space,
      const StateSpaceSystem& sys,
      const typename pp::topology_traits<StateSpace>::point_type& start_point,
      typename pp::topology_traits<StateSpace>::point_type& end_point,
      const InputTrajectory& u_traj,
      double start_time,
      double end_time,
      double time_step) {
    typedef typename pp::topology_traits<StateSpace>::point_type PointType;
    typedef typename pp::topology_traits<StateSpace>::point_difference_type PointDiffType;
    
    typedef typename pp::spatial_trajectory_traits<InputTrajectory>::const_waypoint_descriptor InputWaypoint;
    typedef typename pp::spatial_trajectory_traits<InputTrajectory>::point_type InputType;
    std::pair< InputWaypoint, InputType> u_wp = u_traj.get_waypoint_at_time(start_time);
    
    PointDiffType dp = sys.get_state_derivative(space, start_point, u_wp.second.pt, start_time);
    double t = start_time;
    end_point = start_point;
    
    while(((time_step > 0.0) && (t < end_time)) || 
          ((time_step < 0.0) && (t > end_time))) {
      end_point = space.adjust(end_point, time_step * dp);
      t += time_step;
      u_wp = u_traj.move_time_diff_from(u_wp, time_step);
      dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    };
  };
  
};


/**
 * This class is a factory for integrators that use the Forward-Euler Method. This
 * method is a single-step fixed-step algorithm for numerical integration based on a 1st order
 * finite differencing scheme (forward-Euler rule). This is the simplest of all integrators, 
 * it is extremely cheap computationally, but it is also very inaccurate and guaranteed to rapidly 
 * divergence for all but the most trivial systems (which is the case for most explicit single-step methods).
 * \tparam TemporalSpace The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations), see SSSystemConcept.
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times, see pp::SpatialTrajectoryConcept.
 */
template <typename TemporalSpace, 
          typename StateSpaceSystem,
          typename InputTrajectory>
class euler_integrator_factory : public named_object {
  public:
    typedef euler_integrator_factory<TemporalSpace,StateSpaceSystem,InputTrajectory> self;
    typedef TemporalSpace topology;
    typedef typename pp::topology_traits< TemporalSpace >::point_type point_type;
    
    typedef typename pp::temporal_space_traits< TemporalSpace >::space_topology space_topology;
    typedef typename pp::temporal_space_traits< TemporalSpace >::time_topology time_topology;
    
    BOOST_CONCEPT_ASSERT((pp::TemporalSpaceConcept<TemporalSpace>));
    BOOST_CONCEPT_ASSERT((SSSystemConcept<StateSpaceSystem, space_topology >));
    
  private:
    shared_ptr< const TemporalSpace > m_t_space;
    shared_ptr< const StateSpaceSystem > m_sys;
    shared_ptr< const InputTrajectory > m_input_traj;
    double m_time_step;
    
  public:
    
    /** 
     * This class represents an integration task, from a starting temporal point.
     */
    class extrapolator_type {
      private:
        const self* m_parent;
        const point_type* m_start_point;
        
      public:
        extrapolator_type(const self* aParent, 
                          const point_type* aStartPoint) :
                          m_parent(aParent), m_start_point(aStartPoint) { };
        
        /**
         * Sets the starting point (by a raw-pointer) of the interpolation.
         * \param aStartPoint A raw-pointer to the starting point of the interpolation.
         */ 
        void set_start_point(const point_type* aStartPoint) {
          m_start_point = aStartPoint;
        };
        
        /**
         * Returns the pointer to the starting point of the interpolation.
         * \return The pointer to the starting point of the interpolation.
         */
        const point_type* get_start_point() const { return m_start_point; };
        
        /**
         * Returns the point integrated up to the given time (within the time-step precision).
         * \param end_time The end of the integration period.
         * \return The point resulting from integration from the starting point up to the given end-time.
         */
        point_type get_point_at_time(double end_time) const {
          point_type end_point;
          end_point.time = end_time;
          detail::euler_integrate_impl(
            m_parent->m_t_space->get_space_topology(),
            *(m_parent->m_sys),
            m_start_point->pt,
            end_point.pt,
            *(m_parent->m_input_traj),
            m_start_point->time,
            end_point.time,
            m_parent->m_time_step
          );
          return end_point;
        };
        
    };
    
    /**
     * Parametrized Constructor.
     * \param aName The name of the integrator factory object.
     * \param aTSpace A pointer to the temporal space to use.
     * \param aSystem A pointer to the state-space system to be integrated.
     * \param aInputTraj A pointer to the trajectory object which can deliver input vectors at any given time point.
     * \param aTimeStep The integration time-step to use.
     */
    euler_integrator_factory(
      const std::string& aName = "", 
      const shared_ptr< const TemporalSpace >& aTSpace = shared_ptr< const TemporalSpace >(),
      const shared_ptr< const StateSpaceSystem >& aSystem = shared_ptr< const StateSpaceSystem >(),
      const shared_ptr< const InputTrajectory >& aInputTraj = shared_ptr< const InputTrajectory >(),
      double aTimeStep = 1e-3) :
      named_object(),
      m_t_space(aTSpace),
      m_sys(aSystem),
      m_input_traj(aInputTraj),
      m_time_step(aTimeStep) { };
    
    /**
     * Sets the pointer to the temporal space used by this integrator factory.
     * \param aTSpace A pointer to the temporal space used by this integrator factory.
     */
    void set_temporal_space(const shared_ptr< const TemporalSpace >& aTSpace) {
      m_t_space = aTSpace;
    };
    /**
     * Returns a pointer to the temporal space used by this integrator factory.
     * \return A pointer to the temporal space used by this integrator factory.
     */
    const shared_ptr< const TemporalSpace >& get_temporal_space() const { return m_t_space; };
    
    /**
     * Sets the pointer to the state-space system used by this integrator factory.
     * \param aSystem A pointer to the state-space system used by this integrator factory.
     */
    void set_system(const shared_ptr< const StateSpaceSystem >& aSystem) {
      m_sys = aSystem;
    };
    /**
     * Returns a pointer to the state-space system used by this integrator factory.
     * \return A pointer to the state-space system used by this integrator factory.
     */
    const shared_ptr< const StateSpaceSystem >& get_system() const { return m_sys; };
    
    /**
     * Sets the pointer to the input-vector trajectory used by this integrator factory.
     * \param aInputTraj A pointer to the input-vector trajectory used by this integrator factory.
     */
    void set_input_trajectory(const shared_ptr< const InputTrajectory >& aInputTraj) {
      m_input_traj = aInputTraj;
    };
    /**
     * Returns a pointer to the input-vector trajectory used by this integrator factory.
     * \return A pointer to the input-vector trajectory used by this integrator factory.
     */
    const shared_ptr< const InputTrajectory >& get_input_trajectory() const { return m_input_traj; };
    
    /**
     * Creates an "extrapolator" from a given starting point (by pointer). This creates 
     * an IVP integrator (which is a special kind of extrapolator).
     * \param aStartPoint A raw-pointer to a starting point (raw pointers are used as starting points are usually stored in a waypoint trajectory).
     * \return An IVP integrator from the given starting point.
     */
    extrapolator_type create_extrapolator(const point_type* aStartPoint) const {
      return extrapolator_type(this, aStartPoint);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_t_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_sys)
        & RK_SERIAL_SAVE_WITH_NAME(m_input_traj)
        & RK_SERIAL_SAVE_WITH_NAME(m_time_step);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_t_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_sys)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_traj)
        & RK_SERIAL_LOAD_WITH_NAME(m_time_step);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2431001,1,"euler_integrator_factory",named_object)
};



};


};

#endif



