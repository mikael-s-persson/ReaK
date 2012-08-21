/**
 * \file midpoint_integrator_sys.hpp
 * 
 * 
 * 
 * 
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

#ifndef REAK_MIDPOINT_INTEGRATOR_SYS_HPP
#define REAK_MIDPOINT_INTEGRATOR_SYS_HPP

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
  void midpoint_integrate_impl(
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
      PointType mid_p = space.adjust(end_point, (0.5 * time_step) * dp);
      t += 0.5 * time_step;
      u_wp = u_traj.move_time_diff_from(u_wp, 0.5 * time_step);
      dp = sys.get_state_derivative(space, mid_p, u_wp.second.pt, t);
      
      end_point = space.adjust(end_point, time_step * dp);
      t += 0.5 * time_step;
      u_wp = u_traj.move_time_diff_from(u_wp, 0.5 * time_step);
      dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    };
  };
  
};


template <typename TemporalSpace, 
          typename StateSpaceSystem,
          typename InputTrajectory>
class midpoint_integrator_factory : public named_object {
  public:
    typedef midpoint_integrator_factory<TemporalSpace,StateSpaceSystem,InputTrajectory> self;
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
    
    class extrapolator_type {
      private:
        const self* m_parent;
        const point_type* m_start_point;
        
      public:
        extrapolator_type(const self* aParent, 
                          const point_type* aStartPoint) :
                          m_parent(aParent), m_start_point(aStartPoint) { };
        
        void set_start_point(const point_type* aStartPoint) {
          m_start_point = aStartPoint;
        };
        
        const point_type* get_start_point() const { return m_start_point; };
        
        point_type get_point_at_time(double end_time) const {
          point_type end_point;
          end_point.time = end_time;
          detail::midpoint_integrate_impl(
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
    
    midpoint_integrator_factory(
      const std::string& aName = "", 
      const shared_ptr< const TemporalSpace >& aTSpace = shared_ptr< const TemporalSpace >(),
      const shared_ptr< const StateSpaceSystem >& aSystem = shared_ptr< const StateSpaceSystem >(),
      const shared_ptr< const InputTrajectory >& aInputTraj = shared_ptr< const InputTrajectory >(),
      double aTimeStep = 1e-3) :
      named_object(),
      m_t_space(aTSpace),
      m_sys(aSystem),
      m_input_traj(aInputTraj),
      m_time_step(aTimeStep) { setName(aName); };
    
    void set_temporal_space(const shared_ptr< const TemporalSpace >& aTSpace) {
      m_t_space = aTSpace;
    };
    const shared_ptr< const TemporalSpace >& get_temporal_space() const { return m_t_space; };
    
    void set_system(const shared_ptr< const StateSpaceSystem >& aSystem) {
      m_sys = aSystem;
    };
    const shared_ptr< const StateSpaceSystem >& get_system() const { return m_sys; };
    
    void set_input_trajectory(const shared_ptr< const InputTrajectory >& aInputTraj) {
      m_input_traj = aInputTraj;
    };
    const shared_ptr< const InputTrajectory >& get_input_trajectory() const { return m_input_traj; };
    
    extrapolator_type create_extrapolator(const point_type* aStartPoint) const {
      return extrapolator_type(this, aStartPoint);
    };
    
};



};


};

#endif



