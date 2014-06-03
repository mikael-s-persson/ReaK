/**
 * \file hamming_iter_mod_integrator_sys.hpp
 * 
 * This library implements the Iterated Modified Hamming Method for the numerical integration of a 
 * state-space continuous-time system. Hamming methods are a type of multi-step predictor-corrector 
 * algorithm for numerical integration in which the prediction step is implicitly merged with the 
 * correction step (except for the first steps). These methods show good stability margins even for 
 * high-order stiff problems. These methods are 3rd order method, but show stability for up to 5th or 6th 
 * order systems (and errors are also of that order, roughly speaking).
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

#ifndef REAK_HAMMING_ITER_MOD_INTEGRATOR_SYS_HPP
#define REAK_HAMMING_ITER_MOD_INTEGRATOR_SYS_HPP

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/ctrl/ctrl_sys/state_space_sys_concept.hpp>
#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/temporal_space_concept.hpp>
#include <ReaK/ctrl/interpolation/spatial_trajectory_concept.hpp>

#include <ReaK/core/integrators/integration_exceptions.hpp>
#include <ReaK/core/lin_alg/vect_alg.hpp>
#include <ReaK/core/lin_alg/arithmetic_tuple.hpp>

namespace ReaK {

namespace ctrl {


namespace detail {
  
  
  /*
 P_n+1 = Y_n-3 + dt * 4/3 * (2*Yp_n - Yp_n-1 + 2*Yp_n-2)
 M = P_n+1 - 112/121 * (P_n - C)
 C = 1/8 * (9*Y_n - Y_n-2 + dt * 3 * (Mp + 2*Yp_n - Yp_n-1))
 Y_n+1 = C + 9/121 * (P_n+1 - C)
*/
  
  template <typename StateSpace, 
            typename StateSpaceSystem,
            typename InputTrajectory> 
  void hamming_iter_mod_integrate_impl(
      const StateSpace& space,
      const StateSpaceSystem& sys,
      const typename pp::topology_traits<StateSpace>::point_type& start_point,
      typename pp::topology_traits<StateSpace>::point_type& end_point,
      const InputTrajectory& u_traj,
      double start_time,
      double end_time,
      double time_step,
      double tolerance,
      std::size_t max_iter) {
    if ((time_step == 0.0) ||
        ((time_step > 0.0) && (start_time > end_time)) ||
        ((time_step < 0.0) && (start_time < end_time)))
      throw impossible_integration(start_time, end_time, time_step);
    
    using std::sqrt;
    using ReaK::to_vect;
    typedef typename pp::topology_traits<StateSpace>::point_type PointType;
    typedef typename pp::topology_traits<StateSpace>::point_difference_type PointDiffType;
    
    typedef typename pp::spatial_trajectory_traits<InputTrajectory>::const_waypoint_descriptor InputWaypoint;
    typedef typename pp::spatial_trajectory_traits<InputTrajectory>::point_type InputType;
    std::pair< InputWaypoint, InputType> u_wp = u_traj.get_waypoint_at_time(start_time);
    
    double back_time = start_time;
    PointDiffType dp = sys.get_state_derivative(space, start_point, u_wp.second.pt, start_time);
    PointType w = start_point;
    PointType Y_n[3];
    PointDiffType Yp_n[3];
    
    for(std::size_t i = 0; i < 3; ++i) {
      PointDiffType k1 = time_step * dp;
      Y_n[i] = space.adjust(w, -0.25 * k1);
      
      back_time -= time_step * 0.25;
      u_wp = u_traj.move_time_diff_from(u_wp, -0.25 * time_step);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      PointDiffType k2 = time_step * Yp_n_1;
      Y_n[i] = space.adjust(Y_n[i], (5.0 / 32.0) * k1 - (9.0 / 32.0) * k2);
      
      back_time -= time_step * 0.125;
      u_wp = u_traj.move_time_diff_from(u_wp, -0.125 * time_step);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      PointDiffType k3 = time_step * Yp_n[i];
      Y_n[i] = space.adjust(Y_n[i], (1250865.0 / 351520.0) * k2 - (276165.0 / 351520.0) * k1 - (1167360.0 / 351520.0) * k3);
      
      back_time -= 57.0 * time_step / 104.0;
      u_wp = u_traj.move_time_diff_from(u_wp, -57.0 / 104.0 * time_step);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      PointDiffType k4 = time_step * Yp_n[i];
      Y_n[i] = space.adjust(w, 8.0 * k2 - (439.0 / 216.0) * k1 - (3680.0 / 513.0) * k3 + (845.0 / 4104.0) * k4);
      
      back_time -= time_step / 13.0;
      u_wp = u_traj.move_time_diff_from(u_wp, -time_step / 13.0);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      PointDiffType k5 = time_step * Yp_n[i];
      Y_n[i] = space.adjust(w, (8.0 / 27.0) * k1 - 2.0 * k2 + (3544.0 / 2565.0) * k3 - (1859.0 / 4104.0) * k4 + (11.0 / 40.0) * k5);
      
      back_time += time_step * 0.5;
      u_wp = u_traj.move_time_diff_from(u_wp, 0.5 * time_step);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      Y_n[i] = space.adjust(w, (-16.0 / 135.0) * k1 - (6656.0 / 12825.0) * k3 - (28561.0 / 56430.0) * k4 + (9.0 / 50.0) * k5 - (2.0 * time_step / 55.0) * Yp_n[i]);
      
      back_time -= time_step * 0.5;
      u_wp = u_traj.move_time_diff_from(u_wp, -0.5 * time_step);
      Yp_n[i] = sys.get_state_derivative(space, Y_n[i], u_wp.second.pt, back_time);
      
      w = Y_n[i];
      dp = Yp_n[i];
    };
    
    PointType C = start_point;
    PointType P_n = start_point; // this makes sure the initial error sweep is zero.
    
    double t = start_time;
    end_point = start_point;
    u_wp = u_traj.get_waypoint_at_time(t);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    
    while(((time_step > 0.0) && (t < end_time)) || 
          ((time_step < 0.0) && (t > end_time))) {

      t += time_step;
      u_wp = u_traj.move_time_diff_from(u_wp, time_step);
      
      PointType P = space.adjust(Y_n[2], (time_step * 4.0 / 3.0) * (2.0 * (dp + Yp_n[1]) - Yp_n[0]));
      PointType prevM = space.adjust(P, (112.0 / 121.0) * space.difference(C, P_n));
      
      double ErrorEst = 2.0 * tolerance;
      std::size_t i = 0;
      while((ErrorEst > tolerance) && (i < max_iter)) {
        PointDiffType Mp = sys.get_state_derivative(space, prevM, u_wp.second.pt, t);
        C = space.adjust(end_point, 0.125 * space.difference(end_point, Y_n[1]) + (time_step * 3.0 / 8.0) * (Mp + 2.0 * dp - Yp_n[0]));
        
        PointType M = space.adjust(C, (9.0 / 121.0) * space.difference(P, C));
        vect_n<double> err_vect = to_vect<double>(space.difference(M, prevM));
        ErrorEst = 0.0;
        for(std::size_t j = 0; j < err_vect.size(); ++j)
          ErrorEst += err_vect[j] * err_vect[j];
        ErrorEst = sqrt(ErrorEst);
        prevM = M;
        ++i;
      };
      Y_n[2] = Y_n[1];
      Y_n[1] = Y_n[0];
      Y_n[0] = end_point;
      end_point = prevM;
      
      Yp_n[1] = Yp_n[0];
      Yp_n[0] = dp;
      dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    };
    
  };
  
};


/**
 * This class is a factory for integrators that use the Iterated Modified Hamming Method. Hamming 
 * methods are a type of multi-step predictor-corrector algorithm for numerical integration in which 
 * the prediction step is implicitly merged with the correction step (except for the first steps). 
 * These methods show good stability margins even for high-order stiff problems. These methods are 
 * 3rd order method, but show stability for up to 5th or 6th order systems (and errors are also of 
 * that order, roughly speaking).
 * \tparam TemporalSpace The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations), see SSSystemConcept.
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times, see pp::SpatialTrajectoryConcept.
 */
template <typename TemporalSpace, 
          typename StateSpaceSystem,
          typename InputTrajectory>
class hamming_iter_mod_integrator_factory : public named_object {
  public:
    typedef hamming_iter_mod_integrator_factory<TemporalSpace,StateSpaceSystem,InputTrajectory> self;
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
    double m_tolerance;
    std::size_t m_max_iter;
    
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
          detail::hamming_mod_integrate_impl(
            m_parent->m_t_space->get_space_topology(),
            *(m_parent->m_sys),
            m_start_point->pt,
            end_point.pt,
            *(m_parent->m_input_traj),
            m_start_point->time,
            end_point.time,
            m_parent->m_time_step,
            m_parent->m_tolerance,
            m_parent->m_max_iter
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
     * \param aTolerance The tolerance on the norm of the prediction-correction error to stop the iterations.
     * \param aMaxIter The maximum number of correction iterations.
     */
    hamming_iter_mod_integrator_factory(
      const std::string& aName = "", 
      const shared_ptr< const TemporalSpace >& aTSpace = shared_ptr< const TemporalSpace >(),
      const shared_ptr< const StateSpaceSystem >& aSystem = shared_ptr< const StateSpaceSystem >(),
      const shared_ptr< const InputTrajectory >& aInputTraj = shared_ptr< const InputTrajectory >(),
      double aTimeStep = 1e-3,
      double aTolerance = 1e-6,
      std::size_t aMaxIter = 10) :
      named_object(),
      m_t_space(aTSpace),
      m_sys(aSystem),
      m_input_traj(aInputTraj),
      m_time_step(aTimeStep),
      m_tolerance(aTolerance),
      m_max_iter(aMaxIter) { };
    
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
        & RK_SERIAL_SAVE_WITH_NAME(m_time_step)
        & RK_SERIAL_SAVE_WITH_NAME(m_tolerance)
        & RK_SERIAL_SAVE_WITH_NAME(m_max_iter);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_t_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_sys)
        & RK_SERIAL_LOAD_WITH_NAME(m_input_traj)
        & RK_SERIAL_LOAD_WITH_NAME(m_time_step)
        & RK_SERIAL_LOAD_WITH_NAME(m_tolerance)
        & RK_SERIAL_LOAD_WITH_NAME(m_max_iter);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC243100A,1,"hamming_iter_mod_integrator_factory",named_object)
    
};



};


};

#endif



