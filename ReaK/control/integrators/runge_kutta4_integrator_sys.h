/**
 * \file runge_kutta4_integrator_sys.h
 *
 * This library implements an integrator that uses the 4th Order Runge-Kutta Method. This
 * method is a single-step fixed-step algorithm for numerical integration based on a 4th order
 * finite differencing scheme (Runge-Kutta).
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

#ifndef REAK_CONTROL_INTEGRATORS_RUNGE_KUTTA4_INTEGRATOR_SYS_H_
#define REAK_CONTROL_INTEGRATORS_RUNGE_KUTTA4_INTEGRATOR_SYS_H_

#include "ReaK/control/systems/state_space_sys_concept.h"
#include "ReaK/core/base/named_object.h"
#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

namespace ReaK::ctrl {

namespace detail {

template <typename StateSpace, typename StateSpaceSystem,
          typename InputTrajectory>
void runge_kutta4_integrate_impl(
    const StateSpace& space, const StateSpaceSystem& sys,
    const typename pp::topology_traits<StateSpace>::point_type& start_point,
    typename pp::topology_traits<StateSpace>::point_type& end_point,
    const InputTrajectory& u_traj, double start_time, double end_time,
    double time_step) {
  using PointType = typename pp::topology_traits<StateSpace>::point_type;
  using PointDiffType =
      typename pp::topology_traits<StateSpace>::point_difference_type;

  using InputWaypoint = typename pp::spatial_trajectory_traits<
      InputTrajectory>::const_waypoint_descriptor;
  using InputType =
      typename pp::spatial_trajectory_traits<InputTrajectory>::point_type;
  std::pair<InputWaypoint, InputType> u_wp =
      u_traj.get_waypoint_at_time(start_time);

  PointDiffType dp =
      sys.get_state_derivative(space, start_point, u_wp.second.pt, start_time);
  double t = start_time;
  end_point = start_point;

  while (((time_step > 0.0) && (t < end_time)) ||
         ((time_step < 0.0) && (t > end_time))) {

    PointType w = end_point;
    PointDiffType k1 = time_step * dp;
    end_point = space.adjust(end_point, 0.5 * k1);

    t += time_step * 0.5;
    u_wp = u_traj.move_time_diff_from(u_wp, 0.5 * time_step);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k2 = time_step * dp;
    end_point = space.adjust(w, 0.5 * k2);

    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k3 = time_step * dp;
    end_point = space.adjust(w, k3);

    t += time_step * 0.5;
    u_wp = u_traj.move_time_diff_from(u_wp, 0.5 * time_step);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    end_point =
        space.adjust(end_point, (1.0 / 6.0) * k1 + (2.0 / 6.0) * k2 +
                                    (time_step / 6.0) * dp - (2.0 / 3.0) * k3);

    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
  }
}

}  // namespace detail

/**
 * This class is a factory for integrators that use the 4th Order Runge-Kutta Method. This
 * method is a single-step fixed-step algorithm for numerical integration based on a 4th order
 * finite differencing scheme (Runge-Kutta).
 * \tparam Space The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations).
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times.
 */
template <pp::TemporalSpace Space,
          SSSystem<typename pp::temporal_space_traits<Space>::space_topology>
              StateSpaceSystem,
          pp::SpatialTrajectory InputTrajectory>
class runge_kutta4_integrator_factory : public named_object {
 public:
  using self =
      runge_kutta4_integrator_factory<Space, StateSpaceSystem, InputTrajectory>;
  using topology = Space;
  using point_type = typename pp::topology_traits<Space>::point_type;

  using space_topology =
      typename pp::temporal_space_traits<Space>::space_topology;
  using time_topology =
      typename pp::temporal_space_traits<Space>::time_topology;

 private:
  std::shared_ptr<const Space> m_t_space;
  std::shared_ptr<const StateSpaceSystem> m_sys;
  std::shared_ptr<const InputTrajectory> m_input_traj;
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
    extrapolator_type(const self* aParent, const point_type* aStartPoint)
        : m_parent(aParent), m_start_point(aStartPoint) {}

    /**
     * Sets the starting point (by a raw-pointer) of the interpolation.
     * \param aStartPoint A raw-pointer to the starting point of the interpolation.
     */
    void set_start_point(const point_type* aStartPoint) {
      m_start_point = aStartPoint;
    }

    /**
     * Returns the pointer to the starting point of the interpolation.
     * \return The pointer to the starting point of the interpolation.
     */
    const point_type* get_start_point() const { return m_start_point; }

    /**
     * Returns the point integrated up to the given time (within the time-step precision).
     * \param end_time The end of the integration period.
     * \return The point resulting from integration from the starting point up to the given end-time.
     */
    point_type get_point_at_time(double end_time) const {
      point_type end_point;
      end_point.time = end_time;
      detail::runge_kutta4_integrate_impl(
          m_parent->m_t_space->get_space_topology(), *(m_parent->m_sys),
          m_start_point->pt, end_point.pt, *(m_parent->m_input_traj),
          m_start_point->time, end_point.time, m_parent->m_time_step);
      return end_point;
    }
  };

  /**
   * Parametrized Constructor.
   * \param aName The name of the integrator factory object.
   * \param aTSpace A pointer to the temporal space to use.
   * \param aSystem A pointer to the state-space system to be integrated.
   * \param aInputTraj A pointer to the trajectory object which can deliver input vectors at any given time point.
   * \param aTimeStep The integration time-step to use.
   */
  explicit runge_kutta4_integrator_factory(
      const std::string& aName,
      const std::shared_ptr<const Space>& aTSpace = {},
      const std::shared_ptr<const StateSpaceSystem>& aSystem = {},
      const std::shared_ptr<const InputTrajectory>& aInputTraj = {},
      double aTimeStep = 1e-3)
      : named_object(),
        m_t_space(aTSpace),
        m_sys(aSystem),
        m_input_traj(aInputTraj),
        m_time_step(aTimeStep) {}

  runge_kutta4_integrator_factory() : runge_kutta4_integrator_factory("") {}

  /**
   * Sets the pointer to the temporal space used by this integrator factory.
   * \param aTSpace A pointer to the temporal space used by this integrator factory.
   */
  void set_temporal_space(const std::shared_ptr<const Space>& aTSpace) {
    m_t_space = aTSpace;
  }
  /**
   * Returns a pointer to the temporal space used by this integrator factory.
   * \return A pointer to the temporal space used by this integrator factory.
   */
  const std::shared_ptr<const Space>& get_temporal_space() const {
    return m_t_space;
  }

  /**
   * Sets the pointer to the state-space system used by this integrator factory.
   * \param aSystem A pointer to the state-space system used by this integrator factory.
   */
  void set_system(const std::shared_ptr<const StateSpaceSystem>& aSystem) {
    m_sys = aSystem;
  }
  /**
   * Returns a pointer to the state-space system used by this integrator factory.
   * \return A pointer to the state-space system used by this integrator factory.
   */
  const std::shared_ptr<const StateSpaceSystem>& get_system() const {
    return m_sys;
  }

  /**
   * Sets the pointer to the input-vector trajectory used by this integrator factory.
   * \param aInputTraj A pointer to the input-vector trajectory used by this integrator factory.
   */
  void set_input_trajectory(
      const std::shared_ptr<const InputTrajectory>& aInputTraj) {
    m_input_traj = aInputTraj;
  }
  /**
   * Returns a pointer to the input-vector trajectory used by this integrator factory.
   * \return A pointer to the input-vector trajectory used by this integrator factory.
   */
  const std::shared_ptr<const InputTrajectory>& get_input_trajectory() const {
    return m_input_traj;
  }

  /**
   * Creates an "extrapolator" from a given starting point (by pointer). This creates
   * an IVP integrator (which is a special kind of extrapolator).
   * \param aStartPoint A raw-pointer to a starting point (raw pointers are used as starting points are usually stored
   * in a waypoint trajectory).
   * \return An IVP integrator from the given starting point.
   */
  extrapolator_type create_extrapolator(const point_type* aStartPoint) const {
    return extrapolator_type(this, aStartPoint);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(m_t_space) & RK_SERIAL_SAVE_WITH_NAME(m_sys) &
        RK_SERIAL_SAVE_WITH_NAME(m_input_traj) &
        RK_SERIAL_SAVE_WITH_NAME(m_time_step);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::save(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(m_t_space) & RK_SERIAL_LOAD_WITH_NAME(m_sys) &
        RK_SERIAL_LOAD_WITH_NAME(m_input_traj) &
        RK_SERIAL_LOAD_WITH_NAME(m_time_step);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2431003, 1,
                              "runge_kutta4_integrator_factory", named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_INTEGRATORS_RUNGE_KUTTA4_INTEGRATOR_SYS_H_
