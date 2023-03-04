/**
 * \file dormand_prince45_integrator_sys.hpp
 *
 * This library implements an integrator that uses the 4-5th Order Dormand-Prince Method. This
 * method is of the Runge-Kutta-Fehlberg type which are single-step variable-step algorithm for numerical
 * integration in which the error is estimated by the difference between two integration rules (order 4 and 5)
 * and the time-step is adapted on-the-fly to control the norm of that error estimate.
 * These variable-step integrators can be very accurate due to active control of errors, but they can also
 * lead to very long computational times for the integration when applied to stiff systems.
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

#ifndef REAK_DORMAND_PRINCE45_INTEGRATOR_SYS_HPP
#define REAK_DORMAND_PRINCE45_INTEGRATOR_SYS_HPP

#include <ReaK/control/systems/state_space_sys_concept.hpp>
#include <ReaK/core/base/named_object.hpp>
#include <ReaK/topologies/interpolation/spatial_trajectory_concept.hpp>
#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include <ReaK/math/integrators/integration_exceptions.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

namespace ReaK::ctrl {

namespace detail {

//----------CDORMANDPRINCE45-----------------------------------------------------

/*
     0     |
     1/5   |  1/5
     3/10  |  3/40           9/40
     4/5   |  44/45         -56/15          32/9
     8/9   |  19372/6561    -25360/2187     64448/6561    -212/729
     1     |  9017/3168     -355/33        -46732/5247     49/176        -5103/18656
     1     |  35/384         0              500/1113       125/192       -2187/6784      11/84
     ______________________________________________________________________________________________
     O(h^6)|  5170/57600     0              7571/16695     393/640       -92097/339200   187/2100       1/40
     O(h^5)|  35/384         0              500/1113       125/192       -2187/6784      11/84          0
*/

template <typename StateSpace, typename StateSpaceSystem,
          typename InputTrajectory>
void dormand_prince45_integrate_impl(
    const StateSpace& space, const StateSpaceSystem& sys,
    const typename pp::topology_traits<StateSpace>::point_type& start_point,
    typename pp::topology_traits<StateSpace>::point_type& end_point,
    const InputTrajectory& u_traj, double start_time, double end_time,
    double time_step, double tolerance, double min_step, double max_step) {
  using ReaK::to_vect;
  using std::abs;
  using std::pow;
  using PointType = typename pp::topology_traits<StateSpace>::point_type;
  using PointDiffType =
      typename pp::topology_traits<StateSpace>::point_difference_type;

  if ((time_step == 0.0) || ((time_step > 0.0) && (start_time > end_time)) ||
      ((time_step < 0.0) && (end_time > start_time)) || (tolerance <= 0.0) ||
      (min_step > max_step)) {
    throw impossible_integration(start_time, end_time, time_step);
  }

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

    PointType prevY = end_point;
    PointDiffType k1 = time_step * dp;
    end_point = space.adjust(end_point, 0.2 * k1);

    t += time_step / 5.0;
    u_wp = u_traj.move_time_diff_from(u_wp, time_step / 5.0);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k2 = time_step * dp;
    end_point = space.adjust(prevY, (3.0 / 40.0) * k1 + (9.0 / 40.0) * k2);

    t += time_step / 10.0;
    u_wp = u_traj.move_time_diff_from(u_wp, time_step / 10.0);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k3 = time_step * dp;
    end_point = space.adjust(
        prevY, (44.0 / 45.0) * k1 - (56.0 / 15.0) * k2 + (32.0 / 9.0) * k3);

    t += time_step / 2.0;
    u_wp = u_traj.move_time_diff_from(u_wp, time_step / 2.0);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k4 = time_step * dp;
    end_point =
        space.adjust(prevY, (19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2 +
                                (64448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4);

    t += 4.0 * time_step / 45.0;
    u_wp = u_traj.move_time_diff_from(u_wp, 4.0 * time_step / 45.0);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k5 = time_step * dp;
    end_point =
        space.adjust(prevY, (9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2 -
                                (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 -
                                (5103.0 / 18656.0) * k5);

    t += time_step / 9.0;
    u_wp = u_traj.move_time_diff_from(u_wp, time_step / 9.0);
    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k6 = time_step * dp;
    end_point =
        space.adjust(prevY, (35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3 +
                                (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5 +
                                (11.0 / 84.0) * k6);

    dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
    PointDiffType k7 = time_step * dp;

    vect_n<double> err_vect =
        to_vect<double>((5170.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3 +
                        (393.0 / 640.0) * k4 - (92097.0 / 339200.0) * k5 +
                        (187.0 / 2100.0) * k6 - (39.0 / 40.0) * k7);
    double Rmax = 0.0;
    std::size_t worst_DOF = 0;
    for (std::size_t i = 0; i < err_vect.size(); ++i) {
      double R = abs(err_vect[i] / time_step);
      if (R > Rmax) {
        Rmax = R;
        worst_DOF = i;
      }
    }

    if (Rmax > tolerance) {
      if (abs(time_step) <= min_step) {
        throw untolerable_integration(tolerance, Rmax, worst_DOF, time_step, t);
      }

      t -= time_step;
      u_wp = u_traj.move_time_diff_from(u_wp, -time_step);
      end_point = prevY;
      dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);
      double R = 0.84 * pow(tolerance / Rmax, 0.25);
      if (R < 0.1) {
        time_step *= 0.1;
      } else {
        time_step *= R;
      }
    } else {
      end_point = space.adjust(
          prevY, (5170.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3 +
                     (393.0 / 640.0) * k4 - (92097.0 / 339200.0) * k5 +
                     (187.0 / 2100.0) * k6 + (1.0 / 40.0) * k7);
      dp = sys.get_state_derivative(space, end_point, u_wp.second.pt, t);

      double R = 0.84 * pow(tolerance / Rmax, 0.25);
      if (R >= 4.0) {
        time_step *= 4.0;
      } else if (R > 1.0) {
        time_step *= R;
      }
      // if(((time_step > 0.0) && (t + time_step > end_time)) || ((time_step < 0.0) && (t + time_step < end_time)))
      // time_step = end_time - t;
    }

    if (abs(time_step) < min_step) {
      time_step *= abs(min_step / time_step);
    }
    if (abs(time_step) > max_step) {
      time_step *= abs(max_step / time_step);
    }
  }
}

}  // namespace detail

/**
 * This class is a factory for integrators that use the 4-5th Order Dormand-Prince Method. This
 * method is of the Runge-Kutta-Fehlberg type which are single-step variable-step algorithm for numerical
 * integration in which the error is estimated by the difference between two integration rules (order 4 and 5)
 * and the time-step is adapted on-the-fly to control the norm of that error estimate.
 * These variable-step integrators can be very accurate due to active control of errors, but they can also
 * lead to very long computational times for the integration when applied to stiff systems.
 * \tparam TemporalSpace The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations), see
 * SSSystemConcept.
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times, see
 * pp::SpatialTrajectoryConcept.
 */
template <typename TemporalSpace, typename StateSpaceSystem,
          typename InputTrajectory>
class dormand_prince45_integrator_factory : public named_object {
 public:
  using self =
      dormand_prince45_integrator_factory<TemporalSpace, StateSpaceSystem,
                                          InputTrajectory>;
  using topology = TemporalSpace;
  using point_type = typename pp::topology_traits<TemporalSpace>::point_type;

  using space_topology =
      typename pp::temporal_space_traits<TemporalSpace>::space_topology;
  using time_topology =
      typename pp::temporal_space_traits<TemporalSpace>::time_topology;

  BOOST_CONCEPT_ASSERT((pp::TemporalSpaceConcept<TemporalSpace>));
  BOOST_CONCEPT_ASSERT((SSSystemConcept<StateSpaceSystem, space_topology>));

 private:
  std::shared_ptr<const TemporalSpace> m_t_space;
  std::shared_ptr<const StateSpaceSystem> m_sys;
  std::shared_ptr<const InputTrajectory> m_input_traj;
  double m_time_step;
  double m_tolerance;
  double m_min_step;
  double m_max_step;

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
      detail::dormand_prince45_integrate_impl(
          m_parent->m_t_space->get_space_topology(), *(m_parent->m_sys),
          m_start_point->pt, end_point.pt, *(m_parent->m_input_traj),
          m_start_point->time, end_point.time, m_parent->m_time_step,
          m_parent->m_tolerance, m_parent->m_min_step, m_parent->m_max_step);
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
   * \param aTolerance The tolerance on the norm of the estimated integration error.
   * \param aMinStep The minimum integration time-step to use (will throw an intolerable_integration exception if
   * reached).
   * \param aMaxStep The maximum integration time-step to use.
   */
  explicit dormand_prince45_integrator_factory(
      const std::string& aName,
      const std::shared_ptr<const TemporalSpace>& aTSpace = {},
      const std::shared_ptr<const StateSpaceSystem>& aSystem = {},
      const std::shared_ptr<const InputTrajectory>& aInputTraj = {},
      double aTimeStep = 1e-3, double aTolerance = 1e-6, double aMinStep = 1e-8,
      double aMaxStep = 1e-1)
      : named_object(),
        m_t_space(aTSpace),
        m_sys(aSystem),
        m_input_traj(aInputTraj),
        m_time_step(aTimeStep),
        m_tolerance(aTolerance),
        m_min_step(aMinStep),
        m_max_step(aMaxStep) {}

  dormand_prince45_integrator_factory()
      : dormand_prince45_integrator_factory("") {}

  /**
   * Sets the pointer to the temporal space used by this integrator factory.
   * \param aTSpace A pointer to the temporal space used by this integrator factory.
   */
  void set_temporal_space(const std::shared_ptr<const TemporalSpace>& aTSpace) {
    m_t_space = aTSpace;
  }
  /**
   * Returns a pointer to the temporal space used by this integrator factory.
   * \return A pointer to the temporal space used by this integrator factory.
   */
  const std::shared_ptr<const TemporalSpace>& get_temporal_space() const {
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
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_t_space) & RK_SERIAL_SAVE_WITH_NAME(m_sys) &
        RK_SERIAL_SAVE_WITH_NAME(m_input_traj) &
        RK_SERIAL_SAVE_WITH_NAME(m_time_step) &
        RK_SERIAL_SAVE_WITH_NAME(m_tolerance) &
        RK_SERIAL_SAVE_WITH_NAME(m_min_step) &
        RK_SERIAL_SAVE_WITH_NAME(m_max_step);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_t_space) & RK_SERIAL_LOAD_WITH_NAME(m_sys) &
        RK_SERIAL_LOAD_WITH_NAME(m_input_traj) &
        RK_SERIAL_LOAD_WITH_NAME(m_time_step) &
        RK_SERIAL_LOAD_WITH_NAME(m_tolerance) &
        RK_SERIAL_LOAD_WITH_NAME(m_min_step) &
        RK_SERIAL_LOAD_WITH_NAME(m_max_step);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2431006, 1,
                              "dormand_prince45_integrator_factory",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif
