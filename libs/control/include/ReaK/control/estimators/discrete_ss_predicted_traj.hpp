/**
 * \file discrete_ss_predicted_traj.hpp
 *
 * This library provides a class template which can generate a predicted trajectory of
 * belief-states. This class template relies on several classes to implement its functionality.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2013
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

#ifndef REAK_DISCRETE_SS_PREDICTED_TRAJ_HPP
#define REAK_DISCRETE_SS_PREDICTED_TRAJ_HPP

#include <ReaK/topologies/interpolation/predicted_trajectory_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include <ReaK/topologies/spaces/temporal_space.hpp>
#include <ReaK/topologies/spaces/time_poisson_topology.hpp>

#include <ReaK/topologies/interpolation/waypoint_container.hpp>

#include <ReaK/control/systems/discrete_sss_concept.hpp>

namespace ReaK::ctrl {

/**
 * This class template can generate a predicted trajectory of the states of a discrete-time
 * state-space system. This class template relies on several classes to implement its functionality.
 * Given a space (or topology) whose point-type is a state, and a state-space system, this
 * class predict the states at equal time intervals from a starting state, and can provide
 * state predictions at any given time. The predictions are populated in a Just-In-Time (JIT)
 * fashion (i.e., lazily), and internally maintains a current horizon (time to last pre-computed
 * prediction) which can be explicitely extended or pruned (eliminating predictions beyond a certain time).
 *
 * Models: SpatialTrajectoryConcept, PredictedTrajectoryConcept.
 *
 * \tparam StateSpaceType The topology type in which the state vectors can reside.
 * \tparam DiscreteSSSystem The discrete-time state-space system used to make the predictions.
 * \tparam InputTrajectory The input vector trajectory to provide input vectors at any given time, should model the
 *SpatialTrajectoryConcept over a vector-topology of input vectors.
 */
template <typename StateSpaceType, typename DiscreteSSSystem,
          typename InputTrajectory>
class discrete_ss_predicted_traj
    : public pp::waypoint_container<pp::temporal_space<
          StateSpaceType, pp::time_poisson_topology, pp::time_distance_only>> {
 public:
  using self = discrete_ss_predicted_traj<StateSpaceType, DiscreteSSSystem,
                                          InputTrajectory>;
  using topology = pp::temporal_space<StateSpaceType, pp::time_poisson_topology,
                                      pp::time_distance_only>;
  using base_class_type = pp::waypoint_container<topology>;

  using topology_ptr = std::shared_ptr<topology>;

  using waypoint_descriptor = typename base_class_type::waypoint_descriptor;
  using const_waypoint_descriptor =
      typename base_class_type::const_waypoint_descriptor;
  using const_waypoint_bounds = typename base_class_type::const_waypoint_bounds;
  using distance_metric = typename base_class_type::distance_metric;

  using waypoint_pair = typename base_class_type::waypoint_pair;

  using point_type = typename pp::topology_traits<topology>::point_type;
  using point_difference_type =
      typename pp::topology_traits<topology>::point_difference_type;

  using time_topology =
      typename pp::temporal_space_traits<topology>::time_topology;
  using time_type = typename pp::topology_traits<time_topology>::point_type;
  using time_difference_type =
      typename pp::topology_traits<time_topology>::point_difference_type;

  using space_topology =
      typename pp::temporal_space_traits<topology>::space_topology;
  using belief_state = typename pp::topology_traits<space_topology>::point_type;
  using belief_state_diff =
      typename pp::topology_traits<space_topology>::point_difference_type;

  using state_space_system = DiscreteSSSystem;
  using state_space_system_ptr = std::shared_ptr<state_space_system>;

  using input_type =
      typename discrete_sss_traits<state_space_system>::input_type;
  using output_type =
      typename discrete_sss_traits<state_space_system>::output_type;

  BOOST_CONCEPT_ASSERT((pp::TemporalSpaceConcept<topology>));
  BOOST_CONCEPT_ASSERT(
      (DiscreteSSSConcept<state_space_system, space_topology>));

 protected:
  InputTrajectory input;
  state_space_system_ptr dt_system;

  waypoint_descriptor updated_end;

  virtual double travel_distance_impl(const point_type& a,
                                      const point_type& b) const {
    using std::abs;
    return abs(b.time - a.time);
  }

  waypoint_pair get_point_at_time_impl(double t,
                                       const_waypoint_bounds wpb_a) const {
    if (t > get_current_horizon()) {
      set_minimal_horizon(t);
      wpb_a.second = updated_end;
      --wpb_a.second;
      wpb_a.first = wpb_a.second;
      if (wpb_a.second->first > t)
        --wpb_a.first;
    }

    if ((wpb_a.first == wpb_a.second) ||
        ((t - wpb_a.first->first) <= (wpb_a.second->first - t))) {
      return waypoint_pair(wpb_a.first, point_type(t, wpb_a.first->second.pt));
    } else {
      return waypoint_pair(wpb_a.first, point_type(t, wpb_a.second->second.pt));
    }
  }

  virtual waypoint_pair move_time_diff_from_impl(
      const point_type& a, const const_waypoint_bounds& wpb_a,
      double dt) const {
    if ((a.time + dt >= wpb_a.first->first) &&
        (a.time + dt <= wpb_a.second->first)) {
      return get_point_at_time_impl(a.time + dt, wpb_a);
    } else {
      return get_point_at_time_impl(a.time + dt,
                                    this->get_waypoint_bounds(a.time + dt));
    }
  }

 public:
  /**
   * Default constructor.
   */
  discrete_ss_predicted_traj()
      : base_class_type(),
        input(),
        dt_system(),
        updated_end(this->waypoints.end()) {}

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the trajectory is.
   * \param aDTSystem A shared-pointer to the discrete-time state-space system used to make the predictions.
   * \param aInitialPoint The starting point of the predicted trajectory.
   * \param aInputTrajectory The input-vector trajectory used to compute the input vectors necessary for the prediction.
   */
  explicit discrete_ss_predicted_traj(
      const topology_ptr& aSpace, const state_space_system_ptr& aDTSystem,
      const point_type& aInitialPoint,
      const InputTrajectory& aInputTrajectory = InputTrajectory())
      : base_class_type(aSpace),
        input(aInputTrajectory),
        dt_system(aDTSystem),
        updated_end(this->waypoints.end()) {
    set_initial_point(aInitialPoint);
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;

    // this little tricky piece of code is due to iterators being invalidated across the container-swap.
    time_type lhs_horiz = lhs.get_current_horizon();
    time_type rhs_horiz = rhs.get_current_horizon();
    swap(static_cast<base_class_type&>(lhs),
         static_cast<base_class_type&>(rhs));
    lhs.updated_end = lhs.waypoints.lower_bound(lhs_horiz);
    rhs.updated_end = rhs.waypoints.lower_bound(rhs_horiz);

    swap(lhs.input, rhs.input);
    swap(lhs.dt_system, rhs.dt_system);
  }

  /**
   * Copy-constructor with standard semantics. The current horizon (all predicted) in rhs is also copied or reproduced.
   */
  discrete_ss_predicted_traj(const self& rhs)
      : base_class_type(static_cast<const base_class_type&>(rhs)),
        input(rhs.input),
        dt_system(rhs.dt_system),
        updated_end(this->waypoints.end()) {
    waypoint_descriptor last_valid = this->waypoints.begin();
    for (waypoint_descriptor it = rhs.waypoints.begin(); it != rhs.updated_end;
         ++it, ++last_valid) {
      pred_segments[&(last_valid->second)] = pred_factory.create_predictor(
          this->space->get_space_topology(), &(last_valid->second.pt),
          last_valid->first, input.get_point_at_time(last_valid->first));
    }
    updated_end = last_valid;
  }

  /**
   * Move-constructor with standard semantics. Iterators are not invalidated (iterators of rhs become those of this
   * object).
   */
  discrete_ss_predicted_traj(self&& rhs) noexcept
      : base_class_type(),
        input(std::move(rhs.input)),
        dt_system(std::move(rhs.dt_system)),
        updated_end(this->waypoints.end()) {
    using std::swap;
    // must deal with the waypoint container separately due to completely retarded STL iterator invalidation guarantees:
    bool is_end = (rhs.updated_end == rhs.waypoints.end());
    swap(static_cast<base_class_type&>(*this),
         static_cast<base_class_type&>(rhs));
    if (is_end) {
      updated_end = this->waypoints.end();
    } else {
      updated_end = rhs.updated_end;
    }
    rhs.updated_end = rhs.waypoints.end();
    // NOTE: This move-constructor leave the rhs object in a destructible but crippled state (not even a default-ctor
    // state).
  }

  /**
   * Copy- or Move-and-swap assignment operator.
   * \param rhs The right-hand-side of the assignment.
   * \return A reference to this object.
   */
  self& operator=(self rhs) noexcept {
    swap(*this, rhs);
    return *this;
  }

  /**
   * Gets the predicted trajectory horizon, i.e., the time up to which predictions are currently valid.
   * \return The horizon for the predictions currently up-to-date in this trajectory.
   */
  time_type get_current_horizon() const {
    waypoint_descriptor last_valid = updated_end;
    --last_valid;
    return last_valid->first;
  }

  /**
   * Sets the minimal predicted trajectory horizon, i.e., extending predictions to reach that horizon
   * if they are not already reaching it.
   * \param t The new minimal horizon for the predictions up-to-date in this trajectory.
   */
  void set_minimal_horizon(const time_type& t) {

    waypoint_descriptor it = this->waypoints.lower_bound(t);

    waypoint_descriptor last_valid = updated_end;
    --last_valid;
    while ((updated_end != this->waypoints.end()) && (last_valid->first < t)) {
      updated_end->second.pt = dt_system->get_next_state(
          this->space->get_space_topology(), last_valid->second.pt,
          input.get_point_at_time(last_valid->first), last_valid->first);
      ++last_valid;
      ++updated_end;
    }

    if (it == this->waypoints.end()) {
      // must perform predictions from end to past time t.
      while (last_valid->first < t) {
        input_type u = input.get_point_at_time(last_valid->first);
        this->push_back(point_type(
            last_valid->first + dt_system->get_time_step(),
            dt_system->get_next_state(
                this->space->get_space_topology(), last_valid->second.pt,
                input.get_point_at_time(last_valid->first),
                last_valid->first)));
        ++last_valid;
      }
      updated_end = this->waypoints.end();
    }
  }

  /**
   * Sets the predicted trajectory horizon, i.e., pruning predictions that are
   * "too far ahead", or extending predictions to reach that horizon.
   * \note This function will always preserve the starting point (stem), which can be reset with set_initial_point.
   * \param t The new horizon for the predictions kept .
   */
  void set_exact_horizon(const time_type& t) {
    if (t > get_current_horizon()) {
      set_minimal_horizon(t);
      return;
    }

    waypoint_descriptor it = this->waypoints.lower_bound(t);
    ++it;  // <- this means that the pruning "conservative", and also, will always keep the beginning iterator.
    this->erase(it, this->waypoints.end());
    updated_end = this->waypoints.end();
  }

  /**
   * Resets the initial point of the predicted trajector with a point and an associated time.
   * This function triggers the elimination of all waypoints prior to this time and triggers
   * the recomputation of the belief predictions of all waypoints after this time.
   * \param pt The point that will become the first, initial waypoint.
   */
  void set_initial_point(const point_type& pt) {
    using std::abs;

    waypoint_descriptor it =
        this->waypoints.lower_bound(pt.time - 0.5 * dt_system->get_time_step());
    // check if the time difference is too much:
    if ((it == this->waypoints.end()) ||
        (abs(pt.time - it->first) > 1e-4 * dt_system->get_time_step())) {
      // we have to completely reset the entire container:
      this->waypoints.clear();
      this->push_back(pt);
      updated_end = this->waypoints.begin();
      ++updated_end;
      return;
    }

    // trim away the starting part of the trajectory:
    this->waypoints.erase(this->waypoints.begin(), it);
    // reset the start point and predictor, and reset the horizon:
    it->second = pt;
    updated_end = it;
    ++updated_end;
  }

  /**
   * Resets the initial point of the predicted trajector with a point and an associated time.
   * This function triggers the elimination of all waypoints prior to this time and triggers
   * the recomputation of the belief predictions of all waypoints after this time.
   * \param wpt The point that will become the first, initial waypoint.
   */
  void set_initial_point(const waypoint_pair& wpt) {
    set_initial_point(wpt.second);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int) const override {
    base_class_type::save(
        A, base_class_type::getStaticObjectType()->TypeVersion());
    std::size_t prediction_assumption = pred_assumption;
    A& RK_SERIAL_SAVE_WITH_NAME(input) & RK_SERIAL_SAVE_WITH_NAME(dt_system);
  }

  void load(serialization::iarchive& A, unsigned int) override {
    base_class_type::load(
        A, base_class_type::getStaticObjectType()->TypeVersion());
    std::size_t prediction_assumption = 0;
    A& RK_SERIAL_LOAD_WITH_NAME(input) & RK_SERIAL_LOAD_WITH_NAME(dt_system);
    updated_end = this->waypoints.begin();
    if (updated_end != this->waypoints.end())
      ++updated_end;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300015, 1, "discrete_ss_predicted_traj",
                              base_class_type)
};

}  // namespace ReaK::ctrl

#endif
