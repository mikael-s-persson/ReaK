/**
 * \file transformed_trajectory.h
 *
 * This library provides an implementation of a trajectory which is a map between an underlying
 * trajectory (on its topology) to another topology using a homeomorphism across the topologies.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_TRANSFORMED_TRAJECTORY_H_
#define REAK_TOPOLOGIES_INTERPOLATION_TRANSFORMED_TRAJECTORY_H_

#include "ReaK/core/base/shared_object.h"

#include "ReaK/topologies/interpolation/sequential_trajectory_concept.h"
#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/topological_map_concepts.h"

#include <map>

namespace ReaK::pp {

namespace detail {

template <typename ParentType, typename InputIterator>
struct transformed_point_iterator {
  using self = transformed_point_iterator<ParentType, InputIterator>;

  const ParentType* p_parent;
  InputIterator base_it;

  transformed_point_iterator(const ParentType* aParent, InputIterator aBaseIt)
      : p_parent(aParent), base_it(aBaseIt) {}

  friend self& operator+=(self& lhs, double rhs) {
    lhs.base_it += rhs;
    return lhs;
  }

  friend self& operator-=(self& lhs, double rhs) {
    lhs.base_it -= rhs;
    return lhs;
  }

  friend self operator+(self lhs, double rhs) { return (lhs += rhs); }
  friend self operator+(double lhs, self rhs) { return (rhs += lhs); }
  friend self operator-(self lhs, double rhs) { return (lhs -= rhs); }

  friend bool operator==(const self& lhs, const self& rhs) {
    return (lhs.base_it == rhs.base_it);
  }

  friend bool operator!=(const self& lhs, const self& rhs) {
    return !(lhs == rhs);
  }

  typename ParentType::point_type operator*() const {
    return p_parent->map_point_forward(*base_it);
  }
};

}  // namespace detail

/**
 * This class implements a trajectory which is a map between an underlying trajectory (on its
 * own topology) to another topology using a homeomorphism across the topologies.
 * This class models SpatialTrajectory.
 * \tparam Space The topology type on which the points and the trajectory can reside.
 * \tparam InputTrajectory The underlying trajectory, should model SpatialTrajectory.
 * \tparam Mapping The homeomorphic mapping between the spatial topology of the underlying
 *                 topology (of the InputTrajectory) and the given topology.
 */
template <TemporalSpace Space, typename InputTrajectory, typename Mapping>
requires SpatialTrajectory<InputTrajectory> class transformed_trajectory
    : public shared_object {
 public:
  using input_topology =
      typename spatial_trajectory_traits<InputTrajectory>::topology;

  static_assert(Bijection<Mapping, input_topology, Space>);

  using self = transformed_trajectory<Space, InputTrajectory, Mapping>;

  using waypoint_descriptor =
      typename spatial_trajectory_traits<InputTrajectory>::waypoint_descriptor;
  using const_waypoint_descriptor = typename spatial_trajectory_traits<
      InputTrajectory>::const_waypoint_descriptor;
  using point_type = topology_point_type_t<Space>;
  using point_difference_type = topology_point_difference_type_t<Space>;
  using topology = Space;
  using distance_metric =
      typename spatial_trajectory_traits<InputTrajectory>::distance_metric;

  using waypoint_pair = std::pair<const_waypoint_descriptor, point_type>;

  using point_time_iterator = detail::transformed_point_iterator<
      self, typename sequential_trajectory_traits<
                InputTrajectory>::point_time_iterator>;

  using point_fraction_iterator = detail::transformed_point_iterator<
      self, typename sequential_trajectory_traits<
                InputTrajectory>::point_fraction_iterator>;

  friend struct detail::transformed_point_iterator<
      self, typename sequential_trajectory_traits<
                InputTrajectory>::point_time_iterator>;

  friend struct detail::transformed_point_iterator<
      self, typename sequential_trajectory_traits<
                InputTrajectory>::point_fraction_iterator>;

 private:
  std::shared_ptr<Space> space;
  std::shared_ptr<InputTrajectory> traject;
  Mapping map;

  point_type map_point_forward(
      const typename spatial_trajectory_traits<InputTrajectory>::point_type&
          pt_in) const {
    return map.map_to_space(pt_in, traject->get_temporal_space(), *space);
  }

  typename spatial_trajectory_traits<InputTrajectory>::point_type
  map_point_backward(const point_type& pt_in) const {
    return map.map_to_space(pt_in, *space, traject->get_temporal_space());
  }

 public:
  /**
   * Returns the space on which the path resides.
   * \return The space on which the path resides.
   */
  const topology& get_temporal_space() const noexcept { return *space; }

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the trajectory is.
   * \param aTrajectory The underlying trajectory to use.
   * \param aMap The homeomorphic mapping object to use.
   */
  explicit transformed_trajectory(
      const std::shared_ptr<Space>& aSpace = {},
      const std::shared_ptr<InputTrajectory>& aTrajectory = {},
      const Mapping& aMap = Mapping())
      : space(aSpace), traject(aTrajectory), map(aMap) {}

  /**
   * Computes the travel distance between two points, if traveling along the path.
   * \param a The first point.
   * \param b The second point.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance(const point_type& a, const point_type& b) const {
    return get(ReaK::pp::distance_metric, *space)(a, b, *space);
  }

  /**
   * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
   * \param a The first waypoint-point-pair.
   * \param b The second waypoint-point-pair.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
    return get(ReaK::pp::distance_metric, *space)(a.second, b.second, *space);
  }

  /**
   * Returns the total travel-distance of the trajectory.
   * \return The total travel-distance of the trajectory.
   */
  double get_total_length() const { return traject->get_total_length(); }

  /**
   * Computes the point that is a time-difference away from a point on the trajectory.
   * \param a The point on the trajectory.
   * \param dt The time to move away from the point.
   * \return The point that is a time away from the given point.
   */
  point_type move_time_diff_from(const point_type& a, double dt) const {
    return get_point_at_time(a.time + dt);
  }

  /**
   * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
   * \param a The waypoint-point-pair on the trajectory.
   * \param dt The time to move away from the waypoint-point-pair.
   * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
   */
  waypoint_pair move_time_diff_from(const waypoint_pair& a, double dt) const {
    return get_waypoint_at_time(a.second.time + dt);
  }

  /**
   * Computes the point that is on the trajectory at the given time.
   * \param t The time at which the point is sought.
   * \return The point that is on the trajectory at the given time.
   */
  point_type get_point_at_time(double t) const {
    return map_point_forward(traject->get_point_at_time(t));
  }

  /**
   * Computes the waypoint-point pair that is on the trajectory at the given time.
   * \param t The time at which the waypoint-point pair is sought.
   * \return The waypoint-point pair that is on the trajectory at the given time.
   */
  waypoint_pair get_waypoint_at_time(double t) const {
    auto result = traject->get_waypoint_at_time(t);
    return {result.first, map_point_forward(result.second)};
  }

  /**
   * Returns the starting time of the trajectory.
   * \return The starting time of the trajectory.
   */
  double get_start_time() const { return traject->get_start_time(); }

  /**
   * Returns the end time of the trajectory.
   * \return The end time of the trajectory.
   */
  double get_end_time() const { return traject->get_end_time(); }

  /**
   * Returns the starting point of the waypoints.
   * \return The starting point of the waypoints.
   */
  const point_type& get_start_point() const {
    return map_point_forward(traject->get_start_point());
  }

  /**
   * Returns the starting waypoint-point-pair of the waypoints.
   * \return The starting waypoint-point-pair of the waypoints.
   */
  waypoint_pair get_start_waypoint() const {
    auto result = traject->get_start_waypoint();
    return {result.first, map_point_forward(result.second)};
  }

  /**
   * Returns the end point of the waypoints.
   * \return The end point of the waypoints.
   */
  const point_type& get_end_point() const {
    return map_point_forward(traject->get_end_point());
  }

  /**
   * Returns the end waypoint-point-pair of the waypoints.
   * \return The end waypoint-point-pair of the waypoints.
   */
  waypoint_pair get_end_waypoint() const {
    auto result = traject->get_end_waypoint();
    return {result.first, map_point_forward(result.second)};
  }

  /**
   * Returns the starting time-iterator along the trajectory.
   * \return The starting time-iterator along the trajectory.
   */
  point_time_iterator begin_time_travel() const {
    return {this, traject->begin_time_travel()};
  }

  /**
   * Returns the end time-iterator along the trajectory.
   * \return The end time-iterator along the trajectory.
   */
  point_time_iterator end_time_travel() const {
    return {this, traject->end_time_travel()};
  }

  /**
   * Returns the starting fraction-iterator along the trajectory.
   * \return The starting fraction-iterator along the trajectory.
   */
  point_fraction_iterator begin_fraction_travel() const {
    return {this, traject->begin_fraction_travel()};
  }

  /**
   * Returns the end fraction-iterator along the trajectory.
   * \return The end fraction-iterator along the trajectory.
   */
  point_fraction_iterator end_fraction_travel() const {
    return {this, traject->end_fraction_travel()};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(space) & RK_SERIAL_SAVE_WITH_NAME(traject) &
        RK_SERIAL_SAVE_WITH_NAME(map);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(space) & RK_SERIAL_LOAD_WITH_NAME(traject) &
        RK_SERIAL_LOAD_WITH_NAME(map);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440008, 1, "transformed_trajectory",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_TRANSFORMED_TRAJECTORY_H_
