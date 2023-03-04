/**
 * \file trajectory_base.hpp
 *
 * This library provides the base-class for trajectories within a temporal topology.
 * This is a base-class that stems the object-oriented compatibility of other temporal
 * trajectory classes. Then, this library provides a trajectory-wrapper class template
 * which makes an OOP-compatible trajectory class for a given trajectory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_TRAJECTORY_BASE_HPP
#define REAK_TRAJECTORY_BASE_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include <any>

#include "seq_trajectory_base.hpp"
#include "sequential_trajectory_concept.hpp"
#include "spatial_trajectory_concept.hpp"

#include <boost/concept_check.hpp>

#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include <boost/any.hpp>

namespace ReaK::pp {

/**
 * This class defines the OOP interface for a trajectory in a temporal topology.
 * \tparam Topology The topology type on which the points and the path can reside, should model the
 * TemporalSpaceConcept.
 */
template <typename Topology>
class trajectory_base : public seq_trajectory_base<Topology> {
 public:
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));

  using topology = Topology;
  using distance_metric =
      typename metric_space_traits<topology>::distance_metric_type;
  using point_type = topology_point_type_t<topology>;
  using point_difference_type = topology_point_difference_type_t<topology>;
  using self = trajectory_base<Topology>;
  using base_type = seq_trajectory_base<Topology>;

  using waypoint_descriptor = std::any;
  using const_waypoint_descriptor = std::any;
  using waypoint_pair = std::pair<const_waypoint_descriptor, point_type>;

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aName The name for this object.
   */
  explicit trajectory_base(const std::string& aName) : base_type(aName){};

  virtual const topology& get_temporal_space() const = 0;

  using seq_trajectory_base<Topology>::travel_distance;

  /**
   * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
   * \param a The first waypoint-point-pair.
   * \param b The second waypoint-point-pair.
   * \return The travel distance between two points if traveling along the path.
   */
  virtual double travel_distance(waypoint_pair& a, waypoint_pair& b) const = 0;

  /**
   * Computes the point that is a time-difference away from a point on the trajectory.
   * \param a The point on the trajectory.
   * \param dt The time to move away from the point.
   * \return The point that is a time away from the given point.
   */
  virtual point_type move_time_diff_from(const point_type& a,
                                         double dt) const = 0;

  /**
   * Computes the waypoint-point-pair that is a time-difference away from a waypoint-point-pair on the trajectory.
   * \param a The waypoint-point-pair on the trajectory.
   * \param dt The time to move away from the waypoint-point-pair.
   * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
   */
  virtual waypoint_pair move_time_diff_from(const waypoint_pair& a,
                                            double dt) const = 0;

  /**
   * Computes the point that is on the trajectory at the given time.
   * \param t The time at which the point is sought.
   * \return The point that is on the trajectory at the given time.
   */
  virtual point_type get_point_at_time(double t) const = 0;

  /**
   * Computes the waypoint-point pair that is on the trajectory at the given time.
   * \param t The time at which the waypoint-point pair is sought.
   * \return The waypoint-point pair that is on the trajectory at the given time.
   */
  virtual waypoint_pair get_waypoint_at_time(double t) const = 0;

  /**
   * Returns the starting time of the trajectory.
   * \return The starting time of the trajectory.
   */
  virtual double get_start_time() const = 0;

  /**
   * Returns the end time of the trajectory.
   * \return The end time of the trajectory.
   */
  virtual double get_end_time() const = 0;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2440009, 1, "trajectory_base", base_type)
};

/**
 * This class wraps a generic spatial trajectory class into an OOP interface.
 * It, itself, also models the generic SpatialTrajectoryConcept and SequentialTrajectoryConcept,
 * so this wrapper can be used for both purposes (i.e. this is a type-erasure class).
 * \tparam SpatialTrajectory The trajectory type to be wrapped.
 */
template <typename SpatialTrajectory>
class trajectory_wrapper
    : public trajectory_base<
          typename spatial_trajectory_traits<SpatialTrajectory>::topology> {
 public:
  using base_type = trajectory_base<
      typename spatial_trajectory_traits<SpatialTrajectory>::topology>;
  using seq_base_type = typename base_type::base_type;
  using self = trajectory_wrapper<SpatialTrajectory>;

  using topology = typename base_type::topology;
  using point_type = typename base_type::point_type;
  using point_difference_type = typename base_type::point_difference_type;

  BOOST_CONCEPT_ASSERT(
      (SequentialTrajectoryConcept<SpatialTrajectory, topology>));

  using time_topology = typename temporal_space_traits<topology>::time_topology;
  using space_topology =
      typename temporal_space_traits<topology>::space_topology;
  using distance_metric =
      typename spatial_trajectory_traits<SpatialTrajectory>::distance_metric;

  using wrapped_waypoint_descriptor = typename spatial_trajectory_traits<
      SpatialTrajectory>::waypoint_descriptor;
  using wrapped_const_waypoint_descriptor = typename spatial_trajectory_traits<
      SpatialTrajectory>::const_waypoint_descriptor;
  using wrapped_waypoint_pair =
      std::pair<wrapped_const_waypoint_descriptor, point_type>;

  using waypoint_descriptor = typename base_type::waypoint_descriptor;
  using const_waypoint_descriptor =
      typename base_type::const_waypoint_descriptor;
  using waypoint_pair = typename base_type::waypoint_pair;

  using wrapped_type = SpatialTrajectory;

 protected:
  SpatialTrajectory m_traj;
  mutable wrapped_waypoint_pair m_last_waypoint;

  using base_pt_time_iterator_impl =
      typename seq_base_type::point_time_iterator_impl;
  using gen_pt_time_iterator = typename sequential_trajectory_traits<
      SpatialTrajectory>::point_time_iterator;

  struct point_time_iterator_impl : public base_pt_time_iterator_impl {

    gen_pt_time_iterator base_it;

    explicit point_time_iterator_impl(gen_pt_time_iterator aBaseIt)
        : base_it(aBaseIt) {}

    ~point_time_iterator_impl() override = default;

    void move_by_time(double d) override { base_it += d; }

    bool is_equal_to(const base_pt_time_iterator_impl* rhs) const override {
      return (base_it ==
              static_cast<const point_time_iterator_impl*>(rhs)->base_it);
    }

    point_type get_point() const override { return *base_it; }

    base_pt_time_iterator_impl* clone() const override {
      return new point_time_iterator_impl(base_it);
    }
  };

  using base_pt_frac_iterator_impl =
      typename seq_base_type::point_fraction_iterator_impl;
  using gen_pt_frac_iterator = typename sequential_trajectory_traits<
      SpatialTrajectory>::point_fraction_iterator;

  struct point_fraction_iterator_impl : public base_pt_frac_iterator_impl {

    gen_pt_frac_iterator base_it;

    explicit point_fraction_iterator_impl(gen_pt_frac_iterator aBaseIt)
        : base_it(aBaseIt) {}

    ~point_fraction_iterator_impl() override = default;

    void move_by_fraction(double f) override { base_it += f; }

    bool is_equal_to(const base_pt_frac_iterator_impl* rhs) const override {
      return (base_it ==
              static_cast<const point_fraction_iterator_impl*>(rhs)->base_it);
    }

    point_type get_point() const override { return *base_it; }

    base_pt_frac_iterator_impl* clone() const override {
      return new point_fraction_iterator_impl(base_it);
    }
  };

 public:
  using point_time_iterator = typename base_type::point_time_iterator;
  using point_fraction_iterator = typename base_type::point_fraction_iterator;

  wrapped_type& get_underlying_trajectory() { return m_traj; }
  const wrapped_type& get_underlying_trajectory() const { return m_traj; }

  wrapped_type& get_wrapped_object() { return m_traj; }
  const wrapped_type& get_wrapped_object() const { return m_traj; }

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aName The name for this object.
   * \param aTraj The wrapped trajectory object to use.
   */
  explicit trajectory_wrapper(
      const std::string& aName = "",
      const SpatialTrajectory& aTraj = SpatialTrajectory())
      : base_type(aName), m_traj(aTraj), m_last_waypoint() {}

  const topology& get_temporal_space() const override {
    return m_traj.get_temporal_space();
  }

  double travel_distance(const point_type& a,
                         const point_type& b) const override {
    wrapped_waypoint_pair next_waypoint = m_last_waypoint;
    m_last_waypoint.second = a;
    next_waypoint.second = b;
    return m_traj.travel_distance(m_last_waypoint, next_waypoint);
  }

  double travel_distance(waypoint_pair& a, waypoint_pair& b) const override {
    m_last_waypoint = wrapped_waypoint_pair(
        std::any_cast<wrapped_const_waypoint_descriptor>(a.first), a.second);
    wrapped_waypoint_pair b_real(
        std::any_cast<wrapped_const_waypoint_descriptor>(b.first), b.second);
    double result = m_traj.travel_distance(m_last_waypoint, b_real);
    a.first = m_last_waypoint.first;
    b.first = b_real.first;
    return result;
  }

  point_type move_time_diff_from(const point_type& a,
                                 double dt) const override {
    m_last_waypoint.second = a;
    m_last_waypoint = m_traj.move_time_diff_from(m_last_waypoint, dt);
    return m_last_waypoint.second;
  }

  waypoint_pair move_time_diff_from(const waypoint_pair& a,
                                    double dt) const override {
    m_last_waypoint = m_traj.move_time_diff_from(
        wrapped_waypoint_pair(
            std::any_cast<wrapped_const_waypoint_descriptor>(a.first),
            a.second),
        dt);
    return waypoint_pair(m_last_waypoint.first, m_last_waypoint.second);
  }

  point_type get_point_at_time(double t) const override {
    m_last_waypoint = m_traj.get_waypoint_at_time(t);
    return m_last_waypoint.second;
  }

  waypoint_pair get_waypoint_at_time(double t) const override {
    m_last_waypoint = m_traj.get_waypoint_at_time(t);
    return waypoint_pair(m_last_waypoint.first, m_last_waypoint.second);
  }

  double get_start_time() const override { return m_traj.get_start_time(); }

  double get_end_time() const override { return m_traj.get_end_time(); }

  point_time_iterator begin_time_travel() const override {
    return point_time_iterator(
        new point_time_iterator_impl(m_traj.begin_time_travel()));
  }

  point_time_iterator end_time_travel() const override {
    return point_time_iterator(
        new point_time_iterator_impl(m_traj.end_time_travel()));
  }

  point_fraction_iterator begin_fraction_travel() const override {
    return point_fraction_iterator(
        new point_fraction_iterator_impl(m_traj.begin_fraction_travel()));
  }

  point_fraction_iterator end_fraction_travel() const override {
    return point_fraction_iterator(
        new point_fraction_iterator_impl(m_traj.end_fraction_travel()));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_traj);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_traj);
    m_last_waypoint = wrapped_waypoint_pair();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC244000A, 1, "trajectory_wrapper",
                              base_type)
};

}  // namespace ReaK::pp

#endif
