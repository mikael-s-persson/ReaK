/**
 * \file no_obstacle_space.h
 *
 * This library defines a class template for path-planning problems on a topology without
 * any obstacles. This class can also be used to attach a new distance-metric or random-sampler
 * to a given topology type.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_TOPOLOGIES_SPACES_NO_OBSTACLE_SPACE_H_
#define REAK_TOPOLOGIES_SPACES_NO_OBSTACLE_SPACE_H_

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/proper_metric_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"

namespace ReaK::pp {

/**
 * This class is used for path-planning problems on a topology without
 * any obstacles. This class can also be used to attach a new distance-metric or random-sampler
 * to a given topology type.
 */
template <typename Topology,
          typename DistanceMetric =
              typename metric_space_traits<Topology>::distance_metric_type,
          typename RandomSampler =
              typename point_distribution_traits<Topology>::random_sampler_type>
class no_obstacle_space : public named_object {
 public:
  using self = no_obstacle_space<Topology, DistanceMetric, RandomSampler>;
  using super_space_type = Topology;
  using point_type = topology_point_type_t<super_space_type>;
  using point_difference_type =
      topology_point_difference_type_t<super_space_type>;

  static constexpr std::size_t dimensions =
      topology_traits<super_space_type>::dimensions;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

 private:
  double max_edge_length;

  super_space_type m_space;
  DistanceMetric m_distance;
  RandomSampler m_rand_sampler;

  point_type m_start_point;
  point_type m_goal_point;

 public:
  double get_max_edge_length() const { return max_edge_length; }

  /**
   * Returns a reference to the super-space in which this topology is embedded.
   * \return A reference to the super-space in which this topology is embedded.
   */
  super_space_type& get_super_space() { return m_space; }

  /**
   * Returns a const-reference to the super-space in which this topology is embedded.
   * \return A const-reference to the super-space in which this topology is embedded.
   */
  const super_space_type& get_super_space() const { return m_space; }

  /**
   * Checks if the given point is within the free-space.
   * \param p The point to be checked for being collision-free.
   * \return True if p is collision-free.
   */
  bool is_free(const point_type& pt) const { return m_space.is_in_bounds(pt); }

  // Topology concepts:

  /**
   * Produces a random, collision-free point.
   * \return A random, collision-free point.
   */
  point_type random_point() const { return m_rand_sampler(m_space); }

  /**
   * Computes the distance between two points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return The distance between the two given points.
   */
  double distance(const point_type& p1, const point_type& p2) const {
    if (m_distance(p2, move_position_toward(p1, 1.0, p2), m_space) <
        std::numeric_limits<double>::epsilon()) {
      return m_distance(
          p1, p2,
          m_space);  // if p2 is reachable from p1, use Euclidean distance.
    }
    return std::numeric_limits<
        double>::infinity();  // p2 is not reachable from p1.
  }

  /**
   * Computes the norm of the difference between two points.
   * \param dp The point difference.
   * \return The norm of the difference between the two points.
   */
  double norm(const point_difference_type& dp) const {
    return m_distance(dp, m_space);
  }

  /**
   * Returns the difference between two points (a - b).
   */
  point_difference_type difference(const point_type& p1,
                                   const point_type& p2) const {
    return m_space.difference(p1, p2);
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type origin() const { return m_space.origin(); }

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type adjust(const point_type& p,
                    const point_difference_type& dp) const {
    return move_position_toward(p, 1.0, m_space.adjust(p, dp));
  }

  /**
   * Returns a point which is at a fraction between two points a to b, or as
   * far as it can get before a collision.
   */
  point_type move_position_toward(const point_type& p1, double fraction,
                                  const point_type& p2) const {
    double total_dist = m_distance(p1, p2, m_space);
    if (total_dist * fraction > max_edge_length) {
      return m_space.move_position_toward(p1, max_edge_length / total_dist, p2);
    }
    return m_space.move_position_toward(p1, fraction, p2);
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b, or as
   * far as it can get before a collision.
   */
  point_type move_position_back_to(const point_type& p1, double fraction,
                                   const point_type& p2) const {
    double total_dist = m_distance(p1, p2, m_space);
    if (total_dist * fraction > max_edge_length) {
      return m_space.move_position_back_to(p1, max_edge_length / total_dist,
                                           p2);
    }
    return m_space.move_position_back_to(p1, fraction, p2);
  }

  /**
   * Returns a random point fairly near to the given point.
   */
  std::pair<point_type, bool> random_walk(const point_type& p_u) const {
    return {move_position_toward(
                p_u, 1.0,
                m_space.adjust(p_u, m_space.difference(m_rand_sampler(m_space),
                                                       m_space.origin()))),
            true};
  }

  /**
   * Returns a random point fairly near to the given point.
   */
  std::pair<point_type, bool> random_back_walk(const point_type& p_u) const {
    return {move_position_back_to(
                m_space.adjust(p_u, -m_space.difference(m_rand_sampler(m_space),
                                                        m_space.origin())),
                1.0, p_u),
            true};
  }

  /**
   * Returns the bird-flight distance to the goal from the given point.
   * \param p_u The point from which the bird-flight distance is sought.
   * \return The bird-flight distance to the goal from the given point.
   */
  double bird_fly_to_goal(const point_type& p_u) const {
    return m_distance(p_u, m_goal_point, m_space);
  }

  /**
   * Returns the bird-flight distance to the given point from the start.
   * \param p_u The point from which the bird-flight distance is sought.
   * \return The bird-flight distance to the given point from the start.
   */
  double bird_fly_to_start(const point_type& p_u) const {
    return m_distance(m_start_point, p_u, m_space);
  }

  /**
   * Returns the start point.
   * \return The start point.
   */
  const point_type& get_start_pos() const { return m_start_point; }

  /**
   * Returns the goal point.
   * \return The goal point.
   */
  const point_type& get_goal_pos() const { return m_goal_point; }

  /**
   * Set the start point.
   * \param aStart The new start point.
   */
  void set_start_pos(const point_type& aStart) { m_start_point = aStart; }

  /**
   * Set the goal point.
   * \param aStart The new goal point.
   */
  void set_goal_pos(const point_type& aGoal) { m_goal_point = aGoal; }

  /**
   * Parametrized constructor (this class is a RAII class).
   * \param aName The name of the topology.
   * \param aSpace The super-space topology.
   * \param aMaxEdgeLength The maximum length of an added edge, in pixel-units.
   * \param aDistMetric The distance metric functor to use.
   * \param aRandSampler The random sampler functor to use.
   */
  explicit no_obstacle_space(
      const std::string& aName,
      const super_space_type& aSpace = super_space_type{},
      double aMaxEdgeLength = std::numeric_limits<double>::infinity(),
      DistanceMetric aDistMetric = DistanceMetric{},
      RandomSampler aRandSampler = RandomSampler{})
      : named_object(),
        max_edge_length(aMaxEdgeLength),
        m_space(aSpace),
        m_distance(aDistMetric),
        m_rand_sampler(aRandSampler),
        m_start_point(aRandSampler(aSpace)),
        m_goal_point(aRandSampler(aSpace)) {
    set_name(aName);
  }

  no_obstacle_space() : no_obstacle_space("") {}

  ~no_obstacle_space() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(max_edge_length) &
        RK_SERIAL_SAVE_WITH_NAME(m_space) &
        RK_SERIAL_SAVE_WITH_NAME(m_distance) &
        RK_SERIAL_SAVE_WITH_NAME(m_rand_sampler) &
        RK_SERIAL_SAVE_WITH_NAME(m_start_point) &
        RK_SERIAL_SAVE_WITH_NAME(m_goal_point);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(max_edge_length) &
        RK_SERIAL_LOAD_WITH_NAME(m_space) &
        RK_SERIAL_LOAD_WITH_NAME(m_distance) &
        RK_SERIAL_LOAD_WITH_NAME(m_rand_sampler) &
        RK_SERIAL_LOAD_WITH_NAME(m_start_point) &
        RK_SERIAL_LOAD_WITH_NAME(m_goal_point);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400021, 1, "no_obstacle_space",
                              named_object)
};

template <typename Topology, typename DistanceMetric, typename RandomSampler>
struct is_metric_symmetric<
    no_obstacle_space<Topology, DistanceMetric, RandomSampler>>
    : std::integral_constant<bool, is_metric_symmetric_v<DistanceMetric> &&
                                       is_metric_symmetric_v<Topology>> {};

template <typename Topology, typename DistanceMetric, typename RandomSampler>
struct get_proper_metric<
    no_obstacle_space<Topology, DistanceMetric, RandomSampler>>
    : get_proper_metric_from_metric<DistanceMetric> {};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_NO_OBSTACLE_SPACE_H_
