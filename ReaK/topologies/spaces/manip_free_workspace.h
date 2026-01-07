/**
 * \file manip_free_workspace.h
 *
 * This library defines a class for path-planning for a manipulator moving inside an environment
 * with obstacles and physical limits (joint limits). This class is essentially just an assembly
 * of many of the building blocks in the ReaK path-planning library. This class can also draw the
 * elements of the motion graph (edges) as end-effector trajectories in a Coin3D scene-graph.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_TOPOLOGIES_SPACES_MANIP_FREE_WORKSPACE_H_
#define REAK_TOPOLOGIES_SPACES_MANIP_FREE_WORKSPACE_H_

#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"

#include "ReaK/topologies/spaces/interpolated_topologies.h"
#include "ReaK/topologies/spaces/proxy_model_updater.h"

#include "ReaK/geometry/proximity/proxy_query_model.h"  // for proxy-query class

#include "ReaK/topologies/spaces/default_random_sampler.h"

#include <functional>
#include <utility>

namespace ReaK::pp {

namespace detail {

template <Topology BaseJointSpace>
class manip_dk_proxy_env_impl {
 public:
  using point_type = topology_point_type_t<BaseJointSpace>;

  std::shared_ptr<proxy_model_applicator<BaseJointSpace>> m_applicator;

  std::vector<std::shared_ptr<geom::proxy_query_pair_2D>> m_proxy_env_2D;
  std::vector<std::shared_ptr<geom::proxy_query_pair_3D>> m_proxy_env_3D;

  explicit manip_dk_proxy_env_impl(
      std::shared_ptr<proxy_model_applicator<BaseJointSpace>> aApplicator = {})
      : m_applicator(std::move(aApplicator)) {}

  bool is_free(const point_type& pt, const BaseJointSpace& space) const {
    if (m_applicator) {
      m_applicator->apply_to_model(pt, space);
    }

    // NOTE: it is assumed that the proxy environments are properly linked with the manipulator kinematic model.

    const auto not_penetrating = [](const auto& sp) {
      return sp->findMinimumDistance().mDistance >= 0.0;
    };
    return std::all_of(m_proxy_env_2D.begin(), m_proxy_env_2D.end(),
                       not_penetrating) &&
           std::all_of(m_proxy_env_3D.begin(), m_proxy_env_3D.end(),
                       not_penetrating);
  }
};

}  // namespace detail

/**
 * This class is used to represent the free-space of a manipulator in a quasi-static environment.
 * Here, the term quasi-static refers to the fact that proximity queries are performed against an
 * environment that is assumed to be unchanging (at least, not significantly while this topology is used).
 */
template <Topology BaseJointSpace>
class manip_quasi_static_env : public named_object {
 public:
  using self = manip_quasi_static_env<BaseJointSpace>;
  using super_space_base_type = interpolated_topology_base<BaseJointSpace>;
  using super_space_type = wrapped_interp_topology<BaseJointSpace>;
  using point_type = topology_point_type_t<super_space_type>;
  using point_difference_type =
      topology_point_difference_type_t<super_space_type>;

  using applicator_type = proxy_model_applicator<BaseJointSpace>;
  using dk_proxy_env_type = detail::manip_dk_proxy_env_impl<BaseJointSpace>;

  static constexpr std::size_t dimensions =
      topology_traits<super_space_type>::dimensions;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

 private:
  double min_interval;
  super_space_type m_space;
  dk_proxy_env_type m_prox_env;

 public:
  /**
   * Returns a reference to the super-space in which this test-world is embedded.
   * \return A reference to the super-space in which this test-world is embedded.
   */
  super_space_type& get_super_space() { return m_space; }

  /**
   * Returns a const-reference to the super-space in which this test-world is embedded.
   * \return A const-reference to the super-space in which this test-world is embedded.
   */
  const super_space_type& get_super_space() const { return m_space; }

  /**
   * Checks if the given point is within the free-space.
   * \param p The point to be checked for being collision-free.
   * \return True if p is collision-free.
   */
  bool is_free(const point_type& p) const {
    return m_space.is_in_bounds(p) &&
           m_prox_env.is_free(p, m_space.get_super_space());
  }

  struct is_free_predicate {
    const self* parent;
    explicit is_free_predicate(const self* aParent) : parent(aParent) {}
    bool operator()(const point_type& p) const { return parent->is_free(p); }
  };

  // Topology concepts:

  /**
   * Produces a random, collision-free point.
   * \return A random, collision-free point.
   */
  point_type random_point() const {
    return m_space.random_point(is_free_predicate(this));
  }

  /**
   * Computes the distance between two points. If there is no collision-free line between
   * the two points, the distance is infinite.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return The collision-free distance between the two given points.
   */
  double distance(const point_type& p1, const point_type& p2) const {
    return m_space.distance(p1, p2, min_interval, is_free_predicate(this));
  }

  /**
   * Computes the norm of the difference between two points.
   * \param dp The point difference.
   * \return The norm of the difference between the two points.
   */
  double norm(const point_difference_type& dp) const {
    return m_space.norm(dp);
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
    return m_space.adjust(p, dp);
  }

  /**
   * Returns a point which is at a fraction between two points a to b, or as
   * far as it can get before a collision.
   */
  point_type move_position_toward(const point_type& p1, double fraction,
                                  const point_type& p2) const {
    return m_space.move_position_toward(p1, fraction, p2, min_interval,
                                        is_free_predicate(this));
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b, or as
   * far as it can get before a collision.
   */
  point_type move_position_back_to(const point_type& p1, double fraction,
                                   const point_type& p2) const {
    return m_space.move_position_back_to(p1, fraction, p2, min_interval,
                                         is_free_predicate(this));
  }

  /**
   * Parametrized constructor (this class is a RAII class).
   * \param aMinInterval The minimum length of the travel between collision detection calls.
   */
  explicit manip_quasi_static_env(
      const BaseJointSpace& aSpace = BaseJointSpace{},
      const std::shared_ptr<applicator_type>& aApplicator = {},
      double aMinInterval = 0.1)
      : min_interval(aMinInterval),
        m_space(std::make_shared<super_space_base_type>(aSpace)),
        m_prox_env(aApplicator) {}

  /**
   * Parametrized constructor (this class is a RAII class).
   * \param aMinInterval The minimum length of the travel between collision detection calls.
   */
  template <typename InterpMethodTag>
  explicit manip_quasi_static_env(
      InterpMethodTag aInterpTag,
      const BaseJointSpace& aSpace = BaseJointSpace{},
      const std::shared_ptr<applicator_type>& aApplicator = {},
      double aMinInterval = 0.1)
      : min_interval(aMinInterval),
        m_space(std::make_shared<
                interpolated_topology<BaseJointSpace, InterpMethodTag>>(
            aSpace)),
        m_prox_env(aApplicator) {}

  ~manip_quasi_static_env() override = default;

  /**
   * Add a 2D proxy query pair to the collision environment.
   * \param aProxy The new 2D proxy query pair to add to the collision environment.
   * \return A reference back to 'this'.
   */
  self& operator<<(const std::shared_ptr<geom::proxy_query_pair_2D>& aProxy) {
    m_prox_env.m_proxy_env_2D.push_back(aProxy);
    return *this;
  }

  /**
   * Add a 3D proxy query pair to the collision environment.
   * \param aProxy The new 3D proxy query pair to add to the collision environment.
   * \return A reference back to 'this'.
   */
  self& operator<<(const std::shared_ptr<geom::proxy_query_pair_3D>& aProxy) {
    m_prox_env.m_proxy_env_3D.push_back(aProxy);
    return *this;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(min_interval) &
        RK_SERIAL_SAVE_WITH_NAME(m_space) &
        RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_applicator) &
        RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_2D) &
        RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_3D);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(min_interval) &
        RK_SERIAL_LOAD_WITH_NAME(m_space) &
        RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_applicator) &
        RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_2D) &
        RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400027, 1, "manip_quasi_static_env",
                              named_object)
};

template <typename BaseJointSpace>
struct is_metric_symmetric<manip_quasi_static_env<BaseJointSpace>>
    : is_metric_symmetric<
          typename manip_quasi_static_env<BaseJointSpace>::super_space_type> {};

}  // namespace ReaK::pp

#include "ReaK/topologies/spaces/manip_free_workspace_ext.h"

#endif  // REAK_TOPOLOGIES_SPACES_MANIP_FREE_WORKSPACE_H_
