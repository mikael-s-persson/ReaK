/**
 * \file interpolated_topologies.hpp
 *
 *
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
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

#ifndef REAK_INTERPOLATED_TOPOLOGIES_HPP
#define REAK_INTERPOLATED_TOPOLOGIES_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/topologies/interpolation/generic_interpolator_factory.hpp>
#include <ReaK/topologies/spaces/default_random_sampler.hpp>
#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/proper_metric_concept.hpp>
#include <ReaK/topologies/spaces/rate_limited_space_metamaps.hpp>
#include <ReaK/topologies/spaces/reversible_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include <boost/concept_check.hpp>

#include <functional>
#include <type_traits>
#include <utility>

namespace ReaK::pp {

namespace detail {

template <typename BaseTopology, bool IsTemporal>
struct interp_topo_base_get_temporal_topos {
  using space_topology =
      typename temporal_space_traits<BaseTopology>::space_topology;
  using time_topology =
      typename temporal_space_traits<BaseTopology>::time_topology;
};

template <typename BaseTopology>
struct interp_topo_base_get_temporal_topos<BaseTopology, false> {
  using space_topology = int;
  using time_topology = int;
};

template <typename BaseTopology>
struct interp_topo_base_temporal_topos
    : interp_topo_base_get_temporal_topos<BaseTopology,
                                          is_temporal_space_v<BaseTopology>> {};
};  // namespace detail

/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and
 * sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class interpolated_topology_base : public BaseTopology {
 public:
  BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));

  using self = interpolated_topology_base<BaseTopology>;

  using point_type = topology_point_type_t<BaseTopology>;
  using point_difference_type = topology_point_difference_type_t<BaseTopology>;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  using TemporalToposSelect =
      detail::interp_topo_base_temporal_topos<BaseTopology>;
  using space_topology = typename TemporalToposSelect::space_topology;
  using time_topology = typename TemporalToposSelect::time_topology;

  using super_space_type = BaseTopology;

  static constexpr std::size_t dimensions =
      topology_traits<BaseTopology>::dimensions;

  using validity_predicate_type = std::function<bool(const point_type&)>;

 protected:
  virtual point_type interp_topo_move_position_toward(
      const point_type& a, double fraction, const point_type& b) const {
    return static_cast<const BaseTopology&>(*this).move_position_toward(
        a, fraction, b);
  }

  virtual point_type interp_topo_move_position_back_to(
      const point_type& a, double fraction, const point_type& b) const {
    return static_cast<const BaseTopology&>(*this).move_position_back_to(
        a, fraction, b);
  }

  virtual double interp_topo_get_distance(const point_type& a,
                                          const point_type& b) const {
    const auto& b_space = static_cast<const BaseTopology&>(*this);
    return get(distance_metric, b_space)(a, b, b_space);
  }

  virtual double interp_topo_get_proper_distance(const point_type& a,
                                                 const point_type& b) const {
    const auto& b_space = static_cast<const BaseTopology&>(*this);
    return get(proper_metric, b_space)(a, b, b_space);
  }

  virtual point_type interp_topo_move_position_toward_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval, validity_predicate_type predicate) const {

    double dist_tot = this->interp_topo_get_distance(a, b);
    if (dist_tot == std::numeric_limits<double>::infinity()) {
      return a;
    }
    if (dist_tot < min_dist_interval) {
      return this->interp_topo_move_position_toward(a, fraction, b);
    }

    double dist_inter = dist_tot * fraction;
    double dist_cur = min_dist_interval;
    point_type result = a;
    point_type last_result = a;
    while (dist_cur < dist_inter) {
      result =
          this->interp_topo_move_position_toward(a, dist_cur / dist_tot, b);
      if (!predicate(result)) {
        return last_result;
      }
      dist_cur += min_dist_interval;
      last_result = result;
    }
    if (fraction == 1.0) {
      return b;
    }
    if (fraction == 0.0) {
      return a;
    }

    return this->interp_topo_move_position_toward(a, fraction, b);
  }

  virtual point_type interp_topo_move_position_back_to_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval, validity_predicate_type predicate) const {

    double dist_tot = this->interp_topo_get_distance(a, b);
    if (dist_tot == std::numeric_limits<double>::infinity()) {
      return b;
    }
    if (dist_tot < min_dist_interval) {
      return this->interp_topo_move_position_back_to(a, fraction, b);
    }

    double dist_inter = dist_tot * fraction;
    double dist_cur = min_dist_interval;
    point_type result = b;
    point_type last_result = b;
    while (dist_cur < dist_inter) {
      result =
          this->interp_topo_move_position_back_to(a, dist_cur / dist_tot, b);
      if (!predicate(result)) {
        return last_result;
      }
      dist_cur += min_dist_interval;
      last_result = result;
    }
    if (fraction == 1.0) {
      return b;
    }
    if (fraction == 0.0) {
      return a;
    }

    return this->interp_topo_move_position_back_to(a, fraction, b);
  }

  virtual double interp_topo_get_distance_pred(
      const point_type& a, const point_type& b, double min_dist_interval,
      validity_predicate_type predicate) const {
    point_type b_tmp = this->interp_topo_move_position_toward_pred(
        a, 1.0, b, min_dist_interval, predicate);
    if (this->interp_topo_get_distance(b_tmp, b) <
        std::numeric_limits<double>::epsilon()) {
      return this->interp_topo_get_distance(a, b);  // if b is reachable from a.
    }
    return std::numeric_limits<
        double>::infinity();  // b is not reachable from a.
  }

  virtual double interp_topo_get_norm(const point_difference_type& dp) const {
    const auto& b_space = static_cast<const BaseTopology&>(*this);
    return get(distance_metric, b_space)(dp, b_space);
  }

  virtual double interp_topo_get_proper_norm(
      const point_difference_type& dp) const {
    const auto& b_space = static_cast<const BaseTopology&>(*this);
    return get(proper_metric, b_space)(dp, b_space);
  }

  virtual bool interp_topo_is_in_bounds(const point_type& a) const {
    return static_cast<const BaseTopology&>(*this).is_in_bounds(a);
  }

  virtual point_type interp_topo_random_point() const {
    const auto& b_space = static_cast<const BaseTopology&>(*this);
    return get(random_sampler, b_space)(b_space);
  }

  virtual point_type interp_topo_random_point_pred(
      validity_predicate_type predicate) const {
    point_type result;
    while (!predicate(result = this->interp_topo_random_point())) {
      // output only free C-space points.
    }
    return result;
  }

 public:
  // non-copyable
  interpolated_topology_base(const self& aTopo) = delete;
  self& operator=(const self& aTopo) = delete;

  explicit interpolated_topology_base(const BaseTopology& aTopo)
      : BaseTopology(aTopo) {}

  template <typename... Args>
  explicit interpolated_topology_base(Args&&... args)
      : BaseTopology(std::forward<Args>(args)...) {}

  /**
   * Returns a const-reference to the super-space of this topology.
   * \note This function returns a const-reference to itself since the super-space is also
   *       the base-class of this topology. The base class is not polymorphic, meaning that its
   *       distance metric and random-sampler are not overridden (non-virtual calls).
   */
  const super_space_type& get_super_space() const { return *this; }

  /*************************************************************************
   *                             MetricSpaceConcept
   * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& a, const point_type& b) const {
    return this->interp_topo_get_distance(a, b);
  }

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& a, const point_type& b,
                  double min_dist_interval,
                  validity_predicate_type predicate) const {
    return this->interp_topo_get_distance_pred(a, b, min_dist_interval,
                                               predicate);
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double norm(const point_difference_type& delta) const {
    return this->interp_topo_get_norm(delta);
  }

  /**
   * Returns the distance between two points.
   */
  double proper_distance(const point_type& a, const point_type& b) const {
    return this->interp_topo_get_proper_distance(a, b);
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double proper_norm(const point_difference_type& delta) const {
    return this->interp_topo_get_proper_norm(delta);
  }

  /*************************************************************************
  *                             BoundedSpaceConcept
  * **********************************************************************/

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const {
    return this->interp_topo_is_in_bounds(a);
  }

  /*************************************************************************
   *                         for PointDistributionConcept
   * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed within the reachable space.
   */
  point_type random_point() const { return this->interp_topo_random_point(); }

  /**
   * Generates a random point in the space, uniformly distributed within the reachable space.
   */
  point_type random_point(validity_predicate_type predicate) const {
    return this->interp_topo_random_point_pred(predicate);
  }

  /*************************************************************************
   *                             LieGroupConcept
   * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b) const {
    return this->interp_topo_move_position_toward(a, fraction, b);
  }

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b, double min_dist_interval,
                                  validity_predicate_type predicate) const {
    return this->interp_topo_move_position_toward_pred(
        a, fraction, b, min_dist_interval, predicate);
  }

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b) const {
    return this->interp_topo_move_position_back_to(a, fraction, b);
  }

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b,
                                   double min_dist_interval,
                                   validity_predicate_type predicate) const {
    return this->interp_topo_move_position_back_to_pred(
        a, fraction, b, min_dist_interval, predicate);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    BaseTopology::save(A, BaseTopology::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    BaseTopology::load(A, BaseTopology::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400039, 1, "interpolated_topology_base",
                              BaseTopology)
};

template <typename BaseTopology>
struct is_metric_space<interpolated_topology_base<BaseTopology>>
    : std::true_type {};

template <typename BaseTopology>
struct is_reversible_space<interpolated_topology_base<BaseTopology>>
    : std::true_type {};

template <typename BaseTopology>
struct is_point_distribution<interpolated_topology_base<BaseTopology>>
    : std::true_type {};

template <typename BaseTopology>
struct get_rate_illimited_space<interpolated_topology_base<BaseTopology>>
    : get_rate_illimited_space<BaseTopology> {};

template <typename BaseTopology>
struct is_temporal_space<interpolated_topology_base<BaseTopology>>
    : is_temporal_space<BaseTopology> {};

template <typename BaseTopology>
struct is_metric_symmetric<interpolated_topology_base<BaseTopology>>
    : is_metric_symmetric<BaseTopology> {};

template <typename BaseTopology>
struct get_proper_metric<interpolated_topology_base<BaseTopology>> {
  using type = default_proper_metric;
};

/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and
 * sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class wrapped_interp_topology : public named_object {
 public:
  using self = wrapped_interp_topology<BaseTopology>;

  using wrapped_base_type = interpolated_topology_base<BaseTopology>;

  using point_type = typename topology_traits<wrapped_base_type>::point_type;
  using point_difference_type =
      typename topology_traits<wrapped_base_type>::point_difference_type;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  using super_space_type = BaseTopology;

  using space_topology = typename wrapped_base_type::space_topology;
  using time_topology = typename wrapped_base_type::time_topology;

  using validity_predicate_type =
      typename wrapped_base_type::validity_predicate_type;

  static constexpr std::size_t dimensions =
      topology_traits<wrapped_base_type>::dimensions;

 private:
  std::shared_ptr<wrapped_base_type> m_space;

 public:
  explicit wrapped_interp_topology(std::shared_ptr<wrapped_base_type> aSpace =
                                       std::shared_ptr<wrapped_base_type>())
      : m_space(std::move(aSpace)) {
    setName("wrapped_interp_topology");
  }

  explicit operator wrapped_base_type&() { return *m_space; }
  explicit operator const wrapped_base_type&() const { return *m_space; }

  /**
   * Returns a const-reference to the super-space of this topology.
   * \note This function returns a const-reference to itself since the super-space is also
   *       the base-class of this topology. The base class is not polymorphic, meaning that its
   *       distance metric and random-sampler are not overridden (non-virtual calls).
   */
  const super_space_type& get_super_space() const {
    return m_space->get_super_space();
  }

  /** Returns the underlying space topology. */
  const space_topology& get_space_topology() const {
    return m_space->get_space_topology();
  }
  /** Returns the underlying time topology. */
  const time_topology& get_time_topology() const {
    return m_space->get_time_topology();
  }

  /*************************************************************************
   *                             TopologyConcept
   * **********************************************************************/

  /**
   * Returns the difference between two points (a - b).
   * \param a The first point.
   * \param b The second point.
   * \return The difference between the two points.
   */
  point_difference_type difference(const point_type& a,
                                   const point_type& b) const {
    return m_space->difference(a, b);
  }

  /**
   * Returns the addition of a point-difference to a point.
   * \param a The starting point.
   * \param delta The point-difference.
   * \return The addition of a point-difference to a point.
   */
  point_type adjust(const point_type& a,
                    const point_difference_type& delta) const {
    return m_space->adjust(a, delta);
  }

  /**
   * Returns the origin of the temporal-space.
   * \return The origin of the temporal-space.
   */
  point_type origin() const { return m_space->origin(); }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const {
    return m_space->is_in_bounds(a);
  }

  /*************************************************************************
   *                             MetricSpaceConcept
   * **********************************************************************/

  /**
   * Computes the distance between two points in the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \return The distance between a and b.
   */
  double distance(const point_type& a, const point_type& b) const {
    return m_space->distance(a, b);
  }

  double distance(const point_type& a, const point_type& b,
                  double min_dist_interval,
                  validity_predicate_type predicate) const {
    return m_space->distance(a, b, min_dist_interval, predicate);
  }

  /**
   * Computes the norm of a difference between two points.
   * \param a The difference between two points.
   * \return The norm of a difference between two points.
   */
  double norm(const point_difference_type& a) const { return m_space->norm(a); }

  /**
   * Computes the distance between two points in the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \return The distance between a and b.
   */
  double proper_distance(const point_type& a, const point_type& b) const {
    return m_space->proper_distance(a, b);
  }

  /**
   * Computes the norm of a difference between two points.
   * \param a The difference between two points.
   * \return The norm of a difference between two points.
   */
  double proper_norm(const point_difference_type& a) const {
    return m_space->proper_norm(a);
  }

  /*************************************************************************
   *                             LieGroupConcept
   * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   * \param a The first point.
   * \param fraction The fraction between the two points (0 to 1).
   * \param b The second point.
   * \return The point which is at a fraction between two points.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b) const {
    return m_space->move_position_toward(a, fraction, b);
  }

  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b, double min_dist_interval,
                                  validity_predicate_type predicate) const {
    return m_space->move_position_toward(a, fraction, b, min_dist_interval,
                                         predicate);
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b.
   * \param a The first point.
   * \param fraction The backward fraction between the two points (0 to 1).
   * \param b The second point.
   * \return The point which is at a backward fraction between two points.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b) const {
    return m_space->move_position_back_to(a, fraction, b);
  }

  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b,
                                   double min_dist_interval,
                                   validity_predicate_type predicate) const {
    return m_space->move_position_back_to(a, fraction, b, min_dist_interval,
                                          predicate);
  }

  /*************************************************************************
   *                             PointDistributionConcept
   * **********************************************************************/

  /**
   * Returns a random point within the temporal-space.
   * \return A random point within the temporal-space.
   */
  point_type random_point() const { return m_space->random_point(); }

  point_type random_point(validity_predicate_type predicate) const {
    return m_space->random_point(predicate);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_space);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_space);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240003D, 1, "wrapped_interp_topology",
                              named_object)
};

template <typename BaseTopology>
struct is_metric_space<wrapped_interp_topology<BaseTopology>> : std::true_type {
};

template <typename BaseTopology>
struct is_reversible_space<wrapped_interp_topology<BaseTopology>>
    : std::true_type {};

template <typename BaseTopology>
struct is_point_distribution<wrapped_interp_topology<BaseTopology>>
    : std::true_type {};

template <typename BaseTopology>
struct get_rate_illimited_space<wrapped_interp_topology<BaseTopology>>
    : get_rate_illimited_space<BaseTopology> {};

template <typename BaseTopology>
struct is_temporal_space<wrapped_interp_topology<BaseTopology>>
    : is_temporal_space<BaseTopology> {};

template <typename BaseTopology>
struct is_metric_symmetric<wrapped_interp_topology<BaseTopology>>
    : is_metric_symmetric<BaseTopology> {};

template <typename BaseTopology>
struct get_proper_metric<wrapped_interp_topology<BaseTopology>> {
  using type = default_proper_metric;
};

namespace detail {

template <typename InterpMethodTag, typename BaseTopology>
topology_point_type_t<BaseTopology> interp_topo_move_pt_impl(
    const BaseTopology& b_space, const topology_point_type_t<BaseTopology>& a,
    double fraction, const topology_point_type_t<BaseTopology>& b,
    double dist_tot) {
  if constexpr (is_temporal_space_v<BaseTopology>) {
    if (a.time > b.time) {  // Am I trying to go backwards in time (impossible)?
      return a;             // b is not reachable from a.
    }
  }

  if (dist_tot == std::numeric_limits<double>::infinity()) {
    return a;
  }
  if (fraction == 1.0) {
    return b;
  }
  if (fraction == 0.0) {
    return a;
  }

  auto result = a;
  if constexpr (is_temporal_space_v<BaseTopology>) {
    using SpaceTopoType =
        typename temporal_space_traits<BaseTopology>::space_topology;
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, SpaceTopoType,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type;

    InterpType interp;
    double dt_total =
        (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                      b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * fraction;
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                         b_space.get_time_topology(), dt, dt_total,
                         InterpFactoryType());
    result.time = a.time + dt;
  } else {
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, BaseTopology,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type;

    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(),
                      InterpFactoryType());
    double dist_inter = dist_tot * fraction;
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter,
                         dist_tot, InterpFactoryType());
  }
  return result;
}

template <typename InterpMethodTag, typename BaseTopology,
          typename PredFunction>
topology_point_type_t<BaseTopology> interp_topo_move_pt_impl(
    const BaseTopology& b_space, const topology_point_type_t<BaseTopology>& a,
    double fraction, const topology_point_type_t<BaseTopology>& b,
    double dist_tot, double min_interval, PredFunction predicate) {
  if constexpr (is_temporal_space_v<BaseTopology>) {
    if (a.time > b.time) {  // Am I trying to go backwards in time (impossible)?
      return a;             // b is not reachable from a.
    }
  }

  if (dist_tot == std::numeric_limits<double>::infinity()) {
    return a;
  }
  if (dist_tot < min_interval) {
    return interp_topo_move_pt_impl<InterpMethodTag>(b_space, a, fraction, b,
                                                     dist_tot);
  }
  if (fraction == 1.0) {
    return b;
  }
  if (fraction == 0.0) {
    return a;
  }

  auto result = a;
  if constexpr (is_temporal_space_v<BaseTopology>) {
    using SpaceTopoType =
        typename temporal_space_traits<BaseTopology>::space_topology;
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, SpaceTopoType,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type;

    InterpType interp;
    double dt_total =
        (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                      b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * fraction;
    double d = min_interval;
    auto last_result = a;
    while (d < dt) {
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), d, dt_total,
                           InterpFactoryType());
      result.time = a.time + d;
      if (!predicate(result)) {
        return last_result;
      }
      d += min_interval;
      last_result = result;
    }
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                         b_space.get_time_topology(), dt, dt_total,
                         InterpFactoryType());
    result.time = a.time + dt;
  } else {
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, BaseTopology,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type;

    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(),
                      InterpFactoryType());
    double dist_inter = dist_tot * fraction;
    double dist_cur = min_interval;
    auto last_result = a;
    while (dist_cur < dist_inter) {
      interp.compute_point(result, a, b, b_space, time_topology(), dist_cur,
                           dist_tot, InterpFactoryType());
      if (!predicate(result)) {
        return last_result;
      }
      dist_cur += min_interval;
      last_result = result;
    }
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter,
                         dist_tot, InterpFactoryType());
  }
  return result;
}

template <typename InterpMethodTag, typename BaseTopology>
topology_point_type_t<BaseTopology> interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, const topology_point_type_t<BaseTopology>& a,
    double fraction, const topology_point_type_t<BaseTopology>& b,
    double dist_tot) {
  if constexpr (is_temporal_space_v<BaseTopology>) {
    if (a.time > b.time) {  // Am I trying to go backwards in time (impossible)?
      return a;             // b is not reachable from a.
    }
  }

  if (dist_tot == std::numeric_limits<double>::infinity()) {
    return b;
  }
  if (fraction == 1.0) {
    return a;
  }
  if (fraction == 0.0) {
    return b;
  }

  auto result = a;

  if constexpr (is_temporal_space_v<BaseTopology>) {
    using SpaceTopoType =
        typename temporal_space_traits<BaseTopology>::space_topology;
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, SpaceTopoType,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type;

    InterpType interp;
    double dt_total =
        (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                      b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * (1.0 - fraction);
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                         b_space.get_time_topology(), dt, dt_total,
                         InterpFactoryType());
    result.time = a.time + dt;
  } else {
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, BaseTopology,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type;

    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(),
                      InterpFactoryType());
    double dist_inter = dist_tot * (1.0 - fraction);
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter,
                         dist_tot, InterpFactoryType());
  }
  return result;
}

template <typename InterpMethodTag, typename BaseTopology,
          typename PredFunction>
topology_point_type_t<BaseTopology> interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, const topology_point_type_t<BaseTopology>& a,
    double fraction, const topology_point_type_t<BaseTopology>& b,
    double dist_tot, double min_interval, PredFunction predicate) {
  if constexpr (is_temporal_space_v<BaseTopology>) {
    if (a.time > b.time) {  // Am I trying to go backwards in time (impossible)?
      return a;             // b is not reachable from a.
    }
  }

  if (dist_tot == std::numeric_limits<double>::infinity()) {
    return b;
  }
  if (dist_tot < min_interval) {
    return interp_topo_move_pt_back_impl<InterpMethodTag>(b_space, a, fraction,
                                                          b, dist_tot);
  }
  if (fraction == 1.0) {
    return a;
  }
  if (fraction == 0.0) {
    return b;
  }

  auto result = b;

  if constexpr (is_temporal_space_v<BaseTopology>) {
    using SpaceTopoType =
        typename temporal_space_traits<BaseTopology>::space_topology;
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, SpaceTopoType,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type;

    InterpType interp;
    double dt_total =
        (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                      b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * (1.0 - fraction);
    double d = dt_total - min_interval;
    auto last_result = b;
    while (d > dt) {
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), d, dt_total,
                           InterpFactoryType());
      result.time = a.time + d;
      if (!predicate(result)) {
        return last_result;
      }
      d -= min_interval;
      last_result = result;
    }
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                         b_space.get_time_topology(), dt, dt_total,
                         InterpFactoryType());
    result.time = a.time + dt;
  } else {
    using InterpType =
        get_tagged_spatial_interpolator_t<InterpMethodTag, BaseTopology,
                                          time_topology>;
    using InterpFactoryType = typename get_tagged_spatial_interpolator<
        InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type;

    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(),
                      InterpFactoryType());
    double dist_inter = dist_tot * (1.0 - fraction);
    double dist_cur = dist_tot - min_interval;
    auto last_result = b;
    while (dist_cur > dist_inter) {
      interp.compute_point(result, a, b, b_space, time_topology(), dist_cur,
                           dist_tot, InterpFactoryType());
      if (!predicate(result)) {
        return last_result;
      }
      dist_cur -= min_interval;
      last_result = result;
    }
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter,
                         dist_tot, InterpFactoryType());
  }
  return result;
}

}  // namespace detail

/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and
 * sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology, typename InterpMethodTag>
class interpolated_topology : public interpolated_topology_base<BaseTopology> {
 public:
  using base_type = interpolated_topology_base<BaseTopology>;
  using self = interpolated_topology<BaseTopology, InterpMethodTag>;

  using point_type = typename base_type::point_type;
  using point_difference_type = typename base_type::point_difference_type;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  using space_topology = typename base_type::space_topology;
  using time_topology = typename base_type::time_topology;

  static constexpr std::size_t dimensions = base_type::dimensions;

  using validity_predicate_type = typename base_type::validity_predicate_type;

 protected:
  point_type interp_topo_move_position_toward(
      const point_type& a, double fraction,
      const point_type& b) const override {
    return detail::interp_topo_move_pt_impl<InterpMethodTag, BaseTopology>(
        static_cast<const BaseTopology&>(*this), a, fraction, b,
        this->interp_topo_get_distance(a, b));
  }

  point_type interp_topo_move_position_toward_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval,
      validity_predicate_type predicate) const override {
    return detail::interp_topo_move_pt_impl<InterpMethodTag, BaseTopology,
                                            validity_predicate_type>(
        static_cast<const BaseTopology&>(*this), a, fraction, b,
        this->interp_topo_get_distance(a, b), min_dist_interval, predicate);
  }

  point_type interp_topo_move_position_back_to(
      const point_type& a, double fraction,
      const point_type& b) const override {
    return detail::interp_topo_move_pt_back_impl<InterpMethodTag, BaseTopology>(
        static_cast<const BaseTopology&>(*this), a, fraction, b,
        this->interp_topo_get_distance(a, b));
  }

  point_type interp_topo_move_position_back_to_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval,
      validity_predicate_type predicate) const override {
    return detail::interp_topo_move_pt_back_impl<InterpMethodTag, BaseTopology,
                                                 validity_predicate_type>(
        static_cast<const BaseTopology&>(*this), a, fraction, b,
        this->interp_topo_get_distance(a, b), min_dist_interval, predicate);
  }

 public:
  explicit interpolated_topology(const BaseTopology& aTopo)
      : base_type(aTopo) {}

  template <typename... Args>
  explicit interpolated_topology(Args&&... args)
      : base_type(std::forward<Args>(args)...) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240003A, 1, "interpolated_topology",
                              base_type)
};

template <typename BaseTopology, typename InterpMethodTag>
struct is_metric_space<interpolated_topology<BaseTopology, InterpMethodTag>>
    : std::true_type {};

template <typename BaseTopology, typename InterpMethodTag>
struct is_reversible_space<interpolated_topology<BaseTopology, InterpMethodTag>>
    : std::true_type {};

template <typename BaseTopology, typename InterpMethodTag>
struct is_point_distribution<
    interpolated_topology<BaseTopology, InterpMethodTag>> : std::true_type {};

template <typename BaseTopology, typename InterpMethodTag>
struct get_rate_illimited_space<
    interpolated_topology<BaseTopology, InterpMethodTag>>
    : get_rate_illimited_space<BaseTopology> {};

template <typename BaseTopology, typename InterpMethodTag>
struct is_temporal_space<interpolated_topology<BaseTopology, InterpMethodTag>>
    : is_temporal_space<BaseTopology> {};

template <typename BaseTopology, typename InterpMethodTag>
struct is_metric_symmetric<interpolated_topology<BaseTopology, InterpMethodTag>>
    : std::integral_constant<bool, is_metric_symmetric_v<InterpMethodTag> &&
                                       is_metric_symmetric_v<BaseTopology>> {};

template <typename BaseTopology, typename InterpMethodTag>
struct get_proper_metric<interpolated_topology<BaseTopology, InterpMethodTag>> {
  using type = default_proper_metric;
};

}  // namespace ReaK::pp

namespace ReaK {

/* Specialization, see general template docs. */
template <typename BaseTopology>
struct arithmetic_tuple_size<pp::interpolated_topology_base<BaseTopology>>
    : arithmetic_tuple_size<BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology>
struct arithmetic_tuple_element<Idx,
                                pp::interpolated_topology_base<BaseTopology>>
    : arithmetic_tuple_element<Idx, BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology>
struct arithmetic_tuple_element<
    Idx, const pp::interpolated_topology_base<BaseTopology>>
    : arithmetic_tuple_element<Idx, const BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology>
struct arithmetic_tuple_element<
    Idx, volatile pp::interpolated_topology_base<BaseTopology>>
    : arithmetic_tuple_element<Idx, volatile BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology>
struct arithmetic_tuple_element<
    Idx, const volatile pp::interpolated_topology_base<BaseTopology>>
    : arithmetic_tuple_element<Idx, const volatile BaseTopology> {};

/* Specialization, see general template docs. */
template <typename BaseTopology, typename InterpMethodTag>
struct arithmetic_tuple_size<
    pp::interpolated_topology<BaseTopology, InterpMethodTag>>
    : arithmetic_tuple_size<BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology, typename InterpMethodTag>
struct arithmetic_tuple_element<
    Idx, pp::interpolated_topology<BaseTopology, InterpMethodTag>>
    : arithmetic_tuple_element<Idx, BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology, typename InterpMethodTag>
struct arithmetic_tuple_element<
    Idx, const pp::interpolated_topology<BaseTopology, InterpMethodTag>>
    : arithmetic_tuple_element<Idx, const BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology, typename InterpMethodTag>
struct arithmetic_tuple_element<
    Idx, volatile pp::interpolated_topology<BaseTopology, InterpMethodTag>>
    : arithmetic_tuple_element<Idx, volatile BaseTopology> {};

/* Specialization, see general template docs. */
template <int Idx, typename BaseTopology, typename InterpMethodTag>
struct arithmetic_tuple_element<Idx, const volatile pp::interpolated_topology<
                                         BaseTopology, InterpMethodTag>>
    : arithmetic_tuple_element<Idx, const volatile BaseTopology> {};

}  // namespace ReaK

#endif
