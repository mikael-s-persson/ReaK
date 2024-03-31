/**
 * \file svp_reach_topologies.h
 *
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a rate-limited sustained velocity pulse (SVP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SVP_REACH_TOPOLOGIES_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SVP_REACH_TOPOLOGIES_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/optimization/optim_exceptions.h"
#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"
#include "ReaK/topologies/interpolation/sustained_velocity_pulse.h"
#include "ReaK/topologies/interpolation/svp_metrics.h"
#include "ReaK/topologies/interpolation/svp_samplers.h"
#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/generic_interpolator_factory.h"
#include "ReaK/topologies/spaces/generic_sampler_factory.h"
#include "ReaK/topologies/spaces/interpolated_topologies.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"

namespace ReaK::pp {

namespace detail {

template <typename BaseTopology, bool IsTemporal>
struct svp_reach_topo_impl {

  using point_type = topology_point_type_t<BaseTopology>;
  using point_difference_type = topology_point_difference_type_t<BaseTopology>;

  using rt_metric_type = svp_reach_time_metric<time_topology>;
  using sampler_type =
      generic_sampler<svp_rate_limited_sampler<time_topology>, BaseTopology>;

  static rt_metric_type make_rt_metric(const BaseTopology&) {
    return rt_metric_type();
  }
  static sampler_type make_sampler(const BaseTopology&) {
    return sampler_type();
  }

  using validity_predicate_type = std::function<bool(const point_type&)>;

  static point_type move_pt_toward(const BaseTopology& b_space,
                                   const rt_metric_type& rt_dist,
                                   const point_type& a, double fraction,
                                   const point_type& b) {
    try {
      generic_interpolator_impl<svp_interpolator, BaseTopology, time_topology>
          interp;
      interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
      double dt_min = interp.get_minimum_travel_time();
      if (dt_min == std::numeric_limits<double>::infinity()) {
        return a;
      }
      double dt = dt_min * fraction;
      point_type result = a;
      interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min,
                           rt_dist);
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return a;
    }
  }

  static point_type move_pt_toward(const BaseTopology& b_space,
                                   const rt_metric_type& rt_dist,
                                   const point_type& a, double fraction,
                                   const point_type& b,
                                   double min_dist_interval,
                                   validity_predicate_type predicate) {
    try {
      generic_interpolator_impl<svp_interpolator, BaseTopology, time_topology>
          interp;
      interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
      double dt_min = interp.get_minimum_travel_time();
      if (dt_min == std::numeric_limits<double>::infinity()) {
        return a;
      }
      double dt = dt_min * fraction;
      double d = min_dist_interval;
      point_type result = a;
      point_type last_result = a;
      while (d < dt) {
        interp.compute_point(result, a, b, b_space, time_topology(), d, dt_min,
                             rt_dist);
        if (!predicate(result)) {
          return last_result;
        }
        d += min_dist_interval;
        last_result = result;
      }
      if (fraction == 1.0) {
        return b;
      }
      if (fraction == 0.0) {
        return a;
      }
      interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min,
                           rt_dist);
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return a;
    }
  }

  static point_type move_pt_back_to(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, double fraction,
                                    const point_type& b) {
    try {
      generic_interpolator_impl<svp_interpolator, BaseTopology, time_topology>
          interp;
      interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
      double dt_min = interp.get_minimum_travel_time();
      if (dt_min == std::numeric_limits<double>::infinity()) {
        return b;
      }
      double dt = dt_min * (1.0 - fraction);
      point_type result = b;
      interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min,
                           rt_dist);
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return b;
    }
  }

  static point_type move_pt_back_to(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, double fraction,
                                    const point_type& b,
                                    double min_dist_interval,
                                    validity_predicate_type predicate) {
    try {
      generic_interpolator_impl<svp_interpolator, BaseTopology, time_topology>
          interp;
      interp.initialize(a, b, 0.0, b_space, time_topology(), rt_dist);
      double dt_min = interp.get_minimum_travel_time();
      if (dt_min == std::numeric_limits<double>::infinity()) {
        return b;
      }
      double dt = dt_min * (1.0 - fraction);
      double d = dt_min - min_dist_interval;
      point_type result = b;
      point_type last_result = b;
      while (d > dt) {
        interp.compute_point(result, a, b, b_space, time_topology(), d, dt_min,
                             rt_dist);
        if (!predicate(result)) {
          return last_result;
        }
        d -= min_dist_interval;
        last_result = result;
      }
      if (fraction == 1.0) {
        return a;
      }
      if (fraction == 0.0) {
        return b;
      }
      interp.compute_point(result, a, b, b_space, time_topology(), dt, dt_min,
                           rt_dist);
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return b;
    }
  }

  static double get_distance(const BaseTopology& b_space,
                             const rt_metric_type& rt_dist, const point_type& a,
                             const point_type& b) {
    return rt_dist(a, b, b_space);
  }

  static double get_norm(const BaseTopology& b_space,
                         const rt_metric_type& rt_dist,
                         const point_difference_type& dp) {
    return rt_dist(dp, b_space);
  }

  static double get_proper_distance(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, const point_type& b) {
    svp_reach_time_metric<time_topology, true> p_rt_dist(rt_dist);
    return p_rt_dist(a, b, b_space);
  }

  static double get_proper_norm(const BaseTopology& b_space,
                                const rt_metric_type& rt_dist,
                                const point_difference_type& dp) {
    svp_reach_time_metric<time_topology, true> p_rt_dist(rt_dist);
    return p_rt_dist(dp, b_space);
  }

  static bool is_in_bounds(const BaseTopology& b_space, const point_type& a) {
    return svp_is_in_bounds(a, b_space, time_topology());
  }

  static point_type random_point(const BaseTopology& b_space,
                                 const sampler_type& rl_sampler) {
    return rl_sampler(b_space);
  }
};

// Implementation for the temporal spaces:
template <typename BaseTopology>
struct svp_reach_topo_impl<BaseTopology, true> {

  using point_type = topology_point_type_t<BaseTopology>;
  using point_difference_type = topology_point_difference_type_t<BaseTopology>;

  using base_time_topo =
      typename temporal_space_traits<BaseTopology>::time_topology;
  using base_space_topo =
      typename temporal_space_traits<BaseTopology>::space_topology;

  using rt_metric_type = svp_reach_time_metric<base_time_topo>;
  using sampler_type = generic_sampler<svp_rate_limited_sampler<base_time_topo>,
                                       base_space_topo>;

  static rt_metric_type make_rt_metric(const BaseTopology& b_space) {
    return rt_metric_type(std::shared_ptr<const base_time_topo>(
        &(b_space.get_time_topology()), null_deleter()));
  }

  static sampler_type make_sampler(const BaseTopology& b_space) {
    return sampler_type(std::shared_ptr<const base_time_topo>(
        &(b_space.get_time_topology()), null_deleter()));
  }

  using validity_predicate_type = std::function<bool(const point_type&)>;

  static point_type move_pt_toward(const BaseTopology& b_space,
                                   const rt_metric_type& rt_dist,
                                   const point_type& a, double fraction,
                                   const point_type& b) {

    if (a.time > b.time) {
      return a;
    }

    try {
      generic_interpolator_impl<svp_interpolator, base_space_topo,
                                base_time_topo>
          interp;
      double dt_total =
          (b.time - a.time);  // the free time that I have along the path.
      interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                        b_space.get_time_topology(), rt_dist);
      if (interp.get_minimum_travel_time() ==
          std::numeric_limits<double>::infinity()) {
        return a;
      }
      double dt = dt_total * fraction;
      point_type result = a;
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), dt, dt_total, rt_dist);
      result.time += dt;
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return a;
    }
  }

  static point_type move_pt_toward(const BaseTopology& b_space,
                                   const rt_metric_type& rt_dist,
                                   const point_type& a, double fraction,
                                   const point_type& b,
                                   double min_dist_interval,
                                   validity_predicate_type predicate) {

    if (a.time > b.time) {
      return a;
    }

    try {
      double dt_total =
          (b.time - a.time);  // the free time that I have along the path.
      if (dt_total < min_dist_interval) {
        return move_pt_toward(b_space, rt_dist, a, fraction, b);
      }

      generic_interpolator_impl<svp_interpolator, base_space_topo,
                                base_time_topo>
          interp;
      interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                        b_space.get_time_topology(), rt_dist);
      if (interp.get_minimum_travel_time() ==
          std::numeric_limits<double>::infinity()) {
        return a;
      }
      double dt = dt_total * fraction;
      double d = min_dist_interval;
      point_type result = a;
      point_type last_result = a;
      while (d < dt) {
        interp.compute_point(result.pt, a.pt, b.pt,
                             b_space.get_space_topology(),
                             b_space.get_time_topology(), d, dt_total, rt_dist);
        result.time = a.time + d;
        if (!predicate(result)) {
          return last_result;
        }
        d += min_dist_interval;
        last_result = result;
      };
      if (fraction == 1.0) {
        return b;
      }
      if (fraction == 0.0) {
        return a;
      }
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), dt, dt_total, rt_dist);
      result.time = a.time + dt;
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return a;
    }
  };

  static point_type move_pt_back_to(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, double fraction,
                                    const point_type& b) {

    if (a.time > b.time) {
      return b;
    }

    try {
      generic_interpolator_impl<svp_interpolator, base_space_topo,
                                base_time_topo>
          interp;
      double dt_total =
          (b.time - a.time);  // the free time that I have along the path.
      interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                        b_space.get_time_topology(), rt_dist);
      if (interp.get_minimum_travel_time() ==
          std::numeric_limits<double>::infinity()) {
        return b;
      }
      double dt = dt_total * (1.0 - fraction);
      point_type result = b;
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), dt, dt_total, rt_dist);
      result.time = a.time + dt;
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return b;
    }
  }

  static point_type move_pt_back_to(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, double fraction,
                                    const point_type& b,
                                    double min_dist_interval,
                                    validity_predicate_type predicate) {

    if (a.time > b.time) {
      return b;
    }

    try {
      double dt_total =
          (b.time - a.time);  // the free time that I have along the path.
      if (dt_total < min_dist_interval) {
        return move_pt_back_to(b_space, rt_dist, a, fraction, b);
      }

      generic_interpolator_impl<svp_interpolator, base_space_topo,
                                base_time_topo>
          interp;
      interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(),
                        b_space.get_time_topology(), rt_dist);
      if (interp.get_minimum_travel_time() ==
          std::numeric_limits<double>::infinity()) {
        return b;
      }
      double dt = dt_total * (1.0 - fraction);
      double d = dt_total - min_dist_interval;
      point_type result = b;
      point_type last_result = b;
      while (d > dt) {
        interp.compute_point(result.pt, a.pt, b.pt,
                             b_space.get_space_topology(),
                             b_space.get_time_topology(), d, dt_total, rt_dist);
        result.time = a.time + d;
        if (!predicate(result)) {
          return last_result;
        }
        d -= min_dist_interval;
        last_result = result;
      };
      if (fraction == 1.0) {
        return a;
      }
      if (fraction == 0.0) {
        return b;
      }
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(),
                           b_space.get_time_topology(), dt, dt_total, rt_dist);
      result.time = a.time + dt;
      return result;
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      return b;
    }
  }

  static double get_distance(const BaseTopology& b_space,
                             const rt_metric_type& rt_dist, const point_type& a,
                             const point_type& b) {
    if (a.time > b.time) {
      return std::numeric_limits<double>::infinity();
    }
    double reach_time = rt_dist(a.pt, b.pt, b_space.get_space_topology());
    if ((b.time - a.time) < reach_time) {
      return std::numeric_limits<double>::infinity();
    }
    return (b.time - a.time) + reach_time;
  }

  static double get_norm(const BaseTopology& b_space,
                         const rt_metric_type& rt_dist,
                         const point_difference_type& dp) {
    if (dp.time < 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    double reach_time = rt_dist(dp.pt, b_space.get_space_topology());
    if (dp.time < reach_time) {
      return std::numeric_limits<double>::infinity();
    }
    return dp.time + reach_time;
  }

  static double get_proper_distance(const BaseTopology& b_space,
                                    const rt_metric_type& rt_dist,
                                    const point_type& a, const point_type& b) {
    using std::abs;
    svp_reach_time_metric<base_time_topo, true> p_rt_dist(rt_dist);
    double reach_time = p_rt_dist(a.pt, b.pt, b_space.get_space_topology());
    return abs(b.time - a.time) + reach_time;
  }

  static double get_proper_norm(const BaseTopology& b_space,
                                const rt_metric_type& rt_dist,
                                const point_difference_type& dp) {
    using std::abs;
    svp_reach_time_metric<base_time_topo, true> p_rt_dist(rt_dist);
    double reach_time = p_rt_dist(dp.pt, b_space.get_space_topology());
    return abs(dp.time) + reach_time;
  }

  static bool is_in_bounds(const BaseTopology& b_space, const point_type& a) {
    return svp_is_in_bounds(a.pt, b_space.get_space_topology(),
                            b_space.get_time_topology());
  }

  static point_type random_point(const BaseTopology& b_space,
                                 const sampler_type& rl_sampler) {
    return point_type(get(random_sampler, b_space.get_time_topology())(
                          b_space.get_time_topology()),
                      rl_sampler(b_space.get_space_topology()));
  }
};

template <typename BaseTopology>
struct svp_reach_topo_selector {
  using type =
      svp_reach_topo_impl<BaseTopology, is_temporal_space_v<BaseTopology>>;
};

}  // namespace detail

/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and
 * sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class interpolated_topology<BaseTopology, svp_interpolation_tag>
    : public interpolated_topology_base<BaseTopology> {
 public:
  using base_type = interpolated_topology_base<BaseTopology>;
  using self = interpolated_topology<BaseTopology, svp_interpolation_tag>;

  using point_type = typename base_type::point_type;
  using point_difference_type = typename base_type::point_difference_type;

  using distance_metric_type = typename base_type::distance_metric_type;
  using random_sampler_type = typename base_type::random_sampler_type;

  static constexpr std::size_t dimensions = base_type::dimensions;

  using validity_predicate_type = std::function<bool(const point_type&)>;

 protected:
  using Impl = typename detail::svp_reach_topo_selector<BaseTopology>::type;
  using rt_metric_type = typename Impl::rt_metric_type;
  using sampler_type = typename Impl::sampler_type;

  rt_metric_type rt_dist;
  sampler_type rl_sampler;

  virtual point_type interp_topo_move_position_toward(
      const point_type& a, double fraction, const point_type& b) const {
    return Impl::move_pt_toward(*this, rt_dist, a, fraction, b);
  }

  virtual point_type interp_topo_move_position_toward_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval, validity_predicate_type predicate) const {
    return Impl::move_pt_toward(*this, rt_dist, a, fraction, b,
                                min_dist_interval, predicate);
  }

  virtual point_type interp_topo_move_position_back_to(
      const point_type& a, double fraction, const point_type& b) const {
    return Impl::move_pt_back_to(*this, rt_dist, a, fraction, b);
  }

  virtual point_type interp_topo_move_position_back_to_pred(
      const point_type& a, double fraction, const point_type& b,
      double min_dist_interval, validity_predicate_type predicate) const {
    return Impl::move_pt_back_to(*this, rt_dist, a, fraction, b,
                                 min_dist_interval, predicate);
  }

  virtual double interp_topo_get_distance(const point_type& a,
                                          const point_type& b) const {
    return Impl::get_distance(*this, rt_dist, a, b);
  }
  virtual double interp_topo_get_norm(const point_difference_type& dp) const {
    return Impl::get_norm(*this, rt_dist, dp);
  }
  virtual double interp_topo_get_proper_distance(const point_type& a,
                                                 const point_type& b) const {
    return Impl::get_proper_distance(*this, rt_dist, a, b);
  }
  virtual double interp_topo_get_proper_norm(
      const point_difference_type& dp) const {
    return Impl::get_proper_norm(*this, rt_dist, dp);
  }
  virtual bool interp_topo_is_in_bounds(const point_type& a) const {
    return Impl::is_in_bounds(*this, a);
  }
  virtual point_type interp_topo_random_point() const {
    return Impl::random_point(*this, rl_sampler);
  }

 public:
  const rt_metric_type& get_pseudo_factory() const { return rt_dist; }

  interpolated_topology(const BaseTopology& aTopo)
      : base_type(aTopo),
        rt_dist(Impl::make_rt_metric(*this)),
        rl_sampler(Impl::make_sampler(*this)) {}

  template <typename... Args>
  interpolated_topology(Args&&... args)
      : base_type(std::forward<Args>(args)...),
        rt_dist(Impl::make_rt_metric(*this)),
        rl_sampler(Impl::make_sampler(*this)) {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240003A, 1, "interpolated_topology",
                              base_type)
};

template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator<svp_interpolation_tag, SpaceType,
                                       TimeTopology> {
  using type = detail::generic_interpolator_impl<svp_interpolator, SpaceType,
                                                 TimeTopology>;
  using pseudo_factory_type = svp_reach_time_metric<TimeTopology>;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator<svp_interpolation_tag,
                                        TemporalSpaceType> {
  using type = generic_interpolator<svp_interpolator_factory<TemporalSpaceType>,
                                    svp_interpolator>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SVP_REACH_TOPOLOGIES_H_
