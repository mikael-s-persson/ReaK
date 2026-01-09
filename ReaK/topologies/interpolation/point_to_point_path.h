/**
 * \file point_to_point_path.h
 *
 * This library provides an implementation of a linear-interpolated path within a topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a linear interpolation based on the distance.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_POINT_TO_POINT_PATH_H_
#define REAK_TOPOLOGIES_INTERPOLATION_POINT_TO_POINT_PATH_H_

#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/topologies/interpolation/spatial_path_concept.h"

#include "ReaK/topologies/interpolation/waypoint_container.h"

#include <cmath>
#include <limits>
#include <list>
#include <map>

namespace ReaK::pp {

/**
 * This class implements a point-to-point path within a topology (interpolated by whatever method is used in
 * the topology's move function). The path is represented by a set of waypoints and all intermediate points
 * are computed with the topology's move_position_toward function based on the distance or fraction of travel.
 * This class models the SpatialPathConcept and SequentialPathConcept.
 * \tparam Space The topology type on which the points and the path can reside.
 * \tparam Metric The distance metric used to assess the distance between points in the path.
 */
template <MetricSpace Space,
          DistanceMetric<Space> Metric =
              typename metric_space_traits<Space>::distance_metric_type>
class point_to_point_path : public waypoint_container<Space, Metric> {
 public:
  using self = point_to_point_path<Space, Metric>;
  using base_class_type = waypoint_container<Space, Metric>;

  using const_waypoint_descriptor =
      typename base_class_type::const_waypoint_descriptor;
  using const_waypoint_bounds = typename base_class_type::const_waypoint_bounds;
  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  using waypoint_pair = typename base_class_type::waypoint_pair;

 private:
  using container_type = typename base_class_type::container_type;

  const container_type& get_waypoints() const { return this->waypoints; }
  const distance_metric& get_dist() const { return this->dist; }
  const topology& get_space() const { return *(this->space); }

  double travel_distance_impl(const point_type& a,
                              const const_waypoint_bounds& wpb_a,
                              const point_type& b,
                              const const_waypoint_bounds& wpb_b) const {
    if (std::distance(wpb_a.first, this->waypoints.end()) <
        std::distance(wpb_b.first, this->waypoints.end())) {
      return travel_distance_impl(b, wpb_b, a, wpb_a);
    }

    double sum = 0;
    if (((wpb_a.first == wpb_b.first) && (wpb_a.second == wpb_b.second)) ||
        (wpb_a.second == wpb_b.first)) {
      // this means that a and b are in the same segment.
      return this->dist(a, b, *(this->space));
    }

    sum += this->dist(a, *(wpb_a.second), *(this->space));

    const_waypoint_descriptor it = wpb_a.second;
    auto it_prev = it;
    while (++it != wpb_b.first) {
      sum += this->dist(*it_prev++, *it, *(this->space));
    }

    sum += this->dist(*it_prev, *(wpb_b.first), *(this->space));

    sum += this->dist(*(wpb_b.first), b, *(this->space));

    return sum;
  }

  waypoint_pair get_point_at_distance_impl(
      double d, const const_waypoint_bounds& wpb_a) const {
    const_waypoint_descriptor it_prev = wpb_a.first;
    const_waypoint_descriptor it = it_prev;
    ++it;

    if (it == this->waypoints.end()) {
      return waypoint_pair(it_prev, *it_prev);
    }

    double total_d = this->dist(*it_prev, *it, *(this->space));

    return waypoint_pair(
        it_prev, this->space->move_position_toward(*it_prev, d / total_d, *it));
  }

  waypoint_pair move_away_from_impl(point_type a, const_waypoint_bounds wpb_a,
                                    double d) const {
    while (true) {
      double d_to_second = this->dist(a, *(wpb_a.second), *(this->space));
      if (d_to_second >= d) {
        return waypoint_pair(
            wpb_a.first, this->space->move_position_toward(a, d / d_to_second,
                                                           *(wpb_a.second)));
      }
      const_waypoint_descriptor it = wpb_a.second;
      ++it;
      if (it == this->waypoints.end()) {
        return waypoint_pair(wpb_a.first, *(wpb_a.second));
      }
      d -= d_to_second;
      a = *(wpb_a.second);
      wpb_a.first = wpb_a.second;
      wpb_a.second = it;
    }
  }

 public:
  struct point_distance_iterator {
    const point_to_point_path* parent;
    const_waypoint_bounds current_wpbound;
    point_type current_pt;
    double segment_distance;
    double current_distance;

    point_distance_iterator(const point_to_point_path* aParent,
                            const const_waypoint_bounds& aWPB,
                            const point_type& aCurrentPt)
        : parent(aParent), current_wpbound(aWPB), current_pt(aCurrentPt) {
      if (current_wpbound.first == current_wpbound.second) {
        segment_distance = 0.0;
        current_distance = 0.0;
        return;
      }
      segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      double cur_d_back = parent->get_dist()(*(current_wpbound.first),
                                             current_pt, parent->get_space());
      double cur_d_forth = parent->get_dist()(
          current_pt, *(current_wpbound.second), parent->get_space());
      current_distance = (cur_d_back + (segment_distance - cur_d_forth)) * 0.5;
      if (current_distance < 0.0) {
        current_distance = cur_d_back;
      }
    }

    point_distance_iterator(const point_to_point_path* aParent,
                            const const_waypoint_bounds& aWPB)
        : parent(aParent), current_wpbound(aWPB), current_pt(*(aWPB.first)) {
      segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      current_distance = 0.0;
    }

    point_distance_iterator(const point_to_point_path* aParent,
                            const_waypoint_bounds&& aWPB,
                            point_type&& aCurrentPt)
        : parent(aParent),
          current_wpbound(std::move(aWPB)),
          current_pt(std::move(aCurrentPt)) {
      if (current_wpbound.first == current_wpbound.second) {
        segment_distance = 0.0;
        current_distance = 0.0;
        return;
      }
      segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      double cur_d_back = parent->get_dist()(*(current_wpbound.first),
                                             current_pt, parent->get_space());
      double cur_d_forth = parent->get_dist()(
          current_pt, *(current_wpbound.second), parent->get_space());
      current_distance = (cur_d_back + (segment_distance - cur_d_forth)) * 0.5;
      if (current_distance < 0.0) {
        current_distance = cur_d_back;
      }
    }

    point_distance_iterator(const point_to_point_path* aParent,
                            const_waypoint_bounds&& aWPB)
        : parent(aParent),
          current_wpbound(std::move(aWPB)),
          current_pt(*(current_wpbound.first)) {
      segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      current_distance = 0.0;
    }

    point_distance_iterator& operator+=(double rhs) {
      if (current_wpbound.first == current_wpbound.second) {
        // then we are at some edge of the set of waypoints.
        if (current_wpbound.first == parent->get_waypoints().begin()) {
          // we are at the beginning.
          if (rhs < 0.0) {
            return *this;  // cannot move further back.
          }
          ++(current_wpbound.second);
          if (current_wpbound.second == parent->get_waypoints().end()) {
            current_wpbound.second = current_wpbound.first;
            return *this;  // nothing beyond the first point.
          }
          segment_distance = parent->get_dist()(*(current_wpbound.first),
                                                *(current_wpbound.second),
                                                parent->get_space());
        } else {
          // we are at the end.
          if (rhs > 0.0) {
            return *this;  // cannot move further forward.
          }
          --(current_wpbound.first);
          segment_distance = parent->get_dist()(*(current_wpbound.first),
                                                *(current_wpbound.second),
                                                parent->get_space());
          current_distance = segment_distance;
        }
      }

      if (rhs < 0.0) {
        while (current_distance + rhs <
               0.0) {  // while spilling over before the current bounds.
          if (current_wpbound.first == parent->get_waypoints().begin()) {
            current_pt = *(current_wpbound.first);
            current_wpbound.second = current_wpbound.first;
            segment_distance = 0.0;
            current_distance = 0.0;
            return *this;
          }
          --(current_wpbound.first);
          --(current_wpbound.second);
          rhs += current_distance;
          segment_distance = parent->get_dist()(*(current_wpbound.first),
                                                *(current_wpbound.second),
                                                parent->get_space());
          current_distance = segment_distance;
        }
      } else {
        while (segment_distance <=
               current_distance +
                   rhs) {  // while spilling over beyond the current bounds.
          ++(current_wpbound.first);
          ++(current_wpbound.second);
          rhs -= segment_distance - current_distance;
          current_distance = 0.0;
          if (current_wpbound.second == parent->get_waypoints().end()) {
            current_pt = *(current_wpbound.first);
            current_wpbound.second = current_wpbound.first;
            segment_distance = 0.0;
            current_distance = 0.0;
            return *this;
          }
          segment_distance = parent->get_dist()(*(current_wpbound.first),
                                                *(current_wpbound.second),
                                                parent->get_space());
        }
      }

      double frac = (current_distance + rhs) / segment_distance;
      current_pt = parent->get_space().move_position_toward(
          *(current_wpbound.first), frac, *(current_wpbound.second));
      current_distance += rhs;
      return *this;
    }

    friend point_distance_iterator operator+(point_distance_iterator lhs,
                                             double rhs) {
      return (lhs += rhs);
    }

    friend point_distance_iterator operator+(double lhs,
                                             point_distance_iterator rhs) {
      return (rhs += lhs);
    }

    friend point_distance_iterator operator-(point_distance_iterator lhs,
                                             double rhs) {
      return (lhs += -rhs);
    }

    friend point_distance_iterator& operator-=(point_distance_iterator& lhs,
                                               double rhs) {
      return (lhs += -rhs);
    }

    friend bool operator==(const point_distance_iterator& lhs,
                           const point_distance_iterator& rhs) {
      using std::abs;
      return static_cast<bool>(
          (lhs.parent == rhs.parent) &&
          (lhs.current_wpbound.first == rhs.current_wpbound.first) &&
          (lhs.current_wpbound.second == rhs.current_wpbound.second) &&
          (abs(lhs.current_distance - rhs.current_distance) <=
           std::numeric_limits<double>::epsilon()));
    }

    friend bool operator!=(const point_distance_iterator& lhs,
                           const point_distance_iterator& rhs) {
      return !(lhs == rhs);
    }

    const point_type& operator*() const { return current_pt; }
  };

  struct point_fraction_iterator {
    const point_to_point_path* parent;
    const_waypoint_bounds current_wpbound;
    point_type current_pt;
    double current_fraction;

    point_fraction_iterator(const point_to_point_path* aParent,
                            const const_waypoint_bounds& aWPB,
                            const point_type& aCurrentPt)
        : parent(aParent), current_wpbound(aWPB), current_pt(aCurrentPt) {
      if (current_wpbound.first == current_wpbound.second) {
        current_fraction = 0.0;
        return;
      }
      double segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      double cur_d_back = parent->get_dist()(*(current_wpbound.first),
                                             current_pt, parent->get_space());
      double cur_d_forth = parent->get_dist()(
          current_pt, *(current_wpbound.second), parent->get_space());
      double current_distance =
          (cur_d_back + (segment_distance - cur_d_forth)) * 0.5;
      if (current_distance < 0.0) {
        current_distance = cur_d_back;
      }
      if (current_distance > segment_distance) {
        current_distance = segment_distance;
      }
      current_fraction = current_distance / segment_distance;
    }

    point_fraction_iterator(const point_to_point_path* aParent,
                            const const_waypoint_bounds& aWPB)
        : parent(aParent),
          current_wpbound(aWPB),
          current_pt(*(aWPB.first)),
          current_fraction(0.0) {}

    point_fraction_iterator(const point_to_point_path* aParent,
                            const_waypoint_bounds&& aWPB,
                            point_type&& aCurrentPt)
        : parent(aParent),
          current_wpbound(std::move(aWPB)),
          current_pt(std::move(aCurrentPt)) {
      if (current_wpbound.first == current_wpbound.second) {
        current_fraction = 0.0;
        return;
      }
      double segment_distance =
          parent->get_dist()(*(current_wpbound.first),
                             *(current_wpbound.second), parent->get_space());
      double cur_d_back = parent->get_dist()(*(current_wpbound.first),
                                             current_pt, parent->get_space());
      double cur_d_forth = parent->get_dist()(
          current_pt, *(current_wpbound.second), parent->get_space());
      double current_distance =
          (cur_d_back + (segment_distance - cur_d_forth)) * 0.5;
      if (current_distance < 0.0) {
        current_distance = cur_d_back;
      }
      if (current_distance > segment_distance) {
        current_distance = segment_distance;
      }
      current_fraction = current_distance / segment_distance;
    }

    point_fraction_iterator(const point_to_point_path* aParent,
                            const_waypoint_bounds&& aWPB)
        : parent(aParent),
          current_wpbound(std::move(aWPB)),
          current_pt(*(current_wpbound.first)),
          current_fraction(0.0) {}

    point_fraction_iterator& operator+=(double rhs) {
      if (current_wpbound.first == current_wpbound.second) {
        // then we are at some edge of the set of waypoints.
        if (current_wpbound.first == parent->get_waypoints().begin()) {
          // we are at the beginning.
          if (rhs < 0.0) {
            return *this;  // cannot move further back.
          }
          ++(current_wpbound.second);
          if (current_wpbound.second == parent->get_waypoints().end()) {
            current_wpbound.second = current_wpbound.first;
            return *this;  // nothing beyond the first point.
          }
        } else {
          // we are at the end.
          if (rhs > 0.0) {
            return *this;  // cannot move further forward.
          }
          --(current_wpbound.first);
          current_fraction = 1.0;
        }
      }

      if (rhs < 0.0) {
        while (current_fraction + rhs <
               0.0) {  // while spilling over before the current bounds.
          if (current_wpbound.first == parent->get_waypoints().begin()) {
            current_pt = *(current_wpbound.first);
            current_wpbound.second = current_wpbound.first;
            current_fraction = 0.0;
            return *this;
          }
          --(current_wpbound.first);
          --(current_wpbound.second);
          rhs += current_fraction;
          current_fraction = 1.0;
        }
      } else {
        while (1.0 <=
               current_fraction +
                   rhs) {  // while spilling over beyond the current bounds.
          ++(current_wpbound.first);
          ++(current_wpbound.second);
          rhs -= 1.0 - current_fraction;
          current_fraction = 0.0;
          if (current_wpbound.second == parent->get_waypoints().end()) {
            current_pt = *(current_wpbound.first);
            current_wpbound.second = current_wpbound.first;
            current_fraction = 0.0;
            return *this;
          }
        }
      }

      current_pt = parent->get_space().move_position_toward(
          *(current_wpbound.first), current_fraction + rhs,
          *(current_wpbound.second));
      current_fraction += rhs;
      return *this;
    }

    friend point_fraction_iterator operator+(point_fraction_iterator lhs,
                                             double rhs) {
      return (lhs += rhs);
    }

    friend point_fraction_iterator operator+(double lhs,
                                             point_fraction_iterator rhs) {
      return (rhs += lhs);
    }

    friend point_fraction_iterator operator-(point_fraction_iterator lhs,
                                             double rhs) {
      return (lhs += -rhs);
    }

    friend point_fraction_iterator& operator-=(point_fraction_iterator& lhs,
                                               double rhs) {
      return (lhs += -rhs);
    }

    friend bool operator==(const point_fraction_iterator& lhs,
                           const point_fraction_iterator& rhs) {
      using std::abs;
      return static_cast<bool>(
          (lhs.parent == rhs.parent) &&
          (lhs.current_wpbound.first == rhs.current_wpbound.first) &&
          (lhs.current_wpbound.second == rhs.current_wpbound.second) &&
          (abs(lhs.current_fraction - rhs.current_fraction) <=
           std::numeric_limits<double>::epsilon()));
    }

    friend bool operator!=(const point_fraction_iterator& lhs,
                           const point_fraction_iterator& rhs) {
      return !(lhs == rhs);
    }

    const point_type& operator*() const { return current_pt; }
  };

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit point_to_point_path(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  point_to_point_path(const std::shared_ptr<topology>& aSpace,
                      const point_type& aStart, const point_type& aEnd,
                      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the on-past-last point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  point_to_point_path(ForwardIter aBegin, ForwardIter aEnd,
                      const std::shared_ptr<topology>& aSpace,
                      const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(static_cast<base_class_type&>(lhs),
         static_cast<base_class_type&>(rhs));
  }

  /**
   * Returns the starting distance-iterator along the path.
   * \return The starting distance-iterator along the path.
   */
  point_distance_iterator begin_distance_travel() const {
    return point_distance_iterator(
        this, const_waypoint_bounds(this->waypoints.begin(),
                                    this->waypoints.begin()));
  }

  /**
   * Returns the end distance-iterator along the path.
   * \return The end distance-iterator along the path.
   */
  point_distance_iterator end_distance_travel() const {
    auto it = this->waypoints.end();
    --it;
    return point_distance_iterator(this, const_waypoint_bounds(it, it));
  }

  /**
   * Returns the starting fraction-iterator along the path.
   * \return The starting fraction-iterator along the path.
   */
  point_fraction_iterator begin_fraction_travel() const {
    return point_fraction_iterator(
        this, const_waypoint_bounds(this->waypoints.begin(),
                                    this->waypoints.begin()));
  }

  /**
   * Returns the end fraction-iterator along the path.
   * \return The end fraction-iterator along the path.
   */
  point_fraction_iterator end_fraction_travel() const {
    auto it = this->waypoints.end();
    --it;
    return point_fraction_iterator(this, const_waypoint_bounds(it, it));
  }

  /**
   * Computes the travel distance between two points, if traveling along the path.
   * \param a The first point.
   * \param b The second point.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance(const point_type& a, const point_type& b) const {
    const_waypoint_bounds wpb_a =
        this->get_waypoint_bounds(a, this->waypoints.begin());
    const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b, wpb_a.first);
    return this->travel_distance_impl(a, wpb_a, b, wpb_b);
  }

  /**
   * Computes the travel distance between two waypoint-point-pairs, if traveling along the path.
   * \param a The first waypoint-point-pair.
   * \param b The second waypoint-point-pair.
   * \return The travel distance between two points if traveling along the path.
   */
  double travel_distance(waypoint_pair& a, waypoint_pair& b) const {
    const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
    const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b.second, b.first);
    a.first = wpb_a.first;
    b.first = wpb_b.first;
    return this->travel_distance_impl(a.second, wpb_a, b.second, wpb_b);
  }

  /**
   * Computes the point that is a distance away from a point on the path.
   * \param a The point on the path.
   * \param d The distance to move away from the point.
   * \return The point that is a distance away from the given point.
   */
  point_type move_away_from(const point_type& a, double d) const {
    const_waypoint_bounds wpb_a =
        this->get_waypoint_bounds(a, this->waypoints.begin());
    return move_away_from_impl(a, wpb_a, d).second;
  }

  /**
   * Computes the waypoint-point-pair that is a distance away from a waypoint-point-pair on the path.
   * \param a The waypoint-point-pair on the path.
   * \param dt The distance to move away from the waypoint-point-pair.
   * \return The waypoint-point-pair that is a time away from the given waypoint-point-pair.
   */
  waypoint_pair move_away_from(const waypoint_pair& a, double d) const {
    const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.second, a.first);
    return move_away_from_impl(a.second, wpb_a, d);
  }

  /**
   * Computes the point that is on the path at the given distance from the start.
   * \param d The distance (from start) at which the point is sought.
   * \return The point that is on the path at the given distance from the start.
   */
  point_type get_point_at_distance(double d) const {
    const_waypoint_descriptor start = this->waypoints.begin();
    return move_away_from_impl(*start, const_waypoint_bounds(start, start), d)
        .second;
  }

  /**
   * Computes the waypoint-point pair that is on the path at the given distance from the start.
   * \param t The distance from the start at which the waypoint-point pair is sought.
   * \return The waypoint-point pair that is on the trajectory at the given distance from the start.
   */
  waypoint_pair get_waypoint_at_distance(double d) const {
    const_waypoint_descriptor start = this->waypoints.begin();
    return move_away_from_impl(*start, const_waypoint_bounds(start, start), d);
  }

  /**
   * Returns the total length of the path.
   * \return The total length of the path.
   */
  double get_total_length() const {
    const_waypoint_descriptor start = this->waypoints.begin();
    const_waypoint_descriptor end = this->waypoints.end();
    --end;
    return travel_distance_impl(*start, const_waypoint_bounds(start, start),
                                *end, const_waypoint_bounds(end, end));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_class_type::save(A,
                          base_class_type::get_static_object_type()->version());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_class_type::load(A,
                          base_class_type::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC244000B, 1, "point_to_point_path",
                              base_class_type)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_POINT_TO_POINT_PATH_H_
