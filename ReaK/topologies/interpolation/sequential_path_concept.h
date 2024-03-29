/**
 * \file sequential_path_concept.h
 *
 * This library defines the traits and concepts related to a sequential spatial path. A
 * path is simply a continuous curve in a topology (or space) which can be travelled
 * sequentially via either increments in distance or in fractions between waypoints.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_PATH_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_PATH_CONCEPT_H_

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize a sequential spatial path within a
 * topology.
 * \tparam SequentialPath The spatial path type for which the traits are sought.
 */
template <typename Path>
struct sequential_path_traits {
  /** This type describes a point in the space or topology. */
  using point_type = typename Path::point_type;

  /** This type describes an iterator, corresponding to a point on the path, which can be incremented by distance to
   * travel to the next iterator. */
  using point_distance_iterator = typename Path::point_distance_iterator;
  /** This type describes an iterator, corresponding to a point on the path, which can be incremented by a fraction
   * between waypoints to travel to the next iterator. */
  using point_fraction_iterator = typename Path::point_fraction_iterator;

  /** This type is the topology type in which the path exists. */
  using topology = typename Path::topology;
  /** This type is the distance metric type used on the topology and defining the travel distances along the path. */
  using distance_metric = typename Path::distance_metric;
};

/**
 * This concept class defines the requirements for a type to model a sequential spatial-path
 * as used in ReaK::pp. A sequential spatial path is a continuous curve within a topology
 * which can be travelled sequentially via either increments in distance or in fractions
 * between waypoints.
 *
 * Required concepts:
 *
 * The topology should model the MetricSpace.
 *
 * The distance-metric should model the DistanceMetric.
 *
 * Valid expressions:
 *
 * dit = path.begin_distance_travel();  The start of the distance-iterator range of the sequential path can be obtained.
 *
 * dit = path.end_distance_travel();  The end of the distance-iterator range (one-past-last) of the sequential path can
 *be obtained.
 *
 * pt = *dit;  A point can be obtained from dereferencing a distance-iterator.
 *
 * dit = dit + d;
 * dit = d + dit;
 * dit += d;
 * dit = dit - d;
 * dit -= d;  A distance-iterator can be incremented by a distance (double).
 *
 * b = (dit != dit);
 * b = (dit == dit);  Two distance-iterator can be compared for inequality.
 *
 * fit = path.begin_fraction_travel();  The start of the fraction-iterator range of the sequential path can be obtained.
 *
 * fit = path.end_fraction_travel();  The end of the fraction-iterator range (one-past-last) of the sequential path can
 *be obtained.
 *
 * pt = *fit;  A point can be obtained from dereferencing a fraction-iterator.
 *
 * fit = fit + f;
 * fit = f + fit;
 * fit += f;
 * fit = fit - f;
 * fit -= f;  A fraction-iterator can be incremented by a fraction (double).
 *
 * b = (fit != fit);
 * b = (fit == fit);  Two fraction-iterator can be compared for equality.
 *
 * d = path.travel_distance(pt,pt);  The travel distance, along the path (p), between two points (pt,pt), can be
 *obtained.
 */
template <typename Path, typename Space>
concept SequentialPath = Topology<Space>&& DistanceMetric<
    typename sequential_path_traits<Path>::distance_metric, Space>&&
requires(const Path& path, const topology_point_type_t<Space>& pt) {
  { path.travel_distance(pt, pt) } -> std::convertible_to<double>;
}
&&requires(const Path& path) {
  {
    path.begin_distance_travel()
    } -> std::convertible_to<
        typename sequential_path_traits<Path>::point_distance_iterator>;
  {
    path.end_distance_travel()
    } -> std::convertible_to<
        typename sequential_path_traits<Path>::point_distance_iterator>;
}
&&requires(
    double d,
    typename sequential_path_traits<Path>::point_distance_iterator& dit) {
  { *dit } -> std::convertible_to<topology_point_type_t<Space>>;
  dit = dit + d;
  dit = d + dit;
  dit += d;
  dit = dit - d;
  dit -= d;
  { dit != dit } -> std::convertible_to<bool>;
  { dit == dit } -> std::convertible_to<bool>;
}
&&requires(const Path& path) {
  {
    path.begin_fraction_travel()
    } -> std::convertible_to<
        typename sequential_path_traits<Path>::point_fraction_iterator>;
  {
    path.end_fraction_travel()
    } -> std::convertible_to<
        typename sequential_path_traits<Path>::point_fraction_iterator>;
}
&&requires(
    double d,
    typename sequential_path_traits<Path>::point_fraction_iterator& fit) {
  { *fit } -> std::convertible_to<topology_point_type_t<Space>>;
  fit = fit + d;
  fit = d + fit;
  fit += d;
  fit = fit - d;
  fit -= d;
  { fit != fit } -> std::convertible_to<bool>;
  { fit == fit } -> std::convertible_to<bool>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_PATH_CONCEPT_H_
