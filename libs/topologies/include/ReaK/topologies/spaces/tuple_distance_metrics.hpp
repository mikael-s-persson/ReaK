/**
 * \file tuple_distance_metrics.hpp
 *
 * This library defines the distance metrics classes that work on points of a topology-tuple.
 * All the classes satisfy the DistanceMetricConcept and are suitable for topology-tuple (see metric_space_tuple).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_TUPLE_DISTANCE_METRICS_HPP
#define REAK_TUPLE_DISTANCE_METRICS_HPP

#include <ReaK/core/base/serializable.hpp>

#include <cmath>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

#include "metric_space_concept.hpp"

namespace ReaK {

namespace pp {


namespace detail {
namespace {

struct manhattan_tuple_distance_accum {
  double* p_result;
  manhattan_tuple_distance_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Pt >
  void operator()( const Space& s, const Pt& p1, const Pt& p2 ) const {
    using std::fabs;
    ( *p_result ) += fabs( get( distance_metric, s )( p1, p2, s ) );
  };
};

struct manhattan_tuple_norm_accum {
  double* p_result;
  manhattan_tuple_norm_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Dp >
  void operator()( const Space& s, const Dp& dp ) const {
    using std::fabs;
    ( *p_result ) += fabs( get( distance_metric, s )( dp, s ) );
  };
};


struct euclidean_tuple_distance_accum {
  double* p_result;
  euclidean_tuple_distance_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Pt >
  void operator()( const Space& s, const Pt& p1, const Pt& p2 ) const {
    double r = get( distance_metric, s )( p1, p2, s );
    ( *p_result ) += r * r;
  };
};

struct euclidean_tuple_norm_accum {
  double* p_result;
  euclidean_tuple_norm_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Dp >
  void operator()( const Space& s, const Dp& dp ) const {
    double r = get( distance_metric, s )( dp, s );
    ( *p_result ) += r * r;
  };
};


struct p_norm_tuple_distance_accum {
  int p;
  double* p_result;
  p_norm_tuple_distance_accum( int p_value, double& result ) : p( p_value ), p_result( &result ){};
  template < typename Space, typename Pt >
  void operator()( const Space& s, const Pt& p1, const Pt& p2 ) const {
    using std::pow;
    using std::fabs;
    double r = fabs( get( distance_metric, s )( p1, p2, s ) );
    ( *p_result ) += pow( r, p );
  };
};

struct p_norm_tuple_norm_accum {
  int p;
  double* p_result;
  p_norm_tuple_norm_accum( int p_value, double& result ) : p( p_value ), p_result( &result ){};
  template < typename Space, typename Dp >
  void operator()( const Space& s, const Dp& dp ) const {
    using std::pow;
    using std::fabs;
    double r = fabs( get( distance_metric, s )( dp, s ) );
    ( *p_result ) += pow( r, p );
  };
};


struct inf_norm_tuple_distance_accum {
  double* p_result;
  inf_norm_tuple_distance_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Pt >
  void operator()( const Space& s, const Pt& p1, const Pt& p2 ) const {
    using std::fabs;
    double r = fabs( get( distance_metric, s )( p1, p2, s ) );
    if( ( *p_result ) < r )
      ( *p_result ) = r;
  };
};

struct inf_norm_tuple_norm_accum {
  double* p_result;
  inf_norm_tuple_norm_accum( double& result ) : p_result( &result ){};
  template < typename Space, typename Dp >
  void operator()( const Space& s, const Dp& dp ) const {
    using std::fabs;
    double r = fabs( get( distance_metric, s )( dp, s ) );
    if( ( *p_result ) < r )
      ( *p_result ) = r;
  };
};
};
};


/**
 * This class is a Manhattan distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Manhattan norm to the point-difference vectors in the
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct manhattan_tuple_distance : public serializable {

  manhattan_tuple_distance(){};

  template < typename Topology >
  manhattan_tuple_distance( const Topology& ){};

  /**
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template < typename PointDiff, typename Topology >
  double operator()( const PointDiff& a, const Topology& s ) const {
    double result = 0.0;
    tuple_for_each( s, a, detail::manhattan_tuple_norm_accum( result ) );
    return result;
  };

  /**
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template < typename Point, typename Topology >
  double operator()( const Point& a, const Point& b, const Topology& s ) const {
    double result = 0.0;
    tuple_for_each( s, a, b, detail::manhattan_tuple_distance_accum( result ) );
    return result;
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( manhattan_tuple_distance, 0xC2410005, 1, "manhattan_tuple_distance", serializable )
};


/**
 * This class is a Euclidean distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Euclidean norm to the point-difference vectors in the
 * given topology (assuming it models the MetricSpaceConcept and that it is a tuple).
 */
struct euclidean_tuple_distance : public serializable {

  euclidean_tuple_distance(){};

  template < typename Topology >
  euclidean_tuple_distance( const Topology& ){};

  /**
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology-tuple.
   */
  template < typename PointDiff, typename Topology >
  double operator()( const PointDiff& a, const Topology& s ) const {
    using std::sqrt;
    double result = 0.0;
    tuple_for_each( s, a, detail::euclidean_tuple_norm_accum( result ) );
    return sqrt( result );
  };

  /**
   * This function returns the distance between two points on a topology-tuple.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology-tuple.
   */
  template < typename Point, typename Topology >
  double operator()( const Point& a, const Point& b, const Topology& s ) const {
    using std::sqrt;
    double result = 0.0;
    tuple_for_each( s, a, b, detail::euclidean_tuple_distance_accum( result ) );
    return sqrt( result );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( euclidean_tuple_distance, 0xC2410006, 1, "euclidean_tuple_distance", serializable )
};


/**
 * This class is a Infinity-norm distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Infinity-norm to the point-difference vectors in the
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct inf_norm_tuple_distance : public serializable {

  inf_norm_tuple_distance(){};

  template < typename Topology >
  inf_norm_tuple_distance( const Topology& ){};

  /**
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \return The norm of the difference between two points on a topology-tuple.
   */
  template < typename PointDiff, typename Topology >
  double operator()( const PointDiff& a, const Topology& s ) const {
    double result = 0.0;
    tuple_for_each( s, a, detail::inf_norm_tuple_norm_accum( result ) );
    return result;
  };

  /**
   * This function returns the distance between two points on a topology-tuple.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology-tuple.
   */
  template < typename Point, typename Topology >
  double operator()( const Point& a, const Point& b, const Topology& s ) const {
    double result = 0.0;
    tuple_for_each( s, a, b, detail::inf_norm_tuple_distance_accum( result ) );
    return result;
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_MAKE_ABSTRACT_1BASE( inf_norm_tuple_distance, 0xC2410007, 1, "inf_norm_tuple_distance", serializable )
};


/**
 * This class is a Euclidean distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Euclidean norm to the point-difference vectors in the
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct p_norm_tuple_distance : public serializable {

  int p_value;

  p_norm_tuple_distance( int aP = 2 ) : p_value( aP ){};

  template < typename Topology >
  p_norm_tuple_distance( const Topology&, int aP = 2 )
      : p_value( aP ){};

  /**
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points  on a topology-tuple.
   */
  template < typename PointDiff, typename Topology >
  double operator()( const PointDiff& a, const Topology& s ) const {
    using std::pow;
    double result = 0.0;
    tuple_for_each( s, a, detail::p_norm_tuple_norm_accum( p_value, result ) );
    return pow( result, 1.0 / p_value );
  };

  /**
   * This function returns the distance between two points on a topology-tuple.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology-tuple.
   */
  template < typename Point, typename Topology >
  double operator()( const Point& a, const Point& b, const Topology& s ) const {
    using std::pow;
    double result = 0.0;
    tuple_for_each( s, a, b, detail::p_norm_tuple_distance_accum( p_value, result ) );
    return pow( result, 1.0 / p_value );
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( p_value );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) { A& RK_SERIAL_LOAD_WITH_NAME( p_value ); };

  RK_RTTI_MAKE_ABSTRACT_1BASE( p_norm_tuple_distance, 0xC2410008, 1, "p_norm_tuple_distance", serializable )
};
};
};


#endif
