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
#include <ReaK/core/lin_alg/arithmetic_tuple.hpp>

#include <ReaK/ctrl/path_planning/metric_space_concept.hpp>

namespace ReaK {

namespace pp {
  

namespace detail {

  template <std::size_t Idx, typename SpaceTuple> 
  struct manhattan_tuple_distance_impl {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      double result = manhattan_tuple_distance_impl<Idx-1,SpaceTuple>::distance(s,p1,p2);
      using std::fabs;
      result += fabs(get(distance_metric,get<Idx>(s))(get<Idx>(p1), get<Idx>(p2), get<Idx>(s)));
      return result;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      double result = manhattan_tuple_distance_impl<Idx-1,SpaceTuple>::norm(s,dp);
      using std::fabs;
      result += fabs(get(distance_metric,get<Idx>(s))(get<Idx>(dp), get<Idx>(s)));
      return result;
    };
  };
  
  template <typename SpaceTuple> 
  struct manhattan_tuple_distance_impl<0,SpaceTuple> {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      using std::fabs;
      return fabs(get(distance_metric,get<0>(s))(get<0>(p1), get<0>(p2), get<0>(s)));
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      using std::fabs;
      return fabs(get(distance_metric,get<0>(s))(get<0>(dp), get<0>(s)));
    };
  };
  
  
  template <std::size_t Idx, typename SpaceTuple> 
  struct euclidean_tuple_distance_impl {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      double result = euclidean_tuple_distance_impl<Idx-1,SpaceTuple>::distance(s,p1,p2);
      double r = get(distance_metric,get<Idx>(s))(get<Idx>(p1), get<Idx>(p2), get<Idx>(s));
      result += r * r;
      return result;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      double result = euclidean_tuple_distance_impl<Idx-1,SpaceTuple>::norm(s,dp);
      double r = get(distance_metric,get<Idx>(s))(get<Idx>(dp), get<Idx>(s));
      result += r * r;
      return result;
    };
  };
  
  template <typename SpaceTuple> 
  struct euclidean_tuple_distance_impl<0,SpaceTuple> {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      double r = get(distance_metric,get<0>(s))(get<0>(p1), get<0>(p2), get<0>(s));
      return r * r;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      double r = get(distance_metric,get<0>(s))(get<0>(dp), get<0>(s));
      return r * r;
    };
  };
  
  
  template <std::size_t Idx, typename SpaceTuple> 
  struct p_norm_tuple_distance_impl {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2, int p_value) {
      double result = p_norm_tuple_distance_impl<Idx-1,SpaceTuple>::distance(s,p1,p2,p_value);
      double r = get(distance_metric,get<Idx>(s))(get<Idx>(p1), get<Idx>(p2), get<Idx>(s));
      using std::pow;
      result += pow(r,p_value);
      return result;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp, int p_value) {
      double result = p_norm_tuple_distance_impl<Idx-1,SpaceTuple>::norm(s,dp,p_value);
      double r = get(distance_metric,get<Idx>(s))(get<Idx>(dp), get<Idx>(s));
      using std::pow;
      result += pow(r,p_value);
      return result;
    };
  };
  
  template <typename SpaceTuple> 
  struct p_norm_tuple_distance_impl<0,SpaceTuple> {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2, int p_value) {
      double r = get(distance_metric,get<0>(s))(get<0>(p1), get<0>(p2), get<0>(s));
      using std::pow;
      return pow(r,p_value);
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp, int p_value) {
      double r = get(distance_metric,get<0>(s))(get<0>(dp), get<0>(s));
      using std::pow;
      return pow(r,p_value);
    };
  };
  
  
  template <std::size_t Idx, typename SpaceTuple> 
  struct inf_norm_tuple_distance_impl {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      double result = inf_norm_tuple_distance_impl<Idx-1,SpaceTuple>::distance(s,p1,p2);
      using std::fabs;
      double r = fabs(get(distance_metric,get<Idx>(s))(get<Idx>(p1), get<Idx>(p2), get<Idx>(s)));
      if(result < r)
        return r;
      else
        return result;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      double result = inf_norm_tuple_distance_impl<Idx-1,SpaceTuple>::norm(s,dp);
      using std::fabs;
      double r = fabs(get(distance_metric,get<Idx>(s))(get<Idx>(dp), get<Idx>(s)));
      if(result < r)
        return r;
      else
        return result;
    };
  };
  
  template <typename SpaceTuple> 
  struct inf_norm_tuple_distance_impl<0,SpaceTuple> {
    template <typename PointType>
    static double distance(const SpaceTuple& s, const PointType& p1, const PointType& p2) {
      using std::fabs;
      double r = fabs(get(distance_metric,get<0>(s))(get<0>(p1), get<0>(p2), get<0>(s)));
      return r;
    };
      
    template <typename PointDiff>
    static double norm(const SpaceTuple& s, const PointDiff& dp) {
      using std::fabs;
      double r = fabs(get(distance_metric,get<0>(s))(get<0>(dp), get<0>(s)));
      return r;
    };
  };
  


};


/**
 * This class is a Manhattan distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Manhattan norm to the point-difference vectors in the 
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct manhattan_tuple_distance : public serialization::serializable {
  
  manhattan_tuple_distance() { };
  
  template <typename Topology>
  manhattan_tuple_distance(const Topology&) { };
  
  /** 
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    return detail::manhattan_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::norm(s,a);
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
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    return detail::manhattan_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::distance(s,a,b);
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };

  RK_RTTI_MAKE_ABSTRACT_1BASE(manhattan_tuple_distance,0xC2410005,1,"manhattan_tuple_distance",serialization::serializable)
};



/**
 * This class is a Euclidean distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Euclidean norm to the point-difference vectors in the 
 * given topology (assuming it models the MetricSpaceConcept and that it is a tuple).
 */
struct euclidean_tuple_distance : public serialization::serializable {
  
  euclidean_tuple_distance() { };

  template <typename Topology>
  euclidean_tuple_distance(const Topology&) { };
  
  /** 
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology-tuple.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    using std::sqrt;
    return sqrt(detail::euclidean_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::norm(s,a));
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
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    using std::sqrt;
    return sqrt(detail::euclidean_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::distance(s,a,b));
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };

  RK_RTTI_MAKE_ABSTRACT_1BASE(euclidean_tuple_distance,0xC2410006,1,"euclidean_tuple_distance",serialization::serializable)
};





/**
 * This class is a Infinity-norm distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Infinity-norm to the point-difference vectors in the 
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct inf_norm_tuple_distance : public serialization::serializable {
  
  inf_norm_tuple_distance() { };
  
  template <typename Topology>
  inf_norm_tuple_distance(const Topology&) { };
  
  /** 
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \return The norm of the difference between two points on a topology-tuple.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    return detail::inf_norm_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::norm(s,a);
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
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    return detail::inf_norm_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::distance(s,a,b);
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };

  RK_RTTI_MAKE_ABSTRACT_1BASE(inf_norm_tuple_distance,0xC2410007,1,"inf_norm_tuple_distance",serialization::serializable)
};



/**
 * This class is a Euclidean distance metric functor which models the DistanceMetricConcept on a topology-tuple.
 * This class will simply apply the Euclidean norm to the point-difference vectors in the 
 * given topology (assuming it models the MetricSpaceConcept, and that it is a tuple).
 */
struct p_norm_tuple_distance : public serialization::serializable {
  
  int p_value;
  
  p_norm_tuple_distance(int aP = 2) : p_value(aP) { };
  
  template <typename Topology>
  p_norm_tuple_distance(const Topology&, int aP = 2) : p_value(aP) { };
  
  /** 
   * This function returns the norm of a difference between two points on a topology-tuple.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points  on a topology-tuple.
   */
  template <typename PointDiff, typename Topology>
  double operator()(const PointDiff& a, const Topology& s) const {
    using std::pow;
    return pow(detail::euclidean_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::norm(s,a,p_value), 1.0 / p_value);
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
  template <typename Point, typename Topology>
  double operator()(const Point& a, const Point& b, const Topology& s) const {
    using std::pow;
    return pow(detail::euclidean_tuple_distance_impl< arithmetic_tuple_size<Topology>::type::value - 1, Topology >::distance(s,a,b,p_value), 1.0 / p_value);
  };
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { 
    A & RK_SERIAL_SAVE_WITH_NAME(p_value);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
    A & RK_SERIAL_LOAD_WITH_NAME(p_value);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(p_norm_tuple_distance,0xC2410008,1,"p_norm_tuple_distance",serialization::serializable)
};



};

};


#endif


