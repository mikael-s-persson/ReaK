/**
 * \file sap_Ndof_metrics.hpp
 * 
 * This library provides an implementation of a distance metric within a temporal topology
 * which is based on the reach-time required by a sustained acceleration pulse motion (SAP) 
 * between two points.
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

#ifndef REAK_SAP_NDOF_METRICS_HPP
#define REAK_SAP_NDOF_METRICS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/optimization/optim_exceptions.hpp>

#include <ReaK/topologies/spaces/proper_metric_concept.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>

#include <boost/concept_check.hpp>

#include "generic_interpolator_factory.hpp"
#include "sustained_acceleration_pulse_Ndof.hpp"

namespace ReaK {

namespace pp {



/**
 * This functor class is a distance metric based on the reach-time of a SAP interpolation between
 * two points in a differentiable space.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology, bool MakeProper = false>
struct sap_Ndof_reach_time_metric : public serializable {
  
  typedef sap_Ndof_reach_time_metric<TimeSpaceType, MakeProper> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  sap_Ndof_reach_time_metric(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
                             t_space(aTimeSpace) { };
  
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
    try {
      detail::generic_interpolator_impl<sap_Ndof_interpolator,Topology,TimeSpaceType> interp;
      interp.initialize(a, b, 0.0, s, *t_space, *this);
      return interp.get_minimum_travel_time();
    } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
      return std::numeric_limits<double>::infinity();
    };
  };
  
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
    try {
      detail::generic_interpolator_impl<sap_Ndof_interpolator,Topology,TimeSpaceType> interp;
      interp.initialize(s.origin(), s.adjust(s.origin(),a), 0.0, s, *t_space, *this);
      return interp.get_minimum_travel_time();
    } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
      return std::numeric_limits<double>::infinity();
    };
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(t_space);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(t_space);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC241000D,1,"sap_Ndof_reach_time_metric",serializable)
};



/**
 * This functor class is a distance metric based on the reach-time of a SAP interpolation between
 * two points in a differentiable space.
 * \note This specialization is used for a proper distance metric.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType>
struct sap_Ndof_reach_time_metric<TimeSpaceType, true> : public serializable {
  
  typedef sap_Ndof_reach_time_metric<TimeSpaceType, true> self;
  
  shared_ptr<TimeSpaceType> t_space;
  
  sap_Ndof_reach_time_metric(const sap_Ndof_reach_time_metric<TimeSpaceType, false>& aRHS) : 
                             t_space(aRHS.t_space) { };
  
  sap_Ndof_reach_time_metric(const shared_ptr<TimeSpaceType>& aTimeSpace = shared_ptr<TimeSpaceType>(new TimeSpaceType())) : 
                             t_space(aTimeSpace) { };
  
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
    double d = std::numeric_limits<double>::infinity();
    try {
      // guarantee symmetry by computing reach-times in both directions:
      detail::generic_interpolator_impl<sap_Ndof_interpolator,Topology,TimeSpaceType> interp;
      interp.initialize(a, b, 0.0, s, *t_space, *this);
      d = interp.get_minimum_travel_time();
      interp.initialize(b, a, 0.0, s, *t_space, *this);
      double d2 = interp.get_minimum_travel_time();
      // pick the minimal value:
      if( d2 < d)
        d = d2;
    } catch(optim::infeasible_problem& e) { RK_UNUSED(e); };
    if( d == std::numeric_limits<double>::infinity() )
      return get(proper_metric, s)(a, b, s);
    return d;
  };
  
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
    typedef typename topology_traits<Topology>::point_type PointType;
    double d = std::numeric_limits<double>::infinity();
    try {
      PointType p = s.adjust(s.origin(),a);
      // guarantee symmetry by computing reach-times in both directions:
      detail::generic_interpolator_impl<sap_Ndof_interpolator,Topology,TimeSpaceType> interp;
      interp.initialize(s.origin(), p, 0.0, s, *t_space, *this);
      d = interp.get_minimum_travel_time();
      interp.initialize(p, s.origin(), 0.0, s, *t_space, *this);
      double d2 = interp.get_minimum_travel_time();
      // pick the minimal value:
      if( d2 < d)
        d = d2;
    } catch(optim::infeasible_problem& e) { RK_UNUSED(e); };
    if( d == std::numeric_limits<double>::infinity() )
      return get(proper_metric, s)(a, s);
    return d;
  };
  
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(t_space);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(t_space);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC241000D,1,"sap_Ndof_reach_time_metric",serializable)
};



template <typename TimeSpaceType>
struct is_metric_symmetric< sap_Ndof_reach_time_metric<TimeSpaceType, false> > : boost::mpl::false_ { };

template <typename TimeSpaceType>
struct is_metric_symmetric< sap_Ndof_reach_time_metric<TimeSpaceType, true> > : boost::mpl::true_ { };


template <typename TimeSpaceType, bool MakeProper>
struct get_proper_metric_from_metric< sap_Ndof_reach_time_metric<TimeSpaceType, MakeProper> > {
  typedef sap_Ndof_reach_time_metric<TimeSpaceType, true> type;
};


};

};



#endif









