/**
 * \file svp_Ndof_metrics.hpp
 * 
 * This library provides an implementation of a distance metric within a temporal topology
 * which is based on the reach-time required by a sustained velocity pulse motion (SVP) 
 * between two points in an N-dof differentiable space.
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

#ifndef REAK_SVP_NDOF_METRICS_HPP
#define REAK_SVP_NDOF_METRICS_HPP

#include "base/defs.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "generic_interpolator_factory.hpp"
#include "sustained_velocity_pulse_Ndof.hpp"
#include "topologies/time_topology.hpp"
#include "optimization/optim_exceptions.hpp"

namespace ReaK {

namespace pp {


/**
 * This functor class is a distance metric based on the reach-time of a SVP interpolation between
 * two points in a differentiable space.
 * \tparam TimeSpaceType The time topology type against which the interpolation is done.
 */
template <typename TimeSpaceType = time_topology>
struct svp_Ndof_reach_time_metric : public serialization::serializable {
  
  typedef svp_Ndof_reach_time_metric<TimeSpaceType> self;
  
  shared_ptr<const TimeSpaceType> t_space;
  
  svp_Ndof_reach_time_metric(const shared_ptr<const TimeSpaceType>& aTimeSpace = shared_ptr<const TimeSpaceType>(new TimeSpaceType())) : 
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
      detail::generic_interpolator_impl<svp_Ndof_interpolator,Topology,TimeSpaceType> interp;
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
      detail::generic_interpolator_impl<svp_Ndof_interpolator,Topology,TimeSpaceType> interp;
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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC241000C,1,"svp_Ndof_reach_time_metric",serialization::serializable)
};


};

};

#endif









