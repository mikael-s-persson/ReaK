/**
 * \file sequential_trajectory_concept.hpp
 * 
 * This library defines the traits and concepts related to a sequential spatial trajectory. A 
 * trajectory is simply a continuous curve in a temporal topology (or time-space) which can be travelled 
 * sequentially via either increments in time or in fractions between waypoints.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_SEQUENTIAL_TRAJECTORY_CONCEPT_HPP
#define REAK_SEQUENTIAL_TRAJECTORY_CONCEPT_HPP


#include <boost/concept_check.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/temporal_space_concept.hpp>

namespace ReaK {

namespace pp {

  
/**
 * This traits class defines the traits that characterize a sequential spatial trajectory within a 
 * temporal topology.
 * \tparam SequentialTraj The spatial trajectory type for which the traits are sought.
 */
template <typename SequentialTraj>
struct sequential_trajectory_traits {
  /** This type describes a point in the space or topology. */
  typedef typename SequentialTraj::point_type point_type;
  
  /** This type describes an iterator, corresponding to a point on the trajectory, which can be incremented by time to travel to the next iterator. */
  typedef typename SequentialTraj::point_time_iterator point_time_iterator;
  /** This type describes an iterator, corresponding to a point on the trajectory, which can be incremented by a fraction between waypoints to travel to the next iterator. */
  typedef typename SequentialTraj::point_fraction_iterator point_fraction_iterator;
  
  /** This type is the topology type in which the path exists. */
  typedef typename SequentialTraj::topology topology;
  
};


/**
 * This concept class defines the requirements for a type to model a sequential trajectory 
 * as used in ReaK::pp. A sequential trajectory is a continuous curve within a temporal topology 
 * which can be travelled sequentially via either increments in time or in fractions 
 * between waypoints.
 * 
 * Required concepts:
 * 
 * The topology should model the TemporalSpaceConcept.
 * 
 * Valid expressions:
 * 
 * tit = traj.begin_time_travel();  The start of the time-iterator range of the sequential trajectory can be obtained.
 * 
 * tit = traj.end_time_travel();  The end of the time-iterator range (one-past-last) of the sequential trajectory can be obtained.
 * 
 * pt = *tit;  A point can be obtained from dereferencing a time-iterator.
 * 
 * tit = tit + d; 
 * tit = d + tit; 
 * tit += d;
 * tit = tit - d;
 * tit -= d;  A time-iterator can be incremented by a time (double).
 * 
 * b = (tit != tit);
 * b = (tit == tit);  Two time-iterator can be compared for inequality.
 * 
 * fit = traj.begin_fraction_travel();  The start of the fraction-iterator range of the sequential trajectory can be obtained.
 * 
 * fit = traj.end_fraction_travel();  The end of the fraction-iterator range (one-past-last) of the sequential trajectory can be obtained.
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
 * d = traj.travel_distance(pt,pt);  The travel distance (as of the distance-metric), along the trajectory (p), between two points (pt,pt), can be obtained.
 * 
 * \tparam SequentialTraj The type to be checked for the requirements of this concept.
 * \tparam Topology The topology in which the trajectory should reside.
 */
template <typename SequentialTraj, typename Topology>
struct SequentialTrajectoryConcept {
  
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  
  SequentialTraj* traj;
  typename topology_traits<Topology>::point_type pt;
  double d;
  bool b;
  typename sequential_trajectory_traits<SequentialTraj>::point_time_iterator tit;
  typename sequential_trajectory_traits<SequentialTraj>::point_fraction_iterator fit;
  
  BOOST_CONCEPT_USAGE(SequentialTrajectoryConcept)
  {
    tit = traj->begin_time_travel();
    tit = traj->end_time_travel();
    
    pt = *tit;
    
    tit = tit + d; 
    tit = d + tit; 
    tit += d;
    tit = tit - d;
    tit -= d;
    
    b = (tit != tit);
    b = (tit == tit);
    
    fit = traj->begin_fraction_travel();
    fit = traj->end_fraction_travel();
    
    pt = *fit;
    
    fit = fit + d; 
    fit = d + fit; 
    fit += d;
    fit = fit - d;
    fit -= d;
    
    b = (fit != fit);
    b = (fit == fit);
    
    d = traj->travel_distance(pt,pt);
  };
  
};



};

};


#endif









