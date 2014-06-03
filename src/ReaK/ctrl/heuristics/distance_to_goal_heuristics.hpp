/**
 * \file distance_to_goal_heuristics.hpp
 * 
 * This library defines the distance to goal heuristic that work on points of a topology and a given goal 
 * trajectory.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_DISTANCE_TO_GOAL_HEURISTICS_HPP
#define REAK_DISTANCE_TO_GOAL_HEURISTICS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/serializable.hpp>

#include <ReaK/ctrl/interpolation/spatial_trajectory_concept.hpp>


namespace ReaK {

namespace pp {
  
/**
 * 
 */
template <typename Trajectory, typename DistanceMetric = default_distance_metric>
class distance_to_goal_heuristic : public serialization::serializable {
  public:
    typedef distance_to_goal_heuristic<Trajectory,DistanceMetric> self;
    typedef typename spatial_trajectory_traits<Trajectory>::point_type point_type;
    typedef typename spatial_trajectory_traits<Trajectory>::topology topology;
    typedef typename spatial_trajectory_traits<Trajectory>::space_topology space_topology;
    
    BOOST_CONCEPT_ASSERT((SpatialTrajectoryConcept<Trajectory>));
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,topology>));
    
  private:
    shared_ptr<const Trajectory> traject;
    DistanceMetric dist;
  public:
    
    distance_to_goal_heuristic(const shared_ptr<const Trajectory>& aTraject, 
                               const DistanceMetric& aDist = DistanceMetric()) :
                               traject(aTraject), dist(aDist) { };
  
    /** 
     * This function returns the heuristic value of a point on a given topology.
     * \param a The point to evaluate.
     * \param s The topology or space on which the points lie.
     * \return The distance between two points on a topology.
     */
    double evaluate_point(const point_type& a, const topology& s) const {
      return dist(a, traject->get_point_at_time(a.time), s);
    };
    /** 
     * This function returns the heuristic value of an edge between two points on a topology.
     * \param a The start point of the edge.
     * \param b The end point of the edge.
     * \param s The topology or space on which the points lie.
     * \return The evalutation of the heuristic on the given edge.
     */
    double evaluate_edge(const point_type& a, const point_type& b, const topology& s) const {
      return dist(a, b, s);
    };
    
      
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(traject)
        & RK_SERIAL_SAVE_WITH_NAME(dist);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(traject)
        & RK_SERIAL_LOAD_WITH_NAME(dist);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2410021,1,"distance_to_goal_heuristic",serialization::serializable)
};


};

};


#endif


