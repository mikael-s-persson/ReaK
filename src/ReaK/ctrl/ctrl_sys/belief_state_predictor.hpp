
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

#ifndef BELIEF_STATE_PREDICTOR_HPP
#define BELIEF_STATE_PREDICTOR_HPP

#include "belief_state_concept.hpp"
#include "state_vector_concept.hpp"
#include "path_planning/predicted_trajectory_concept.hpp"
#include "path_planning/temporal_space.hpp"

#include <deque>

namespace ReaK {

namespace ctrl {


template <typename BeliefTopology, typename BeliefPredictor>
class belief_predicted_trajectory {
  public:
    typedef belief_predicted_trajectory<BeliefTopology,BeliefPredictor> self;
    typedef BeliefPredictor predictor_type;
    
    typedef pp::temporal_space< BeliefTopology > topology;
    typedef typename pp::temporal_topology_traits<topology>::time_topology time_topology;
    typedef typename pp::temporal_topology_traits<topology>::space_topology space_topology;
    typedef typename topology::distance_metric;
    
    typedef typename pp::temporal_topology_traits< topology >::point_type point_type;
    typedef typename pp::temporal_topology_traits< topology >::point_difference_type point_difference_type;
    
    typedef typename time_topology::point_type time_type;
    typedef typename time_topology::point_difference_type time_difference_type;
    
    struct waypoint {
      point_type point;
      predictor_type predictor;
      waypoint(const point_type& aPoint = point_type(),
	       const point_difference_type& aPredictor = point_difference_type()) :
	       point(aPoint),
	       predictor(aPredictor) { };
    };
    
    typedef std::deque<waypoint> container_type;
    typedef typename container_type::iterator waypoint_descriptor;
    typedef typename container_type::const_iterator const_waypoint_descriptor;
    typedef std::pair< const_waypoint_descriptor, point_type > waypoint_point_pair;
    
  private:
    container_type waypoints;
    
  public:
    
    // pt = p.move_away_from(pt, dt);
    point_type move_away_from(const point_type& pt, const time_difference_type& dt) const {
      const_waypoint_descriptor it = waypoints.begin();
      const_waypoint_descriptor pit = it++;
      for(; it != waypoints.end(); ++it) {
	if(pt.time < it->point.time)
	  break;
      };
      //...TODO
    };
    // d = p.travel_distance(pt, pt);
    double travel_distance(const point_type& p1, const point_type& p2) const {
      
    };
    // pt = p.get_point(t);
    point_type get_point(const time_type& t) const {
      const_waypoint_descriptor it = waypoints.begin();
      for(; it != waypoints.end(); ++it) {
	if(t < it->point.time)
	  break;
      };
      if(it == waypoints.end())
	it = (++(waypoints.rbegin()).base();
      else
	--it;
      return it->point;
    };
    
    // w_p = p.move_away_from(w_p, dt);
    waypoint_point_pair move_away_from(const waypoint_point_pair& wp, const time_difference_type& dt) const {
      
    };
    // d = p.travel_distance(w_p, w_p);
    double travel_distance(const waypoint_point_pair& wp1, const waypoint_point_pair& wp2) const {
      
    };
    // w_p = p.get_point(t);
    waypoint_point_pair get_waypoint(const time_type& t) const {
      const_waypoint_descriptor it = waypoints.begin();
      for(; it != waypoints.end(); ++it) {
	if(t < it->point.time)
	  break;
      };
      if(it == waypoints.end())
	it = (++(waypoints.rbegin()).base();
      else
	--it;
      return waypoint_point_pair(it, it->point);
    };
    
    // p.set_initial_point(pt, t);
    void set_initial_point(const waypoint_point_pair& wp, const time_type& t) {
      
    };
    // p.set_initial_point(w_p, t);
    void set_initial_point(const waypoint_point_pair& wp, const time_type& t) {
      
    };

};


};


};


#endif









