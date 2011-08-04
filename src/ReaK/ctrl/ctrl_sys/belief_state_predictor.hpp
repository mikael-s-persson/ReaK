
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
#include <iterator>

namespace ReaK {

namespace ctrl {


template <typename BeliefTopology, 
          typename BeliefPredictor,
	  typename InputTrajectory>
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
    
    typedef typename pp::metric_topology_traits<space_topology>::point_type belief_state;
    typedef typename pp::metric_topology_traits<space_topology>::point_difference_type belief_state_diff;
    
    struct waypoint {
      point_type point;
      predictor_type predictor;
      waypoint(const point_type& aPoint,
	       const predictor_type& aPredictor) :
	       point(aPoint),
	       predictor(aPredictor) { };
    };
    
    typedef std::deque<waypoint> container_type;
    typedef typename container_type::iterator waypoint_descriptor;
    typedef typename container_type::const_iterator const_waypoint_descriptor;
    typedef std::pair< const_waypoint_descriptor, point_type > waypoint_point_pair;
    typedef typename container_type::size_type size_type;
    
  private:
    container_type waypoints;
    topology space;
    InputTrajectory input;
    
  public:
    
    belief_predicted_trajectory(const point_type& aInitialPoint,
                                const predictor_type& aInitialPredictor,
				const topology& aSpace = topology(),
				const InputTrajectory& aInputTrajectory = InputTrajectory()) :
				waypoints(1,waypoint(aInitialPoint,aInitialPredictor)),
				space(aSpace),
				input(aInputTrajectory) { };
    
    const waypoint& front() const { return waypoints.front(); };
    const waypoint& back() const { return waypoints.back(); };
    
    size_type size() const { return waypoints.size(); };
    bool empty() const { return false; };
    size_type max_size() const { return waypoints.max_size(); };
    void resize( size_type sz, const waypoint& value ) {
      if(sz > waypoints.size()) {
	waypoint_descriptor it = (++(waypoints.rbegin())).base();
	waypoints.resize(sz,value);
	belief_state b = it->point.pt;
	time_type t = it->point.time;
        for(;it != waypoints.end();++it) {
	  it->point = point_type(t, b);
	  it->predictor.predict_belief(b, t, input.get_point(t));
	  t = t + it->predictor.get_time_step();
        };
      } else if(sz != waypoints.size())
	waypoints.resize((sz ? sz : 1),value);
    };
    void resize( size_type sz ) { resize(sz, waypoints.back()); };
    
    const_waypoint_descriptor begin() const { return waypoints.begin(); };
    const_waypoint_descriptor end() const { return waypoints.end(); };
    typename container_type::const_reverse_iterator rbegin() const { return waypoints.rbegin(); };
    typename container_type::const_reverse_iterator rend() const { return waypoints.rend(); };
    
    typename container_type::allocator_type get_allocator() const { return waypoints.get_allocator(); };
    
    const waypoint& operator[](time_type t) const { return *(get_waypoint(t).first); };
    
    void pop_back() { if(waypoints.size() > 1) waypoints.pop_back(); };
    void pop_front() { if(waypoints.size() > 1) waypoints.pop_front(); };
    void push_back(const waypoint& value) {
      waypoint_descriptor it = (++(waypoints.rbegin())).base();
      waypoints.push_back(value);
      belief_state b = it->point.pt;
      time_type t = it->point.time;
      b = it->predictor.predict_belief(b, t, input.get_point(t));
      waypoints.back().point = point_type(t + it->predictor.get_time_step(),b);
    };
    void push_front(const waypoint& value) {
      waypoints.push_front(value);
      waypoint_descriptor it = waypoints.begin();
      belief_state b = it->point.pt;
      time_type t = it->point.time;
      for(;it != waypoints.end();++it) {
	it->point = point_type(t, b);
	b = it->predictor.predict_belief(b, t, input.get_point(t));
	t = t + it->predictor.get_time_step();
      };
    };
    
    // w_p = p.move_away_from(w_p, dt);
    waypoint_point_pair move_away_from(const waypoint_point_pair& wp, const time_difference_type& dt) const {
      const_waypoint_descriptor cit = wp.first;
      time_difference_type dt_tot = cit->point.time - wp.second.time;
      while(dt_tot < dt) {
        ++cit;
        if(cit == waypoints.end())
	  return waypoint_point_pair((++(waypoints.rbegin())).base(),point_type(wp.second.time + dt, wp.second.pt));
        dt_tot = cit->point.time - wp.second.time;
      };
      const point_type& tmp = cit->point;
      return waypoint_point_pair(--cit,space.move_position_toward(wp.second,dt / dt_tot,tmp));
    };
    // d = p.travel_distance(w_p, w_p);
    double travel_distance(const waypoint_point_pair& wp1, const waypoint_point_pair& wp2) const {
      if(wp1.first == wp2.first)
	return space.distance(wp1.second,wp2.second);
      ++(wp1.first);
      if(wp1.first == waypoints.end())
	return space.distance(wp1.second,wp2.second);
      double sum = space.distance(wp1.second,wp1.first->point);
      while( wp1.first != wp2.first ) {
	const point_type& tmp = wp1.first->point;
	++(wp1.first);
	sum += space.distance(tmp, wp1.first->point);
      };
      sum += space.distance(wp2.first->point, wp2.second);
      return sum;
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
    
    // p.set_initial_point(w_p, t);
    void set_initial_point(const waypoint_point_pair& wp, time_type t) {
      const_waypoint_descriptor cit = wp.first;
      waypoint_descriptor it = waypoints.begin();
      std::advance(it, std::distance<const_waypoint_descriptor>(it,cit));
      it->point = wp.second;
      waypoints.erase(waypoints.begin(),it);
      belief_state b = it->point.pt;
      for(;it != waypoints.end();++it) {
	it->point = point_type(t, b);
	b = it->predictor.predict_belief(b, t, input.get_point(t));
	t = t + it->predictor.get_time_step();
      };
    };
    
    // pt = p.move_away_from(pt, dt);
    point_type move_away_from(const point_type& pt, const time_difference_type& dt) const {
      waypoint_point_pair wp = get_waypoint(pt.time);
      wp.second = pt;
      return move_away_from(wp,dt).second;
    };
    // d = p.travel_distance(pt, pt);
    double travel_distance(const point_type& p1, const point_type& p2) const {
      waypoint_point_pair wp1 = get_waypoint(p1.time); wp1.second = p1;
      waypoint_point_pair wp2 = get_waypoint(p2.time); wp2.second = p2;
      return travel_distance(wp1,wp2);
    };
    // pt = p.get_point(t);
    point_type get_point(const time_type& t) const {
      return get_waypoint(t).second;
    };

    // p.set_initial_point(pt, t);
    void set_initial_point(const point_type& pt, time_type t) {
      waypoint_point_pair wp = get_waypoint(t);
      wp.second = pt;
      set_initial_point(wp,t);
    };
    
};


};


};


#endif









