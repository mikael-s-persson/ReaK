/**
 * \file belief_state_predictor.hpp
 * 
 * This library provides a class template which can generate a predicted trajectory of 
 * belief-states. This class template relies on several classes to implement its functionality.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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



/**
 * This class template can generate a predicted trajectory of 
 * belief-states. This class template relies on several classes to implement its functionality. 
 * Given a temporal space (or topology) whose point-type is a belief-state and a belief-state
 * predictor function type, this class predict the beliefs at equal time intervals from a 
 * starting belief-state, and can provide belief-state predictions at any given time.
 * 
 * Models: SpatialTrajectoryConcept, PredictedTrajectoryConcept, and provides an interface similar to STL containers (a sort-of hybrid of std::map and std::list).
 * 
 * \tparam BeliefTopology The topology of the belief-space, should model the TemporalSpaceConcept.
 * \tparam BeliefPredictor The belief-state predictor function type, should model BeliefPredictorConcept.
 * \tparam InputTrajectory The input vector trajectory to provide input vectors at any given time, should model the SpatialTrajectoryConcept.
 */
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
   
    /**
     * This nested class defines a waypoint in a predicted belief trajectory. A waypoint 
     * associates a belief point and a predictor functor associated to that belief point.
     * This class acts as the value-type of the container-like interface.
     */
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
    
    /**
     * Parametrized constructor.
     * \param aInitialPoint The starting belief-point of the predicted trajectory.
     * \param aInitialPredictor The predictor functor associated to the starting belief-point.
     * \param aSpace The belief-state temporal topology object.
     * \param aInputTrajectory The input-vector trajectory used to compute the input vectors necessary for the belief prediction.
     */
    belief_predicted_trajectory(const point_type& aInitialPoint,
                                const predictor_type& aInitialPredictor,
				const topology& aSpace = topology(),
				const InputTrajectory& aInputTrajectory = InputTrajectory()) :
				waypoints(1,waypoint(aInitialPoint,aInitialPredictor)),
				space(aSpace),
				input(aInputTrajectory) { };
    
    /**
     * Returns the first waypoint in the underlying waypoint container.
     * \return The first waypoint in the underlying waypoint container.
     */
    const waypoint& front() const { return waypoints.front(); };
    /**
     * Returns the last waypoint in the underlying waypoint container.
     * \return The last waypoint in the underlying waypoint container.
     */
    const waypoint& back() const { return waypoints.back(); };
    
    /**
     * Returns the number of waypoints in the underlying container.
     * \return The number of waypoints in the underlying container.
     */
    size_type size() const { return waypoints.size(); };
    /**
     * Determines if the waypoint container is empty.
     * \return Always false, the underlying container is not allowed to be emptied.
     */
    bool empty() const { return false; };
    /**
     * Returns the maximum number of waypoints the underlying container will allow.
     * \return The maximum number of waypoints the underlying container will allow.
     */
    size_type max_size() const { return waypoints.max_size(); };
    /**
     * Resizes the underlying container to contain sz number of waypoints. Any new 
     * waypoint will first be a copy of value but will then be predicted from the 
     * existing waypoints. If sz is zero, the size after the resize will be 1.
     * \param sz The new size for the waypoint container.
     * \param value The value of the waypoints to be appended to the container if the new size is larger than the current size.
     */
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
    /**
     * Resizes the underlying container to contain sz number of waypoints. Any new 
     * waypoint will first be a copy of the last waypoint but will then be predicted from the 
     * existing waypoints. If sz is zero, the size after the resize will be 1.
     * \param sz The new size for the waypoint container.
     */
    void resize( size_type sz ) { resize(sz, waypoints.back()); };
    
    /**
     * Returns a const-iterator to the first waypoint in the underlying container.
     * \return A const-iterator to the first waypoint in the underlying container.
     */
    const_waypoint_descriptor begin() const { return waypoints.begin(); };
    /**
     * Returns a const-iterator to the one-past-last waypoint in the underlying container.
     * \return A const-iterator to the one-past-last waypoint in the underlying container.
     */
    const_waypoint_descriptor end() const { return waypoints.end(); };
    /**
     * Returns a const-reverse-iterator to the first waypoint in the underlying container.
     * \return A const-reverse-iterator to the first waypoint in the underlying container.
     */
    typename container_type::const_reverse_iterator rbegin() const { return waypoints.rbegin(); };
    /**
     * Returns a const-reverse-iterator to the one-past-last waypoint in the underlying container.
     * \return A const-reverse-iterator to the one-past-last waypoint in the underlying container.
     */
    typename container_type::const_reverse_iterator rend() const { return waypoints.rend(); };
    
    /**
     * Returns the allocator object of the underlying container.
     * \return The allocator object of the underlying container.
     */
    typename container_type::allocator_type get_allocator() const { return waypoints.get_allocator(); };
    
    /**
     * Gets the waypoint that is associated to a time just before the given time.
     * \param t The time at which the waypoint is sought.
     * \return The waypoint that is associated to a time just before the given time.
     */
    const waypoint& operator[](time_type t) const { return *(get_waypoint(t).first); };
    
    /**
     * Pops the last element off the underlying container, unless there is just one waypoint remaining.
     */
    void pop_back() { if(waypoints.size() > 1) waypoints.pop_back(); };
    /**
     * Pops the first element off the underlying container, unless there is just one waypoint remaining.
     */
    void pop_front() { if(waypoints.size() > 1) waypoints.pop_front(); };
    /**
     * Adds a waypoint at the end of the underlying container.
     * \param value The new element to be added to the end of the container.
     */
    void push_back(const waypoint& value) {
      waypoint_descriptor it = (++(waypoints.rbegin())).base();
      waypoints.push_back(value);
      belief_state b = it->point.pt;
      time_type t = it->point.time;
      b = it->predictor.predict_belief(b, t, input.get_point(t));
      waypoints.back().point = point_type(t + it->predictor.get_time_step(),b);
    };
    /**
     * Adds a waypoint at the begining of the underlying container.
     * \param value The new element to be added to the begining of the container.
     */
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
    
    /**
     * Moves the belief-point of a waypoint forward by a time-difference (or backward if negative).
     * \param wp The waypoint from which to move.
     * \param dt The time difference to move from the given waypoint.
     * \return The adjusted waypoint.
     */
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
    /**
     * Returns the travel distance, along the trajectory, from one waypoint to another.
     * \param wp1 The first waypoint.
     * \param wp2 The second waypoint.
     * \return The travel distance, along the trajectory, from wp1 to wp2.
     */
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
    /**
     * Gets the waypoint that is associated to a time just before the given time.
     * \param t The time at which the waypoint is sought.
     * \return The waypoint that is associated to a time just before the given time.
     */
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
    
    /**
     * Resets the initial point of the predicted trajector with a waypoint and an associated time. 
     * This function triggers the elimination of all waypoints prior to this time and triggers 
     * the recomputation of the belief predictions of all waypoints after this time.
     * \param wp The waypoint that will become the first, initial waypoint.
     * \param t The time of the initial waypoint.
     */
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
    
    /**
     * Moves the belief-point forward by a time-difference (or backward if negative). 
     * \note It is more efficient to use the overload of this function that uses a waypoint, if the waypoint is already known.
     * \param wp The waypoint from which to move.
     * \param dt The time difference to move from the given waypoint.
     * \return The adjusted waypoint.
     */
    point_type move_away_from(const point_type& pt, const time_difference_type& dt) const {
      waypoint_point_pair wp = get_waypoint(pt.time);
      wp.second = pt;
      return move_away_from(wp,dt).second;
    };
    /**
     * Returns the travel distance, along the trajectory, from one belief-point to another.
     * \note It is more efficient to use the overload of this function that uses a waypoint, if the waypoint is already known.
     * \param p1 The first belief-point.
     * \param p2 The second belief-point.
     * \return The travel distance, along the trajectory, from p1 to p2.
     */
    double travel_distance(const point_type& p1, const point_type& p2) const {
      waypoint_point_pair wp1 = get_waypoint(p1.time); wp1.second = p1;
      waypoint_point_pair wp2 = get_waypoint(p2.time); wp2.second = p2;
      return travel_distance(wp1,wp2);
    };
    /**
     * Gets the point that is associated to a time just before the given time.
     * \param t The time at which the point is sought.
     * \return The point that is associated to a time just before the given time.
     */
    point_type get_point(const time_type& t) const {
      return get_waypoint(t).second;
    };

    /**
     * Resets the initial point of the predicted trajector with a point and an associated time. 
     * This function triggers the elimination of all waypoints prior to this time and triggers 
     * the recomputation of the belief predictions of all waypoints after this time.
     * \note It is more efficient to use the overload of this function that uses a waypoint, if the waypoint is already known.
     * \param pt The point that will become the first, initial waypoint.
     * \param t The time of the initial point.
     */
    void set_initial_point(const point_type& pt, time_type t) {
      waypoint_point_pair wp = get_waypoint(t);
      wp.second = pt;
      set_initial_point(wp,t);
    };
    
};


};


};


#endif









