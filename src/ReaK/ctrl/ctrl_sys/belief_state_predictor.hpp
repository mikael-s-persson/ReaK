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

#ifndef REAK_BELIEF_STATE_PREDICTOR_HPP
#define REAK_BELIEF_STATE_PREDICTOR_HPP

#include "belief_state_concept.hpp"
#include "state_vector_concept.hpp"
#include "discrete_sss_concept.hpp"

#include "path_planning/predicted_trajectory_concept.hpp"
#include "path_planning/temporal_space_concept.hpp"

#include "topologies/temporal_space.hpp"
#include "topologies/vector_topology.hpp"

#include <deque>
#include <iterator>

namespace ReaK {

namespace ctrl {



/**
 * This class template can generate a predicted trajectory of 
 * belief-states. This class template relies on several classes to implement its functionality. 
 * Given a space (or topology) whose point-type is a belief-state, and a belief-state
 * predictor function type, this class predict the beliefs at equal time intervals from a 
 * starting belief-state, and can provide belief-state predictions at any given time.
 * The predictions are populated in a Just-In-Time (JIT) fashion (i.e., lazily), and 
 * internally maintains a current horizon (time to last pre-computed prediction) which can
 * be explicitely extended or pruned (eliminating predictions beyond a certain time).
 * 
 * Models: SpatialTrajectoryConcept, PredictedTrajectoryConcept.
 * 
 * \tparam BeliefTopology The topology of the belief-space, should model the BeliefSpaceConcept.
 * \tparam BeliefPredictor The belief-state predictor function type, should model BeliefPredictorConcept.
 * \tparam InputTrajectory The input vector trajectory to provide input vectors at any given time, should model the SpatialTrajectoryConcept over a vector-topology of input vectors.
 */
template <typename BeliefTopology, 
          typename BeliefPredictor,
	  typename InputTrajectory>
class belief_predicted_trajectory {
  public:
    
    typedef belief_predicted_trajectory<BeliefTopology,BeliefPredictor,InputTrajectory> self;
    typedef BeliefPredictor predictor_type;
    
    typedef pp::temporal_space< BeliefTopology > topology;
    typedef shared_ptr<topology> topology_ptr;
    typedef typename pp::temporal_space_traits<topology>::time_topology time_topology;
    typedef typename pp::temporal_space_traits<topology>::space_topology space_topology;
    
    typedef typename pp::temporal_space_traits< topology >::point_type point_type;
    typedef typename pp::temporal_space_traits< topology >::point_difference_type point_difference_type;
    
    typedef typename pp::topology_traits<time_topology>::point_type time_type;
    typedef typename pp::topology_traits<time_topology>::point_difference_type time_difference_type;
    
    typedef typename pp::topology_traits<space_topology>::point_type belief_state;
    typedef typename pp::topology_traits<space_topology>::point_difference_type belief_state_diff;
    
    typedef typename discrete_sss_traits< predictor_type >::input_type input_type;
    
    BOOST_CONCEPT_ASSERT((pp::TemporalSpaceConcept<topology>));
    BOOST_CONCEPT_ASSERT((BeliefSpaceConcept<space_topology>));
    BOOST_CONCEPT_ASSERT((BeliefPredictorConcept<BeliefPredictor, space_topology>));
    BOOST_CONCEPT_ASSERT((SpatialTrajectoryConcept<InputTrajectory, pp::vector_topology<input_type> >));
   
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
    
    typedef std::map<time_type, waypoint> container_type;
    typedef typename container_type::value_type container_value;
    typedef typename container_type::iterator waypoint_descriptor;
    typedef typename container_type::const_iterator const_waypoint_descriptor;
    typedef std::pair< const_waypoint_descriptor, point_type > waypoint_point_pair;
    typedef typename container_type::size_type size_type;
    
  private:
    container_type waypoints;
    topology_ptr space;
    InputTrajectory input;
    
    waypoint_descriptor updated_end;
    
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
                                const topology_ptr& aSpace = topology_ptr(new topology()),
                                const InputTrajectory& aInputTrajectory = InputTrajectory()) :
                                waypoints(1,container_value(aInitialPoint.time, waypoint(aInitialPoint,aInitialPredictor))),
                                space(aSpace),
                                input(aInputTrajectory),
                                updated_end(waypoints.end()) { };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      
      // this little tricky piece of code is due to iterators being invalidated across the container-swap.
      time_type lhs_horiz = lhs.get_current_horizon();
      time_type rhs_horiz = rhs.get_current_horizon();
      lhs.waypoints.swap(rhs.waypoints);
      lhs.updated_end = lhs.waypoints.lower_bound(lhs_horiz);
      rhs.updated_end = rhs.waypoints.lower_bound(rhs_horiz);
      
      swap(lhs.space,rhs.space);
      swap(lhs.input,rhs.input);
    };
    
    
    /**
     * Gets the predicted trajectory horizon, i.e., the time up to which predictions are currently valid.
     * \return The horizon for the predictions currently up-to-date in this trajectory.
     */
    time_type get_current_horizon() const {
      waypoint_descriptor last_valid = updated_end; --last_valid;
      return last_valid->first;
    };
    
    /**
     * Sets the minimal predicted trajectory horizon, i.e., extending predictions to reach that horizon 
     * if they are not already reaching it.
     * \param t The new minimal horizon for the predictions up-to-date in this trajectory.
     */
    void set_minimal_horizon(const time_type& t) {
      
      waypoint_descriptor it = waypoints.lower_bound(t);
      
      waypoint_descriptor last_valid = updated_end; --last_valid;
      while( ( updated_end != waypoints.end() ) && ( last_valid->first < t ) ) {
        updated_end->second.point.pt = last_valid->second.predictor.predict_belief(
          space->get_space_topology(), 
          last_valid->second.point.pt, last_valid->second.point.time, 
          input.get_point_at_time(last_valid->second.point.time));
        ++last_valid; ++updated_end;
      };
      
      if(it == waypoints.end()) {
        // must perform predictions from end to past time t.
        while( last_valid->first < t ) {
          updated_end = waypoints.insert(waypoints.end(), container_value(
            last_valid->first + it->predictor.get_time_step(), 
            last_valid->second));
          updated_end->second.point.pt = last_valid->second.predictor.predict_belief(
            space->get_space_topology(), 
            last_valid->second.point.pt, last_valid->first, 
            input.get_point_at_time(last_valid->first));
          updated_end->second.point.time = updated_end->first;
          last_valid = updated_end;
        };
        updated_end = waypoints.end();
      };
      
    };
    
    /**
     * Sets the predicted trajectory horizon, i.e., pruning predictions that are 
     * "too far ahead", or extending predictions to reach that horizon.
     * \note This function will always preserve the starting point (stem), which can be reset with set_start_point.
     * \param t The new horizon for the predictions kept .
     */
    void set_exact_horizon(const time_type& t) {
      if( t > get_current_horizon() ) {
        set_minimal_horizon(t);
        return;
      };
      
      waypoint_descriptor it = waypoints.lower_bound(t);
      ++it; // <- this means that the pruning "conservative", and also, will always keep the beginning iterator.
      waypoints.erase(it, waypoints.end());
      updated_end = waypoints.end();
      
    };
    
    
    /**
     * Returns the travel distance, along the trajectory, from one waypoint to another.
     * \param wp1 The first waypoint.
     * \param wp2 The second waypoint.
     * \return The travel distance, along the trajectory, from wp1 to wp2.
     */
    double travel_distance(waypoint_point_pair wp1, const waypoint_point_pair& wp2) const {
      if(wp1.second.time > wp2.second.time)
        return travel_distance(wp2, wp1);
      
      if(wp1.first == wp2.first)
        return get(pp::distance_metric, *space)(wp1.second, wp2.second, *space);
      ++(wp1.first);
      if(wp1.first == waypoints.end())
        return get(pp::distance_metric, *space)(wp1.second, wp2.second, *space);
      double sum = get(pp::distance_metric, *space)(wp1.second, wp1.first->second.point, *space);
      while( wp1.first != wp2.first ) {
        const point_type& tmp = wp1.first->second.point;
        ++(wp1.first);
        sum += get(pp::distance_metric, *space)(tmp, wp1.first->second.point, *space);
      };
      sum += get(pp::distance_metric, *space)(wp2.first->second.point, wp2.second, *space);
      return sum;
    };
    
    /**
     * Gets the waypoint that is associated to a time just before the given time.
     * \note This function might add the necessary amount of waypoint predictions to provide a belief at the given time.
     * \param t The time at which the waypoint is sought.
     * \return The waypoint that is associated to a time just before the given time, and an interpolated approximation of the belief-state at the given time.
     */
    waypoint_point_pair get_waypoint_at_time(const time_type& t) {
      
      waypoint_descriptor it = waypoints.lower_bound(t);
      
      waypoint_descriptor last_valid = updated_end; --last_valid;
      while( ( updated_end != waypoints.end() ) && ( last_valid->first < t ) ) {
        updated_end->second.point.pt = last_valid->second.predictor.predict_belief(
          space->get_space_topology(), 
          last_valid->second.point.pt, last_valid->second.point.time, 
          input.get_point_at_time(last_valid->second.point.time));
        ++last_valid; ++updated_end;
      };
      
      if(it == waypoints.end()) {
        // must perform predictions from end to past time t.
        while( last_valid->first < t ) {
          updated_end = waypoints.insert(waypoints.end(), container_value(
            last_valid->first + it->predictor.get_time_step(), 
            last_valid->second));
          updated_end->second.point.pt = last_valid->second.predictor.predict_belief(
            space->get_space_topology(), 
            last_valid->second.point.pt, last_valid->first, 
            input.get_point_at_time(last_valid->first));
          updated_end->second.point.time = updated_end->first;
          last_valid = updated_end;
        };
        updated_end = waypoints.end();
        it = last_valid;
      };
      
      waypoint_point_pair result(it, it->second.point);
      --(result.first);
      result.second.time = t;
      
      time_difference_type dt_total = it->first - result.first->first;
      time_difference_type dt       = t - result.first->first;
      result.second.pt = space->get_space_topology().move_position_toward(result.first->second.point.pt, dt / dt_total, it->second.point.pt);
      
      return result;
    };
    
    /**
     * Gets the point that is associated to a time just before the given time.
     * \param t The time at which the point is sought.
     * \return The point that is associated to a time just before the given time.
     */
    point_type get_point_at_time(const time_type& t) {
      return get_waypoint_at_time(t).second;
    };
    
    /**
     * Returns the travel distance, along the trajectory, from one belief-point to another.
     * \note It is more efficient to use the overload of this function that uses a waypoint, if the waypoint is already known.
     * \param p1 The first belief-point.
     * \param p2 The second belief-point.
     * \return The travel distance, along the trajectory, from p1 to p2.
     */
    double travel_distance(const point_type& p1, const point_type& p2) const {
      waypoint_point_pair wp1 = get_waypoint_at_time(p1.time); wp1.second = p1;
      waypoint_point_pair wp2 = get_waypoint_at_time(p2.time); wp2.second = p2;
      return travel_distance(wp1,wp2);
    };
    
    /**
     * Moves the belief-point of a waypoint forward by a time-difference (or backward if negative).
     * \note This function might add the necessary amount of waypoint predictions to provide a belief at the given time.
     * \param wp The waypoint from which to move.
     * \param dt The time difference to move from the given waypoint.
     * \return The adjusted waypoint.
     */
    waypoint_point_pair move_time_diff_from(const waypoint_point_pair& wp, const time_difference_type& dt) {
      return get_waypoint_at_time(wp.second.time + dt);
    };
    
    /**
     * Moves the belief-point forward by a time-difference (or backward if negative). 
     * \note This function might add the necessary amount of waypoint predictions to provide a belief at the given time.
     * \param wp The waypoint from which to move.
     * \param dt The time difference to move from the given waypoint.
     * \return The adjusted waypoint.
     */
    point_type move_time_diff_from(const point_type& pt, const time_difference_type& dt) {
      return get_waypoint_at_time(pt.time + dt).second;
    };
    
    
    
    /**
     * Returns the starting point of the waypoints.
     * \return The starting point of the waypoints.
     */
    const point_type& get_start_point() const {
      return waypoints.begin()->second;
    };
    
    /**
     * Returns the starting waypoint-point-pair of the waypoints.
     * \return The starting waypoint-point-pair of the waypoints.
     */
    waypoint_pair get_start_waypoint() const {
      const_waypoint_descriptor start = waypoints.begin();
      return waypoint_pair(start, start->second);
    };
    
    /**
     * Returns the end point of the waypoints.
     * \return The end point of the waypoints.
     */
    const point_type& get_end_point() const {
      return waypoints.rbegin()->second;
    };
    
    /**
     * Returns the end waypoint-point-pair of the waypoints.
     * \return The end waypoint-point-pair of the waypoints.
     */
    waypoint_pair get_end_waypoint() const {
      const_waypoint_descriptor end = (++waypoints.rbegin()).base();
      return waypoint_pair(end, end->second);
    };
    
    
    /**
     * Resets the initial point of the predicted trajector with a point and an associated time. 
     * This function triggers the elimination of all waypoints prior to this time and triggers 
     * the recomputation of the belief predictions of all waypoints after this time.
     * \param pt The point that will become the first, initial waypoint.
     */
    void set_start_point(const point_type& pt, const predictor_type& pred) {
      using std::fabs;
      
      waypoint_descriptor it = waypoints.lower_bound(pt.time - 0.5 * dt);
      // check if the time difference is too much:
      if( fabs(pt.time - it->first) > 1e-4 * dt ) {
        // we have to completely reset the entire container:
        waypoints.clear();
        waypoints.insert(waypoints.end(), container_value(pt.time, waypoint(pt, pred)));
        updated_end = waypoints.end();
        return;
      };
      
      waypoints.erase(waypoints.begin(), it);
      
      it->second.point.time = it->first;
      it->second.point.pt   = pt.pt;
      it->second.predictor  = pred;
      updated_end = it; ++updated_end;
      
    };
    
    
    /**
     * Returns the starting time of the predictions.
     * \return The starting time of the predictions.
     */
    time_type get_start_time() const {
      return waypoints.begin()->first;
    };
    
    /**
     * Returns the end time of the predictions.
     * \return The end time of the predictions.
     */
    time_type get_end_time() const {
      return std::numeric_limits< time_type >::infinity();
    };
    
    
    
    
    
    struct point_time_iterator {
      const self* parent;
      time_type current_time;
      time_type end_time;
      
      point_time_iterator(const self* aParent, time_type aCurrentTime) :
                          parent(aParent), current_time(aCurrentTime), end_time(parent->get_current_horizon()) { };
      
      explicit point_time_iterator(const self* aParent) :
                                   parent(aParent), current_time(parent->get_current_horizon()), end_time(current_time) { };
      
      friend point_time_iterator& operator+=(point_time_iterator& lhs, double rhs) {
        lhs.current_time += rhs;
        if(current_time > end_time)
          current_time = end_time;
        return *this;
      };
      
      friend point_time_iterator operator+(point_time_iterator lhs, double rhs) { return (lhs += rhs); };
      friend point_time_iterator operator+(double lhs, point_time_iterator rhs) { return (rhs += lhs); };
      friend point_time_iterator operator-(point_time_iterator lhs, double rhs) { return (lhs += -rhs); };
      friend point_time_iterator& operator-=(point_time_iterator& lhs, double rhs) { return (lhs += -rhs); };
      
      friend bool operator==(const point_time_iterator& lhs, const point_time_iterator& rhs) {
        return ( ( lhs.parent == rhs.parent ) && ( lhs.current_time == rhs.current_time ) );
      };
      friend bool operator!=(const point_time_iterator& lhs, const point_time_iterator& rhs) { return !(lhs == rhs); };
      
      point_type operator*() const {
        return parent->get_point_at_time(current_time);
      };
      
    };
    
    /**
     * Returns the starting time-iterator along the trajectory.
     * \return The starting time-iterator along the trajectory.
     */
    point_time_iterator begin_time_travel() const {
      return point_time_iterator(this, this->waypoints.begin()->first);
    };
    
    /**
     * Returns the end time-iterator along the trajectory.
     * \return The end time-iterator along the trajectory.
     */
    point_time_iterator end_time_travel() const {
      return point_time_iterator(this);
    };
    
    
    
    struct point_fraction_iterator {
      const self* parent;
      time_difference_type interval_time;
      time_type current_time;
      time_type end_time;
      
      point_fraction_iterator(const self* aParent, time_difference_type aDeltaTime, time_type aCurrentTime) :
                              parent(aParent), interval_time(aDeltaTime), 
                              current_time(aCurrentTime), end_time(parent->get_current_horizon()) { };
      
      point_fraction_iterator(const self* aParent, time_difference_type aDeltaTime) :
                              parent(aParent), interval_time(aDeltaTime), 
                              current_time(parent->get_current_horizon()), end_time(current_time) { };
      
      friend point_fraction_iterator& operator+=(point_fraction_iterator& lhs, double rhs) {
        current_time += interval_time * rhs;
        if(current_time > end_time)
          current_time = end_time;
        return *this;
      };
      
      friend point_fraction_iterator operator+(point_fraction_iterator lhs, double rhs) { return (lhs += rhs); };
      friend point_fraction_iterator operator+(double lhs, point_fraction_iterator rhs) { return (rhs += lhs); };
      friend point_fraction_iterator operator-(point_fraction_iterator lhs, double rhs) { return (lhs += -rhs); };
      friend point_fraction_iterator& operator-=(point_fraction_iterator& lhs, double rhs) { return (lhs += -rhs); };
      
      friend bool operator==(const point_fraction_iterator& lhs, const point_fraction_iterator& rhs) {
        return ( ( lhs.parent == rhs.parent ) && ( lhs.current_time == rhs.current_time ) );
      };
      friend bool operator!=(const point_fraction_iterator& lhs, const point_fraction_iterator& rhs) { return !(lhs == rhs); };
      
      point_type operator*() const {
        return parent->get_point_at_time(current_time);
      };
      
    };
    
    /**
     * Returns the starting fraction-iterator along the trajectory.
     * \return The starting fraction-iterator along the trajectory.
     */
    point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(this, this->waypoints.begin()->second.predictor.get_time_step(), this->waypoints.begin()->first);
    };
    
    /**
     * Returns the end fraction-iterator along the trajectory.
     * \return The end fraction-iterator along the trajectory.
     */
    point_fraction_iterator end_fraction_travel() const {
      return point_fraction_iterator(this, this->waypoints.begin()->second.predictor.get_time_step());
    };
    
    
    
    
    
};


};


};


#endif









