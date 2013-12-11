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
#include "topologies/time_poisson_topology.hpp"

#include "interpolation/waypoint_container.hpp"

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
class belief_predicted_trajectory : public pp::waypoint_container< pp::temporal_space<BeliefTopology, pp::time_poisson_topology, pp::time_distance_only> > {
  public:
    
    typedef belief_predicted_trajectory<BeliefTopology,BeliefPredictor,InputTrajectory> self;
    typedef pp::temporal_space<BeliefTopology, pp::time_poisson_topology, pp::time_distance_only> topology;
    typedef pp::waypoint_container< topology > base_class_type;
    typedef BeliefPredictor predictor_type;
    
    typedef shared_ptr<topology> topology_ptr;
    
    typedef typename base_class_type::waypoint_descriptor waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::const_waypoint_bounds const_waypoint_bounds;
    typedef typename base_class_type::distance_metric distance_metric;
    
    typedef typename base_class_type::waypoint_pair waypoint_pair;
    
    typedef typename pp::topology_traits< topology >::point_type point_type;
    typedef typename pp::topology_traits< topology >::point_difference_type point_difference_type;
    
    typedef typename pp::temporal_space_traits<topology>::time_topology time_topology;
    typedef typename pp::topology_traits<time_topology>::point_type time_type;
    typedef typename pp::topology_traits<time_topology>::point_difference_type time_difference_type;
    
    typedef typename pp::temporal_space_traits<topology>::space_topology space_topology;
    typedef typename pp::topology_traits<space_topology>::point_type belief_state;
    typedef typename pp::topology_traits<space_topology>::point_difference_type belief_state_diff;
    
    typedef typename discrete_sss_traits< predictor_type >::input_type input_type;
    
    BOOST_CONCEPT_ASSERT((pp::TemporalSpaceConcept<topology>));
    BOOST_CONCEPT_ASSERT((BeliefSpaceConcept<space_topology>));
    BOOST_CONCEPT_ASSERT((BeliefPredictorConcept<BeliefPredictor, space_topology>));
    BOOST_CONCEPT_ASSERT((SpatialTrajectoryConcept<InputTrajectory, pp::vector_topology<input_type> >));
    
    typedef std::map< const point_type*, predictor_type > predictor_map_type;
    
  private:
    InputTrajectory input;
    mutable predictor_map_type pred_segments;
    
    waypoint_descriptor updated_end;
    
    virtual double travel_distance_impl(const point_type& a, const point_type& b) const {
      using std::fabs;
      return fabs(b.time - a.time);
    };
    
    waypoint_pair get_point_at_time_impl(double t, const_waypoint_bounds wpb_a) const {
      if( t > get_current_horizon() ) {
        set_minimal_horizon(t);
        wpb_a.second = updated_end;
        --wpb_a.second;
        wpb_a.first = wpb_a.second;
        if( wpb_a.second->first > t )
          --wpb_a.first;
      };
      
      if( ( wpb_a.first == wpb_a.second ) || ( (t - wpb_a.first->first) <= (wpb_a.second->first - t) ) )
        return waypoint_pair(wpb_a.first, point_type(t, wpb_a.first->second.pt));
      else
        return waypoint_pair(wpb_a.first, point_type(t, wpb_a.second->second.pt));
    };
    
    virtual waypoint_pair move_time_diff_from_impl(const point_type& a, const const_waypoint_bounds& wpb_a, double dt) const {
      if( ( a.time + dt >= wpb_a.first->first ) && ( a.time + dt <= wpb_a.second->first ) )
        return get_point_at_time_impl(a.time + dt, wpb_a);
      else
        return get_point_at_time_impl(a.time + dt, this->get_waypoint_bounds(a.time + dt));
    };
    
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
      swap(static_cast<base_class_type&>(lhs),static_cast<base_class_type&>(rhs));
      lhs.updated_end = lhs.waypoints.lower_bound(lhs_horiz);
      rhs.updated_end = rhs.waypoints.lower_bound(rhs_horiz);
      
      swap(lhs.pred_segments, rhs.pred_segments);
      swap(lhs.input, rhs.input);
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
      
      waypoint_descriptor it = this->waypoints.lower_bound(t);
      
      waypoint_descriptor last_valid = updated_end; --last_valid;
      while( ( updated_end != this->waypoints.end() ) && ( last_valid->first < t ) ) {
        updated_end->second.pt = pred_segments[&(last_valid->second)].predict_belief(
          this->space->get_space_topology(), 
          last_valid->second.pt, last_valid->first, 
          input.get_point_at_time(last_valid->first));
        ++last_valid; ++updated_end;
      };
      
      if(it == this->waypoints.end()) {
        // must perform predictions from end to past time t.
        while( last_valid->first < t ) {
          const predictor_type& pred = pred_segments[&(last_valid->second)];
          this->push_back( point_type( 
            last_valid->first + pred.get_time_step(), 
            pred.predict_belief(
              space->get_space_topology(), 
              last_valid->second.pt, last_valid->first, 
              input.get_point_at_time(last_valid->first))));
          ++last_valid;
          pred_segments[&(last_valid->second)] = pred;
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
      
      waypoint_descriptor it = this->waypoints.lower_bound(t);
      ++it; // <- this means that the pruning "conservative", and also, will always keep the beginning iterator.
      for(waypoint_descriptor it2 = it; it2 != this->waypoints.end(); ++it2)
        pred_segments.erase(&(it2->second));
      this->erase(it, this->waypoints.end());
      updated_end = this->waypoints.end();
      
    };
    
    
    /**
     * Resets the initial point of the predicted trajector with a point and an associated time. 
     * This function triggers the elimination of all waypoints prior to this time and triggers 
     * the recomputation of the belief predictions of all waypoints after this time.
     * \param pt The point that will become the first, initial waypoint.
     */
    void set_start_point(const point_type& pt, const predictor_type& pred) {
      using std::fabs;
      
      waypoint_descriptor it = this->waypoints.lower_bound(pt.time - 0.5 * pred.get_time_step());
      // check if the time difference is too much:
      if( fabs(pt.time - it->first) > 1e-4 * pred.get_time_step() ) {
        // we have to completely reset the entire container:
        pred_segments.clear();
        this->waypoints.clear();
        this->push_back(pt);
        updated_end = this->waypoints.begin();
        pred_segments[&(updated_end->second)] = pred;
        ++updated_end;
        return;
      };
      
      // trim away the starting part of the trajectory:
      for(waypoint_descriptor it2 = this->waypoints.begin(); it2 != it; ++it2)
        pred_segments.erase(&(it2->second));
      this->waypoints.erase(waypoints.begin(), it);
      // reset the start point and predictor, and reset the horizon:
      pred_segments[&(it->second)] = pred;
      it->second = pt;
      updated_end = it; ++updated_end;
    };
    
    
    
    
    
    
    
};


};


};


#endif









