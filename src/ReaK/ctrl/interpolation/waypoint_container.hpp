/**
 * \file waypoint_container.hpp
 * 
 * This library provides an implementation of a simple container of temporal waypoints.
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

#ifndef REAK_WAYPOINT_CONTAINER_HPP
#define REAK_WAYPOINT_CONTAINER_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include "path_planning/spatial_trajectory_concept.hpp"

#include "topologies/temporal_space.hpp"

#include "topologies/basic_distance_metrics.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>

namespace ReaK {

namespace pp {

  
/**
 * This class implements a simple container of waypoints.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points.
 */
template <typename Topology, typename DistanceMetric = typename metric_space_traits<Topology>::distance_metric_type >
class waypoint_container_base : public shared_object {
  public:
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,Topology>));
    
    typedef waypoint_container_base<Topology,DistanceMetric> self;
    typedef Topology topology;
    typedef DistanceMetric distance_metric;
    typedef typename topology_traits<Topology>::point_type point_type;
    typedef typename topology_traits<Topology>::point_difference_type point_difference_type;
    
    typedef std::list<point_type> container_type;
    
    typedef typename container_type::iterator waypoint_descriptor;
    typedef typename container_type::const_iterator const_waypoint_descriptor;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  protected:
    
    typedef typename container_type::value_type waypoint_value;
    
    static const point_type& make_value(const point_type& p) { return p; };
#ifdef RK_ENABLE_CXX11_FEATURES
    static point_type&& make_value(point_type&& p) { return std::move(p); };
#endif
    static const point_type& extract_point(const point_type& v) { return v; };
    static point_type& extract_point(point_type& v) { return v; };
    
    shared_ptr<topology> space;
    distance_metric dist;
    
    container_type waypoints;
    
    typedef std::pair<const_waypoint_descriptor, const_waypoint_descriptor> const_waypoint_bounds;
        
    const_waypoint_bounds get_waypoint_bounds(const point_type& p, const_waypoint_descriptor start) const {
      if(!space)
        return std::make_pair(waypoints.begin(),waypoints.begin());
      const_waypoint_descriptor it_end = waypoints.end();
      if(start == it_end)
        throw invalid_path("Waypoints exhausted during waypoint query!");
      double d1 = dist(*start,p,*space);
      const_waypoint_descriptor it2 = start; ++it2;
      if(it2 == it_end)
        return std::make_pair(start,start);
      double d12 = dist(*start,*it2,*space);
      while(d1 > d12) {
        d1 = dist(*it2, p, *space);
        if(d1 == std::numeric_limits<double>::infinity())
          return std::make_pair(start, it2);
        ++start; ++it2;
        if(it2 == it_end)
          return std::make_pair(start, start);
        d12 = dist(*start,*it2,*space);
      };
      return std::make_pair(start,it2);
      
//       double d2 = dist(p,*it2,*space);
//       double c12 = (d1 * d1 + d2 * d2 - d12 * d12) / (2.0 * d1 * d2); //cosine of the summet angle of the triangle.
//       const_waypoint_descriptor it3 = it2; ++it3;
// //       NOTE basically this loop will select the interval with the triangle (start,p,it2) with the largest angle at p.
//       while(it3 != it_end) {
//         double d3 = dist(*it3,p,*space);
//         double d23 = dist(*it2, *it3, *space);
//         double c23 = (d2 * d2 + d3 * d3 - d23 * d23) / (2.0 * d2 * d3);
//         if(c23 > c12) {
//           break;
//         } else {
//           start = it2; d1 = d2;
//           it2 = it3; d2 = d3;
//           c12 = c23;
//           ++it3;
//         };
//       };
//       return std::make_pair(start,it2);
    };
    
  public:
    /**
     * Constructs the waypoint-container from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    explicit waypoint_container_base(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), const distance_metric& aDist = distance_metric()) : 
                                     space(aSpace), 
                                     dist(aDist),
                                     waypoints() { 
    };
    
    /**
     * Constructs the waypoint-container from a space, the start and end points.
     * \param aSpace The space on which the waypoints are.
     * \param aStart The start point of the waypoints.
     * \param aEnd The end-point of the waypoints.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    waypoint_container_base(const shared_ptr<topology>& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                            space(aSpace), dist(aDist), waypoints() {
      waypoints.push_back(aStart);
      waypoints.push_back(aEnd);
    };
    
    /**
     * Constructs the waypoint-container from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the waypoints with.
     * \param aBegin An iterator to the first point of the waypoints.
     * \param aEnd An iterator to the last point of the waypoints.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    template <typename ForwardIter>
    waypoint_container_base(ForwardIter aBegin, ForwardIter aEnd, const shared_ptr<topology>& aSpace, const distance_metric& aDist = distance_metric()) : 
                            space(aSpace), dist(aDist), waypoints(aBegin,aEnd) {
      if(aBegin == aEnd)
        throw invalid_path("Empty list of waypoints!");
    };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& getSpace() const throw() { return *space; };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& get_temporal_space() const throw() { return *space; };
    
    /**
     * Returns the distance metric that the path uses.
     * \return The distance metric that the path uses.
     */
    const distance_metric& getDistanceMetric() const throw() { return dist; };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.space,rhs.space);
      swap(lhs.dist,rhs.dist);
      lhs.waypoints.swap(rhs.waypoints);
    };
    
    
    /**
     * Returns the starting point of the waypoints.
     * \return The starting point of the waypoints.
     */
    point_type get_start_point() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return *start;
    };
    
    /**
     * Returns the starting waypoint-point-pair of the waypoints.
     * \return The starting waypoint-point-pair of the waypoints.
     */
    waypoint_pair get_start_waypoint() const {
      const_waypoint_descriptor start = this->waypoints.begin();
      return waypoint_pair(start,*start);
    };
    
    /**
     * Returns the end point of the waypoints.
     * \return The end point of the waypoints.
     */
    point_type get_end_point() const {
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return *end;
    };
    
    /**
     * Returns the end waypoint-point-pair of the waypoints.
     * \return The end waypoint-point-pair of the waypoints.
     */
    waypoint_pair get_end_waypoint() const {
      const_waypoint_descriptor end = this->waypoints.end(); --end;
      return waypoint_pair(end,*end);
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(space)
        & RK_SERIAL_SAVE_WITH_NAME(dist)
        & RK_SERIAL_SAVE_WITH_NAME(waypoints);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(space)
        & RK_SERIAL_LOAD_WITH_NAME(dist)
        & RK_SERIAL_LOAD_WITH_NAME(waypoints);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440000,1,"waypoint_container_base",shared_object)
};



/**
 * This class implements a simple container of waypoints.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points.
 */
template <typename SpaceTopology, typename TimeTopology, typename DistanceMetric, typename DistanceMetricBase>
class waypoint_container_base< temporal_space<SpaceTopology, TimeTopology, DistanceMetric>, DistanceMetricBase > : public shared_object {
  public:
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept< temporal_space<SpaceTopology, TimeTopology, DistanceMetric> >));
    
    typedef waypoint_container_base<temporal_space<SpaceTopology, TimeTopology, DistanceMetric>,DistanceMetricBase> self;
    typedef temporal_space<SpaceTopology, TimeTopology, DistanceMetric> topology;
    typedef DistanceMetricBase distance_metric;
    typedef typename topology_traits< temporal_space<SpaceTopology, TimeTopology, DistanceMetric> >::point_type point_type;
    typedef typename topology_traits< temporal_space<SpaceTopology, TimeTopology, DistanceMetric> >::point_difference_type point_difference_type;
    
    typedef std::map<double, point_type> container_type;
    
    typedef typename container_type::iterator waypoint_descriptor;
    typedef typename container_type::const_iterator const_waypoint_descriptor;
    
    typedef std::pair<const_waypoint_descriptor, point_type> waypoint_pair;
    
  protected:
    
    typedef typename container_type::value_type waypoint_value;
    
    static waypoint_value make_value(const point_type& p) { return waypoint_value(p.time, p); };
#ifdef RK_ENABLE_CXX11_FEATURES
    static waypoint_value make_value(point_type&& p) { return waypoint_value(p.time, std::move(p)); };
#endif
    static const point_type& extract_point(const waypoint_value& v) { return v.second; };
    static point_type& extract_point(waypoint_value& v) { return v.second; };
    
    shared_ptr<topology> space;
    distance_metric dist;
    
    container_type waypoints;
    
    typedef std::pair<const_waypoint_descriptor, const_waypoint_descriptor> const_waypoint_bounds;
    
    const_waypoint_bounds get_waypoint_bounds(double t) const {
      const_waypoint_descriptor it2 = waypoints.lower_bound(t);
      if(it2 == waypoints.begin())
        return const_waypoint_bounds(it2,it2);
      if(it2 == waypoints.end())
        return const_waypoint_bounds((++waypoints.rbegin()).base(),(++waypoints.rbegin()).base());
      const_waypoint_descriptor it1 = it2; --it1;
      return const_waypoint_bounds(it1,it2);
    };
    
    virtual double travel_distance_impl(const point_type& a, const point_type& b) const {
      return 0.0;
    };
    
  public:
    /**
     * Constructs the waypoint-container from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    explicit waypoint_container_base(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), const distance_metric& aDist = distance_metric()) : 
                                     space(aSpace), 
                                     dist(aDist),
                                     waypoints() { 
    };
    
    /**
     * Constructs the waypoint-container from a space, the start and end points.
     * \param aSpace The space on which the waypoints are.
     * \param aStart The start point of the waypoints.
     * \param aEnd The end-point of the waypoints.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    waypoint_container_base(const shared_ptr<topology>& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                            space(aSpace), dist(aDist), waypoints() {
      waypoints.insert(waypoint_value(aStart.time, aStart));
      waypoints.insert(waypoints.end(), waypoint_value(aEnd.time, aEnd));
    };
    
    /**
     * Constructs the waypoint-container from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the waypoints with.
     * \param aBegin An iterator to the first point of the waypoints.
     * \param aEnd An iterator to the last point of the waypoints.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    template <typename ForwardIter>
    waypoint_container_base(ForwardIter aBegin, ForwardIter aEnd, const shared_ptr<topology>& aSpace, const distance_metric& aDist = distance_metric()) : 
                            space(aSpace), dist(aDist), waypoints() {
      if(aBegin == aEnd)
        throw invalid_path("Empty list of waypoints!");
      while(aBegin != aEnd) {
        waypoints.insert(waypoints.end(), waypoint_value(aBegin->time, *aBegin));
        ++aBegin;
      };
    };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& getSpace() const throw() { return *space; };
    
    /**
     * Returns the space on which the path resides.
     * \return The space on which the path resides.
     */
    const topology& get_temporal_space() const throw() { return *space; };
    
    /**
     * Returns the distance metric that the path uses.
     * \return The distance metric that the path uses.
     */
    const distance_metric& getDistanceMetric() const throw() { return dist; };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.space,rhs.space);
      swap(lhs.dist,rhs.dist);
      lhs.waypoints.swap(rhs.waypoints);
    };
    
    
    /**
     * Returns the starting time of the waypoints.
     * \return The starting time of the waypoints.
     */
    double get_start_time() const {
      return waypoints.begin()->first;
    };
    
    /**
     * Returns the end time of the waypoints.
     * \return The end time of the waypoints.
     */
    double get_end_time() const {
      const_waypoint_descriptor end = waypoints.end(); --end;
      return end->first;
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
     * Computes the travel distance between two points, if traveling along the trajectory.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the trajectory.
     */
    double travel_distance(const point_type& a, const point_type& b) const {
      return this->travel_distance_impl(a, b);
    };
    
    /**
     * Computes the travel distance between two waypoint-point-pairs, if traveling along the trajectory.
     * \param a The first waypoint-point-pair.
     * \param b The second waypoint-point-pair.
     * \return The travel distance between two points if traveling along the trajectory.
     */
    double travel_distance(const waypoint_pair& a, const waypoint_pair& b) const {
      return this->travel_distance_impl(a.second, b.second);
    };
    
    /**
     * Returns the total travel-distance of the trajectory.
     * \return The total travel-distance of the trajectory.
     */
    double get_total_length() const {
      return this->travel_distance_impl(this->get_start_point(), this->get_end_point());
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(space)
        & RK_SERIAL_SAVE_WITH_NAME(dist)
        & RK_SERIAL_SAVE_WITH_NAME(waypoints);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(space)
        & RK_SERIAL_LOAD_WITH_NAME(dist)
        & RK_SERIAL_LOAD_WITH_NAME(waypoints);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440000,1,"waypoint_container_base",shared_object)
    
};







/**
 * This class implements a simple container of waypoints.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class waypoint_container : public waypoint_container_base<Topology,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,Topology>));
    
    typedef waypoint_container<Topology,DistanceMetric> self;
    typedef waypoint_container_base<Topology,DistanceMetric> base_class_type;
    
    typedef typename base_class_type::container_type container_type;
    typedef typename base_class_type::topology topology;
    typedef typename base_class_type::distance_metric distance_metric;
    typedef typename base_class_type::point_type point_type;
    typedef typename base_class_type::const_waypoint_descriptor const_waypoint_descriptor;
    typedef typename base_class_type::waypoint_descriptor waypoint_descriptor;
    
    typedef typename base_class_type::waypoint_pair waypoint_pair;
    
    typedef typename container_type::size_type size_type;
    typedef typename container_type::value_type value_type;
    
  public:
    /**
     * Constructs the waypoint-container from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    explicit waypoint_container(const shared_ptr<topology>& aSpace = shared_ptr<topology>(new topology()), const distance_metric& aDist = distance_metric()) : 
                                base_class_type(aSpace,aDist) { };
    
    /**
     * Constructs the waypoint-container from a space, the start and end points.
     * \param aSpace The space on which the waypoints are.
     * \param aStart The start point of the waypoints.
     * \param aEnd The end-point of the waypoints.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    waypoint_container(const shared_ptr<topology>& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                       base_class_type(aSpace,aStart,aEnd,aDist) { };

    /**
     * Constructs the waypoint-container from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the waypoints with.
     * \param aBegin An iterator to the first point of the waypoints.
     * \param aEnd An iterator to the last point of the waypoints.
     * \param aSpace The space on which the waypoints are.
     * \param aDist The distance metric functor that the waypoint-container should use.
     */
    template <typename ForwardIter>
    waypoint_container(ForwardIter aBegin, ForwardIter aEnd, const shared_ptr<topology>& aSpace, const distance_metric& aDist = distance_metric()) : 
                       base_class_type(aBegin, aEnd, aSpace, aDist) { };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(static_cast<base_class_type&>(lhs),static_cast<base_class_type&>(rhs));
    };
    
    /* **************************************************************
     *                   STL container interface
     * ************************************************************** */
    
    const_waypoint_descriptor begin() const {
      return this->waypoints.begin();
    };
    
    const_waypoint_descriptor end() const {
      return this->waypoints.end();
    };
    
    typename container_type::const_reverse_iterator rbegin() const {
      return this->waypoints.rbegin();
    };
    
    typename container_type::const_reverse_iterator rend() const {
      return this->waypoints.rend();
    };
    
    waypoint_descriptor begin() {
      return this->waypoints.begin();
    };
    
    waypoint_descriptor end() {
      return this->waypoints.end();
    };
    
    typename container_type::reverse_iterator rbegin() {
      return this->waypoints.rbegin();
    };
    
    typename container_type::reverse_iterator rend() {
      return this->waypoints.rend();
    };
    
    size_type size() const { 
      return this->waypoints.size();
    };
    
    bool empty() const { return this->waypoints.empty(); };
    
    size_type max_size() const { 
      return this->waypoints.max_size();
    };
    
    const point_type& front() const {
      return base_class_type::extract_point( *(this->waypoints.begin()) );
    };
    
    const point_type& back() const {
      return base_class_type::extract_point( *(this->waypoints.rbegin()) );
    };
    
    point_type& front() {
      return base_class_type::extract_point( *(this->waypoints.begin()) );
    };
    
    point_type& back() {
      return base_class_type::extract_point( *(this->waypoints.rbegin()) );
    };
    
    template <typename InputIterator>
    void assign( InputIterator first, InputIterator last ) {
      waypoint_container tmp(first, last, this->space, this->dist);
      swap(tmp, *this);
    };
    
    void push_front(const point_type& p) {
      this->waypoints.insert(this->waypoints.begin(), base_class_type::make_value(p));
    };
    
    void push_back(const point_type& p) {
      this->waypoints.insert(this->waypoints.end(), base_class_type::make_value(p));
    };
    
    void pop_front() {
      this->waypoints.erase(this->waypoints.begin());
    };
    
    void pop_back() {
      this->waypoints.erase((++(this->waypoints.rbegin())).base());
    };
    
    waypoint_descriptor insert( waypoint_descriptor position, const point_type& p) {
      return this->waypoints.insert(position, base_class_type::make_value(p));
    };
    
    template <typename InputIterator>
    void insert( waypoint_descriptor position, InputIterator first, InputIterator last) {
      for(; first != last; ++first)
        this->waypoints.insert(position, base_class_type::make_value(*first));
    };
    
    void erase( waypoint_descriptor position) {
      if((position == this->waypoints.begin()) && (this->waypoints.size == 1))
        throw invalid_path("Cannot empty the list of waypoints!");
      this->waypoints.erase(position);
    };
    
    void erase( waypoint_descriptor first, waypoint_descriptor last) {
      if((first == this->waypoints.begin()) && (last == this->waypoints.end()))
        throw invalid_path("Cannot empty the list of waypoints!");
      this->waypoints.erase(first, last);
    };
    
    typename container_type::allocator_type get_allocator() const { return this->waypoints.get_allocator(); };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_class_type::save(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_class_type::load(A,base_class_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440001,1,"waypoint_container",base_class_type)
    
    
};




};

};

#endif









