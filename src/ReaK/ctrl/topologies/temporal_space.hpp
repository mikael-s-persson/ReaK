/**
 * \file temporal_space.hpp
 * 
 * This library provides an implementation of a temporal-space which augments a 
 * topology with a temporal dimension (time-stamp).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_TEMPORAL_SPACE_HPP
#define REAK_TEMPORAL_SPACE_HPP

#include <boost/config.hpp>
#include <cmath>
#include <boost/random/uniform_01.hpp>
#include <boost/random/linear_congruential.hpp>

#include "time_topology.hpp"
#include "path_planning/temporal_space_concept.hpp"

#include "base/named_object.hpp"

namespace ReaK {

namespace pp {


/**
 * This class implementats a temporal-space which augments a 
 * topology with a temporal dimension (time-stamp). The time-dimension resides on a line-segment 
 * topology (see line_segment_topology), while the spatial topology and distance-metric 
 * is provided by the user. Models the TemporalSpaceConcept.
 * \tparam Topology The topology type which represents the spatial dimensions, should model MetricSpaceConcept.
 * \tparam TimeTopology The topology type which represents the time dimension, should model MetricSpaceConcept.
 * \tparam TemporalDistanceMetric The distance metric type for the temporal-space, should model the TemporalDistMetricConcept.
 */
template <typename Topology, typename TimeTopology, typename TemporalDistanceMetric = spatial_distance_only>
class temporal_space : public named_object {
  public:
    typedef TimeTopology time_topology;
    typedef Topology space_topology;
    typedef TemporalDistanceMetric distance_metric;
    
    typedef temporal_space<Topology,TimeTopology,TemporalDistanceMetric> self;
    
    /**
     * This nested type represents the points of the temporal-space.
     */
    struct point {
      typename metric_topology_traits<time_topology>::point_type time; ///< Holds the time associated to the space-time point.
      typename metric_topology_traits<space_topology>::point_type pt; ///< Holds the spatial-point associated to the space-time point.
    
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = metric_topology_traits<space_topology>::dimensions);
      
      /**
       * Default constructor.
       */
      point() { };
      /**
       * Constructor from a time and spatial point.
       * \param aTime The time associated to the space-time point.
       * \param aPt The spatial-point associated to the space-time point.
       */
      point(const typename metric_topology_traits<time_topology>::point_type& aTime, 
	    const typename metric_topology_traits<space_topology>::point_type& aPt) : 
	    time(aTime), pt(aPt) { };
     
      /* This should not be there, because of the principle of minimum requirements.
      double& operator[](std::size_t i) { return pt[i]; };
      const double& operator[](std::size_t i) const { return pt[i]; };
      */
    };
  
    /**
     * This nested type represents the difference between two points of the temporal-space.
     */
    struct point_difference {
      typename metric_topology_traits<time_topology>::point_difference_type time; ///< Holds the time-difference.
      typename metric_topology_traits<space_topology>::point_difference_type pt; ///< Holds the spatial-difference.
      
      BOOST_STATIC_CONSTANT(std::size_t, dimensions = metric_topology_traits<space_topology>::dimensions);
      /**
       * Default constructor.
       */
      point_difference() : time(), pt() { };
      /**
       * Constructor from a time and space difference.
       * \param aTime The time difference.
       * \param aPt The spatial-difference.
       */
      point_difference(const typename metric_topology_traits<time_topology>::point_difference_type& aTime,
	               const typename metric_topology_traits<space_topology>::point_difference_type& aPt) : 
	               time(aTime), pt(aPt) { };

      point_difference operator-() const {
        return point_difference(-time,-pt);
      };

      friend point_difference operator*(const point_difference& a, double b) {
        return point_difference(a.time * b, a.pt * b);
      };

      friend point_difference operator*(double a, const point_difference& b) {
        return point_difference(a * b.time, a * b.pt);
      };
      
      point_difference& operator+=(const point_difference& b) {
        time += b.time;
	pt += b.pt;
        return *this;
      };

      point_difference& operator-=(const point_difference& b) {
        time -= b.time;
	pt -= b.pt;
        return *this;
      };
      
      friend point_difference operator+(const point_difference& a, const point_difference& b) {
        point_difference result;
	result.time = a.time + b.time;
	result.pt = a.pt + b.pt;
        return result;
      };

      friend point_difference operator-(const point_difference& a, const point_difference& b) {
        point_difference result;
	result.time = a.time - b.time;
	result.pt = a.pt - b.pt;
        return result;
      };

      /* This should not be there, because of the principle of minimum requirements.
      double& operator[](std::size_t i) { return pt[i]; };
      const double& operator[](std::size_t i) const { return pt[i]; };
      
      friend double dot(const point_difference& a, const point_difference& b) {
        return dot(a.pt, b.pt);
      };*/

    };
  protected:
    space_topology space;
    time_topology time;
    distance_metric dist;
        
    temporal_space(const temporal_space&); //non-copyable.
    temporal_space& operator=(const temporal_space&);
    
  public:
    /**
     * Parametrized constructor.
     * \param aSpace The space topology to be used.
     * \param aTime The time topology to be used.
     * \param aDist The temporal distance metric functor to use.
     */
    explicit temporal_space(const std::string& aName = "",
			    const space_topology& aSpace = space_topology(), 
			    const time_topology& aTime = time_topology(),
			    const distance_metric& aDist = distance_metric()) :
                            named_object(), space(aSpace), time(aTime), dist(aDist) { this->setName(aName); };
    
    typedef point point_type;
    typedef point_difference point_difference_type;
    
    /** Returns the underlying space topology. */
    const space_topology& get_space_topology() const { return space; };
    /** Returns the underlying time topology. */
    const time_topology& get_time_topology() const { return time; };
    /** Returns the temporal distance metric functor used. */
    const distance_metric& get_distance_metric() const { return dist; };
    
    /** Returns the underlying space topology. */
    space_topology& get_space_topology() { return space; };
    /** Returns the underlying time topology. */
    time_topology& get_time_topology() { return time; };
    /** Returns the temporal distance metric functor used. */
    distance_metric& get_distance_metric() { return dist; };

    /**
     * Returns a random point within the temporal-space.
     * \return A random point within the temporal-space.
     */
    point random_point() const 
    {
      return point(time.random_point(), space.random_point());
    }
    
    /**
     * Computes the distance between two points in the temporal-space.
     * \param a The first point.
     * \param b The second point.
     * \return The distance between a and b.
     */
    double distance(const point& a, const point& b) const {
      return dist(a, b, time, space);
    };

    /**
     * Returns a point which is at a fraction between two points a to b.
     * \param a The first point.
     * \param fraction The fraction between the two points (0 to 1).
     * \param b The second point.
     * \return The point which is at a fraction between two points.
     */
    point move_position_toward(const point& a, double fraction, const point& b) const {
      return point(time.move_position_toward(a.time, fraction, b.time),
	           space.move_position_toward(a.pt, fraction, b.pt));
    };

    /**
     * Returns the difference between two points (a - b).
     * \param a The first point.
     * \param b The second point.
     * \return The difference between the two points.
     */
    point_difference difference(const point& a, const point& b) const {
      return point_difference(time.difference(a.time, b.time), space.difference(a.pt, b.pt));
    };

    /**
     * Returns the addition of a point-difference to a point.
     * \param a The starting point.
     * \param delta The point-difference.
     * \return The addition of a point-difference to a point.
     */
    point adjust(const point& a, const point_difference& delta) const {
      return point(time.adjust(a.time, delta.time), space.adjust(a.pt, delta.pt));
    };
  
    /**
     * Returns the origin of the temporal-space.
     * \return The origin of the temporal-space.
     */
    point origin() const {
      return point(time.origin(), space.origin());
    };

    /**
     * Computes the norm of a difference between two points.
     * \param a The difference between two points.
     * \return The norm of a difference between two points.
     */
    double norm(const point_difference& a) const {
      return dist(a, time, space);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(space)
        & RK_SERIAL_SAVE_WITH_NAME(time)
	& RK_SERIAL_SAVE_WITH_NAME(dist);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(space)
        & RK_SERIAL_LOAD_WITH_NAME(time)
	& RK_SERIAL_LOAD_WITH_NAME(dist);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400004,1,"temporal_space",named_object)
    
};



};

};

#endif
















