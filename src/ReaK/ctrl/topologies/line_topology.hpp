/**
 * \file line_topology.hpp
 * 
 * This library provides classes that define a line-topology. A line-topology is 
 * a simple metric-space where the points are real values (doubles) along a 1D 
 * space (line-segment).
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

#ifndef REAK_LINE_TOPOLOGY_HPP
#define REAK_LINE_TOPOLOGY_HPP


#include "base/defs.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "path_planning/global_rng.hpp"

#include <cmath>
#include "base/named_object.hpp"

#include "basic_distance_metrics.hpp"
#include "default_random_sampler.hpp"

namespace ReaK {

namespace pp {

/**
 * This class implements an infinite line topology. This class models the TopologyConcept,
 * the LieGroupConcept, and the MetricSpaceConcept.
 * \tparam T The value-type for the topology (should be an arithmetic type that is implicitly convertable to double).
 */
template <typename T = double>
class line_topology : public named_object
{
  public:
    typedef line_topology<T> self;
    
    typedef T point_type;
    typedef T point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 1);
    
    line_topology(const std::string& aName = "line_topology") : named_object() {
      setName(aName);
    };

    
   /*************************************************************************
    *                             TopologyConcept
    * **********************************************************************/
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& a, const point_type& b) const {
      return a - b;
    }

    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const {
      return a + delta;
    }

    /**
     * Returns the origin of the space (the lower-limit).
     */
    virtual point_type origin() const {
      return point_type(0);
    };
    
    /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      using std::fabs;
      return fabs(b - a);
    }

    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      using std::fabs;
      return fabs(delta);
    }

    /**
     * Returns the volume of the difference between two points.
     */
    double volume(const point_difference_type& delta) const {
      using std::fabs;
      return fabs(delta);
    }

    /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const 
    {
      return a + (b - a) * fraction;
    }

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400001,1,"line_topology",named_object)
    
};

/**
 * This class implements a line-segment topology. The space extends from the minimum value up to some 
 * maximum value. Models the TopologyConcept, LieGroupConcept, MetricSpaceConcept, 
 * BoundedSpaceConcept, and SphereBoundedSpaceConcept, and also provides a distance-metric and a 
 * random-sampler (default, uniform sampler).
 * 
 * \tparam T A value-type (scalar value).
 */
template<typename T = double>
class line_segment_topology : public line_topology<T>
{
  typedef boost::uniform_01<global_rng_type&, T> rand_t;

  public:
    typedef line_segment_topology<T> self;
    
    typedef typename line_topology<T>::point_type point_type;
    typedef typename line_topology<T>::point_difference_type point_difference_type;
    
    typedef typename line_topology<T>::distance_metric_type distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = line_topology<T>::dimensions);

    /**
     * Default constructor.
     * \param aStart The minimum bound of the line-segment.
     * \param aEnd The maximum bound of the line-segment.
     */
    explicit line_segment_topology(const std::string& aName = "line_segment_topology", 
				   point_type aStart = point_type(0.0), 
				   point_type aEnd = point_type(1.0)) 
      : line_topology<T>(aName), start_pt(aStart), end_pt(aEnd) { };
   
    
    /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      return rand_t(get_global_rng())() * (end_pt - start_pt) + start_pt;
    };
    
    /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/

    /**
     * Takes a point and clips it to within this line-segment space.
     */
    point_type bound(point_type a) const {
      if(end_pt > start_pt) {
        if(a > end_pt)
  	  return end_pt;
        else if(a < start_pt) 
 	  return start_pt;
        else
	  return a;
      } else {
        if(a < end_pt)
	  return end_pt;
        else if(a > start_pt)
	  return start_pt;
        else
	  return a;
      };
    };

    /**
     * Returns the distance to the boundary of the space.
     */
    double distance_from_boundary(point_type a) const {
      using std::fabs;
      double dist = fabs(end_pt - a);
      if(fabs(a - start_pt) < dist)
        return fabs(a - start_pt);
      else
        return dist;
    };
    
    /**
     * Returns the difference to the closest boundary.
     */
    point_difference_type get_diff_to_boundary(point_type a) const {
      using std::fabs;
      double dist = fabs(end_pt - a);
      if(fabs(a - start_pt) < dist)
        return start_pt - a;
      else
        return end_pt - a;
    };
      
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      if(end_pt > start_pt) {
        if((a > end_pt) || (a < start_pt))
  	  return false;
      } else {
        if((a < end_pt) || (a > start_pt))
	  return false;
      };
      return true;
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return (end_pt - start_pt) * 0.5 + start_pt;
    };
    
    /*************************************************************************
    *                             SphereBoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the radius of the space.
     */
    double get_radius() const {
      using std::fabs;
      return fabs(end_pt - start_pt) * 0.5;
    };

  private:
    point_type start_pt;
    point_type end_pt;
    
  public:
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      line_topology<T>::save(A,line_topology<T>::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(start_pt)
        & RK_SERIAL_SAVE_WITH_NAME(end_pt);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      line_topology<T>::load(A,line_topology<T>::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(start_pt)
        & RK_SERIAL_LOAD_WITH_NAME(end_pt);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400006,1,"line_segment_topology",line_topology<T>)
};


};

};

#endif








