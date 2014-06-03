/**
 * \file path_base.hpp
 * 
 * This library provides the base-class for paths within a topology.
 * This is a base-class that stems the object-oriented compatibility of other 
 * path classes. Then, this library provides a path-wrapper class template 
 * which makes an OOP-compatible path class for a given path.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PATH_BASE_HPP
#define REAK_PATH_BASE_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>

#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  

/**
 * This class defines the OOP interface for a path in a topology.
 * \tparam Topology The topology type on which the points and the path can reside, should model the MetricSpaceConcept.
 */
template <typename Topology>
class path_base : public named_object {
  public:
    
    BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Topology>));
    
    typedef Topology topology;
    typedef typename topology_traits<topology>::point_type point_type;
    typedef path_base<Topology> self;
    
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     */
    explicit path_base(const std::string& aName = "") : 
                       named_object() { 
      setName(aName);
    };
    
    /**
     * Computes the travel distance between two points, if traveling along the path.
     * \param a The first point.
     * \param b The second point.
     * \return The travel distance between two points if traveling along the path.
     */
    virtual double travel_distance(const point_type& a, const point_type& b) const = 0;
    
    /**
     * Computes the point that is a distance away from a point on the path.
     * \param a The point on the path.
     * \param d The distance to move away from the point.
     * \return The point that is a distance away from the given point.
     */
    virtual point_type move_away_from(const point_type& a, double d) const = 0;
    
    /**
     * Returns the starting point of the trajectory.
     * \return The starting point of the trajectory.
     */
    virtual point_type get_start_point() const = 0;
    
    /**
     * Returns the end point of the trajectory.
     * \return The end point of the trajectory.
     */
    virtual point_type get_end_point() const = 0;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC244000C,1,"path_base",named_object)
};



};

};

#endif









