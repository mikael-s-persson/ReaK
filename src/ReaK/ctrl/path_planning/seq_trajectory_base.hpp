/**
 * \file seq_trajectory_base.hpp
 * 
 * This library provides the base-class for sequential trajectories within a temporal topology.
 * This is a base-class that stems the object-oriented compatibility of other 
 * sequential trajectory classes. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_SEQ_TRAJECTORY_BASE_HPP
#define REAK_SEQ_TRAJECTORY_BASE_HPP

#include "base/defs.hpp"

#include "base/named_object.hpp"

#include "temporal_space_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  

/**
 * This class defines the OOP interface for a sequential trajectory in a temporal topology.
 * \tparam Topology The temporal topology type on which the points and the sequential trajectory can reside, should model the TemporalSpaceConcept.
 */
template <typename Topology>
class seq_trajectory_base : public named_object {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    
    typedef Topology topology;
    typedef typename topology_traits<topology>::point_type point_type;
    typedef seq_trajectory_base<Topology> self;
    
  protected:
    
    struct point_time_iterator_impl {
      
      virtual ~point_time_iterator_impl() { };
      
      virtual void move_by_time(double t) = 0;
      
      virtual bool is_equal_to(const point_time_iterator_impl* rhs) const = 0;
      
      virtual const point_type& get_point() const = 0;
      
      virtual point_time_iterator_impl* clone() const = 0;
      
    };
    
    
    struct point_fraction_iterator_impl {
      
      virtual ~point_fraction_iterator_impl() { };
      
      virtual void move_by_fraction(double f) = 0;
      
      virtual bool is_equal_to(const point_fraction_iterator_impl* rhs) const = 0;
      
      virtual const point_type& get_point() const = 0;
      
      virtual point_fraction_iterator_impl* clone() const = 0;
      
    };
    
    
  public:
    
    class point_time_iterator {
      private:
        point_time_iterator_impl* p_impl;
        
      public:
        
        point_time_iterator(point_time_iterator_impl* aPImpl) : p_impl(aPImpl) { };
        
        point_time_iterator(const point_time_iterator& rhs) : p_impl(rhs.p_impl->clone()) { };
#ifdef RK_ENABLE_CXX11_FEATURES
        point_time_iterator(point_time_iterator&& rhs) : p_impl(rhs.p_impl) { rhs.p_impl = NULL; };
#endif
        friend 
        void swap(point_time_iterator& rhs, point_time_iterator& lhs) {
          std::swap(rhs.p_impl, lhs.p_impl);
        };
        point_time_iterator& operator=(point_time_iterator rhs) {
          swap(*this,rhs);
          return *this;
        };
        ~point_time_iterator() { delete p_impl; };
        
        friend 
        point_time_iterator operator+(point_time_iterator lhs, double rhs) {
          lhs.p_impl->move_by_time(rhs);
          return lhs;
        };
        
        friend 
        point_time_iterator operator+(double lhs, point_time_iterator rhs) {
          rhs.p_impl->move_by_time(lhs);
          return rhs;
        };
        
        friend
        point_time_iterator& operator+=(point_time_iterator& lhs, double rhs) {
          lhs.p_impl->move_by_time(rhs);
          return lhs;
        };
        
        friend 
        point_time_iterator operator-(point_time_iterator lhs, double rhs) {
          lhs.p_impl->move_by_time(-rhs);
          return lhs;
        };
        
        friend
        point_time_iterator& operator-=(point_time_iterator& lhs, double rhs) {
          lhs.p_impl->move_by_time(-rhs);
          return lhs;
        };
        
        friend 
        bool operator==(const point_time_iterator& lhs, 
                        const point_time_iterator& rhs) {
          return lhs.p_impl->is_equal_to(rhs.p_impl);
        };
        
        friend 
        bool operator!=(const point_time_iterator& lhs, 
                        const point_time_iterator& rhs) {
          return !(lhs.p_impl->is_equal_to(rhs.p_impl));
        };
        
        const point_type& operator*() const {
          return p_impl->get_point();
        };
        
    };
    
    
    class point_fraction_iterator {
      private:
        point_fraction_iterator_impl* p_impl;
        
      public:
        
        point_fraction_iterator(point_fraction_iterator_impl* aPImpl) : p_impl(aPImpl) { };
        
        point_fraction_iterator(const point_fraction_iterator& rhs) : p_impl(rhs.p_impl->clone()) { };
#ifdef RK_ENABLE_CXX11_FEATURES
        point_fraction_iterator(point_fraction_iterator&& rhs) : p_impl(rhs.p_impl) { rhs.p_impl = NULL; };
#endif
        friend 
        void swap(point_fraction_iterator& rhs, point_fraction_iterator& lhs) {
          std::swap(rhs.p_impl, lhs.p_impl);
        };
        point_fraction_iterator& operator=(point_fraction_iterator rhs) {
          swap(*this,rhs);
          return *this;
        };
        ~point_fraction_iterator() { delete p_impl; };
        
        friend 
        point_fraction_iterator operator+(point_fraction_iterator lhs, double rhs) {
          lhs.p_impl->move_by_fraction(rhs);
          return lhs;
        };
        
        friend 
        point_fraction_iterator operator+(double lhs, point_fraction_iterator rhs) {
          rhs.p_impl->move_by_fraction(lhs);
          return rhs;
        };
        
        friend
        point_fraction_iterator& operator+=(point_fraction_iterator& lhs, double rhs) {
          lhs.p_impl->move_by_fraction(rhs);
          return lhs;
        };
        
        friend 
        point_fraction_iterator operator-(point_fraction_iterator lhs, double rhs) {
          lhs.p_impl->move_by_fraction(-rhs);
          return lhs;
        };
        
        friend
        point_fraction_iterator& operator-=(point_fraction_iterator& lhs, double rhs) {
          lhs.p_impl->move_by_fraction(-rhs);
          return lhs;
        };
        
        friend 
        bool operator==(const point_fraction_iterator& lhs, 
                        const point_fraction_iterator& rhs) {
          return lhs.p_impl->is_equal_to(rhs.p_impl);
        };
        
        friend 
        bool operator!=(const point_fraction_iterator& lhs, 
                        const point_fraction_iterator& rhs) {
          return !(lhs.p_impl->is_equal_to(rhs.p_impl));
        };
        
        const point_type& operator*() const {
          return p_impl->get_point();
        };
      
    };
    
    
    
  public:
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     */
    explicit seq_trajectory_base(const std::string& aName = "") : named_object() { 
      setName(aName);
    };
    
    virtual ~seq_trajectory_base() { };
    
    /**
     * Returns the starting time-iterator along the path.
     * \return The starting time-iterator along the path.
     */
    virtual point_time_iterator begin_time_travel() const = 0;
    
    /**
     * Returns the end time-iterator along the path.
     * \return The end time-iterator along the path.
     */
    virtual point_time_iterator end_time_travel() const = 0;
    
    /**
     * Returns the starting fraction-iterator along the path.
     * \return The starting fraction-iterator along the path.
     */
    virtual point_fraction_iterator begin_fraction_travel() const = 0;
    
    /**
     * Returns the end fraction-iterator along the path.
     * \return The end fraction-iterator along the path.
     */
    virtual point_fraction_iterator end_fraction_travel() const = 0;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2440014,1,"seq_trajectory_base",named_object)
};



};

};

#endif









