/**
 * \file seq_trajectory_wrapper.hpp
 * 
 * This library provides a sequential trajectory-wrapper class template which makes an OOP-compatible trajectory class 
 * for a given sequential trajectory (in the generic sense).
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

#ifndef REAK_SEQ_TRAJECTORY_WRAPPER_HPP
#define REAK_SEQ_TRAJECTORY_WRAPPER_HPP

#include "base/defs.hpp"

#include "seq_trajectory_base.hpp"

#include "sequential_trajectory_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  


/**
 * This class wraps a generic trajectory class into an OOP interface.
 * It, itself, also models the generic SpatialTrajectoryConcept, so this wrapper can 
 * be used for both purposes.
 * \tparam SequentialTraj The trajectory type to be wrapped.
 */
template <typename SequentialTraj>
class seq_trajectory_wrapper : public seq_trajectory_base< typename sequential_trajectory_traits<SequentialTraj>::topology > {
  public:
    typedef seq_trajectory_base< typename sequential_trajectory_traits<SequentialTraj>::topology > base_type;
    typedef seq_trajectory_wrapper<SequentialTraj> self;
    
    typedef typename base_type::topology topology;
    typedef typename base_type::point_type point_type;
    
    BOOST_CONCEPT_ASSERT((SequentialTrajectoryConcept<SequentialTraj,topology>));
    
    typedef SequentialTraj wrapped_type;
    
  protected:
    SequentialTraj m_traj;
    
    typedef typename base_type::point_time_iterator_impl base_pt_time_iterator_impl;
    typedef typename sequential_trajectory_traits<SequentialTraj>::point_time_iterator gen_pt_time_iterator;
    
    struct point_time_iterator_impl : public base_pt_time_iterator_impl {
      
      gen_pt_time_iterator base_it;
      
      point_time_iterator_impl(gen_pt_time_iterator aBaseIt) : base_it(aBaseIt) { };
      
      virtual ~point_time_iterator_impl() { };
      
      virtual void move_by_distance(double d) {
        base_it += d;
      };
      
      virtual bool is_equal_to(const base_pt_time_iterator_impl* rhs) const {
        return (base_it == static_cast<const point_time_iterator_impl*>(rhs)->base_it);
      };
      
      virtual const point_type& get_point() const {
        return *base_it;
      };
      
      virtual base_pt_time_iterator_impl* clone() const {
        return new point_time_iterator_impl(base_it);
      };
      
    };
    
    typedef typename base_type::point_fraction_iterator_impl base_pt_frac_iterator_impl;
    typedef typename sequential_trajectory_traits<SequentialTraj>::point_fraction_iterator gen_pt_frac_iterator;
    
    struct point_fraction_iterator_impl : public base_pt_frac_iterator_impl {
      
      gen_pt_frac_iterator base_it;
      
      point_fraction_iterator_impl(gen_pt_frac_iterator aBaseIt) : base_it(aBaseIt) { };
      
      virtual ~point_fraction_iterator_impl() { };
      
      virtual void move_by_fraction(double f) {
        base_it += f;
      };
      
      virtual bool is_equal_to(const base_pt_frac_iterator_impl* rhs) const {
        return (base_it == static_cast<const point_fraction_iterator_impl*>(rhs)->base_it);
      };
      
      virtual const point_type& get_point() const {
        return *base_it;
      };
      
      virtual base_pt_frac_iterator_impl* clone() const {
        return new point_fraction_iterator_impl(base_it);
      };
      
    };
    
  public:
    
    typedef typename base_type::point_time_iterator point_time_iterator;
    typedef typename base_type::point_fraction_iterator point_fraction_iterator;
    
    wrapped_type& get_underlying_trajectory() { return m_traj; };
    const wrapped_type& get_underlying_trajectory() const { return m_traj; };
    
    wrapped_type& get_wrapped_object() { return m_traj; };
    const wrapped_type& get_wrapped_object() const { return m_traj; };
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     * \param aTraj The wrapped trajectory object to use.
     */
    explicit seq_trajectory_wrapper(const std::string& aName = "",
                                    const SequentialTraj& aTraj = SequentialTraj()) : 
                                    base_type(aName), m_traj(aTraj) { };
    
    virtual ~seq_trajectory_wrapper() { };
    
    /**
     * Returns the starting time-iterator along the trajectory.
     * \return The starting time-iterator along the trajectory.
     */
    virtual point_time_iterator begin_time_travel() const {
      return point_time_iterator(new point_time_iterator_impl(m_traj.begin_time_travel()));
    };
    
    /**
     * Returns the end time-iterator along the trajectory.
     * \return The end time-iterator along the trajectory.
     */
    virtual point_time_iterator end_time_travel() const {
      return point_time_iterator(new point_time_iterator_impl(m_traj.end_time_travel()));
    };
    
    /**
     * Returns the starting fraction-iterator along the trajectory.
     * \return The starting fraction-iterator along the trajectory.
     */
    virtual point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.begin_fraction_travel()));
    };
    
    /**
     * Returns the end fraction-iterator along the trajectory.
     * \return The end fraction-iterator along the trajectory.
     */
    virtual point_fraction_iterator end_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.end_fraction_travel()));
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_traj);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_traj);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440015,1,"seq_trajectory_wrapper",base_type)
  
};


};

};

#endif









