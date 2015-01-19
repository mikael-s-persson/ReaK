/**
 * \file seq_path_wrapper.hpp
 * 
 * This library provides a sequential path-wrapper class template which makes an OOP-compatible path class 
 * for a given sequential path (in the generic sense).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_SEQ_PATH_WRAPPER_HPP
#define REAK_SEQ_PATH_WRAPPER_HPP

#include <ReaK/core/base/defs.hpp>

#include "seq_path_base.hpp"
#include "sequential_path_concept.hpp"

#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {
  


/**
 * This class wraps a generic spatial path class into an OOP interface.
 * It, itself, also models the generic SpatialPathConcept, so this wrapper can 
 * be used for both purposes.
 * \tparam SequentialPath The path type to be wrapped.
 */
template <typename SequentialPath>
class seq_path_wrapper : public seq_path_base< typename sequential_path_traits<SequentialPath>::topology > {
  public:
    typedef seq_path_base< typename sequential_path_traits<SequentialPath>::topology > base_type;
    typedef seq_path_wrapper<SequentialPath> self;
    
    typedef typename base_type::topology topology;
    typedef typename base_type::point_type point_type;
    
    BOOST_CONCEPT_ASSERT((SequentialPathConcept<SequentialPath,topology>));
    
    typedef SequentialPath wrapped_type;
    
  protected:
    SequentialPath m_traj;
    
    typedef typename base_type::point_distance_iterator_impl base_pt_dist_iterator_impl;
    typedef typename sequential_path_traits<SequentialPath>::point_distance_iterator gen_pt_dist_iterator;
    
    struct point_distance_iterator_impl : public base_pt_dist_iterator_impl {
      
      gen_pt_dist_iterator base_it;
      
      point_distance_iterator_impl(gen_pt_dist_iterator aBaseIt) : base_it(aBaseIt) { };
      
      virtual ~point_distance_iterator_impl() { };
      
      virtual void move_by_distance(double d) {
        base_it += d;
      };
      
      virtual bool is_equal_to(const base_pt_dist_iterator_impl* rhs) const {
        return (base_it == static_cast<const point_distance_iterator_impl*>(rhs)->base_it);
      };
      
      virtual const point_type& get_point() const {
        return *base_it;
      };
      
      virtual base_pt_dist_iterator_impl* clone() const {
        return new point_distance_iterator_impl(base_it);
      };
      
    };
    
    typedef typename base_type::point_fraction_iterator_impl base_pt_frac_iterator_impl;
    typedef typename sequential_path_traits<SequentialPath>::point_fraction_iterator gen_pt_frac_iterator;
    
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
    
    typedef typename base_type::point_distance_iterator point_distance_iterator;
    typedef typename base_type::point_fraction_iterator point_fraction_iterator;
    
    wrapped_type& get_underlying_path() { return m_traj; };
    const wrapped_type& get_underlying_path() const { return m_traj; };
    
    wrapped_type& get_wrapped_object() { return m_traj; };
    const wrapped_type& get_wrapped_object() const { return m_traj; };
    
    /**
     * Constructs the trajectory from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aName The name for this object.
     * \param aTraj The wrapped path object to use.
     */
    explicit seq_path_wrapper(const std::string& aName = "",
                              const SequentialPath& aTraj = SequentialPath()) : 
                              base_type(aName),
                              m_traj(aTraj) { };
    
    virtual ~seq_path_wrapper() { };
    
    /**
     * Returns the starting distance-iterator along the path.
     * \return The starting distance-iterator along the path.
     */
    virtual point_distance_iterator begin_distance_travel() const {
      return point_distance_iterator(new point_distance_iterator_impl(m_traj.begin_distance_travel()));
    };
    
    /**
     * Returns the end distance-iterator along the path.
     * \return The end distance-iterator along the path.
     */
    virtual point_distance_iterator end_distance_travel() const {
      return point_distance_iterator(new point_distance_iterator_impl(m_traj.end_distance_travel()));
    };
    
    /**
     * Returns the starting fraction-iterator along the path.
     * \return The starting fraction-iterator along the path.
     */
    virtual point_fraction_iterator begin_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.begin_fraction_travel()));
    };
    
    /**
     * Returns the end fraction-iterator along the path.
     * \return The end fraction-iterator along the path.
     */
    virtual point_fraction_iterator end_fraction_travel() const {
      return point_fraction_iterator(new point_fraction_iterator_impl(m_traj.end_fraction_travel()));
    };
    
    virtual double travel_distance(const point_type& a, const point_type& b) const {
      return m_traj.travel_distance(a, b);
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2440012,1,"seq_path_wrapper",base_type)
  
};


};

};

#endif









