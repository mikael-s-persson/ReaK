/**
 * \file planning_space_options.hpp
 * 
 * This library defines the options available when creating the configuration-space in which a path-planner is run.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PLANNING_SPACE_OPTIONS_HPP
#define REAK_PLANNING_SPACE_OPTIONS_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"


namespace ReaK {
  
namespace pp {


/// This mask indicates the planning-space's order for the planning problem.
const std::size_t PLANNING_SPACE_ORDER_MASK = 0x07;

/// This indicates the planning-space's order is 0 (position only).
const std::size_t PLANNING_SPACE_ORDER_ZERO  = 0;
/// This indicates the planning-space's order is 1 (position + velocity).
const std::size_t PLANNING_SPACE_ORDER_ONE   = 1;
/// This indicates the planning-space's order is 2 (position + velocity + acceleration).
const std::size_t PLANNING_SPACE_ORDER_TWO   = 2;
/// This indicates the planning-space's order is 3 (position + velocity + acceleration + jerk).
const std::size_t PLANNING_SPACE_ORDER_THREE = 3;

/// This mask indicates the planning-space's interpolation method for the planning problem.
const std::size_t PLANNING_SPACE_INTERPOLATOR_MASK = 0x1F << 3;

/// This indicates the planning-space uses linear interpolation.
const std::size_t PLANNING_SPACE_LINEAR_INTERP  = 0 << 3;
/// This indicates the planning-space uses cubic interpolation.
const std::size_t PLANNING_SPACE_CUBIC_INTERP   = 1 << 3;
/// This indicates the planning-space uses quintic interpolation.
const std::size_t PLANNING_SPACE_QUINTIC_INTERP = 2 << 3;
/// This indicates the planning-space uses sustained velocity pulse interpolation (acceleration-bounded).
const std::size_t PLANNING_SPACE_SVP_INTERP     = 3 << 3;
/// This indicates the planning-space uses sustained acceleration pulse interpolation (jerk-bounded).
const std::size_t PLANNING_SPACE_SAP_INTERP     = 4 << 3;


/// This flag indicates that the planning is done in a temporal-space (space-time, or dynamic).
const std::size_t PLAN_IN_TEMPORAL_SPACE = 0x01 << 8;

/// This flag indicates that the spatial components are normalized to their reach-time distances.
const std::size_t PLAN_IN_RATE_LIMITED_SPACE = 0x02 << 8;



class planning_space_options : public shared_object {
  public:
    std::size_t space_options;
    std::size_t output_space_options;
    double min_travel;
    double max_travel;
    
    planning_space_options() : 
      space_options(0), output_space_options(0),      // <- 0-th order, linear-interp, non-temporal, not rate-limited.
      min_travel(0.0), max_travel(1.0) { };
    
    
    std::size_t get_space_order() const {
      return space_options & PLANNING_SPACE_ORDER_MASK;
    };
    void set_space_order(std::size_t aOrder) {
      space_options &= ~PLANNING_SPACE_ORDER_MASK;
      space_options |= aOrder & PLANNING_SPACE_ORDER_MASK;
    };
    
    std::size_t get_interp_id() const {
      return ( space_options & PLANNING_SPACE_INTERPOLATOR_MASK ) >> 3;
    };
    void set_interp_id(std::size_t aID) {
      space_options &= ~PLANNING_SPACE_INTERPOLATOR_MASK;
      space_options |= ( aID << 3 ) & PLANNING_SPACE_INTERPOLATOR_MASK;
    };
    
    bool is_temporal_space() const {
      return space_options & PLAN_IN_TEMPORAL_SPACE;
    };
    void set_temporal_space(bool aIsTemporal = true) {
      if(aIsTemporal)
        space_options |= PLAN_IN_TEMPORAL_SPACE;
      else
        space_options &= ~PLAN_IN_TEMPORAL_SPACE;
    };
    
    bool is_rate_limited() const {
      return space_options & PLAN_IN_RATE_LIMITED_SPACE;
    };
    void set_rate_limited(bool aIsRateLimited = true) {
      if(aIsRateLimited)
        space_options |= PLAN_IN_RATE_LIMITED_SPACE;
      else
        space_options &= ~PLAN_IN_RATE_LIMITED_SPACE;
    };
    
    
    std::size_t get_output_space_order() const {
      return output_space_options & PLANNING_SPACE_ORDER_MASK;
    };
    void set_output_space_order(std::size_t aOrder) {
      output_space_options &= ~PLANNING_SPACE_ORDER_MASK;
      output_space_options |= aOrder & PLANNING_SPACE_ORDER_MASK;
    };
    
    std::size_t get_output_interp_id() const {
      return ( output_space_options & PLANNING_SPACE_INTERPOLATOR_MASK ) >> 3;
    };
    void set_output_interp_id(std::size_t aID) {
      output_space_options &= ~PLANNING_SPACE_INTERPOLATOR_MASK;
      output_space_options |= ( aID << 3 ) & PLANNING_SPACE_INTERPOLATOR_MASK;
    };
    
    bool is_temporal_output_space() const {
      return output_space_options & PLAN_IN_TEMPORAL_SPACE;
    };
    void set_temporal_output_space(bool aIsTemporal = true) {
      if(aIsTemporal)
        output_space_options |= PLAN_IN_TEMPORAL_SPACE;
      else
        output_space_options &= ~PLAN_IN_TEMPORAL_SPACE;
    };
    
    bool is_output_rate_limited() const {
      return output_space_options & PLAN_IN_RATE_LIMITED_SPACE;
    };
    void set_output_rate_limited(bool aIsRateLimited = true) {
      if(aIsRateLimited)
        output_space_options |= PLAN_IN_RATE_LIMITED_SPACE;
      else
        output_space_options &= ~PLAN_IN_RATE_LIMITED_SPACE;
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(space_options)
        & RK_SERIAL_SAVE_WITH_NAME(output_space_options)
        & RK_SERIAL_SAVE_WITH_NAME(min_travel)
        & RK_SERIAL_SAVE_WITH_NAME(max_travel);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(space_options)
        & RK_SERIAL_LOAD_WITH_NAME(output_space_options)
        & RK_SERIAL_LOAD_WITH_NAME(min_travel)
        & RK_SERIAL_LOAD_WITH_NAME(max_travel);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(planning_space_options,0xC246001A,1,"planning_space_options",shared_object)
    
    
};




};

};

#endif

