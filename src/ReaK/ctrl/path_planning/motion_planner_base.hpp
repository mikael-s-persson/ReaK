/**
 * \file motion_planner_base.hpp
 * 
 * This library defines a class
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

#ifndef REAK_MOTION_PLANNER_BASE_HPP
#define REAK_MOTION_PLANNER_BASE_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "metric_space_concept.hpp"
#include "random_sampler_concept.hpp"
#include "subspace_concept.hpp"

namespace ReaK {
  
namespace pp {

/**
 * This class is the basic OOP interface for a motion planner.
 */
template <typename FreeSpaceType>
class motion_planner_base : public named_object {
  public:
    typedef motion_planner_base<FreeSpaceType> self;
    typedef FreeSpaceType space_type;
    typedef subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    

  protected:
    
    shared_ptr< space_type > m_space;
    
  public:
    
    /**
     * This function computes a valid trajectory in the C-free. If it cannot 
     * achieve a valid trajectory, an exception will be thrown. This algorithmic
     * motion solver class is such that any settings that ought to be set for the 
     * motion planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The trajectory object that can be used to map out the motion in time.
     */
    virtual shared_ptr< trajectory_base< super_space_type > > solve_motion() = 0;
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    motion_planner_base(const std::string& aName,
                        const shared_ptr< space_type >& aWorld) :
                        named_object()
                        m_space(aWorld) { setName(aName); };
    
    virtual ~motion_planner_base() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_space);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460000,1,"motion_planner_base",named_object)
};


/**
 * This class is the basic OOP interface for a path planner.
 */
template <typename FreeSpaceType>
class path_planner_base : public named_object {
  public:
    typedef path_planner_base<FreeSpaceType> self;
    typedef FreeSpaceType space_type;
    typedef subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    

  protected:
    
    shared_ptr< space_type > m_space;
    
  public:
    
    /**
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The path object that can be used to map out the path.
     */
    virtual shared_ptr< path_base< super_space_type > > solve_path() = 0;
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    path_planner_base(const std::string& aName,
                        const shared_ptr< space_type >& aWorld) :
                        named_object()
                        m_space(aWorld) { setName(aName); };
    
    virtual ~path_planner_base() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_space);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460001,1,"path_planner_base",named_object)
};


/**
 * This class is the basic OOP interface for a sampling-based motion planner.
 */
template <typename BaseType>
class sample_based_planner : public BaseType {
  protected:
    typedef BaseType base_type;
    
    std::size_t m_max_vertex_count;
    std::size_t m_progress_interval;
    
    
  public:
    
    /**
     * Returns the maximum number of samples to generate during the motion planning.
     * \return the maximum number of samples to generate during the motion planning.
     */
    std::size_t get_max_vertex_count() const { return m_max_vertex_count; };
    
    /**
     * Sets the maximum number of samples to generate during the motion planning.
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     */
    virtual void set_max_vertex_count(std::size_t aMaxVertexCount) { 
      m_max_vertex_count = aMaxVertexCount;
    };
    
    /**
     * Returns the number of new samples between each "progress report".
     * \return the number of new samples between each "progress report".
     */
    std::size_t get_progress_interval() const { return m_progress_interval; };
    
    /**
     * Sets the number of new samples between each "progress report".
     * \param aProgressInterval The number of new samples between each "progress report".
     */
    virtual void set_progress_interval(std::size_t aProgressInterval) { 
      m_progress_interval = aProgressInterval;
    };
    
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     */
    sample_based_planner(const std::string& aName,
                         const shared_ptr< space_type >& aWorld, 
                         std::size_t aMaxVertexCount = 0, 
                         std::size_t aProgressInterval = 0) :
                         base_type(aName,aWorld), 
                         m_max_vertex_count(aMaxVertexCount), 
                         m_progress_interval(aProgressInterval) { };
    
    virtual ~sample_based_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_max_vertex_count)
        & RK_SERIAL_SAVE_WITH_NAME(m_progress_interval);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_max_vertex_count)
        & RK_SERIAL_LOAD_WITH_NAME(m_progress_interval);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460002,1,"sample_based_planner",base_type)
};

};

};

#endif

