/**
 * \file pp_tester_base.hpp
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

#ifndef REAK_PP_TESTER_BASE_HPP
#define REAK_PP_TESTER_BASE_HPP


#include "metric_space_concept.hpp"
#include "random_sampler_concept.hpp"

namespace ReaK {
  
namespace pp {

/**
 * This class is
 */
template <typename FreeSpaceType>
class pp_tester_base {
  public:
    typedef FreeSpaceType space_type;
    typedef topology_traits< FreeSpaceType >::point_type point_type;
    typedef topology_traits< FreeSpaceType >::point_difference_type point_difference_type;
    
    typedef metric_space_traits< FreeSpaceType >::distance_metric_type distance_metric_type;
    typedef point_distribution_traits< FreeSpaceType >::random_sampler_type random_sampler_type;
    

  protected:
    
    shared_ptr< space_type > m_space;
    
  public:
    
    virtual void draw_motion_graph() = 0;
    virtual void draw_best_solution() = 0;
    
    virtual void get_best_solution(std::list< point_type >& path) = 0;

    //PRM Visitor concepts:
    
    virtual bool keep_going() const {
      return true;
    };
    
    virtual double heuristic(const point_type& p_u) {
      return m_space->bird_fly_to_goal(p_u);
    };

    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    pp_tester_base(const shared_ptr< space_type >& aWorld) :
                   m_space(aWorld) { };

    //This is a RAII class, so the destructor is meaningless. Hurray for RAII!!
    ~pp_tester_base() { };
    
};

};

};

#endif

