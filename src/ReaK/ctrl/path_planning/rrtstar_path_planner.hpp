/**
 * \file rrtstar_path_planner.hpp
 * 
 * This library defines a class to solve path planning problems using the 
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_RRTSTAR_PATH_PLANNER_HPP
#define REAK_RRTSTAR_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"

#include "metric_space_concept.hpp"
#include "any_sbmp_reporter.hpp"


namespace ReaK {
  
namespace pp {


/**
 * This class solves path planning problems using the 
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename FreeSpaceType>
class rrtstar_planner : public sample_based_planner<FreeSpaceType> {
  public:
    typedef sample_based_planner<FreeSpaceType> base_type;
    typedef rrtstar_planner<FreeSpaceType> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    
    template <typename RRTStarFactory>
    void solve_planning_query_impl(planning_query<FreeSpaceType>& aQuery);
    
  public:
    
    virtual std::size_t get_motion_graph_kind() const { 
      if((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) == UNIDIRECTIONAL_PLANNING)
        return OPTIMAL_MOTION_GRAPH_KIND;
      else 
        return BIDIR_OPTIMAL_MOTION_GRAPH_KIND;
    };
    
    /**
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \param aQuery The query object that defines as input the parameters of the query, 
     *               and as output, the recorded solutions.
     */
    virtual void solve_planning_query(planning_query<FreeSpaceType>& aQuery);
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     * \param aDataStructureFlags An integer flags representing the kind of motion graph data-structure to use in the 
     *                            planning algorithm. Can be ADJ_LIST_MOTION_GRAPH or DVP_ADJ_LIST_MOTION_GRAPH.
     *                            Any combination of those two and of KNN method flags to use for nearest
     *                            neighbor queries in the graph. KNN method flags can be LINEAR_SEARCH_KNN, 
     *                            DVP_BF2_TREE_KNN, DVP_BF4_TREE_KNN, DVP_COB2_TREE_KNN, or DVP_COB4_TREE_KNN.
     *                            See path_planner_options.hpp documentation.
     * \param aPlanningMethodFlags The integer flags that identify various options to use with this planner.
     *                             The options available include EAGER_COLLISION_CHECKING or LAZY_COLLISION_CHECKING, 
     *                             NOMINAL_PLANNER_ONLY or any combination of PLAN_WITH_VORONOI_PULL, 
     *                             PLAN_WITH_NARROW_PASSAGE_PUSH and PLAN_WITH_ANYTIME_HEURISTIC, UNIDIRECTIONAL_PLANNING 
     *                             or BIDIRECTIONAL_PLANNING, and USE_BRANCH_AND_BOUND_PRUNING_FLAG. 
     * \param aSteerProgressTolerance The steer progress tolerance to be used by this planner when making connections.
     * \param aConnectionTolerance The connection tolerance to be used by this planner when making connections.
     * \param aSpaceDimensionality The dimensionality of the space used by this planner.
     * \param aReporter The path-planning reporter to be used by this planner.
     */
    rrtstar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                    std::size_t aMaxVertexCount = 5000, 
                    std::size_t aProgressInterval = 100,
                    std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                    std::size_t aPlanningMethodFlags = BIDIRECTIONAL_PLANNING,
                    double aSteerProgressTolerance = 0.1,
                    double aConnectionTolerance = 0.1,
                    std::size_t aSpaceDimensionality = 1,
                    const any_sbmp_reporter_chain<space_type>& aReporter = any_sbmp_reporter_chain<space_type>()) :
                    base_type("rrtstar_planner", aWorld, aMaxVertexCount, aProgressInterval,
                              aDataStructureFlags, aPlanningMethodFlags,
                              aSteerProgressTolerance, aConnectionTolerance, 1.0, 
                              aSpaceDimensionality, aReporter) { };
    
    virtual ~rrtstar_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460009,1,"rrtstar_planner",base_type)
};

};

};



#endif

