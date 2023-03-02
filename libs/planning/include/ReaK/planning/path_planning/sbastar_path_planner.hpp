/**
 * \file sbastar_path_planner.hpp
 *
 * This library defines a class to solve path planning problems using the
 * Sampling-based A* algorithm (or one of its variants). Given a C_free (configuration space
 * restricted to non-colliding points) and a result reporting policy, this class
 * will probabilistically construct a motion-graph that will connect a starting point
 * and a goal point with a path through C-free that is as close as possible to the
 * optimal path in terms of distance. The planner uses a selectable variant of the
 * Sampling-based A* (SBA*) algorithm, including the basic version, the SBA*-RRT*
 * alternating algorithm, and the Anytime SBA* algorithm. In all cases, collision
 * checking and connectivity can be either full or lazy (and pruned) to either construct
 * a full-connectivity graph containing only collision-free edges, or a single-query motion-tree
 * that includes only optimal edges (whose collisions are checked lazily).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_SBASTAR_PATH_PLANNER_HPP
#define REAK_SBASTAR_PATH_PLANNER_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "motion_planner_base.hpp"

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include "any_sbmp_reporter.hpp"

namespace ReaK::pp {

/**
 * This class solves path planning problems using the
 * Sampling-based A* algorithm (or one of its variants). Given a C_free (configuration space
 * restricted to non-colliding points) and a result reporting policy, this class
 * will probabilistically construct a motion-graph that will connect a starting point
 * and a goal point with a path through C-free that is as close as possible to the
 * optimal path in terms of distance. The planner uses a selectable variant of the
 * Sampling-based A* (SBA*) algorithm, including the basic version, the SBA*-RRT*
 * alternating algorithm, and the Anytime SBA* algorithm. In all cases, collision
 * checking and connectivity can be either full or lazy (and pruned) to either construct
 * a full-connectivity graph containing only collision-free edges, or a single-query motion-tree
 * that includes only optimal edges (whose collisions are checked lazily).
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger
 * configuration space.
 */
template <typename FreeSpaceType>
class sbastar_planner : public sample_based_planner<FreeSpaceType> {
 public:
  using base_type = sample_based_planner<FreeSpaceType>;
  using self = sbastar_planner<FreeSpaceType>;

  using space_type = FreeSpaceType;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));

  using point_type = topology_point_type_t<super_space_type>;
  using point_difference_type =
      topology_point_difference_type_t<super_space_type>;

 protected:
  double m_init_dens_threshold;
  double m_init_relaxation;
  double m_SA_init_temperature;

  template <typename SBAStarFactory>
  void solve_planning_query_impl(planning_query<FreeSpaceType>& aQuery);

 public:
  std::size_t get_motion_graph_kind() const override {
    if ((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) ==
        UNIDIRECTIONAL_PLANNING) {
      return ASTAR_MOTION_GRAPH_KIND | RECURSIVE_DENSE_MOTION_GRAPH_KIND;
    }
    return BIDIR_ASTAR_MOTION_GRAPH_KIND | RECURSIVE_DENSE_MOTION_GRAPH_KIND;
  }

  /**
   * This function computes a valid path in the C-free. If it cannot
   * achieve a valid path, an exception will be thrown. This algorithmic
   * path solver class is such that any settings that ought to be set for the
   * path planning algorithm should be set before calling this function, otherwise
   * the function is likely to fail.
   * \param aQuery The query object that defines as input the parameters of the query,
   *               and as output, the recorded solutions.
   */
  void solve_planning_query(planning_query<FreeSpaceType>& aQuery) override;

  /**
   * Returns the initial density-value threshold used by this planner.
   * \return The initial density-value threshold used by this planner.
   */
  double get_initial_density_threshold() const { return m_init_dens_threshold; }
  /**
   * Sets the initial density-value threshold to be used by this planner.
   * \param aInitialThreshold The initial density-value threshold to be used by this planner.
   */
  void set_initial_density_threshold(double aInitialThreshold) {
    m_init_dens_threshold = aInitialThreshold;
  }

  /**
   * Returns the initial relaxation factor used by this planner.
   * \return The initial relaxation factor used by this planner.
   */
  double get_initial_relaxation() const { return m_init_relaxation; }
  /**
   * Sets the initial relaxation factor to be used by this planner.
   * \param aInitialRelaxation The initial relaxation factor to be used by this planner.
   */
  void set_initial_relaxation(double aInitialRelaxation) {
    m_init_relaxation = aInitialRelaxation;
  }

  /**
   * Returns the initial Simulated Annealing temperature use by this planner, if the
   * added-bias is set to an exploratory bias (e.g., PLAN_WITH_VORONOI_PULL). If negative,
   * then simulated annealing is not used, and the exploratory bias (if any) is applied
   * only when SBA* seaching stalls (isn't progressing anymore).
   * \return The initial Simulated Annealing temperature used by this planner.
   */
  double get_initial_SA_temperature() const { return m_SA_init_temperature; }
  /**
   * Sets the initial Simulated Annealing temperature use by this planner, if the
   * added-bias is set to an exploratory bias (e.g., PLAN_WITH_VORONOI_PULL). If negative,
   * then simulated annealing is not used, and the exploratory bias (if any) is applied
   * only when SBA* seaching stalls (isn't progressing anymore).
   * \param aInitialSATemperature The initial Simulated Annealing temperature to be used by this planner.
   */
  void set_initial_SA_temperature(double aInitialSATemperature) {
    m_SA_init_temperature = aInitialSATemperature;
  }

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
   * \param aSamplingRadius The sampling radius to be used by this planner when doing random walks.
   * \param aSpaceDimensionality The dimensionality of the space used by this planner.
   * \param aReporter The path-planning reporter to be used by this planner.
   */
  explicit sbastar_planner(
      const std::shared_ptr<space_type>& aWorld,
      std::size_t aMaxVertexCount = 5000, std::size_t aProgressInterval = 100,
      std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH |
                                        DVP_BF2_TREE_KNN,
      std::size_t aPlanningMethodFlags = LAZY_COLLISION_CHECKING |
                                         NOMINAL_PLANNER_ONLY,
      double aSteerProgressTolerance = 0.1, double aConnectionTolerance = 0.1,
      double aSamplingRadius = 1.0, std::size_t aSpaceDimensionality = 1,
      const any_sbmp_reporter_chain<space_type>& aReporter =
          any_sbmp_reporter_chain<space_type>())
      : base_type("sbastar_planner", aWorld, aMaxVertexCount, aProgressInterval,
                  aDataStructureFlags, aPlanningMethodFlags,
                  aSteerProgressTolerance, aConnectionTolerance,
                  aSamplingRadius, aSpaceDimensionality, aReporter),
        m_init_dens_threshold(0.8),
        m_init_relaxation(0.0),
        m_SA_init_temperature(-1.0) {}

  sbastar_planner() : sbastar_planner(std::shared_ptr<space_type>()) {}

  ~sbastar_planner() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_init_dens_threshold) &
        RK_SERIAL_SAVE_WITH_NAME(m_init_relaxation) &
        RK_SERIAL_SAVE_WITH_NAME(m_SA_init_temperature);
  }

  void load(serialization::iarchive& A, unsigned int) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_init_dens_threshold) &
        RK_SERIAL_LOAD_WITH_NAME(m_init_relaxation) &
        RK_SERIAL_LOAD_WITH_NAME(m_SA_init_temperature);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC246000C, 1, "sbastar_planner", base_type)
};

}  // namespace ReaK::pp

#endif
