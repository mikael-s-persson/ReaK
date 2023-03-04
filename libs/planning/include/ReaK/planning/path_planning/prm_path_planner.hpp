/**
 * \file prm_path_planner.hpp
 *
 * This library defines a class to solve path planning problems using the
 * Probabilistic Road-map (PRM) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
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

#ifndef REAK_PRM_PATH_PLANNER_HPP
#define REAK_PRM_PATH_PLANNER_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "motion_planner_base.hpp"

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include "any_sbmp_reporter.hpp"

namespace ReaK::pp {

/**
 * This class solves path planning problems using the
 * Probabilistic Road-map (PRM) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger
 * configuration space.
 */
template <typename FreeSpaceType>
class prm_planner : public sample_based_planner<FreeSpaceType> {
 public:
  using base_type = sample_based_planner<FreeSpaceType>;
  using self = prm_planner<FreeSpaceType>;

  using space_type = FreeSpaceType;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));

  using point_type = topology_point_type_t<super_space_type>;
  using point_difference_type =
      topology_point_difference_type_t<super_space_type>;

  std::size_t get_motion_graph_kind() const override {
    return ASTAR_MOTION_GRAPH_KIND | DENSE_MOTION_GRAPH_KIND;
  }

  /**
   * This function is called to reset the internal state of the planner.
   */
  void reset_internal_state() override { base_type::reset_internal_state(); }

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
   * \param aSteerProgressTolerance The steer progress tolerance to be used by this planner when making connections.
   * \param aConnectionTolerance The connection tolerance to be used by this planner when making connections.
   * \param aSamplingRadius The sampling radius to be used by this planner when doing random walks.
   * \param aSpaceDimensionality The dimensionality of the space used by this planner.
   * \param aReporter The path-planning reporter to be used by this planner.
   */
  explicit prm_planner(const std::shared_ptr<space_type>& aWorld,
                       std::size_t aMaxVertexCount = 5000,
                       std::size_t aProgressInterval = 100,
                       std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH |
                                                         DVP_BF2_TREE_KNN,
                       double aSteerProgressTolerance = 0.1,
                       double aConnectionTolerance = 0.1,
                       double aSamplingRadius = 1.0,
                       std::size_t aSpaceDimensionality = 1,
                       const any_sbmp_reporter_chain<space_type>& aReporter =
                           any_sbmp_reporter_chain<space_type>())
      : base_type("prm_planner", aWorld, aMaxVertexCount, aProgressInterval,
                  aDataStructureFlags, 0, aSteerProgressTolerance,
                  aConnectionTolerance, aSamplingRadius, aSpaceDimensionality,
                  aReporter) {}

  prm_planner() : prm_planner(std::shared_ptr<space_type>()) {}

  ~prm_planner() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2460008, 1, "prm_planner", base_type)
};

}  // namespace ReaK::pp

#endif
