/**
 * \file sbmp_reporter_concept.hpp
 * 
 * This library defines the traits and concepts that pertain to what can be considered 
 * a sampling-based motion planning reporter. Such a reporter report on the progress of 
 * a SBMP method by drawing the complete motion graph and/or a trajectory obtained from
 * the motion planner.
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

#ifndef REAK_SBMP_REPORTER_CONCEPT_HPP
#define REAK_SBMP_REPORTER_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <boost/concept_check.hpp>

#include <ReaK/topologies/interpolation/seq_trajectory_base.hpp>
#include <ReaK/topologies/interpolation/seq_path_base.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {

  
  

/**
 * This concept defines the requirements to fulfill in order to model 
 * a Sampling-based Path-planner (SBPP) reporter as used in ReaK::pp.
 * 
 * Valid expressions:
 * 
 * where:
 * 
 * SBPPReporter reporter;
 * 
 * FreeSpaceType free_space;
 * 
 * MotionGraph g;
 * 
 * PositionMap pos_map;
 * 
 * shared_ptr< path_base< super_space_type > > path;
 * 
 * 
 * reporter.draw_motion_graph(free_space, g, pos_map);  The reporter can be asked to draw the current motion-graph.
 * 
 * reporter.draw_solution(free_space, path);  The reporter can be asked to draw a solution path.
 * 
 * \tparam SBMPReporter The reporter type to be checked for this concept.
 * \tparam FreeSpaceType The topology type that represents the C-free sub-space, should model SubSpaceConcept.
 * \tparam MotionGraph The motion-graph type that represents the motion samples.
 * \tparam PositionMap The property-map type to fetch positions associated to vertices of the motion graph.
 */
template <typename SBMPReporter, 
          typename FreeSpaceType, 
          typename MotionGraph, 
          typename PositionMap>
struct SBPPReporterConcept {
  SBMPReporter reporter;
  FreeSpaceType free_space;
  MotionGraph g;
  PositionMap pos_map;
  shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > > path;
  
  BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
  
  BOOST_CONCEPT_USAGE(SBPPReporterConcept) 
  {
    reporter.reset_internal_state();
    reporter.draw_motion_graph(free_space, g, pos_map);
    reporter.draw_solution(free_space, path);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model 
 * a Sampling-based Motion-planner (SBMP) reporter as used in ReaK::pp.
 * 
 * Valid expressions:
 * 
 * where:
 * 
 * SBMPReporter reporter;
 * 
 * FreeSpaceType free_space;
 * 
 * MotionGraph g;
 * 
 * PositionMap pos_map;
 * 
 * shared_ptr< trajectory_base< super_space_type > > traj;
 * 
 * 
 * reporter.draw_motion_graph(free_space, g, pos_map);  The reporter can be asked to draw the current motion-graph.
 * 
 * reporter.draw_solution(free_space, traj);  The reporter can be asked to draw a solution trajectory.
 * 
 * \tparam SBMPReporter The reporter type to be checked for this concept.
 * \tparam FreeSpaceType The topology type that represents the C-free sub-space, should model SubSpaceConcept.
 * \tparam MotionGraph The motion-graph type that represents the motion samples.
 * \tparam PositionMap The property-map type to fetch positions associated to vertices of the motion graph.
 */
template <typename SBMPReporter, 
          typename FreeSpaceType, 
          typename MotionGraph, 
          typename PositionMap>
struct SBMPReporterConcept {
  SBMPReporter reporter;
  FreeSpaceType free_space;
  MotionGraph g;
  PositionMap pos_map;
  shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > > traj;
  
  BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
  
  BOOST_CONCEPT_USAGE(SBMPReporterConcept) 
  {
    reporter.reset_internal_state();
    reporter.draw_motion_graph(free_space, g, pos_map);
    reporter.draw_solution(free_space, traj);
  };
  
};



};

};


#endif


