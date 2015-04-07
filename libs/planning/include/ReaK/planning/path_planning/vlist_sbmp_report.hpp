/**
 * \file vlist_sbmp_report.hpp
 *
 * This library defines a sampling-based motion/path planning reporter that prints a list of
 * vertex and their properties.
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

#ifndef REAK_VLIST_SBMP_REPORT_HPP
#define REAK_VLIST_SBMP_REPORT_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include <boost/concept_check.hpp>

#include <ReaK/topologies/spaces/subspace_concept.hpp>
#include <ReaK/topologies/interpolation/seq_trajectory_base.hpp>
#include <ReaK/topologies/interpolation/seq_path_base.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <fstream>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and prints all the vertices (the properties) of the motion graph into files.
 * \tparam VertexPrinter A functor type that can print, in a single line, all the properties of a vertex of the
 * motion-graph.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template < typename VertexPrinter, typename NextReporter = no_sbmp_report >
struct vlist_sbmp_report : public shared_object {
  typedef vlist_sbmp_report< VertexPrinter, NextReporter > self;

  /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
  NextReporter next_reporter;
  /// Holds the functor that can print, in a single line, all the properties of a vertex of the motion-graph.
  VertexPrinter print_to_stream;
  /// Holds the file-path where to output the reports.
  std::string file_path;

  /**
   * Parametrized constructor.
   * \param aFilePath The path where to create the output files.
   * \param aPrinter The functor that can print, in a single line, all the properties of a vertex of the motion-graph.
   * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
   */
  explicit vlist_sbmp_report( const std::string& aFilePath = "", VertexPrinter aPrinter = VertexPrinter(),
                              NextReporter aNextReporter = NextReporter() )
      : next_reporter( aNextReporter ), print_to_stream( aPrinter ), file_path( aFilePath ){};

  void reset_internal_state() { next_reporter.reset_internal_state(); };

  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   * \param free_space The C-free topology.
   * \param g The motion-graph.
   * \param pos The position-map to obtain positions of the motion-graph vertices.
   */
  template < typename FreeSpaceType, typename MotionGraph, typename PositionMap >
  void draw_motion_graph( const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos ) const {
    typedef typename boost::graph_traits< MotionGraph >::vertex_iterator VIter;

    std::stringstream ss;
    ss << std::setw( 6 ) << std::setfill( '0' ) << num_vertices( g );
    std::ofstream file_out( file_path + "vlist_" + ss.str() );

    VIter vi, vi_end;
    for( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
      print_to_stream( file_out, *vi, g );

    next_reporter.draw_motion_graph( free_space, g, pos );
  };


  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template < typename FreeSpaceType >
  void draw_solution(
    const FreeSpaceType& free_space,
    typename boost::enable_if< is_temporal_space< FreeSpaceType >,
    const shared_ptr< seq_trajectory_base< typename subspace_traits< FreeSpaceType >::super_space_type > >& >::type traj )
    const {
    next_reporter.draw_solution( free_space, traj );
  };

  /**
   * Draws the solution path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param p The solution path.
   */
  template < typename FreeSpaceType >
  void draw_solution(
    const FreeSpaceType& free_space,
    typename boost::disable_if< is_temporal_space< FreeSpaceType >,
    const shared_ptr< seq_path_base< typename subspace_traits< FreeSpaceType >::super_space_type > >& >::type p ) const {
    next_reporter.draw_solution( free_space, p );
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    shared_object::save( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( next_reporter ) & RK_SERIAL_SAVE_WITH_NAME( print_to_stream )
      & RK_SERIAL_SAVE_WITH_NAME( file_path );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    shared_object::load( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( next_reporter ) & RK_SERIAL_LOAD_WITH_NAME( print_to_stream )
      & RK_SERIAL_LOAD_WITH_NAME( file_path );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2460010, 1, "vlist_sbmp_report", shared_object )
};
};
};


#endif
