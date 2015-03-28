/**
 * \file frame_tracer_coin3d.hpp
 *
 * This library defines a sampling-based motion/path planning reporter for tracing out the path of a frame
 * linked to a DK kinematic model.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_FRAME_TRACER_COIN3D_HPP
#define REAK_FRAME_TRACER_COIN3D_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include <boost/concept_check.hpp>

#include <ReaK/topologies/interpolation/trajectory_base.hpp>
#include "basic_sbmp_reporters.hpp"
#include <boost/graph/graph_concepts.hpp>

#include <ReaK/topologies/spaces/direct_kinematics_topomap.hpp>
#include <ReaK/topologies/spaces/topological_map_concepts.hpp>

#include <ReaK/topologies/spaces/proxy_model_updater.hpp>

#include <ReaK/mbd/coin3D/frame_tracer_coin3d_impl.hpp>

class SoSeparator; // forward-declare


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept)
 * and uses the Coin3D library to trace out the position of a number of 3D poses linked to a model
 * linked with the given DirectKinMapper.
 * \tparam JointStateSpace A joint-state space type.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template < typename JointStateSpace, typename NextReporter = no_sbmp_report >
class frame_tracer_3D : public shared_object {
public:
  typedef frame_tracer_3D< JointStateSpace, NextReporter > self;

  /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
  NextReporter next_reporter;

  typedef proxy_model_applicator< JointStateSpace > applicator_type;

protected:
  /// A shared-pointer to the joint configuration applicator (e.g., direct-kinematics calculator).
  shared_ptr< applicator_type > mdl_applicator;
  /// Holds the interval-size between output points of the solution trajectory/path.
  double interval_size;

  /// Holds the list of frames (poses) to trace out.
  std::vector< shared_ptr< pose_3D< double > > > traced_frames;
  /// Holds the list of traces that represent the current motion-graph.
  std::vector< geom::tracer_coin3d_impl > motion_graph_traces;
  /// Holds the list of traces that represent each solution recorded.
  mutable std::map< double, std::vector< geom::tracer_coin3d_impl > > solution_traces;

  bool enable_mg_traces;

public:
  /**
   * Parametrized constructor.
   * \param aMdlApplicator The shared-pointer to the joint configuration applicator (e.g., direct-kinematics
   * calculator).
   * \param aIntervalSize The interval-size between output points of the solution trajectory/path.
   * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
   */
  explicit frame_tracer_3D( const shared_ptr< applicator_type >& aMdlApplicator = shared_ptr< applicator_type >(),
                            double aIntervalSize = 0.1, bool aTraceMotionGraph = true,
                            NextReporter aNextReporter = NextReporter() )
      : next_reporter( aNextReporter ), mdl_applicator( aMdlApplicator ), interval_size( aIntervalSize ),
        enable_mg_traces( aTraceMotionGraph ){};

  /**
   * This function adds a frame to the list of traced frames, i.e., the frames whose
   * position will be traced out in the Coin3D scene-graph.
   * \param aPose The frame to add to the list of traced frames.
   * \return A reference to this tracer.
   */
  self& add_traced_frame( const shared_ptr< pose_3D< double > >& aPose ) {
    traced_frames.push_back( aPose );
    motion_graph_traces.push_back( geom::tracer_coin3d_impl( false ) );
    return *this;
  };

  /**
   * This function returns the motion-graph trace associated to the given frame (pose).
   * \param aPose The frame of which the motion-graph trace is sought.
   * \return A const-reference to the motion-graph trace as a Coin3d scene-graph.
   */
  const geom::tracer_coin3d_impl& get_motion_graph_tracer( const shared_ptr< pose_3D< double > >& aPose ) const {
    std::vector< shared_ptr< pose_3D< double > > >::const_iterator it
      = std::find( traced_frames.begin(), traced_frames.end(), aPose );
    if( it != traced_frames.end() ) {
      return motion_graph_traces[it - traced_frames.begin()];
    } else
      throw std::range_error( "The given pose is not being traced by this frame-tracer!" );
  };

  /**
   * This function returns the number of solutions registers by this reporter.
   * \return The number of solutions registers by this reporter.
   */
  std::size_t get_solution_count() const { return solution_traces.size(); };

  /**
   * This function returns the cost of the best solution registered by this reporter.
   * \return The cost of the best solution registered by this reporter.
   */
  double get_best_solution_value() const {
    if( solution_traces.empty() )
      return std::numeric_limits< double >::infinity();
    return solution_traces.begin()->first;
  };

  /**
   * This function returns the solution trace associated to the given frame (pose).
   * \param aPose The frame of which the solution trace is sought.
   * \param aSolutionId The index of the solution trace, 0 is the best solution.
   * \return A const-reference to the solution trace as a Coin3d scene-graph.
   */
  const geom::tracer_coin3d_impl& get_solution_tracer( const shared_ptr< pose_3D< double > >& aPose,
                                                       std::size_t aSolutionId = 0 ) const {
    if( aSolutionId >= solution_traces.size() )
      aSolutionId = 0;
    std::vector< shared_ptr< pose_3D< double > > >::const_iterator it
      = std::find( traced_frames.begin(), traced_frames.end(), aPose );
    if( it != traced_frames.end() ) {
      std::map< double, std::vector< geom::tracer_coin3d_impl > >::const_iterator itm = solution_traces.begin();
      std::advance( itm, aSolutionId );
      return itm->second[it - traced_frames.begin()];
    } else
      throw std::range_error( "The given pose is not being traced by this frame-tracer!" );
  };

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
  typename boost::disable_if< is_steerable_space< FreeSpaceType >, void >::type
    draw_motion_graph( const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos ) const {

    if( !enable_mg_traces ) {
      next_reporter.draw_motion_graph( free_space, g, pos );
      return;
    };

    typedef typename boost::graph_traits< MotionGraph >::vertex_iterator VIter;
    typedef typename boost::graph_traits< MotionGraph >::out_edge_iterator EIter;
    typedef typename topology_traits< FreeSpaceType >::point_type PointType;

    VIter vi, vi_end;
    for( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi ) {
      EIter ei, ei_end;
      PointType p_u = get( pos, *vi );
      for( boost::tie( ei, ei_end ) = out_edges( *vi, g ); ei != ei_end; ++ei ) {
        mdl_applicator->apply_to_model( p_u, free_space.get_super_space() );
        for( std::size_t i = 0; i < traced_frames.size(); ++i )
          motion_graph_traces[i].begin_edge( traced_frames[i]->getGlobalPose().Position );
        PointType p_v = get( pos, target( *ei, g ) );
        for( double j = 0.1; j <= 1.01; j += 0.1 ) {
          PointType p_new = free_space.get_super_space().move_position_toward( p_u, j, p_v );
          mdl_applicator->apply_to_model( p_new, free_space.get_super_space() );
          for( std::size_t i = 0; i < traced_frames.size(); ++i )
            motion_graph_traces[i].add_point( traced_frames[i]->getGlobalPose().Position );
        };
        for( std::size_t i = 0; i < traced_frames.size(); ++i )
          motion_graph_traces[i].end_edge();
      };
    };

    next_reporter.draw_motion_graph( free_space, g, pos );
  };


  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam SteerRecMap The property-map type that can map motion-graph edge descriptors into steer-records.
   * \param free_space The C-free topology.
   * \param g The motion-graph.
   * \param steer_rec The steer-record-map to obtain steer-records of the motion-graph edges.
   */
  template < typename FreeSpaceType, typename MotionGraph, typename SteerRecMap >
  typename boost::enable_if< is_steerable_space< FreeSpaceType >, void >::type
    draw_motion_graph( const FreeSpaceType& free_space, const MotionGraph& g, SteerRecMap steer_rec ) const {

    if( !enable_mg_traces ) {
      next_reporter.draw_motion_graph( free_space, g, steer_rec );
      return;
    };

    typedef typename boost::graph_traits< MotionGraph >::vertex_iterator VIter;
    typedef typename boost::graph_traits< MotionGraph >::out_edge_iterator EIter;
    typedef typename topology_traits< FreeSpaceType >::point_type PointType;
    typedef typename boost::property_traits< SteerRecMap >::value_type SteerRecType;
    typedef typename SteerRecType::point_fraction_iterator SteerIter;

    VIter vi, vi_end;
    for( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi ) {
      EIter ei, ei_end;
      for( boost::tie( ei, ei_end ) = out_edges( *vi, g ); ei != ei_end; ++ei ) {
        const SteerRecType& st_rec = get( steer_rec, *ei );
        SteerIter it = st_rec.begin_fraction_travel();
        mdl_applicator->apply_to_model( PointType( *it ), free_space.get_super_space() );
        for( std::size_t i = 0; i < traced_frames.size(); ++i )
          motion_graph_traces[i].begin_edge( traced_frames[i]->getGlobalPose().Position );

        for( ; it != st_rec.end_fraction_travel(); it += 0.1 ) {
          mdl_applicator->apply_to_model( PointType( *it ), free_space.get_super_space() );
          for( std::size_t i = 0; i < traced_frames.size(); ++i )
            motion_graph_traces[i].add_point( traced_frames[i]->getGlobalPose().Position );
        };

        for( std::size_t i = 0; i < traced_frames.size(); ++i )
          motion_graph_traces[i].end_edge();
      };
    };

    next_reporter.draw_motion_graph( free_space, g, steer_rec );
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
    const shared_ptr< seq_trajectory_base< typename subspace_traits< FreeSpaceType >::super_space_type > >& traj )
    const {
    typedef typename topology_traits< FreeSpaceType >::point_type PointType;
    typedef typename seq_trajectory_base< typename subspace_traits< FreeSpaceType >::super_space_type >::
      point_time_iterator PtIter;

    std::vector< geom::tracer_coin3d_impl > current_trace;
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace.push_back( geom::tracer_coin3d_impl( true ) );
    std::cout << "Starting to register solution..." << std::endl;
    double t = 0.0;
    PtIter it = traj->begin_time_travel();
    PointType last_pt = *it;
    mdl_applicator->apply_to_model( last_pt, free_space.get_super_space() );
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace[i].begin_edge( traced_frames[i]->getGlobalPose().Position );
    while( it != traj->end_time_travel() ) {
      t += interval_size;
      it += interval_size;
      last_pt = *it;
      mdl_applicator->apply_to_model( last_pt, free_space.get_super_space() );
      for( std::size_t i = 0; i < traced_frames.size(); ++i )
        current_trace[i].add_point( traced_frames[i]->getGlobalPose().Position );
    };
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace[i].end_edge();

    if( solution_traces.empty() || ( t <= solution_traces.begin()->first ) )
      solution_traces[t].swap( current_trace );

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
    const shared_ptr< seq_path_base< typename subspace_traits< FreeSpaceType >::super_space_type > >& p ) const {
    typedef typename topology_traits< FreeSpaceType >::point_type PointType;
    typedef typename seq_path_base< typename subspace_traits< FreeSpaceType >::super_space_type >::
      point_fraction_iterator PtIter;

    std::vector< geom::tracer_coin3d_impl > current_trace;
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace.push_back( geom::tracer_coin3d_impl( true ) );

    double t = 0.0;
    PtIter it = p->begin_fraction_travel();
    PointType last_pt = *it;
    mdl_applicator->apply_to_model( last_pt, free_space.get_super_space() );
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace[i].begin_edge( traced_frames[i]->getGlobalPose().Position );
    while( it != p->end_fraction_travel() ) {
      it += 0.1;
      t += get( distance_metric, free_space.get_super_space() )( last_pt, *it, free_space.get_super_space() );
      last_pt = *it;
      mdl_applicator->apply_to_model( last_pt, free_space.get_super_space() );
      for( std::size_t i = 0; i < traced_frames.size(); ++i )
        current_trace[i].add_point( traced_frames[i]->getGlobalPose().Position );
    };
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      current_trace[i].end_edge();

    if( solution_traces.empty() || ( t <= solution_traces.begin()->first ) )
      solution_traces[t].swap( current_trace );

    next_reporter.draw_solution( free_space, p );
  };

  void reset_internal_state() {
    motion_graph_traces.clear();
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      motion_graph_traces.push_back( geom::tracer_coin3d_impl( false ) );
    solution_traces.clear();

    next_reporter.reset_internal_state();
  };

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    shared_object::save( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( next_reporter ) & RK_SERIAL_SAVE_WITH_NAME( mdl_applicator )
      & RK_SERIAL_SAVE_WITH_NAME( interval_size ) & RK_SERIAL_SAVE_WITH_NAME( traced_frames );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    shared_object::load( A, shared_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( next_reporter ) & RK_SERIAL_LOAD_WITH_NAME( mdl_applicator )
      & RK_SERIAL_LOAD_WITH_NAME( interval_size ) & RK_SERIAL_LOAD_WITH_NAME( traced_frames );
    motion_graph_traces.clear();
    for( std::size_t i = 0; i < traced_frames.size(); ++i )
      motion_graph_traces.push_back( geom::tracer_coin3d_impl( false ) );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC246000A, 1, "frame_tracer_3D", shared_object )
};
};
};


#endif
