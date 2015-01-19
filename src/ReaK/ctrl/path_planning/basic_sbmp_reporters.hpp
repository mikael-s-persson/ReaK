/**
 * \file basic_sbmp_reporters.hpp
 * 
 * This library defines simple sampling-based motion/path planning reporters. 
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

#ifndef REAK_SIMPLE_SBMP_REPORTERS_HPP
#define REAK_SIMPLE_SBMP_REPORTERS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include <ReaK/ctrl/topologies/subspace_concept.hpp>
#include <ReaK/ctrl/topologies/steerable_space_concept.hpp>
#include <ReaK/ctrl/interpolation/seq_trajectory_base.hpp>
#include <ReaK/ctrl/interpolation/seq_path_base.hpp>

#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and reports nothing.
 */
struct no_sbmp_report : public shared_object {
  
  void reset_internal_state() { };
  
  /**
   * Draws the entire motion-graph.
   * \tparam FreeSpaceType The C-free topology type.
   * \tparam MotionGraph The graph structure type representing the motion-graph.
   * \tparam PositionMap The property-map type that can map motion-graph vertex descriptors into point values.
   */
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename PositionMap>
  void draw_motion_graph(const FreeSpaceType&, const MotionGraph&, PositionMap) const { };
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType&, 
                     const shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >&) const { };
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType&, 
                     const shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >&) const { };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
  };
  
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
  };
  
  RK_RTTI_MAKE_CONCRETE_1BASE(no_sbmp_report,0xC2460002,1,"no_sbmp_report",shared_object)
  
  
};


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the underlying C-free (free_space) to draw individual edges of the motion graph or 
 * solution trajectory. The underlying space should have the functions:
 * 
 * void reset_output() const;
 * 
 * void draw_edge(const point_type& a, const point_type& b, bool is_solution_path) const;
 * 
 * void save_output(const std::string& filename) const;
 */
template <typename NextReporter = no_sbmp_report>
struct differ_sbmp_report_to_space : public shared_object {
  typedef differ_sbmp_report_to_space<NextReporter> self;
  
  NextReporter next_reporter;
  /// Holds the interval-size between output points of the solution trajectory/path.
  double interval_size;
  /// Holds the file-path where to output the reports.
  std::string file_path;
  
  mutable std::size_t progress_count;
  mutable std::size_t solution_count;
  
  explicit differ_sbmp_report_to_space(const std::string& aFilePath = "", 
                                       double aIntervalSize = 0.1,
                                       NextReporter aNextReporter = NextReporter()) : 
                                       next_reporter(aNextReporter),
                                       interval_size(aIntervalSize),
                                       file_path(aFilePath),
                                       progress_count(0),
                                       solution_count(0) { };
  
  void reset_internal_state() { 
    progress_count = 0;
    solution_count = 0;
    
    next_reporter.reset_internal_state();
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
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename PositionMap>
  typename boost::disable_if< is_steerable_space<FreeSpaceType>,
  void >::type draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos) const {
    typedef typename boost::graph_traits<MotionGraph>::vertex_iterator VIter;
    typedef typename boost::graph_traits<MotionGraph>::out_edge_iterator EIter;
    
    free_space.reset_output();
    
    VIter vi, vi_end;
    for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
      EIter ei, ei_end;
      for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei)
        free_space.draw_edge(get(pos, *vi), get(pos, target(*ei, g)), false);
    };
    
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << (progress_count++);
    free_space.save_output(file_path + "progress_" + ss.str());
    
    next_reporter.draw_motion_graph(free_space, g, pos);
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
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename SteerRecMap>
  typename boost::enable_if< is_steerable_space<FreeSpaceType>,
  void >::type draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g, SteerRecMap steer_rec) const {
    typedef typename boost::graph_traits<MotionGraph>::vertex_iterator VIter;
    typedef typename boost::graph_traits<MotionGraph>::out_edge_iterator EIter;
    typedef typename topology_traits<FreeSpaceType>::point_type PointType;
    typedef typename boost::property_traits<SteerRecMap>::value_type SteerRecType;
    typedef typename SteerRecType::point_fraction_iterator SteerIter;
    
    free_space.reset_output();
    
    VIter vi, vi_end;
    for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
      EIter ei, ei_end;
      for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
        const SteerRecType& st_rec = get(steer_rec, *ei);
        SteerIter it = st_rec.begin_fraction_travel();
        SteerIter prev_it = it; it += 0.1;
        for(; prev_it != st_rec.end_fraction_travel(); it += 0.1) {
          free_space.draw_edge(PointType(*prev_it), PointType(*it), false);
          prev_it = it;
        };
      };
    };
    
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << (progress_count++);
    free_space.save_output(file_path + "progress_" + ss.str());
    
    next_reporter.draw_motion_graph(free_space, g, steer_rec);
  };
  
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
    typedef typename topology_traits<FreeSpaceType>::point_type PointType;
    
    free_space.reset_output();
    
    double t = traj->get_start_time();
    PointType u_pt = traj->get_point_at_time(t);
    PointType v_pt;
    while(t < traj->get_end_time()) {
      t += interval_size;
      v_pt = traj->get_point_at_time(t);
      free_space.draw_edge(u_pt, v_pt, true);
      u_pt = v_pt;
    };
    
    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << (solution_count++) << "_" << (traj->get_end_time() - traj->get_start_time());
    free_space.save_output(file_path + "solution_" + ss.str());
    
    next_reporter.draw_solution(free_space, traj);
  };
  
  /**
   * Draws the solution path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param p The solution path.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
    typedef typename seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_distance_iterator PtIter;
    
    free_space.reset_output();
    
    double total_dist = 0.0;
    
    PtIter u_it = p->begin_distance_travel();
    PtIter v_it = u_it; v_it += interval_size;
    while(u_it != p->end_distance_travel()) {
      total_dist += interval_size;
      free_space.draw_edge(*u_it, *v_it, true);
      u_it = v_it;
      v_it += interval_size;
    };
    
    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << (solution_count++) << "_" << total_dist;
    free_space.save_output(file_path + "solution_" + ss.str());
    
    next_reporter.draw_solution(free_space, p);
  };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_SAVE_WITH_NAME(next_reporter)
      & RK_SERIAL_SAVE_WITH_NAME(interval_size)
      & RK_SERIAL_SAVE_WITH_NAME(file_path);
  };
  
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_LOAD_WITH_NAME(next_reporter)
      & RK_SERIAL_LOAD_WITH_NAME(interval_size)
      & RK_SERIAL_LOAD_WITH_NAME(file_path);
  };
  
  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460003,1,"differ_sbmp_report_to_space",shared_object)
  
};


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and simply times the progress of the planner, which it outputs to a given output stream.
 */
template <typename NextReporter = no_sbmp_report>
struct timing_sbmp_report : public shared_object {
  typedef timing_sbmp_report<NextReporter> self;
  
  NextReporter next_reporter;
  /// Holds the interval-size between output points of the solution trajectory/path.
  mutable boost::posix_time::ptime last_time;
  std::ostream* p_out;
  
  explicit timing_sbmp_report(std::ostream& aOutStream = std::cout,
                              NextReporter aNextReporter = NextReporter()) : 
                              next_reporter(aNextReporter),
                              last_time(boost::posix_time::microsec_clock::local_time()),
                              p_out(&aOutStream) { };
  
  void reset_internal_state() {
    last_time = boost::posix_time::microsec_clock::local_time();
    
    next_reporter.reset_internal_state();
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
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos) const {
    
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - last_time;
    (*p_out) << num_vertices(g) << " " << dt.total_microseconds() << std::endl;
    
    next_reporter.draw_motion_graph(free_space, g, pos);
    
    last_time = boost::posix_time::microsec_clock::local_time();
  };
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
    next_reporter.draw_solution(free_space, traj);
  };
  
  /**
   * Draws the solution path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param p The solution path.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
    next_reporter.draw_solution(free_space, p);
  };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  };
  
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_LOAD_WITH_NAME(next_reporter);
    last_time = boost::posix_time::microsec_clock::local_time();
  };
  
  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460004,1,"timing_sbmp_report",shared_object)
  
};


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and simply prints the progress of the planner.
 */
template <typename NextReporter = no_sbmp_report>
struct print_sbmp_progress : public shared_object {
  typedef print_sbmp_progress<NextReporter> self;
  
  NextReporter next_reporter;
  
  explicit print_sbmp_progress(NextReporter aNextReporter = NextReporter()) : 
                               next_reporter(aNextReporter) { };
  
  void reset_internal_state() { 
    next_reporter.reset_internal_state();
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
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos) const {
    std::cout << "\r" << std::setw(15) << num_vertices(g);
    std::cout.flush();
     
    next_reporter.draw_motion_graph(free_space, g, pos);
  };
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
    std::cout << "Solution Found!" << std::endl;
    
    next_reporter.draw_solution(free_space, traj);
  };
  
  /**
   * Draws the solution path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param p The solution path.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
    std::cout << "Solution Found!" << std::endl;
    
    next_reporter.draw_solution(free_space, p);
  };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  };
  
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_LOAD_WITH_NAME(next_reporter);
  };
  
  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460005,1,"print_sbmp_progress",shared_object)
  
};



/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and records the current best solution for each progress interval, which it outputs to a given output stream.
 */
template <typename NextReporter = no_sbmp_report>
struct least_cost_sbmp_report : public shared_object {
  typedef least_cost_sbmp_report<NextReporter> self;
  
  NextReporter next_reporter;
  std::ostream* p_out;
  std::ostream* p_sol;
  mutable double current_best;
  mutable std::size_t last_node_count;
  
  explicit least_cost_sbmp_report(std::ostream& aOutStream = std::cout,
                                  std::ostream* aPSolutionOutput = NULL,
                                  NextReporter aNextReporter = NextReporter()) : 
                                  next_reporter(aNextReporter),
                                  p_out(&aOutStream),
                                  p_sol(aPSolutionOutput),
                                  current_best(1e10), 
                                  last_node_count(0) { };
  
  void reset_internal_state() { 
    current_best = 1e10;
    last_node_count = 0;
    
    next_reporter.reset_internal_state();
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
  template <typename FreeSpaceType,
            typename MotionGraph,
            typename PositionMap>
  void draw_motion_graph(const FreeSpaceType& free_space, const MotionGraph& g, PositionMap pos) const {
    last_node_count = num_vertices(g);
    (*p_out) << num_vertices(g) << " " << current_best << std::endl;
    
    next_reporter.draw_motion_graph(free_space, g, pos);
  };
  
  /**
   * Draws the solution trajectory.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param traj The solution trajectory.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
//     typedef typename seq_trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_time_iterator TIter;
//     double total_cost = 0.0;
//     TIter it_prev = traj->begin_time_travel();
//     TIter it = traj->begin_time_travel(); it += 0.1;
//     for(; it != traj->end_time_travel(); it += 0.1) {
//       total_cost += get(distance_metric, free_space.get_super_space())(*it_prev, *it, free_space.get_super_space());
//     };
//     double total_cost = get(distance_metric, free_space.get_super_space())(
//       *(traj->begin_time_travel()), *(traj->end_time_travel()), free_space.get_super_space());
    double total_cost = traj->travel_distance(*(traj->begin_time_travel()), *(traj->end_time_travel()));
    if(total_cost < current_best)
      current_best = total_cost;
    
    if(p_sol)
      (*p_sol) << last_node_count << " " << current_best << std::endl;
    
    next_reporter.draw_solution(free_space, traj);
  };
  
  /**
   * Draws the solution path.
   * \tparam FreeSpaceType The C-free topology type.
   * \param free_space The C-free topology.
   * \param p The solution path.
   */
  template <typename FreeSpaceType>
  void draw_solution(const FreeSpaceType& free_space, 
                     const shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
    typedef typename seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_fraction_iterator FIter;
    double total_cost = 0.0;
    for(FIter it = p->begin_fraction_travel(); it != p->end_fraction_travel(); ) {
      FIter it_next = it; it_next += 1.0;
      total_cost += get(distance_metric, free_space.get_super_space())(*it, *it_next, free_space.get_super_space());
      it = it_next;
    };
    if(total_cost < current_best)
      current_best = total_cost;
    
    if(p_sol)
      (*p_sol) << last_node_count << " " << current_best << std::endl;
    
    next_reporter.draw_solution(free_space, p);
  };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
    shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_SAVE_WITH_NAME(next_reporter);
  };
  
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
    shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    A & RK_SERIAL_LOAD_WITH_NAME(next_reporter);
    current_best = 1e10;
    last_node_count = 0;
  };
  
  RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460006,1,"least_cost_sbmp_report",shared_object)
  
};



};

};


#endif


