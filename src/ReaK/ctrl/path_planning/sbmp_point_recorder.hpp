/**
 * \file sbmp_point_recorder.hpp
 * 
 * This library defines a sampling-based motion/path planning reporter for tracing out the points along a motion graph
 * and solution paths.
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

#ifndef REAK_SBMP_POINT_RECORDER_HPP
#define REAK_SBMP_POINT_RECORDER_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "trajectory_base.hpp"
#include "basic_sbmp_reporters.hpp"
#include <boost/graph/graph_concepts.hpp>

#include "topologies/direct_kinematics_topomap.hpp"
#include "topological_map_concepts.hpp"

#include "recorders/data_record.hpp"
#include "recorders/ssv_recorder.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the recorder library to trace out the points of a motion-graph and solution paths
 * as space-separated files (i.e., loadable into Matlab / Octave).
 */
template <typename JointStateSpace, typename JointStateMapping = identity_topo_map, typename NextReporter = no_sbmp_report>
class sbmp_point_recorder : public shared_object {
  public:
    typedef sbmp_point_recorder<JointStateSpace, JointStateMapping, NextReporter> self;
    
    NextReporter next_reporter;
    
  protected:
    shared_ptr<JointStateSpace> jt_space;
    JointStateMapping map_to_jt_space;
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    /// Holds the file-path where to output the reports.
    std::string file_path;
    
  public:
    
    explicit sbmp_point_recorder(const shared_ptr<JointStateSpace>& aJointSpace = shared_ptr<JointStateSpace>(),
                                 const JointStateMapping& aMapToJtSpace = JointStateMapping(),
                                 const std::string& aFilePath = "", 
                                 double aIntervalSize = 0.1,
                                 NextReporter aNextReporter = NextReporter()) : 
                                 next_reporter(aNextReporter),
                                 jt_space(aJointSpace),
                                 map_to_jt_space(aMapToJtSpace),
                                 interval_size(aIntervalSize),
                                 file_path(aFilePath) { };
    
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
      typedef typename boost::graph_traits<MotionGraph>::vertex_iterator VIter;
      typedef typename boost::graph_traits<MotionGraph>::out_edge_iterator EIter;
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      using ReaK::to_vect;
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
      recorder::ssv_recorder rec_out(file_path + "progress_" + ss.str());
      bool not_initialized_yet = true;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        PointType p_u = get(pos, *vi);
        JointState s_u = map_to_jt_space.map_to_space(p_u, free_space.get_super_space(), *jt_space);
        vect_n<double> v_s_u = to_vect<double>(s_u);
        
        if(not_initialized_yet) {
          for(std::size_t i = 0; i < v_s_u.size(); ++i) {
            std::stringstream ss2;
            ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
            rec_out << ss2.str();
          };
          rec_out << recorder::data_recorder::end_name_row;
          not_initialized_yet = false;
        };
        for(std::size_t i = 0; i < v_s_u.size(); ++i)
          rec_out << v_s_u[i];
        rec_out << recorder::data_recorder::end_value_row;
        
        for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
          PointType p_v = get(pos, target(*ei, g));
          JointState s_v = map_to_jt_space.map_to_space(p_v, free_space.get_super_space(), *jt_space);
          
          dk_map.apply_to_model(s_u, *jt_space);
          for(std::size_t i = 0; i < traced_frames.size(); ++i)
            motion_graph_traces[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
          for(double j = 0.1; j <= 1.01; j += 0.1) {
            PointType p_new = free_space.get_super_space().move_position_toward(p_u, j, p_v);
            JointState s_new = map_to_jt_space.map_to_space(p_new, free_space.get_super_space(), *jt_space);
            vect_n<double> v_s_new = to_vect<double>(s_new);
            
            dk_map.apply_to_model(s_new, *jt_space);
            for(std::size_t i = 0; i < traced_frames.size(); ++i)
              motion_graph_traces[i].add_point(traced_frames[i]->getGlobalPose().Position);
          };
          for(std::size_t i = 0; i < traced_frames.size(); ++i)
            motion_graph_traces[i].end_edge();
          
        };
      };
      
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
                       const shared_ptr< trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      
      double t_total = traj->get_end_time() - traj->get_start_time();
      if(!(solution_traces.empty()) && t_total > solution_traces.begin()->first) {
        next_reporter.draw_solution(free_space, traj);
        return;
      };
      std::vector< geom::tracer_coin3d_impl >& current_trace = solution_traces[t_total];
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = traj->get_start_time();
      PointType u_pt = traj->get_point_at_time(t);
      JointState s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
      dk_map.apply_to_model(s_u, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(t < traj->get_end_time()) {
        t += interval_size;
        u_pt = traj->get_point_at_time(t);
        s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
        
        dk_map.apply_to_model(s_u, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
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
      typedef typename topology_traits<FreeSpaceType>::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      typedef typename seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_fraction_iterator PtIter;
      
      std::vector< geom::tracer_coin3d_impl > current_trace;
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = 0.0;
      PtIter it = p->begin_fraction_travel();
      PointType last_pt = *it;
      JointState s_u = map_to_jt_space.map_to_space(last_pt, free_space.get_super_space(), *jt_space);
      dk_map.apply_to_model(s_u, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(it != p->end_fraction_travel()) {
        it += 0.1;
        t += get(distance_metric, free_space.get_super_space())(last_pt, *it, free_space.get_super_space());
        last_pt = *it;
        s_u = map_to_jt_space.map_to_space(last_pt, free_space.get_super_space(), *jt_space);
        dk_map.apply_to_model(s_u, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
      if(solution_traces.empty() || (t <= solution_traces.begin()->first))
        solution_traces[t].swap(current_trace);
      
      next_reporter.draw_solution(free_space, p);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(next_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(map_to_jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(interval_size)
        & RK_SERIAL_SAVE_WITH_NAME(file_path);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(next_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(map_to_jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(interval_size)
        & RK_SERIAL_LOAD_WITH_NAME(file_path);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000F,1,"sbmp_point_recorder",shared_object)
    
};




/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the recorder library to trace out the points of a motion-graph and solution paths
 * as space-separated files (i.e., loadable into Matlab / Octave).
 */
template <typename JointStateSpace, typename NextReporter>
class sbmp_point_recorder<JointStateSpace, identity_topo_map, NextReporter> : public shared_object {
  public:
    typedef sbmp_point_recorder<JointStateSpace, identity_topo_map, NextReporter> self;
    
    NextReporter next_reporter;
    
  protected:
    shared_ptr<JointStateSpace> jt_space;
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    /// Holds the file-path where to output the reports.
    std::string file_path;
  
  public:
    
    explicit sbmp_point_recorder(const DirectKinMapper& aDKMap = DirectKinMapper(),
                                 const shared_ptr<JointStateSpace>& aJointSpace = shared_ptr<JointStateSpace>(),
                                 const std::string& aFilePath = "", 
                                 double aIntervalSize = 0.1,
                                 NextReporter aNextReporter = NextReporter()) : 
                                 next_reporter(aNextReporter),
                                 jt_space(aJointSpace),
                                 interval_size(aIntervalSize),
                                 file_path(aFilePath) { };
    
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
      typedef typename boost::graph_traits<MotionGraph>::vertex_iterator VIter;
      typedef typename boost::graph_traits<MotionGraph>::out_edge_iterator EIter;
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        PointType p_u = get(pos, *vi);
        for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
          PointType p_v = get(pos, target(*ei, g));
          dk_map.apply_to_model(p_u, *jt_space);
          for(std::size_t i = 0; i < traced_frames.size(); ++i)
            motion_graph_traces[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
          for(double j = 0.1; j <= 1.01; j += 0.1) {
            PointType p_new = free_space.get_super_space().move_position_toward(p_u, j, p_v);
            dk_map.apply_to_model(p_new, *jt_space);
            for(std::size_t i = 0; i < traced_frames.size(); ++i)
              motion_graph_traces[i].add_point(traced_frames[i]->getGlobalPose().Position);
          };
          for(std::size_t i = 0; i < traced_frames.size(); ++i)
            motion_graph_traces[i].end_edge();
        };
      };
      
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
                       const shared_ptr< trajectory_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& traj) const {
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      
      double t_total = traj->get_end_time() - traj->get_start_time();
      if(!(solution_traces.empty()) && t_total > solution_traces.begin()->first) {
        next_reporter.draw_solution(free_space, traj);
        return;
      };
      std::vector< geom::tracer_coin3d_impl >& current_trace = solution_traces[t_total];
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = traj->get_start_time();
      PointType u_pt = traj->get_point_at_time(t);
      dk_map.apply_to_model(u_pt, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(t < traj->get_end_time()) {
        t += interval_size;
        u_pt = traj->get_point_at_time(t);
        dk_map.apply_to_model(u_pt, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
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
      typedef typename topology_traits<FreeSpaceType>::point_type PointType;
      typedef typename seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_fraction_iterator PtIter;
      
      std::vector< geom::tracer_coin3d_impl > current_trace;
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = 0.0;
      PtIter it = p->begin_fraction_travel();
      PointType last_pt = *it;
      dk_map.apply_to_model(last_pt, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(it != p->end_fraction_travel()) {
        it += 0.1;
        t += get(distance_metric, free_space.get_super_space())(last_pt, *it, free_space.get_super_space());
        last_pt = *it;
        dk_map.apply_to_model(last_pt, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
      if(solution_traces.empty() || (t <= solution_traces.begin()->first))
        solution_traces[t].swap(current_trace);
      
      next_reporter.draw_solution(free_space, p);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(next_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(interval_size)
        & RK_SERIAL_SAVE_WITH_NAME(file_path);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(next_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(interval_size)
        & RK_SERIAL_LOAD_WITH_NAME(file_path);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000F,1,"sbmp_point_recorder",shared_object)
    
};




};

};


#endif


