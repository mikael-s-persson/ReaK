/**
 * \file sbmp_point_recorder.hpp
 * 
 * This library defines a sampling-based motion/path planning reporter for recording the 
 * points along a motion graph and solution paths.
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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

#include <boost/concept_check.hpp>

#include <ReaK/ctrl/interpolation/seq_trajectory_base.hpp>
#include "basic_sbmp_reporters.hpp"
#include <boost/graph/graph_concepts.hpp>

#include <ReaK/ctrl/topologies/direct_kinematics_topomap.hpp>
#include <ReaK/ctrl/topologies/topological_map_concepts.hpp>

#include <ReaK/core/recorders/data_record.hpp>
#include <ReaK/core/recorders/ascii_recorder.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the recorder library to print out the points of a motion-graph and solution paths
 * as space-separated files (i.e., loadable into Matlab / Octave).
 * \tparam JointStateSpace A joint-state space type.
 * \tparam JointStateMapping A map type that can map the points of the path-planning topology into points of the joint-state space.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template <typename JointStateSpace, typename JointStateMapping = identity_topo_map, typename NextReporter = no_sbmp_report>
class sbmp_point_recorder : public shared_object {
  public:
    typedef sbmp_point_recorder<JointStateSpace, JointStateMapping, NextReporter> self;
    
    /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
    NextReporter next_reporter;
    
  protected:
    /// A shared-pointer to the joint-state space.
    shared_ptr<JointStateSpace> jt_space;
    /// A map that can map the points of the path-planning topology into points of the joint-state space.
    JointStateMapping map_to_jt_space;
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    /// Holds the file-path where to output the reports.
    std::string file_path;
    mutable std::size_t solution_count;
    
  public:
    
    /**
     * Parametrized constructor.
     * \param aJointSpace The shared-pointer to the joint-state space.
     * \param aMapToJtSpace The map that can map the points of the path-planning topology into points of the joint-state space.
     * \param aFilePath The path where to create the output files.
     * \param aIntervalSize The interval-size between output points of the solution trajectory/path.
     * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
     */
    explicit sbmp_point_recorder(const shared_ptr<JointStateSpace>& aJointSpace = shared_ptr<JointStateSpace>(),
                                 const JointStateMapping& aMapToJtSpace = JointStateMapping(),
                                 const std::string& aFilePath = "", 
                                 double aIntervalSize = 0.1,
                                 NextReporter aNextReporter = NextReporter()) : 
                                 next_reporter(aNextReporter),
                                 jt_space(aJointSpace),
                                 map_to_jt_space(aMapToJtSpace),
                                 interval_size(aIntervalSize),
                                 file_path(aFilePath), 
                                 solution_count(0) { };
    
    void reset_internal_state() { 
      solution_count = 0;
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      using ReaK::to_vect;
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "progress_" + ss.str());
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
          
          for(std::size_t i = 0; i < v_s_u.size(); ++i)
            rec_out << v_s_u[i];
          rec_out << recorder::data_recorder::end_value_row;
          
          for(double j = 0.1; j <= 1.01; j += 0.1) {
            PointType p_new = free_space.get_super_space().move_position_toward(p_u, j, p_v);
            JointState s_new = map_to_jt_space.map_to_space(p_new, free_space.get_super_space(), *jt_space);
            vect_n<double> v_s_new = to_vect<double>(s_new);
            for(std::size_t i = 0; i < v_s_new.size(); ++i)
              rec_out << v_s_new[i];
            rec_out << recorder::data_recorder::end_value_row;
          };
          
        };
      };
      rec_out << recorder::data_recorder::flush;
      rec_out << recorder::data_recorder::close;
      
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      typedef typename boost::property_traits<SteerRecMap>::value_type SteerRecType;
      typedef typename SteerRecType::point_fraction_iterator SteerIter;
      using ReaK::to_vect;
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "progress_" + ss.str());
      bool not_initialized_yet = true;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
          const SteerRecType& st_rec = get(steer_rec, *ei);
          for(SteerIter it = st_rec.begin_fraction_travel(); it != st_rec.end_fraction_travel(); it += 0.1) { 
            JointState s_new = map_to_jt_space.map_to_space(PointType(*it), free_space.get_super_space(), *jt_space);
            vect_n<double> v_s_new = to_vect<double>(s_new);
            
            if(not_initialized_yet) {
              for(std::size_t i = 0; i < v_s_new.size(); ++i) {
                std::stringstream ss2;
                ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
                rec_out << ss2.str();
              };
              rec_out << recorder::data_recorder::end_name_row;
              not_initialized_yet = false;
            };
            
            for(std::size_t i = 0; i < v_s_new.size(); ++i)
              rec_out << v_s_new[i];
            rec_out << recorder::data_recorder::end_value_row;
          };
        };
      };
      rec_out << recorder::data_recorder::flush;
      rec_out << recorder::data_recorder::close;
      
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      
      std::stringstream ss;
      ss << std::setw(3) << std::setfill('0') << (solution_count++) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "solution_" + ss.str());
      
      double t = traj->get_start_time();
      PointType u_pt = traj->get_point_at_time(t);
      JointState s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
      vect_n<double> v_s_u = to_vect<double>(s_u);
      
      for(std::size_t i = 0; i < v_s_u.size(); ++i) {
        std::stringstream ss2;
        ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
        rec_out << ss2.str();
      };
      rec_out << recorder::data_recorder::end_name_row;
      for(std::size_t i = 0; i < v_s_u.size(); ++i)
        rec_out << v_s_u[i];
      rec_out << recorder::data_recorder::end_value_row;
      
      while(t < traj->get_end_time()) {
        t += interval_size;
        u_pt = traj->get_point_at_time(t);
        s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
        v_s_u = to_vect<double>(s_u);
        for(std::size_t i = 0; i < v_s_u.size(); ++i)
          rec_out << v_s_u[i];
        rec_out << recorder::data_recorder::end_value_row;
      };
      rec_out << recorder::data_recorder::flush;
      rec_out << recorder::data_recorder::close;
      
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
      
      std::stringstream ss;
      ss << std::setw(3) << std::setfill('0') << (solution_count++) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "solution_" + ss.str());
      
      PtIter it = p->begin_fraction_travel();
      PointType last_pt = *it;
      JointState s_u = map_to_jt_space.map_to_space(last_pt, free_space.get_super_space(), *jt_space);
      vect_n<double> v_s_u = to_vect<double>(s_u);
      
      for(std::size_t i = 0; i < v_s_u.size(); ++i) {
        std::stringstream ss2;
        ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
        rec_out << ss2.str();
      };
      rec_out << recorder::data_recorder::end_name_row;
      for(std::size_t i = 0; i < v_s_u.size(); ++i)
        rec_out << v_s_u[i];
      rec_out << recorder::data_recorder::end_value_row;
      
      while(it != p->end_fraction_travel()) {
        it += 0.1;
        last_pt = *it;
        s_u = map_to_jt_space.map_to_space(last_pt, free_space.get_super_space(), *jt_space);
        v_s_u = to_vect<double>(s_u);
        for(std::size_t i = 0; i < v_s_u.size(); ++i)
          rec_out << v_s_u[i];
        rec_out << recorder::data_recorder::end_value_row;
      };
      rec_out << recorder::data_recorder::flush;
      rec_out << recorder::data_recorder::close;
      
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
      solution_count = 0;
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000F,1,"sbmp_point_recorder",shared_object)
    
};




/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the recorder library to print out the points of a motion-graph and solution paths
 * as space-separated files (i.e., loadable into Matlab / Octave).
 * \tparam JointStateSpace A joint-state space type.
 * \tparam NextReporter A SBMP/SBPP reporter type to chain to this reporter.
 */
template <typename JointStateSpace, typename NextReporter>
class sbmp_point_recorder<JointStateSpace, identity_topo_map, NextReporter> : public shared_object {
  public:
    typedef sbmp_point_recorder<JointStateSpace, identity_topo_map, NextReporter> self;
    
    /// Holds the instance of the SBMP/SBPP reporter to which calls are forwarded to.
    NextReporter next_reporter;
    
  protected:
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    /// Holds the file-path where to output the reports.
    std::string file_path;
    mutable std::size_t solution_count;
  
  public:
    
    /**
     * Parametrized constructor.
     * \param aFilePath The path where to create the output files.
     * \param aIntervalSize The interval-size between output points of the solution trajectory/path.
     * \param aNextReporter The instance of the SBMP/SBPP reporter to which calls are forwarded to.
     */
    explicit sbmp_point_recorder(const std::string& aFilePath = "", 
                                 double aIntervalSize = 0.1,
                                 NextReporter aNextReporter = NextReporter()) : 
                                 next_reporter(aNextReporter),
                                 interval_size(aIntervalSize),
                                 file_path(aFilePath),
                                 solution_count(0) { };
    
    void reset_internal_state() { 
      solution_count = 0;
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      using ReaK::to_vect;
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "progress_" + ss.str());
      bool not_initialized_yet = true;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        PointType p_u = get(pos, *vi);
        vect_n<double> v_s_u = to_vect<double>(p_u);
        
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
          
          for(std::size_t i = 0; i < v_s_u.size(); ++i)
            rec_out << v_s_u[i];
          rec_out << recorder::data_recorder::end_value_row;
          
          for(double j = 0.1; j <= 1.01; j += 0.1) {
            PointType p_new = free_space.get_super_space().move_position_toward(p_u, j, p_v);
            vect_n<double> v_s_new = to_vect<double>(p_new);
            for(std::size_t i = 0; i < v_s_new.size(); ++i)
              rec_out << v_s_new[i];
            rec_out << recorder::data_recorder::end_value_row;
          };
          
        };
      };
      rec_out << recorder::data_recorder::flush;
      
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename boost::property_traits<SteerRecMap>::value_type SteerRecType;
      typedef typename SteerRecType::point_fraction_iterator SteerIter;
      using ReaK::to_vect;
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << num_vertices(g) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "progress_" + ss.str());
      bool not_initialized_yet = true;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
          const SteerRecType& st_rec = get(steer_rec, *ei);
          for(SteerIter it = st_rec.begin_fraction_travel(); it != st_rec.end_fraction_travel(); it += 0.1) { 
            vect_n<double> v_s_new = to_vect<double>(PointType(*it));
            
            if(not_initialized_yet) {
              for(std::size_t i = 0; i < v_s_new.size(); ++i) {
                std::stringstream ss2;
                ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
                rec_out << ss2.str();
              };
              rec_out << recorder::data_recorder::end_name_row;
              not_initialized_yet = false;
            };
            
            for(std::size_t i = 0; i < v_s_new.size(); ++i)
              rec_out << v_s_new[i];
            rec_out << recorder::data_recorder::end_value_row;
          };
        };
      };
      rec_out << recorder::data_recorder::flush;
      
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
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      
      double t_total = traj->get_end_time() - traj->get_start_time();
      std::stringstream ss;
      ss << std::setw(3) << std::setfill('0') << (solution_count++) << "_" << t_total << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "solution_" + ss.str());
      
      double t = traj->get_start_time();
      PointType u_pt = traj->get_point_at_time(t);
      vect_n<double> v_s_u = to_vect<double>(u_pt);
      
      for(std::size_t i = 0; i < v_s_u.size(); ++i) {
        std::stringstream ss2;
        ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
        rec_out << ss2.str();
      };
      rec_out << recorder::data_recorder::end_name_row;
      for(std::size_t i = 0; i < v_s_u.size(); ++i)
        rec_out << v_s_u[i];
      rec_out << recorder::data_recorder::end_value_row;
      
      while(t < traj->get_end_time()) {
        t += interval_size;
        u_pt = traj->get_point_at_time(t);
        v_s_u = to_vect<double>(u_pt);
        for(std::size_t i = 0; i < v_s_u.size(); ++i)
          rec_out << v_s_u[i];
        rec_out << recorder::data_recorder::end_value_row;
      };
      rec_out << recorder::data_recorder::flush;
      
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
      typedef typename seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type >::point_fraction_iterator PtIter;
      
      std::stringstream ss;
      ss << std::setw(3) << std::setfill('0') << (solution_count++) << ".ssv";
      recorder::ascii_recorder rec_out(file_path + "solution_" + ss.str());
      
      PtIter it = p->begin_fraction_travel();
      vect_n<double> v_s_u = to_vect<double>(*it);
      
      for(std::size_t i = 0; i < v_s_u.size(); ++i) {
        std::stringstream ss2;
        ss2 << "state_" << std::setw(2) << std::setfill('0') << i;
        rec_out << ss2.str();
      };
      rec_out << recorder::data_recorder::end_name_row;
      for(std::size_t i = 0; i < v_s_u.size(); ++i)
        rec_out << v_s_u[i];
      rec_out << recorder::data_recorder::end_value_row;
      
      while(it != p->end_fraction_travel()) {
        it += 0.1;
        v_s_u = to_vect<double>(*it);
        for(std::size_t i = 0; i < v_s_u.size(); ++i)
          rec_out << v_s_u[i];
        rec_out << recorder::data_recorder::end_value_row;
      };
      rec_out << recorder::data_recorder::flush;
      
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
      solution_count = 0;
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000F,1,"sbmp_point_recorder",shared_object)
    
};




};

};


#endif


