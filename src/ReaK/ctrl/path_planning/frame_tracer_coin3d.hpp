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

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "trajectory_base.hpp"
#include "basic_sbmp_reporters.hpp"
#include <boost/graph/graph_concepts.hpp>

#include "topologies/direct_kinematics_topomap.hpp"
#include "topological_map_concepts.hpp"

#include "shapes/frame_tracer_coin3d_impl.hpp"

class SoSeparator;   // forward-declare


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the Coin3D library to trace out the position of a number of 3D poses linked to a model
 * linked with the given DirectKinMapper.
 */
template <typename DirectKinMapper, typename JointStateSpace, typename JointStateMapping = identity_topo_map, typename NextReporter = no_sbmp_report>
class frame_tracer_3D : public shared_object {
  public:
    typedef frame_tracer_3D<DirectKinMapper, JointStateSpace, JointStateMapping, NextReporter> self;
    
    NextReporter next_reporter;
    
  protected:
    DirectKinMapper dk_map;
    shared_ptr<JointStateSpace> jt_space;
    JointStateMapping map_to_jt_space;
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    
    std::vector< shared_ptr< pose_3D<double> > > traced_frames;
    std::vector< geom::tracer_coin3d_impl > motion_graph_traces;
    mutable std::map< double, std::vector< geom::tracer_coin3d_impl > > solution_traces;
  
  public:
    
    explicit frame_tracer_3D(const DirectKinMapper& aDKMap = DirectKinMapper(),
                             const shared_ptr<JointStateSpace>& aJointSpace = shared_ptr<JointStateSpace>(),
                             const JointStateMapping& aMapToJtSpace = JointStateMapping(),
                             double aIntervalSize = 0.1,
                             NextReporter aNextReporter = NextReporter()) : 
                             next_reporter(aNextReporter),
                             dk_map(aDKMap),
                             jt_space(aJointSpace),
                             map_to_jt_space(aMapToJtSpace),
                             interval_size(aIntervalSize) { };
    
    self& add_traced_frame(const shared_ptr< pose_3D<double> >& aPose) {
      traced_frames.push_back(aPose);
      motion_graph_traces.push_back(geom::tracer_coin3d_impl(false));
      return *this;
    };
    
    const geom::tracer_coin3d_impl& get_motion_graph_tracer(const shared_ptr< pose_3D<double> >& aPose) const {
      std::vector< shared_ptr< pose_3D<double> > >::const_iterator it = std::find(traced_frames.begin(),traced_frames.end(),aPose);
      if(it != traced_frames.end()) {
        return motion_graph_traces[it - traced_frames.begin()];
      } else
        throw std::range_error("The given pose is not being traced by this frame-tracer!");
    };
    
    std::size_t get_solution_count() const {
      return solution_traces.size();
    };
    
    double get_best_solution_value() const {
      if(solution_traces.empty())
        return std::numeric_limits<double>::infinity();
      return solution_traces.begin()->first;
    };
    
    const geom::tracer_coin3d_impl& get_solution_tracer(const shared_ptr< pose_3D<double> >& aPose, std::size_t aSolutionId = 0) const {
      if(aSolutionId >= solution_traces.size())
        aSolutionId = 0;
      std::vector< shared_ptr< pose_3D<double> > >::const_iterator it = std::find(traced_frames.begin(),traced_frames.end(),aPose);
      if(it != traced_frames.end()) {
        std::map< double, std::vector< geom::tracer_coin3d_impl > >::const_iterator itm = solution_traces.begin();
        std::advance(itm, aSolutionId);
        return itm->second[it - traced_frames.begin()];
      } else
        throw std::range_error("The given pose is not being traced by this frame-tracer!");
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
      typedef typename boost::graph_traits<MotionGraph>::vertex_iterator VIter;
      typedef typename boost::graph_traits<MotionGraph>::out_edge_iterator EIter;
      typedef typename topology_traits< FreeSpaceType >::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      
      VIter vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi) {
        EIter ei, ei_end;
        PointType p_u = get(pos, *vi);
        JointState s_u = map_to_jt_space.map_to_space(p_u, free_space.get_super_space(), *jt_space);
        for(boost::tie(ei,ei_end) = out_edges(*vi,g); ei != ei_end; ++ei) {
          PointType p_v = get(pos, target(*ei, g));
          JointState s_v = map_to_jt_space.map_to_space(p_v, free_space.get_super_space(), *jt_space);
          dk_map.apply_to_model(s_u, *jt_space);
          for(std::size_t i = 0; i < traced_frames.size(); ++i)
            motion_graph_traces[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
          for(double j = 0.1; j <= 1.01; j += 0.1) {
            PointType p_new = free_space.get_super_space().move_position_toward(p_u, j, p_v);
            JointState s_new = map_to_jt_space.map_to_space(p_new, free_space.get_super_space(), *jt_space);
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
                       const shared_ptr< path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
      typedef typename topology_traits<FreeSpaceType>::point_type PointType;
      typedef typename topology_traits< JointStateSpace >::point_type JointState;
      
      double t_total = p->travel_distance(p->get_start_point(), p->get_end_point());
      if(!(solution_traces.empty()) && t_total > solution_traces.begin()->first) {
        next_reporter.draw_solution(free_space, p);
        return;
      };
      std::vector< geom::tracer_coin3d_impl >& current_trace = solution_traces[t_total];
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = 0.0;
      PointType u_pt = p->get_start_point();
      JointState s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
      JointState s_prev = map_to_jt_space.map_to_space(p->get_end_point(), free_space.get_super_space(), *jt_space);
      double total_dist = get(distance_metric, *jt_space)(s_u, s_prev, *jt_space);
      dk_map.apply_to_model(s_u, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(get(distance_metric, *jt_space)(s_prev, s_u, *jt_space) > 0.01 * total_dist) {
      //while(get(distance_metric, free_space.get_super_space())(u_pt, p->get_end_point(), free_space.get_super_space()) > interval_size) {
        s_prev = s_u;
        t += interval_size;
        u_pt = p->move_away_from(u_pt, interval_size);
        s_u = map_to_jt_space.map_to_space(u_pt, free_space.get_super_space(), *jt_space);
        std::cout << "s_u = " << s_u << std::endl;
        std::cout << "s_prev = " << s_prev << std::endl;
        dk_map.apply_to_model(s_u, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
      next_reporter.draw_solution(free_space, p);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(next_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(dk_map)
        & RK_SERIAL_SAVE_WITH_NAME(jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(map_to_jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(interval_size)
        & RK_SERIAL_SAVE_WITH_NAME(traced_frames);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(next_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(dk_map)
        & RK_SERIAL_LOAD_WITH_NAME(jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(map_to_jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(interval_size)
        & RK_SERIAL_LOAD_WITH_NAME(traced_frames);
      motion_graph_traces.clear();
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        motion_graph_traces.push_back(geom::tracer_coin3d_impl(false));
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000A,1,"frame_tracer_3D",shared_object)
    
};




/**
 * This class can be used as a SBMP/SBPP Reporter (SBMPReporterConcept and SBPPReporterConcept) 
 * and uses the Coin3D library to trace out the position of a number of 3D poses linked to a model
 * linked with the given DirectKinMapper.
 */
template <typename DirectKinMapper, typename JointStateSpace, typename NextReporter>
class frame_tracer_3D<DirectKinMapper, JointStateSpace, identity_topo_map, NextReporter> : public shared_object {
  public:
    typedef frame_tracer_3D<DirectKinMapper, JointStateSpace, identity_topo_map, NextReporter> self;
    
    NextReporter next_reporter;
    
  protected:
    DirectKinMapper dk_map;
    shared_ptr<JointStateSpace> jt_space;
    /// Holds the interval-size between output points of the solution trajectory/path.
    double interval_size;
    
    std::vector< shared_ptr< pose_3D<double> > > traced_frames;
    std::vector< geom::tracer_coin3d_impl > motion_graph_traces;
    mutable std::map< double, std::vector< geom::tracer_coin3d_impl > > solution_traces;
  
  public:
    
    explicit frame_tracer_3D(const DirectKinMapper& aDKMap = DirectKinMapper(),
                             const shared_ptr<JointStateSpace>& aJointSpace = shared_ptr<JointStateSpace>(),
                             double aIntervalSize = 0.1,
                             NextReporter aNextReporter = NextReporter()) : 
                             next_reporter(aNextReporter),
                             dk_map(aDKMap),
                             jt_space(aJointSpace),
                             interval_size(aIntervalSize) { };
    
    self& add_traced_frame(const shared_ptr< pose_3D<double> >& aPose) {
      traced_frames.push_back(aPose);
      motion_graph_traces.push_back(geom::tracer_coin3d_impl(false));
      return *this;
    };
    
    const geom::tracer_coin3d_impl& get_motion_graph_tracer(const shared_ptr< pose_3D<double> >& aPose) const {
      std::vector< shared_ptr< pose_3D<double> > >::const_iterator it = std::find(traced_frames.begin(),traced_frames.end(),aPose);
      if(it != traced_frames.end()) {
        return motion_graph_traces[it - traced_frames.begin()];
      } else
        throw std::range_error("The given pose is not being traced by this frame-tracer!");
    };
    
    std::size_t get_solution_count() const {
      return solution_traces.size();
    };
    
    double get_best_solution_value() const {
      if(solution_traces.empty())
        return std::numeric_limits<double>::infinity();
      return solution_traces.begin()->first;
    };
    
    const geom::tracer_coin3d_impl& get_solution_tracer(const shared_ptr< pose_3D<double> >& aPose, std::size_t aSolutionId = 0) const {
      if(aSolutionId >= solution_traces.size())
        aSolutionId = 0;
      std::vector< shared_ptr< pose_3D<double> > >::const_iterator it = std::find(traced_frames.begin(),traced_frames.end(),aPose);
      if(it != traced_frames.end()) {
        std::map< double, std::vector< geom::tracer_coin3d_impl > >::const_iterator itm = solution_traces.begin();
        std::advance(itm, aSolutionId);
        return itm->second[it - traced_frames.begin()];
      } else
        throw std::range_error("The given pose is not being traced by this frame-tracer!");
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
                       const shared_ptr< path_base< typename subspace_traits<FreeSpaceType>::super_space_type > >& p) const {
      typedef typename topology_traits<FreeSpaceType>::point_type PointType;
      
      double t_total = p->travel_distance(p->get_start_point(), p->get_end_point());
      if(!(solution_traces.empty()) && t_total > solution_traces.begin()->first) {
        next_reporter.draw_solution(free_space, p);
        return;
      };
      std::vector< geom::tracer_coin3d_impl >& current_trace = solution_traces[t_total];
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace.push_back(geom::tracer_coin3d_impl(true));
      
      double t = 0.0;
      PointType u_pt = p->get_start_point();
      dk_map.apply_to_model(u_pt, *jt_space);
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].begin_edge(traced_frames[i]->getGlobalPose().Position);
      while(u_pt != p->get_end_point()) {
        t += interval_size;
        u_pt = p->move_away_from(u_pt, interval_size);
        dk_map.apply_to_model(u_pt, *jt_space);
        for(std::size_t i = 0; i < traced_frames.size(); ++i)
          current_trace[i].add_point(traced_frames[i]->getGlobalPose().Position);
      };
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        current_trace[i].end_edge();
      
      next_reporter.draw_solution(free_space, p);
    };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(next_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(dk_map)
        & RK_SERIAL_SAVE_WITH_NAME(jt_space)
        & RK_SERIAL_SAVE_WITH_NAME(interval_size)
        & RK_SERIAL_SAVE_WITH_NAME(traced_frames);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(next_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(dk_map)
        & RK_SERIAL_LOAD_WITH_NAME(jt_space)
        & RK_SERIAL_LOAD_WITH_NAME(interval_size)
        & RK_SERIAL_LOAD_WITH_NAME(traced_frames);
      motion_graph_traces.clear();
      for(std::size_t i = 0; i < traced_frames.size(); ++i)
        motion_graph_traces.push_back(geom::tracer_coin3d_impl(false));
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000A,1,"frame_tracer_3D",shared_object)
    
};




};

};


#endif


