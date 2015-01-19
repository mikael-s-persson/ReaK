/**
 * \file any_motion_graphs.hpp
 *
 * This library provides a type-erased classes for motion graphs.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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
 *    If not, sMOTION_ee <http://www.gnu.org/licenses/>.
 */

#ifndef RK_ANY_MOTION_GRAPHS_HPP
#define RK_ANY_MOTION_GRAPHS_HPP

#include <ReaK/core/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/ctrl/graph_alg/any_graph.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/steerable_space_concept.hpp>

#include <iostream>
#include <iomanip>

namespace ReaK {

namespace pp {


namespace detail {
  
  template <typename Topology, typename Graph>
  typename boost::enable_if< is_steerable_space<Topology>,
  steerable_space_traits<Topology> >::type::steer_record_type& 
    try_get_steer_record(Graph& g, typename boost::graph_traits< Graph >::edge_descriptor e) {
    return g[e].steer_record;
  };
  
  template <typename Topology, typename Graph>
  typename boost::disable_if< is_steerable_space<Topology>,
  int >::type& try_get_steer_record(Graph&, typename boost::graph_traits< Graph >::edge_descriptor) {
    throw std::invalid_argument("Required property 'edge_steer_record' on a non-steerable space!");
  };
  
  
  template <typename Topology, bool IsSteerable>
  struct mg_edge_data_base { 
    mg_edge_data_base() { };
  };
  
  template <typename Topology>
  struct mg_edge_data_base<Topology,true> { 
    typedef typename steerable_space_traits<Topology>::steer_record_type steer_record_type;
    steer_record_type steer_record;
    
    mg_edge_data_base() : steer_record() { };
    mg_edge_data_base(const steer_record_type& aRec) : steer_record(aRec) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    mg_edge_data_base(steer_record_type&& aRec) : steer_record(std::move(aRec)) { };
#endif
  };
  
  
};



/**
 * This struct contains the data required on a per-vertex basis for any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct mg_vertex_data {
  /// The position associated to the vertex.
  typename topology_traits<Topology>::point_type position;
};

/**
 * This struct contains the data required on a per-edge basis for any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct mg_edge_data : detail::mg_edge_data_base<Topology, is_steerable_space<Topology>::value > { 
  typedef detail::mg_edge_data_base<Topology, is_steerable_space<Topology>::value > base_type;
  
  mg_edge_data() : base_type() { };
  
  template <typename SteerRec>
  mg_edge_data(const SteerRec& aRec) : base_type(aRec) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  template <typename SteerRec>
  mg_edge_data(SteerRec&& aRec) : base_type(std::move(aRec)) { };
#endif
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const mg_vertex_data<Topology>& vp) {
  using ReaK::to_vect;
  vect_n<double> v_pos = to_vect<double>(vp.position);
  for(std::size_t i = 0; i < v_pos.size(); ++i)
    out << " " << std::setw(10) << v_pos[i];
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any basic path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_motion_graph : public graph::type_erased_graph<Graph> {
  protected:
    typedef any_motion_graph<Topology, Graph> self;
    typedef graph::type_erased_graph<Graph> base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_position")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].position));
      if(aProperty == "edge_steer_record")
        return static_cast<void*>( &( detail::try_get_steer_record<Topology>(*(this->p_graph), boost::any_cast<real_edge_desc>(aElement)) ) );
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_position")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].position);
      if(aProperty == "edge_steer_record")
        return boost::any( detail::try_get_steer_record<Topology>(*(this->p_graph), boost::any_cast<real_edge_desc>(aElement)) );
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};



/**
 * This struct contains the data required on a per-vertex basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct optimal_mg_vertex : mg_vertex_data<Topology> {
  /// The travel-distance accumulated in the vertex, i.e., the travel-distance from the root vertex to this vertex.
  double distance_accum;
  /// The predecessor associated to the vertex, e.g., following the predecessor links starting at the goal node yields a backward trace of the optimal path.
  std::size_t predecessor;
};

/**
 * This struct contains the data required on a per-edge basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct optimal_mg_edge : mg_edge_data<Topology> {
  /// The travel-distance associated to the edge (from source to target).
  double weight;
  
  explicit optimal_mg_edge(double aWeight = 0.0) : mg_edge_data<Topology>(), weight(aWeight) { };
  
  template <typename SteerRec>
  optimal_mg_edge(double aWeight, const SteerRec& aRec) : mg_edge_data<Topology>(aRec), weight(aWeight) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  template <typename SteerRec>
  optimal_mg_edge(double aWeight, SteerRec&& aRec) : mg_edge_data<Topology>(std::move(aRec)), weight(aWeight) { };
#endif
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const optimal_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  vect_n<double> v_pos = to_vect<double>(vp.position);
  for(std::size_t i = 0; i < v_pos.size(); ++i)
    out << " " << std::setw(10) << v_pos[i];
  out << " " << std::setw(10) << vp.distance_accum;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_optimal_motion_graph : public any_motion_graph<Topology, Graph> {
  protected:
    typedef any_optimal_motion_graph<Topology, Graph> self;
    typedef any_motion_graph<Topology, Graph> base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_distance_accum")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].distance_accum));
      if(aProperty == "vertex_predecessor")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].predecessor));
      if(aProperty == "edge_weight")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_edge_desc>(aElement)].weight));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_distance_accum")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].distance_accum);
      if(aProperty == "vertex_predecessor")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].predecessor);
      if(aProperty == "edge_weight")
        return boost::any((*(this->p_graph))[boost::any_cast<real_edge_desc>(aElement)].weight);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_optimal_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};



/**
 * This struct contains the data required on a per-vertex basis for any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct bidir_optimal_mg_vertex : optimal_mg_vertex<Topology> {
  /// The foward travel-distance accumulated in the vertex, i.e., the travel-distance to the goal vertex from this vertex.
  double fwd_distance_accum;
  /// The successor associated to the vertex, e.g., following the successor links starting at the start node yields a trace of the optimal path.
  std::size_t successor;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const bidir_optimal_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  vect_n<double> v_pos = to_vect<double>(vp.position);
  for(std::size_t i = 0; i < v_pos.size(); ++i)
    out << " " << std::setw(10) << v_pos[i];
  out << " " << std::setw(10) << vp.distance_accum;
  out << " " << std::setw(10) << vp.fwd_distance_accum;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any optimal path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_bidir_optimal_motion_graph : public any_optimal_motion_graph<Topology, Graph> {
  protected:
    typedef any_bidir_optimal_motion_graph<Topology, Graph> self;
    typedef any_optimal_motion_graph<Topology, Graph> base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_fwd_distance_accum")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].fwd_distance_accum));
      if(aProperty == "vertex_successor")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].successor));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_fwd_distance_accum")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].fwd_distance_accum);
      if(aProperty == "vertex_successor")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].successor);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_bidir_optimal_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};



/**
 * This struct contains the data required on a per-vertex basis for any A*-like path-planning algorithm (heuristically driven).
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct astar_mg_vertex : optimal_mg_vertex<Topology> {
  /// The heuristic-value associated to the vertex, i.e., the bird-flight distance to the goal.
  double heuristic_value;
  /// The key-value associated to the vertex, computed by the algorithm (usually a combination of accumulated and heuristic distances).
  double key_value;
  /// The color-value associated to the vertex, computed by the algorithm.
  boost::default_color_type astar_color;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const astar_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  vect_n<double> v_pos = to_vect<double>(vp.position);
  for(std::size_t i = 0; i < v_pos.size(); ++i)
    out << " " << std::setw(10) << v_pos[i];
  out << " " << std::setw(10) << vp.distance_accum
      << " " << std::setw(10) << vp.heuristic_value 
      << " " << std::setw(10) << vp.key_value;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any A*-like path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_astar_motion_graph : public any_optimal_motion_graph<Topology, Graph> {
  protected:
    typedef any_astar_motion_graph<Topology, Graph> self;
    typedef any_optimal_motion_graph<Topology, Graph> base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_heuristic_value")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].heuristic_value));
      if(aProperty == "vertex_key_value")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].key_value));
      if(aProperty == "vertex_astar_color")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].astar_color));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_heuristic_value")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].heuristic_value);
      if(aProperty == "vertex_key_value")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].key_value);
      if(aProperty == "vertex_astar_color")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].astar_color);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_astar_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};



/**
 * This struct contains the data required on a per-vertex basis for any A*-like path-planning algorithm (heuristically driven).
 * \tparam Topology The topology type on which the planning is performed.
 */
template <typename Topology>
struct bidir_astar_mg_vertex : bidir_optimal_mg_vertex<Topology> {
  /// The key-value associated to the vertex, computed by the algorithm (usually a combination of accumulated and heuristic distances).
  double key_value;
  /// The color-value associated to the vertex, computed by the algorithm.
  boost::default_color_type astar_color;
};

template <typename Topology>
void print_mg_vertex(std::ostream& out, const bidir_astar_mg_vertex<Topology>& vp) {
  using ReaK::to_vect;
  vect_n<double> v_pos = to_vect<double>(vp.position);
  for(std::size_t i = 0; i < v_pos.size(); ++i)
    out << " " << std::setw(10) << v_pos[i];
  out << " " << std::setw(10) << vp.distance_accum
      << " " << std::setw(10) << vp.fwd_distance_accum 
      << " " << std::setw(10) << vp.key_value;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used by any A*-like path-planning algorithm.
 * \tparam Topology The topology type on which the planning is performed.
 * \tparam Graph The graph type used by the path-planning algorithm.
 */
template <typename Topology, typename Graph>
class any_bidir_astar_motion_graph : public any_bidir_optimal_motion_graph<Topology, Graph> {
  protected:
    typedef any_bidir_astar_motion_graph<Topology, Graph> self;
    typedef any_bidir_optimal_motion_graph<Topology, Graph> base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_key_value")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].key_value));
      if(aProperty == "vertex_astar_color")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].astar_color));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_key_value")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].key_value);
      if(aProperty == "vertex_astar_color")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].astar_color);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_bidir_astar_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};




/**
 * This struct contains the data required on a per-vertex basis for calculating the density
 * from a survey (non-recursive) of the neighborhood.
 * \tparam BaseVertex The type of the underlying vertex.
 */
template <typename BaseVertex>
struct dense_mg_vertex : BaseVertex {
  /// The density-value associated to the vertex, calculated by a survey of the neighborhood.
  double density;
};

template <typename BaseVertex>
void print_mg_vertex(std::ostream& out, const dense_mg_vertex<BaseVertex>& vp) {
  print_mg_vertex(out, static_cast<const BaseVertex&>(vp));
  out << " " << std::setw(10) << vp.density;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used 
 * by a path-planning algorithm that uses a (non-recursive) density metric.
 * \tparam BaseMotionGraph The type-erased motion-graph base-type used.
 */
template <typename BaseMotionGraph>
class any_dense_motion_graph : public BaseMotionGraph {
  protected:
    typedef any_dense_motion_graph<BaseMotionGraph> self;
    typedef BaseMotionGraph base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_density")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].density));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_density")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].density);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_dense_motion_graph(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};


/**
 * This struct contains the data required on a per-vertex basis for calculating the density
 * from a recursive accumulation of the neighborhood statistics (e.g., recursive KL-divergence).
 * \tparam BaseVertex The type of the underlying vertex.
 */
template <typename BaseVertex>
struct recursive_dense_mg_vertex : BaseVertex {
  /// The density associated to the vertex, which represents the probability that a sample drawn from the neighborhood of the vertex will not yield any information gain.
  double density;
  /// Keeps track of the number of neighbors of the vertex.
  std::size_t expansion_trials;
  /// The constriction associated to the vertex, which represents the probability that a sample drawn from the neighborhood of the vertex will be colliding (or unreachable by a collision-free path).
  double constriction;
  /// Keeps track of the number of neighbors of the vertex that could not be connected to it due to a collision.
  std::size_t collision_count;
};

template <typename BaseVertex>
void print_mg_vertex(std::ostream& out, const recursive_dense_mg_vertex<BaseVertex>& vp) {
  print_mg_vertex(out, static_cast<const BaseVertex&>(vp));
  out << " " << std::setw(10) << vp.density
      << " " << std::setw(10) << vp.expansion_trials
      << " " << std::setw(10) << vp.constriction
      << " " << std::setw(10) << vp.collision_count;
};

/**
 * This class template can be used as a type-erased encapsulation of a graph used 
 * by a path-planning algorithm that uses a recursive density metric.
 * \tparam BaseMotionGraph The type-erased motion-graph base-type used.
 */
template <typename BaseMotionGraph>
class any_recursive_dense_mg : public BaseMotionGraph {
  protected:
    typedef any_recursive_dense_mg<BaseMotionGraph> self;
    typedef BaseMotionGraph base_type;
    typedef typename base_type::original_graph_type original_graph_type;
    typedef typename base_type::real_vertex_desc real_vertex_desc;
    typedef typename base_type::real_edge_desc real_edge_desc;
    
    virtual void* get_property_by_ptr(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_density")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].density));
      if(aProperty == "vertex_constriction")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].constriction));
      if(aProperty == "vertex_expansion_trials")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].expansion_trials));
      if(aProperty == "vertex_collision_count")
        return static_cast<void*>(&((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].collision_count));
      
      return base_type::get_property_by_ptr(aProperty, aElement);
    };
    
    virtual boost::any get_property_by_any(const std::string& aProperty, const boost::any& aElement) const {
      
      if(aProperty == "vertex_density")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].density);
      if(aProperty == "vertex_constriction")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].constriction);
      if(aProperty == "vertex_expansion_trials")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].expansion_trials);
      if(aProperty == "vertex_collision_count")
        return boost::any((*(this->p_graph))[boost::any_cast<real_vertex_desc>(aElement)].collision_count);
      
      return base_type::get_property_by_any(aProperty, aElement);
    };
    
  public:
    
    any_recursive_dense_mg(original_graph_type* aPGraph) : base_type(aPGraph) { };
    
};


template <typename Topology, typename Graph>
struct te_mg_selector {
  typedef typename Graph::vertex_bundled VertexProp;
  
  typedef boost::is_convertible< VertexProp*, mg_vertex_data<Topology>* > IsBasicMG;
  typedef boost::is_convertible< VertexProp*, optimal_mg_vertex<Topology>* > IsOptimMG;
  typedef boost::is_convertible< VertexProp*, astar_mg_vertex<Topology>* > IsAStarMG;
  
  typedef 
  typename boost::mpl::if_< IsBasicMG,
    typename boost::mpl::if_< IsOptimMG,
      typename boost::mpl::if_< IsAStarMG,
        astar_mg_vertex<Topology>,
        optimal_mg_vertex<Topology> >::type,
      mg_vertex_data<Topology> >::type,
    void >::type BaseVertexProp;
  
  typedef 
  typename boost::mpl::if_< IsBasicMG,
    typename boost::mpl::if_< IsOptimMG,
      typename boost::mpl::if_< IsAStarMG,
        any_astar_motion_graph<Topology,Graph>,
        any_optimal_motion_graph<Topology,Graph> >::type,
      any_motion_graph<Topology,Graph> >::type,
    graph::type_erased_graph<Graph> >::type BaseMG;
  
  typedef boost::is_convertible< VertexProp*, dense_mg_vertex<BaseVertexProp>* > IsDenseMG;
  typedef boost::is_convertible< VertexProp*, recursive_dense_mg_vertex<BaseVertexProp>* > IsRecDenseMG;
  
  typedef
  typename boost::mpl::if_< IsDenseMG,
    any_dense_motion_graph<BaseMG>,
    typename boost::mpl::if_< IsRecDenseMG,
      any_recursive_dense_mg<BaseMG>,
      BaseMG >::type >::type type;
  
};



/**
 * This stateless functor type can be used to print out the information about an A*-like motion-graph vertex.
 * This is a printing policy type for the vlist_sbmp_report class.
 */
struct mg_vertex_printer : serializable {
  
  /**
   * This call operator prints all the information about a given vertex to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the planning algorithm.
   * \param out The output-stream to which to print the information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  template <typename Vertex, typename Graph>
  void operator()(std::ostream& out, Vertex u, const Graph& g) const {
    print_mg_vertex(out, g[u]);
    out << std::endl;
  };
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };
  
  RK_RTTI_MAKE_ABSTRACT_1BASE(mg_vertex_printer,0xC2460011,1,"mg_vertex_printer",serializable)
};


const std::size_t BASE_MOTION_GRAPH_KIND_MASK         = 0x0F;
const std::size_t BASIC_MOTION_GRAPH_KIND             = 0x00;
const std::size_t OPTIMAL_MOTION_GRAPH_KIND           = 0x01;
const std::size_t ASTAR_MOTION_GRAPH_KIND             = 0x03;
const std::size_t BIDIR_MOTION_GRAPH_KIND             = 0x04;
const std::size_t BIDIR_OPTIMAL_MOTION_GRAPH_KIND     = 0x05;
const std::size_t BIDIR_ASTAR_MOTION_GRAPH_KIND       = 0x07;

const std::size_t DENSITY_MOTION_GRAPH_KIND_MASK      = 0xF0;
const std::size_t DENSE_MOTION_GRAPH_KIND             = 0x10;
const std::size_t RECURSIVE_DENSE_MOTION_GRAPH_KIND   = 0x20;


/**
 * This stateless functor type can be used to print out the information about an A*-like motion-graph vertex.
 * This is a printing policy type for the vlist_sbmp_report class.
 */
template <typename Topology>
struct any_mg_vertex_printer : serializable {
  typedef any_mg_vertex_printer<Topology> self;
  typedef typename topology_traits<Topology>::point_type PointType;
  
  std::size_t graph_kind;
  
  /**
   * This call operator prints all the information about a given vertex to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the planning algorithm.
   * \param out The output-stream to which to print the information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  void operator()(std::ostream& out, graph::any_graph::vertex_descriptor u, const graph::any_graph& g) const {
    using ReaK::to_vect;
    using ReaK::graph::get_dyn_prop;
    
    vect_n<double> v_pos = to_vect<double>( get_dyn_prop<const PointType&>("vertex_position", u, g) );
    for(std::size_t i = 0; i < v_pos.size(); ++i)
      out << " " << std::setw(10) << v_pos[i];
    
    if( graph_kind & OPTIMAL_MOTION_GRAPH_KIND ) {
      out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_distance_accum", u, g);
      if( graph_kind & BIDIR_MOTION_GRAPH_KIND ) {
        out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_fwd_distance_accum", u, g);
      };
    };
    
    if( graph_kind & ASTAR_MOTION_GRAPH_KIND ) {
      if( graph_kind & BIDIR_MOTION_GRAPH_KIND ) {
        out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_key_value", u, g);
      } else {
        out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_heuristic_value", u, g) 
            << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_key_value", u, g);
      };
    };
    
    if( graph_kind & DENSE_MOTION_GRAPH_KIND ) {
      out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_density", u, g);
    } else if( graph_kind & RECURSIVE_DENSE_MOTION_GRAPH_KIND ) {
      out << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_density", u, g)
          << " " << std::setw(10) << get_dyn_prop<const std::size_t&>("vertex_expansion_trials", u, g)
          << " " << std::setw(10) << get_dyn_prop<const double&>("vertex_constriction", u, g)
          << " " << std::setw(10) << get_dyn_prop<const std::size_t&>("vertex_collision_count", u, g);
    };
    
    out << std::endl;
  };
  
  explicit any_mg_vertex_printer(std::size_t aGraphKind = 0) : graph_kind(aGraphKind) { };
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { 
    A & RK_SERIAL_SAVE_WITH_NAME(graph_kind);
  };
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
    A & RK_SERIAL_LOAD_WITH_NAME(graph_kind);
  };
  
  RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460012,1,"any_mg_vertex_printer",serializable)
};



};

};


#endif



