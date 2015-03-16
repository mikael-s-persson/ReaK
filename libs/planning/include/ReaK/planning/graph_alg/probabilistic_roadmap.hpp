/**
 * \file probabilistic_roadmap.hpp
 *
 * This library provides function templates and concepts that implement the probabilistic roadmap (PRM) 
 * algorithm (as of "Geraerts and Overmars, 2002"). A PRM is generated in two phases. First, a number of 
 * vertices are generated at random from the free-space (configuration-space which is not occupied by 
 * obstacles) and bi-directional (or undirected) connections between those vertices are established, as
 * much as possible (if they can be directly connected through free-space). Second, difficult areas 
 * of the roadmap are determined using some density function (heuristic) and a priority-queue of vertices
 * with high density (whatever that means for the heuristic chosen). Then, high-density vertices are 
 * selected for expansion, in which case a number of neighboring vertices are generated coming out of 
 * the selected vertex into free-space. The first phase is called the construction phase, and the 
 * second is called the expansion phase. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, how to measure density, how to compare 
 * densities, when to stop the algorithm, what proportion of constructed vs. expanded vertices should 
 * be generated, etc. All these customization points are left to the user to implement, some are 
 * defined by the PRMVisitorConcept (expand-vertex and compute-density) while others are provided 
 * as functors to the function template that generates the PRM (generate_prm).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_PROBABILISTIC_ROADMAP_HPP
#define REAK_PROBABILISTIC_ROADMAP_HPP

#include <utility>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>
#include <boost/property_map/property_map.hpp>

#include <ReaK/core/base/global_rng.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/random_sampler_concept.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>

#include "prm_connector.hpp"
#include "sbmp_visitor_concepts.hpp"
#include <set>


namespace ReaK {

namespace graph {


namespace detail { namespace {
  
  
  
  template <typename PRMVisitor,
            typename UpdatableQueue, 
            typename IndexInHeapMap,
            typename PositionMap, 
            typename CCRootMap>
  struct prm_conn_visitor
  {
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::property_traits<CCRootMap>::value_type CCRootValue;
    
    prm_conn_visitor(PRMVisitor vis, UpdatableQueue& Q, IndexInHeapMap index_in_heap, 
                     PositionMap pos, CCRootMap cc_root, std::set<CCRootValue>& cc_set) : 
                     m_vis(vis), m_Q(Q), m_index_in_heap(index_in_heap), 
                     m_position(pos), m_cc_root(cc_root), m_cc_set(cc_set) { };
    
    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor create_vertex(const PositionValue& p, Graph& g) const {
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename Graph::vertex_bundled VertexProp;
      
      VertexProp up;
      put(m_position, up, p);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      put(m_cc_root, u, u);
      m_cc_set.insert(u);
      m_vis.vertex_added(u,g);
      put(m_index_in_heap, u, static_cast<std::size_t>(-1));
      
      return u;
    };
    
    template <typename Vertex, typename Graph>
    void shortcut_cc_root(Vertex u, Graph& g) const {
      std::stack<Vertex> u_trace;
      u_trace.push(u);
      while( get(m_cc_root, u) != u ) {
        u = get(m_cc_root, u);
        u_trace.push(u);
      };
      while(!u_trace.empty()) {
        put(m_cc_root, u_trace.top(), u);
        u_trace.pop();
      };
    };
    
    template <typename Edge, typename Graph>
    void edge_added(Edge e, Graph& g) const { 
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      
      m_vis.edge_added(e,g); 
      
      Vertex u = source(e, g);
      Vertex v = target(e, g);
      shortcut_cc_root(u, g);
      shortcut_cc_root(v, g);
      if(get(m_cc_root, v) != get(m_cc_root, u)) {
        Vertex r1 = get(m_cc_root, u);
        Vertex r2 = get(m_cc_root, v);
        put(m_cc_root, r2, r1);
        put(m_cc_root, v, r1);
        m_cc_set.erase(r2);
        
        if(m_cc_set.size() < 2)
          m_vis.publish_path(g);
      };
    };
    
    template <typename Vertex, typename Graph>
    void travel_explored(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_explored(u, v, g);
    };
    
    template <typename Vertex, typename Graph>
    void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_succeeded(u, v, g);
    };
    
    template <typename Vertex, typename Graph>
    void travel_failed(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_failed(u, v, g);
    };
    
    template <typename Vertex, typename Graph>
    void requeue_vertex(Vertex u, Graph& g) const {
      m_vis.affected_vertex(u, g);
      m_Q.push_or_update(u);
    };
    template <typename Vertex, typename Graph>
    void affected_vertex(Vertex u, Graph& g) const { requeue_vertex(u,g); }; // same function, different name.
    
    bool keep_going() const { return m_vis.keep_going(); };
    
    template <typename Vertex, typename Graph>
    std::pair< bool, typename Graph::edge_bundled > can_be_connected(Vertex u, Vertex v, Graph& g) const {
      return m_vis.can_be_connected(u, v, g);
    };
    
    template <typename Vertex, typename Graph>
    boost::tuple< PositionValue, bool, typename Graph::edge_bundled > random_walk(Vertex u, Graph& g) const {
      return m_vis.random_walk(u, g);
    };
    
    bool is_position_free(const PositionValue& p) const {
      return m_vis.is_position_free(p);
    };
    
    PRMVisitor m_vis;
    UpdatableQueue& m_Q; 
    IndexInHeapMap m_index_in_heap;
    
    PositionMap m_position;
    CCRootMap m_cc_root;
    
    std::set<CCRootValue>& m_cc_set;
  };
  
  
  
  
  
  template <typename Graph,
            typename Topology,
            typename PRMConnVisitor,
            typename PositionMap,
            typename RandomSampler,
            typename MutableQueue,
            typename NodeConnector,
            typename NcSelector>
  inline void generate_prm_impl(Graph& g,
                                const Topology& super_space,
                                PRMConnVisitor& vis,
                                PositionMap position,
                                RandomSampler get_sample,
                                MutableQueue Q,
                                NodeConnector connect_vertex,
                                const NcSelector& select_neighborhood,
                                double expand_probability) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    
    while(vis.keep_going()) {
      
      double rand_value = boost::uniform_01<ReaK::global_rng_type&,double>(ReaK::get_global_rng())(); // generate random-number between 0 and 1.
      
      if(rand_value > expand_probability) {
        //Construction node:
        EdgeProp ep;
        PositionValue p_rnd = get_sample(super_space);
        while(!vis.is_position_free(p_rnd))
          p_rnd = get_sample(super_space);
        
        connect_vertex(p_rnd, boost::graph_traits<Graph>::null_vertex(), ep, g, super_space, vis, position, select_neighborhood);
        
      } else {
        //Expansion node:
        //use the priority queue to get the vertices that need expansion.
        Vertex v = Q.top();
        PositionValue p_rnd; bool expanding_worked; EdgeProp ep;
        boost::tie(p_rnd, expanding_worked, ep) = vis.random_walk(v, g);
        
        if(expanding_worked)
          connect_vertex(p_rnd, v, ep, g, super_space, vis, position, select_neighborhood);
        else
          Q.pop(); // if one cannot expand from this vertex, then leave it out of future attempts.
      };
      
    };

  };
  
  
}; }; // detail


/**
 * This function is the basic PRM algorithm (as of "Geraerts and Overmars, 2002"). Note that 
 * all the design decisions (how to generate vertices, how to select the neighborhood, how to 
 * check the collision) have been externalized via the provided free_space topology, visitor 
 * object, and neighborhood selection functor.
 * 
 * \tparam Graph A mutable graph type that can store the roadmap, should model boost::MutableGraphConcept 
 *         and boost::VertexListGraphConcept (either bidirectional or undirected graph, the algorithm 
 *         will deal with either cases as appropriate).
 * \tparam Topology A topology type on which the vertex positions lie, should model the TopologyConcept,
 *         the MetricSpaceConcept (has an attached distance-metric) and the PointDistributionConcept (has 
 *         an attached random-sampler functor).
 * \tparam PRMVisitor A PRM visitor type, should model the PRMVisitorConcept.
 * \tparam PositionMap A property-map type that can store the position of each vertex. 
 * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
 * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
 *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
 *         See classes in the topological_search.hpp header-file.
 * \tparam CompareFunction A functor type that can be used to compare density values (strict weak-ordering).
 * \param g A mutable graph that should initially store the starting 
 *        vertex (if not it will be randomly generated) and will store 
 *        the generated graph once the algorithm has finished.
 * \param super_space A topology (as defined by ReaK) that represents the configuration-space
 *        used for the planning (a space with no obstacles).
 * \param vis A PRM visitor implementing the PRMVisitorConcept. This is the 
 *        main point of customization and recording of results that the 
 *        user can implement.
 * \param position A mapping that implements the MutablePropertyMap Concept. Also,
 *        the value_type of this map should be the same type as the topology's 
 *        value_type.
 * \param get_sample A random sampler of positions in the super-space (obstacle-free topology).
 * \param density A property-map that provides the density values assiciated to each vertex.
 * \param select_neighborhood A callable object (functor) that can select a list of 
 *        vertices of the graph that ought to be connected to a new 
 *        vertex. The list should be sorted in order of increasing "distance".
 * \param expand_probability The probability (between 0 and 1) that an expansion will be 
 *        performed as opposed to a general sampling from the free-space.
 */
template <typename Graph,
          typename Topology,
          typename PRMVisitor,
          typename PositionMap,
          typename RandomSampler,
          typename DensityMap,
          typename NcSelector>
inline void generate_prm(Graph& g,
                         const Topology& super_space,
                         PRMVisitor vis,
                         PositionMap position,
                         RandomSampler get_sample,
                         DensityMap density,
                         const NcSelector& select_neighborhood,
                         double expand_probability) {
  BOOST_CONCEPT_ASSERT((PRMVisitorConcept<PRMVisitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>)); // for the distance-metric.
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
  
  typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename Graph::vertex_bundled VertexProp;
  
  typedef boost::vector_property_map<Vertex> CCRootMap;
  CCRootMap cc_root;
  std::set<Vertex> cc_set;
  
  typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
  IndexInHeapMap index_in_heap;
  
  typedef boost::composite_property_map< DensityMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphDensityMap;
  GraphDensityMap g_density = GraphDensityMap(density, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
  
  typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, GraphDensityMap, std::less<double> > MutableQueue;
  MutableQueue Q(g_density, index_in_heap, std::less<double>()); //priority queue holding the "expandable" vertices.
  
  if(num_vertices(g) == 0) {
    PositionValue p = get_sample(super_space);
    while(!vis.is_position_free(p))
      p = get_sample(super_space);
    VertexProp up;
    put(position, up, p);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    Vertex u = add_vertex(std::move(up), g);
#else
    Vertex u = add_vertex(up, g);
#endif
    put(cc_root, u, u);
    vis.vertex_added(u, g);
    cc_set.insert(u);
  } else {
    
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      vis.affected_vertex(*ui,g);
      Q.push(*ui);
      put(cc_root, *ui, *ui);
    };
    
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      if( get(cc_root, *ui) != *ui )
        continue;
      cc_set.insert(*ui);
      std::stack<Vertex> u_trace;
      u_trace.push(*ui);
      while(!u_trace.empty()) {
        Vertex v = u_trace.top(); u_trace.pop();
        if(get(cc_root, v) == *ui)
          continue;
        put(cc_root, v, *ui);
        typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = out_edges(v,g); ei != ei_end; ++ei) {
          if(target(*ei, g) == v)
            u_trace.push(source(*ei, g));
          else
            u_trace.push(target(*ei, g));
        };
      };
    };
    
  };
  
  detail::prm_conn_visitor<PRMVisitor, MutableQueue, IndexInHeapMap, PositionMap, CCRootMap>
    prm_conn_vis(vis, Q, index_in_heap, position, cc_root, cc_set);
  
  detail::generate_prm_impl( g, super_space, prm_conn_vis, position, get_sample, Q, prm_node_connector(), select_neighborhood, expand_probability);
  
};





};


};

#endif










