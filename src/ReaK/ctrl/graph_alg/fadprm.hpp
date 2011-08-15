/**
 * \file fadprm.hpp
 *
 * This library provides function templates and concepts that implement the flexible anytime dynamic 
 * probabilistic roadmap (FADPRM) algorithm (as of "Belghith et al., 2006"). A FADPRM is uses the  
 * AD* search algorithm to drive the expansion of a probabilistic roadmap into the free-space 
 * in order to connect a start and goal location. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, how to measure density, how to compare 
 * densities, when to stop the algorithm, etc. All these customization points are left to the user 
 * to implement, some are defined by the FADPRMVisitorConcept (expand-vertex and compute-density) 
 * while others are provided as functors to the function template that generates the FADPRM (generate_fadprm).
 *
 * The FADPRM algorithm is so closely related to the AD* algorithm that it is actually implemented 
 * using the AD* loop (a function in the detail namespace of the AD* functions) and a customized 
 * AD* BFS visitor that wraps the FADPRM customizations along with the FADPRM visitor for the 
 * user's customizations.
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

#ifndef FADPRM_HPP
#define FADPRM_HPP

#include <functional>
#include <boost/utility/enable_if.hpp>

#include "adstar_search.hpp"
#include "probabilistic_roadmap.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {

  /**
   * This concept class defines the valid expressions required of a class to be used as a visitor 
   * class for the FADPRM algorithm. A visitor class is essentially a class that regroups a number of 
   * callback functions that can be used to inject customization into the FADPRM algorithm. In other 
   * words, the visitor pattern in generic programming is an implementation of IoC 
   * (Inversion of Control), since the FADPRM algorithm is in control of execution, but custom behavior can
   * be injected in several places, even blocking the algorithm if needed.
   * 
   * Required concepts:
   * 
   * The visitor class should model the boost::CopyConstructibleConcept.
   * 
   * The visitor class should model the PRMVisitorConcept.
   * 
   * The visitor class should model the ADStarVisitorConcept.
   * 
   * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   */
  template <typename Visitor, typename Graph>
  struct FADPRMVisitorConcept {
    void constraints()
    {
      boost::function_requires< boost::CopyConstructibleConcept<Visitor> >();
      boost::function_requires< PRMVisitorConcept<Visitor,Graph> >();
      boost::function_requires< ADStarVisitorConcept<Visitor,Graph> >();
    }
  };
  
  /**
   * This class is simply a "null" visitor for the FADPRM algorithm. It is null in the sense that it
   * will do nothing on all accounts.
   * \tparam Topology The topology type that represents the free-space.
   * \tparam PositionMap The property-map type which can store the position associated to each vertex.
   */
  template <typename Topology, typename PositionMap>
  class default_fadprm_visitor {
    public:
      default_fadprm_visitor(const Topology& free_space, PositionMap position) : m_free_space(free_space), m_position(position) {};
      default_fadprm_visitor(const default_fadprm_visitor<Topology,PositionMap>& aVis) : m_free_space(aVis.m_free_space), m_position(aVis.m_position) { };

      template <class Vertex, class Graph>
      void initialize_vertex(Vertex u, Graph& g) { };
      template <class Vertex, class Graph>
      void discover_vertex(Vertex u, Graph& g) { };
      template <class Vertex, class Graph>
      void inconsistent_vertex(Vertex u, Graph& g) { };
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) { };
      template <class Edge, class Graph>
      void examine_edge(Edge e, Graph& g) { };
      template <class Edge, class Graph>
      void edge_relaxed(Edge e, Graph& g) { };
      template <class Vertex, class Graph>
      void forget_vertex(Vertex u, Graph& g) { };
      template <class Vertex, class Graph>
      void finish_vertex(Vertex u, Graph& g) { };
      template <class Vertex, class Graph>
      void recycle_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void vertex_added(Vertex,Graph&) { };
      template <typename Edge, typename Graph>
      void edge_added(Edge,Graph&) { };
      template <typename Vertex, typename Graph>
      void expand_vertex(Vertex,Graph&,std::vector<Vertex>&) { };
      template <typename Vertex, typename Graph>
      void update_density(Vertex, Graph& g) { };

      const Topology& m_free_space;
      PositionMap m_position;
  };

  /**
   * This class is a composite visitor class template. It can be used to glue together two classes
   * which model the ADStarVisitorConcept and PRMVisitorConcept, respectively. Note that it 
   * is preferred to use the make_composite_fadprm_visitor function template to let the 
   * compiler deduce the type instead of providing the template arguments manually.
   * \tparam ADStarVisitor The visitor type which models the ADStarVisitorConcept.
   * \tparam PRMVisitor The visitor type which models the PRMVisitorConcept.
   */
  template <typename ADStarVisitor, typename PRMVisitor>
  struct composite_fadprm_visitor : public ADStarVisitor, public PRMVisitor {
    composite_fadprm_visitor(ADStarVisitor aVis_adstar, PRMVisitor aVis_prm) : ADStarVisitor(aVis_adstar), PRMVisitor(aVis_prm) { };
    composite_fadprm_visitor(const composite_fadprm_visitor<ADStarVisitor,PRMVisitor>& aVis) : ADStarVisitor(aVis), PRMVisitor(aVis) { };
  };
  
  
  /**
   * This function template creates a composite visitor class. It can be used to glue together two classes
   * which model the ADStarVisitorConcept and PRMVisitorConcept, respectively.
   * \tparam ADStarVisitor The visitor type which models the ADStarVisitorConcept.
   * \tparam PRMVisitor The visitor type which models the PRMVisitorConcept.
   * \param aVis_adstar The visitor object for the AD* functions.
   * \param aVis_prm The visitor object for the PRM functions.
   */
  template <typename ADStarVisitor, typename PRMVisitor>
  inline composite_fadprm_visitor<ADStarVisitor,PRMVisitor> make_composite_fadprm_visitor(ADStarVisitor aVis_adstar, PRMVisitor aVis_prm) {
    return composite_fadprm_visitor<ADStarVisitor,PRMVisitor>(aVis_adstar,aVis_prm);
  };
  
  
  
  
  
  
  
  namespace detail {
  
    template <typename Topology,
              typename AStarHeuristicMap, 
              typename UniformCostVisitor,
              typename UpdatableQueue, 
	      typename List,
	      typename IndexInHeapMap,
	      typename PredecessorMap,
	      typename KeyMap, 
	      typename DistanceMap, 
	      typename RHSMap, 
	      typename WeightMap,
	      typename DensityMap, 
	      typename PositionMap, 
	      typename NcSelector,
	      typename ColorMap, 
	      typename ScalarType,
	      typename CompareFunction,
	      typename EqualCompareFunction, 
	      typename CombineFunction,
	      typename ComposeFunction>
    struct fadprm_bfs_visitor
    {

      typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
      typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
      typedef boost::color_traits<ColorValue> Color;
      typedef typename boost::property_traits<DistanceMap>::value_type distance_type;
      typedef typename boost::property_traits<WeightMap>::value_type weight_type;
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;

      fadprm_bfs_visitor(const Topology& free_space, AStarHeuristicMap h, UniformCostVisitor vis,
                         UpdatableQueue& Q, List& I, IndexInHeapMap index_in_heap, PredecessorMap p,
                         KeyMap k, DistanceMap d, RHSMap rhs, WeightMap w,
                         DensityMap dens, PositionMap pos, NcSelector select_neighborhood, ColorMap col, const ScalarType& beta,
                         CompareFunction compare, EqualCompareFunction equal_compare,
			 CombineFunction combine, ComposeFunction compose,
			 distance_type inf, distance_type zero)
        : m_free_space(free_space), m_h(h), m_vis(vis), m_Q(Q), m_I(I), 
          m_index_in_heap(index_in_heap), m_predecessor(p), m_key(k),
          m_distance(d), m_rhs(rhs), m_weight(w), m_density(dens), m_position(pos), 
          m_select_neighborhood(select_neighborhood), m_color(col), m_beta(beta),
          m_compare(compare), m_equal_compare(equal_compare), m_combine(combine),
          m_compose(compose), m_inf(inf), m_zero(zero) {};

      fadprm_bfs_visitor(fadprm_bfs_visitor<Topology, AStarHeuristicMap, UniformCostVisitor,
                                            UpdatableQueue, List, IndexInHeapMap, PredecessorMap,
					    KeyMap, DistanceMap, RHSMap, WeightMap,
					    DensityMap, PositionMap, NcSelector, 
					    ColorMap, ScalarType, CompareFunction,
					    EqualCompareFunction, CombineFunction,
					    ComposeFunction>& aVis)
        : m_free_space(aVis.m_free_space), m_h(aVis.m_h), m_vis(aVis.m_vis), m_Q(aVis.m_Q), m_I(aVis.m_I),
          m_index_in_heap(aVis.m_index_in_heap), m_predecessor(aVis.m_predecessor),
          m_key(aVis.m_key), m_distance(aVis.m_distance), m_rhs(aVis.m_rhs), m_weight(aVis.m_weight),
          m_density(aVis.m_density), m_position(aVis.m_position), m_select_neighborhood(aVis.m_select_neighborhood), m_color(aVis.m_color), m_beta(aVis.m_beta), m_compare(aVis.m_compare),
          m_equal_compare(aVis.m_equal_compare), m_combine(aVis.m_combine),
          m_compose(aVis.m_compose), m_inf(aVis.m_inf), m_zero(aVis.m_zero) {};

      template <class Vertex, class Graph>
      void initialize_vertex(Vertex u, Graph& g) {
        m_vis.initialize_vertex(u, g);
      };
      template <class Vertex, class Graph>
      void discover_vertex(Vertex u, Graph& g) {
        m_vis.discover_vertex(u, g);
      };
      template <class Vertex, class Graph>
      void inconsistent_vertex(Vertex u, Graph& g) {
	m_vis.inconsistent_vertex(u, g);
      };
      template <class Vertex, class Graph>
      typename boost::enable_if< boost::is_undirected_graph<Graph> >::type connect_vertex(Vertex u, Graph& g) {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	
	PositionValue p = get(m_position, u); 
        std::vector<Vertex> Nc;
        m_select_neighborhood(p,Nc,g,m_free_space,m_position); 
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
	  if((u != *it) && (m_free_space.distance(boost::get(m_position,*it),p) != std::numeric_limits<distance_type>::infinity())) {
	    //this means that u is reachable from *it.
	    std::pair<Edge, bool> ep = boost::add_edge(*it,u,g); 
	    if(ep.second) {
	      m_vis.edge_added(ep.first, g); 
	      update_key(*it,g); 
	      update_vertex(*it,g);
	      //m_Q.update(*it);
	    };
	  };
        }; 
	put(m_index_in_heap, u, (std::size_t)(-1));
	update_key(u,g); 
	//m_Q.push(u); put(m_color, u, Color::gray());
	update_vertex(u,g);
      };
      template <class Vertex, class Graph>
      typename boost::enable_if< boost::is_directed_graph<Graph> >::type connect_vertex(Vertex u, Graph& g) {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	
	PositionValue p = get(m_position, u); 
        std::vector<Vertex> Pred;
	std::vector<Vertex> Succ;
        m_select_neighborhood(p,Pred,Succ,g,m_free_space,m_position); 
        for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
	  if((u != *it) && (m_free_space.distance(boost::get(m_position,*it),p) != std::numeric_limits<distance_type>::infinity())) {
	    //this means that u is reachable from *it.
	    std::pair<Edge, bool> ep = boost::add_edge(*it,u,g); 
	    if(ep.second) {
	      m_vis.edge_added(ep.first, g); 
	      update_key(*it,g); 
	      update_vertex(*it,g);
	      //m_Q.update(*it);
	    };
	  };
        };
	put(m_index_in_heap, u, (std::size_t)(-1));
	update_key(u,g); 
	//m_Q.push(u); put(m_color, u, Color::gray());
	update_vertex(u,g);
	for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
	  if((u != *it) && (m_free_space.distance(p,boost::get(m_position,*it)) != std::numeric_limits<distance_type>::infinity())) {
	    //this means that u is reachable from *it.
	    std::pair<Edge, bool> ep = boost::add_edge(u,*it,g); 
	    if(ep.second) {
	      m_vis.edge_added(ep.first, g); 
	      update_key(*it,g); 
	      update_vertex(*it,g);
	      //m_Q.update(*it);
	    };
	  };
        }; 
      };
      
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        
	m_vis.examine_vertex(u, g);
	
	std::vector<Vertex> v_list;
	m_vis.expand_vertex(u, g, v_list);
	
	for(typename std::vector<Vertex>::iterator vi = v_list.begin(); vi != v_list.end(); ++vi) {
	  Vertex v = *vi;
	
          m_vis.vertex_added(v,g); 
	  put(m_color, v, Color::white());
          put(m_distance, v, m_inf);
          put(m_rhs, v, m_inf);
          put(m_key, v, KeyValue(m_inf,m_inf,m_compare,m_equal_compare));
          put(m_predecessor, v, v);
          
	  connect_vertex(v,g);
	};
	update_key(u,g); 
	//m_Q.update(u);
	
      };
      template <class Edge, class Graph>
      void examine_edge(Edge e, Graph& g) {
        if (m_compare(get(m_weight, e), m_zero))
          throw boost::negative_edge();
        m_vis.examine_edge(e, g);
      };
      template <class Vertex, class Graph>
      void forget_vertex(Vertex u, Graph& g) {
        m_vis.forget_vertex(u, g);
      };
      template <class Vertex, class Graph>
      void finish_vertex(Vertex u, Graph& g) {
	m_vis.finish_vertex(u, g);
      };
      template <class Vertex, class Graph>
      void recycle_vertex(Vertex u, Graph& g) {
	m_vis.recycle_vertex(u, g);
      };

      template <class Vertex, typename Graph>
      void update_key(Vertex u, Graph& g) {
	m_vis.update_density(u, g);
	distance_type g_u = get(m_distance, u);
	distance_type rhs_u = get(m_rhs, u);
	distance_type h_u = get(m_h, u);
	if( m_compare( rhs_u, g_u ) ) {
	  distance_type f_u = m_combine( rhs_u, h_u ); //m_compose(0.5, m_combine( rhs_u, h_u ));
	  put(m_key, u, KeyValue( m_combine(m_compose(1.0 - m_beta, get(m_density, u)), m_compose(m_beta, f_u)), rhs_u, m_compare, m_equal_compare));
	} else {
	  distance_type f_u = m_combine( g_u, h_u ); // m_compose(0.5, m_combine( g_u, h_u ));
	  put(m_key, u, KeyValue( f_u, rhs_u, m_compare, m_equal_compare));
	};
      };

      template <class Vertex, class BidirectionalGraph>
      void update_vertex(Vertex u, BidirectionalGraph& g) {
	boost::function_requires< boost::BidirectionalGraphConcept<BidirectionalGraph> >();
	typedef boost::graph_traits<BidirectionalGraph> GTraits;
	typename GTraits::in_edge_iterator ei, ei_end;
        
	ColorValue col_u = get(m_color, u);
        
	if(col_u == Color::white()) {
	  put(m_distance, u, m_inf);
	  col_u = Color::green();
	  put(m_color, u, col_u);
	};
         
	distance_type g_u = get(m_distance, u); 
	distance_type rhs_u = get(m_rhs, u);

	if(!m_equal_compare(rhs_u,m_zero)) {
	  rhs_u = m_inf; 
	  Vertex pred_u = get(m_predecessor, u);
	  typename GTraits::edge_descriptor pred_e; 
	  for(boost::tie(ei,ei_end) = boost::in_edges(u,g); ei != ei_end; ++ei) {
	    distance_type rhs_tmp = m_combine(get(m_weight, *ei), get(m_distance, boost::source(*ei,g))); 
	    if(m_compare(rhs_tmp, rhs_u)) {
	      rhs_u = rhs_tmp; 
	      put(m_rhs, u, rhs_u);
	      put(m_predecessor, u, boost::source(*ei,g)); 
	      pred_e = *ei;
	    };
	  };
	  rhs_u = get(m_rhs, u); 
	  if(pred_u != get(m_predecessor, u))                          m_vis.edge_relaxed(pred_e, g);
	};

	if(!m_equal_compare(rhs_u,g_u)) { 
	  if((col_u != Color::black()) && (col_u != Color::red())) {
	    update_key(u,g); 
	    m_Q.push_or_update(u); 
	    put(m_color, u, Color::gray());                            m_vis.discover_vertex(u, g);
	  } else if(col_u == Color::black()) {
	    m_I.push_back(u); 
	    put(m_color, u, Color::red());                             m_vis.inconsistent_vertex(u,g);
	  };
	} else if(m_Q.contains(u)) {
	  put(m_key, u, KeyValue(m_zero,m_zero,m_compare,m_equal_compare)); 
	  m_Q.update(u); 
	  m_Q.pop(); //remove from OPEN set
	  update_key(u,g); 
	  put(m_color, u, Color::green());                             m_vis.forget_vertex(u, g);
	}; 
      };

      const Topology& m_free_space;
      AStarHeuristicMap m_h;
      UniformCostVisitor m_vis;
      UpdatableQueue& m_Q;
      List& m_I;
      IndexInHeapMap m_index_in_heap;
      PredecessorMap m_predecessor;
      KeyMap m_key;
      DistanceMap m_distance;
      RHSMap m_rhs;
      WeightMap m_weight;
      DensityMap m_density;
      PositionMap m_position;
      NcSelector m_select_neighborhood;
      ColorMap m_color;
      const ScalarType& m_beta;
      CompareFunction m_compare;
      EqualCompareFunction m_equal_compare;
      CombineFunction m_combine;
      ComposeFunction m_compose;
      distance_type m_inf;
      distance_type m_zero;
    };

  
  
  }; //end of detail namespace.
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the FADPRM algorithm (as of Belghith et al., 2006), without initialization of the 
   * existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values 
   *         for each vertex in the graph.
   * \tparam FADPRMVisitor The type of the FADPRM visitor to be used, should model the FADPRMVisitorConcept.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of 
   *         each vertex to the goal (internal use to AD*).
   * \tparam KeyMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam PositionMap A property-map type that can store the position of each vertex. 
   * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
   *         See classes in the topological_search.hpp header-file.
   * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors 
   *         are used to mark vertices by their status in the AD* algorithm (white = not visited, 
   *         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED, 
   *         not OPEN), red = inconsistent (in INCONS)).
   * \tparam AdjustFunction This functor type represents the callback function that will adjust the 
   *         epsilon value (sloppyness of the A* search).
   * \tparam RunningPredicate This function type represents the callback function that will evaluate 
   *         if it is still worth continuing the AD* search passes.
   * \tparam GetStartNodeFunction This function type represents the callback function that will give 
   *         the starting point for the search (usually the goal vertex).
   * \tparam ChangeEventFunction This function type represents the callback that can wait for a change 
   *         in edge weights, and provide the list of changed edges.
   * \tparam PublishPathFunction This function type represents the callback that is called when a 
   *         path has been obtained (so that it can be buffered and executed).
   * \tparam ScalarType The type of the scalar value epsilon (amplification factor, controls the 
   *         anytime nature of this algorithm).
   * \tparam CompareFunction A binary comparison functor type that returns true if the first operand 
   *         is strictly better (less-than) than the second operand.
   * \tparam EqualCompareFunction A binary comparison functor type that returns true if both operands 
   *         are equal to each other.
   * \tparam CombineFunction A binary combination functor type that returns the sum of its operands, 
   *         semantically-speaking.
   * \tparam ComposeFunction A binary composition functor type that amplifies a heuristic distance 
   *         metric by a scalar value (i.e. epsilon x h(u) ).
   * 
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated graph once the algorithm has finished.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *        that it is required to generate only random points in 
   *        the free-space and to only allow interpolation within the free-space.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param vis A FADPRM visitor implementing the FADPRMVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param predecessor The property-map which will store the resulting path by connecting 
   *        vertices together with their optimal predecessor (follow in reverse to discover the 
   *        complete path).
   * \param distance The property-map which stores the estimated distance of each vertex to the goal.
   * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the 
   *        goal (for internal use).
   * \param key The property-map which stores the AD* key-values associated to each vertex.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
   * \param density A property-map that provides the density values assiciated to each vertex.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param color The property-map which stores the color-value of the vertices, colors are used to mark
   *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), 
   *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
   * \param adj_epsilon The callback functor that will adjust the epsilon value (sloppyness of the sub-optimal
   *        A* search) depending (usually) on the level of change in the weights that was recorded after a callback
   *        using the edge_change functor. Used as: epsilon = adj_epsilon(epsilon, max_w_change); (takes the current 
   *        epsilon value and the maximum recorder weight change, and outputs the new epsilon value).
   * \param keep_going The callback functor that will tell whether the AD* search should keep on going (running 
   *        predicate). Should be callable with no parameter and output a bool value (true to keep going, false
   *        to stop).
   * \param get_start The functor that is called to obtain the starting point for the AD* search (usually the 
   *        goal vertex). Should be callable with no parameter and output a vertex-descriptor of the graph (i.e. 
   *        a vertex of the graph g).
   * \param edge_change The callback functor that is called to make a check whether the weights of the graph
   *        have changed do to a dynamic environment. Used as: max_w_change = edge_change(affected_edges); 
   *        (takes an STL container (std::vector) of affected edges and populates it all the edges that 
   *        should be updated because their weight has changed, also returns the maximum weight change 
   *        measured (or any other measure of how bad the change in the environment is, this value will 
   *        be forwarded unchanged to the adj_epsilon functor)). The container of affected edges is subject 
   *        to change, so it is recommended to define this functor with a templated call-operator that 
   *        can take a generic type of STL container, by non-const reference.
   * \param publish_path The callback functor that is called when the AD* search has ended, and presumably,
   *        found a path. Used as: publish_path(); (takes no parameters, returns nothing). The path is retrieved
   *        by walking the predecessors of the vertices (using the given predecessor property-map).
   * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime 
   *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
   * \param inf The quantity that represents infinity (either a very large value or the infinity value for 
   *        the underlying value-type).
   * \param zero The quantity that represents zero with the given value-type.
   * \param compare A binary comparison functor that returns true if the first operand is strictly 
   *        better (less-than) than the second operand.
   * \param equal_compare A binary comparison functor that returns true if both operands are equal 
   *        to each other.
   * \param combine A binary combination functor that returns the sum of its operands, 
   *        semantically-speaking.
   * \param compose A binary composition functor that amplifies a heuristic distance metric by a 
   *        scalar value (i.e. epsilon x h(u) ).
   */
  template <typename Graph,
	    typename Topology,
            typename AStarHeuristicMap,
            typename FADPRMVisitor,
	    typename PredecessorMap,
            typename DistanceMap,
	    typename RHSMap,
	    typename KeyMap,
            typename WeightMap,
            typename PositionMap,
            typename DensityMap,
	    typename NcSelector,
	    typename ColorMap,
	    typename AdjustFunction,
	    typename RunningPredicate,
	    typename GetStartNodeFunction,
	    typename ChangeEventFunction,
	    typename PublishPathFunction,
	    typename ScalarType,
            typename CompareFunction,
	    typename EqualCompareFunction,
            typename CombineFunction,
	    typename ComposeFunction>
  inline void
  generate_fadprm_no_init
    (Graph &g, const Topology& free_space,
     AStarHeuristicMap hval, FADPRMVisitor vis,
     PredecessorMap predecessor, DistanceMap distance,
     RHSMap rhs, KeyMap key, WeightMap weight, DensityMap density, PositionMap position, 
     NcSelector select_neighborhood, ColorMap color, AdjustFunction adj_epsilon, 
     RunningPredicate keep_going, GetStartNodeFunction get_start,
     ChangeEventFunction edge_change, PublishPathFunction publish_path,
     ScalarType epsilon, typename boost::property_traits<DistanceMap>::value_type inf,
     typename boost::property_traits<DistanceMap>::value_type zero = typename boost::property_traits<DistanceMap>::value_type(0),
     CompareFunction compare = CompareFunction(), EqualCompareFunction equal_compare = EqualCompareFunction(),
     CombineFunction combine = CombineFunction(), ComposeFunction compose = ComposeFunction())
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef typename boost::property_traits<DistanceMap>::value_type DistanceValue;
    typedef typename boost::property_traits<WeightMap>::value_type WeightValue;
    typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
    typedef typename adstar_key_traits<KeyValue>::compare_type KeyCompareType;
    typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
    IndexInHeapMap index_in_heap;
    {
      typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
      for (boost::tie(ui, ui_end) = boost::vertices(g); ui != ui_end; ++ui) {
        put(index_in_heap,*ui, (std::size_t)(-1)); //this ugly C-style cast is required to match the boost::d_ary_heap_indirect implementation.
      };
    };

    typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, KeyMap, KeyCompareType> MutableQueue;
    MutableQueue Q(key, index_in_heap, KeyCompareType()); //priority queue holding the OPEN set.
    std::vector<Vertex> I; //list holding the INCONS set (inconsistent nodes).


    detail::fadprm_bfs_visitor<Topology, AStarHeuristicMap, FADPRMVisitor,
        MutableQueue, std::vector<Vertex>, IndexInHeapMap, PredecessorMap, KeyMap, DistanceMap, RHSMap,
        WeightMap, DensityMap, PositionMap, NcSelector, ColorMap, ScalarType, CompareFunction,
	EqualCompareFunction, CombineFunction, ComposeFunction>
      bfs_vis(free_space, hval, vis, Q, I, index_in_heap, predecessor, key, distance, rhs, weight, density, position, select_neighborhood,
              color, epsilon, compare, equal_compare, combine, compose, inf, zero);

    std::vector<Edge> affected_edges;

    detail::adstar_search_loop(g, hval, bfs_vis, predecessor, distance, rhs, key, weight, 
			       color, Q, I, affected_edges, adj_epsilon, keep_going,
			       get_start, edge_change, publish_path, epsilon, 
			       inf, zero, compare, equal_compare, combine, compose);
    
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the FADPRM algorithm (as of Belghith et al., 2006), with initialization of the 
   * existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values 
   *         for each vertex in the graph.
   * \tparam FADPRMVisitor The type of the FADPRM visitor to be used, should model the FADPRMVisitorConcept.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of 
   *         each vertex to the goal (internal use to AD*).
   * \tparam KeyMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam PositionMap A property-map type that can store the position of each vertex. 
   * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
   *         See classes in the topological_search.hpp header-file.
   * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors 
   *         are used to mark vertices by their status in the AD* algorithm (white = not visited, 
   *         gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED, 
   *         not OPEN), red = inconsistent (in INCONS)).
   * \tparam AdjustFunction This functor type represents the callback function that will adjust the 
   *         epsilon value (sloppyness of the A* search).
   * \tparam RunningPredicate This function type represents the callback function that will evaluate 
   *         if it is still worth continuing the AD* search passes.
   * \tparam GetStartNodeFunction This function type represents the callback function that will give 
   *         the starting point for the search (usually the goal vertex).
   * \tparam ChangeEventFunction This function type represents the callback that can wait for a change 
   *         in edge weights, and provide the list of changed edges.
   * \tparam PublishPathFunction This function type represents the callback that is called when a 
   *         path has been obtained (so that it can be buffered and executed).
   * \tparam ScalarType The type of the scalar value epsilon (amplification factor, controls the 
   *         anytime nature of this algorithm).
   * \tparam CompareFunction A binary comparison functor type that returns true if the first operand 
   *         is strictly better (less-than) than the second operand.
   * \tparam EqualCompareFunction A binary comparison functor type that returns true if both operands 
   *         are equal to each other.
   * \tparam CombineFunction A binary combination functor type that returns the sum of its operands, 
   *         semantically-speaking.
   * \tparam ComposeFunction A binary composition functor type that amplifies a heuristic distance 
   *         metric by a scalar value (i.e. epsilon x h(u) ).
   * 
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated graph once the algorithm has finished.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *        that it is required to generate only random points in 
   *        the free-space and to only allow interpolation within the free-space.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param vis A FADPRM visitor implementing the FADPRMVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param predecessor The property-map which will store the resulting path by connecting 
   *        vertices together with their optimal predecessor (follow in reverse to discover the 
   *        complete path).
   * \param distance The property-map which stores the estimated distance of each vertex to the goal.
   * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the 
   *        goal (for internal use).
   * \param key The property-map which stores the AD* key-values associated to each vertex.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
   * \param density A property-map that provides the density values assiciated to each vertex.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param color The property-map which stores the color-value of the vertices, colors are used to mark
   *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), 
   *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
   * \param adj_epsilon The callback functor that will adjust the epsilon value (sloppyness of the sub-optimal
   *        A* search) depending (usually) on the level of change in the weights that was recorded after a callback
   *        using the edge_change functor. Used as: epsilon = adj_epsilon(epsilon, max_w_change); (takes the current 
   *        epsilon value and the maximum recorder weight change, and outputs the new epsilon value).
   * \param keep_going The callback functor that will tell whether the AD* search should keep on going (running 
   *        predicate). Should be callable with no parameter and output a bool value (true to keep going, false
   *        to stop).
   * \param get_start The functor that is called to obtain the starting point for the AD* search (usually the 
   *        goal vertex). Should be callable with no parameter and output a vertex-descriptor of the graph (i.e. 
   *        a vertex of the graph g).
   * \param edge_change The callback functor that is called to make a check whether the weights of the graph
   *        have changed do to a dynamic environment. Used as: max_w_change = edge_change(affected_edges); 
   *        (takes an STL container (std::vector) of affected edges and populates it all the edges that 
   *        should be updated because their weight has changed, also returns the maximum weight change 
   *        measured (or any other measure of how bad the change in the environment is, this value will 
   *        be forwarded unchanged to the adj_epsilon functor)). The container of affected edges is subject 
   *        to change, so it is recommended to define this functor with a templated call-operator that 
   *        can take a generic type of STL container, by non-const reference.
   * \param publish_path The callback functor that is called when the AD* search has ended, and presumably,
   *        found a path. Used as: publish_path(); (takes no parameters, returns nothing). The path is retrieved
   *        by walking the predecessors of the vertices (using the given predecessor property-map).
   * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime 
   *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
   * \param inf The quantity that represents infinity (either a very large value or the infinity value for 
   *        the underlying value-type).
   * \param zero The quantity that represents zero with the given value-type.
   * \param compare A binary comparison functor that returns true if the first operand is strictly 
   *        better (less-than) than the second operand.
   * \param equal_compare A binary comparison functor that returns true if both operands are equal 
   *        to each other.
   * \param combine A binary combination functor that returns the sum of its operands, 
   *        semantically-speaking.
   * \param compose A binary composition functor that amplifies a heuristic distance metric by a 
   *        scalar value (i.e. epsilon x h(u) ).
   */
  template <typename Graph, //this is the actual graph, should comply to BidirectionalMutableGraphConcept.
	    typename Topology,
            typename AStarHeuristicMap, //this the map of heuristic function value for each vertex.
            typename FADPRMVisitor, //this is a visitor class that can perform special operations at event points.
	    typename PredecessorMap, //this is the map that stores the preceeding edge for each vertex.
            typename DistanceMap, //this is the map of distance values associated with each vertex.
	    typename RHSMap,
	    typename KeyMap, //this is the map of key values associated to each vertex.
            typename WeightMap, //this is the map of edge weight (or cost) associated to each edge of the graph.
            typename DensityMap,
            typename PositionMap,
	    typename NcSelector,
	    typename ColorMap, //this is a color map for each vertex, i.e. white=not visited, gray=discovered, black=expanded.
	    typename AdjustFunction, //this is the function object type that adjusts scalar implification factor epsilon for the current max change in edge weight.
	    typename RunningPredicate, //this is a function object type that immediately returns true if the algorithm should keep going.
	    typename GetStartNodeFunction, //this is a function object type that returns the current start node of the graph (to start the bfs).
	    typename ChangeEventFunction, //this is a function object type that can wait for a change in edge weights, and provide the list of changed edges.
	    typename PublishPathFunction, //this is a function object type that is called when a path has been obtained (so that it can be buffered and executed).
	    typename ScalarType/* = typename property_traits<DistanceMap>::value_type*/, //the type of the scalar value epsilon (amplification factor, controls the anytime nature of this algorithm).
            typename CompareFunction /*= std::less<typename property_traits<DistanceMap>::value_type>*/, //a binary comparison function object that returns true if the first operand is strictly better (less-than) than the second operand.
	    typename EqualCompareFunction /*= std::equal_to<typename property_traits<DistanceMap>::value_type >*/, //a binary comparison function object that returns true if both operands are equal to each other.
            typename CombineFunction /*= std::plus<typename property_traits<DistanceMap>::value_type>*/, //a binary combination function object that returns the sum of its operands (sum in the broad sense).
	    typename ComposeFunction /*= std::multiplies<typename property_traits<DistanceMap>::value_type>*/ > //a binary composition function object that amplifies a heuristic distance metric by a scalar value (i.e. epsilon x h(u) ).
  inline void
  generate_fadprm
    (Graph &g, const Topology& free_space,
     AStarHeuristicMap hval, FADPRMVisitor vis,
     PredecessorMap predecessor, DistanceMap distance,
     RHSMap rhs, KeyMap key, WeightMap weight, DensityMap density, PositionMap position, NcSelector select_neighborhood,
     ColorMap color, AdjustFunction adj_epsilon, RunningPredicate keep_going, GetStartNodeFunction get_start,
     ChangeEventFunction edge_change, PublishPathFunction publish_path,
     ScalarType epsilon, typename boost::property_traits<DistanceMap>::value_type inf,
     typename boost::property_traits<DistanceMap>::value_type zero = typename boost::property_traits<DistanceMap>::value_type(0),
     CompareFunction compare = CompareFunction(), EqualCompareFunction equal_compare = EqualCompareFunction(),
     CombineFunction combine = CombineFunction(), ComposeFunction compose = ComposeFunction())
  {

    typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    
    if(boost::num_vertices(g) == 0) {
      Vertex u = boost::add_vertex(g);
      PositionValue p = free_space.random_point();
      put(position, u, p);
      vis.vertex_added(u, g);
    };
    
    for (boost::tie(ui, ui_end) = boost::vertices(g); ui != ui_end; ++ui) {
      put(color, *ui, Color::white());
      put(distance, *ui, inf);
      put(rhs, *ui, inf);
      put(key, *ui, KeyValue(inf,inf,compare,equal_compare));
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    generate_fadprm_no_init
      (g, free_space, hval, vis, predecessor, distance, rhs, key, weight, density, position, select_neighborhood,
       color, adj_epsilon, keep_going, get_start, edge_change, publish_path,
       epsilon, inf, zero, compare, equal_compare, combine, compose);

  };
  
  
  

};

};

#endif
















