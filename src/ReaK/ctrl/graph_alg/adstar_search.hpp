/**
 * \file adstar_search.hpp
 *
 * This library implements the Anytime Dynamic A* (AD*) algorithm according to the algorithmic description given
 * in the original article:
 *
 * M. Likhachev, D. Ferguson, G. Gordon, A. Stentz and S. Thrun (2005). "Anytime Dynamic A*: An Anytime, Replanning
 * Algorithm". Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS), June, 2005.
 *
 * The algorithm assumes a fixed starting point for the search (usually the fixed goal of a path-planner) and an A*
 * heuristic function that can be computed for any vertex that is being explored. The underlying graph search is an
 * A* algorithm that was modified for sub-optimal, anytime computation. It will start at a given initial
 * epsilon value and rapidly find a epsilon-optimal path. Afterwards, the espilon can be adjusted between each run
 * of the graph search as a function of the amount of change in the environment (usually decreasing while no important
 * changes are reported, and a reset to the initial value if important changes occur). The AD* will keep in memory the
 * previously computed solution to improve upon it at future graph searches. AD* uses two callback functions: one to
 * publish a path once it has been resolved by the graph search; and another to query for the list of weights of the
 * graph that have undergone some (significant) change and a scalar value that is representative of the amount of
 * change. For all other aspects related to observing the changes in the graph during the graph search, the A* visitor
 * concept is used, as documented in the Boost.Graph library's homepage.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2011
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


#ifndef REAK_ADSTAR_SEARCH_HPP
#define REAK_ADSTAR_SEARCH_HPP

#include <functional>
#include <vector>
#include <boost/limits.hpp>
#include <boost/graph/exception.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>
#include <boost/property_map/property_map.hpp>


namespace boost {

  enum vertex_heuristic_t { vertex_heuristic };
  enum vertex_rhs_t { vertex_rhs };
  enum vertex_key_t { vertex_key };

  BOOST_INSTALL_PROPERTY(vertex, heuristic);
  BOOST_INSTALL_PROPERTY(vertex, rhs);
  BOOST_INSTALL_PROPERTY(vertex, key);

};

namespace ReaK {
  
namespace graph {

  /**
   * This concept class defines the valid expressions required of a class to be used as a visitor 
   * class for the AD* algorithm. A visitor class is essentially a class that regroups a number of 
   * callback functions that can be used to inject customization into the AD* algorithm (see fadprm.hpp for 
   * example). In other words, the visitor pattern in generic programming is an implementation of IoC 
   * (Inversion of Control), since the AD* algorithm is in control of execution, but custom behavior can
   * be injected in several places, even blocking the algorithm if needed.
   * 
   * Required concepts:
   * 
   * The visitor class should model the boost::CopyConstructibleConcept.
   * 
   * Valid expressions:
   * 
   * vis.initialize_vertex(u,g);  A function that gets called whenever a vertex (u) is first initialized before the search.
   * 
   * vis.finish_vertex(u,g);  A function that gets called whenever a vertex (u) is put into the CLOSED set (after being explored by the current search pass).
   * 
   * vis.recycle_vertex(u,g);  A function that gets called whenever a vertex (u) is taken out of the CLOSED set to be recycled in the search.
   * 
   * vis.discover_vertex(u,g);  A function that gets called whenever a vertex (u) is added to the OPEN set (or updated in the OPEN set).
   * 
   * vis.examine_vertex(u,g);  A function that gets called whenever a vertex (u) is taken out of the OPEN set to be examined, this is called before it gets expanded.
   * 
   * vis.examine_edge(e,g);  A function that gets called whenever an edge (e) is being looked at, as it comes out of the vertex that is currently being examined in the search.
   * 
   * vis.edge_relaxed(e,g);  A function that gets called whenever an edge (e) has been newly declared as a better alternative than the current surrounding edges (i.e. the edge is added to the optimal path).
   * 
   * vis.forget_vertex(u,g);  A function that gets called whenever a vertex (u) is deemed uninteresting for the search and is taken out of the OPEN set and not examined.
   * 
   * vis.inconsistent_vertex(u,g);  A function that gets called whenever a CLOSEd vertex becomes inconsistent (added to the INCONS set) due to changes in the environment or sub-optimal search.
   * 
   * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   */
  template <typename Visitor, typename Graph>
  struct ADStarVisitorConcept {
    BOOST_CONCEPT_USAGE(ADStarVisitorConcept)
    {
      boost::function_requires< boost::CopyConstructibleConcept<Visitor> >();
      vis.initialize_vertex(u, g);   //whenever the vertex is first initialized.
      vis.finish_vertex(u, g);       //whenever a vertex is added to the CLOSED set.
      vis.recycle_vertex(u, g);      //whenever a vertex is taken out of the CLOSED set.
      vis.discover_vertex(u, g);     //whenever a vertex is added to the OPEN set (or updated in OPEN).
      vis.examine_vertex(u, g);      //whenever a vertex is taken out of OPEN, before it gets "expanded".
      vis.examine_edge(e, g);        //whenever an edge is being looked at (an out_edge of the vertex under examination).
      vis.edge_relaxed(e, g);        //whenever it is newly decided that an edge is relaxed (has improved the distance for its target)
      vis.forget_vertex(u, g);       //whenever a vertex is deemed uninteresting and is taken out of OPEN, but not yet expanded.
      vis.inconsistent_vertex(u, g); //whenever a closed vertex becomes INCONS.
    }
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
  };

  /**
   * This class is the default implementation of an AD* visitor (see ADStarVisitorConcept). 
   * Basically, this implementation models the concept required by AD*, but does nothing at all 
   * (all functions are empty).
   */
  class default_adstar_visitor {
    public:
      default_adstar_visitor() {}

      template <typename Vertex, typename Graph>
      void initialize_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void discover_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void inconsistent_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void examine_vertex(Vertex u, Graph& g) { };
      template <typename Edge, typename Graph>
      void examine_edge(Edge e, Graph& g) { };
      template <typename Edge, typename Graph>
      void edge_relaxed(Edge e, Graph& g) { };
      template <typename Vertex, typename Graph>
      void forget_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void finish_vertex(Vertex u, Graph& g) { };
      template <typename Vertex, typename Graph>
      void recycle_vertex(Vertex u, Graph& g) { };

  };



  /**
   * This class template is used by the AD* algorithm to constitute the key-values which 
   * drive the ordering in the priority-queue that the AD* algorithm uses to choose the next 
   * vertex to examine. This class simply implements the key-value described in the original
   * paper on AD* (ICAPS 2005).
   * \tparam DistanceValueType The scalar type that describes the distance values.
   * \tparam CompareFunction The strict weak-ordering function that can sort the distance values.
   * \tparam EqualCompareFunction The equal-comparison function that can compare two distance values to be equal.
   */
  template <typename DistanceValueType,
            typename CompareFunction = std::less<DistanceValueType>,
	    typename EqualCompareFunction = std::equal_to<DistanceValueType> >
  struct adstar_key_value {
    /**
     * Default constructor.
     */
    adstar_key_value() : m_k1(0), m_k2(0), m_compare(), m_equal() { };
    /**
     * Parametrized constructor.
     * \param k1 The first value of the key-value.
     * \param k2 The second value of the key-value.
     * \param compare The functor of type CompareFunction.
     * \param equalCompare The functor of type EqualCompareFunction.
     */
    adstar_key_value(DistanceValueType k1, DistanceValueType k2,
                     CompareFunction compare, EqualCompareFunction equalCompare) :
                     m_k1(k1), m_k2(k2), m_compare(compare), m_equal(equalCompare) { };

    /**
     * The less-than operator for strict weak-ordering.
     */
    bool operator <(const adstar_key_value<DistanceValueType,CompareFunction,EqualCompareFunction>& aKey) const {
      return m_compare(m_k1, aKey.m_k1) || (m_equal(m_k1, aKey.m_k1) && m_compare(m_k2, aKey.m_k2));
    };

    DistanceValueType m_k1, m_k2;
    CompareFunction m_compare;
    EqualCompareFunction m_equal;
  };

  /**
   * This traits class defines that traits that an AD* key-value should have.
   * \tparam ADStarKeyType The key-value type of which the traits are sought.
   */
  template <typename ADStarKeyType>
  struct adstar_key_traits {
    /** The type of comparison to use for strict weak-ordering of the key-values. */
    typedef std::less< ADStarKeyType > compare_type;
  };





  namespace detail {

    template <typename AStarHeuristicMap, typename UniformCostVisitor,
              typename UpdatableQueue, typename List, typename PredecessorMap,
	      typename KeyMap, typename DistanceMap, typename RHSMap, typename WeightMap,
	      typename ColorMap, typename ScalarType, typename CompareFunction,
	      typename EqualCompareFunction, typename CombineFunction,
	      typename ComposeFunction>
    struct adstar_bfs_visitor
    {

      typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
      typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
      typedef boost::color_traits<ColorValue> Color;
      typedef typename boost::property_traits<DistanceMap>::value_type distance_type;
      typedef typename boost::property_traits<WeightMap>::value_type weight_type;

      adstar_bfs_visitor(AStarHeuristicMap h, UniformCostVisitor vis,
                         UpdatableQueue& Q, List& I, PredecessorMap p,
                         KeyMap k, DistanceMap d, RHSMap rhs, WeightMap w,
                         ColorMap col, const ScalarType& epsilon,
                         CompareFunction compare, EqualCompareFunction equal_compare,
			 CombineFunction combine, ComposeFunction compose,
			 distance_type inf, distance_type zero)
        : m_h(h), m_vis(vis), m_Q(Q), m_I(I), m_predecessor(p), m_key(k),
          m_distance(d), m_rhs(rhs), m_weight(w), m_color(col), m_epsilon(epsilon),
          m_compare(compare), m_equal_compare(equal_compare), m_combine(combine),
          m_compose(compose), m_inf(inf), m_zero(zero) {};

      adstar_bfs_visitor(const adstar_bfs_visitor<AStarHeuristicMap, UniformCostVisitor,
                                                  UpdatableQueue, List, PredecessorMap,
						  KeyMap, DistanceMap, RHSMap, WeightMap,
						  ColorMap, ScalarType, CompareFunction,
						  EqualCompareFunction, CombineFunction,
						  ComposeFunction>& aVis)
        : m_h(aVis.m_h), m_vis(aVis.m_vis), m_Q(aVis.m_Q), m_I(aVis.m_I), m_predecessor(aVis.m_predecessor),
          m_key(aVis.m_key), m_distance(aVis.m_distance), m_rhs(aVis.m_rhs), m_weight(aVis.m_weight),
          m_color(aVis.m_color), m_epsilon(aVis.m_epsilon), m_compare(aVis.m_compare),
          m_equal_compare(aVis.m_equal_compare), m_combine(aVis.m_combine),
          m_compose(aVis.m_compose), m_inf(aVis.m_inf), m_zero(aVis.m_zero) {};

      template <typename Vertex, typename Graph>
      void initialize_vertex(Vertex u, Graph& g) {
        m_vis.initialize_vertex(u, g);
      };
      template <typename Vertex, typename Graph>
      void discover_vertex(Vertex u, Graph& g) {
        m_vis.discover_vertex(u, g);
      };
      template <typename Vertex, typename Graph>
      void inconsistent_vertex(Vertex u, Graph& g) {
	m_vis.inconsistent_vertex(u, g);
      };
      template <typename Vertex, typename Graph>
      void examine_vertex(Vertex u, Graph& g) {
        m_vis.examine_vertex(u, g);
      };
      template <typename Edge, typename Graph>
      void examine_edge(Edge e, Graph& g) {
        if (m_compare(get(m_weight, e), m_zero))
          throw boost::negative_edge();
        m_vis.examine_edge(e, g);
      };
      template <typename Vertex, typename Graph>
      void forget_vertex(Vertex u, Graph& g) {
        m_vis.forget_vertex(u, g);
      };
      template <typename Vertex, typename Graph>
      void finish_vertex(Vertex u, Graph& g) {
	m_vis.finish_vertex(u, g);
      };
      template <typename Vertex, typename Graph>
      void recycle_vertex(Vertex u, Graph& g) {
	m_vis.recycle_vertex(u, g);
      };

      template <typename Vertex, typename Graph>
      void update_key(Vertex u, Graph&) {
	distance_type g_u = get(m_distance, u);
	distance_type rhs_u = get(m_rhs, u);
	if( m_compare( rhs_u, g_u ) )
	  put(m_key, u, KeyValue( m_combine(rhs_u, m_compose(m_epsilon, get(m_h, u))), rhs_u, m_compare, m_equal_compare));
	else
	  put(m_key, u, KeyValue( m_combine(g_u, get(m_h, u)), g_u, m_compare, m_equal_compare));
      };

      template <typename Vertex, typename BidirectionalGraph>
      void update_vertex(Vertex u, BidirectionalGraph& g) {
	BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<BidirectionalGraph>));
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

	if(!m_equal_compare(rhs_u,m_zero)) { //if u is not the start node.
	  rhs_u = m_inf; //This was in the original code!
	  //Vertex pred_u = u;
	  //Vertex pred_u = get(m_predecessor, u); //this alone, works. Because it doesn't make anything inconsistent!
	  typename GTraits::edge_descriptor pred_e;
	  for(tie(ei,ei_end) = in_edges(u,g); ei != ei_end; ++ei) {
	    distance_type rhs_tmp = m_combine(get(m_weight, *ei), get(m_distance, source(*ei,g)));
	    if(m_compare(rhs_tmp, rhs_u)) {
	      rhs_u = rhs_tmp;
	      //pred_u = source(*ei,g);
	      put(m_rhs, u, rhs_u);  //this was the original code!
	      put(m_predecessor, u, source(*ei,g));
	      pred_e = *ei;
	    };
	  };
	  rhs_u = get(m_rhs, u); //put(m_rhs, u, rhs_u);
	  //if(pred_u != get(m_predecessor, u))                          m_vis.edge_relaxed(pred_e, g);
	  //put(m_predecessor, u, pred_u);
	};

	if(!m_equal_compare(rhs_u,g_u)) {
	  if((col_u != Color::black()) && (col_u != Color::red())) { //if not in CLOSED set (i.e. either just closed (black) or inconsistent (red)).
	    update_key(u,g);
	    m_Q.push_or_update(u);
	    put(m_color, u, Color::gray());                            m_vis.discover_vertex(u, g);
	  } else if(col_u == Color::black()) {
	    m_I.push_back(u);
	    put(m_color, u, Color::red());                             m_vis.inconsistent_vertex(u,g);
	  };
	} else if(m_Q.contains(u)) { //if u is in the OPEN set, then remove it.
	  //KeyValue k = get(m_key, u);
	  put(m_key, u, KeyValue(-m_inf,-m_inf,m_compare,m_equal_compare));
	  m_Q.update(u);
	  m_Q.pop(); //remove from OPEN set
	  //put(m_key, u, k);
	  update_key(u, g); //this was the original code!
	  put(m_color, u, Color::green());                             m_vis.forget_vertex(u, g);
	};
      };


      AStarHeuristicMap m_h;
      UniformCostVisitor m_vis;
      UpdatableQueue& m_Q;
      List& m_I;
      PredecessorMap m_predecessor;
      KeyMap m_key;
      DistanceMap m_distance;
      RHSMap m_rhs;
      WeightMap m_weight;
      ColorMap m_color;
      const ScalarType& m_epsilon;
      CompareFunction m_compare;
      EqualCompareFunction m_equal_compare;
      CombineFunction m_combine;
      ComposeFunction m_compose;
      distance_type m_inf;
      distance_type m_zero;
    };
    
    
    template <typename VertexListGraph, //this is the actual graph, should comply to BidirectionalGraphConcept.
              typename AStarHeuristicMap, //this the map of heuristic function value for each vertex.
              typename ADStarBFSVisitor, //this is a visitor class that can perform special operations at event points.
	      typename PredecessorMap, //this is the map that stores the preceeding edge for each vertex.
              typename DistanceMap, //this is the map of distance values associated with each vertex.
	      typename RHSMap,
	      typename KeyMap, //this is the map of key values associated to each vertex.
              typename WeightMap, //this is the map of edge weight (or cost) associated to each edge of the graph.
	      typename ColorMap, //this is a color map for each vertex, i.e. white=not visited, gray=discovered, black=expanded.
	      typename MutableQueue,
	      typename InconsList,
	      typename AffectedEdgeList,
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
    adstar_search_loop
      (VertexListGraph &g,
       AStarHeuristicMap hval, ADStarBFSVisitor& bfs_vis,
       PredecessorMap predecessor, DistanceMap distance,
       RHSMap rhs, KeyMap key, WeightMap weight,
       ColorMap color, MutableQueue& Q, InconsList& I, AffectedEdgeList& affected_edges,
       AdjustFunction adj_epsilon, RunningPredicate keep_going,
       GetStartNodeFunction get_start,
       ChangeEventFunction edge_change, PublishPathFunction publish_path,
       ScalarType epsilon, typename boost::property_traits<DistanceMap>::value_type inf,
       typename boost::property_traits<DistanceMap>::value_type zero = typename boost::property_traits<DistanceMap>::value_type(0),
       CompareFunction compare = CompareFunction(), EqualCompareFunction equal_compare = EqualCompareFunction(),
       CombineFunction combine = CombineFunction(), ComposeFunction compose = ComposeFunction())
    {
      using namespace boost;
      typedef typename graph_traits<VertexListGraph>::vertex_descriptor Vertex;
      typedef typename graph_traits<VertexListGraph>::edge_descriptor Edge;
      typedef typename property_traits<ColorMap>::value_type ColorValue;
      typedef color_traits<ColorValue> Color;
      typedef typename property_traits<DistanceMap>::value_type DistanceValue;
      typedef typename property_traits<WeightMap>::value_type WeightValue;
      typedef typename property_traits<KeyMap>::value_type KeyValue;
      typedef typename adstar_key_traits<KeyValue>::compare_type KeyCompareType;
      
      Vertex s = get_start();
      put(distance, s, inf);
      put(rhs, s, zero);
      put(predecessor, s, s);
      bfs_vis.update_key(s,g);
    
      put(color, s, Color::gray());                          bfs_vis.discover_vertex(s, g);
      Q.push(s);
      
      while(keep_going()) {

        typename graph_traits<VertexListGraph>::out_edge_iterator eig, eig_end;
        typename graph_traits<VertexListGraph>::in_edge_iterator eii, eii_end;

        //put(color, s, Color::gray());                          bfs_vis.discover_vertex(s, g);
        //Q.push(s);              //this was in the original code!
        while (! Q.empty()) { 
          Vertex u = Q.top(); Q.pop();  bfs_vis.examine_vertex(u, g);
	  //put(color, u, Color::green());   /*this was an addition from working version*/   
	  DistanceValue g_u = get(distance, u);
	  DistanceValue rhs_u = get(rhs,u);
	  //if( KeyCompareType()(get(key,s), get(key,u)) && equal_compare(get(rhs,s),get(distance,s)) ) {
	  if( equal_compare(get(hval, u), zero) && equal_compare(g_u,rhs_u) ) { //if we have a consistent node at the goal
	    break;
	  };
	  if(compare(rhs_u, g_u)) { //if g_u is greater than rhs_u, then make u consistent and close it.
	    put(distance, u, rhs_u);
	    g_u = rhs_u; 
	    put(color, u, Color::black());                     bfs_vis.finish_vertex(u, g);
	  } else { 
	    put(distance, u, inf);
	    bfs_vis.update_vertex(u,g); 
 	  };
	  for (tie(eig, eig_end) = out_edges(u, g); eig != eig_end; ++eig) {
	    bfs_vis.examine_edge(*eig, g); 
	    bfs_vis.update_vertex(target(*eig, g),g); 
	  };
        }; // end while
        
        publish_path();
	
        affected_edges.clear();
        WeightValue max_w_change = edge_change(affected_edges);

        //s = get_start();          //this was in the original code, but the way it was, it had no effect (I mean none at all).
        //put(distance, s, inf);
        //put(rhs, s, zero);
        //put(predecessor, s, s);
        //bfs_vis.update_key(s,g);

        //update all nodes that were affected.
        for(typename AffectedEdgeList::iterator ei = affected_edges.begin(); ei != affected_edges.end(); ++ei) {
	  if( get(color, source(*ei,g)) == Color::black() )
	    put(color, source(*ei,g), Color::green());
	  bfs_vis.update_vertex(source(*ei, g),g);
	  if( get(color, target(*ei,g)) == Color::black() )
	    put(color, target(*ei,g), Color::green());
	  bfs_vis.update_vertex(target(*ei, g),g);
        };

        epsilon = adj_epsilon(epsilon, max_w_change);

        //merge the OPEN and INCONS sets
        for(typename InconsList::iterator ui = I.begin(); ui != I.end(); ++ui) {
	  bfs_vis.update_key(*ui,g);
	  Q.push_or_update(*ui);
	  put(color, *ui, Color::gray());
        };
        I.clear();

        //update keys for all OPEN nodes, and change all black nodes to green (empty the CLOSED set).
        {
          typename graph_traits<VertexListGraph>::vertex_iterator ui, ui_end;
          for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
	    ColorValue u_color = get(color, *ui);
	    if( Q.contains(*ui) ) {
	      bfs_vis.update_key(*ui,g);
	      Q.update(*ui);
	      put(color, *ui, Color::gray());                  bfs_vis.discover_vertex(*ui, g);
	    } else if( u_color == Color::black() ) {
	      put(color, *ui, Color::green());                 bfs_vis.recycle_vertex(*ui, g);
            };
          };
        };
      };
    };
    
    
    
    
    
    

  }; // namespace detail



  /**
   * This function template performs an AD* search over a graph, without initialization. The AD* search 
   * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
   * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then 
   * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it 
   * uses callback functions to allow the user to check for changes in the environment, triggering an update
   * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
   * \tparam VertexListGraph The type of the graph on which the search is performed, should model 
   *         the BGL's BidirectionalGraphConcept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values 
   *         for each vertex in the graph.
   * \tparam ADStarVisitor The type of the AD* visitor to be used, should model the ADStarVisitorConcept.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each 
   *         vertex to the goal.
   * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of 
   *         each vertex to the goal (internal use to AD*).
   * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
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
   * \param g The graph on which to apply the AD* algorithm.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param vis The AD* visitor object, should model ADStarVisitorConcept.
   * \param predecessor The property-map which will store the resulting path by connecting 
   *        vertices together with their optimal predecessor (follow in reverse to discover the 
   *        complete path).
   * \param distance The property-map which stores the estimated distance of each vertex to the goal.
   * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the 
   *        goal (for internal use).
   * \param key The property-map which stores the AD* key-values associated to each vertex.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
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
  template <typename VertexListGraph,
            typename AStarHeuristicMap,
            typename ADStarVisitor,
	    typename PredecessorMap,
            typename DistanceMap,
	    typename RHSMap,
	    typename KeyMap,
            typename WeightMap,
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
  adstar_search_no_init
    (VertexListGraph &g,
     AStarHeuristicMap hval, ADStarVisitor vis,
     PredecessorMap predecessor, DistanceMap distance,
     RHSMap rhs, KeyMap key, WeightMap weight,
     ColorMap color, AdjustFunction adj_epsilon, RunningPredicate keep_going,
     GetStartNodeFunction get_start,
     ChangeEventFunction edge_change, PublishPathFunction publish_path,
     ScalarType epsilon, typename boost::property_traits<DistanceMap>::value_type inf,
     typename boost::property_traits<DistanceMap>::value_type zero = typename boost::property_traits<DistanceMap>::value_type(0),
     CompareFunction compare = CompareFunction(), EqualCompareFunction equal_compare = EqualCompareFunction(),
     CombineFunction combine = CombineFunction(), ComposeFunction compose = ComposeFunction())
  {
    typedef typename boost::graph_traits<VertexListGraph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<VertexListGraph>::edge_descriptor Edge;
    typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef typename boost::property_traits<DistanceMap>::value_type DistanceValue;
    typedef typename boost::property_traits<WeightMap>::value_type WeightValue;
    typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
    typedef typename adstar_key_traits<KeyValue>::compare_type KeyCompareType;
    typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
    IndexInHeapMap index_in_heap;
    {
      typename boost::graph_traits<VertexListGraph>::vertex_iterator ui, ui_end;
      for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
        put(index_in_heap,*ui, static_cast<std::size_t>(-1)); //this ugly C-style cast is required to match the boost::d_ary_heap_indirect implementation.
      };
    };

    typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, KeyMap, KeyCompareType> MutableQueue;
    MutableQueue Q(key, index_in_heap, KeyCompareType()); //priority queue holding the OPEN set.
    std::vector<Vertex> I; //list holding the INCONS set (inconsistent nodes).


    detail::adstar_bfs_visitor<AStarHeuristicMap, ADStarVisitor,
        MutableQueue, std::vector<Vertex>, PredecessorMap, KeyMap, DistanceMap, RHSMap,
        WeightMap, ColorMap, ScalarType, CompareFunction,
	EqualCompareFunction, CombineFunction, ComposeFunction>
      bfs_vis(hval, vis, Q, I, predecessor, key, distance, rhs, weight,
              color, epsilon, compare, equal_compare, combine, compose, inf, zero);

    std::vector<Edge> affected_edges;

    detail::adstar_search_loop(g, hval, bfs_vis, predecessor, distance, rhs, key, weight, 
			       color, Q, I, affected_edges, adj_epsilon, keep_going,
			       get_start, edge_change, publish_path, epsilon, 
			       inf, zero, compare, equal_compare, combine, compose);
  };



  /**
   * This function template performs an AD* search over a graph. The AD* search 
   * uses a sequence of sub-optimal A* searches of increasing level of optimality (i.e. it performs sloppy
   * or relaxed A* searches which yields results quicker by searching a subset of the graph only, and then 
   * it tightens the sloppiness of the search, improving on the results of previous searches). Then, it 
   * uses callback functions to allow the user to check for changes in the environment, triggering an update
   * of the affected edges and vertices, and also, possibly relaxing the A* search if changes are too significant.
   * \tparam VertexListGraph The type of the graph on which the search is performed, should model the BGL's BidirectionalGraphConcept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values for each vertex in the graph.
   * \tparam ADStarVisitor The type of the AD* visitor to be used, should model the ADStarVisitorConcept.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting vertex together with its optimal predecessor.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex to the goal.
   * \tparam RHSMap This property-map type is used to store the inconsistent estimated distance of each vertex to the goal (internal use to AD*).
   * \tparam KeyMap This property-map type is used to store the AD* key-values associated to each vertex.
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
   * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors are used to mark vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
   * \tparam AdjustFunction This functor type represents the callback function that will adjust the epsilon value (sloppyness of the A* search).
   * \tparam RunningPredicate This function type represents the callback function that will evaluate if it is still worth continuing the AD* search passes.
   * \tparam GetStartNodeFunction This function type represents the callback function that will give the starting point for the search (usually the goal vertex).
   * \tparam ChangeEventFunction This function type represents the callback that can wait for a change in edge weights, and provide the list of changed edges.
   * \tparam PublishPathFunction This function type represents the callback that is called when a path has been obtained (so that it can be buffered and executed).
   * \tparam ScalarType The type of the scalar value epsilon (amplification factor, controls the anytime nature of this algorithm).
   * \tparam CompareFunction A binary comparison functor type that returns true if the first operand is strictly better (less-than) than the second operand.
   * \tparam EqualCompareFunction A binary comparison functor type that returns true if both operands are equal to each other.
   * \tparam CombineFunction A binary combination functor type that returns the sum of its operands, semantically-speaking.
   * \tparam ComposeFunction A binary composition functor type that amplifies a heuristic distance metric by a scalar value (i.e. epsilon x h(u) ).
   * 
   * \param g The graph on which to apply the AD* algorithm.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param vis The AD* visitor object, should model ADStarVisitorConcept.
   * \param predecessor The property-map which will store the resulting path by connecting 
   *        vertices together with their optimal predecessor (follow in reverse to discover the 
   *        complete path).
   * \param distance The property-map which stores the estimated distance of each vertex to the goal.
   * \param rhs The property-map which stores the inconsistent estimated distance of each vertex to the 
   *        goal (for internal use).
   * \param key The property-map which stores the AD* key-values associated to each vertex.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
   *        along the edge).
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
   * \param compare A binary comparison functor that returns true if the first operand is strictly better (less-than) than the second operand.
   * \param equal_compare A binary comparison functor that returns true if both operands are equal to each other.
   * \param combine A binary combination functor that returns the sum of its operands, semantically-speaking.
   * \param compose A binary composition functor that amplifies a heuristic distance metric by a scalar value (i.e. epsilon x h(u) ).
   */
  template <typename VertexListGraph,
            typename AStarHeuristicMap,
            typename ADStarVisitor,
	    typename PredecessorMap,
            typename DistanceMap,
	    typename RHSMap,
	    typename KeyMap,
            typename WeightMap,
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
  adstar_search
    (VertexListGraph &g,
     AStarHeuristicMap hval, ADStarVisitor vis,
     PredecessorMap predecessor, DistanceMap distance,
     RHSMap rhs, KeyMap key, WeightMap weight,
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
    typename boost::graph_traits<VertexListGraph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(color, *ui, Color::white());
      put(distance, *ui, inf);
      put(rhs, *ui, inf);
      put(key, *ui, KeyValue(inf,inf,compare,equal_compare));
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    adstar_search_no_init
      (g, hval, vis, predecessor, distance, rhs, key, weight,
       color, adj_epsilon, keep_going, get_start, edge_change, publish_path,
       epsilon, inf, zero, compare, equal_compare, combine, compose);

  };



};

};


#endif










