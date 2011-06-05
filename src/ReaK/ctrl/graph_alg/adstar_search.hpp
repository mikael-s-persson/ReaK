/**
 * \file adstar_search.hpp
 *
 * \author S. Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 9th, 2011
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


#ifndef ADSTAR_SEARCH_HPP
#define ADSTAR_SEARCH_HPP

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

  template <typename Visitor, typename Graph>
  struct ADStarVisitorConcept {
    void constraints()
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




  template <typename DistanceValueType,
            typename CompareFunction = std::less<DistanceValueType>,
	    typename EqualCompareFunction = std::equal_to<DistanceValueType> >
  struct adstar_key_value {
    adstar_key_value() : m_k1(0), m_k2(0), m_compare(), m_equal() { };
    adstar_key_value(DistanceValueType k1, DistanceValueType k2,
                     CompareFunction compare, EqualCompareFunction equalCompare) :
                     m_k1(k1), m_k2(k2), m_compare(compare), m_equal(equalCompare) { };

    bool operator <(const adstar_key_value<DistanceValueType,CompareFunction,EqualCompareFunction>& aKey) const {
      return m_compare(m_k1, aKey.m_k1) || (m_equal(m_k1, aKey.m_k1) && m_compare(m_k2, aKey.m_k2));
    };

    DistanceValueType m_k1, m_k2;
    CompareFunction m_compare;
    EqualCompareFunction m_equal;
  };

  template <typename ADStarKeyType>
  struct adstar_key_traits {
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



  template <typename VertexListGraph, //this is the actual graph, should comply to BidirectionalGraphConcept.
            typename AStarHeuristicMap, //this the map of heuristic function value for each vertex.
            typename ADStarVisitor, //this is a visitor class that can perform special operations at event points.
	    typename PredecessorMap, //this is the map that stores the preceeding edge for each vertex.
            typename DistanceMap, //this is the map of distance values associated with each vertex.
	    typename RHSMap,
	    typename KeyMap, //this is the map of key values associated to each vertex.
            typename WeightMap, //this is the map of edge weight (or cost) associated to each edge of the graph.
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
      for (boost::tie(ui, ui_end) = boost::vertices(g); ui != ui_end; ++ui) {
        put(index_in_heap,*ui, (std::size_t)(-1)); //this ugly C-style cast is required to match the boost::d_ary_heap_indirect implementation.
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


  // Non-named parameter interface
  template <typename VertexListGraph, //this is the actual graph, should comply to BidirectionalGraphConcept.
            typename AStarHeuristicMap, //this the map of heuristic function value for each vertex.
            typename ADStarVisitor, //this is a visitor class that can perform special operations at event points.
	    typename PredecessorMap, //this is the map that stores the preceeding edge for each vertex.
            typename DistanceMap, //this is the map of distance values associated with each vertex.
	    typename RHSMap,
	    typename KeyMap, //this is the map of key values associated to each vertex.
            typename WeightMap, //this is the map of edge weight (or cost) associated to each edge of the graph.
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
    for (boost::tie(ui, ui_end) = boost::vertices(g); ui != ui_end; ++ui) {
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










