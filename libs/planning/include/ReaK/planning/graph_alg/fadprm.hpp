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

#ifndef REAK_FADPRM_HPP
#define REAK_FADPRM_HPP

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
template < typename Visitor, typename Graph, typename Topology >
struct FADPRMVisitorConcept {
  BOOST_CONCEPT_USAGE( FADPRMVisitorConcept ) {
    BOOST_CONCEPT_ASSERT( (boost::CopyConstructibleConcept< Visitor >));
    BOOST_CONCEPT_ASSERT( (PRMVisitorConcept< Visitor, Graph, Topology >));
    BOOST_CONCEPT_ASSERT( (ADStarVisitorConcept< Visitor, Graph >));
  }
};

/**
  * This class is simply a "null" visitor for the FADPRM algorithm. It is null in the sense that it
  * will do nothing on all accounts.
  * \tparam Topology The topology type that represents the free-space.
  * \tparam PositionMap The property-map type which can store the position associated to each vertex.
  */
template < typename Topology, typename PositionMap >
class default_fadprm_visitor {
public:
  typedef typename ReaK::pp::topology_traits< Topology >::point_type PointType;

  default_fadprm_visitor( const Topology& free_space, PositionMap position )
      : m_free_space( free_space ), m_position( position ){};
  default_fadprm_visitor( const default_fadprm_visitor< Topology, PositionMap >& aVis )
      : m_free_space( aVis.m_free_space ), m_position( aVis.m_position ){};

  // AD* visitor functions:
  template < typename Vertex, typename Graph >
  void initialize_vertex( Vertex u, const Graph& g ) const {};
  template < typename Vertex, typename Graph >
  void discover_vertex( Vertex u, const Graph& g ) const {};
  template < typename Vertex, typename Graph >
  void examine_vertex( Vertex u, const Graph& g ) const {};
  template < typename Edge, typename Graph >
  void examine_edge( Edge e, const Graph& g ) const {};
  template < typename Graph >
  void publish_path( const Graph& g ) const {};
  template < typename EdgeIter, typename Graph >
  std::pair< double, EdgeIter > detect_edge_change( EdgeIter ei, const Graph& g ) const {
    RK_UNUSED( g );
    return std::pair< double, EdgeIter >( 0.0, ei );
  };
  template < typename Graph >
  double adjust_relaxation( double old_eps, const Graph& ) const {
    return old_eps;
  };


  // PRM visitor functions:
  template < typename Vertex, typename Graph >
  void vertex_added( Vertex, const Graph& ) const {};
  template < typename Edge, typename Graph >
  void edge_added( Edge, const Graph& ) const {};
  template < typename Vertex, typename Graph >
  std::pair< PointType, bool > random_walk( Vertex, const Graph& ) const {
    return std::make_pair( PointType(), false );
  };
  template < typename Vertex, typename Graph >
  void update_density( Vertex, const Graph& g ) const {};

  // Common to both PRM and AD* visitors:
  bool keep_going() const { return true; };

  const Topology& m_free_space;
  PositionMap m_position;
};


namespace detail {
namespace {

template < typename Topology, typename AStarHeuristicMap, typename UniformCostVisitor, typename UpdatableQueue,
           typename List, typename IndexInHeapMap, typename PredecessorMap, typename KeyMap, typename DistanceMap,
           typename RHSMap, typename WeightMap, typename DensityMap, typename PositionMap, typename NcSelector,
           typename ColorMap >
struct fadprm_bfs_visitor {

  typedef typename boost::property_traits< KeyMap >::value_type KeyValue;
  typedef typename boost::property_traits< ColorMap >::value_type ColorValue;
  typedef boost::color_traits< ColorValue > Color;
  typedef typename boost::property_traits< DistanceMap >::value_type distance_type;
  typedef typename boost::property_traits< WeightMap >::value_type weight_type;
  typedef typename boost::property_traits< PositionMap >::value_type PositionValue;

  fadprm_bfs_visitor( const Topology& super_space, AStarHeuristicMap h, UniformCostVisitor vis, UpdatableQueue& Q,
                      List& I, IndexInHeapMap index_in_heap, PredecessorMap p, KeyMap k, DistanceMap d, RHSMap rhs,
                      WeightMap w, DensityMap dens, PositionMap pos, NcSelector select_neighborhood, ColorMap col,
                      double& beta )
      : m_super_space( super_space ), m_h( h ), m_vis( vis ), m_Q( Q ), m_I( I ), m_index_in_heap( index_in_heap ),
        m_predecessor( p ), m_key( k ), m_distance( d ), m_rhs( rhs ), m_weight( w ), m_density( dens ),
        m_position( pos ), m_select_neighborhood( select_neighborhood ), m_color( col ), m_beta( beta ){};

  template < class Vertex, class Graph >
  void initialize_vertex( Vertex u, Graph& g ) const {
    m_vis.initialize_vertex( u, g );
  };
  template < class Vertex, class Graph >
  void discover_vertex( Vertex u, Graph& g ) const {
    m_vis.discover_vertex( u, g );
  };
  template < class Vertex, class Graph >
  void inconsistent_vertex( Vertex u, Graph& g ) const {};


  template < class Graph >
  typename boost::graph_traits< Graph >::vertex_descriptor create_vertex( const PositionValue& p, Graph& g ) const {
    typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
    typedef typename Graph::vertex_bundled VertexProp;

    VertexProp up;
    put( m_position, up, p );
    Vertex u = add_vertex( std::move( up ), g );
    m_vis.vertex_added( u, g );
    put( m_color, u, Color::white() );
    put( m_index_in_heap, u, static_cast< std::size_t >( -1 ) );
    put( m_distance, u, std::numeric_limits< double >::infinity() );
    put( m_rhs, u, std::numeric_limits< double >::infinity() );
    put( m_key, u, KeyValue( std::numeric_limits< double >::infinity(), std::numeric_limits< double >::infinity() ) );
    put( m_predecessor, u, u );

    return u;
  };

  template < typename Edge, typename Graph >
  void edge_added( Edge e, Graph& g ) const {
    m_vis.edge_added( e, g );
  };

  template < class Vertex, class Graph >
  void examine_vertex( Vertex u, Graph& g ) {
    typedef typename Graph::edge_bundled EdgeProp;

    m_vis.examine_vertex( u, g );

    prm_node_connector connect_vertex;

    std::size_t max_node_degree = 10;
    if( out_degree( u, g ) >= max_node_degree )
      return;
    max_node_degree = ( max_node_degree - out_degree( u, g ) ) / 2;
    for( std::size_t i = 0; i < max_node_degree; ++i ) {
      PositionValue p_rnd;
      bool expanding_worked;
      EdgeProp ep;
      boost::tie( p_rnd, expanding_worked, ep ) = m_vis.random_walk( u, g );
      if( expanding_worked )
        connect_vertex( p_rnd, u, ep, g, m_super_space, *this, m_position, m_select_neighborhood );
      else
        break;
    };
  };

  template < class Edge, class Graph >
  void examine_edge( Edge e, Graph& g ) const {
    if( get( m_weight, e ) < 0.0 )
      throw boost::negative_edge();
    m_vis.examine_edge( e, g );
  };
  template < class Vertex, class Graph >
  void forget_vertex( Vertex u, Graph& g ) const {};
  template < class Vertex, class Graph >
  void finish_vertex( Vertex u, Graph& g ) const {};
  template < class Vertex, class Graph >
  void recycle_vertex( Vertex u, Graph& g ) const {};

  template < typename Graph >
  void publish_path( const Graph& g ) const {
    m_vis.publish_path( g );
  };

  bool keep_going() const { return m_vis.keep_going(); };

  template < typename EdgeIter, typename Graph >
  std::pair< double, EdgeIter > detect_edge_change( EdgeIter ei, const Graph& g ) const {
    return m_vis.detect_edge_change( ei, g );
  };

  template < typename Graph >
  double adjust_epsilon( double old_eps, double, const Graph& g ) const {
    return m_vis.adjust_relaxation( old_eps, g );
  };

  template < class Vertex, typename Graph >
  void update_key( Vertex u, Graph& g ) const {
    m_vis.affected_vertex( u, g );
    distance_type g_u = get( m_distance, u );
    distance_type rhs_u = get( m_rhs, u );
    distance_type h_u = get( m_h, u );
    if( rhs_u < g_u ) {
      distance_type f_u = rhs_u + h_u; // 0.5 * ( rhs_u + h_u );
      put( m_key, u, KeyValue( ( 1.0 - m_beta ) * get( m_density, u ) + m_beta * f_u, rhs_u ) );
    } else {
      distance_type f_u = g_u + h_u; // 0.5 * ( g_u + h_u );
      put( m_key, u, KeyValue( f_u, rhs_u ) );
    };
  };

  template < class Vertex, class Graph >
  void update_vertex( Vertex u, Graph& g ) {
    BOOST_CONCEPT_ASSERT( (boost::BidirectionalGraphConcept< Graph >));
    typedef boost::graph_traits< Graph > GTraits;
    typename GTraits::in_edge_iterator ei, ei_end;

    ColorValue col_u = get( m_color, u );

    if( col_u == Color::white() ) {
      put( m_distance, u, std::numeric_limits< double >::infinity() );
      col_u = Color::green();
      put( m_color, u, col_u );
    };

    distance_type g_u = get( m_distance, u );
    distance_type rhs_u = get( m_rhs, u );

    if( rhs_u != 0.0 ) {
      rhs_u = std::numeric_limits< double >::infinity();
      typename GTraits::edge_descriptor pred_e;
      for( boost::tie( ei, ei_end ) = in_edges( u, g ); ei != ei_end; ++ei ) {
        distance_type rhs_tmp = get( m_weight, *ei ) + get( m_distance, source( *ei, g ) );
        if( rhs_tmp < rhs_u ) {
          rhs_u = rhs_tmp;
          put( m_rhs, u, rhs_u );
          put( m_predecessor, u, source( *ei, g ) );
          pred_e = *ei;
        };
      };
      rhs_u = get( m_rhs, u );
    };

    if( rhs_u != g_u ) {
      if( ( col_u != Color::black() ) && ( col_u != Color::red() ) ) {
        update_key( u, g );
        m_Q.push_or_update( u );
        put( m_color, u, Color::gray() );
        m_vis.discover_vertex( u, g );
      } else if( col_u == Color::black() ) {
        m_I.push_back( u );
        put( m_color, u, Color::red() );
      };
    } else if( m_Q.contains( u ) ) {
      put( m_key, u, KeyValue( 0.0, 0.0 ) );
      m_Q.update( u );
      m_Q.pop(); // remove from OPEN set
      update_key( u, g );
      put( m_color, u, Color::green() );
    };
  };


  template < typename Vertex, typename Graph >
  void travel_explored( Vertex u, Vertex v, Graph& g ) const {
    m_vis.travel_explored( u, v, g );
  };

  template < typename Vertex, typename Graph >
  void travel_succeeded( Vertex u, Vertex v, Graph& g ) const {
    m_vis.travel_succeeded( u, v, g );
  };

  template < typename Vertex, typename Graph >
  void travel_failed( Vertex u, Vertex v, Graph& g ) const {
    m_vis.travel_failed( u, v, g );
  };

  template < typename Vertex, typename Graph >
  void affected_vertex( Vertex u, Graph& g ) {
    update_key( u, g );
    update_vertex( u, g );
  };


  template < typename Vertex, typename Graph >
  std::pair< bool, typename Graph::edge_bundled > can_be_connected( Vertex u, Vertex v, Graph& g ) const {
    return m_vis.can_be_connected( u, v, g );
  };

  template < typename Vertex, typename Graph >
  boost::tuple< PositionValue, bool, typename Graph::edge_bundled > random_walk( Vertex u, Graph& g ) const {
    return m_vis.random_walk( u, g );
  };

  bool is_position_free( const PositionValue& p ) const { return m_vis.is_position_free( p ); };


  const Topology& m_super_space;
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
  distance_type& m_beta;
};
};
}; // detail


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the FADPRM algorithm (as of Belghith et al., 2006), without initialization of the
  * existing graph.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         BidirectionalGraphConcept and MutableGraphConcept.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
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
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
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
  * \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
  *        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
  */
template < typename Graph, typename Vertex, typename Topology, typename AStarHeuristicMap, typename FADPRMVisitor,
           typename PredecessorMap, typename DistanceMap, typename RHSMap, typename WeightMap, typename PositionMap,
           typename DensityMap, typename NcSelector, typename ColorMap >
inline void generate_fadprm_no_init( Graph& g, Vertex start_vertex, const Topology& free_space, AStarHeuristicMap hval,
                                     FADPRMVisitor vis, PredecessorMap predecessor, DistanceMap distance, RHSMap rhs,
                                     WeightMap weight, DensityMap density, PositionMap position,
                                     NcSelector select_neighborhood, ColorMap color, double epsilon ) {
  typedef adstar_key_value< double > KeyValue;
  typedef typename adstar_key_traits< KeyValue >::compare_type KeyCompareType;
  typedef boost::vector_property_map< KeyValue > KeyValueMap;
  KeyValueMap key_map;

  typedef boost::vector_property_map< std::size_t > IndexInHeapMap;
  IndexInHeapMap index_in_heap;
  {
    typename boost::graph_traits< Graph >::vertex_iterator ui, ui_end;
    for( boost::tie( ui, ui_end ) = vertices( g ); ui != ui_end; ++ui ) {
      put( index_in_heap, *ui, static_cast< std::size_t >( -1 ) );
      put( key_map, *ui,
           KeyValue( std::numeric_limits< double >::infinity(), std::numeric_limits< double >::infinity() ) );
    };
  };

  typedef boost::d_ary_heap_indirect< Vertex, 4, IndexInHeapMap, KeyValueMap, KeyCompareType > MutableQueue;
  MutableQueue Q( key_map, index_in_heap, KeyCompareType() ); // priority queue holding the OPEN set.
  std::vector< Vertex > I;                                    // list holding the INCONS set (inconsistent nodes).

  detail::fadprm_bfs_visitor< Topology, AStarHeuristicMap, FADPRMVisitor, MutableQueue, std::vector< Vertex >,
                              IndexInHeapMap, PredecessorMap, KeyValueMap, DistanceMap, RHSMap, WeightMap, DensityMap,
                              PositionMap, NcSelector, ColorMap > bfs_vis( free_space, hval, vis, Q, I, index_in_heap,
                                                                           predecessor, key_map, distance, rhs, weight,
                                                                           density, position, select_neighborhood,
                                                                           color, epsilon );

  detail::adstar_search_loop( g, start_vertex, hval, bfs_vis, predecessor, distance, rhs, key_map, weight, color,
                              index_in_heap, Q, I, epsilon, std::numeric_limits< double >::infinity(), 0.0,
                              std::less< double >(), std::equal_to< double >(), std::plus< double >(),
                              std::multiplies< double >() );
};


/**
* This function template generates a roadmap to connect a goal location to a start location
* using the FADPRM algorithm (as of Belghith et al., 2006), with initialization of the
* existing graph to (re)start the search.
* \tparam Graph The graph type that can store the generated roadmap, should model
*         BidirectionalGraphConcept and MutableGraphConcept.
* \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
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
*
* \param g A mutable graph that should initially store the starting
*        vertex (if not it will be randomly generated) and will store
*        the generated graph once the algorithm has finished.
* \param start_vertex The starting point of the algorithm, on the graph.
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
* \param epsilon The initial epsilon value that relaxes the A* search to give the AD* its anytime
*        characteristic. Epsilon values usually range from 1 to 10 (theoretically, the range is 1 to infinity).
*/
template < typename Graph, // this is the actual graph, should comply to BidirectionalMutableGraphConcept.
           typename Vertex, typename Topology,
           typename AStarHeuristicMap, // this the map of heuristic function value for each vertex.
           typename FADPRMVisitor,     // this is a visitor class that can perform special operations at event points.
           typename PredecessorMap,    // this is the map that stores the preceeding edge for each vertex.
           typename DistanceMap,       // this is the map of distance values associated with each vertex.
           typename RHSMap,
           typename WeightMap, // this is the map of edge weight (or cost) associated to each edge of the graph.
           typename DensityMap, typename PositionMap, typename NcSelector,
           typename ColorMap > // this is a color map for each vertex, i.e. white=not visited, gray=discovered,
// black=expanded.
inline void
  generate_fadprm( Graph& g, Vertex start_vertex, const Topology& free_space, AStarHeuristicMap hval, FADPRMVisitor vis,
                   PredecessorMap predecessor, DistanceMap distance, RHSMap rhs, WeightMap weight, DensityMap density,
                   PositionMap position, NcSelector select_neighborhood, ColorMap color, double epsilon ) {
  BOOST_CONCEPT_ASSERT( (boost::VertexListGraphConcept< Graph >));
  // BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT( (ReaK::pp::MetricSpaceConcept< Topology >));
  BOOST_CONCEPT_ASSERT( (ReaK::pp::PointDistributionConcept< Topology >));
  //     BOOST_CONCEPT_ASSERT((FADPRMVisitorConcept<FADPRMVisitor,Graph,Topology>));

  typedef typename boost::property_traits< ColorMap >::value_type ColorValue;
  typedef boost::color_traits< ColorValue > Color;
  typedef typename boost::property_traits< PositionMap >::value_type PositionValue;
  typename boost::graph_traits< Graph >::vertex_iterator ui, ui_end;
  typedef typename Graph::vertex_bundled VertexProp;

  if( num_vertices( g ) == 0 ) {
    VertexProp up;
    PositionValue p = get( ReaK::pp::random_sampler, free_space )( free_space );
    put( position, up, p );
    Vertex u = add_vertex( std::move( up ), g );
    vis.vertex_added( u, g );
    start_vertex = u;
  };

  for( boost::tie( ui, ui_end ) = vertices( g ); ui != ui_end; ++ui ) {
    put( color, *ui, Color::white() );
    put( distance, *ui, std::numeric_limits< double >::infinity() );
    put( rhs, *ui, std::numeric_limits< double >::infinity() );
    put( predecessor, *ui, *ui );
    vis.initialize_vertex( *ui, g );
  };

  generate_fadprm_no_init( g, start_vertex, free_space, hval, vis, predecessor, distance, rhs, weight, density,
                           position, select_neighborhood, color, epsilon );
};
};
};

#endif
