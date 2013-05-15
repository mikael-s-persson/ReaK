/**
 * \file prm_path_planner.hpp
 * 
 * This library defines a class to solve path planning problems using the 
 * Probabilistic Road-map (PRM) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
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

#ifndef REAK_PRM_PATH_PLANNER_HPP
#define REAK_PRM_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"
#include "sbmp_reporter_concept.hpp"

#include "metric_space_concept.hpp"
#include "seq_path_wrapper.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "basic_sbmp_reporters.hpp"

#include "graph_alg/probabilistic_roadmap.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "dvp_layout_adjacency_list.hpp"
#include "metric_space_search.hpp"
#include "topological_search.hpp"
#include "path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/graph/astar_search.hpp>

#include <stack>

namespace ReaK {
  
namespace pp {
  


/**
 * This POD type contains the data required on a per-vertex basis for the PRM path-planning algorithm.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
struct prm_vertex_data {
  /// The position associated to the vertex.
  typename topology_traits<FreeSpaceType>::point_type position;  //for PRM
  /// The density associated to the vertex.
  double density;                                                //for PRM
  /// Keeps track of the root of the connected component this node belongs to.
  std::size_t cc_root;                                           //for PRM
  /// The travel-distance accumulated in the vertex, i.e., the travel-distance from the start vertex to this vertex.
  double distance_accum;                 //for A*
  /// The key-value associated to the vertex, computed by the A* algorithm.
  double astar_rhs_value;                //for A*
  /// The color-value associated to the vertex, computed by the A* algorithm.
  boost::default_color_type astar_color; //for A*
  /// The predecessor associated to the vertex, i.e., following the predecessor links starting at the goal node yields a backward trace of the optimal path.
  std::size_t predecessor;       //for A*
  
  /**
   * Default constructor.
   */
  prm_vertex_data() : position(typename topology_traits<FreeSpaceType>::point_type()),
                      density(0.0), distance_accum(0.0), astar_rhs_value(0.0), 
                      astar_color(), predecessor(0) { };
};

/**
 * This POD type contains the data required on a per-edge basis for the PRM path-planning algorithm.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
struct prm_edge_data { 
  /// The travel-distance associated to the edge (from source to target).
  double astar_weight; //for A*
  
  /**
   * Default constructor.
   * \param aWeight The travel-distance to be associated to this edge.
   */
  prm_edge_data(double aWeight = 0.0) : astar_weight(aWeight) { };
};


/**
 * This stateless functor type can be used to print out the information about a given PRM vertex.
 * This functor type can be used as a printing policy type for the vlist_sbmp_report class 
 * template that prints the list of vertices to a file.
 * \note This is mostly useful for debugging purposes (recording all information about the 
 *       motion-graph), it should not be used as the "output" of the path-planner.
 */
struct prm_vprinter : serialization::serializable {
  
  /**
   * This call operator prints all the PRM information about a given vertex 
   * to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the PRM planning algorithm.
   * \param out The output-stream to which to print the PRM information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  template <typename Vertex, typename Graph>
  void operator()(std::ostream& out, Vertex u, const Graph& g) const {
    using ReaK::to_vect;
    vect_n<double> v_pos = to_vect<double>(g[u].position);
    for(std::size_t i = 0; i < v_pos.size(); ++i)
      out << " " << std::setw(10) << v_pos[i];
    out << " " << std::setw(10) << g[u].density 
        << " " << std::setw(10) << g[u].cc_root 
        << " " << std::setw(10) << g[u].distance_accum 
        << " " << std::setw(10) << g[u].astar_rhs_value << std::endl;
  };
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };
  
  RK_RTTI_MAKE_ABSTRACT_1BASE(prm_vprinter,0xC2460015,1,"prm_vprinter",serialization::serializable)
};




/**
 * This class solves path planning problems using the 
 * Probabilistic Road-map (PRM) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename FreeSpaceType, 
          typename SBPPReporter = no_sbmp_report>
class prm_path_planner : public sample_based_planner< path_planner_base<FreeSpaceType> > {
  public:
    typedef sample_based_planner< path_planner_base<FreeSpaceType> > base_type;
    typedef prm_path_planner<FreeSpaceType, SBPPReporter> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    SBPPReporter m_reporter;
    point_type m_start_pos;
    point_type m_goal_pos;
    std::size_t max_num_results;
    bool has_reached_max_vertices;
    
    bool m_start_goal_connected;
    std::size_t m_v_count_at_connect;
    std::map<double, shared_ptr< seq_path_base< super_space_type > > > m_solutions;
    
  public:
    
    /**
     * Returns the best solution distance obtained after solve_path() has been called.
     * \return The best solution distance obtained after solve_path() has been called.
     */
    double get_best_solution_distance() const {
      if(m_solutions.size() == 0)
        return std::numeric_limits<double>::infinity();
      else
        return m_solutions.begin()->first;
    };
    
    /**
     * Returns true if the solver should keep on going trying to solve the path-planning problem.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \return True if the solver should keep on going trying to solve the path-planning problem.
     */
    bool keep_going() const {
      return (max_num_results > m_solutions.size()) && !has_reached_max_vertices;
    };
    
    /**
     * This function invokes the path-planning reporter to report on the progress of the path-planning
     * solver.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param g The current motion-graph.
     */
    template <typename Graph>
    void report_progress(Graph& g) {
      if(num_vertices(g) % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, get(&prm_vertex_data<FreeSpaceType>::position,g));
      has_reached_max_vertices = (num_vertices(g) >= this->m_max_vertex_count);
    };
    
    /**
     * This function computes the heuristic distance value of a given node, i.e., the bird-flight 
     * distance to the goal position.
     * \note This function is used internally by the path-planning algorithm (a visitor callback).
     * \param u The node for which the heuristic value is sought.
     * \param g The current motion-graph.
     * \return The heuristic distance value of the given node.
     */
    template <typename Graph>
    double heuristic(typename boost::graph_traits<Graph>::vertex_descriptor u, const Graph& g) const {
      return get(distance_metric, this->m_space->get_super_space())(g[u].position, this->m_goal_pos, this->m_space->get_super_space());
    };
    
    /**
     * This function constructs a solution path (if one is found) and invokes the path-planning 
     * reporter to report on that solution path.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param start_node The start node in the motion-graph.
     * \param goal_node The goal node in the motion-graph.
     * \param g The current motion-graph.
     */
    template <typename Vertex, typename Graph>
    void create_solution_path(Vertex start_node, Vertex goal_node, Graph& g) {
      
      if(!(m_start_goal_connected && ( (num_vertices(g) - m_v_count_at_connect) % (this->m_max_vertex_count / (2 * math::highest_set_bit(this->m_max_vertex_count))) == 0)))
        return;
      
      boost::astar_search(
        g, start_node,
        boost::bind(&prm_path_planner<FreeSpaceType,SBPPReporter>::heuristic<Graph>,this,_1,boost::cref(g)),
        boost::default_astar_visitor(),
        get(&prm_vertex_data<FreeSpaceType>::predecessor, g),
        get(&prm_vertex_data<FreeSpaceType>::astar_rhs_value, g),
        get(&prm_vertex_data<FreeSpaceType>::distance_accum, g),
        get(&prm_edge_data<FreeSpaceType>::astar_weight, g),
        boost::identity_property_map(),
        get(&prm_vertex_data<FreeSpaceType>::astar_color,g),
        std::less<double>(), std::plus<double>(),
        std::numeric_limits< double >::infinity(),
        double(0.0)); 
      
      double goal_distance = g[goal_node].distance_accum;
      
      if(goal_distance < std::numeric_limits<double>::infinity()) {
        //Draw the edges of the current best solution:
        
        shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
        shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("prm_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
        point_to_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
        std::set<Vertex> path;
        
        Vertex v = goal_node;
        point_type p_v = this->m_goal_pos;
        Vertex u = g[v].predecessor;
        point_type p_u = g[u].position;
        
        waypoints.push_front(m_goal_pos);
        waypoints.push_front(p_u);
        path.insert(v);
      
        while((u != start_node) && (path.insert(u).second)) {
          v = u; p_v = p_u;
          u = g[v].predecessor; 
          p_u = g[u].position; 
          waypoints.push_front(p_u);
        };
        
        if(u == start_node) {
          m_solutions[goal_distance] = new_sol;
          m_reporter.draw_solution(*(this->m_space), m_solutions[goal_distance]);
        };
      };
    };
    
    
    /**
     * This function constructs a solution path (if one is found) and invokes the path-planning 
     * reporter to report on that solution path.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param p_u The latest added position.
     * \param u The latest added node in the motion-graph.
     * \param g The current motion-graph.
     */
    template <typename Vertex, typename EdgeType, typename Graph>
    void check_goal_connection(Vertex start_node, Vertex goal_node, EdgeType e, Graph& g) {
      typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIter;
      
      // merge the connected components.
      std::stack<Vertex> trace_src;
      Vertex u = source(e,g);
      while(g[u].cc_root != u) {
        trace_src.push(u);
        u = g[u].cc_root;
      };
      trace_src.push(u);
      Vertex src_root = u;
      
      std::stack<Vertex> trace_dst;
      Vertex v = target(e,g);
      while(g[v].cc_root != v) {
        trace_dst.push(v);
        v = g[v].cc_root;
      };
      trace_dst.push(v);
      Vertex dst_root = v;
      
      Vertex common_root = src_root;
      if((dst_root == goal_node) || (src_root == goal_node))
        common_root = goal_node;
      else if((dst_root == start_node) || (src_root == start_node))
        common_root = start_node;
      else if(in_degree(target(e,g),g) > in_degree(source(e,g),g))
        common_root = dst_root;
      
      while( !trace_src.empty() ) {
        u = trace_src.top();
        g[u].cc_root = common_root;
        trace_src.pop();
      };
      
      while( !trace_dst.empty() ) {
        u = trace_dst.top();
        g[u].cc_root = common_root;
        trace_dst.pop();
      };
      
      //check if a start-goal connection was found:
      if( ((dst_root == start_node) && (src_root == goal_node)) ||
          ((dst_root == goal_node) && (src_root == start_node)) ) {
        m_start_goal_connected = true;
        m_v_count_at_connect = num_vertices(g);
      };
      
    };
    
    /**
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The path object that can be used to map out the path.
     */
    virtual shared_ptr< seq_path_base< super_space_type > > solve_path();
    
    /**
     * Returns a const-reference to the path-planning reporter used by this planner.
     * \return A const-reference to the path-planning reporter used by this planner.
     */
    const SBPPReporter& get_reporter() const { return m_reporter; };
    /**
     * Sets the path-planning reporter to be used by this planner.
     * \param aNewReporter The path-planning reporter to be used by this planner.
     */
    void set_reporter(const SBPPReporter& aNewReporter) { m_reporter = aNewReporter; };
    
    /**
     * Returns a const-reference to the start position used by this planner.
     * \return A const-reference to the start position used by this planner.
     */
    const point_type& get_start_pos() const { return m_start_pos; };
    /**
     * Sets the start position to be used by this planner.
     * \param aNewReporter The start position to be used by this planner.
     */
    void set_start_pos(const point_type& aStartPos) { m_start_pos = aStartPos; };
    
    /**
     * Returns a const-reference to the goal position used by this planner.
     * \return A const-reference to the goal position used by this planner.
     */
    const point_type& get_goal_pos() const { return m_goal_pos; };
    /**
     * Sets the goal position to be used by this planner.
     * \param aNewReporter The goal position to be used by this planner.
     */
    void set_goal_pos(const point_type& aGoalPos) { m_goal_pos = aGoalPos; };
    
    /**
     * Returns the maximum number of solutions that this planner should register.
     * \note Most probabilistic path-planners produce an initial solution that isn't perfect, and so, 
     *       more solutions should be sought to eventually have a more optimal one (shorter distance).
     * \return The maximum number of solutions that this planner should register.
     */
    std::size_t get_max_result_count() const { return max_num_results; };
    /**
     * Sets the maximum number of solutions that this planner should register.
     * \note Most probabilistic path-planners produce an initial solution that isn't perfect, and so, 
     *       more solutions should be sought to eventually have a more optimal one (shorter distance).
     * \param aMaxResultCount The maximum number of solutions that this planner should register.
     */
    void set_max_result_count(std::size_t aMaxResultCount) { max_num_results = aMaxResultCount; };
    
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aStartPos The position value of the starting location.
     * \param aGoalPos The position value of the goal location.
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     * \param aDataStructureFlags An integer flags representing the kind of motion graph data-structure to use in the 
     *                            planning algorithm. Can be ADJ_LIST_MOTION_GRAPH or DVP_ADJ_LIST_MOTION_GRAPH.
     *                            Any combination of those two and of KNN method flags to use for nearest
     *                            neighbor queries in the graph. KNN method flags can be LINEAR_SEARCH_KNN, 
     *                            DVP_BF2_TREE_KNN, DVP_BF4_TREE_KNN, DVP_COB2_TREE_KNN, or DVP_COB4_TREE_KNN.
     *                            See path_planner_options.hpp documentation.
     * \param aReporter The SBPP reporter object to use to report results and progress.
     * \param aMaxResultCount The maximum number of successful start-goal connections to make before 
     *                        stopping the path planner (the higher the number the more likely that a 
     *                        good path will be found, however, running time can become much longer).
     */
    prm_path_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                     const point_type& aStartPos = point_type(),
                     const point_type& aGoalPos = point_type(),
                     std::size_t aMaxVertexCount = 5000, 
                     std::size_t aProgressInterval = 100,
                     std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                     SBPPReporter aReporter = SBPPReporter(),
                     std::size_t aMaxResultCount = 50) :
                     base_type("prm_planner", aWorld, aMaxVertexCount, aProgressInterval, aDataStructureFlags, 0),
                     m_reporter(aReporter),
                     m_start_pos(aStartPos),
                     m_goal_pos(aGoalPos),
                     max_num_results(aMaxResultCount),
                     has_reached_max_vertices(false),
                     m_start_goal_connected(false),
                     m_v_count_at_connect(0) { };
    
    virtual ~prm_path_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(m_start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(m_start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results);
      has_reached_max_vertices = false;
      m_start_goal_connected = false;
      m_v_count_at_connect = 0;
      m_solutions.clear();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460008,1,"prm_path_planner",base_type)
};




/**
 * This class template is used by the PRM path-planner as the visitor object needed to 
 * collaborate with the PRM algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the PRMVisitorConcept and AStarVisitorConcept.
 * As with most planning algorithms in ReaK, the algorithm is really made up of a high-level 
 * algorithmic logic in the form of function templates, and a number of customization points 
 * collected as member functions of an algorithm visitor class that implement the problem-specific 
 * behaviors (random-walks / local-planning, heuristic computation, progress reporting, 
 * completion criteria, etc.). 
 */
template <typename FreeSpaceType, typename MotionGraph, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct prm_planner_visitor {
  typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
  
  shared_ptr< FreeSpaceType > m_space;
  prm_path_planner<FreeSpaceType,SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  Vertex m_start_node;
  Vertex m_goal_node;
  bool m_start_goal_connected;
  std::size_t m_v_count_at_connect;
  
  prm_planner_visitor(const shared_ptr< FreeSpaceType >& aSpace, 
                      prm_path_planner<FreeSpaceType,SBPPReporter>* aPlanner,
                      NNFinderSynchro aNNSynchro,
                      Vertex aStartNode, Vertex aGoalNode) : 
                      m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro),
                      m_start_node(aStartNode), m_goal_node(aGoalNode), 
                      m_start_goal_connected(false), m_v_count_at_connect(0) { };
  
  typedef typename topology_traits<FreeSpaceType>::point_type PointType;
  typedef prm_edge_data<FreeSpaceType> EdgeProp;
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    
    g[u].cc_root = u;
    
    g[u].astar_color = boost::color_traits<boost::default_color_type>::white();
    g[u].distance_accum = std::numeric_limits<double>::infinity();
    g[u].astar_rhs_value = std::numeric_limits<double>::infinity();
    g[u].predecessor = u;
    
    // Call progress reporter...
    m_planner->report_progress(g);
    
    // try to create a goal connection path
    m_planner->create_solution_path(m_start_node, m_goal_node, g);
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const {
    
    g[e].astar_weight = get(ReaK::pp::distance_metric, m_space->get_super_space())(
      g[source(e,g)].position,
      g[target(e,g)].position,
      m_space->get_super_space());
    
    m_planner->check_goal_connection(m_start_node, m_goal_node, e, g);
    
  };
  
  bool keep_going() const {
    return m_planner->keep_going();
  };
  
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, EdgeProp > random_walk(Vertex u, Graph& g) const {
    std::pair<PointType, bool> result = m_space->random_walk(g[u].position);
    if(result.second) {
      double dist = get(distance_metric, m_space->get_super_space())(g[u].position, result.first, m_space->get_super_space());
      return boost::tuple<PointType, bool, EdgeProp >(result.first, result.second, EdgeProp(dist));
    } else 
      return boost::tuple<PointType, bool, EdgeProp >(result.first, result.second, EdgeProp());
  };
  
  template <typename Vertex, typename Graph>
  std::pair<bool, EdgeProp > can_be_connected(Vertex u, Vertex v, const Graph& g) const {
    double dist = get(distance_metric, *m_space)(g[u].position, g[v].position, *m_space);
    return std::pair<bool, EdgeProp>((dist < std::numeric_limits<double>::infinity()), EdgeProp(dist));
  };
  
  bool is_position_free(const PointType& p) const {
    return m_space->is_free(p);
  };
  
  
  
  template <typename Vertex, typename Graph>
  void travel_explored(Vertex, Vertex, Graph&) const { };
  
  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex, Vertex, Graph&) const { };
  
  template <typename Vertex, typename Graph>
  void travel_failed(Vertex, Vertex, Graph&) const { };
  
  
  template <typename Vertex, typename Graph>
#if 0
  // Old connection strategy:
  void update_density(Vertex u, Graph& g) const {
#endif
  void affected_vertex(Vertex u, Graph& g) const {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    //take the sum of all weights of outgoing edges.
    if(out_degree(u,g) == 0) {
      g[u].density = 0.0;
      return;
    };
    double sum = 0.0;
    OutEdgeIter ei, ei_end;
    for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei) {
      sum += g[*ei].astar_weight;
    };
    sum /= out_degree(u,g) * out_degree(u,g); //this gives the average edge distances divided by the number of adjacent nodes.
    g[u].density = std::exp(-sum*sum);
  };
  
  
};






template <typename FreeSpaceType, 
          typename SBPPReporter>
shared_ptr< seq_path_base< typename prm_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > > 
  prm_path_planner<FreeSpaceType,SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  this->has_reached_max_vertices = false;
  this->m_start_goal_connected = false;
  this->m_v_count_at_connect = 0;
  this->m_solutions.clear();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef boost::data_member_property_map<PointType, prm_vertex_data<FreeSpaceType> > PositionMap;
  PositionMap pos_map = PositionMap(&prm_vertex_data<FreeSpaceType>::position);
  
  typedef boost::data_member_property_map<double, prm_vertex_data<FreeSpaceType> > DensityMap;
  DensityMap dens_map = DensityMap(&prm_vertex_data<FreeSpaceType>::density);
  
  typedef boost::data_member_property_map<std::size_t, prm_vertex_data<FreeSpaceType> > CCRootMap;
  CCRootMap cc_root_map = CCRootMap(&prm_vertex_data<FreeSpaceType>::cc_root);
  
  double space_dim = double((to_vect<double>(this->m_space->get_super_space().difference(this->m_goal_pos,this->m_start_pos))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  

#ifdef RK_ENABLE_CXX0X_FEATURES

#define RK_PRM_INITIALIZE_START_AND_GOAL \
    prm_vertex_data<FreeSpaceType> vs_p, vg_p; \
    vs_p.position = this->m_start_pos; \
    vg_p.position = this->m_goal_pos; \
    Vertex start_node = add_vertex(std::move(vs_p), motion_graph); \
    Vertex goal_node = add_vertex(std::move(vg_p), motion_graph); \
    motion_graph[start_node].cc_root = start_node; \
    motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white(); \
    motion_graph[start_node].distance_accum = 0.0; \
    motion_graph[start_node].astar_rhs_value = this->heuristic(start_node, motion_graph); \
    motion_graph[start_node].predecessor = start_node; \
    motion_graph[goal_node].cc_root = goal_node; \
    motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white(); \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].astar_rhs_value = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node;
    
#else
    
#define RK_PRM_INITIALIZE_START_AND_GOAL \
    prm_vertex_data<FreeSpaceType> vs_p, vg_p; \
    vs_p.position = this->m_start_pos; \
    vg_p.position = this->m_goal_pos; \
    Vertex start_node = add_vertex(vs_p, motion_graph); \
    Vertex goal_node = add_vertex(vg_p, motion_graph); \
    motion_graph[start_node].cc_root = start_node; \
    motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white(); \
    motion_graph[start_node].distance_accum = 0.0; \
    motion_graph[start_node].astar_rhs_value = this->heuristic(start_node, motion_graph); \
    motion_graph[start_node].predecessor = start_node; \
    motion_graph[goal_node].cc_root = goal_node; \
    motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white(); \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].astar_rhs_value = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node;
    
#endif
  
  
#define RK_PRM_MAKE_GENERATE_PRM_CALL \
  ReaK::graph::generate_prm(motion_graph, this->m_space->get_super_space(), \
                            vis, pos_map, get(random_sampler, this->m_space->get_super_space()), \
                            dens_map, cc_root_map, \
                            ReaK::graph::star_neighborhood< NNFinderType >( \
                              nn_finder, \
                              space_dim, 3.0 * space_Lc), \
                            this->m_max_vertex_count, \
                            0.2);
  
  
  
  if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
    
    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                           prm_vertex_data<FreeSpaceType>,
                           prm_edge_data<FreeSpaceType>,
                           boost::listS> MotionGraphType;
    typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
    typedef typename MotionGraphType::vertex_property_type VertexProp;
    typedef boost::composite_property_map< 
    PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
    
    MotionGraphType motion_graph;
    GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
    
    RK_PRM_INITIALIZE_START_AND_GOAL
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
      prm_planner_visitor<FreeSpaceType, MotionGraphType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node);
      
      typedef linear_neighbor_search<> NNFinderType;
      NNFinderType nn_finder;
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    };
    
  } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        prm_vertex_data<FreeSpaceType>,
        prm_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_PRM_INITIALIZE_START_AND_GOAL
      
      typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        prm_vertex_data<FreeSpaceType>,
        prm_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_PRM_INITIALIZE_START_AND_GOAL
      
      typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        prm_vertex_data<FreeSpaceType>,
        prm_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_PRM_INITIALIZE_START_AND_GOAL
      
      typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        prm_vertex_data<FreeSpaceType>,
        prm_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_PRM_INITIALIZE_START_AND_GOAL
      
      typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
      NNFinderType nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      prm_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node);
      
      RK_PRM_MAKE_GENERATE_PRM_CALL
      
    };
    
  };
  
#undef RK_PRM_INITIALIZE_START_AND_GOAL
#undef RK_PRM_MAKE_GENERATE_PRM_CALL
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< seq_path_base< SuperSpace > >();
};


};

};

#endif

