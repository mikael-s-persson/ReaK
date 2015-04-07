/**
 * \file path_planner_options_po.hpp
 *
 * This library defines functions to load path-planner options from Boost.Program-Options.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PATH_PLANNER_OPTIONS_PO_HPP
#define REAK_PATH_PLANNER_OPTIONS_PO_HPP

#include "path_planner_options.hpp"

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <boost/program_options.hpp>


namespace ReaK {

namespace pp {


boost::program_options::options_description get_planning_option_po_desc() {
  namespace po = boost::program_options;
  po::options_description planner_alg_options( "Planning algorithm options" );
  planner_alg_options.add_options()( "planner-options", po::value< std::string >(),
                                     "specify the file containing the planner-options data." )

    ( "planner-alg", po::value< std::string >(),
      "specify the planner algorithm to use, can be any of (rrt, rrt_star, prm, sba_star, fadprm)." )

    ( "max-vertices", po::value< std::size_t >()->default_value( 5000 ),
      "maximum number of vertices during runs (default is 5000)" )(
      "max-results", po::value< std::size_t >()->default_value( 50 ),
      "maximum number of result-paths during runs (default is 50)" )(
      "prog-interval", po::value< std::size_t >()->default_value( 10 ),
      "number of vertices between progress reports during runs (default is 10)" )

    ( "max-random-walk", po::value< double >()->default_value( 1.0 ),
      "specify the maximum random-walk distance allowed in the algorithm (default: 1.0). Only meaningful for expanding "
      "algorithms (SBA*, FADPRM, PRM)." )

    ( "no-lazy-connect", "if set, disable lazy connection strategy during planning." )(
      "bi-directional", "specify whether to use a bi-directional algorithm or not during planning. Only supported for "
                        "some algorithms (RRT, RRT*, SBA*)." )(
      "with-bnb", "specify whether to use a Branch-and-bound or not during planning to prune useless nodes from the "
                  "motion-graph. Only supported for optimizing algorithms (RRT*, SBA*)." )(
      "relaxation-factor", po::value< double >()->default_value( 0.0 ), "specify the initial relaxation factor for the "
                                                                        "algorithm (default: 0.0). Only supported for "
                                                                        "heuristic-driven algorithms." )(
      "start-delay", po::value< double >()->default_value( 0.0 ),
      "specify the starting time delay for the algorithm (default: 0.0). Only for dynamic problems (temporal space)." )
    //     ("density-cutoff", po::value< double >()->default_value(0.0), "specify the density cutoff (default: 0.0).
    //     Only supported for density-driven algorithms.")
    ( "with-voronoi-pull",
      "specify whether to use a Voronoi pull or not to add an exploratory bias to the search (default: not)." )(
      "sa-temperature", po::value< double >()->default_value( -1.0 ),
      "specify the initial Simulated Annealing temperature for algorithms that work on a exploration-exploitation "
      "schedule (e.g., SA-SBA*)." )

    ( "knn-method", po::value< std::string >(),
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
      "specify the KNN method to use (supported options: linear, bf2, bf4, cob2, cob4) (default: bf2)"
#else
      "specify the KNN method to use (supported options: linear, bf2, bf4) (default: bf2)"
#endif
      )

    ( "mg-storage", po::value< std::string >(),
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
      "specify the KNN method to use (supported options: adj-list, dvp-adj-list) (default: adj-list)"
#else
      "specify the KNN method to use (supported options: adj-list) (default: adj-list)"
#endif
      );
  return planner_alg_options;
};


planning_option_collection get_planning_option_from_po( boost::program_options::variables_map& vm ) {
  planning_option_collection plan_options;

  if( vm.count( "planner-options" ) ) {
    try {
      ( *serialization::open_iarchive( vm["planner-options"].as< std::string >() ) ) >> plan_options;
    } catch( std::exception& e ) {
      RK_UNUSED( e );
    };
  };

  if( vm.count( "planner-alg" ) ) {
    plan_options.planning_algo = 0; // RRT (default)
    if( vm["planner-alg"].as< std::string >() == "rrt_star" )
      plan_options.planning_algo = 1;
    else if( vm["planner-alg"].as< std::string >() == "prm" )
      plan_options.planning_algo = 2;
    else if( vm["planner-alg"].as< std::string >() == "sba_star" )
      plan_options.planning_algo = 3;
    else if( vm["planner-alg"].as< std::string >() == "fadprm" )
      plan_options.planning_algo = 4;
  };

  if( vm["max-vertices"].as< std::size_t >() != 5000 )
    plan_options.max_vertices = vm["max-vertices"].as< std::size_t >();
  if( vm["max-results"].as< std::size_t >() != 50 )
    plan_options.max_results = vm["max-results"].as< std::size_t >();
  if( vm["prog-interval"].as< std::size_t >() != 10 )
    plan_options.prog_interval = vm["prog-interval"].as< std::size_t >();

  if( vm.count( "knn-method" ) ) {
    plan_options.knn_method = 0;
    if( ( vm["knn-method"].as< std::string >() == "linear" ) && ( vm["mg-storage"].as< std::string >() == "adj-list" ) )
      plan_options.knn_method |= LINEAR_SEARCH_KNN;
    else if( vm["knn-method"].as< std::string >() == "bf4" )
      plan_options.knn_method |= DVP_BF4_TREE_KNN;
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
    else if( vm["knn-method"].as< std::string >() == "cob2" )
      plan_options.knn_method |= DVP_COB2_TREE_KNN;
    else if( vm["knn-method"].as< std::string >() == "cob4" )
      plan_options.knn_method |= DVP_COB4_TREE_KNN;
#endif
    else
      plan_options.knn_method |= DVP_BF2_TREE_KNN;
  };

  if( vm.count( "mg-storage" ) ) {
    plan_options.store_policy = 0;
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
    if( vm["mg-storage"].as< std::string >() == "dvp-adj-list" )
      plan_options.store_policy |= DVP_ADJ_LIST_MOTION_GRAPH;
    else
#endif
      plan_options.store_policy |= ADJ_LIST_MOTION_GRAPH;
  };

  if( !vm.count( "no-lazy-connect" ) )
    plan_options.planning_options |= LAZY_COLLISION_CHECKING; // use lazy, if supported.

  if( vm.count( "bi-directional" ) )
    plan_options.planning_options |= BIDIRECTIONAL_PLANNING;

  if( vm.count( "with-bnb" ) )
    plan_options.planning_options |= USE_BRANCH_AND_BOUND_PRUNING_FLAG;

  if( vm["relaxation-factor"].as< double >() > 1e-6 ) {
    plan_options.init_relax = vm["relaxation-factor"].as< double >();
    plan_options.planning_options |= PLAN_WITH_ANYTIME_HEURISTIC;
  };

  if( vm["sa-temperature"].as< double >() > -1.0 ) // default has been overriden
    plan_options.init_SA_temp = vm["sa-temperature"].as< double >();

  if( vm["max-random-walk"].as< double >() != 1.0 )
    plan_options.max_random_walk = vm["max-random-walk"].as< double >();

  if( vm.count( "with-voronoi-pull" ) )
    plan_options.planning_options |= PLAN_WITH_VORONOI_PULL;

  //   ("density-cutoff", po::value< double >()->default_value(0.0), "specify the density cutoff (default: 0.0). Only
  //   supported for density-driven algorithms.")

  return plan_options;
};
};
};

#endif
