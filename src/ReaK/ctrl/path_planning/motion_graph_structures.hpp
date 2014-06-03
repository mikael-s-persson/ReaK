/**
 * \file motion_graph_structures.hpp
 * 
 * This library defines a class to solve path planning problems using the 
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_MOTION_GRAPH_STRUCTURES_HPP
#define REAK_MOTION_GRAPH_STRUCTURES_HPP


// #define RK_PLANNERS_ENABLE_VEBL_TREE
// #define RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT


#include <ReaK/core/base/defs.hpp>

// BGL-Extra includes:
#include <boost/graph/bfl_d_ary_tree.hpp>
#include <boost/graph/adjacency_list_BC.hpp>

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
#include <boost/graph/vebl_d_ary_tree.hpp>
#endif

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
#include "dvp_layout_adjacency_list.hpp"
#endif



#include <ReaK/ctrl/topologies/metric_space_concept.hpp>

#include <boost/mpl/if.hpp>


namespace ReaK {
  
namespace pp {

template <typename FreeSpaceType>
struct motion_segment_directionality {
  typedef typename boost::mpl::if_<
    is_metric_symmetric< FreeSpaceType >,
    boost::undirectedS,
    boost::bidirectionalS >::type type;
};

};

};



#endif

