/**
 * \file path_planner_options.hpp
 * 
 * This library defines the options available when creating a path-planner.
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

#ifndef REAK_PATH_PLANNER_OPTIONS_HPP
#define REAK_PATH_PLANNER_OPTIONS_HPP

namespace ReaK {
  
namespace pp {


const std::size_t ADJ_LIST_MOTION_GRAPH = 0;
const std::size_t DVP_ADJ_LIST_MOTION_GRAPH = 1;
const std::size_t LINKED_TREE_MOTION_GRAPH = 2;

const std::size_t LINEAR_SEARCH_KNN = 0;
const std::size_t APPROX_LINEAR_SEARCH_KNN = 1;
const std::size_t DVP_LINKED_TREE_KNN = 2;
const std::size_t DVP_BF2_TREE_KNN = 3;
const std::size_t DVP_BF4_TREE_KNN = 4;
const std::size_t DVP_COB2_TREE_KNN = 5;
const std::size_t DVP_COB4_TREE_KNN = 6;
const std::size_t DVP_ALT_BF2_KNN = 7;
const std::size_t DVP_ALT_BF4_KNN = 8;
const std::size_t DVP_ALT_COB2_KNN = 9;
const std::size_t DVP_ALT_COB4_KNN = 10;


};

};

#endif

