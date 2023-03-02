/**
 * \file simple_graph_traits.hpp
 *
 * This library contains the
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2023
 */

/*
 *    Copyright 2023 Sven Mikael Persson
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

#ifndef REAK_SIMPLE_GRAPH_TRAITS_HPP
#define REAK_SIMPLE_GRAPH_TRAITS_HPP

#include <boost/graph/graph_traits.hpp>

namespace ReaK::graph {

template <typename Graph>
using graph_vertex_t = typename boost::graph_traits<Graph>::vertex_descriptor;

template <typename Graph>
using graph_edge_t = typename boost::graph_traits<Graph>::edge_descriptor;

template <typename Graph>
using graph_vertex_bundle_t = std::decay_t<
    decltype(std::declval<Graph>()[std::declval<graph_vertex_t<Graph>>()])>;

template <typename Graph>
using graph_edge_bundle_t = std::decay_t<
    decltype(std::declval<Graph>()[std::declval<graph_edge_t<Graph>>()])>;

template <typename Graph>
using graph_vertex_property_t = typename Graph::vertex_property_type;

template <typename Graph>
using graph_edge_property_t = typename Graph::edge_property_type;

template <typename PropertyMap>
using property_value_t =
    typename boost::property_traits<PropertyMap>::value_type;

}  // namespace ReaK::graph

#endif  // REAK_SIMPLE_GRAPH_TRAITS_HPP
