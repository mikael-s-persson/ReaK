/**
 * \file bgl_raw_property_graph.hpp
 *
 * This library provides traits and concepts definitions for a RawPropertyGraphConcept.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
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

#ifndef REAK_BGL_RAW_PROPERTY_GRAPH_HPP
#define REAK_BGL_RAW_PROPERTY_GRAPH_HPP

#include "boost/graph/adjacency_list_BC.hpp"
#include "boost/graph/graph_concepts.hpp"
#include "boost/graph/properties.hpp"

// BGL-Extra includes:
#include "boost/graph/more_property_maps.hpp"

#include <type_traits>

namespace boost {

enum vertex_raw_prop_to_bundle_t { vertex_raw_prop_to_bundle };
enum vertex_raw_property_t { vertex_raw_property };

BOOST_INSTALL_PROPERTY(vertex, raw_prop_to_bundle);
BOOST_INSTALL_PROPERTY(vertex, raw_property);

enum edge_raw_prop_to_bundle_t { edge_raw_prop_to_bundle };
enum edge_raw_property_t { edge_raw_property };

BOOST_INSTALL_PROPERTY(edge, raw_prop_to_bundle);
BOOST_INSTALL_PROPERTY(edge, raw_property);

template <typename Graph>
struct raw_vertex_property_map {
  using type = whole_bundle_property_map<Graph, vertex_bundle_t>;
  using const_type = whole_bundle_property_map<const Graph, vertex_bundle_t>;
};

template <typename Graph>
struct raw_edge_property_map {
  using type = whole_bundle_property_map<Graph, edge_bundle_t>;
  using const_type = whole_bundle_property_map<const Graph, edge_bundle_t>;
};

template <typename Graph>
struct raw_vertex_to_bundle_map {
  using type = self_property_map<typename Graph::vertex_bundled>;
  using const_type = self_property_map<const typename Graph::vertex_bundled>;
};

template <typename Graph>
struct raw_edge_to_bundle_map {
  using type = self_property_map<typename Graph::edge_bundled>;
  using const_type = self_property_map<const typename Graph::edge_bundled>;
};

template <typename Graph>
struct RawPropertyGraphConcept {
  typename graph_traits<Graph>::vertex_descriptor v;
  typename graph_traits<Graph>::edge_descriptor e;
  Graph g;

  using vertex_raw_prop_to_bundle_map =
      typename raw_vertex_to_bundle_map<Graph>::type;
  using vertex_raw_prop_map = typename raw_vertex_property_map<Graph>::type;
  using edge_raw_prop_to_bundle_map =
      typename raw_edge_to_bundle_map<Graph>::type;
  using edge_raw_prop_map = typename raw_edge_property_map<Graph>::type;

  vertex_raw_prop_to_bundle_map vrp_to_bundle;
  vertex_raw_prop_map vrp_map;
  edge_raw_prop_to_bundle_map erp_to_bundle;
  edge_raw_prop_map erp_map;

  using v_prop_type =
      typename property_traits<vertex_raw_prop_to_bundle_map>::key_type;
  using v_bundle_type =
      typename property_traits<vertex_raw_prop_to_bundle_map>::value_type;
  using e_prop_type =
      typename property_traits<vertex_raw_prop_to_bundle_map>::key_type;
  using e_bundle_type =
      typename property_traits<vertex_raw_prop_to_bundle_map>::value_type;

  v_prop_type vrp;
  v_bundle_type vb;
  e_prop_type erp;
  e_bundle_type eb;

  BOOST_CONCEPT_ASSERT((boost::GraphConcept<Graph>));

  BOOST_CONCEPT_USAGE(RawPropertyGraphConcept) {
    vrp_to_bundle = get_raw_vertex_to_bundle_map(g);
    vrp_map = get_raw_vertex_property_map(g);
    erp_to_bundle = get_raw_edge_to_bundle_map(g);
    erp_map = get_raw_edge_property_map(g);

    vb = get_raw_vertex_to_bundle(g, vrp);
    vrp = get_raw_vertex_property(g, v);
    eb = get_raw_edge_to_bundle(g, erp);
    erp = get_raw_edge_property(g, e);

    put_raw_vertex_to_bundle(g, vrp, vb);
    put_raw_vertex_property(g, v, vrp);
    put_raw_edge_to_bundle(g, erp, eb);
    put_raw_edge_property(g, e, erp);

    vb = get(vrp_to_bundle, vrp);
    vrp = get(vrp_map, v);
    eb = get(erp_to_bundle, erp);
    erp = get(erp_map, e);

    put(vrp_to_bundle, vrp, vb);
    put(vrp_map, v, vrp);
    put(erp_to_bundle, erp, eb);
    put(erp_map, e, erp);
  }
};

template <typename Graph>
auto get_raw_vertex_property_map(Graph& g) {
  return typename raw_vertex_property_map<Graph>::type(&g);
}

template <typename Graph>
auto get_raw_vertex_property_map(const Graph& g) {
  return typename raw_vertex_property_map<Graph>::const_type(&g);
}

template <typename Graph>
auto get_raw_edge_property_map(Graph& g) {
  return typename raw_edge_property_map<Graph>::type(&g);
}

template <typename Graph>
auto get_raw_edge_property_map(const Graph& g) {
  return typename raw_edge_property_map<Graph>::const_type(&g);
}

template <typename Graph>
auto get_raw_vertex_to_bundle_map(Graph& g) {
  return typename raw_vertex_to_bundle_map<Graph>::type();
}

template <typename Graph>
auto get_raw_vertex_to_bundle_map(const Graph& g) {
  return typename raw_vertex_to_bundle_map<Graph>::const_type();
}

template <typename Graph>
auto get_raw_edge_to_bundle_map(Graph& g) {
  return typename raw_edge_to_bundle_map<Graph>::type();
}

template <typename Graph>
auto get_raw_edge_to_bundle_map(const Graph& g) {
  return typename raw_edge_to_bundle_map<Graph>::const_type();
}

template <typename Graph, typename Vertex>
auto& get_raw_vertex_property(Graph& g, Vertex v) {
  return g[v];
}

template <typename Graph, typename Vertex>
const auto& get_raw_vertex_property(const Graph& g, Vertex v) {
  return g[v];
}

template <typename Graph, typename Edge>
auto& get_raw_edge_property(Graph& g, Edge e) {
  return g[e];
}

template <typename Graph, typename Edge>
const auto& get_raw_edge_property(const Graph& g, Edge e) {
  return g[e];
}

template <typename Graph, typename Bundle>
Bundle& get_raw_vertex_to_bundle(const Graph& /*unused*/, Bundle& b) {
  return b;
}

template <typename Graph, typename Bundle>
Bundle& get_raw_edge_to_bundle(const Graph& /*unused*/, Bundle& b) {
  return b;
}

template <typename Graph, typename Vertex, typename Bundle>
void put_raw_vertex_property(Graph& g, Vertex v, Bundle&& value) {
  g[v] = std::forward<Bundle>(value);
}

template <typename Graph, typename Edge, typename Bundle>
void put_raw_edge_property(Graph& g, Edge e, Bundle&& value) {
  g[e] = std::forward<Bundle>(value);
}

template <typename Graph, typename Bundle1, typename Bundle2>
void put_raw_vertex_to_bundle(const Graph& /*unused*/, Bundle1& b1,
                              Bundle2&& b2) {
  b1 = std::forward<Bundle2>(b2);
}

template <typename Graph, typename Bundle1, typename Bundle2>
void put_raw_edge_to_bundle(const Graph& /*unused*/, Bundle1& b1,
                            Bundle2&& b2) {
  b1 = std::forward<Bundle2>(b2);
}

/**
 * This property-map uses a graph's "get_raw_vertex_property" function to obtain the
 * property value associated to a given descriptor.
 * \tparam T The value-type of the property.
 * \tparam Graph The graph type.
 */
template <typename T, typename Graph>
struct raw_vertex_propgraph_map
    : public put_get_helper<T&, raw_vertex_propgraph_map<T, Graph>> {
 private:
  Graph* pg;

 public:
  using is_vertex_prop = std::true_type;
  using is_const_graph = is_const<Graph>;
  using value_type = T;
  using reference = T&;

  using key_type = typename graph_traits<Graph>::vertex_descriptor;
  using category =
      std::conditional_t<std::is_const_v<T>, readable_property_map_tag,
                         lvalue_property_map_tag>;

  explicit raw_vertex_propgraph_map(Graph* aPG = nullptr) : pg(aPG) {}
  reference operator[](key_type k) const {
    return get_raw_vertex_property(*pg, k);
  }
};

/**
 * This property-map uses a graph's "get_raw_edge_property" function to obtain the
 * property value associated to a given descriptor.
 * \tparam T The value-type of the property.
 * \tparam Graph The graph type.
 */
template <typename T, typename Graph>
struct raw_edge_propgraph_map
    : public put_get_helper<T&, raw_edge_propgraph_map<T, Graph>> {
 private:
  Graph* pg;

 public:
  using is_vertex_prop = std::false_type;
  using is_const_graph = is_const<Graph>;
  using value_type = T;
  using reference = T&;

  using key_type = typename graph_traits<Graph>::edge_descriptor;
  using category =
      std::conditional_t<std::is_const_v<T>, readable_property_map_tag,
                         lvalue_property_map_tag>;

  explicit raw_edge_propgraph_map(Graph* aPG = nullptr) : pg(aPG) {}
  reference operator[](key_type k) const {
    return get_raw_edge_property(*pg, k);
  }
};

}  // namespace boost

#endif
