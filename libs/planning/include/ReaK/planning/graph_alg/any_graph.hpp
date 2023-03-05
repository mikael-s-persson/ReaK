/**
 * \file any_graph.hpp
 *
 * This library provides a type-erasure base-class to represent any graph.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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
 *    If not, sMOTION_ee <http://www.gnu.org/licenses/>.
 */

#ifndef RK_ANY_GRAPH_HPP
#define RK_ANY_GRAPH_HPP

#include <stdexcept>
#include <string>

#include "boost/graph/graph_concepts.hpp"
#include "boost/graph/graph_mutability_traits.hpp"
#include "boost/graph/properties.hpp"
#include "boost/iterator/iterator_facade.hpp"

#include <any>
#include "boost/range/any_range.hpp"

#include <tuple>
#include <type_traits>

#include "ReaK/planning/graph_alg/simple_graph_traits.hpp"

// BGL-Extra includes:
#include "boost/graph/more_property_tags.hpp"

namespace ReaK::graph {

class any_graph {
 public:
  using any_descriptor = std::any;

  using directed_category = void;
  using edge_parallel_category = void;
  using traversal_category = void;

  struct vertex_descriptor {
    any_descriptor base;
    explicit vertex_descriptor(const any_descriptor& aBase) : base(aBase) {}
    vertex_descriptor() : vertex_descriptor(any_descriptor()) {}
    operator any_descriptor&() { return base; }              // NOLINT
    operator const any_descriptor&() const { return base; }  // NOLINT
  };

  struct edge_descriptor {
    any_descriptor base;
    explicit edge_descriptor(const any_descriptor& aBase) : base(aBase) {}
    edge_descriptor() : edge_descriptor(any_descriptor()) {}
    operator any_descriptor&() { return base; }              // NOLINT
    operator const any_descriptor&() const { return base; }  // NOLINT
  };

  using any_bundled = std::any;
  using vertex_bundled = any_bundled;
  using edge_bundled = any_bundled;

  using edge_range =
      boost::any_range<edge_descriptor, boost::forward_traversal_tag,
                       edge_descriptor, std::ptrdiff_t>;
  using out_edge_range = edge_range;
  using in_edge_range = edge_range;

  using edge_iterator = edge_range::iterator;
  using out_edge_iterator = edge_iterator;
  using in_edge_iterator = edge_iterator;

  using vertex_range =
      boost::any_range<vertex_descriptor, boost::forward_traversal_tag,
                       vertex_descriptor, std::ptrdiff_t>;

  using vertex_iterator = vertex_range::iterator;
  using adjacency_iterator = vertex_range::iterator;

  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;

 protected:
  virtual void* get_property_by_ptr(std::string_view aProperty,
                                    const std::any& aElement) const = 0;
  virtual std::any get_property_by_any(std::string_view aProperty,
                                       const std::any& aElement) const = 0;

  virtual vertex_bundled get_vertex_bundled(
      const vertex_descriptor& aVertex) const = 0;
  virtual edge_bundled get_edge_bundled(const edge_descriptor& aEdge) const = 0;

  virtual bool equal_vertex_descriptors(const vertex_descriptor& aU,
                                        const vertex_descriptor& aV) const = 0;
  virtual bool equal_edge_descriptors(const edge_descriptor& aE,
                                      const edge_descriptor& aF) const = 0;

  virtual edge_range edges_mem() const = 0;
  virtual edges_size_type num_edges_mem() const = 0;

  virtual vertex_descriptor source_mem(const edge_descriptor& aEdge) const = 0;
  virtual vertex_descriptor target_mem(const edge_descriptor& aEdge) const = 0;

  virtual std::pair<edge_descriptor, bool> edge_mem(
      const vertex_descriptor& aU, const vertex_descriptor& aV) const = 0;

  virtual out_edge_range out_edges_mem(
      const vertex_descriptor& aVertex) const = 0;
  virtual degree_size_type out_degree_mem(
      const vertex_descriptor& aVertex) const = 0;

  virtual in_edge_range in_edges_mem(
      const vertex_descriptor& aVertex) const = 0;
  virtual degree_size_type in_degree_mem(
      const vertex_descriptor& aVertex) const = 0;
  virtual degree_size_type degree_mem(
      const vertex_descriptor& aVertex) const = 0;

  virtual vertex_range vertices_mem() const = 0;
  virtual vertices_size_type num_vertices_mem() const = 0;

  virtual vertex_descriptor add_vertex_mem() = 0;
  virtual void clear_vertex_mem(const vertex_descriptor& aVertex) = 0;
  virtual void remove_vertex_mem(const vertex_descriptor& aVertex) = 0;

  virtual std::pair<edge_descriptor, bool> add_edge_mem(
      const vertex_descriptor& aU, const vertex_descriptor& aV) = 0;
  virtual void remove_edge_mem(const edge_descriptor& aEdge) = 0;

 public:
  vertex_bundled operator[](const vertex_descriptor& aVertex) const {
    return get_vertex_bundled(aVertex);
  }
  edge_bundled operator[](const edge_descriptor& aEdge) const {
    return get_edge_bundled(aEdge);
  }

  bool equal_descriptors(const vertex_descriptor& aU,
                         const vertex_descriptor& aV) const {
    return equal_vertex_descriptors(aU, aV);
  }
  bool equal_descriptors(const edge_descriptor& aE,
                         const edge_descriptor& aF) const {
    return equal_edge_descriptors(aE, aF);
  }

  virtual ~any_graph() = default;

  /*******************************************************************************************
   *                  IncidenceGraph concept
   ******************************************************************************************/

  friend inline auto out_edges(vertex_descriptor v, const any_graph& g) {
    using std::begin;
    using std::end;
    out_edge_range er = g.out_edges_mem(v);
    return std::pair(begin(er), end(er));
  }

  friend inline auto source(edge_descriptor e, const any_graph& g) {
    return g.source_mem(e);
  }

  friend inline auto target(edge_descriptor e, const any_graph& g) {
    return g.target_mem(e);
  }

  friend inline auto out_degree(vertex_descriptor v, const any_graph& g) {
    return g.out_degree_mem(v);
  }

  /*******************************************************************************************
   *                  BidirectionalGraph concept
   ******************************************************************************************/

  friend inline auto in_edges(vertex_descriptor v, const any_graph& g) {
    using std::begin;
    using std::end;
    in_edge_range er = g.in_edges_mem(v);
    return std::pair(begin(er), end(er));
  }

  friend inline auto in_degree(vertex_descriptor v, const any_graph& g) {
    return g.in_degree_mem(v);
  }

  friend inline auto degree(vertex_descriptor v, const any_graph& g) {
    return g.degree_mem(v);
  }

  /*******************************************************************************************
   *                  VertexListGraph concept
   ******************************************************************************************/

  friend inline auto vertices(const any_graph& g) {
    using std::begin;
    using std::end;
    vertex_range vr = g.vertices_mem();
    return std::pair(begin(vr), end(vr));
  }

  friend inline auto num_vertices(const any_graph& g) {
    return g.num_vertices_mem();
  }

  /*******************************************************************************************
   *                  EdgeListGraph concept
   ******************************************************************************************/

  friend inline auto edges(const any_graph& g) {
    using std::begin;
    using std::end;
    edge_range er = g.edges_mem();
    return std::pair(begin(er), end(er));
  }

  friend inline auto num_edges(const any_graph& g) { return g.num_edges_mem(); }

  /*******************************************************************************************
   *                  AdjacencyMatrix concept
   ******************************************************************************************/

  friend inline auto edge(vertex_descriptor u, vertex_descriptor v,
                          const any_graph& g) {
    return g.edge_mem(u, v);
  }

  /*******************************************************************************************
   *                  MutableGraph concept
   ******************************************************************************************/

  friend inline auto add_vertex(any_graph& g) { return g.add_vertex_mem(); }

  friend inline void clear_vertex(vertex_descriptor v, any_graph& g) {
    g.clear_vertex_mem(v);
  }

  friend inline void remove_vertex(vertex_descriptor v, any_graph& g) {
    g.remove_vertex_mem(v);
  }

  friend inline auto add_edge(vertex_descriptor u, vertex_descriptor v,
                              any_graph& g) {
    return g.add_edge_mem(u, v);
  }

  friend inline void remove_edge(vertex_descriptor u, vertex_descriptor v,
                                 any_graph& g) {
    auto [e, e_exists] = g.edge_mem(u, v);
    if (e_exists) {
      g.remove_edge_mem(e);
    }
  }

  friend inline void remove_edge(edge_descriptor e, any_graph& g) {
    g.remove_edge_mem(e);
  }

  friend inline void remove_edge(edge_iterator e_iter, any_graph& g) {
    g.remove_edge_mem(*e_iter);
  }

  template <typename ValueType>
  struct property_map_by_ptr
      : public boost::put_get_helper<ValueType&,
                                     property_map_by_ptr<ValueType>> {
    using self = property_map_by_ptr<ValueType>;

    using value_type = ValueType;
    using reference = ValueType&;
    using key_type = std::any;
    using category = std::conditional_t<std::is_const_v<ValueType>,
                                        boost::readable_property_map_tag,
                                        boost::lvalue_property_map_tag>;

    const any_graph* p_parent;
    std::string_view prop_name;

    property_map_by_ptr(const any_graph* aParentPtr, std::string_view aProperty)
        : p_parent(aParentPtr), prop_name(aProperty) {}
    reference operator[](const key_type& k) const {
      return *(
          static_cast<ValueType*>(p_parent->get_property_by_ptr(prop_name, k)));
    }
  };

  template <typename ValueType>
  struct property_map_by_any
      : public boost::put_get_helper<ValueType,
                                     property_map_by_any<ValueType>> {
    using self = property_map_by_any<ValueType>;

    using value_type = ValueType;
    using reference = ValueType;
    using key_type = std::any;
    using category = boost::readable_property_map_tag;

    const any_graph* p_parent;
    std::string_view prop_name;

    property_map_by_any(const any_graph* aParentPtr, std::string_view aProperty)
        : p_parent(aParentPtr), prop_name(aProperty) {}
    reference operator[](const key_type& k) const {
      return std::any_cast<ValueType>(
          p_parent->get_property_by_any(prop_name, k));
    }
  };
};

template <typename ValueType>
auto get_dyn_prop(const std::string& aProperty, const any_graph& aGraph) {
  if constexpr (std::is_reference_v<ValueType>) {
    return any_graph::property_map_by_ptr<std::remove_reference_t<ValueType>>(
        &aGraph, aProperty);
  } else {
    return any_graph::property_map_by_any<ValueType>(&aGraph, aProperty);
  }
}

template <typename ValueType>
ValueType get_dyn_prop(const std::string& aProperty, std::any aDesc,
                       const any_graph& aGraph) {
  return get_dyn_prop<ValueType>(aProperty, aGraph)[aDesc];
}

namespace detail {
namespace {

template <typename ErasedDesc, typename Iterator>
class type_erased_graph_iterator
    : public boost::iterator_facade<
          type_erased_graph_iterator<ErasedDesc, Iterator>, ErasedDesc,
          typename std::iterator_traits<Iterator>::iterator_category,
          ErasedDesc> {
 public:
  using self = type_erased_graph_iterator<ErasedDesc, Iterator>;

  explicit type_erased_graph_iterator(Iterator aIt) : it(aIt) {}
  type_erased_graph_iterator() : type_erased_graph_iterator(Iterator()) {}

  // private:
  friend class boost::iterator_core_access;

  void increment() { ++it; }
  void decrement() { --it; }
  bool equal(const self& rhs) const { return this->it == rhs.it; }
  ErasedDesc dereference() const { return ErasedDesc(std::any(*it)); }

  void advance(std::ptrdiff_t i) { it += i; }
  std::ptrdiff_t distance_to(const self& rhs) const {
    return rhs.it - this->it;
  }

  Iterator it;
};

template <typename ErasedDesc, typename Iterator>
auto type_erase_graph_iterator(Iterator it) {
  return type_erased_graph_iterator<ErasedDesc, Iterator>(it);
}

}  // namespace
}  // namespace detail

template <typename Graph>
class type_erased_graph : public any_graph {
 public:
  using any_descriptor = any_graph::any_descriptor;
  using vertex_descriptor = any_graph::vertex_descriptor;
  using edge_descriptor = any_graph::edge_descriptor;

  using any_bundled = any_graph::any_bundled;
  using vertex_bundled = any_graph::vertex_bundled;
  using edge_bundled = any_graph::edge_bundled;

  using edge_range = any_graph::edge_range;
  using out_edge_range = any_graph::out_edge_range;
  using in_edge_range = any_graph::in_edge_range;

  using vertex_range = any_graph::vertex_range;

  using degree_size_type = any_graph::degree_size_type;
  using vertices_size_type = any_graph::vertices_size_type;
  using edges_size_type = any_graph::edges_size_type;

 protected:
  using original_graph_type = Graph;
  using real_graph_type = std::remove_cv_t<Graph>;
  using real_vertex_desc = graph_vertex_t<real_graph_type>;
  using real_edge_desc = graph_edge_t<real_graph_type>;

  real_graph_type* p_graph;

  void* get_property_by_ptr(std::string_view aProperty,
                            const std::any& /*aElement*/) const override {
    std::string err_message = "Unknown property: ";
    err_message += aProperty;
    throw std::invalid_argument(err_message.c_str());
  }

  std::any get_property_by_any(std::string_view aProperty,
                               const std::any& aElement) const override {
    std::string err_message = "Unknown property: ";
    err_message += aProperty;
    throw std::invalid_argument(err_message.c_str());
  }

  vertex_bundled get_vertex_bundled(
      const vertex_descriptor& aVertex) const override {
    return vertex_bundled(
        (*p_graph)[std::any_cast<real_vertex_desc>(aVertex.base)]);
  }

  edge_bundled get_edge_bundled(const edge_descriptor& aEdge) const override {
    return edge_bundled((*p_graph)[std::any_cast<real_edge_desc>(aEdge.base)]);
  }

  bool equal_vertex_descriptors(const vertex_descriptor& aU,
                                const vertex_descriptor& aV) const override {
    return (std::any_cast<real_vertex_desc>(aU.base) ==
            std::any_cast<real_vertex_desc>(aV.base));
  }

  bool equal_edge_descriptors(const edge_descriptor& aE,
                              const edge_descriptor& aF) const override {
    return (std::any_cast<real_edge_desc>(aE.base) ==
            std::any_cast<real_edge_desc>(aF.base));
  }

  edge_range edges_mem() const override {
    if constexpr (boost::is_edge_list_graph<real_graph_type>::value) {
      auto [first, last] = edges(*p_graph);
      return edge_range(
          detail::type_erase_graph_iterator<edge_descriptor>(first),
          detail::type_erase_graph_iterator<edge_descriptor>(last));
    }
    return {};
  }
  edges_size_type num_edges_mem() const override {
    if constexpr (boost::is_edge_list_graph<real_graph_type>::value) {
      return num_edges(*p_graph);
    }
    return 0;
  }

  vertex_descriptor source_mem(const edge_descriptor& aEdge) const override {
    if constexpr (boost::is_edge_list_graph<real_graph_type>::value ||
                  boost::is_incidence_graph<real_graph_type>::value) {
      return vertex_descriptor(
          std::any(source(std::any_cast<real_edge_desc>(aEdge), *p_graph)));
    }
    return vertex_descriptor(
        std::any(boost::graph_traits<real_graph_type>::null_vertex()));
  }
  vertex_descriptor target_mem(const edge_descriptor& aEdge) const override {
    if constexpr (boost::is_edge_list_graph<real_graph_type>::value ||
                  boost::is_incidence_graph<real_graph_type>::value) {
      return vertex_descriptor(
          std::any(target(std::any_cast<real_edge_desc>(aEdge), *p_graph)));
    }
    return vertex_descriptor(
        std::any(boost::graph_traits<real_graph_type>::null_vertex()));
  }

  std::pair<edge_descriptor, bool> edge_mem(
      const vertex_descriptor& aU, const vertex_descriptor& aV) const override {
    if constexpr (boost::is_adjacency_matrix<real_graph_type>::value) {
      auto [e, e_exists] = edge(std::any_cast<real_vertex_desc>(aU),
                                std::any_cast<real_vertex_desc>(aV), *p_graph);
      return {edge_descriptor(std::any(e)), e_exists};
    }
    return {edge_descriptor(), false};
  }

  out_edge_range out_edges_mem(
      const vertex_descriptor& aVertex) const override {
    if constexpr (boost::is_incidence_graph<real_graph_type>::value) {
      auto [first, last] =
          out_edges(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
      return out_edge_range(
          detail::type_erase_graph_iterator<edge_descriptor>(first),
          detail::type_erase_graph_iterator<edge_descriptor>(last));
    }
    return {};
  }

  degree_size_type out_degree_mem(
      const vertex_descriptor& aVertex) const override {
    if constexpr (boost::is_incidence_graph<real_graph_type>::value) {
      return out_degree(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
    }
    return 0;
  }

  in_edge_range in_edges_mem(const vertex_descriptor& aVertex) const override {
    if constexpr (boost::is_bidirectional_graph<real_graph_type>::value) {
      auto [first, last] =
          in_edges(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
      return in_edge_range(
          detail::type_erase_graph_iterator<edge_descriptor>(first),
          detail::type_erase_graph_iterator<edge_descriptor>(last));
    }
    return {};
  }

  degree_size_type in_degree_mem(
      const vertex_descriptor& aVertex) const override {
    if constexpr (boost::is_bidirectional_graph<real_graph_type>::value) {
      return in_degree(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
    }
    return 0;
  }

  degree_size_type degree_mem(const vertex_descriptor& aVertex) const override {
    if constexpr (boost::is_bidirectional_graph<real_graph_type>::value) {
      return degree(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
    }
    return 0;
  }

  vertex_range vertices_mem() const override {
    if constexpr (boost::is_vertex_list_graph<real_graph_type>::value) {
      auto [first, last] = vertices(*p_graph);
      return vertex_range(
          detail::type_erase_graph_iterator<vertex_descriptor>(first),
          detail::type_erase_graph_iterator<vertex_descriptor>(last));
    }
    return {};
  }

  vertices_size_type num_vertices_mem() const override {
    if constexpr (boost::is_vertex_list_graph<real_graph_type>::value) {
      return num_vertices(*p_graph);
    }
    return 0;
  }

  vertex_descriptor add_vertex_mem() override {
    if constexpr (boost::is_mutable_vertex_graph<real_graph_type>::value) {
      return vertex_descriptor(std::any(add_vertex(*p_graph)));
    }
    return vertex_descriptor(
        std::any(boost::graph_traits<real_graph_type>::null_vertex()));
  }

  void clear_vertex_mem(const vertex_descriptor& aVertex) override {
    if constexpr (boost::is_mutable_vertex_graph<real_graph_type>::value) {
      clear_vertex(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
    }
  }

  void remove_vertex_mem(const vertex_descriptor& aVertex) override {
    if constexpr (boost::is_mutable_vertex_graph<real_graph_type>::value) {
      remove_vertex(std::any_cast<real_vertex_desc>(aVertex), *p_graph);
    }
  }

  std::pair<edge_descriptor, bool> add_edge_mem(
      const vertex_descriptor& aU, const vertex_descriptor& aV) override {
    if constexpr (boost::is_mutable_edge_graph<real_graph_type>::value) {
      auto [e, e_exists] =
          add_edge(std::any_cast<real_vertex_desc>(aU),
                   std::any_cast<real_vertex_desc>(aV), *p_graph);
      return {edge_descriptor(std::any(e)), e_exists};
    }
    return {edge_descriptor(std::any(real_edge_desc())), false};
  }

  void remove_edge_mem(const edge_descriptor& aEdge) override {
    if constexpr (boost::is_mutable_edge_graph<real_graph_type>::value) {
      remove_edge(std::any_cast<real_edge_desc>(aEdge), *p_graph);
    }
  }

 public:
  real_graph_type& base() { return *p_graph; }
  const real_graph_type& base() const { return *p_graph; }

  explicit type_erased_graph(original_graph_type* aPGraph)
      : p_graph(const_cast<real_graph_type*>(aPGraph)) {}
};

}  // namespace ReaK::graph

#endif
